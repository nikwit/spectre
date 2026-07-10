// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <random>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/UpdateAcceleration.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {
static constexpr size_t Dim = 3;

using variables_tag = ::Tags::Variables<
    tmpl::list<Tags::EvolvedPosition<Dim>, Tags::EvolvedVelocity<Dim>,
               Tags::Psi0, Tags::dtPsi0>>;
using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;

// wraps the Psi0 ODE part of `UpdateAcceleration` for comparison with the
// python implementation. The evolved variables Psi0 and dtPsi0 are packed
// into a tensor so they can be tested on the python side.
void time_derivative_wrapper(
    const gsl::not_null<Scalar<double>*> dt_psi0,
    const gsl::not_null<Scalar<double>*> dt2_psi0,
    const Scalar<double>& psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
    const tnsr::ii<double, Dim, Frame::Grid>& psi_quadrupole,
    const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole,
    const tnsr::AA<double, Dim, Frame::Grid>& inverse_spacetime_metric,
    const tnsr::A<double, Dim, Frame::Grid>& trace_christoffel,
    const tnsr::i<double, Dim, Frame::Grid>& evolved_vars_tnsr,
    const Scalar<double> wt_radius) {
  typename variables_tag::type evolved_vars(1, 0.);
  get(get<Tags::Psi0>(evolved_vars))[0] = get<0>(evolved_vars_tnsr);
  get(get<Tags::dtPsi0>(evolved_vars))[0] = get<1>(evolved_vars_tnsr);
  typename dt_variables_tag::type dt_evolved_vars(1);
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{tnsr::I<double, Dim>(0.),
                                                    tnsr::I<double, Dim>(0.)};
  const tnsr::I<double, Dim, Frame::Inertial> geodesic_acc(0.);
  UpdateAcceleration::apply(
      make_not_null(&dt_evolved_vars), evolved_vars, pos_vel, geodesic_acc,
      psi_monopole, psi_dipole, psi_quadrupole, dt_psi_dipole,
      inverse_spacetime_metric, trace_christoffel, 2, get(wt_radius), 0);
  get(*dt_psi0) = get(get<::Tags::dt<Tags::Psi0>>(dt_evolved_vars))[0];
  get(*dt2_psi0) = get(get<::Tags::dt<Tags::dtPsi0>>(dt_evolved_vars))[0];
}

// checks the Psi0 ODE at expansion order 2 against the python implementation
void test_psi0_ode() {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions"};
  pypp::check_with_random_values<1>(&time_derivative_wrapper, "TimeDerivatives",
                                    {"dt_psi0", "dt2_psi0"}, {{{0.1, 1.}}}, 1,
                                    1.e-12);
}

// checks the geodesic evolution of position and velocity and that Psi0 is
// not evolved below expansion order 2. Also checks that requesting
// acceleration iterations (i.e. the scalar self-force) is an error.
void test_geodesic_and_errors() {
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(-1., 1.);
  const DataVector used_for_size(1);
  auto dt_evolved_vars = make_with_random_values<dt_variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  auto evolved_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  const auto geodesic_acc = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto vel = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto pos = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{{pos, vel}};
  const auto psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto psi_dipole =
      make_with_random_values<tnsr::i<double, Dim, Frame::Grid>>(
          make_not_null(&gen), dist, 1);
  const auto psi_quadrupole =
      make_with_random_values<tnsr::ii<double, Dim, Frame::Grid>>(
          make_not_null(&gen), dist, 1);
  const auto dt_psi_dipole =
      make_with_random_values<tnsr::i<double, Dim, Frame::Grid>>(
          make_not_null(&gen), dist, 1);
  const auto inverse_spacetime_metric =
      make_with_random_values<tnsr::AA<double, Dim, Frame::Grid>>(
          make_not_null(&gen), dist, 1);
  const auto trace_christoffel =
      make_with_random_values<tnsr::A<double, Dim, Frame::Grid>>(
          make_not_null(&gen), dist, 1);
  const double wt_radius = 1.6;

  for (size_t expansion_order = 0; expansion_order < 2; ++expansion_order) {
    CAPTURE(expansion_order);
    auto box = db::create<db::AddSimpleTags<
        dt_variables_tag, variables_tag, Tags::ParticlePositionVelocity<Dim>,
        Tags::GeodesicAcceleration<Dim>,
        Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>,
        Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>,
        Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Grid>,
        Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                             Frame::Grid>,
        gr::Tags::InverseSpacetimeMetric<double, Dim, Frame::Grid>,
        gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim, Frame::Grid>,
        Tags::ExpansionOrder, Tags::WorldtubeRadius, Tags::MaxIterations>>(
        dt_evolved_vars, evolved_vars, pos_vel, geodesic_acc, psi_monopole,
        psi_dipole, psi_quadrupole, dt_psi_dipole, inverse_spacetime_metric,
        trace_christoffel, expansion_order, wt_radius, size_t(0));

    db::mutate_apply<UpdateAcceleration>(make_not_null(&box));
    const auto& dt_vars = db::get<dt_variables_tag>(box);
    for (size_t i = 0; i < Dim; ++i) {
      CHECK(get<::Tags::dt<Tags::EvolvedPosition<Dim>>>(dt_vars).get(i)[0] ==
            vel.get(i));
      CHECK(get<::Tags::dt<Tags::EvolvedVelocity<Dim>>>(dt_vars).get(i)[0] ==
            geodesic_acc.get(i));
    }
    // Psi0 is only evolved at expansion order 2
    CHECK(get(get<::Tags::dt<Tags::Psi0>>(dt_vars))[0] == 0.);
    CHECK(get(get<::Tags::dt<Tags::dtPsi0>>(dt_vars))[0] == 0.);

    db::mutate<Tags::MaxIterations>(
        [](const gsl::not_null<size_t*> max_iterations) {
          *max_iterations = 1;
        },
        make_not_null(&box));
    CHECK_THROWS_WITH(db::mutate_apply<UpdateAcceleration>(make_not_null(&box)),
                      Catch::Matchers::ContainsSubstring(
                          "The scalar self-force is not supported"));
  }
}

SPECTRE_TEST_CASE("Unit.Evolution.Systems.CSW.Worldtube.UpdateAcceleration",
                  "[Unit][Evolution]") {
  test_psi0_ode();
  test_geodesic_and_errors();
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
