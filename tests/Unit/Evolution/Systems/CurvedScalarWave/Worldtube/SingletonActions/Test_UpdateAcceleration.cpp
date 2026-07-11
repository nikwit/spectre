// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <random>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SelfForce.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/UpdateAcceleration.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {
constexpr size_t Dim = 3;

using variables_tag = UpdateAcceleration::variables_tag;
using dt_variables_tag = UpdateAcceleration::dt_variables_tag;

using worldtube_tags = db::AddSimpleTags<
    dt_variables_tag, variables_tag, Tags::ParticlePositionVelocity<Dim>,
    Tags::BackgroundQuantities<Dim>, Tags::GeodesicAcceleration<Dim>,
    Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Inertial>,
    Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim,
                         Frame::Inertial>,
    Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Inertial>,
    Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                         Frame::Inertial>,
    Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Inertial>,
    Tags::Charge, Tags::Mass, Tags::MaxIterations, ::Tags::Time,
    Tags::SelfForceTurnOnTime, Tags::SelfForceTurnOnInterval,
    Tags::ExpansionOrder, Tags::WorldtubeRadius>;

// Wraps the Psi0 ODE part of `UpdateAcceleration` for comparison with the
// python implementation. With `max_iterations` set to 0 the self-force
// vanishes, so the full particle acceleration entering the ODE is the
// geodesic acceleration which corresponds to the `particle_acceleration`
// argument on the python side. The evolved variables Psi0 and dtPsi0 are
// packed into a rank-1 tensor so they can be passed to the python side.
void time_derivative_wrapper(
    const gsl::not_null<Scalar<double>*> dt_psi0_result,
    const gsl::not_null<Scalar<double>*> dt2_psi0_result,
    const Scalar<double>& psi_monopole,
    const tnsr::i<double, Dim, Frame::Inertial>& psi_dipole,
    const tnsr::ii<double, Dim, Frame::Inertial>& psi_quadrupole,
    const tnsr::i<double, Dim, Frame::Inertial>& dt_psi_dipole,
    const tnsr::AA<double, Dim>& inverse_spacetime_metric,
    const tnsr::A<double, Dim>& trace_christoffel,
    const tnsr::i<double, Dim, Frame::Inertial>& evolved_vars_tnsr,
    const Scalar<double>& wt_radius,
    const tnsr::I<double, Dim, Frame::Inertial>& particle_velocity,
    const tnsr::I<double, Dim, Frame::Inertial>& particle_acceleration) {
  typename variables_tag::type evolved_vars(1, 0.);
  get(get<Tags::Psi0>(evolved_vars))[0] = get<0>(evolved_vars_tnsr);
  get(get<Tags::dtPsi0>(evolved_vars))[0] = get<1>(evolved_vars_tnsr);
  typename dt_variables_tag::type dt_evolved_vars(1);
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{tnsr::I<double, Dim>(0.),
                                                    particle_velocity};
  // only the inverse metric and the trace of the Christoffel symbol enter
  // the ODE
  const typename Tags::BackgroundQuantities<Dim>::type background{
      tnsr::aa<double, Dim>(0.), inverse_spacetime_metric,
      tnsr::Abb<double, Dim>(0.), trace_christoffel, Scalar<double>(1.)};
  const Scalar<double> dt_psi_monopole(0.);
  UpdateAcceleration::apply(make_not_null(&dt_evolved_vars), evolved_vars,
                            pos_vel, background, particle_acceleration,
                            psi_monopole, dt_psi_monopole, psi_dipole,
                            dt_psi_dipole, psi_quadrupole, 0.1, std::nullopt, 0,
                            0., std::nullopt, std::nullopt, 2, get(wt_radius));
  get(*dt_psi0_result) = get(get<::Tags::dt<Tags::Psi0>>(dt_evolved_vars))[0];
  get(*dt2_psi0_result) =
      get(get<::Tags::dt<Tags::dtPsi0>>(dt_evolved_vars))[0];
}

// checks the Psi0 ODE at expansion order 2 against the python implementation
void test_psi0_ode() {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions"};
  pypp::check_with_random_values<1>(&time_derivative_wrapper, "TimeDerivatives",
                                    {"dt_psi0", "dt2_psi0"}, {{{0.1, 1.}}}, 1,
                                    1.e-12);
}

// checks that position and velocity are evolved geodesically when the
// self-force is disabled and that Psi0 is not evolved below expansion order 2
void test_geodesic_acceleration() {
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(-1., 1.);
  const DataVector used_for_size(1);
  auto dt_evolved_vars = make_with_random_values<dt_variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  const auto evolved_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  const auto geodesic_acc = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto vel = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto pos = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{{pos, vel}};
  const auto metric = make_with_random_values<tnsr::aa<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto inverse_metric = make_with_random_values<tnsr::AA<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto christoffel = make_with_random_values<tnsr::Abb<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto trace_christoffel = make_with_random_values<tnsr::A<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dilation =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const typename Tags::BackgroundQuantities<Dim>::type background_quantities{
      metric, inverse_metric, christoffel, trace_christoffel, dilation};
  const auto psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto dt_psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dt_psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto psi_quadrupole = make_with_random_values<tnsr::ii<double, Dim>>(
      make_not_null(&gen), dist, 1);

  const double charge = 0.2;
  const double mass = 0.1;
  const double time = 105.;
  const double turn_on_time = 100.;
  const double turn_on_interval = 10.;
  const double wt_radius = 1.6;

  for (size_t expansion_order = 0; expansion_order < 2; ++expansion_order) {
    CAPTURE(expansion_order);
    auto box = db::create<worldtube_tags>(
        dt_evolved_vars, evolved_vars, pos_vel, background_quantities,
        geodesic_acc, psi_monopole, dt_psi_monopole, psi_dipole, dt_psi_dipole,
        psi_quadrupole, charge, std::make_optional(mass), size_t{0}, time,
        std::make_optional(turn_on_time), std::make_optional(turn_on_interval),
        expansion_order, wt_radius);

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
  }
}

// checks that the acceleration includes the rolled-on self-force when
// iterations are requested and that the same full acceleration enters the
// Psi0 ODE at expansion order 2
void test_self_force_acceleration() {
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(-1., 1.);
  const DataVector used_for_size(1);
  auto dt_evolved_vars = make_with_random_values<dt_variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  const auto evolved_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  const auto geodesic_acc = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto vel = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto pos = make_with_random_values<tnsr::I<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{{pos, vel}};
  const auto metric = make_with_random_values<tnsr::aa<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto inverse_metric = make_with_random_values<tnsr::AA<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto christoffel = make_with_random_values<tnsr::Abb<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto trace_christoffel = make_with_random_values<tnsr::A<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dilation =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const typename Tags::BackgroundQuantities<Dim>::type background_quantities{
      metric, inverse_metric, christoffel, trace_christoffel, dilation};
  const auto psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto dt_psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dt_psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto psi_quadrupole = make_with_random_values<tnsr::ii<double, Dim>>(
      make_not_null(&gen), dist, 1);

  const double charge = 0.2;
  const double mass = 0.1;
  const double time = 105.;
  const double turn_on_time = 100.;
  const double turn_on_interval = 10.;
  const double wt_radius = 1.6;
  const size_t expansion_order = 2;
  const size_t max_iterations = 1;

  auto box = db::create<worldtube_tags>(
      dt_evolved_vars, evolved_vars, pos_vel, background_quantities,
      geodesic_acc, psi_monopole, dt_psi_monopole, psi_dipole, dt_psi_dipole,
      psi_quadrupole, charge, std::make_optional(mass), max_iterations, time,
      std::make_optional(turn_on_time), std::make_optional(turn_on_interval),
      expansion_order, wt_radius);

  db::mutate_apply<UpdateAcceleration>(make_not_null(&box));

  const double roll_on =
      turn_on_function(time - turn_on_time, turn_on_interval);
  const double evolved_mass = mass - charge * get(psi_monopole);
  const auto self_force_acc =
      self_force_acceleration(dt_psi_monopole, psi_dipole, vel, charge,
                              evolved_mass, inverse_metric, dilation);

  const Approx local_approx = Approx::custom().epsilon(1e-12);
  const auto& dt_vars = db::get<dt_variables_tag>(box);
  tnsr::I<double, Dim> full_acceleration{};
  for (size_t i = 0; i < Dim; ++i) {
    CHECK(get<::Tags::dt<Tags::EvolvedPosition<Dim>>>(dt_vars).get(i)[0] ==
          local_approx(vel.get(i)));
    CHECK(get<::Tags::dt<Tags::EvolvedVelocity<Dim>>>(dt_vars).get(i)[0] ==
          local_approx(geodesic_acc.get(i) + roll_on * self_force_acc.get(i)));
    full_acceleration.get(i) =
        geodesic_acc.get(i) + roll_on * self_force_acc.get(i);
  }

  // The Psi0 ODE must use the full particle acceleration including the
  // rolled-on self-force: applying the mutator with the self-force disabled
  // but the geodesic acceleration replaced by the full acceleration has to
  // give the same time derivatives of Psi0 and dtPsi0.
  typename dt_variables_tag::type expected_dt_vars(1);
  UpdateAcceleration::apply(
      make_not_null(&expected_dt_vars), evolved_vars, pos_vel,
      background_quantities, full_acceleration, psi_monopole, dt_psi_monopole,
      psi_dipole, dt_psi_dipole, psi_quadrupole, charge,
      std::make_optional(mass), 0, time, std::make_optional(turn_on_time),
      std::make_optional(turn_on_interval), expansion_order, wt_radius);
  CHECK(get(get<::Tags::dt<Tags::Psi0>>(dt_vars))[0] ==
        local_approx(get(get<::Tags::dt<Tags::Psi0>>(expected_dt_vars))[0]));
  CHECK(get(get<::Tags::dt<Tags::dtPsi0>>(dt_vars))[0] ==
        local_approx(get(get<::Tags::dt<Tags::dtPsi0>>(expected_dt_vars))[0]));
}

// checks that for a static particle the Psi0 ODE at expansion order 2
// reduces to the static limit which we compute inline here
void test_static_ode() {
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<> dist(-1., 1.);
  std::uniform_real_distribution<> dist_positive(0.5, 2.);
  const DataVector used_for_size(1);
  const auto evolved_vars = make_with_random_values<variables_tag::type>(
      make_not_null(&gen), dist, used_for_size);
  typename dt_variables_tag::type dt_evolved_vars(1);
  // static particle: zero velocity and zero acceleration
  const std::array<tnsr::I<double, Dim>, 2> pos_vel{
      make_with_random_values<tnsr::I<double, Dim>>(make_not_null(&gen), dist,
                                                    1),
      tnsr::I<double, Dim>(0.)};
  const tnsr::I<double, Dim> geodesic_acc(0.);
  const auto metric = make_with_random_values<tnsr::aa<double, Dim>>(
      make_not_null(&gen), dist, 1);
  // ensure the 00-component is not close to zero as we divide by it
  const auto inverse_metric = make_with_random_values<tnsr::AA<double, Dim>>(
      make_not_null(&gen), dist_positive, 1);
  const auto christoffel = make_with_random_values<tnsr::Abb<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto trace_christoffel = make_with_random_values<tnsr::A<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dilation =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const typename Tags::BackgroundQuantities<Dim>::type background_quantities{
      metric, inverse_metric, christoffel, trace_christoffel, dilation};
  const auto psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto dt_psi_monopole =
      make_with_random_values<Scalar<double>>(make_not_null(&gen), dist, 1);
  const auto psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto dt_psi_dipole = make_with_random_values<tnsr::i<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const auto psi_quadrupole = make_with_random_values<tnsr::ii<double, Dim>>(
      make_not_null(&gen), dist, 1);
  const double charge = 0.2;
  const double mass = 0.1;
  const double wt_radius = 1.6;

  UpdateAcceleration::apply(make_not_null(&dt_evolved_vars), evolved_vars,
                            pos_vel, background_quantities, geodesic_acc,
                            psi_monopole, dt_psi_monopole, psi_dipole,
                            dt_psi_dipole, psi_quadrupole, charge,
                            std::make_optional(mass), 0, 0., std::nullopt,
                            std::nullopt, 2, wt_radius);

  const double psi0 = get(get<Tags::Psi0>(evolved_vars))[0];
  const double dt_psi0 = get(get<Tags::dtPsi0>(evolved_vars))[0];

  // full second-order coefficient: STF part plus trace
  tnsr::ii<double, Dim> psi_ij_full = psi_quadrupole;
  const double trace_part = (get(psi_monopole) - psi0) / square(wt_radius);
  for (size_t i = 0; i < Dim; ++i) {
    psi_ij_full.get(i, i) += trace_part;
  }
  // static limit of the ODE: with zero velocity and acceleration only the
  // terms independent of the particle motion remain
  double expected_dt2_psi0 = get<0>(trace_christoffel) * dt_psi0;
  for (size_t i = 0; i < Dim; ++i) {
    expected_dt2_psi0 -=
        2. * inverse_metric.get(0, i + 1) * dt_psi_dipole.get(i);
    expected_dt2_psi0 += trace_christoffel.get(i + 1) * psi_dipole.get(i);
    for (size_t j = 0; j < Dim; ++j) {
      expected_dt2_psi0 -=
          2. * inverse_metric.get(i + 1, j + 1) * psi_ij_full.get(i, j);
    }
  }
  expected_dt2_psi0 /= inverse_metric.get(0, 0);

  const Approx local_approx = Approx::custom().epsilon(1e-12);
  for (size_t i = 0; i < Dim; ++i) {
    CHECK(get<::Tags::dt<Tags::EvolvedPosition<Dim>>>(dt_evolved_vars)
              .get(i)[0] == 0.);
    CHECK(get<::Tags::dt<Tags::EvolvedVelocity<Dim>>>(dt_evolved_vars)
              .get(i)[0] == 0.);
  }
  CHECK(get(get<::Tags::dt<Tags::Psi0>>(dt_evolved_vars))[0] == dt_psi0);
  CHECK(get(get<::Tags::dt<Tags::dtPsi0>>(dt_evolved_vars))[0] ==
        local_approx(expected_dt2_psi0));
}

SPECTRE_TEST_CASE("Unit.Evolution.Systems.CSW.Worldtube.UpdateAcceleration",
                  "[Unit][Evolution]") {
  test_psi0_ode();
  test_geodesic_acceleration();
  test_self_force_acceleration();
  test_static_ode();
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
