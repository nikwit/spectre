// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/KerrSchildGrVars.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave {
namespace {

void test_against_generic_solution(const double mass,
                                   const std::array<double, 3>& center) {
  MAKE_GENERATOR(gen);
  UniformCustomDistribution<double> dist{-10.0, 10.0};
  const size_t num_points = 50;
  // Random coordinates around the center, excluding the singularity by
  // offsetting one component
  auto coords = make_with_random_values<tnsr::I<DataVector, 3>>(
      make_not_null(&gen), make_not_null(&dist), DataVector(num_points));
  for (size_t i = 0; i < 3; ++i) {
    coords.get(i) += gsl::at(center, i);
  }
  get<0>(coords) += 25.0;

  using tags = typename System<3>::spacetime_tag_list;
  const gr::Solutions::KerrSchild solution{
      mass, {{0.0, 0.0, 0.0}}, center, {{0.0, 0.0, 0.0}}};
  const auto expected = solution.variables(coords, 0.0, tags{});

  Scalar<DataVector> lapse{};
  tnsr::i<DataVector, 3> deriv_lapse{};
  tnsr::I<DataVector, 3> shift{};
  tnsr::iJ<DataVector, 3> deriv_shift{};
  tnsr::ii<DataVector, 3> spatial_metric{};
  tnsr::II<DataVector, 3> inverse_spatial_metric{};
  tnsr::I<DataVector, 3> trace_spatial_christoffel{};
  Scalar<DataVector> trace_extrinsic_curvature{};
  zero_spin_kerr_schild_gr_vars(
      make_not_null(&lapse), make_not_null(&deriv_lapse), make_not_null(&shift),
      make_not_null(&deriv_shift), make_not_null(&spatial_metric),
      make_not_null(&inverse_spatial_metric),
      make_not_null(&trace_spatial_christoffel),
      make_not_null(&trace_extrinsic_curvature), coords, mass, center);

  Approx custom_approx = Approx::custom().epsilon(1e-12).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(
      lapse, get<gr::Tags::Lapse<DataVector>>(expected), custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      deriv_lapse,
      (get<::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                         Frame::Inertial>>(expected)),
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      shift, (get<gr::Tags::Shift<DataVector, 3>>(expected)), custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      deriv_shift,
      (get<::Tags::deriv<gr::Tags::Shift<DataVector, 3>, tmpl::size_t<3>,
                         Frame::Inertial>>(expected)),
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      spatial_metric, (get<gr::Tags::SpatialMetric<DataVector, 3>>(expected)),
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      inverse_spatial_metric,
      (get<gr::Tags::InverseSpatialMetric<DataVector, 3>>(expected)),
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      trace_spatial_christoffel,
      (get<gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, 3>>(
          expected)),
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      trace_extrinsic_curvature,
      get<gr::Tags::TraceExtrinsicCurvature<DataVector>>(expected),
      custom_approx);
}

SPECTRE_TEST_CASE("Unit.Evolution.Systems.CurvedScalarWave.KerrSchildGrVars",
                  "[Unit][Evolution]") {
  test_against_generic_solution(1.0, {{0.0, 0.0, 0.0}});
  test_against_generic_solution(0.5, {{0.0, 0.0, 0.0}});
  test_against_generic_solution(2.3, {{-1.5, 4.2, 0.7}});
}

}  // namespace
}  // namespace CurvedScalarWave
