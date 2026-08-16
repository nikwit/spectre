// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace {
using matcher_options = gh::Worldtube::MatcherConfig::options;
using expected_matcher_options =
    tmpl::list<gh::Worldtube::MatcherConfig::Mass,
               gh::Worldtube::MatcherConfig::FitInterval,
               gh::Worldtube::MatcherConfig::OrderOne,
               gh::Worldtube::MatcherConfig::ExcisionSphereName>;
static_assert(std::is_same_v<matcher_options, expected_matcher_options>);

struct SphereFields {
  tnsr::aa<DataVector, 3> metric;
  tnsr::aa<DataVector, 3> pi;
  tnsr::iaa<DataVector, 3> phi;
  Scalar<DataVector> gamma2;
  tnsr::I<DataVector, 3> coords;
  std::array<double, 3> center;
};

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.MatcherOptions",
    "[Unit][Evolution]") {
  const auto config = TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
      "Mass: 0.125\n"
      "FitInterval: 0.5\n"
      "OrderOne: Apply\n"
      "ExcisionSphereName: ExcisionSphereB\n");
  CHECK(config.mass == 0.125);
  CHECK(config.fit_interval == 0.5);
  CHECK(config.fit_l_max == 4);
  CHECK(config.fit_uplus);
  CHECK(config.fit_exact_frame);
  CHECK(config.fit_radial_derivative);
  CHECK(config.radial_derivative_weight == 1.);
  CHECK(config.order_one == gh::Worldtube::OrderOneMode::Apply);
  CHECK(config.fit_radial_index == 0);
  CHECK(config.excision_sphere_name == "ExcisionSphereB");
  test_serialization(config);
  const std::array<double, 15> unit_weights{
      {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}};
  CHECK(config.uplus_block_weights == unit_weights);

  CHECK(TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
            "Mass: 0.125\n"
            "FitInterval: 0.5\n"
            "OrderOne: Shadow\n"
            "ExcisionSphereName: ExcisionSphereB\n")
            .order_one == gh::Worldtube::OrderOneMode::Shadow);
  CHECK(TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
            "Mass: 0.125\n"
            "FitInterval: 0.5\n"
            "OrderOne: Off\n"
            "ExcisionSphereName: ExcisionSphereB\n")
            .order_one == gh::Worldtube::OrderOneMode::Off);

  CHECK_THROWS_WITH(
      TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
          "Mass: 0.125\n"
          "FitInterval: 0.5\n"
          "OrderOne: First\n"
          "ExcisionSphereName: ExcisionSphereB\n"),
      Catch::Matchers::ContainsSubstring("Must be Off, Shadow, or Apply"));

  CHECK_THROWS_WITH(TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
                        "Mass: 0.125\n"
                        "FitInterval: 0.5\n"
                        "OrderOne: Apply\n"
                        "ExcisionSphereName: ExcisionSphereB\n"
                        "RateOde: false\n"),
                    Catch::Matchers::ContainsSubstring("RateOde"));

  CHECK_THROWS_WITH(TestHelpers::test_creation<gh::Worldtube::MatcherConfig>(
                        "Mass: 0.125\n"
                        "FitInterval: 0.5\n"
                        "FitOrderOneShadow: true\n"
                        "ExcisionSphereName: ExcisionSphereB\n"),
                    Catch::Matchers::ContainsSubstring("FitOrderOneShadow"));
}

SphereFields fields_at(const ylm::Spherepack& ylm, const double time) {
  const auto& theta_phi = ylm.theta_phi_points();
  const size_t n_points = theta_phi[0].size();
  const std::array<double, 3> center{{0.2, -0.1, 0.3}};
  constexpr double radius = 3.2;

  SphereFields result{tnsr::aa<DataVector, 3>(n_points, 0.),
                      tnsr::aa<DataVector, 3>(n_points, 0.),
                      tnsr::iaa<DataVector, 3>(n_points, 0.),
                      Scalar<DataVector>(DataVector(n_points, 0.)),
                      tnsr::I<DataVector, 3>(n_points),
                      center};
  get<0>(result.coords) =
      center[0] + radius * sin(theta_phi[0]) * cos(theta_phi[1]);
  get<1>(result.coords) =
      center[1] + radius * sin(theta_phi[0]) * sin(theta_phi[1]);
  get<2>(result.coords) = center[2] + radius * cos(theta_phi[0]);

  const DataVector profile =
      0.4 * sin(0.6 * time) + 0.03 * get<0>(result.coords) -
      0.02 * get<1>(result.coords) + 0.01 * get<2>(result.coords);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      result.metric.get(a, b) =
          (a == b ? (a == 0 ? -1. : 1.) : 0.) +
          2.e-3 * static_cast<double>((a + 1) * (b + 1)) * profile;
      result.pi.get(a, b) = 0.01 * static_cast<double>(1 + a + 2 * b) *
                            (cos(0.3 * time) + 0.02 * get<0>(result.coords) *
                                                   get<2>(result.coords));
      for (size_t k = 0; k < 3; ++k) {
        result.phi.get(k, a, b) =
            0.004 * static_cast<double>((k + 1) * (a + 1) + b) *
            (sin(0.2 * time) + 0.01 * result.coords.get((k + 1) % 3));
      }
    }
  }
  get(result.gamma2) = 0.35 + 0.025 * time +
                       0.003 * get<0>(result.coords) * get<1>(result.coords);
  return result;
}

SphereFields boosted_fields_at(const ylm::Spherepack& ylm, const double time,
                               const std::array<double, 3>& velocity) {
  const auto& theta_phi = ylm.theta_phi_points();
  const size_t n_points = theta_phi[0].size();
  constexpr double radius = 2.;
  const std::array<double, 3> center{{0., 0., 0.}};
  SphereFields result{tnsr::aa<DataVector, 3>{},
                      tnsr::aa<DataVector, 3>{},
                      tnsr::iaa<DataVector, 3>{},
                      Scalar<DataVector>{DataVector(n_points, 0.1)},
                      tnsr::I<DataVector, 3>{n_points},
                      center};
  get<0>(result.coords) = radius * sin(theta_phi[0]) * cos(theta_phi[1]);
  get<1>(result.coords) = radius * sin(theta_phi[0]) * sin(theta_phi[1]);
  get<2>(result.coords) = radius * cos(theta_phi[0]);

  const std::array<double, gh::Worldtube::num_map_parameters> zero{};
  gh::Solutions::affine_map_model::boosted_evolved_variables(
      make_not_null(&result.metric), make_not_null(&result.pi),
      make_not_null(&result.phi), result.coords, time, 1., center, zero, zero,
      velocity);
  return result;
}

SphereFields boosted_fields_at(const ylm::Spherepack& ylm, const double time,
                               const double velocity) {
  return boosted_fields_at(ylm, time, {{0., 0., velocity}});
}

// Manufactured data: the exact pushforward of harmonic Schwarzschild through
// the finite frame map L = B(rapidity) S, evaluated on a sphere.
SphereFields exact_frame_fields_at(
    const ylm::Spherepack& ylm,
    const std::array<double, gh::Worldtube::num_map_parameters>& theta) {
  const auto& theta_phi = ylm.theta_phi_points();
  const size_t n_points = theta_phi[0].size();
  constexpr double radius = 2.;
  const std::array<double, 3> center{{0., 0., 0.}};
  SphereFields result{tnsr::aa<DataVector, 3>{},
                      tnsr::aa<DataVector, 3>{},
                      tnsr::iaa<DataVector, 3>{},
                      Scalar<DataVector>{DataVector(n_points, 0.1)},
                      tnsr::I<DataVector, 3>{n_points},
                      center};
  get<0>(result.coords) = radius * sin(theta_phi[0]) * cos(theta_phi[1]);
  get<1>(result.coords) = radius * sin(theta_phi[0]) * sin(theta_phi[1]);
  get<2>(result.coords) = radius * cos(theta_phi[0]);

  gh::Solutions::exact_frame::evolved_variables(
      make_not_null(&result.metric), make_not_null(&result.pi),
      make_not_null(&result.phi), result.coords, 0., 1., center, theta);
  return result;
}

SphereFields order_one_fields_at(
    const ylm::Spherepack& ylm, const double radius,
    const std::array<double, gh::Worldtube::num_map_parameters>& theta,
    const gh::Solutions::order_by_order_worldtube::AffineRates& rates) {
  const auto& theta_phi = ylm.theta_phi_points();
  const size_t n_points = theta_phi[0].size();
  const std::array<double, 3> center{{0., 0., 0.}};
  SphereFields result{tnsr::aa<DataVector, 3>{},
                      tnsr::aa<DataVector, 3>{},
                      tnsr::iaa<DataVector, 3>{},
                      Scalar<DataVector>{DataVector(n_points, 0.1)},
                      tnsr::I<DataVector, 3>{n_points},
                      center};
  get<0>(result.coords) = radius * sin(theta_phi[0]) * cos(theta_phi[1]);
  get<1>(result.coords) = radius * sin(theta_phi[0]) * sin(theta_phi[1]);
  get<2>(result.coords) = radius * cos(theta_phi[0]);
  gh::Solutions::order_by_order_worldtube::evolved_variables(
      make_not_null(&result.metric), make_not_null(&result.pi),
      make_not_null(&result.phi), result.coords, 0., 1., center,
      gh::Solutions::exact_frame::frame_map(theta), rates);
  return result;
}

gh::Solutions::exact_frame::FrameMatrix frame_tangent_from_rates(
    const gh::Solutions::exact_frame::FrameMatrix& frame,
    const gh::Solutions::order_by_order_worldtube::AffineRates& rates) {
  gh::Solutions::exact_frame::FrameMatrix generator{};
  generator[0][0] = rates[0];
  for (size_t i = 0; i < 3; ++i) {
    generator[0][i + 1] = rates[i + 1];
    generator[i + 1][0] = rates[i + 4];
  }
  static constexpr std::array<std::array<size_t, 2>, 6> symmetric{
      {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 2}}}};
  for (size_t parameter = 0; parameter < symmetric.size(); ++parameter) {
    const auto ij = symmetric[parameter];
    generator[ij[0] + 1][ij[1] + 1] = rates[parameter + 7];
    generator[ij[1] + 1][ij[0] + 1] = rates[parameter + 7];
  }
  gh::Solutions::exact_frame::FrameMatrix tangent{};
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        tangent[a][b] += frame[a][c] * generator[c][b];
      }
      tangent[a][b] /= frame[0][0];  // manufactured mass is one
    }
  }
  return tangent;
}

SphereFields evolving_order_one_fields_at(
    const ylm::Spherepack& ylm, const double radius, const double time,
    const std::array<double, gh::Worldtube::num_map_parameters>& theta,
    const gh::Solutions::order_by_order_worldtube::AffineRates& rates) {
  SphereFields result = order_one_fields_at(ylm, radius, theta, {});
  const auto frame = gh::Solutions::exact_frame::frame_map(theta);
  const auto tangent = frame_tangent_from_rates(frame, rates);
  auto time_frame = frame;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      time_frame[a][b] += time * tangent[a][b];
    }
  }
  gh::Solutions::exact_frame::evolved_variables(
      make_not_null(&result.metric), make_not_null(&result.pi),
      make_not_null(&result.phi), result.coords, time, 1., result.center,
      time_frame);
  tnsr::aa<DataVector, 3> metric_response{};
  tnsr::aa<DataVector, 3> pi_response{};
  tnsr::iaa<DataVector, 3> phi_response{};
  gh::Solutions::order_by_order_worldtube::
      affine_rate_evolved_variables_response(
          make_not_null(&metric_response), make_not_null(&pi_response),
          make_not_null(&phi_response), result.coords, time, 1., result.center,
          frame, rates);
  for (size_t storage = 0; storage < result.metric.size(); ++storage) {
    result.metric[storage] += metric_response[storage];
    result.pi[storage] += pi_response[storage];
  }
  for (size_t storage = 0; storage < result.phi.size(); ++storage) {
    result.phi[storage] += phi_response[storage];
  }
  return result;
}

gh::Worldtube::MatcherConfig exact_frame_config(
    const std::array<double, 3>& center) {
  gh::Worldtube::MatcherConfig config{};
  config.mass = 1.;
  config.center = center;
  config.fit_l_max = 4;
  config.fit_uplus = true;
  config.fit_exact_frame = true;
  return config;
}

template <typename TensorType>
TensorType centered_derivative(const TensorType& upper, const TensorType& lower,
                               const double step) {
  TensorType result{};
  for (size_t storage = 0; storage < upper.size(); ++storage) {
    result[storage] = (upper[storage] - lower[storage]) / (2. * step);
  }
  return result;
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.MovingCenter",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  const SphereFields sphere = fields_at(ylm, 0.);
  CHECK_ITERABLE_APPROX(
      gh::Worldtube::detail::worldtube_center(sphere.coords, ylm),
      sphere.center);

  gh::Worldtube::MatcherConfig config{};
  config.center = {{-8., 7., 6.}};  // ignored once live geometry is available
  config.centre_advection = true;
  gh::Worldtube::MapParameterData state{};
  state.valid = true;
  state.last_fit_time = 2.;
  state.p[4] = 0.3;
  state.p[5] = -0.1;
  state.p[6] = 0.2;
  state.center_offset = {{0.01, 0.02, -0.03}};
  state.worldtube_center_at_last_fit = {{1., 2., 3.}};
  state.worldtube_center = {{1.4, 1.8, 3.3}};
  state.worldtube_center_valid = true;

  // q advances with the independently fitted inertial hole velocity and
  // subtracts the numerical domain-map displacement.
  const std::array<double, 3> expected_offset{{-0.24, 0.17, -0.23}};
  CHECK_ITERABLE_APPROX(
      gh::Worldtube::detail::current_center_offset(config, state, 2.5),
      expected_offset);
  // Consequently the absolute model center is independent of how the
  // excision sphere moved: c_old + q_old + dt V.
  CHECK_ITERABLE_APPROX(gh::Worldtube::detail::model_center(config, state, 2.5),
                        (std::array<double, 3>{{1.16, 1.97, 3.07}}));

  config.centre_advection = false;
  CHECK_ITERABLE_APPROX(
      gh::Worldtube::detail::current_center_offset(config, state, 2.5),
      state.center_offset);
  CHECK_ITERABLE_APPROX(gh::Worldtube::detail::model_center(config, state, 2.5),
                        (std::array<double, 3>{{1.41, 1.82, 3.27}}));

  // In finite-background mode the exact boost, rather than the epsilon-order
  // qdot coefficient, owns inertial center advection.
  config.fit_bulk_boost = true;
  state.bulk_velocity = {{-0.2, 0.4, 0.1}};
  CHECK_ITERABLE_APPROX(
      gh::Worldtube::detail::current_center_offset(config, state, 2.5),
      (std::array<double, 3>{{-0.49, 0.42, -0.28}}));
  CHECK_ITERABLE_APPROX(gh::Worldtube::detail::model_center(config, state, 2.5),
                        (std::array<double, 3>{{0.91, 2.22, 3.02}}));

  // The exact-frame mode advects with the fitted coordinate centre velocity
  // V_c = L^i_0/L^0_0 (spec Eq. Z4), never with the boost parameter V.
  config.fit_bulk_boost = false;
  config.fit_exact_frame = true;
  state.exact_frame_valid = true;
  state.exact_frame_center_velocity = {{0.1, -0.2, 0.3}};
  CHECK_ITERABLE_APPROX(
      gh::Worldtube::detail::current_center_offset(config, state, 2.5),
      (std::array<double, 3>{{-0.34, 0.12, -0.18}}));
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.OrderOneBoundary",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  const SphereFields sphere = fields_at(ylm, 0.);
  gh::Worldtube::MatcherConfig config = exact_frame_config(sphere.center);
  gh::Worldtube::MapParameterData state{};
  state.exact_frame_theta = {{0.02, -0.01, 0.015, 0.003, -0.002, 0.001, 0.004,
                              0.002, -0.001, 0.0005, -0.003, 0.0015, 0.0025}};
  state.order_one_rates = {{1.e-5, -2.e-5, 1.5e-5, 0.5e-5, -1.e-5, 2.e-5,
                            -0.5e-5, 3.e-5, -2.e-5, 1.e-5, 2.5e-5, -1.5e-5,
                            0.5e-5, 0., 0., 0.}};
  state.order_one_valid = true;

  const auto frame =
      gh::Solutions::exact_frame::frame_map(state.exact_frame_theta);
  tnsr::aa<DataVector, 3> expected_metric{};
  tnsr::aa<DataVector, 3> expected_pi{};
  tnsr::iaa<DataVector, 3> expected_phi{};
  tnsr::aa<DataVector, 3> actual_metric{};
  tnsr::aa<DataVector, 3> actual_pi{};
  tnsr::iaa<DataVector, 3> actual_phi{};
  const auto check_tensors = [&expected_metric, &expected_pi, &expected_phi,
                              &actual_metric, &actual_pi, &actual_phi]() {
    for (size_t storage = 0; storage < expected_metric.size(); ++storage) {
      CHECK_ITERABLE_APPROX(actual_metric[storage], expected_metric[storage]);
      CHECK_ITERABLE_APPROX(actual_pi[storage], expected_pi[storage]);
    }
    for (size_t storage = 0; storage < expected_phi.size(); ++storage) {
      CHECK_ITERABLE_APPROX(actual_phi[storage], expected_phi[storage]);
    }
  };

  config.order_one = gh::Worldtube::OrderOneMode::Apply;
  gh::Solutions::order_by_order_worldtube::evolved_variables(
      make_not_null(&expected_metric), make_not_null(&expected_pi),
      make_not_null(&expected_phi), sphere.coords, 0., config.mass,
      sphere.center, frame, state.order_one_rates);
  gh::Worldtube::detail::exact_frame_boundary_evolved_variables(
      make_not_null(&actual_metric), make_not_null(&actual_pi),
      make_not_null(&actual_phi), sphere.coords, config, state, sphere.center);
  check_tensors();

  // Shadow and an invalid Apply solve both preserve the frame-only boundary.
  gh::Solutions::exact_frame::evolved_variables(
      make_not_null(&expected_metric), make_not_null(&expected_pi),
      make_not_null(&expected_phi), sphere.coords, 0., config.mass,
      sphere.center, frame);
  config.order_one = gh::Worldtube::OrderOneMode::Shadow;
  gh::Worldtube::detail::exact_frame_boundary_evolved_variables(
      make_not_null(&actual_metric), make_not_null(&actual_pi),
      make_not_null(&actual_phi), sphere.coords, config, state, sphere.center);
  check_tensors();
  config.order_one = gh::Worldtube::OrderOneMode::Apply;
  state.order_one_valid = false;
  gh::Worldtube::detail::exact_frame_boundary_evolved_variables(
      make_not_null(&actual_metric), make_not_null(&actual_pi),
      make_not_null(&actual_phi), sphere.coords, config, state, sphere.center);
  check_tensors();
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.WeightedModeFit",
    "[Unit][Evolution]") {
  const std::vector<std::vector<double>> columns{{1., 0., 1., 0.},
                                                 {0., 1., 0., 1.}};
  const std::vector<double> target{12., -3., 2., -3.};
  const std::vector<double> held_out_target{2., -3., 4., -6.};
  const std::vector<size_t> row_blocks{0, 1, 2, 3};
  std::array<double, 15> weights{};
  weights.fill(1.);

  const auto equal_weight_fit = gh::Worldtube::detail::fit_weighted_modes(
      target, held_out_target, columns, row_blocks, weights);
  CHECK(equal_weight_fit.coefficients[0] == approx(7.));
  CHECK(equal_weight_fit.coefficients[1] == approx(-3.));

  // Excluding only the corrupted first block recovers the known coefficients.
  // Its residual remains visible, and the same coefficients are independently
  // applied to the held-out target.
  weights[0] = 0.;
  const auto excluded_block_fit = gh::Worldtube::detail::fit_weighted_modes(
      target, held_out_target, columns, row_blocks, weights);
  CHECK_ITERABLE_APPROX(excluded_block_fit.coefficients,
                        (std::vector<double>{2., -3.}));
  CHECK_ITERABLE_APPROX(excluded_block_fit.fitted_residual,
                        (std::vector<double>{10., 0., 0., 0.}));
  CHECK_ITERABLE_APPROX(excluded_block_fit.held_out_residual,
                        (std::vector<double>{0., 0., 2., -3.}));
  CHECK(excluded_block_fit.condition_number == approx(std::sqrt(2.)));
}

// The perturbative truncation belongs to the model, not to the supplied
// fields. Verify that the value matcher fits the strict first-order response
// directly to a fully nonlinear, exactly Lorentz-boosted black-hole field.
SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "FirstOrderFitToNonlinearBoost",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  constexpr double velocity = 1.e-3;
  const SphereFields data = boosted_fields_at(ylm, 0., velocity);

  gh::Worldtube::MatcherConfig config{};
  config.mass = 1.;
  config.center = data.center;
  config.fit_l_max = 4;
  config.fit_trace_strain = true;
  config.fit_velocity = true;
  config.fit_center_offset = false;
  config.centre_advection = true;
  const std::array<double, gh::Worldtube::num_map_parameters> p_start{};
  const std::array<double, 3> center_offset_start{};
  const auto fit = gh::Worldtube::fit_map_parameters(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      p_start, center_offset_start, -1., 0.);

  CAPTURE(fit.residual_initial, fit.residual_final, fit.iterations, fit.p);
  CHECK(fit.residual_final < fit.residual_initial);
  CHECK(fit.residual_final < fit.baseline_residual);
  CHECK(fit.minus_residual_final < fit.minus_baseline_residual);
  CHECK(fit.metric_residual_final < fit.metric_baseline_residual);
  CHECK(fit.phi_residual_final < fit.phi_baseline_residual);
  CHECK(std::isfinite(fit.condition_number));
  CHECK(fit.p[3] == approx(velocity).epsilon(2.e-2));
  CHECK(fit.p[6] == approx(velocity).epsilon(2.e-2));
  // The exact target contains gamma-1 and longitudinal-contraction terms at
  // O(v^2). They are allowed to remain in the residual (or be partly absorbed
  // by first-order coefficients); the input field is never pre-linearized.
  CHECK(fit.residual_final < 2.e-5);

  // A known velocity is a first-order pin by itself. It must not acquire the
  // former (1 + qdot0) resummation, whose product is second order.
  config.fit_velocity = false;
  config.center_velocity = {{0., 0., velocity}};
  const auto pinned_fit = gh::Worldtube::fit_map_parameters(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      p_start, center_offset_start, -1., 0.);
  CHECK(pinned_fit.p[4] == 0.);
  CHECK(pinned_fit.p[5] == 0.);
  CHECK(pinned_fit.p[6] == velocity);

  // The value solve must actually apply the configured block mask while
  // continuing to report the omitted block without that solve weight.
  config.uplus_block_weights[6] = 0.;  // C, ell=1
  const auto no_c1_fit = gh::Worldtube::fit_map_parameters(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      p_start, center_offset_start, -1., 0.);
  CAPTURE(no_c1_fit.condition_number, pinned_fit.condition_number,
          no_c1_fit.block_closure[6], no_c1_fit.block_minus_closure[6]);
  CHECK(std::isfinite(no_c1_fit.condition_number));
  CHECK(std::abs(no_c1_fit.condition_number - pinned_fit.condition_number) >
        1.e-6);
  CHECK(std::isfinite(no_c1_fit.block_closure[6]));
  CHECK(std::isfinite(no_c1_fit.block_minus_closure[6]));
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "FiniteBulkBoostFit",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  const std::array<double, 3> expected_velocity{{0.08, -0.04, 0.06}};
  const SphereFields data = boosted_fields_at(ylm, 0., expected_velocity);

  gh::Worldtube::MatcherConfig config{};
  config.mass = 1.;
  config.center = data.center;
  config.fit_l_max = 4;
  config.fit_bulk_boost = true;
  config.fit_velocity = false;
  config.center_velocity = {{0., 0., 0.}};
  config.fit_center_offset = false;
  config.centre_advection = true;
  const std::array<double, gh::Worldtube::num_map_parameters> p_start{};
  const std::array<double, 3> center_offset_start{};
  const auto fit = gh::Worldtube::fit_map_parameters(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      p_start, center_offset_start, -1., 0.);

  CAPTURE(fit.bulk_velocity, fit.bulk_residual_initial, fit.bulk_residual_final,
          fit.residual_final, fit.p);
  CHECK_ITERABLE_APPROX(fit.bulk_velocity, expected_velocity);
  CHECK(fit.bulk_residual_final < 1.e-10);
  CHECK(fit.bulk_residual_final < fit.bulk_residual_initial);
  CHECK(fit.residual_final < 1.e-10);
  for (const double parameter : fit.p) {
    CHECK(std::abs(parameter) < 1.e-9);
  }
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "CharacteristicDerivativeConvention",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  constexpr double time = 1.3;
  constexpr double step = 1.e-5;
  constexpr size_t fit_l_max = 4;
  const SphereFields lower = fields_at(ylm, time - step);
  const SphereFields middle = fields_at(ylm, time);
  const SphereFields upper = fields_at(ylm, time + step);

  const auto dt_metric = centered_derivative(upper.metric, lower.metric, step);
  const auto dt_pi = centered_derivative(upper.pi, lower.pi, step);
  const auto dt_phi = centered_derivative(upper.phi, lower.phi, step);
  const auto dt_gamma2 = centered_derivative(upper.gamma2, lower.gamma2, step);

  for (const double normal_sign : {-1., 1.}) {
    const std::vector<double> lower_modes =
        gh::Worldtube::detail::characteristic_gauge_modes(
            lower.metric, lower.pi, lower.phi, lower.gamma2, lower.coords,
            lower.center, ylm, fit_l_max, normal_sign);
    const std::vector<double> upper_modes =
        gh::Worldtube::detail::characteristic_gauge_modes(
            upper.metric, upper.pi, upper.phi, upper.gamma2, upper.coords,
            upper.center, ylm, fit_l_max, normal_sign);
    std::vector<double> expected(lower_modes.size());
    for (size_t i = 0; i < expected.size(); ++i) {
      expected[i] = (upper_modes[i] - lower_modes[i]) / (2. * step);
    }

    const std::vector<double> result =
        gh::Worldtube::detail::characteristic_gauge_time_derivative_modes(
            middle.metric, middle.pi, middle.phi, dt_metric, dt_pi, dt_phi,
            middle.gamma2, dt_gamma2, middle.coords, middle.center, ylm,
            fit_l_max, normal_sign);
    CHECK_ITERABLE_CUSTOM_APPROX(result, expected,
                                 Approx::custom().epsilon(2.e-8).scale(1.));

    // Negative control: the fixed-sphere history has a varying gamma2, so its
    // derivative may not be omitted.
    const Scalar<DataVector> zero_dt_gamma2(
        DataVector(get(dt_gamma2).size(), 0.));
    const auto check_is_detectably_different =
        [&expected](const std::vector<double>& incomplete) {
          double error_squared = 0.;
          for (size_t i = 0; i < expected.size(); ++i) {
            error_squared += square(incomplete[i] - expected[i]);
          }
          CHECK(std::sqrt(error_squared) > 1.e-5);
        };
    check_is_detectably_different(
        gh::Worldtube::detail::characteristic_gauge_time_derivative_modes(
            middle.metric, middle.pi, middle.phi, dt_metric, dt_pi, dt_phi,
            middle.gamma2, zero_dt_gamma2, middle.coords, middle.center, ylm,
            fit_l_max, normal_sign));
  }
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "BoostedAccelerationTransport",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  constexpr double velocity = 1.e-3;
  constexpr double step = 1.e-4;
  const SphereFields lower = boosted_fields_at(ylm, -step, velocity);
  const SphereFields middle = boosted_fields_at(ylm, 0., velocity);
  const SphereFields upper = boosted_fields_at(ylm, step, velocity);
  const auto dt_metric = centered_derivative(upper.metric, lower.metric, step);
  const auto dt_pi = centered_derivative(upper.pi, lower.pi, step);
  const auto dt_phi = centered_derivative(upper.phi, lower.phi, step);
  const Scalar<DataVector> dt_gamma2(DataVector(get(middle.gamma2).size(), 0.));

  gh::Worldtube::MatcherConfig config{};
  config.mass = 1.;
  config.center = middle.center;
  config.center_velocity = {{0., 0., velocity}};
  config.fit_l_max = 4;
  config.centre_advection = true;
  config.uplus_block_weights.fill(1.);
  std::array<double, gh::Worldtube::num_map_parameters> p{};
  std::array<double, gh::Worldtube::num_map_parameters> pdot{};
  p[3] = velocity;
  p[6] = velocity;

  const auto with_transport =
      gh::Worldtube::fit_map_parameter_accelerations_uplus(
          middle.metric, middle.pi, middle.phi, dt_metric, dt_pi, dt_phi,
          middle.gamma2, dt_gamma2, p, pdot, middle.coords, ylm, config, true);
  config.centre_advection = false;
  const auto without_transport =
      gh::Worldtube::fit_map_parameter_accelerations_uplus(
          middle.metric, middle.pi, middle.phi, dt_metric, dt_pi, dt_phi,
          middle.gamma2, dt_gamma2, p, pdot, middle.coords, ylm, config, true);

  CAPTURE(with_transport.residual_initial, with_transport.residual_final,
          with_transport.minus_residual_initial,
          with_transport.minus_residual_final,
          with_transport.parameter_derivative_norm,
          without_transport.residual_initial,
          without_transport.parameter_derivative_norm);
  CHECK(with_transport.residual_initial <
        0.01 * without_transport.residual_initial);
  CHECK(with_transport.parameter_derivative_norm <
        0.01 * without_transport.parameter_derivative_norm);
  CHECK(with_transport.residual_final <= with_transport.residual_initial);
  CHECK(with_transport.minus_residual_initial < 1.e-5);
  CHECK(with_transport.minus_residual_final < 1.e-5);

  // Every single-block exclusion must remain a well-defined solve.  Running
  // this test with Catch's -s flag also provides a compact manufactured-data
  // leave-one-block-out audit, indexed [A0..A4,C0..C4,V0..V4].
  config.centre_advection = true;
  for (size_t block = 0; block < config.uplus_block_weights.size(); ++block) {
    config.uplus_block_weights.fill(1.);
    gsl::at(config.uplus_block_weights, block) = 0.;
    const auto excluded = gh::Worldtube::fit_map_parameter_accelerations_uplus(
        middle.metric, middle.pi, middle.phi, dt_metric, dt_pi, dt_phi,
        middle.gamma2, dt_gamma2, p, pdot, middle.coords, ylm, config, true);
    CAPTURE(block, excluded.residual_final, excluded.minus_residual_final,
            excluded.condition_number, excluded.parameter_derivative_norm);
    CHECK((std::isfinite(excluded.residual_final) and
           std::isfinite(excluded.minus_residual_final) and
           std::isfinite(excluded.condition_number)));
  }
}

// Contract C4 of the zeroth-order brief (spec test Z-1, exact recovery):
// manufactured data = the exact pushforward with generic theta; the
// exact-frame fit recovers every parameter family, separately and jointly,
// to fit precision from a cold start, and drives both characteristics --
// including the held-out one -- to roundoff.
SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.ExactFrameRecovery",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  using ThetaArray = std::array<double, gh::Worldtube::num_map_parameters>;
  // the generic frame of the oracle fixtures
  const ThetaArray generic{{0.08031150429463316, -0.04015575214731658,
                            0.06023362822097486, 0.05, 0.004, -0.002, 0.003,
                            -0.05, 0.012, -0.007, -0.04, 0.009, -0.055}};
  ThetaArray boost_only{};
  ThetaArray clock_only{};
  ThetaArray simultaneity_only{};
  ThetaArray strain_only{};
  for (size_t a = 0; a < 3; ++a) {
    gsl::at(boost_only, a) = gsl::at(generic, a);
  }
  clock_only[3] = generic[3];
  for (size_t a = 4; a < 7; ++a) {
    gsl::at(simultaneity_only, a) = gsl::at(generic, a);
  }
  for (size_t a = 7; a < 13; ++a) {
    gsl::at(strain_only, a) = gsl::at(generic, a);
  }

  for (const auto& expected :
       {boost_only, clock_only, simultaneity_only, strain_only, generic}) {
    const SphereFields data = exact_frame_fields_at(ylm, expected);
    const auto config = exact_frame_config(data.center);
    const ThetaArray cold_start{};
    const auto fit = gh::Worldtube::fit_exact_frame_parameters(
        data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
        cold_start, {{0., 0., 0.}});
    CAPTURE(expected, fit.exact_frame_theta, fit.residual_initial,
            fit.residual_final, fit.baseline_residual, fit.iterations);
    CHECK(fit.residual_final < 1.e-10);
    for (size_t a = 0; a < gh::Worldtube::num_map_parameters; ++a) {
      CHECK(std::abs(gsl::at(fit.exact_frame_theta, a) - gsl::at(expected, a)) <
            1.e-7);
    }
    // one frame closes the held-out channel too: nothing was traded
    CHECK(fit.minus_residual_final < 1.e-8);
    CHECK(std::isfinite(fit.condition_number));
  }
}

// Contract C5 (spec test Z-2, the strain ladder at fixed physical
// amplitude): data with s_ij = -Phi delta_ij, s0 = +Phi. The old linear
// matcher floors at the quadratic remainder ~Phi^2 -- doubling Phi
// quadruples its converged residual -- while the exact-frame fit sits at
// roundoff across the ladder. This is the decisive test of the bookkeeping
// claim, the symmetric-sector twin of the boost ladder.
SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "ExactFrameStrainLadder",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  const std::array<double, 3> amplitudes{{0.0125, 0.025, 0.05}};
  std::array<double, 3> linear_residuals{};
  std::array<double, 3> linear_baselines{};
  for (size_t rung = 0; rung < 3; ++rung) {
    const double amplitude = gsl::at(amplitudes, rung);
    std::array<double, gh::Worldtube::num_map_parameters> theta{};
    theta[3] = amplitude;
    theta[7] = -amplitude;
    theta[10] = -amplitude;
    theta[12] = -amplitude;
    const SphereFields data = exact_frame_fields_at(ylm, theta);

    const auto config = exact_frame_config(data.center);
    const auto fit = gh::Worldtube::fit_exact_frame_parameters(
        data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
        {}, {{0., 0., 0.}});
    CAPTURE(amplitude, fit.residual_final, fit.baseline_residual,
            fit.iterations, fit.exact_frame_theta);
    CHECK(fit.residual_final < 1.e-10);
    CHECK(std::abs(fit.exact_frame_theta[3] - amplitude) < 1.e-8);
    CHECK(std::abs(fit.exact_frame_theta[12] + amplitude) < 1.e-8);

    // The old linear matcher on the same data and channel, all 13 free. On
    // this spherically symmetric single-sphere data the sampled u+ gauge
    // modes have fewer independent rows than the 13 coefficients, so the
    // linear fit can zero its own fitted channel exactly -- the strain
    // aliasing the exact-frame model removes. The quadratic floor therefore
    // shows where no coefficient choice can hide it: the full-field
    // collocation metric closure.
    gh::Worldtube::MatcherConfig linear_config{};
    linear_config.mass = 1.;
    linear_config.center = data.center;
    linear_config.fit_l_max = 4;
    linear_config.fit_uplus = true;
    linear_config.fit_trace_strain = true;
    linear_config.fit_velocity = true;
    const std::array<double, gh::Worldtube::num_map_parameters> p_start{};
    const auto linear_fit = gh::Worldtube::fit_map_parameters(
        data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm,
        linear_config, p_start, {{0., 0., 0.}}, -1., 0.);
    CAPTURE(linear_fit.residual_final, linear_fit.baseline_residual,
            linear_fit.metric_residual_final,
            linear_fit.metric_baseline_residual, linear_fit.p);
    gsl::at(linear_residuals, rung) = linear_fit.metric_residual_final;
    gsl::at(linear_baselines, rung) = linear_fit.metric_baseline_residual;
    CHECK(fit.metric_residual_final < 1.e-10);
    CHECK(linear_fit.metric_residual_final >
          1.e3 * std::max(fit.metric_residual_final, 1.e-14));
  }
  // the linear metric floor is quadratic in the amplitude ...
  const double first_ratio = linear_residuals[1] / linear_residuals[0];
  const double second_ratio = linear_residuals[2] / linear_residuals[1];
  CAPTURE(linear_residuals, linear_baselines, first_ratio, second_ratio);
  CHECK(first_ratio > 3.);
  CHECK(first_ratio < 5.5);
  CHECK(second_ratio > 3.);
  CHECK(second_ratio < 5.5);
  // ... while the physical content itself only doubles
  CHECK(linear_baselines[1] / linear_baselines[0] < 2.5);
  CHECK(linear_baselines[2] / linear_baselines[1] < 2.5);
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube."
    "OrderOneAffineRateFit",
    "[Unit][Evolution]") {
  const ylm::Spherepack ylm{5, 5};
  const std::array<double, gh::Worldtube::num_map_parameters> theta{
      {0.03, -0.02, 0.01, 0.012, -0.008, 0.005, -0.004, 0.006, 0.002, -0.003,
       -0.005, 0.004, 0.007}};
  const gh::Solutions::order_by_order_worldtube::AffineRates expected{
      {2.e-4, -1.e-4, 1.5e-4, -0.8e-4, 1.2e-4, -1.1e-4, 0.9e-4, 0.7e-4, -0.6e-4,
       0.5e-4, -0.4e-4, 0.3e-4, -0.2e-4, 0., 0., 0.}};
  constexpr double derivative_step = 2.e-6;
  gh::Worldtube::RadialDerivativeStencil radial_stencil{};
  radial_stencil.fit_shell = 1;
  for (const double radius : {3., 3.2, 3.4}) {
    const SphereFields middle =
        evolving_order_one_fields_at(ylm, radius, 0., theta, expected);
    const SphereFields upper = evolving_order_one_fields_at(
        ylm, radius, derivative_step, theta, expected);
    const SphereFields lower = evolving_order_one_fields_at(
        ylm, radius, -derivative_step, theta, expected);
    radial_stencil.metric.push_back(middle.metric);
    radial_stencil.pi.push_back(middle.pi);
    radial_stencil.phi.push_back(middle.phi);
    radial_stencil.dt_metric.push_back(
        centered_derivative(upper.metric, lower.metric, derivative_step));
    radial_stencil.dt_pi.push_back(
        centered_derivative(upper.pi, lower.pi, derivative_step));
    radial_stencil.dt_phi.push_back(
        centered_derivative(upper.phi, lower.phi, derivative_step));
    radial_stencil.gamma2.push_back(middle.gamma2);
    radial_stencil.coords.push_back(middle.coords);
  }
  const SphereFields data =
      evolving_order_one_fields_at(ylm, 3.2, 0., theta, expected);
  const auto config = exact_frame_config(data.center);
  const auto fit = gh::Worldtube::fit_order_one_affine_rates(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      theta, {{0., 0., 0.}}, radial_stencil);
  CAPTURE(fit.condition_number, fit.time_residual_initial,
          fit.time_residual_final, fit.radial_time_residual_initial,
          fit.radial_time_residual_final, fit.minus_residual_initial,
          fit.minus_residual_final);
  CHECK(fit.valid);
  CHECK(fit.condition_number < 1.e9);
  CHECK(fit.time_residual_final < 2.e-5 * fit.time_residual_initial);
  CHECK(fit.radial_time_residual_final <
        2.e-5 * fit.radial_time_residual_initial);
  for (size_t parameter = 0; parameter < 13; ++parameter) {
    CAPTURE(parameter, expected[parameter], fit.rates[parameter]);
    CHECK(fit.rates[parameter] ==
          Approx::custom().epsilon(2.e-4).scale(1.e-10)(expected[parameter]));
  }
  for (size_t parameter = 13; parameter < 16; ++parameter) {
    CHECK(fit.rates[parameter] == 0.);
  }

  // The full wrapper must feed the projected rate response back into the
  // clean order-zero sectors twice and then re-solve the rates on the final
  // frame. Starting away from the manufactured frame exercises both blocks.
  auto theta_start = theta;
  theta_start[0] += 2.e-3;
  theta_start[7] -= 1.e-3;
  const auto iterated = gh::Worldtube::fit_iterated_order_zero_one(
      data.metric, data.pi, data.phi, data.gamma2, data.coords, ylm, config,
      theta_start, {{0., 0., 0.}}, radial_stencil);
  CAPTURE(iterated.order_zero.exact_frame_theta, iterated.order_one.rates,
          iterated.order_one.final_frame_step_norm);
  CHECK(iterated.order_one.valid);
  CHECK(iterated.order_one.alternations == 2);
  double initial_frame_error = 0.;
  double final_frame_error = 0.;
  for (size_t parameter = 0; parameter < theta.size(); ++parameter) {
    initial_frame_error =
        std::max(initial_frame_error,
                 std::abs(theta_start[parameter] - theta[parameter]));
    final_frame_error =
        std::max(final_frame_error,
                 std::abs(iterated.order_zero.exact_frame_theta[parameter] -
                          theta[parameter]));
    CHECK(std::abs(iterated.order_zero.exact_frame_theta[parameter] -
                   theta[parameter]) < 2.e-5);
  }
  CHECK(final_frame_error < 0.01 * initial_frame_error);
  for (size_t parameter = 0; parameter < 13; ++parameter) {
    CHECK(std::abs(iterated.order_one.rates[parameter] - expected[parameter]) <
          1.e-6);
  }
  for (size_t parameter = 13; parameter < 16; ++parameter) {
    CHECK(iterated.order_one.rates[parameter] == 0.);
  }
}
