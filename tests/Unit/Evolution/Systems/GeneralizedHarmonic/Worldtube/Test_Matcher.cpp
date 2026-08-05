// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace {
struct SphereFields {
  tnsr::aa<DataVector, 3> metric;
  tnsr::aa<DataVector, 3> pi;
  tnsr::iaa<DataVector, 3> phi;
  Scalar<DataVector> gamma2;
  tnsr::I<DataVector, 3> coords;
  std::array<double, 3> center;
};

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
