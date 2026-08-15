// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/HarmonicWorldtubeModel.hpp"
#include "Utilities/Gsl.hpp"

namespace {
using gh::Solutions::order_by_order_worldtube::AffineRates;

const AffineRates rates{{0.07, -0.03, 0.02, 0.04, 0.05, -0.06, 0.08, 0.01,
                         -0.02, 0.03, -0.04, 0.05, -0.01, 0.025, -0.035,
                         0.045}};

void check_against_matrix(
    const tnsr::AA<DataVector, 3>& actual,
    const std::array<std::array<double, 4>, 4>& expected) {
  const Approx custom = Approx::custom().epsilon(2.e-13).scale(1.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      CHECK(actual.get(a, b)[0] == custom(gsl::at(gsl::at(expected, a), b)));
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "HarmonicWorldtubeModel.LocalRateResponse",
    "[PointwiseFunctions][Unit]") {
  tnsr::I<DataVector, 3> point(size_t{1});
  get<0>(point)[0] = 3.1;
  get<1>(point)[0] = -2.2;
  get<2>(point)[0] = 1.7;
  tnsr::AA<DataVector, 3> response{};
  gh::Solutions::order_by_order_worldtube::
      local_affine_rate_inverse_metric_response(make_not_null(&response), point,
                                                1., rates);

  // Independent reference values from the audited q=8 Python response
  // functions (SciPy's spence implementation), not from this C++ evaluator.
  check_against_matrix(
      response, {{{{0.05646917206113843, 0.297296865666648, -0.4070909897705419,
                    -0.1382463933990438}},
                  {{0.297296865666648, 0.10472264628104741,
                    -0.19101910457249496, 0.21433705481505505}},
                  {{-0.4070909897705419, -0.19101910457249496,
                    0.14902859730177065, -0.2819198337553173}},
                  {{-0.1382463933990438, 0.21433705481505505,
                    -0.2819198337553173, 0.12344500694812144}}}});

  AffineRates minus_rates{};
  for (size_t i = 0; i < rates.size(); ++i) {
    gsl::at(minus_rates, i) = -gsl::at(rates, i);
  }
  tnsr::AA<DataVector, 3> minus_response{};
  gh::Solutions::order_by_order_worldtube::
      local_affine_rate_inverse_metric_response(make_not_null(&minus_response),
                                                point, 1., minus_rates);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      CHECK_ITERABLE_APPROX(response.get(a, b), -minus_response.get(a, b));
    }
  }
}

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "HarmonicWorldtubeModel.EvolvedVariableResponse",
    "[PointwiseFunctions][Unit]") {
  const gh::Solutions::exact_frame::FrameMatrix frame{{
      {{1.03, 0.08, -0.03, 0.02}},
      {{0.06, 1.02, 0.01, -0.02}},
      {{-0.01, 0.03, 0.98, 0.04}},
      {{0.04, -0.02, 0.05, 1.01}},
  }};
  const std::array<double, 3> center{{0.2, -0.1, 0.3}};
  constexpr double time = 0.37;
  tnsr::I<DataVector, 3> point(size_t{1});
  get<0>(point)[0] = 4.4;
  get<1>(point)[0] = -1.4;
  get<2>(point)[0] = 2.8;

  const auto metric_from_inverse =
      [&frame, &center](const tnsr::I<DataVector, 3>& local_point,
                        const double local_time,
                        const AffineRates& local_rates) {
        tnsr::AA<DataVector, 3> inverse{};
        gh::Solutions::order_by_order_worldtube::inverse_metric(
            make_not_null(&inverse), local_point, local_time, 1., center, frame,
            local_rates);
        return determinant_and_inverse(inverse).second;
      };

  tnsr::aa<DataVector, 3> metric_response{};
  tnsr::aa<DataVector, 3> pi_response{};
  tnsr::iaa<DataVector, 3> phi_response{};
  gh::Solutions::order_by_order_worldtube::
      affine_rate_evolved_variables_response(
          make_not_null(&metric_response), make_not_null(&pi_response),
          make_not_null(&phi_response), point, time, 1., center, frame, rates);

  // Differentiate the exact matrix inverse of G_0 +/- epsilon delta G. This
  // is independent of the production covariant-sandwich linearization.
  constexpr double rate_step = 1.e-5;
  AffineRates plus_rates{};
  AffineRates minus_rates{};
  for (size_t i = 0; i < rates.size(); ++i) {
    plus_rates[i] = rate_step * rates[i];
    minus_rates[i] = -rate_step * rates[i];
  }
  const auto plus_metric = metric_from_inverse(point, time, plus_rates);
  const auto minus_metric = metric_from_inverse(point, time, minus_rates);
  const Approx tight = Approx::custom().epsilon(2.e-10).scale(1.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      const DataVector numerical_response =
          (plus_metric.get(a, b) - minus_metric.get(a, b)) / (2. * rate_step);
      CHECK_ITERABLE_CUSTOM_APPROX(metric_response.get(a, b),
                                   numerical_response, tight);
    }
  }

  // Differentiate the value-level metric response in space and time. The
  // value response is already checked against the independent Python model,
  // so this directly audits the forward analytic derivatives.
  const auto metric_response_at = [&frame, &center](
                                      const tnsr::I<DataVector, 3>& local_point,
                                      const double local_time) {
    tnsr::aa<DataVector, 3> local_metric_response{};
    tnsr::aa<DataVector, 3> unused_pi{};
    tnsr::iaa<DataVector, 3> unused_phi{};
    gh::Solutions::order_by_order_worldtube::
        affine_rate_evolved_variables_response(
            make_not_null(&local_metric_response), make_not_null(&unused_pi),
            make_not_null(&unused_phi), local_point, local_time, 1., center,
            frame, rates);
    return local_metric_response;
  };
  constexpr double coordinate_step = 2.e-5;
  tnsr::aa<DataVector, 3> numerical_dt_metric_response(size_t{1}, 0.);
  const auto time_plus = metric_response_at(point, time + coordinate_step);
  const auto time_minus = metric_response_at(point, time - coordinate_step);
  const Approx derivative_tolerance = Approx::custom().epsilon(3.e-8).scale(1.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      numerical_dt_metric_response.get(a, b) =
          (time_plus.get(a, b) - time_minus.get(a, b)) / (2. * coordinate_step);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    auto point_plus = point;
    auto point_minus = point;
    point_plus.get(i)[0] += coordinate_step;
    point_minus.get(i)[0] -= coordinate_step;
    const auto spatial_plus = metric_response_at(point_plus, time);
    const auto spatial_minus = metric_response_at(point_minus, time);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        const DataVector numerical_phi =
            (spatial_plus.get(a, b) - spatial_minus.get(a, b)) /
            (2. * coordinate_step);
        CHECK_ITERABLE_CUSTOM_APPROX(phi_response.get(i, a, b), numerical_phi,
                                     derivative_tolerance);
      }
    }
  }

  // Reconstruct the linearized Pi from the independently differentiated
  // delta g and the exact background 3+1 fields.
  tnsr::aa<DataVector, 3> background_metric{};
  tnsr::aa<DataVector, 3> background_pi{};
  tnsr::iaa<DataVector, 3> background_phi{};
  gh::Solutions::exact_frame::evolved_variables(
      make_not_null(&background_metric), make_not_null(&background_pi),
      make_not_null(&background_phi), point, time, 1., center, frame);
  const auto background_inverse =
      determinant_and_inverse(background_metric).second;
  tnsr::AA<DataVector, 3> inverse_response{};
  gh::Solutions::order_by_order_worldtube::affine_rate_inverse_metric_response(
      make_not_null(&inverse_response), point, time, 1., center, frame, rates);
  const DataVector lapse = 1. / sqrt(-get<0, 0>(background_inverse));
  const DataVector lapse_response =
      0.5 * cube(lapse) * get<0, 0>(inverse_response);
  std::array<DataVector, 3> shift{};
  std::array<DataVector, 3> shift_response{};
  for (size_t i = 0; i < 3; ++i) {
    shift[i] =
        -background_inverse.get(0, i + 1) / get<0, 0>(background_inverse);
    shift_response[i] =
        -inverse_response.get(0, i + 1) / get<0, 0>(background_inverse) +
        background_inverse.get(0, i + 1) * get<0, 0>(inverse_response) /
            square(get<0, 0>(background_inverse));
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      DataVector expected_pi = -numerical_dt_metric_response.get(a, b) / lapse -
                               lapse_response / lapse * background_pi.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        expected_pi += (shift_response[i] * background_phi.get(i, a, b) +
                        shift[i] * phi_response.get(i, a, b)) /
                       lapse;
      }
      CHECK_ITERABLE_CUSTOM_APPROX(pi_response.get(a, b), expected_pi,
                                   derivative_tolerance);
    }
  }
}

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "HarmonicWorldtubeModel.FiniteFrameComposition",
    "[PointwiseFunctions][Unit]") {
  const gh::Solutions::exact_frame::FrameMatrix frame{{
      {{1.03, 0.08, -0.03, 0.02}},
      {{0.06, 1.02, 0.01, -0.02}},
      {{-0.01, 0.03, 0.98, 0.04}},
      {{0.04, -0.02, 0.05, 1.01}},
  }};
  tnsr::I<DataVector, 3> point(size_t{1});
  get<0>(point)[0] = 4.2;
  get<1>(point)[0] = -1.3;
  get<2>(point)[0] = 2.5;
  tnsr::AA<DataVector, 3> response{};
  gh::Solutions::order_by_order_worldtube::affine_rate_inverse_metric_response(
      make_not_null(&response), point, 0.37, 1., {{0., 0., 0.}}, frame, rates);
  check_against_matrix(response,
                       {{{{0.08098971473817321, 0.3336723612281542,
                           -0.34001490474724183, -0.1901333566361708}},
                         {{0.33367236122815425, 0.190562811517908,
                           -0.18859680664906187, 0.25089467683413935}},
                         {{-0.3400149047472418, -0.18859680664906184,
                           0.07641818830283627, -0.26400978466961783}},
                         {{-0.19013335663617084, 0.25089467683413935,
                           -0.26400978466961783, 0.16399743958879046}}}});

  tnsr::AA<DataVector, 3> background{};
  gh::Solutions::exact_frame::inverse_metric(make_not_null(&background), point,
                                             0.37, 1., {{0., 0., 0.}}, frame);
  tnsr::AA<DataVector, 3> full_model{};
  gh::Solutions::order_by_order_worldtube::inverse_metric(
      make_not_null(&full_model), point, 0.37, 1., {{0., 0., 0.}}, frame,
      rates);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      CHECK_ITERABLE_APPROX(full_model.get(a, b),
                            background.get(a, b) + response.get(a, b));
    }
  }
}
