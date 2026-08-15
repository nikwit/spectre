// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/HarmonicWorldtubeModel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace gh::Solutions::order_by_order_worldtube {
namespace {

// Li_2(y), for 0 <= y < 1. The power series is only used at y <= 1/2;
// reflection maps the remaining interval back to the same rapidly convergent
// series. This is sufficient for Spence(x) = Li_2(1-x), because the model is
// evaluated outside the harmonic horizon and hence x=(rho+M)/(2M) > 1.
double dilogarithm_on_unit_interval(const double y) {
  ASSERT(y >= 0. and y < 1.,
         "Dilogarithm argument must be in [0,1). Got " << y);
  if (y == 0.) {
    return 0.;
  }
  if (y > 0.5) {
    constexpr double pi_squared_over_six =
        1.64493406684822643647241516664602519;
    return pi_squared_over_six - log(y) * log1p(-y) -
           dilogarithm_on_unit_interval(1. - y);
  }
  double sum = 0.;
  double y_to_k = y;
  for (size_t k = 1; k < 10000; ++k) {
    const double term = y_to_k / square(static_cast<double>(k));
    sum += term;
    if (abs(term) <= 2.e-16 * std::max(1., abs(sum))) {
      break;
    }
    y_to_k *= y;
  }
  return sum;
}

// scipy.special.spence(x) = Li_2(1-x). For x >= 1, transform the negative
// argument z=1-x to y=z/(z-1)=1-1/x in [0,1):
// Li_2(z) = -Li_2(y) - log(x)^2/2.
double spence_for_x_at_least_one(const double x) {
  ASSERT(x >= 1., "The exterior rate model requires x >= 1. Got " << x);
  const double log_x = log(x);
  return -dilogarithm_on_unit_interval(1. - 1. / x) - 0.5 * square(log_x);
}

template <size_t NumberOfDerivatives>
struct FirstDerivative {
  double value = 0.;
  std::array<double, NumberOfDerivatives> derivative{};

  FirstDerivative() = default;
  // NOLINTNEXTLINE(google-explicit-constructor)
  FirstDerivative(const double local_value) : value(local_value) {}

  FirstDerivative& operator+=(const FirstDerivative& rhs) {
    value += rhs.value;
    for (size_t i = 0; i < NumberOfDerivatives; ++i) {
      derivative[i] += rhs.derivative[i];
    }
    return *this;
  }

  friend FirstDerivative operator+(FirstDerivative lhs,
                                   const FirstDerivative& rhs) {
    lhs += rhs;
    return lhs;
  }

  friend FirstDerivative operator-(const FirstDerivative& lhs,
                                   const FirstDerivative& rhs) {
    FirstDerivative result{lhs.value - rhs.value};
    for (size_t i = 0; i < NumberOfDerivatives; ++i) {
      result.derivative[i] = lhs.derivative[i] - rhs.derivative[i];
    }
    return result;
  }

  friend FirstDerivative operator-(const FirstDerivative& operand) {
    FirstDerivative result{-operand.value};
    for (size_t i = 0; i < NumberOfDerivatives; ++i) {
      result.derivative[i] = -operand.derivative[i];
    }
    return result;
  }

  friend FirstDerivative operator*(const FirstDerivative& lhs,
                                   const FirstDerivative& rhs) {
    FirstDerivative result{lhs.value * rhs.value};
    for (size_t i = 0; i < NumberOfDerivatives; ++i) {
      result.derivative[i] =
          lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
    }
    return result;
  }

  friend FirstDerivative operator/(const FirstDerivative& lhs,
                                   const FirstDerivative& rhs) {
    FirstDerivative result{lhs.value / rhs.value};
    for (size_t i = 0; i < NumberOfDerivatives; ++i) {
      result.derivative[i] =
          (lhs.derivative[i] * rhs.value - lhs.value * rhs.derivative[i]) /
          square(rhs.value);
    }
    return result;
  }
};

double scalar_sqrt(const double value) { return sqrt(value); }

template <size_t NumberOfDerivatives>
FirstDerivative<NumberOfDerivatives> scalar_sqrt(
    const FirstDerivative<NumberOfDerivatives>& operand) {
  FirstDerivative<NumberOfDerivatives> result{sqrt(operand.value)};
  for (size_t i = 0; i < NumberOfDerivatives; ++i) {
    result.derivative[i] = 0.5 * operand.derivative[i] / result.value;
  }
  return result;
}

double scalar_log(const double value) { return log(value); }

template <size_t NumberOfDerivatives>
FirstDerivative<NumberOfDerivatives> scalar_log(
    const FirstDerivative<NumberOfDerivatives>& operand) {
  FirstDerivative<NumberOfDerivatives> result{log(operand.value)};
  for (size_t i = 0; i < NumberOfDerivatives; ++i) {
    result.derivative[i] = operand.derivative[i] / operand.value;
  }
  return result;
}

template <size_t NumberOfDerivatives>
FirstDerivative<NumberOfDerivatives> spence_for_x_at_least_one(
    const FirstDerivative<NumberOfDerivatives>& x) {
  FirstDerivative<NumberOfDerivatives> result{
      spence_for_x_at_least_one(x.value)};
  // d Li_2(1-x) / dx = log(x) / (1-x). The exterior chart has x>1.
  const double epsilon = x.value - 1.;
  const double radial_derivative =
      std::abs(epsilon) < 1.e-8 ? -1. + 0.5 * epsilon - square(epsilon) / 3.
                                : log(x.value) / (1. - x.value);
  for (size_t i = 0; i < NumberOfDerivatives; ++i) {
    result.derivative[i] = radial_derivative * x.derivative[i];
  }
  return result;
}

template <typename T>
std::array<std::array<T, 4>, 4> local_response_at_point(
    const std::array<T, 3>& local_coords, const double mass,
    const AffineRates& rates) {
  T radius_squared{0.};
  for (size_t i = 0; i < 3; ++i) {
    radius_squared += local_coords[i] * local_coords[i];
  }
  const T radius = scalar_sqrt(radius_squared);
  std::array<T, 3> direction{};
  for (size_t i = 0; i < 3; ++i) {
    direction[i] = local_coords[i] / radius;
  }

  const T areal_radius = radius + mass;
  const T f = 1. - 2. * mass / areal_radius;
  const T big_f = 1. + 2. * mass / areal_radius;
  const T h = 4. * square(mass) / (areal_radius * areal_radius);
  const T h_prime =
      -8. * square(mass) / (areal_radius * areal_radius * areal_radius);
  const T g_perpendicular = radius * radius / (areal_radius * areal_radius);
  const T g_perpendicular_prime =
      2. * mass * radius / (areal_radius * areal_radius * areal_radius);
  const T g_parallel = -square(mass) / (areal_radius * areal_radius);
  const T g_parallel_prime =
      2. * square(mass) / (areal_radius * areal_radius * areal_radius);
  const T background_tt = f * big_f * big_f - 2. * big_f;
  const T background_tt_prime =
      2. * mass / (areal_radius * areal_radius) +
      8. * square(mass) / (areal_radius * areal_radius * areal_radius) +
      24. * cube(mass) /
          (areal_radius * areal_radius * areal_radius * areal_radius);

  const T x = areal_radius / (2. * mass);
  const T tortoise_radius = areal_radius + 2. * mass * scalar_log(x);
  const T tortoise_prime = 1. + 2. * mass / areal_radius;
  const T dhesi_profile =
      -4. * square(mass) * spence_for_x_at_least_one(x) +
      ((3. * radius + 14. * mass) * tortoise_radius -
       (radius * radius + 8. * mass * radius + 25. * square(mass))) /
          3.;
  const T dhesi_profile_prime =
      4. * square(mass) * scalar_log(x) / (radius - mass) +
      (3. * tortoise_radius + (3. * radius + 14. * mass) * tortoise_prime -
       (2. * radius + 8. * mass)) /
          3.;
  const T acceleration_profile =
      dhesi_profile - 0.5 * tortoise_radius * tortoise_radius;
  const T acceleration_profile_prime =
      dhesi_profile_prime - tortoise_radius * tortoise_prime;

  const double time_acceleration = rates[0];
  const std::array<double, 3> time_gradient_rate{
      {rates[1], rates[2], rates[3]}};
  const std::array<double, 3> spatial_acceleration{
      {rates[4], rates[5], rates[6]}};
  const std::array<std::array<double, 3>, 3> spatial_strain_rate{{
      {{rates[7], rates[8] + rates[13], rates[9] + rates[14]}},
      {{rates[8] - rates[13], rates[10], rates[11] + rates[15]}},
      {{rates[9] - rates[14], rates[11] - rates[15], rates[12]}},
  }};

  T radial_acceleration{0.};
  T time_gradient_radial{0.};
  std::array<T, 3> strain_direction{};
  for (size_t i = 0; i < 3; ++i) {
    radial_acceleration += spatial_acceleration[i] * direction[i];
    time_gradient_radial += time_gradient_rate[i] * direction[i];
    strain_direction[i] = T{0.};
    for (size_t j = 0; j < 3; ++j) {
      strain_direction[i] += spatial_strain_rate[i][j] * direction[j];
    }
  }
  T strain_radial{0.};
  for (size_t i = 0; i < 3; ++i) {
    strain_radial += direction[i] * strain_direction[i];
  }
  const T strain_profile = -2. * (radius - 2. * mass);
  constexpr double strain_profile_prime = -2.;

  std::array<std::array<T, 4>, 4> result{};
  result[0][0] =
      2. * h * acceleration_profile_prime * time_acceleration -
      acceleration_profile * background_tt_prime * radial_acceleration +
      2. * background_tt * (radius / mass) * time_gradient_radial +
      2. * h * strain_profile_prime * time_gradient_radial -
      strain_profile * background_tt_prime * strain_radial;
  for (size_t i = 0; i < 3; ++i) {
    const T acceleration_ti =
        h * acceleration_profile_prime * spatial_acceleration[i] +
        f * acceleration_profile_prime * time_acceleration * direction[i] -
        acceleration_profile * h_prime * radial_acceleration * direction[i] -
        acceleration_profile * h / radius *
            (spatial_acceleration[i] - radial_acceleration * direction[i]);
    const T time_gradient_ti =
        h * (radius / mass) * time_gradient_radial * direction[i] +
        f * strain_profile_prime * time_gradient_radial * direction[i] +
        g_perpendicular * strain_profile / radius *
            (time_gradient_rate[i] - time_gradient_radial * direction[i]);
    const T strain_ti =
        background_tt * (radius / mass) * strain_direction[i] +
        h * strain_profile_prime * strain_direction[i] -
        strain_profile * h_prime * strain_radial * direction[i] -
        strain_profile * h / radius *
            (strain_direction[i] - strain_radial * direction[i]);
    result[0][i + 1] = acceleration_ti + time_gradient_ti + strain_ti;
    result[i + 1][0] = result[0][i + 1];

    for (size_t j = i; j < 3; ++j) {
      const T symmetric_direction_acceleration =
          0.5 * (direction[i] * spatial_acceleration[j] +
                 spatial_acceleration[i] * direction[j]);
      const T symmetric_direction_strain =
          0.5 * (direction[i] * strain_direction[j] +
                 strain_direction[i] * direction[j]);
      const double symmetric_strain_rate =
          0.5 * (spatial_strain_rate[i][j] + spatial_strain_rate[j][i]);
      const double identity = i == j ? 1. : 0.;
      const T acceleration_ij =
          2. * f * acceleration_profile_prime *
              symmetric_direction_acceleration -
          acceleration_profile * g_perpendicular_prime * radial_acceleration *
              identity -
          acceleration_profile * g_parallel_prime * radial_acceleration *
              direction[i] * direction[j] -
          2. * acceleration_profile * g_parallel / radius *
              (symmetric_direction_acceleration -
               radial_acceleration * direction[i] * direction[j]);
      const T strain_ij =
          2. * h * (radius / mass) * symmetric_direction_strain +
          2. * f * strain_profile_prime * symmetric_direction_strain +
          2. * g_perpendicular * strain_profile / radius *
              (symmetric_strain_rate - symmetric_direction_strain) -
          strain_profile * g_perpendicular_prime * strain_radial * identity -
          strain_profile * g_parallel_prime * strain_radial * direction[i] *
              direction[j] -
          2. * strain_profile * g_parallel / radius *
              (symmetric_direction_strain -
               strain_radial * direction[i] * direction[j]);
      result[i + 1][j + 1] = acceleration_ij + strain_ij;
      result[j + 1][i + 1] = result[i + 1][j + 1];
    }
  }
  return result;
}

}  // namespace

void local_affine_rate_inverse_metric_response(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& local_coords, const double mass,
    const AffineRates& rates) {
  ASSERT(mass > 0., "The black-hole mass must be positive. Got " << mass);
  const size_t number_of_points = get<0>(local_coords).size();
  set_number_of_grid_points(result, number_of_points);
  for (size_t s = 0; s < number_of_points; ++s) {
    const std::array<double, 3> point{{local_coords.get(0)[s],
                                       local_coords.get(1)[s],
                                       local_coords.get(2)[s]}};
    const double radius =
        sqrt(square(point[0]) + square(point[1]) + square(point[2]));
    ASSERT(radius > mass,
           "The affine-rate response is restricted to the exterior harmonic "
           "chart, rho > M. Got rho="
               << radius << " and M=" << mass);
    const auto response = local_response_at_point(point, mass, rates);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        result->get(a, b)[s] = response[a][b];
      }
    }
  }
}

namespace {
using InverseMetricDerivatives = std::array<tnsr::AA<DataVector, 3>, 4>;

void affine_rate_inverse_metric_response_and_derivatives(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> result,
    const gsl::not_null<InverseMetricDerivatives*> derivatives,
    const tnsr::I<DataVector, 3>& x, const double time, const double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix,
    const AffineRates& rates) {
  ASSERT(mass > 0., "The black-hole mass must be positive. Got " << mass);
  const size_t number_of_points = get<0>(x).size();
  const auto inverse_frame = exact_frame::inverse(frame_map_matrix);
  set_number_of_grid_points(result, number_of_points);
  for (auto& derivative : *derivatives) {
    set_number_of_grid_points(make_not_null(&derivative), number_of_points);
  }

  for (size_t s = 0; s < number_of_points; ++s) {
    std::array<FirstDerivative<3>, 3> local_point{};
    for (size_t i = 0; i < 3; ++i) {
      local_point[i].value = inverse_frame[i + 1][0] * time;
      local_point[i].derivative[i] = 1.;
      for (size_t j = 0; j < 3; ++j) {
        local_point[i].value +=
            inverse_frame[i + 1][j + 1] * (x.get(j)[s] - center[j]);
      }
    }
    double radius_squared = 0.;
    for (size_t i = 0; i < 3; ++i) {
      radius_squared += square(local_point[i].value);
    }
    ASSERT(sqrt(radius_squared) > mass,
           "The affine-rate response is restricted to the exterior harmonic "
           "chart, rho > M. Got rho="
               << sqrt(radius_squared) << " and M=" << mass);
    const auto local_response =
        local_response_at_point(local_point, mass, rates);

    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        result->get(a, b)[s] = 0.;
        for (size_t c = 0; c < 4; ++c) {
          derivatives->at(c).get(a, b)[s] = 0.;
        }
        for (size_t mu = 0; mu < 4; ++mu) {
          for (size_t nu = 0; nu < 4; ++nu) {
            const double tensor_transform =
                frame_map_matrix[a][mu] * frame_map_matrix[b][nu];
            result->get(a, b)[s] +=
                tensor_transform * local_response[mu][nu].value;
            for (size_t c = 0; c < 4; ++c) {
              for (size_t i = 0; i < 3; ++i) {
                derivatives->at(c).get(a, b)[s] +=
                    tensor_transform * local_response[mu][nu].derivative[i] *
                    inverse_frame[i + 1][c];
              }
            }
          }
        }
      }
    }
  }
}
}  // namespace

void affine_rate_inverse_metric_response(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& x, const double time, const double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix,
    const AffineRates& rates) {
  InverseMetricDerivatives unused_derivatives{};
  affine_rate_inverse_metric_response_and_derivatives(
      result, make_not_null(&unused_derivatives), x, time, mass, center,
      frame_map_matrix, rates);
}

void inverse_metric(const gsl::not_null<tnsr::AA<DataVector, 3>*> result,
                    const tnsr::I<DataVector, 3>& x, const double time,
                    const double mass, const std::array<double, 3>& center,
                    const exact_frame::FrameMatrix& frame_map_matrix,
                    const AffineRates& rates) {
  exact_frame::inverse_metric(result, x, time, mass, center, frame_map_matrix);
  tnsr::AA<DataVector, 3> rate_response{};
  affine_rate_inverse_metric_response(make_not_null(&rate_response), x, time,
                                      mass, center, frame_map_matrix, rates);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      result->get(a, b) += rate_response.get(a, b);
    }
  }
}

void affine_rate_evolved_variables_response(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric_response,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi_response,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi_response,
    const tnsr::I<DataVector, 3>& x, const double time, const double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix,
    const AffineRates& rates) {
  const size_t number_of_points = get<0>(x).size();
  tnsr::aa<DataVector, 3> background_metric{};
  tnsr::aa<DataVector, 3> background_pi{};
  tnsr::iaa<DataVector, 3> background_phi{};
  exact_frame::evolved_variables(
      make_not_null(&background_metric), make_not_null(&background_pi),
      make_not_null(&background_phi), x, time, mass, center, frame_map_matrix);
  const tnsr::AA<DataVector, 3> background_inverse_metric =
      determinant_and_inverse(background_metric).second;

  tnsr::AA<DataVector, 3> inverse_metric_response{};
  InverseMetricDerivatives inverse_metric_response_derivatives{};
  affine_rate_inverse_metric_response_and_derivatives(
      make_not_null(&inverse_metric_response),
      make_not_null(&inverse_metric_response_derivatives), x, time, mass,
      center, frame_map_matrix, rates);

  set_number_of_grid_points(spacetime_metric_response, number_of_points);
  set_number_of_grid_points(pi_response, number_of_points);
  set_number_of_grid_points(phi_response, number_of_points);
  // delta g = -g_0 delta G g_0.
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      spacetime_metric_response->get(a, b) = 0.;
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          spacetime_metric_response->get(a, b) -=
              background_metric.get(a, c) * inverse_metric_response.get(c, d) *
              background_metric.get(d, b);
        }
      }
    }
  }

  const DataVector background_lapse =
      1. / sqrt(-get<0, 0>(background_inverse_metric));
  std::array<DataVector, 3> background_shift{};
  for (size_t i = 0; i < 3; ++i) {
    background_shift[i] = -background_inverse_metric.get(0, i + 1) /
                          get<0, 0>(background_inverse_metric);
  }
  tnsr::aa<DataVector, 3> dt_background_metric(number_of_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_background_metric.get(a, b) =
          -background_lapse * background_pi.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        dt_background_metric.get(a, b) +=
            background_shift[i] * background_phi.get(i, a, b);
      }
    }
  }

  std::array<tnsr::aa<DataVector, 3>, 4> metric_response_derivatives{};
  for (auto& derivative : metric_response_derivatives) {
    derivative = tnsr::aa<DataVector, 3>(number_of_points, 0.);
  }
  // Differentiate delta g = -g_0 delta G g_0 once. All three terms are
  // linear in the affine rates.
  for (size_t derivative_index = 0; derivative_index < 4; ++derivative_index) {
    const auto background_derivative =
        [&background_phi, &dt_background_metric, derivative_index](
            const size_t a, const size_t b) -> const DataVector& {
      return derivative_index == 0
                 ? dt_background_metric.get(a, b)
                 : background_phi.get(derivative_index - 1, a, b);
    };
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t c = 0; c < 4; ++c) {
          for (size_t d = 0; d < 4; ++d) {
            metric_response_derivatives[derivative_index].get(a, b) -=
                background_derivative(a, c) *
                    inverse_metric_response.get(c, d) *
                    background_metric.get(d, b) +
                background_metric.get(a, c) *
                    inverse_metric_response_derivatives[derivative_index].get(
                        c, d) *
                    background_metric.get(d, b) +
                background_metric.get(a, c) *
                    inverse_metric_response.get(c, d) *
                    background_derivative(d, b);
          }
        }
        if (derivative_index > 0) {
          phi_response->get(derivative_index - 1, a, b) =
              metric_response_derivatives[derivative_index].get(a, b);
        }
      }
    }
  }

  const DataVector lapse_response =
      0.5 * cube(background_lapse) * get<0, 0>(inverse_metric_response);
  std::array<DataVector, 3> shift_response{};
  for (size_t i = 0; i < 3; ++i) {
    shift_response[i] = -inverse_metric_response.get(0, i + 1) /
                            get<0, 0>(background_inverse_metric) +
                        background_inverse_metric.get(0, i + 1) *
                            get<0, 0>(inverse_metric_response) /
                            square(get<0, 0>(background_inverse_metric));
  }
  // delta Pi = delta[(beta^i Phi_i - dt g) / alpha].
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      pi_response->get(a, b) =
          -metric_response_derivatives[0].get(a, b) / background_lapse -
          lapse_response / background_lapse * background_pi.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        pi_response->get(a, b) +=
            (shift_response[i] * background_phi.get(i, a, b) +
             background_shift[i] * phi_response->get(i, a, b)) /
            background_lapse;
      }
    }
  }
}

void evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, const double time, const double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix,
    const AffineRates& rates) {
  exact_frame::evolved_variables(spacetime_metric, pi, phi, x, time, mass,
                                 center, frame_map_matrix);
  tnsr::aa<DataVector, 3> metric_response{};
  tnsr::aa<DataVector, 3> local_pi_response{};
  tnsr::iaa<DataVector, 3> local_phi_response{};
  affine_rate_evolved_variables_response(
      make_not_null(&metric_response), make_not_null(&local_pi_response),
      make_not_null(&local_phi_response), x, time, mass, center,
      frame_map_matrix, rates);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      spacetime_metric->get(a, b) += metric_response.get(a, b);
      pi->get(a, b) += local_pi_response.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        phi->get(i, a, b) += local_phi_response.get(i, a, b);
      }
    }
  }
}

}  // namespace gh::Solutions::order_by_order_worldtube
