// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/Chalmers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/StripeIterator.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace Chalmers {
namespace {
// The tail divided by the squared amplitude of the first omitted mode:
// sum_{j=0}^infinity 2 exp(-2 sigma j)/(2(p+1+j)+1).
double tail_factor(const double sigma, const size_t degree) {
  if (sigma <= 0.1) {
    // Sum_{k=0}^infinity q^k/(2k+1) = atanh(sqrt(q))/sqrt(q).
    // Only use this form for q near 1, where subtraction of the first p+1
    // terms is well conditioned for the supported Legendre degrees. expm1
    // avoids rounding sqrt(q) to 1 for very small positive sigma.
    const long double s = sigma;
    const long double t = std::exp(-s);
    const long double q = t * t;
    long double finite_sum = 0.0L;
    long double power = 1.0L;
    for (size_t k = 0; k <= degree; ++k) {
      finite_sum += power / static_cast<long double>(2 * k + 1);
      power *= q;
    }
    const long double full_sum =
        0.5L * (std::log1p(t) - std::log(-std::expm1(-s))) / t;
    return static_cast<double>(2.0L * std::exp(2.0L * s * (degree + 1)) *
                               (full_sum - finite_sum));
  }
  const double q = std::exp(-2.0 * sigma);
  double sum = 0.0;
  double power = 1.0;
  for (size_t j = 0;; ++j) {
    sum += 2.0 * power / static_cast<double>(2 * (degree + 1 + j) + 1);
    power *= q;
    const double remainder_bound =
        2.0 * power /
        (static_cast<double>(2 * (degree + 2 + j) + 1) * (1.0 - q));
    if (remainder_bound <= 4.0 * std::numeric_limits<double>::epsilon() * sum) {
      return sum;
    }
  }
}
// Infinite tail of one parity, divided by its first omitted amplitude squared.
double parity_tail_factor(const double sigma, const size_t last_degree) {
  const size_t first = last_degree + 2;
  if (sigma <= 0.1) {
    const long double s = sigma;
    const long double t = std::exp(-s);
    const long double atanh_t =
        0.5L * (std::log1p(t) - std::log(-std::expm1(-s)));
    const long double full_sum =
        (atanh_t + (last_degree % 2 == 0 ? 1.0L : -1.0L) * std::atan(t)) /
        (2.0L * t);
    long double finite_sum = 0.0L;
    for (size_t k = last_degree % 2; k <= last_degree; k += 2) {
      finite_sum += std::exp(-2.0L * s * k) / (2.0L * k + 1.0L);
    }
    return static_cast<double>(2.0L * std::exp(2.0L * s * first) *
                               (full_sum - finite_sum));
  }
  const double q = std::exp(-4.0 * sigma);
  double sum = 0.0;
  double power = 1.0;
  for (size_t j = 0;; ++j) {
    sum += 2.0 * power / static_cast<double>(2 * (first + 2 * j) + 1);
    power *= q;
    const double remainder =
        2.0 * power /
        (static_cast<double>(2 * (first + 2 * (j + 1)) + 1) * (1.0 - q));
    if (remainder <= 4.0 * std::numeric_limits<double>::epsilon() * sum) {
      return sum;
    }
  }
}

template <size_t Dim>
DataVector directional_spectrum(const ModalVector& coefficients,
                                const Mesh<Dim>& mesh, const size_t d) {
  DataVector spectrum(mesh.extents(d), 0.0);
  for (StripeIterator it(mesh.extents(), d); it; ++it) {
    double transverse_weight = 1.0;
    size_t index = it.offset();
    for (size_t transverse = 0; transverse < Dim; ++transverse) {
      const size_t mode = index % mesh.extents(transverse);
      index /= mesh.extents(transverse);
      if (transverse != d) {
        transverse_weight *= 2.0 / (2.0 * mode + 1.0);
      }
    }
    const double norm = std::sqrt(transverse_weight);
    for (size_t k = 0; k < spectrum.size(); ++k) {
      spectrum[k] = std::hypot(
          spectrum[k], norm * coefficients[it.offset() + k * it.stride()]);
    }
  }
  return spectrum;
}
}  // namespace

SequenceEstimate estimate_sequence(
    const DataVector& coefficients,
    const std::optional<double> coefficient_floor) {
  if (coefficients.size() < 2) {
    ERROR("Chalmers estimator requires at least two modes.");
  }
  double scale = 0.0;
  for (const double coefficient : coefficients) {
    if (not std::isfinite(coefficient)) {
      ERROR("Chalmers estimator received non-finite coefficients.");
    }
    scale = std::max(scale, std::abs(coefficient));
  }
  const double floor = coefficient_floor.value_or(
      100.0 * std::numeric_limits<double>::epsilon() * scale);
  if (not std::isfinite(floor) or floor < 0.0) {
    ERROR("Chalmers coefficient floor must be finite and nonnegative.");
  }
  SequenceEstimate result{};
  size_t last_significant = 0;
  double mean_k = 0.0;
  double mean_log = 0.0;
  for (size_t k = 0; k < coefficients.size(); ++k) {
    const double amplitude = std::abs(coefficients[k]);
    if (amplitude > floor) {
      ++result.significant_modes;
      last_significant = k;
      mean_k += static_cast<double>(k);
      mean_log += std::log(amplitude) - std::log(scale);
    }
  }
  if (result.significant_modes == 0) {
    return result;
  }
  if (last_significant < coefficients.size() - 2) {
    result.status = FitStatus::TerminatedTail;
    return result;
  }
  result.error = std::numeric_limits<double>::infinity();
  if (result.significant_modes < 3) {
    result.status = FitStatus::InsufficientModes;
    return result;
  }
  mean_k /= static_cast<double>(result.significant_modes);
  mean_log /= static_cast<double>(result.significant_modes);
  double covariance = 0.0;
  double variance = 0.0;
  for (size_t k = 0; k < coefficients.size(); ++k) {
    if (std::abs(coefficients[k]) > floor) {
      const double dk = static_cast<double>(k) - mean_k;
      covariance += dk * (std::log(std::abs(coefficients[k])) -
                          std::log(scale) - mean_log);
      variance += dk * dk;
    }
  }
  result.sigma = -covariance / variance;
  const double intercept = mean_log + result.sigma * mean_k;
  double residual = 0.0;
  for (size_t k = 0; k < coefficients.size(); ++k) {
    if (std::abs(coefficients[k]) > floor) {
      const double difference = std::log(std::abs(coefficients[k])) -
                                std::log(scale) - intercept +
                                result.sigma * static_cast<double>(k);
      residual += difference * difference;
    }
  }
  result.log_fit_rms =
      std::sqrt(residual / static_cast<double>(result.significant_modes));
  if (result.sigma <= 0.0) {
    result.status = FitStatus::Nondecaying;
    return result;
  }
  const size_t degree = coefficients.size() - 1;
  const double first_omitted =
      scale * std::exp(intercept - result.sigma * (degree + 1));
  result.error =
      std::hypot(std::sqrt(2.0 / static_cast<double>(2 * degree + 1)) *
                     std::abs(coefficients[degree]),
                 first_omitted * std::sqrt(tail_factor(result.sigma, degree)));
  result.status = FitStatus::Fit;
  return result;
}

DirectionEstimate estimate_sequence_parity(const DataVector& coefficients,
                                           const double variation_budget) {
  if (coefficients.size() < 2 or not std::isfinite(variation_budget) or
      variation_budget < 0.0) {
    ERROR(
        "Parity fit requires at least two modes and a finite nonnegative "
        "budget.");
  }
  double scale = 0.0;
  for (const double coefficient : coefficients) {
    if (not std::isfinite(coefficient)) {
      ERROR("Chalmers estimator received non-finite coefficients.");
    }
    scale = std::max(scale, std::abs(coefficient));
  }
  const double floor = 100.0 * std::numeric_limits<double>::epsilon() * scale;
  DirectionEstimate result{};
  for (size_t parity = 0; parity < 2; ++parity) {
    const size_t last_degree =
        coefficients.size() - 1 - (coefficients.size() - 1 - parity) % 2;
    double variation = 0.0;
    size_t count = 0;
    size_t last_significant = parity;
    double mean_k = 0.0;
    double mean_log = 0.0;
    for (size_t k = parity; k < coefficients.size(); k += 2) {
      const double amplitude = std::abs(coefficients[k]);
      if (k != 0) {
        variation =
            std::hypot(variation, std::sqrt(2.0 / (2.0 * k + 1.0)) * amplitude);
      }
      if (amplitude > floor) {
        ++count;
        last_significant = k;
        mean_k += static_cast<double>(k);
        mean_log += std::log(amplitude / scale);
      }
    }
    if (variation_budget > 0.0 and
        variation <= 0.01 * variation_budget / std::sqrt(2.0)) {
      ++result.ignored_sequences;
      result.error = std::hypot(result.error, variation);
      continue;
    }
    if (count == 0 or last_significant + 4 <= last_degree) {
      continue;
    }
    if (count < 3) {
      ++result.insufficient_sequences;
      result.error = std::numeric_limits<double>::infinity();
      continue;
    }
    mean_k /= static_cast<double>(count);
    mean_log /= static_cast<double>(count);
    double covariance = 0.0;
    double variance = 0.0;
    for (size_t k = parity; k < coefficients.size(); k += 2) {
      if (std::abs(coefficients[k]) > floor) {
        const double dk = static_cast<double>(k) - mean_k;
        covariance +=
            dk * (std::log(std::abs(coefficients[k]) / scale) - mean_log);
        variance += dk * dk;
      }
    }
    const double sigma = -covariance / variance;
    const double intercept = mean_log + sigma * mean_k;
    double residual = 0.0;
    for (size_t k = parity; k < coefficients.size(); k += 2) {
      if (std::abs(coefficients[k]) > floor) {
        const double difference = std::log(std::abs(coefficients[k]) / scale) -
                                  intercept + sigma * static_cast<double>(k);
        residual += difference * difference;
      }
    }
    result.max_log_fit_rms =
        std::max(result.max_log_fit_rms, std::sqrt(residual / count));
    result.sigma =
        std::isfinite(result.sigma) ? std::min(result.sigma, sigma) : sigma;
    if (sigma <= 0.0) {
      result.error = std::numeric_limits<double>::infinity();
      continue;
    }
    const double first_omitted =
        scale * std::exp(intercept - sigma * (last_degree + 2));
    result.error = std::hypot(
        result.error, std::hypot(std::sqrt(2.0 / (2.0 * last_degree + 1.0)) *
                                     coefficients[last_degree],
                                 first_omitted * std::sqrt(parity_tail_factor(
                                                     sigma, last_degree))));
  }
  if (result.insufficient_sequences != 0) {
    result.sigma = std::numeric_limits<double>::quiet_NaN();
  }
  return result;
}

template <size_t Dim>
std::array<DirectionEstimate, Dim> estimate(
    const DataVector& values, const Mesh<Dim>& mesh,
    const double absolute_stripe_budget, const size_t ignored_highest_modes) {
  if (not std::isfinite(absolute_stripe_budget) or
      absolute_stripe_budget < 0.0) {
    ERROR("Chalmers stripe budget must be finite and nonnegative.");
  }
  for (size_t d = 0; d < Dim; ++d) {
    if (mesh.basis(d) != Spectral::Basis::Legendre or
        mesh.extents(d) < 2 + ignored_highest_modes) {
      ERROR(
          "Chalmers estimator requires tensor-product Legendre meshes with "
          "at least two points per dimension.");
    }
  }
  const auto coefficients = to_modal_coefficients(values, mesh);
  for (const double coefficient : coefficients) {
    if (not std::isfinite(coefficient)) {
      ERROR("Chalmers estimator received non-finite coefficients.");
    }
  }
  const double floor =
      100.0 * std::numeric_limits<double>::epsilon() * max(abs(coefficients));
  std::array<DirectionEstimate, Dim> result{};
  for (size_t d = 0; d < Dim; ++d) {
    auto& direction = gsl::at(result, d);
    DataVector stripe(mesh.extents(d) - ignored_highest_modes);
    const double stripe_floor =
        absolute_stripe_budget /
        std::sqrt(static_cast<double>(
            (size_t{1} << (Dim - 1)) *
            (mesh.number_of_grid_points() / mesh.extents(d))));
    for (StripeIterator it(mesh.extents(), d); it; ++it) {
      for (size_t k = 0; k < stripe.size(); ++k) {
        stripe[k] = coefficients[it.offset() + k * it.stride()];
      }
      if (absolute_stripe_budget > 0.0) {
        double norm = 0.0;
        // The k=0 coefficient is constant in this direction and cannot
        // require directional refinement. Bound only the varying part.
        for (size_t k = 1; k < stripe.size(); ++k) {
          norm = std::hypot(norm, std::sqrt(2.0 / (2.0 * k + 1.0)) * stripe[k]);
        }
        if (norm <= stripe_floor) {
          ++direction.ignored_sequences;
          continue;
        }
      }
      const auto fit = estimate_sequence(stripe, floor);
      direction.error = std::max(direction.error, fit.error);
      if (std::isfinite(fit.sigma)) {
        direction.sigma = std::isfinite(direction.sigma)
                              ? std::min(direction.sigma, fit.sigma)
                              : fit.sigma;
      }
      if (std::isfinite(fit.log_fit_rms)) {
        direction.max_log_fit_rms =
            std::max(direction.max_log_fit_rms, fit.log_fit_rms);
      }
      if (fit.status == FitStatus::InsufficientModes) {
        ++direction.insufficient_sequences;
      }
    }
    if (direction.insufficient_sequences != 0) {
      direction.sigma = std::numeric_limits<double>::quiet_NaN();
    }
  }
  return result;
}

template <size_t Dim>
std::array<DirectionEstimate, Dim> estimate_directional_energy(
    const DataVector& values, const Mesh<Dim>& mesh,
    const double absolute_variation_budget, const size_t ignored_highest_modes,
    const bool separate_parity) {
  if (not std::isfinite(absolute_variation_budget) or
      absolute_variation_budget < 0.0) {
    ERROR("Chalmers variation budget must be finite and nonnegative.");
  }
  for (size_t d = 0; d < Dim; ++d) {
    if (mesh.basis(d) != Spectral::Basis::Legendre or
        mesh.extents(d) < 2 + ignored_highest_modes) {
      ERROR(
          "Chalmers estimator requires tensor-product Legendre meshes with "
          "at least two points per dimension.");
    }
  }
  const auto coefficients = to_modal_coefficients(values, mesh);
  for (const double coefficient : coefficients) {
    if (not std::isfinite(coefficient)) {
      ERROR("Chalmers estimator received non-finite coefficients.");
    }
  }
  std::array<DirectionEstimate, Dim> result{};
  for (size_t d = 0; d < Dim; ++d) {
    const auto full_spectrum = directional_spectrum(coefficients, mesh, d);
    DataVector spectrum(full_spectrum.size() - ignored_highest_modes);
    std::copy_n(full_spectrum.begin(), spectrum.size(), spectrum.begin());
    auto& direction = gsl::at(result, d);
    double variation = 0.0;
    for (size_t k = 1; k < spectrum.size(); ++k) {
      variation =
          std::hypot(variation, std::sqrt(2.0 / (2.0 * k + 1.0)) * spectrum[k]);
    }
    if (absolute_variation_budget > 0.0 and
        variation <= absolute_variation_budget) {
      direction.ignored_sequences = 1;
      continue;
    }
    if (separate_parity) {
      direction = estimate_sequence_parity(spectrum, absolute_variation_budget);
      continue;
    }
    const auto fit = estimate_sequence(spectrum);
    direction.error = fit.error;
    direction.sigma = fit.sigma;
    direction.max_log_fit_rms =
        std::isfinite(fit.log_fit_rms) ? fit.log_fit_rms : 0.0;
    direction.insufficient_sequences =
        fit.status == FitStatus::InsufficientModes ? 1 : 0;
  }
  return result;
}

template <size_t Dim>
std::array<double, Dim> p_coarsening_error(
    const DataVector& values, const Mesh<Dim>& mesh,
    const double absolute_variation_budget, const size_t ignored_highest_modes,
    const bool separate_parity) {
  if (not std::isfinite(absolute_variation_budget) or
      absolute_variation_budget < 0.0) {
    ERROR("Chalmers variation budget must be finite and nonnegative.");
  }
  for (size_t d = 0; d < Dim; ++d) {
    if (mesh.basis(d) != Spectral::Basis::Legendre or
        mesh.extents(d) < 2 + ignored_highest_modes) {
      ERROR(
          "Chalmers coarsening requires Legendre meshes with at least two "
          "modes.");
    }
  }
  const auto coefficients = to_modal_coefficients(values, mesh);
  for (const double coefficient : coefficients) {
    if (not std::isfinite(coefficient)) {
      ERROR("Chalmers estimator received non-finite coefficients.");
    }
  }
  std::array<double, Dim> result{};
  for (size_t d = 0; d < Dim; ++d) {
    auto& error = gsl::at(result, d);
    error = std::numeric_limits<double>::infinity();
    if (mesh.extents(d) <= 6 + ignored_highest_modes) {
      continue;
    }
    const auto spectrum = directional_spectrum(coefficients, mesh, d);
    double variation = 0.0;
    for (size_t k = 1; k < spectrum.size(); ++k) {
      variation =
          std::hypot(variation, std::sqrt(2.0 / (2.0 * k + 1.0)) * spectrum[k]);
    }
    if (variation <= absolute_variation_budget) {
      error = variation;
      continue;
    }
    const size_t degree = spectrum.size() - 1 - ignored_highest_modes;
    DataVector retained(degree);
    std::copy_n(spectrum.begin(), degree, retained.begin());
    double discarded = 0.0;
    for (size_t k = degree; k < spectrum.size(); ++k) {
      discarded =
          std::hypot(discarded, std::sqrt(2.0 / (2.0 * k + 1.0)) * spectrum[k]);
    }
    error = discarded +
            (separate_parity
                 ? estimate_sequence_parity(retained, absolute_variation_budget)
                       .error
                 : estimate_sequence(retained).error);
  }
  return result;
}

template <size_t Dim>
JoinEstimate h_coarsening_error(
    const DataVector& lower_values, const Mesh<Dim>& lower_mesh,
    const DataVector& upper_values, const Mesh<Dim>& upper_mesh,
    const size_t dimension, const double absolute_variation_budget,
    const size_t ignored_highest_modes, const bool separate_parity) {
  if (dimension >= Dim or lower_mesh.basis() != upper_mesh.basis() or
      lower_mesh.quadrature() != upper_mesh.quadrature()) {
    ERROR("Chalmers join requires compatible sibling meshes and a valid axis.");
  }
  auto extents = lower_mesh.extents().indices();
  for (size_t d = 0; d < Dim; ++d) {
    if (lower_mesh.basis(d) != Spectral::Basis::Legendre) {
      ERROR("Chalmers join requires Legendre meshes.");
    }
    extents[d] = std::max(extents[d], upper_mesh.extents(d));
  }
  const Mesh<Dim> parent_mesh{extents, lower_mesh.basis(),
                              lower_mesh.quadrature()};
  const auto full = make_array<Dim>(Spectral::SegmentSize::Full);
  auto lower = full;
  auto upper = full;
  lower[dimension] = Spectral::SegmentSize::LowerHalf;
  upper[dimension] = Spectral::SegmentSize::UpperHalf;
  const DataVector parent =
      Spectral::project(lower_values, lower_mesh, parent_mesh, lower, full) +
      Spectral::project(upper_values, upper_mesh, parent_mesh, upper, full);
  double loss_squared = 0.0;
  for (const bool is_lower : {true, false}) {
    const DataVector difference =
        Spectral::project(parent, parent_mesh, parent_mesh, full,
                          is_lower ? lower : upper) -
        Spectral::project(is_lower ? lower_values : upper_values,
                          is_lower ? lower_mesh : upper_mesh, parent_mesh, full,
                          full);
    const auto modes = to_modal_coefficients(difference, parent_mesh);
    for (size_t i = 0; i < modes.size(); ++i) {
      // Each child occupies one half of the reference parent.
      double weight = 0.5;
      size_t index = i;
      for (size_t d = 0; d < Dim; ++d) {
        weight *= 2.0 / (2.0 * (index % extents[d]) + 1.0);
        index /= extents[d];
      }
      loss_squared += weight * modes[i] * modes[i];
    }
  }
  const auto estimates = estimate_directional_energy(
      parent, parent_mesh, absolute_variation_budget / Dim,
      ignored_highest_modes, separate_parity);
  JoinEstimate result{std::sqrt(loss_squared), std::sqrt(loss_squared),
                      max(abs(parent))};
  // Joining is followed by filtering on the parent. Budget the union of
  // discarded parent modal slabs as well as native projection loss.
  if (ignored_highest_modes > 0) {
    const auto modes = to_modal_coefficients(parent, parent_mesh);
    double filter_loss_squared = 0.0;
    for (size_t i = 0; i < modes.size(); ++i) {
      size_t index = i;
      double weight = 1.0;
      bool discarded = false;
      for (size_t d = 0; d < Dim; ++d) {
        const size_t k = index % extents[d];
        index /= extents[d];
        discarded = discarded or k + ignored_highest_modes >= extents[d];
        weight *= 2.0 / (2.0 * k + 1.0);
      }
      if (discarded) {
        filter_loss_squared += weight * modes[i] * modes[i];
      }
    }
    result.error += std::sqrt(filter_loss_squared);
  }
  for (const auto& estimate : estimates) {
    result.error += estimate.error;
  }
  return result;
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define INSTANTIATE(_, data)                                  \
  template std::array<DirectionEstimate, DIM(data)> estimate( \
      const DataVector& values, const Mesh<DIM(data)>& mesh, double, size_t);
GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))
#undef INSTANTIATE
#define INSTANTIATE(_, data)                                               \
  template std::array<DirectionEstimate, DIM(data)>                        \
  estimate_directional_energy(const DataVector& values,                    \
                              const Mesh<DIM(data)>& mesh, double, size_t, \
                              bool);
GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))
#undef INSTANTIATE
#define INSTANTIATE(_, data)                                                 \
  template std::array<double, DIM(data)> p_coarsening_error(                 \
      const DataVector& values, const Mesh<DIM(data)>& mesh, double, size_t, \
      bool);
GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))
#undef INSTANTIATE
#define INSTANTIATE(_, data)                                        \
  template JoinEstimate h_coarsening_error(                         \
      const DataVector&, const Mesh<DIM(data)>&, const DataVector&, \
      const Mesh<DIM(data)>&, size_t, double, size_t, bool);
GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))
#undef INSTANTIATE
#undef DIM
}  // namespace Chalmers
