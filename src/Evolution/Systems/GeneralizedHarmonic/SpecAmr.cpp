// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "NumericalAlgorithms/Interpolation/LinearRegression.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Constraints.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Math.hpp"

namespace gh::spec_amr {
namespace {
constexpr double roundoff = 100.0 * std::numeric_limits<double>::epsilon();
double tail(const DataVector& power, const size_t count, const double floor) {
  if (count < 2) {
    return 1.0;
  }
  if (power[count - 1] <= floor and power[count - 2] <= floor) {
    return 0.0;
  }
  // A missing leading pair cannot demonstrate a small relative error.
  const double leading = std::max(power[0], power[1]);
  if (leading <= floor) {
    return 1.0;
  }
  double weight_sum = 0.0;
  double log_sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    if (power[k] <= floor) {
      continue;
    }
    const double weight =
        exp(-square(static_cast<double>(count - 1 - k) - 0.5));
    weight_sum += weight;
    log_sum += weight * log10(power[k] / leading);
  }
  return weight_sum > 0.0 ? pow(10.0, log_sum / weight_sum) : 1.0;
}
double predict(const double current, const double old) {
  return old <= 1.e-4 * current ? 2.0 * current : current * (current / old);
}
}  // namespace

Spectrum diagnose_power(const DataVector& power, const size_t filtered_modes) {
  ASSERT(power.size() >= filtered_modes + 3, "Insufficient retained modes");
  Spectrum result{};
  const size_t retained = power.size() - filtered_modes;
  double scale = 0.0;
  for (size_t k = 0; k < retained; ++k) {
    ASSERT(std::isfinite(power[k]) and power[k] >= 0.0, "Invalid spectrum");
    scale = std::max(scale, power[k]);
  }
  if (scale == 0.0) {
    return result;
  }
  const double floor = roundoff * scale;
  result.error = tail(power, retained, floor);
  result.next_error = std::max(result.error, tail(power, retained - 1, floor));
  std::vector<double> x{}, y{};
  for (size_t k = 0; k < retained; ++k) {
    if (power[k] > floor) {
      x.push_back(static_cast<double>(k));
      y.push_back(log10(power[k] / scale));
    }
  }
  if (x.size() < 2) {
    result.degenerate = result.next_error == 0.0;
    return result;
  }
  result.degenerate = false;
  if (x.size() < 4) {
    return result;  // insufficient evidence for a pile-up override
  }
  const auto fit = [&x, &y](const size_t begin, const size_t end) {
    if (end - begin == 2) {
      const double slope =
          (y[begin + 1] - y[begin]) / (x[begin + 1] - x[begin]);
      return intrp::LinearRegressionResult{y[begin] - slope * x[begin], slope,
                                           0.0, 0.0};
    }
    const auto first = static_cast<std::ptrdiff_t>(begin);
    const auto last = static_cast<std::ptrdiff_t>(end);
    return intrp::linear_regression(
        std::vector<double>(x.begin() + first, x.begin() + last),
        std::vector<double>(y.begin() + first, y.begin() + last));
  };
  std::vector<intrp::LinearRegressionResult> fits{
      fit(0, std::min(x.size(), std::max(x.size() / 3, size_t{4})))};
  for (size_t start = 0; start < 3; ++start) {
    for (size_t end = start + 5; end <= x.size(); ++end) {
      fits.push_back(fit(start, end));
    }
  }
  double max_error = 0.0;
  for (const auto& f : fits) {
    max_error = std::max(max_error, std::abs(f.delta_slope));
  }
  const double regularizer = std::max(1.e-15, 1.e-3 * max_error);
  double weight_sum = 0.0;
  for (const auto& f : fits) {
    const double weight = 1.0 / (regularizer + std::abs(f.delta_slope));
    result.convergence -= weight * f.slope;
    weight_sum += weight;
  }
  result.convergence /= weight_sum;
  // A flat or growing entire spectrum is unresolved, not excess resolution.
  // Do not count roundoff tails of a resolved polynomial as physical pile-up.
  if (result.convergence <= 1.e-10 or result.error == 0.0) {
    return result;
  }
  for (size_t start = 1; start + 1 < x.size(); ++start) {
    const double local =
        std::max(0.0, -fit(start, std::min(x.size(), start + 5)).slope);
    const double ratio = square(local / result.convergence);
    result.pile_up += ratio < 16.0 ? exp(-32.0 * ratio) : 0.0;
  }
  return result;
}

Diagnostics diagnose(const tnsr::aa<DataVector, 3>& metric,
                     const tnsr::aa<DataVector, 3>& pi,
                     const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
                     const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                                           Frame::Inertial>& inverse_jacobian,
                     const tnsr::I<DataVector, 3>& coordinates,
                     const double time, const gauges::GaugeCondition& gauge,
                     const Scalar<DataVector>& gamma2) {
  Diagnostics result{};
  std::array<std::array<DataVector, 3>, 2> power{};
  for (auto& group : power) {
    for (size_t d = 0; d < 3; ++d) {
      group[d] = DataVector(mesh.extents(d), 0.0);
    }
  }
  const auto accumulate = [&mesh, &power](const auto& tensor,
                                          const size_t group) {
    for (const auto& component : tensor) {
      const auto modes = PowerMonitors::power_monitors(component, mesh);
      for (size_t d = 0; d < 3; ++d) {
        power[group][d] += square(modes[d]);
      }
    }
  };
  accumulate(metric, 0);
  accumulate(pi, 1);
  accumulate(phi, 1);
  for (size_t group = 0; group < 2; ++group) {
    for (size_t d = 0; d < 3; ++d) {
      result.spectra[group][d] = diagnose_power(
          DataVector(sqrt(power[group][d] / (group == 0 ? 10.0 : 40.0))));
    }
  }

  const auto spatial = gr::spatial_metric(metric);
  const auto [det_spatial, inverse_spatial] = determinant_and_inverse(spatial);
  const auto shift = gr::shift(metric, inverse_spatial);
  const auto lapse = gr::lapse(shift, metric);
  const auto normal = gr::spacetime_normal_vector(lapse, shift);
  const auto normal_form =
      gr::spacetime_normal_one_form<DataVector, 3, Frame::Inertial>(lapse);
  const auto inverse_metric =
      gr::inverse_spacetime_metric(lapse, shift, inverse_spatial);
  using GaugeCompute =
      gauges::Tags::GaugeAndDerivativeCompute<3,
                                              gh::Solutions::all_solutions<3>>;
  GaugeCompute::return_type gauge_data{};
  GaugeCompute::function(make_not_null(&gauge_data), lapse, shift, normal,
                         Scalar<DataVector>{sqrt(get(det_spatial))},
                         inverse_spatial, metric, pi, phi, mesh, time,
                         coordinates, inverse_jacobian, gauge);
  const auto c3 = gh::three_index_constraint(
      partial_derivative(metric, mesh, inverse_jacobian), phi);
  const auto c2 = gh::two_index_constraint(
      get<gh::Tags::SpacetimeDerivGaugeH<DataVector, 3>>(gauge_data),
      normal_form, normal, inverse_spatial, inverse_metric, pi, phi,
      partial_derivative(pi, mesh, inverse_jacobian),
      partial_derivative(phi, mesh, inverse_jacobian), gamma2, c3);
  const auto jacobian = determinant_and_inverse(inverse_jacobian).second;
  tnsr::i<DataVector, 3, Frame::ElementLogical> normalization{}, projected{};
  DataVector buffer{};
  amr::Criteria::Constraints_detail::normalization_factor_square(
      make_not_null(&normalization), jacobian);
  const auto project = [&](const auto& constraint, const size_t group) {
    amr::Criteria::Constraints_detail::logical_constraints(
        make_not_null(&projected), make_not_null(&buffer), constraint, jacobian,
        normalization);
    for (size_t d = 0; d < 3; ++d) {
      result.constraints[group][d] =
          sqrt(sum(square(projected.get(d))) /
               static_cast<double>(mesh.number_of_grid_points()));
    }
  };
  project(c2, 0);
  project(c3, 1);
  return result;
}

std::array<double, 3> updated_target_shifts(const History& history,
                                            const Diagnostics& diagnostics,
                                            const Settings& settings) {
  ASSERT(std::isfinite(settings.adjustment_limit) and
             settings.adjustment_limit <= 0.0,
         "The target adjustment limit must be finite and nonpositive");
  std::array<double, 3> result{};
  if (settings.adjustment_limit == 0.0) {
    return result;
  }
  for (size_t d = 0; d < 3; ++d) {
    const double shift = history.log10_target_shift[d];
    result[d] = shift;
    double convergence = std::numeric_limits<double>::infinity();
    double error = 0.0;
    for (const auto& group : diagnostics.spectra) {
      if (not group[d].degenerate) {
        convergence = std::min(convergence, group[d].convergence);
        error = std::max(error, group[d].error);
      }
    }
    // A negative slope cannot define a positive number of spectral modes to
    // tighten by. Retain the previous target until a usable fit is available.
    if (not std::isfinite(convergence) or convergence <= 1.e-10) {
      continue;
    }
    const double difference =
        error > 0.0 ? log10(error) - log10(settings.spectral_target)
                    : -std::numeric_limits<double>::infinity();
    double most_restrictive = std::numeric_limits<double>::infinity();
    for (size_t c = 0; c < diagnostics.constraints.size(); ++c) {
      const double current = diagnostics.constraints[c][d];
      const double previous = history.constraints[c][d];
      const bool have_previous = history.calls > 0;
      // Unlike the prospective p-removal forecast, SpEC's target controller
      // retains the previous constraint even when a p change was realized.
      const double predicted =
          have_previous ? predict(current, previous) : current;
      double next_shift = shift;
      if (predicted > settings.constraint_targets[c]) {
        if (shift == 0.0 and difference > 0.0) {
          next_shift = std::min(0.0, difference - 0.8 * convergence);
        } else if (have_previous and current >= previous and
                   shift > settings.adjustment_limit) {
          next_shift = std::min(shift - 0.5 * convergence,
                                difference - 0.8 * convergence);
        } else if (shift > settings.adjustment_limit) {
          next_shift = shift - 0.5 * convergence;
        }
      } else if (have_previous and current < previous and
                 std::max(current, predicted) <
                     0.5 * settings.constraint_targets[c] and
                 shift < 0.0) {
        next_shift = shift + 0.1 * convergence;
      }
      most_restrictive = std::min(most_restrictive, next_shift);
    }
    result[d] = std::clamp(most_restrictive, settings.adjustment_limit, 0.0);
  }
  return result;
}

namespace {
std::array<amr::Flag, 3> decide_impl(const gsl::not_null<History*> history,
                                     const Diagnostics& diagnostics,
                                     const std::array<size_t, 3>& extents,
                                     const double time,
                                     std::array<amr::Flag, 3> flags,
                                     const Settings& settings,
                                     const bool standalone) {
  // Repeated observations/evaluations at a fixed time must be idempotent.
  if (history->calls > 0 and time == history->time and
      extents == history->extents) {
    return history->decision;
  }
  ASSERT(time >= history->time, "AMR history ran backwards");
  const size_t call = history->calls;
  const bool same_mesh = call > 0 and extents == history->extents;
  const bool joining =
      std::any_of(flags.begin(), flags.end(),
                  [](const auto f) { return f == amr::Flag::Join; });
  bool join_allowed = true;
  for (size_t group = 0; group < 2; ++group) {
    for (size_t d = 0; d < 3; ++d) {
      const double c = diagnostics.constraints[group][d];
      if (not std::isfinite(c)) {
        ERROR("Nonfinite projected GH constraint");
      }
      join_allowed &= 1.25 * c <= settings.constraint_targets[group];
    }
  }
  history->log10_target_shift =
      updated_target_shifts(*history, diagnostics, settings);
  for (size_t d = 0; d < 3; ++d) {
    const double target =
        settings.spectral_target * pow(10.0, history->log10_target_shift[d]);
    if (call > 0 and extents[d] != history->extents[d]) {
      // Change took place following the previous evaluation, which is when
      // SpEC starts its cooldown. Use realized extents, not requested flags.
      history->last_change[d] = call - 1;
    }
    if (flags[d] == amr::Flag::Split) {
      continue;
    }
    if (flags[d] == amr::Flag::Join) {
      if (not join_allowed) {
        flags[d] = amr::Flag::DoNothing;
      }
      continue;
    }
    if (joining) {
      // The sibling criterion checked a parent at the current child extents.
      // Simultaneous p-down would invalidate that prospective-parent check.
      continue;
    }
    if (call > 0 and call < history->last_change[d] + settings.cooldown) {
      flags[d] = amr::Flag::DoNothing;
      continue;
    }
    double error = 0.0, next = 0.0;
    double pile_up = std::numeric_limits<double>::infinity();
    bool has_spectrum = false;
    bool constraints_ok = true;
    bool constraint_excess = false;
    // The derivative spectrum predicts the two-index constraint and the
    // metric spectrum predicts the three-index constraint, as in SpEC.
    for (size_t group = 0; group < 2; ++group) {
      const auto& spectrum = diagnostics.spectra[group][d];
      error = std::max(error, spectrum.error);
      next = std::max(next, spectrum.next_error);
      if (not spectrum.degenerate) {
        has_spectrum = true;
        pile_up = std::min(pile_up, spectrum.pile_up);
      }
      const size_t constraint = 1 - group;
      const double current = diagnostics.constraints[constraint][d];
      const double predicted =
          same_mesh ? predict(current, history->constraints[constraint][d])
                    : current;
      const double addition_forecast =
          (standalone or settings.adjustment_limit < 0.0) and call > 0
              ? predict(current, history->constraints[constraint][d])
              : predicted;
      constraint_excess |=
          addition_forecast >= settings.constraint_targets[constraint];
      const double ratio =
          spectrum.error > 0.0
              ? std::max(1.0, spectrum.next_error / spectrum.error)
              : (spectrum.next_error > 0.0
                     ? std::numeric_limits<double>::infinity()
                     : 1.0);
      const double predicted_next = predicted == 0.0 ? 0.0 : predicted * ratio;
      constraints_ok &=
          1.25 * predicted_next <= settings.constraint_targets[constraint];
    }
    if (not has_spectrum) {
      pile_up = 0.0;
    }
    const bool ordinary =
        pow(10.0, 0.1) * next <= settings.spectral_target and constraints_ok;
    const bool forced = settings.force_pile_up and pile_up > 3.0 and
                        (pile_up - 1.1 > 3.0 or next < target);
    const bool excessive = pile_up > 3.0;
    if (standalone or flags[d] == amr::Flag::IncreaseResolution) {
      // Forced pile-up has priority over p-addition. Otherwise apply SpEC's
      // pile-up hysteresis before adding points.
      if (not forced) {
        // Retain SpEC's base-target/constraint permission gate, and also
        // require LargeTrunc against the shifted target. Only the legacy
        // unadjusted hybrid relies on Chalmers to supply that accuracy check.
        if (pile_up <= (settings.force_pile_up ? 1.5 : 3.0) and
            (error >= settings.spectral_target or constraint_excess) and
            ((not standalone and settings.adjustment_limit == 0.0) or
             error > target)) {
          flags[d] = amr::Flag::IncreaseResolution;
          continue;
        }
        flags[d] = amr::Flag::DoNothing;
      }
    }
    const bool small = error < target and next < target;
    const bool pile_drop =
        excessive and ordinary and (pile_up - 1.1 > 3.0 or next < target);
    flags[d] = call > 0 and extents[d] > settings.minimum_modes and
                       (forced or pile_drop or (ordinary and small))
                   ? amr::Flag::DecreaseResolution
                   : amr::Flag::DoNothing;
  }
  ++history->calls;
  history->time = time;
  history->extents = extents;
  history->constraints = diagnostics.constraints;
  history->decision = flags;
  return flags;
}
}  // namespace

std::array<amr::Flag, 3> decide(const gsl::not_null<History*> history,
                                const Diagnostics& diagnostics,
                                const std::array<size_t, 3>& extents,
                                const double time,
                                const std::array<amr::Flag, 3> field_decision,
                                const Settings& settings) {
  return decide_impl(history, diagnostics, extents, time, field_decision,
                     settings, false);
}

std::array<amr::Flag, 3> decide_p(const gsl::not_null<History*> history,
                                  const Diagnostics& diagnostics,
                                  const std::array<size_t, 3>& extents,
                                  const double time, const Settings& settings) {
  return decide_impl(
      history, diagnostics, extents, time,
      {{amr::Flag::DoNothing, amr::Flag::DoNothing, amr::Flag::DoNothing}},
      settings, true);
}
}  // namespace gh::spec_amr
