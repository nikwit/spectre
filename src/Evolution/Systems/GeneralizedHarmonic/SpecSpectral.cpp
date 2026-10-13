// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"

namespace gh {
PUP::able::PUP_ID SpecSpectral::my_PUP_ID = 0;  // NOLINT

SpecSpectral::SpecSpectral(const double spectral_target,
                           const std::array<double, 2>& constraint_targets,
                           const bool force, const size_t minimum_modes,
                           const size_t cooldown, const double adjustment_limit,
                           const double start_time, const bool verbose,
                           const Options::Context& context)
    : settings_{spectral_target, constraint_targets, force,
                minimum_modes,   cooldown,           adjustment_limit},
      start_time_(start_time),
      verbose_(verbose) {
  if (not std::isfinite(spectral_target) or spectral_target <= 0.0 or
      not std::isfinite(start_time) or minimum_modes < 5 or cooldown < 1) {
    PARSE_ERROR(context,
                "Invalid spectral target, start time, minimum modes "
                "or cooldown for SpecSpectral.");
  }
  for (const double target : constraint_targets) {
    if (not std::isfinite(target) or target <= 0.0) {
      PARSE_ERROR(context, "Constraint targets must be finite and positive.");
    }
  }
  const double minimum_log = log10(std::numeric_limits<double>::min());
  if (not std::isfinite(adjustment_limit) or adjustment_limit > 0.0 or
      adjustment_limit < minimum_log or
      log10(spectral_target) + adjustment_limit < minimum_log) {
    PARSE_ERROR(context,
                "TruncationErrorAdjustmentLimit must be finite and "
                "nonpositive and leave a representable positive target.");
  }
}

void SpecSpectral::pup(PUP::er& p) {
  amr::Criterion::pup(p);
  p | settings_.spectral_target;
  p | settings_.constraint_targets;
  p | settings_.force_pile_up;
  p | settings_.minimum_modes;
  p | settings_.cooldown;
  p | settings_.adjustment_limit;
  p | start_time_;
  p | verbose_;
}
}  // namespace gh
