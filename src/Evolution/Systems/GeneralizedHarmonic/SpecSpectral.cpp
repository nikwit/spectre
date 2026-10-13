// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"

#include "Utilities/Overloader.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gh {
PUP::able::PUP_ID SpecSpectral::my_PUP_ID = 0;  // NOLINT

SpecSpectral::SpecSpectral(const double spectral_target,
                           const std::array<double, 2>& constraint_targets,
                           const bool force, const size_t minimum_modes,
                           const size_t cooldown, const double adjustment_limit,
                           const double start_time, const bool verbose,
                           const Options::Context& context)
    : SpecSpectral(spectral_target, constraint_targets, force,
                   MinimumModes::type{minimum_modes}, cooldown,
                   adjustment_limit, start_time, verbose, context) {}

SpecSpectral::SpecSpectral(const double spectral_target,
                           const std::array<double, 2>& constraint_targets,
                           const bool force, MinimumModes::type minimum_modes,
                           const size_t cooldown, const double adjustment_limit,
                           const double start_time, const bool verbose,
                           const Options::Context& context)
    : settings_{spectral_target, constraint_targets, force, 5,
                cooldown,        adjustment_limit},
      minimum_modes_(std::move(minimum_modes)),
      start_time_(start_time),
      verbose_(verbose) {
  if (not std::isfinite(spectral_target) or spectral_target <= 0.0 or
      not std::isfinite(start_time) or cooldown < 1) {
    PARSE_ERROR(context,
                "Invalid spectral target, start time, minimum modes "
                "or cooldown for SpecSpectral.");
  }
  const auto validate = [&context](const auto& values) {
    for (const size_t value : values) {
      if (value < 5) {
        PARSE_ERROR(context, "MinimumModes must be at least 5.");
      }
    }
  };
  std::visit(
      Overloader{
          [&validate](const size_t value) { validate(make_array<3>(value)); },
          [&validate](const std::array<size_t, 3>& values) {
            validate(values);
          },
          [&validate,
           &context](const std::vector<std::array<size_t, 3>>& values) {
            if (values.empty()) {
              PARSE_ERROR(context,
                          "MinimumModes block list must not be empty.");
            }
            for (const auto& block : values) {
              validate(block);
            }
          }},
      minimum_modes_);
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

std::array<size_t, 3> SpecSpectral::minimum_modes(const size_t block_id) const {
  return std::visit(
      Overloader{[](const size_t value) { return make_array<3>(value); },
                 [](const std::array<size_t, 3>& values) { return values; },
                 [block_id](const std::vector<std::array<size_t, 3>>& values) {
                   if (block_id >= values.size()) {
                     ERROR("MinimumModes has no entry for block " << block_id);
                   }
                   return values[block_id];
                 }},
      minimum_modes_);
}

void SpecSpectral::pup(PUP::er& p) {
  amr::Criterion::pup(p);
  p | settings_.spectral_target;
  p | settings_.constraint_targets;
  p | settings_.force_pile_up;
  p | settings_.minimum_modes;
  p | settings_.cooldown;
  p | settings_.adjustment_limit;
  p | minimum_modes_;
  p | start_time_;
  p | verbose_;
}
}  // namespace gh
