// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"

#include <algorithm>

#include "Domain/Structure/BlockGroups.hpp"

#include "Utilities/Math.hpp"
#include "Utilities/Overloader.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gh {
PUP::able::PUP_ID SpecSpectral::my_PUP_ID = 0;  // NOLINT

SpecSpectral::GaussianTarget::GaussianTarget(std::array<double, 3> center_in,
                                             const double width_in,
                                             const double amplitude_in,
                                             const double target_in,
                                             const Options::Context& context)
    : center(center_in),
      width(width_in),
      amplitude(amplitude_in),
      target(target_in) {
  if (not std::isfinite(width) or width <= 0.0 or
      not std::isfinite(amplitude) or amplitude < 0.0 or
      not std::isfinite(target) or target <= 0.0) {
    PARSE_ERROR(context,
                "Gaussian width and target must be finite and positive; "
                "amplitude must be finite and nonnegative.");
  }
  for (const double x : center) {
    if (not std::isfinite(x)) {
      PARSE_ERROR(context, "Gaussian centers must be finite grid coordinates.");
    }
  }
}

void SpecSpectral::GaussianTarget::pup(PUP::er& p) {
  p | center;
  p | width;
  p | amplitude;
  p | target;
}

SpecSpectral::SpatialTarget::SpatialTarget() = default;

SpecSpectral::SpatialTarget::SpatialTarget(
    const double background_in, std::vector<GaussianTarget> gaussians_in,
    const Options::Context& context)
    : background(background_in), gaussians(std::move(gaussians_in)) {
  if (not std::isfinite(background) or background <= 0.0) {
    PARSE_ERROR(context, "Background target must be finite and positive.");
  }
}

double SpecSpectral::SpatialTarget::operator()(
    const tnsr::I<DataVector, 3, Frame::Grid>& coordinates) const {
  double log_target = log10(background);
  double weight_sum = 1.0;
  for (const auto& gaussian : gaussians) {
    DataVector r_squared(coordinates.get(0).size(), 0.0);
    for (size_t d = 0; d < 3; ++d) {
      r_squared +=
          square((coordinates.get(d) - gaussian.center[d]) / gaussian.width);
    }
    // SpEC takes the maximum profile value on the element, not the value
    // at its center or an average of pointwise blended targets.
    const double weight = gaussian.amplitude * exp(-min(r_squared));
    const double new_sum = weight_sum + weight;
    log_target += (weight / new_sum) * (log10(gaussian.target) - log_target);
    weight_sum = new_sum;
  }
  return pow(10.0, log_target);
}

void SpecSpectral::SpatialTarget::pup(PUP::er& p) {
  p | background;
  p | gaussians;
}

SpecSpectral::RadiusDependentStart::RadiusDependentStart() = default;

SpecSpectral::RadiusDependentStart::RadiusDependentStart(
    const double time_in, std::vector<std::string> blocks_in,
    const double radius_plus_time_in, const Options::Context& context)
    : time(time_in),
      blocks(std::move(blocks_in)),
      radius_plus_time(radius_plus_time_in) {
  if (not std::isfinite(time) or not std::isfinite(radius_plus_time) or
      blocks.empty()) {
    PARSE_ERROR(context,
                "Junk timing requires finite times and nonempty Blocks.");
  }
}

double SpecSpectral::RadiusDependentStart::release_time(
    const Domain<3>& domain, const size_t block_id,
    const tnsr::I<DataVector, 3, Frame::Grid>& coordinates) const {
  std::vector<std::string> names{};
  for (const auto& block : domain.blocks()) {
    names.push_back(block.name());
  }
  const auto ids =
      domain::block_ids_from_names(blocks, names, domain.block_groups());
  if (std::find(ids.begin(), ids.end(), block_id) == ids.end()) {
    return time;
  }
  DataVector radius_squared(coordinates[0].size(), 0.0);
  for (size_t d = 0; d < 3; ++d) {
    radius_squared += square(coordinates[d]);
  }
  return std::max(time, sqrt(max(radius_squared)) + radius_plus_time);
}

void SpecSpectral::RadiusDependentStart::pup(PUP::er& p) {
  p | time;
  p | blocks;
  p | radius_plus_time;
}

SpecSpectral::SpecSpectral(RelativeTruncationTarget::type spectral_target,
                           const std::array<double, 2>& constraint_targets,
                           const bool force, const size_t minimum_modes,
                           const size_t cooldown, const double adjustment_limit,
                           StartTime::type start_time, const bool verbose,
                           const double time_between_changes,
                           const Options::Context& context)
    : SpecSpectral(std::move(spectral_target), constraint_targets, force,
                   MinimumModes::type{minimum_modes}, cooldown,
                   adjustment_limit, std::move(start_time), verbose,
                   time_between_changes, context) {}

SpecSpectral::SpecSpectral(RelativeTruncationTarget::type spectral_target,
                           const std::array<double, 2>& constraint_targets,
                           const bool force, MinimumModes::type minimum_modes,
                           const size_t cooldown, const double adjustment_limit,
                           StartTime::type start_time, const bool verbose,
                           const double time_between_changes,
                           const Options::Context& context)
    : settings_{0.0,
                constraint_targets,
                force,
                5,
                cooldown,
                adjustment_limit,
                time_between_changes},
      minimum_modes_(std::move(minimum_modes)),
      verbose_(verbose) {
  if (std::holds_alternative<double>(start_time)) {
    start_time_ = std::get<double>(start_time);
  } else {
    radius_dependent_start_ =
        std::get<RadiusDependentStart>(std::move(start_time));
    start_time_ = radius_dependent_start_->time;
  }
  if (std::holds_alternative<double>(spectral_target)) {
    settings_.spectral_target = std::get<double>(spectral_target);
  } else {
    spatial_target_ = std::get<SpatialTarget>(std::move(spectral_target));
    settings_.spectral_target = spatial_target_->background;
  }
  if (not std::isfinite(settings_.spectral_target) or
      settings_.spectral_target <= 0.0 or not std::isfinite(start_time_) or
      cooldown < 1 or not std::isfinite(time_between_changes) or
      time_between_changes < 0.0) {
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
      log10(settings_.spectral_target) + adjustment_limit < minimum_log) {
    PARSE_ERROR(context,
                "TruncationErrorAdjustmentLimit must be finite and "
                "nonpositive and leave a representable positive target.");
  }
  if (spatial_target_.has_value()) {
    for (const auto& gaussian : spatial_target_->gaussians) {
      if (log10(gaussian.target) + adjustment_limit < minimum_log) {
        PARSE_ERROR(context,
                    "Gaussian target and adjustment limit must leave "
                    "a representable positive target.");
      }
    }
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
  p | settings_.minimum_time_between_changes;
  p | spatial_target_;
  p | minimum_modes_;
  p | start_time_;
  p | radius_dependent_start_;
  p | verbose_;
}
}  // namespace gh
