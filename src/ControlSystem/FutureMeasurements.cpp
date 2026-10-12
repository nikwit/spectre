// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/FutureMeasurements.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <pup_stl.h>

#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/StdHelpers.hpp"

namespace control_system {
FutureMeasurements::FutureMeasurements(const size_t measurements_per_update,
                                       const double first_measurement_time)
    : measurements_{std::numeric_limits<double>::signaling_NaN(),
                    first_measurement_time},
      measurements_until_update_(measurements_per_update - 1),
      measurements_per_update_(measurements_per_update) {
  ASSERT(measurements_per_update > 0, "Cannot update without measurements.");
}

std::optional<double> FutureMeasurements::next_measurement() const {
  if (measurements_.size() <= 1) {
    return std::nullopt;
  }
  // First entry is an old measurement.
  return measurements_[1];
}

std::optional<double> FutureMeasurements::next_update() const {
  if (measurements_.size() <= measurements_until_update_ + 1) {
    return std::nullopt;
  }
  return {measurements_[measurements_until_update_ + 1]};
}

void FutureMeasurements::pop_front() {
  ASSERT(measurements_.size() > 1, "Popped an empty container.");
  measurements_.pop_front();
  if (measurements_until_update_ == 0) {
    measurements_until_update_ = measurements_per_update_ - 1;
  } else {
    --measurements_until_update_;
  }
}

void FutureMeasurements::update(
    const domain::FunctionsOfTime::FunctionOfTime& measurement_timescale) {
  ASSERT(not measurements_.empty(), "Don't know current measurement time.");
  // Limit the number of stored measurements to a few cycles.  We
  // shouldn't need more than one cycle, but tuning the exact size
  // isn't necessary.
  while (measurement_timescale.time_bounds()[1] >= measurements_.back() and
         measurements_.size() < 3 * measurements_per_update_) {
    measurements_.emplace_back(
        measurements_.back() +
        measurement_timescale.func(measurements_.back())[0][0]);
  }
}

void FutureMeasurements::pup(PUP::er& p) {
  p | measurements_;
  p | measurements_until_update_;
  p | measurements_per_update_;
}

bool operator==(const FutureMeasurements& lhs, const FutureMeasurements& rhs) {
  if (lhs.measurements_.size() != rhs.measurements_.size() or
      lhs.measurements_until_update_ != rhs.measurements_until_update_ or
      lhs.measurements_per_update_ != rhs.measurements_per_update_) {
    return false;
  }
  if (lhs.measurements_.empty()) {
    return true;
  }
  // Before the first measurement, the old-measurement slot contains a NaN
  // sentinel. This persists for inactive control systems. Treat matching
  // sentinels as equal so AMR can verify that joining children agree.
  return ((std::isnan(lhs.measurements_.front()) and
           std::isnan(rhs.measurements_.front())) or
          lhs.measurements_.front() == rhs.measurements_.front()) and
         std::equal(std::next(lhs.measurements_.begin()),
                    lhs.measurements_.end(),
                    std::next(rhs.measurements_.begin()));
}

bool operator!=(const FutureMeasurements& lhs, const FutureMeasurements& rhs) {
  return not(lhs == rhs);
}
}  // namespace control_system
