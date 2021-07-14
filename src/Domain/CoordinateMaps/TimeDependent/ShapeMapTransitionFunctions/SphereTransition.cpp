// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <optional>
#include <pup.h>

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"

namespace domain::CoordinateMaps {

SphereTransition::SphereTransition(double r_min, double r_max)
    : r_min_(r_min),
      r_max_(r_max),
      a_(-1. / (r_max - r_min)),
      b_(r_max / (r_max - r_min)) {
  if (r_min <= 0.) {
    ERROR("The minimum radius must be greater than 0 but is " << r_min);
  }
  if (r_max <= r_min) {
    ERROR(
        "The maximum radius must be greater than the minimum radius but "
        "r_max =  "
        << r_max << ", and r_min = " << r_min);
  }
}

double SphereTransition::operator()(
    const std::array<double, 3>& source_coords) const {
  return call_impl<double>(source_coords);
}
DataVector SphereTransition::operator()(
    const std::array<DataVector, 3>& source_coords) const {
  return call_impl<DataVector>(source_coords);
}

std::optional<double> SphereTransition::original_radius_over_radius(
    const std::array<double, 3>& target_coords, double distorted_radius) const {
  const double mag = magnitude(target_coords);
  const double denom = 1. - distorted_radius * a_;
  // prevent zero division
  if (mag == 0. or denom == 0.) {
    return std::nullopt;
  }
  const double original_radius = (mag + distorted_radius * b_) / denom;

  return original_radius >= r_min_ and original_radius <= r_max_
             ? std::optional<double>{original_radius / mag}
             : std::nullopt;
}

double SphereTransition::map_over_radius(
    const std::array<double, 3>& source_coords) const {
  return map_over_radius_impl<double>(source_coords);
}
DataVector SphereTransition::map_over_radius(
    const std::array<DataVector, 3>& source_coords) const {
  return map_over_radius_impl<DataVector>(source_coords);
}

std::array<double, 3> SphereTransition::gradient(
    const std::array<double, 3>& source_coords) const {
  return gradient_impl<double>(source_coords);
}
std::array<DataVector, 3> SphereTransition::gradient(
    const std::array<DataVector, 3>& source_coords) const {
  return gradient_impl<DataVector>(source_coords);
}

template <typename T>
T SphereTransition::call_impl(const std::array<T, 3>& source_coords) const {
  const T mag = magnitude(source_coords);
  check_magnitudes(mag);
  return a_ + b_ / mag;
}

template <typename T>
T SphereTransition::map_over_radius_impl(
    const std::array<T, 3>& source_coords) const {
  const T mag = magnitude(source_coords);
  check_magnitudes(mag);
  return a_ / mag + b_ / square(mag);
}

template <typename T>
std::array<T, 3> SphereTransition::gradient_impl(
    const std::array<T, 3>& source_coords) const {
  const T mag = magnitude(source_coords);
  check_magnitudes(mag);
  return -b_ * source_coords / cube(mag);
}

bool SphereTransition::operator==(
    const ShapeMapTransitionFunction& other) const {
  try {
    const auto& derived = dynamic_cast<const SphereTransition&>(other);
    // no need to check `a_` and `b_` as they are uniquely determined by
    // `r_min_` and `r_max_`.
    return this->r_min_ == derived.r_min_ and this->r_max_ == derived.r_max_;
  } catch (const std::bad_cast& ex) {
    return false;
  }
}

bool SphereTransition::operator!=(
    const ShapeMapTransitionFunction& other) const {
  return not(*this == other);
}

// checks that the magnitudes are all between `r_min_` and `r_max_`
template <typename T>
void SphereTransition::check_magnitudes(const T& mag) const {
  for (size_t i = 0; i < get_size(mag); ++i) {
    if (get_element(mag, i) < r_min_ or get_element(mag, i) > r_max_) {
      ERROR(
          "The sphere transition map was called with coordinates outside the "
          "set minimum and maximum radius. The minimum radius is "
          << r_min_ << ", the maximum radius is " << r_max_
          << ". The requested point has magnitude: " << get_element(mag, i));
    }
  }
}

void SphereTransition::pup(PUP::er& p) {
  ShapeMapTransitionFunction::pup(p);
  p | r_min_;
  p | r_max_;
  p | a_;
  p | b_;
};

SphereTransition::SphereTransition(CkMigrateMessage* const msg)
    : ShapeMapTransitionFunction(msg) {}

PUP::able::PUP_ID SphereTransition::my_PUP_ID = 0;

}  // namespace domain::CoordinateMaps
