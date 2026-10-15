// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Domain/FunctionsOfTime/FixedInertialRadius.hpp"

#include <cmath>
#include <numbers>
#include <utility>
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/PupStlCpp11.hpp"

namespace domain::FunctionsOfTime {
FixedInertialRadius::FixedInertialRadius(
    const double grid_radius, const double inertial_radius,
    std::unique_ptr<FunctionOfTime> expansion)
    : grid_radius_(grid_radius),
      inertial_radius_(inertial_radius),
      expansion_(std::move(expansion)) {
  if (not std::isfinite(grid_radius_) or grid_radius_ <= 0. or
      not std::isfinite(inertial_radius_) or inertial_radius_ <= 0. or
      expansion_ == nullptr) {
    ERROR(
        "FixedInertialRadius requires positive finite radii and an expansion "
        "history.");
  }
}
std::unique_ptr<FunctionOfTime> FixedInertialRadius::get_clone() const {
  return std::make_unique<FixedInertialRadius>(grid_radius_, inertial_radius_,
                                               expansion_->get_clone());
}
std::array<double, 2> FixedInertialRadius::time_bounds() const {
  return expansion_->time_bounds();
}
double FixedInertialRadius::expiration_after(const double time) const {
  return expansion_->expiration_after(time);
}
void FixedInertialRadius::truncate_at_time(const double time) {
  expansion_->truncate_at_time(time);
}
std::array<DataVector, 3> FixedInertialRadius::func_and_2_derivs(
    const double time) const {
  const auto a = expansion_->func_and_2_derivs(time);
  if (a[0].size() != 1 or not std::isfinite(a[0][0]) or a[0][0] <= 0.) {
    ERROR("FixedInertialRadius requires a positive scalar expansion.");
  }
  const double inverse_a = 1. / a[0][0];
  const double factor = sqrt(4. * std::numbers::pi);
  return {
      DataVector{factor * (grid_radius_ - inertial_radius_ * inverse_a)},
      DataVector{factor * inertial_radius_ * a[1][0] * inverse_a * inverse_a},
      DataVector{factor * inertial_radius_ * inverse_a * inverse_a *
                 (a[2][0] - 2. * a[1][0] * a[1][0] * inverse_a)}};
}
std::array<DataVector, 1> FixedInertialRadius::func(const double time) const {
  return {func_and_2_derivs(time)[0]};
}
std::array<DataVector, 2> FixedInertialRadius::func_and_deriv(
    const double time) const {
  const auto values = func_and_2_derivs(time);
  return {values[0], values[1]};
}
void FixedInertialRadius::pup(PUP::er& p) {
  FunctionOfTime::pup(p);
  p | grid_radius_;
  p | inertial_radius_;
  p | expansion_;
}
PUP::able::PUP_ID FixedInertialRadius::my_PUP_ID = 0;  // NOLINT
}  // namespace domain::FunctionsOfTime
