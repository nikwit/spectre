// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <memory>
#include <pup.h>
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"

namespace domain::FunctionsOfTime {
/// Prescribed shape-map monopole sqrt(4 pi) (R_grid - R_inertial / a(t)).
/// Owns an immutable expansion history; no controller updates are permitted.
class FixedInertialRadius final : public FunctionOfTime {
 public:
  FixedInertialRadius() = default;
  FixedInertialRadius(double grid_radius, double inertial_radius,
                      std::unique_ptr<FunctionOfTime> expansion);
  explicit FixedInertialRadius(CkMigrateMessage*) {}
  WRAPPED_PUPable_decl_template(FixedInertialRadius);
  std::unique_ptr<FunctionOfTime> get_clone() const override;
  std::array<double, 2> time_bounds() const override;
  double expiration_after(double time) const override;
  void truncate_at_time(double time) override;
  std::array<DataVector, 1> func(double time) const override;
  std::array<DataVector, 2> func_and_deriv(double time) const override;
  std::array<DataVector, 3> func_and_2_derivs(double time) const override;
  void pup(PUP::er& p) override;

 private:
  double grid_radius_{0.};
  double inertial_radius_{0.};
  std::unique_ptr<FunctionOfTime> expansion_{};
};
}  // namespace domain::FunctionsOfTime
