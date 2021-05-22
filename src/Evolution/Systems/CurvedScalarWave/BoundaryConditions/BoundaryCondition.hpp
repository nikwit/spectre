// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <pup.h>

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace CurvedScalarWave::BoundaryConditions {
template <size_t Dim>
class ConstraintPreserving;
template <size_t Dim>
class ConstraintPreservingBaylissTurkel;
template <size_t Dim>
class Freezing;
template <size_t Dim>
class Outflowing;

}  // namespace CurvedScalarWave::BoundaryConditions
/// \endcond

/// \brief Boundary conditions for the scalar wave system
namespace CurvedScalarWave::BoundaryConditions {
/// \brief The base class off of which all boundary conditions must inherit
template <size_t Dim>
class BoundaryCondition : public domain::BoundaryConditions::BoundaryCondition {
 public:
  using creatable_classes =
      tmpl::list<ConstraintPreserving<Dim>, Freezing<Dim>,
                 ConstraintPreservingBaylissTurkel<Dim>, Outflowing<Dim>,
                 domain::BoundaryConditions::Periodic<BoundaryCondition<Dim>>>;

  BoundaryCondition() = default;
  BoundaryCondition(BoundaryCondition&&) noexcept = default;
  BoundaryCondition& operator=(BoundaryCondition&&) noexcept = default;
  BoundaryCondition(const BoundaryCondition&) = default;
  BoundaryCondition& operator=(const BoundaryCondition&) = default;
  ~BoundaryCondition() override = default;
  explicit BoundaryCondition(CkMigrateMessage* msg) noexcept;

  void pup(PUP::er& p) override;
};
}  // namespace CurvedScalarWave::BoundaryConditions
