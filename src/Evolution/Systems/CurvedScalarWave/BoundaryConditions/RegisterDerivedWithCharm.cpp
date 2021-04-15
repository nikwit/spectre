// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/RegisterDerivedWithCharm.hpp"

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/BoundaryCondition.hpp"
//#include
//"Evolution/Systems/CurvedScalarWave/
// BoundaryConditions/ConstraintPreservingSphericalRadiation.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/DirichletAnalytic.hpp"
//#include
//"Evolution/Systems/CurvedScalarWave/
// BoundaryConditions/SphericalRadiation.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"

namespace CurvedScalarWave::BoundaryConditions {
void register_derived_with_charm() noexcept {
  Parallel::register_classes_in_list<
      typename BoundaryCondition<1>::creatable_classes>();
  Parallel::register_classes_in_list<
      typename BoundaryCondition<2>::creatable_classes>();
  Parallel::register_classes_in_list<
      typename BoundaryCondition<3>::creatable_classes>();
}
}  // namespace CurvedScalarWave::BoundaryConditions
