// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/RegisterDerivedWithCharm.hpp"

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/ConstraintPreserving.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/ConstraintPreservingBaylissTurkel.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Freezing.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Outflowing.hpp"


//#include
//"Evolution/Systems/CurvedScalarWave/BoundaryConditions/DirichletAnalytic.hpp"

//#include
//"Evolution/Systems/CurvedScalarWave/
// BoundaryConditions/SphericalRadiation.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"

namespace CurvedScalarWave::BoundaryConditions {
void register_derived_with_charm() noexcept {
  Parallel::register_derived_classes_with_charm<BoundaryCondition<1>>();
  Parallel::register_derived_classes_with_charm<BoundaryCondition<2>>();
  Parallel::register_derived_classes_with_charm<BoundaryCondition<3>>();
}
}  // namespace CurvedScalarWave::BoundaryConditions
