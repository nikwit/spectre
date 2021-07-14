// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/RegisterDerivedWithCharm.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"

namespace domain::CoordinateMaps {
void register_derived_with_charm() noexcept {
  Parallel::register_derived_classes_with_charm<ShapeMapTransitionFunction>();
}
}  // namespace domain::CoordinateMaps
