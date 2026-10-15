// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/BoundaryConditions/Cartoon.hpp"
#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DemandOutgoingCharSpeeds.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletMinkowski.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubePhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplay.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplayTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::BoundaryConditions {
/// Typelist of standard BoundaryConditions
template <size_t Dim>
using standard_boundary_conditions = tmpl::append<
    tmpl::list<domain::BoundaryConditions::Cartoon<BoundaryCondition<Dim>>,
               ConstraintPreservingBjorhus<Dim>, DemandOutgoingCharSpeeds<Dim>,
               DirichletAnalytic<Dim>, DirichletMinkowski<Dim>,
               domain::BoundaryConditions::Periodic<BoundaryCondition<Dim>>,
               WorldtubeTypeD<Dim>>,
    tmpl::conditional_t<Dim == 3,
                        tmpl::list<WorldtubeReplay, WorldtubeReplayTypeD,
                                   WorldtubePhysicalGhost>,
                        tmpl::list<>>>;
}  // namespace gh::BoundaryConditions
