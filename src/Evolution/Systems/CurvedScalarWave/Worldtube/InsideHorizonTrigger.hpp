// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Domain/Structure/ExcisionSphere.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Worldtube.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct TimeStepId;
}  // namespace Tags
/// \endcond

namespace Triggers {
/// \ingroup EventsAndTriggersGroup
/// \ingroup TimeGroup
/// Trigger based on a comparison with the slab number.

class InsideHorizon : public Trigger {
 public:
  /// \cond
  InsideHorizon() = default;
  explicit InsideHorizon(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(InsideHorizon);  // NOLINT
  /// \endcond

  static constexpr Options::String help =
      "Trigger in intervals surrounding particular times.";
  static constexpr size_t Dim = 3;
  using options = tmpl::list<>;
  using argument_tags = tmpl::list<
      CurvedScalarWave::Worldtube::Tags::ParticlePositionVelocity<Dim>,
      CurvedScalarWave::Worldtube::Tags::ExcisionSphere<Dim>>;

  bool operator()(const std::array<tnsr::I<double, Dim, Frame::Inertial>, 2>&
                      position_and_velocity,
                  const ExcisionSphere<Dim>& excision_sphere) const {
    const double orbit_radius = get(magnitude(position_and_velocity[0]));

    const double original_orbit_radius =
        get(magnitude(excision_sphere.center()));
    const double worldtube_radius_factor =
        CurvedScalarWave::Worldtube::worldtube_shrink_factor(
            orbit_radius, original_orbit_radius, 3., 2., 0.5);
    return orbit_radius + excision_sphere.radius() * worldtube_radius_factor <
                   1.99
               ? true
               : false;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {}
};
}  // namespace Triggers
