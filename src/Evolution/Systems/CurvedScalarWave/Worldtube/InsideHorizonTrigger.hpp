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
      CurvedScalarWave::Worldtube::Tags::ExcisionSphere<Dim>,
      CurvedScalarWave::Worldtube::Tags::PowerLawParams>;

  bool operator()(const std::array<tnsr::I<double, Dim, Frame::Inertial>, 2>&
                      position_and_velocity,
                  const ExcisionSphere<Dim>& excision_sphere,
                  const std::array<double, 2>& power_law_params) const {
    const double orbit_radius = get(magnitude(position_and_velocity[0]));
    const double worldtube_radius =
        CurvedScalarWave::Worldtube::broken_power_shrink(
            orbit_radius, power_law_params.at(0), power_law_params.at(1));
    return orbit_radius + worldtube_radius < 1.99 ? true : false;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override {}
};
}  // namespace Triggers
