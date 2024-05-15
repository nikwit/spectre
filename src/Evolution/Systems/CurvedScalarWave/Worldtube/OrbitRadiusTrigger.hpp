// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
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

class OrbitRadius : public Trigger {
 public:
  /// \cond
  OrbitRadius() = default;
  explicit OrbitRadius(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(OrbitRadius);  // NOLINT
  /// \endcond

  struct OptionTags {
    struct Radii {
      using type = std::vector<double>;
      static constexpr Options::String help = "Orbit radii to trigger at";
    };
  };

  static constexpr Options::String help =
      "Trigger in intervals surrounding particular times.";
  using options = tmpl::list<typename OptionTags::Radii>;

  OrbitRadius(const std::vector<double>& radii) : radii_(radii) {}

  using argument_tags =
      tmpl::list<CurvedScalarWave::Worldtube::Tags::ParticlePositionVelocity<3>,
                 Tags::TimeStep>;

  bool operator()(const std::array<tnsr::I<double, 3, Frame::Inertial>, 2>&
                      position_and_velocity,
                  const TimeDelta& time_step) const {
    const auto& position = position_and_velocity[0];
    const auto& velocity = position_and_velocity[1];
    const double current_radius = get(magnitude(position));
    const double radial_velocity = (get<0>(position) * get<0>(velocity) +
                                    get<1>(position) * get<1>(velocity)) /
                                   current_radius;
    // ASSERT(radial_velocity < 0., "Particle should be inspiralling!");
    const double last_radius =
        current_radius - 1.2 * radial_velocity * time_step.value();
    for (double radius_ : radii_) {
      // factor 1.2 is for safety because the approximation is just linear
      // (Euler), triggering it multiple times doesn't matter
      if ((current_radius - radius_) * (last_radius - radius_) < 0.) {
        return true;
      }
    }
    return false;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override { p | radii_; }

 private:
  std::vector<double> radii_{};
};
}  // namespace Triggers
