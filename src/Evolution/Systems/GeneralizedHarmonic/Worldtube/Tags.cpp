// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"

#include <pup.h>
#include <pup_stl.h>

namespace gh::Worldtube {
void MatcherConfig::pup(PUP::er& p) {
  p | mass;
  p | center;
  p | center_velocity;
  p | trace_strain_pin;
  p | fit_l_max;
  p | fit_interval;
}

bool operator==(const MatcherConfig& lhs, const MatcherConfig& rhs) {
  return lhs.mass == rhs.mass and lhs.center == rhs.center and
         lhs.center_velocity == rhs.center_velocity and
         lhs.trace_strain_pin == rhs.trace_strain_pin and
         lhs.fit_l_max == rhs.fit_l_max and
         lhs.fit_interval == rhs.fit_interval;
}

bool operator!=(const MatcherConfig& lhs, const MatcherConfig& rhs) {
  return not(lhs == rhs);
}

void MapParameterData::pup(PUP::er& pupper) {
  pupper | last_fit_time;
  pupper | previous_fit_time;
  pupper | p;
  pupper | p_previous;
  pupper | pdot;
  pupper | valid;
}
}  // namespace gh::Worldtube
