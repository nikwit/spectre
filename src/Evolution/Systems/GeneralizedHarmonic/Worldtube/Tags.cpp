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
  p | fit_center_offset;
  p | rate_ode;
  p | second_order_ode;
  p | stepper_ode;
  p | gauge_damping;
  p | uplus_anchor;
  p | fit_uplus;
  p | kretschmann_trace_pin;
  p | trace_pin_interval;
  p | fit_trace_strain;
  p | fit_velocity;
  p | fit_bulk_boost;
  p | fit_exact_frame;
  p | fit_velocity_separately;
  p | pin_symmetric_factor;
  p | fit_radial_derivative;
  p | radial_derivative_weight;
  p | centre_advection;
  p | spatial_monopole_weight;
  p | uplus_block_weights;
  p | fit_radial_index;
  p | excision_sphere_name;
}

bool operator==(const MatcherConfig& lhs, const MatcherConfig& rhs) {
  return lhs.mass == rhs.mass and lhs.center == rhs.center and
         lhs.center_velocity == rhs.center_velocity and
         lhs.trace_strain_pin == rhs.trace_strain_pin and
         lhs.fit_l_max == rhs.fit_l_max and
         lhs.fit_interval == rhs.fit_interval and
         lhs.fit_center_offset == rhs.fit_center_offset and
         lhs.rate_ode == rhs.rate_ode and
         lhs.second_order_ode == rhs.second_order_ode and
         lhs.stepper_ode == rhs.stepper_ode and
         lhs.gauge_damping == rhs.gauge_damping and
         lhs.uplus_anchor == rhs.uplus_anchor and
         lhs.fit_uplus == rhs.fit_uplus and
         lhs.kretschmann_trace_pin == rhs.kretschmann_trace_pin and
         lhs.trace_pin_interval == rhs.trace_pin_interval and
         lhs.fit_trace_strain == rhs.fit_trace_strain and
         lhs.fit_velocity == rhs.fit_velocity and
         lhs.fit_bulk_boost == rhs.fit_bulk_boost and
         lhs.fit_exact_frame == rhs.fit_exact_frame and
         lhs.fit_velocity_separately == rhs.fit_velocity_separately and
         lhs.pin_symmetric_factor == rhs.pin_symmetric_factor and
         lhs.fit_radial_derivative == rhs.fit_radial_derivative and
         lhs.radial_derivative_weight == rhs.radial_derivative_weight and
         lhs.centre_advection == rhs.centre_advection and
         lhs.spatial_monopole_weight == rhs.spatial_monopole_weight and
         lhs.uplus_block_weights == rhs.uplus_block_weights and
         lhs.fit_radial_index == rhs.fit_radial_index and
         lhs.excision_sphere_name == rhs.excision_sphere_name;
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
  pupper | pddot;
  pupper | bulk_velocity;
  pupper | worldtube_center;
  pupper | worldtube_center_at_last_fit;
  pupper | worldtube_center_valid;
  pupper | center_offset;
  pupper | exact_frame_theta;
  pupper | exact_frame_center_velocity;
  pupper | exact_frame_valid;
  pupper | ode_state;
  pupper | ode_history;
  pupper | ode_step_start;
  pupper | ode_step_id;
  pupper | trace_pin_value;
  pupper | trace_pin_time;
  pupper | gb_dipole;
  pupper | gb_dipole_previous;
  pupper | gb_dipole_velocity;
  pupper | gb_dipole_time;
  pupper | gb_dipole_time_previous;
  pupper | anchor_p;
  pupper | anchor_p_previous;
  pupper | anchor_time;
  pupper | anchor_time_previous;
  pupper | anchor_valid;
  pupper | valid;
}
}  // namespace gh::Worldtube
