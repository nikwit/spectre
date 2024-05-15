// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Worldtube.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube::Actions {

/*!
 * \brief Waits for the data from all neighboring elements and changes the slab
 * size if a change in the global time step is detected.
 * \details We check the slab size of the time step id sent by the elements. If
 * this is different from the slab size currently used by the worldtube
 * singleton, we assume a global slab size change has occurred in the elements
 * and adjust the worldtube slab size accordingly.
 */
struct UpdateFunctionsOfTime {
  static constexpr size_t Dim = 3;
  // using compute_tags =
  // tmpl::list<Tags::InertialParticlePositionCompute<Dim>>;
  using inbox_tags = tmpl::list<>;
  using simple_tags =
      tmpl::list<Tags::ExpirationTime, Tags::WorldtubeRadiusAndVelocity>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& time = db::get<::Tags::Time>(box);
    const std::string rot_function_of_time_name = "Rotation";
    const std::string expansion_fot_name = "Expansion";
    const std::string size_a_fot_name = "SizeA";
    const std::string size_b_fot_name = "SizeB";

    const double current_fot_expiration_time =
        db::get<Tags::ExpirationTime>(box);
    if (time > current_fot_expiration_time) {
      const auto& inertial_particle_position = db::get<Tags::Position>(box);
      tnsr::I<double, Dim> particle_pos_double{};
      particle_pos_double.get(0) = inertial_particle_position.get(0)[0];
      particle_pos_double.get(1) = inertial_particle_position.get(1)[0];
      particle_pos_double.get(2) = inertial_particle_position.get(2)[0];

      const auto& particle_velocity = db::get<Tags::Velocity>(box);

      const double& x = get<0>(inertial_particle_position)[0];
      const double& y = get<1>(inertial_particle_position)[0];
      const double& xdot = get<0>(particle_velocity)[0];
      const double& ydot = get<1>(particle_velocity)[0];
      const double r = hypot(x, y);
      const double angle = atan2(y, x);
      const double radial_vel = (xdot * x + ydot * y) / r;
      const double angular_vel = (x * ydot - y * xdot) / square(r);

      const auto [envelope_radius, object_a_radius, object_b_radius] =
          db::get<Tags::EnvelopeAndObjectRadii>(box);

      const double grid_radius_particle =
          get(magnitude(db::get<Tags::ExcisionSphere<3>>(box).center()));
      DataVector angular_update(3, 0.);
      DataVector expansion_update(3, 0.);
      DataVector compression_update_a(3, 0.);
      DataVector compression_update_b(3, 0.);

      const double sqrt_4_pi = sqrt(4. * M_PI);
      angular_update.at(0) = angle;
      angular_update.at(1) = angular_vel;
      expansion_update.at(0) =
          (1 - r / grid_radius_particle) * sqrt_4_pi * envelope_radius;
      expansion_update.at(1) =
          -radial_vel / grid_radius_particle * sqrt_4_pi * envelope_radius;
      const auto& [amp, rb] = db::get<Tags::PowerLawParams>(box);
      const double worldtube_radius_factor =
          broken_power_shrink(r, amp, rb) / object_a_radius;
      const double worldtube_radius_factor_derivative =
          broken_power_shrink_derivative(r, amp, rb, radial_vel) /
          object_a_radius;

      const double bh_radius_factor = sqrt(r / grid_radius_particle);
      const double bh_radius_factor_derivative =
          0.5 * sqrt(grid_radius_particle / r) * radial_vel /
          grid_radius_particle;

      const double factor =
          1. / (1. - expansion_update.at(0) / (sqrt_4_pi * envelope_radius));

      compression_update_a.at(0) =
          sqrt_4_pi * object_a_radius * (1. - worldtube_radius_factor * factor);
      compression_update_a.at(1) =
          -object_a_radius * factor *
          (sqrt_4_pi * worldtube_radius_factor_derivative +
           worldtube_radius_factor * expansion_update.at(1) * factor /
               envelope_radius);

      compression_update_b.at(0) =
          sqrt_4_pi * object_b_radius * (1. - bh_radius_factor * factor);
      compression_update_b.at(1) = -object_b_radius * factor *
                                   (sqrt_4_pi * bh_radius_factor_derivative +
                                    bh_radius_factor * expansion_update.at(1) *
                                        factor / envelope_radius);

      const double new_fot_expiration_time =
          time +
          0.5 * (db::get<::Tags::Next<::Tags::TimeStepId>>(box).substep_time() -
                 time);
      ::Initialization::mutate_assign<simple_tags>(
          make_not_null(&box), new_fot_expiration_time,
          std::array<double, 2>{
              {object_a_radius * worldtube_radius_factor,
               object_a_radius * worldtube_radius_factor_derivative}});

      /*Parallel::printf(MakeString{} << "Mutating Time from "
                                    << current_fot_expiration_time << " to "
                                    << new_fot_expiration_time << "\n");*/
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, rot_function_of_time_name, current_fot_expiration_time,
          angular_update, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, expansion_fot_name, current_fot_expiration_time,
          expansion_update, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, size_a_fot_name, current_fot_expiration_time,
          compression_update_a, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, size_b_fot_name, current_fot_expiration_time,
          compression_update_b, new_fot_expiration_time);
    } else {
      /*Parallel::printf(MakeString{} << "Not mutating Time at " << time
                                    << " with expiration time "
                                    << current_fot_expiration_time << "\n");*/
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
