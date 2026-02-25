// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/UpdateUKokkos3D.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct UpdateUBatched {
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags = tmpl::list<::Tags::TimeStepper<TimeStepper>,
                                   ::Tags::TimeStepId, ::Tags::TimeStep>;

  static void apply(const gsl::not_null<
                        ScalarWave::Batched::Tags::PackedEvolutionState::type*>
                        packed_evolution_state,
                    const TimeStepper& time_stepper,
                    const TimeStepId& time_step_id,
                    const TimeDelta& time_step) {
    auto& device_vars = packed_evolution_state->device_variables;
    ScalarWave::Actions::UpdateUKokkos3D::apply(
        make_not_null(&device_vars), time_stepper, time_step_id, time_step,
        packed_evolution_state->device_step_start,
        packed_evolution_state->device_derivative_history);
  }
};

}  // namespace ScalarWave::Actions
