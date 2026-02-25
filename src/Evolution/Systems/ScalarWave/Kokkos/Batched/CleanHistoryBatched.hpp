// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/CleanHistoryKokkos.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct CleanHistoryBatched {
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags =
      tmpl::list<::Tags::TimeStepper<TimeStepper>, ::Tags::TimeStepId>;

  static void apply(const gsl::not_null<
                        ScalarWave::Batched::Tags::PackedEvolutionState::type*>
                        packed_evolution_state,
                    const TimeStepper& time_stepper,
                    const TimeStepId& time_step_id) {
    auto& device_step_start = packed_evolution_state->device_step_start;
    auto& device_vars = packed_evolution_state->device_variables;
    ScalarWave::Actions::CleanHistoryKokkos<ScalarWave::System<3>>::apply(
        make_not_null(&device_step_start), make_not_null(&device_vars),
        time_stepper, time_step_id);
  }
};

}  // namespace ScalarWave::Actions
