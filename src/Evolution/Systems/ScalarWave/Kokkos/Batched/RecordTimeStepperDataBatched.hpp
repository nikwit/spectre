// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/RecordTimeStepperDataKokkos3D.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct RecordTimeStepperDataBatched {
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags = tmpl::list<::Tags::TimeStepId>;

  static void apply(const gsl::not_null<
                        ScalarWave::Batched::Tags::PackedEvolutionState::type*>
                        packed_evolution_state,
                    const TimeStepId& time_step_id) {
    auto& deriv_history = packed_evolution_state->device_derivative_history;
    ScalarWave::Actions::RecordTimeStepperDataKokkos3D::apply(
        make_not_null(&deriv_history), time_step_id,
        packed_evolution_state->device_dt_variables);
  }
};

}  // namespace ScalarWave::Actions
