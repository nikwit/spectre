// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Kokkos/CleanHistory.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace evolution::Actions::Kokkos {

template <typename System>
void CleanHistory<System>::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const TimeStepper& time_stepper, const TimeStepId& time_step_id) {
  const auto* runge_kutta =
      dynamic_cast<const TimeSteppers::RungeKutta*>(&time_stepper);
  ASSERT(runge_kutta != nullptr,
         "Kokkos CleanHistory currently supports only Runge-Kutta steppers.");

  if (time_step_id.substep() == runge_kutta->number_of_substeps() - 1) {
    auto& device_step_start = packed_evolution_state->device_step_start;
    auto& device_vars = packed_evolution_state->device_variables;
    ::Kokkos::deep_copy(device_step_start.view(), device_vars.view());
  }
}

template struct CleanHistory<ScalarWave::System<3>>;
template struct CleanHistory<gh::System<3>>;

}  // namespace evolution::Actions::Kokkos
