// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Kokkos/RecordTimeStepperData.hpp"

#include <array>
#include <cstddef>

#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace evolution::Actions::Kokkos {

template <typename System>
void RecordTimeStepperData<System>::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const TimeStepId& time_step_id) {
  auto& derivative_history = packed_evolution_state->device_derivative_history;
  const auto& device_dt = packed_evolution_state->device_dt_variables;

  const size_t substep = time_step_id.substep();
  ASSERT(substep < static_cast<size_t>(derivative_history.extent(0)),
         "Substep " << substep << " exceeds derivative history size "
                    << derivative_history.extent(0));

  const auto derivative_history_view = derivative_history;
  const auto dt_view = device_dt.view();
  const size_t num_points = dt_view.extent(0);
  const size_t num_components = dt_view.extent(1);
  ASSERT(derivative_history_view.extent(1) == num_points,
         "Derivative history point extent mismatch: history="
             << derivative_history_view.extent(1) << " dt=" << num_points);
  ASSERT(derivative_history_view.extent(2) == num_components,
         "Derivative history component extent mismatch: history="
             << derivative_history_view.extent(2) << " dt=" << num_components);
  std::array<size_t, 8> derivative_history_strides{};
  derivative_history_view.stride(derivative_history_strides.data());
  std::array<size_t, 8> dt_strides{};
  dt_view.stride(dt_strides.data());

  ::Actions::Kokkos::detail::record_time_stepper_data_impl(
      derivative_history_view.data(), derivative_history_strides[0],
      derivative_history_strides[1], derivative_history_strides[2],
      dt_view.data(), dt_strides[0], dt_strides[1], substep, num_points,
      num_components);
}

template struct RecordTimeStepperData<ScalarWave::System<3>>;

}  // namespace evolution::Actions::Kokkos
