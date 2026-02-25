// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Kokkos/UpdateU.hpp"

#include <array>
#include <cstddef>
#include <vector>

#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace evolution::Actions::Kokkos {

template <typename System>
void UpdateU<System>::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const TimeStepper& time_stepper, const TimeStepId& time_step_id,
    const TimeDelta& time_step) {
  const auto* runge_kutta =
      dynamic_cast<const TimeSteppers::RungeKutta*>(&time_stepper);
  ASSERT(runge_kutta != nullptr,
         "Kokkos UpdateU currently supports only Runge-Kutta time steppers.");

  const size_t substep = time_step_id.substep();
  const size_t number_of_substeps = runge_kutta->number_of_substeps();
  ASSERT(substep < number_of_substeps, "Substep should be less than "
                                           << number_of_substeps << ", not "
                                           << substep);

  const std::vector<double>* coefficients = nullptr;
  if (substep == number_of_substeps - 1) {
    coefficients = &runge_kutta->butcher_tableau().result_coefficients;
  } else {
    coefficients =
        &runge_kutta->butcher_tableau().substep_coefficients[substep];
  }
  ASSERT(coefficients != nullptr, "Missing RK coefficients.");
  ASSERT(coefficients->size() <=
             static_cast<size_t>(
                 packed_evolution_state->device_derivative_history.extent(0)),
         "Not enough entries in device derivative history: coefficients="
             << coefficients->size() << " history="
             << packed_evolution_state->device_derivative_history.extent(0));

  constexpr size_t max_supported_coefficients = 8;
  ASSERT(coefficients->size() <= max_supported_coefficients,
         "Kokkos UpdateU currently supports up to "
             << max_supported_coefficients << " RK coefficients, not "
             << coefficients->size());

  ::Kokkos::Array<double, 8> coefficients_array{};
  const size_t num_coefficients = coefficients->size();
  for (size_t coeff_index = 0; coeff_index < num_coefficients; ++coeff_index) {
    coefficients_array[coeff_index] = (*coefficients)[coeff_index];
  }

  auto& device_vars = packed_evolution_state->device_variables;
  const auto u_view = device_vars.view();
  const auto u0_view = packed_evolution_state->device_step_start.view();
  const auto derivative_history =
      packed_evolution_state->device_derivative_history;
  const size_t num_points = u_view.extent(0);
  const size_t num_components = u_view.extent(1);
  ASSERT(u0_view.extent(0) == num_points,
         "Step-start point extent mismatch: u0=" << u0_view.extent(0)
                                                 << " u=" << num_points);
  ASSERT(u0_view.extent(1) == num_components,
         "Step-start component extent mismatch: u0="
             << u0_view.extent(1) << " u=" << num_components);
  ASSERT(derivative_history.extent(1) == num_points,
         "Derivative history point extent mismatch: history="
             << derivative_history.extent(1) << " u=" << num_points);
  ASSERT(derivative_history.extent(2) == num_components,
         "Derivative history component extent mismatch: history="
             << derivative_history.extent(2) << " u=" << num_components);

  std::array<size_t, 8> u_strides{};
  u_view.stride(u_strides.data());
  std::array<size_t, 8> u0_strides{};
  u0_view.stride(u0_strides.data());
  std::array<size_t, 8> derivative_history_strides{};
  derivative_history.stride(derivative_history_strides.data());

  ::Actions::Kokkos::detail::update_u_impl(
      u_view.data(), u_strides[0], u_strides[1], u0_view.data(), u0_strides[0],
      u0_strides[1], derivative_history.data(), derivative_history_strides[0],
      derivative_history_strides[1], derivative_history_strides[2],
      coefficients_array.data(), num_coefficients, time_step.value(),
      num_points, num_components);
}

template struct UpdateU<ScalarWave::System<3>>;

}  // namespace evolution::Actions::Kokkos
