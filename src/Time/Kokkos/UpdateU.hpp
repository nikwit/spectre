// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace Actions::Kokkos {

template <typename System, template <typename> class DeviceVariablesTag,
          template <typename> class DeviceStepStartTag,
          template <typename> class DeviceDerivativeHistoryTag>
struct UpdateU {
 private:
  using device_variables_tag = DeviceVariablesTag<System>;
  using device_step_start_tag = DeviceStepStartTag<System>;
  using device_derivative_history_tag = DeviceDerivativeHistoryTag<System>;

 public:
  using return_tags = tmpl::list<device_variables_tag>;
  using argument_tags =
      tmpl::list<::Tags::TimeStepper<TimeStepper>, ::Tags::TimeStepId,
                 ::Tags::TimeStep, device_step_start_tag,
                 device_derivative_history_tag>;

  static void apply(
      const gsl::not_null<typename device_variables_tag::type*> device_vars,
      const TimeStepper& time_stepper, const TimeStepId& time_step_id,
      const TimeDelta& time_step,
      const typename device_step_start_tag::type& device_step_start,
      const typename device_derivative_history_tag::type&
          device_derivative_history) {
    const auto* runge_kutta =
        dynamic_cast<const TimeSteppers::RungeKutta*>(&time_stepper);
    ASSERT(runge_kutta != nullptr,
           "Kokkos UpdateU currently supports only Runge-Kutta time "
           "steppers.");

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
               static_cast<size_t>(device_derivative_history.extent(0)),
           "Not enough entries in device derivative history: coefficients="
               << coefficients->size()
               << " history=" << device_derivative_history.extent(0));

    constexpr size_t max_supported_coefficients = 8;
    ASSERT(coefficients->size() <= max_supported_coefficients,
           "Kokkos UpdateU currently supports up to "
               << max_supported_coefficients << " RK coefficients, not "
               << coefficients->size());

    constexpr size_t number_of_components =
        device_variables_tag::type::number_of_independent_components;
    ::Kokkos::Array<double, max_supported_coefficients> coefficients_array{};
    const size_t num_coefficients = coefficients->size();
    for (size_t coeff_index = 0; coeff_index < num_coefficients;
         ++coeff_index) {
      coefficients_array[coeff_index] = (*coefficients)[coeff_index];
    }

    const double dt = time_step.value();
    const size_t num_points = device_vars->number_of_grid_points();
    auto u_view = device_vars->view();
    const auto u0_view = device_step_start.view();
    const auto deriv_history = device_derivative_history;
    ::Kokkos::parallel_for(
        "KokkosUpdateUFused", num_points, KOKKOS_LAMBDA(const int i) {
          for (size_t c = 0; c < number_of_components; ++c) {
            double weighted_sum = 0.0;
            for (size_t coeff_index = 0; coeff_index < num_coefficients;
                 ++coeff_index) {
              weighted_sum += coefficients_array[coeff_index] *
                              deriv_history(coeff_index, i, c);
            }
            u_view(i, c) = u0_view(i, c) + dt * weighted_sum;
          }
        });
  }
};

}  // namespace Actions::Kokkos
