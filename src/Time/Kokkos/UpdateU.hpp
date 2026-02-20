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

namespace detail {

void update_u_impl(double* u_data, size_t u_stride_0, size_t u_stride_1,
                   const double* u0_data, size_t u0_stride_0,
                   size_t u0_stride_1, const double* deriv_history_data,
                   size_t deriv_history_stride_0, size_t deriv_history_stride_1,
                   size_t deriv_history_stride_2, const double* coefficients,
                   size_t num_coefficients, double dt, size_t num_points,
                   size_t num_components);

}  // namespace detail

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

    ::Kokkos::Array<double, 8> coefficients_array{};
    const size_t num_coefficients = coefficients->size();
    for (size_t coeff_index = 0; coeff_index < num_coefficients;
         ++coeff_index) {
      coefficients_array[coeff_index] = (*coefficients)[coeff_index];
    }

    const double dt = time_step.value();
    const auto u_view = device_vars->view();
    const auto u0_view = device_step_start.view();
    const auto deriv_history = device_derivative_history;
    const size_t num_points = u_view.extent(0);
    const size_t num_components = u_view.extent(1);
    ASSERT(u0_view.extent(0) == num_points,
           "Step-start point extent mismatch: u0=" << u0_view.extent(0)
                                                   << " u=" << num_points);
    ASSERT(u0_view.extent(1) == num_components,
           "Step-start component extent mismatch: u0="
               << u0_view.extent(1) << " u=" << num_components);
    ASSERT(deriv_history.extent(1) == num_points,
           "Derivative history point extent mismatch: history="
               << deriv_history.extent(1) << " u=" << num_points);
    ASSERT(deriv_history.extent(2) == num_components,
           "Derivative history component extent mismatch: history="
               << deriv_history.extent(2) << " u=" << num_components);
    std::array<size_t, 8> u_strides{};
    u_view.stride(u_strides.data());
    std::array<size_t, 8> u0_strides{};
    u0_view.stride(u0_strides.data());
    std::array<size_t, 8> deriv_history_strides{};
    deriv_history.stride(deriv_history_strides.data());

    detail::update_u_impl(u_view.data(), u_strides[0], u_strides[1],
                          u0_view.data(), u0_strides[0], u0_strides[1],
                          deriv_history.data(), deriv_history_strides[0],
                          deriv_history_strides[1], deriv_history_strides[2],
                          coefficients_array.data(), num_coefficients, dt,
                          num_points, num_components);
  }
};

}  // namespace Actions::Kokkos
