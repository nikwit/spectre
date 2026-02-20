// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace Actions::Kokkos {

namespace detail {

void record_time_stepper_data_impl(double* deriv_history_data,
                                   size_t deriv_history_stride_0,
                                   size_t deriv_history_stride_1,
                                   size_t deriv_history_stride_2,
                                   const double* dt_data, size_t dt_stride_0,
                                   size_t dt_stride_1, size_t substep,
                                   size_t num_points, size_t num_components);

}  // namespace detail

template <typename System, template <typename> class DeviceDtVariablesTag,
          template <typename> class DeviceDerivativeHistoryTag>
struct RecordTimeStepperData {
 private:
  using device_dt_variables_tag = DeviceDtVariablesTag<System>;
  using device_derivative_history_tag = DeviceDerivativeHistoryTag<System>;

 public:
  using return_tags = tmpl::list<device_derivative_history_tag>;
  using argument_tags = tmpl::list<::Tags::TimeStepId, device_dt_variables_tag>;

  static void apply(
      const gsl::not_null<typename device_derivative_history_tag::type*>
          device_derivative_history,
      const TimeStepId& time_step_id,
      const typename device_dt_variables_tag::type& device_dt) {
    const size_t substep = time_step_id.substep();
    ASSERT(substep < static_cast<size_t>(device_derivative_history->extent(0)),
           "Substep " << substep << " exceeds derivative history size "
                      << device_derivative_history->extent(0));

    const auto deriv_history = *device_derivative_history;
    const auto dt_view = device_dt.view();
    const size_t num_points = dt_view.extent(0);
    const size_t num_components = dt_view.extent(1);
    ASSERT(deriv_history.extent(1) == num_points,
           "Derivative history point extent mismatch: history="
               << deriv_history.extent(1) << " dt=" << num_points);
    ASSERT(deriv_history.extent(2) == num_components,
           "Derivative history component extent mismatch: history="
               << deriv_history.extent(2) << " dt=" << num_components);
    std::array<size_t, 8> deriv_history_strides{};
    deriv_history.stride(deriv_history_strides.data());
    std::array<size_t, 8> dt_strides{};
    dt_view.stride(dt_strides.data());

    detail::record_time_stepper_data_impl(
        deriv_history.data(), deriv_history_strides[0],
        deriv_history_strides[1], deriv_history_strides[2], dt_view.data(),
        dt_strides[0], dt_strides[1], substep, num_points, num_components);
  }
};

}  // namespace Actions::Kokkos
