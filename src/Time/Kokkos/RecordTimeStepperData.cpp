// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/Kokkos/RecordTimeStepperData.hpp"

#include <cstddef>

#include "Utilities/Kokkos/KokkosCore.hpp"

namespace Actions::Kokkos::detail {

void record_time_stepper_data_impl(
    double* deriv_history_data, const size_t deriv_history_stride_0,
    const size_t deriv_history_stride_1, const size_t deriv_history_stride_2,
    const double* dt_data, const size_t dt_stride_0, const size_t dt_stride_1,
    const size_t substep, const size_t num_points,
    const size_t num_components) {
  if (num_points == 0 or num_components == 0) {
    return;
  }
  const size_t substep_offset = substep * deriv_history_stride_0;
  ::Kokkos::parallel_for(
      "KokkosRecordTimeStepperData",
      ::Kokkos::RangePolicy<size_t>{0, num_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        const size_t point = linear_index / num_components;
        const size_t component = linear_index % num_components;
        deriv_history_data[substep_offset + point * deriv_history_stride_1 +
                           component * deriv_history_stride_2] =
            dt_data[point * dt_stride_0 + component * dt_stride_1];
      });
}

}  // namespace Actions::Kokkos::detail
