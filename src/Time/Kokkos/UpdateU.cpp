// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/Kokkos/UpdateU.hpp"

#include <cstddef>

#include "Utilities/Kokkos/KokkosCore.hpp"

namespace Actions::Kokkos::detail {

void update_u_impl(
    double* u_data, const size_t u_stride_0, const size_t u_stride_1,
    const double* u0_data, const size_t u0_stride_0, const size_t u0_stride_1,
    const double* deriv_history_data, const size_t deriv_history_stride_0,
    const size_t deriv_history_stride_1, const size_t deriv_history_stride_2,
    const double* coefficients, const size_t num_coefficients, const double dt,
    const size_t num_points, const size_t num_components) {
  using memory_space = typename ::Kokkos::DefaultExecutionSpace::memory_space;
  using unmanaged_view_2d =
      ::Kokkos::View<double**, ::Kokkos::LayoutStride, memory_space,
                     ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;
  using unmanaged_const_view_2d =
      ::Kokkos::View<const double**, ::Kokkos::LayoutStride, memory_space,
                     ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;
  using unmanaged_const_view_3d =
      ::Kokkos::View<const double***, ::Kokkos::LayoutStride, memory_space,
                     ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>;

  unmanaged_view_2d u_view(
      u_data, ::Kokkos::LayoutStride(num_points, u_stride_0, num_components,
                                     u_stride_1));
  unmanaged_const_view_2d u0_view(
      u0_data, ::Kokkos::LayoutStride(num_points, u0_stride_0, num_components,
                                      u0_stride_1));
  unmanaged_const_view_3d deriv_history_view(
      deriv_history_data,
      ::Kokkos::LayoutStride(num_coefficients, deriv_history_stride_0,
                             num_points, deriv_history_stride_1, num_components,
                             deriv_history_stride_2));

  ::Kokkos::Array<double, 8> coefficients_array{};
  for (size_t coeff_index = 0; coeff_index < num_coefficients; ++coeff_index) {
    coefficients_array[coeff_index] = coefficients[coeff_index];
  }

  ::Kokkos::parallel_for(
      "KokkosUpdateUFused",
      ::Kokkos::RangePolicy<size_t>{0, num_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        const size_t point = linear_index / num_components;
        const size_t component = linear_index % num_components;
        double weighted_sum = 0.0;
        for (size_t coeff_index = 0; coeff_index < num_coefficients;
             ++coeff_index) {
          weighted_sum += coefficients_array[coeff_index] *
                          deriv_history_view(coeff_index, point, component);
        }
        u_view(point, component) =
            u0_view(point, component) + dt * weighted_sum;
      });
}

}  // namespace Actions::Kokkos::detail
