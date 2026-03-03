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
  // using Policy = ::Kokkos::MDRangePolicy<::Kokkos::Rank<2>>;
  // Policy pol({0, 0}, {num_components, num_points}, {2, 256});
  ::Kokkos::parallel_for(
      "KokkosUpdateUFused",
      ::Kokkos::RangePolicy<::Kokkos::Cuda, ::Kokkos::IndexType<int>>(
          0, num_points),
      KOKKOS_LAMBDA(const int point) {
#pragma unroll
        for (int i = 0; i < num_components; ++i) {
          double weighted_sum = 0.0;
          for (int coeff_index = 0; coeff_index < num_coefficients;
               ++coeff_index) {
            weighted_sum += coefficients_array[coeff_index] *
                            deriv_history_view(coeff_index, point, i);
          }
          u_view(point, i) = u0_view(point, i) + dt * weighted_sum;
        }
      });
}

}  // namespace Actions::Kokkos::detail
