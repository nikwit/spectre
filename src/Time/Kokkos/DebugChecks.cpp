// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Time/Kokkos/DebugChecks.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>

#include "Utilities/Kokkos/KokkosCore.hpp"

namespace Actions::Kokkos::detail {

size_t count_non_finite_2d(const double* const data, const size_t stride_0,
                           const size_t stride_1, const size_t extent_0,
                           const size_t extent_1) {
  if (extent_0 == 0 or extent_1 == 0) {
    return 0;
  }
  size_t num_non_finite = 0;
  ::Kokkos::parallel_reduce(
      "KokkosCountNonFinite2D",
      ::Kokkos::RangePolicy<size_t>{0, extent_0 * extent_1},
      KOKKOS_LAMBDA(const size_t linear_index, size_t& local_non_finite) {
        const size_t i = linear_index / extent_1;
        const size_t j = linear_index % extent_1;
        const double value = data[i * stride_0 + j * stride_1];
        // Detect both NaN and +/-Inf on device without relying on math library
        // overloads.
        if ((value - value) != 0.0) {
          ++local_non_finite;
        }
      },
      num_non_finite);
  return num_non_finite;
}

void throw_if_non_finite_2d(const double* const data, const size_t stride_0,
                            const size_t stride_1, const size_t extent_0,
                            const size_t extent_1, const char* const label) {
  const size_t num_non_finite =
      count_non_finite_2d(data, stride_0, stride_1, extent_0, extent_1);
  if (num_non_finite > 0) {
    throw std::runtime_error(std::string{"Detected "} +
                             std::to_string(num_non_finite) +
                             " non-finite entries in " + label + ".");
  }
}

}  // namespace Actions::Kokkos::detail
