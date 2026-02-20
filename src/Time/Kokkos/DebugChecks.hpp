// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "Utilities/Kokkos/KokkosCore.hpp"

namespace Actions::Kokkos::detail {

size_t count_non_finite_2d(const double* data, size_t stride_0, size_t stride_1,
                           size_t extent_0, size_t extent_1);

void throw_if_non_finite_2d(const double* data, size_t stride_0,
                            size_t stride_1, size_t extent_0, size_t extent_1,
                            const char* label);

template <typename View2D>
size_t count_non_finite_2d(const View2D& view) {
  std::array<size_t, 8> strides{};
  view.stride(strides.data());
  return count_non_finite_2d(view.data(), strides[0], strides[1],
                             view.extent(0), view.extent(1));
}

template <typename View2D>
void throw_if_non_finite_2d(const View2D& view, const char* const label) {
  std::array<size_t, 8> strides{};
  view.stride(strides.data());
  throw_if_non_finite_2d(view.data(), strides[0], strides[1], view.extent(0),
                         view.extent(1), label);
}

}  // namespace Actions::Kokkos::detail
