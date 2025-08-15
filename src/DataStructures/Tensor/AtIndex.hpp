// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Metafunctions.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

/// Get the `Tensor` at a specific grid point index.
template <typename TensorType, typename... Is>
KOKKOS_FUNCTION TensorMetafunctions::swap_type<double, TensorType>
make_at_index(const TensorType& tensor, const Is&... i) {
  TensorMetafunctions::swap_type<double, TensorType> result{};
  for (size_t component = 0; component < TensorType::size(); ++component) {
    result[component] = tensor[component](i...);
  }
  return result;
}

/// Set the `Tensor` at a specific grid point index.
template <typename TensorType, typename... Is>
KOKKOS_FUNCTION void set_at_index(
    const gsl::not_null<TensorType*> tensor,
    const TensorMetafunctions::swap_type<double, TensorType>& value,
    const Is&... i) {
  for (size_t component = 0; component < TensorType::size(); ++component) {
    (*tensor)[component](i...) = value[component];
  }
}

namespace Tags {

/*!
 * \brief Tag representing a specific grid point index of a tensor
 *
 * This tag replaces the `Tensor` data type with `double` so that it can be
 * used with pointwise operations, e.g. in a `Kokkos::parallel_for` kernel.
 */
template <typename Tag>
struct AtIndex : db::PrefixTag {
  using tag = Tag;
  using type = TensorMetafunctions::swap_type<double, typename Tag::type>;
};

}  // namespace Tags
