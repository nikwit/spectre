// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Metafunctions.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

/*!
 * \brief Get the `Tensor` at a specific grid point index.
 *
 * Works with a `Tensor` of `Kokkos::View`s. Can also be made to work with
 * `DataVector` if needed (use `operator[]` instead of `operator()`).
 */
template <
    typename TensorType, std::integral... Is,
    typename ValueType = typename TensorType::type::value_type,
    typename ResultType = TensorMetafunctions::swap_type<ValueType, TensorType>>
KOKKOS_FUNCTION ResultType make_at_index(const TensorType& tensor,
                                         const Is&... i) {
  ResultType result{};
  for (size_t component = 0; component < TensorType::size(); ++component) {
    result[component] = tensor[component](i...);
  }
  return result;
}

template <typename TensorType, size_t NumIndices>
struct AtIndexView {
  const TensorType* tensor = nullptr;
  std::array<size_t, NumIndices> indices{};

  KOKKOS_FUNCTION static constexpr size_t size() { return TensorType::size(); }

  template <typename... Ns>
  KOKKOS_FUNCTION decltype(auto) get(const Ns&... n) const {
    return (*this)[TensorType::get_storage_index(n...)];
  }

  KOKKOS_FUNCTION decltype(auto) operator[](const size_t component) const {
    return component_at(component, std::make_index_sequence<NumIndices>{});
  }

 private:
  template <size_t... I>
  KOKKOS_FUNCTION decltype(auto) component_at(
      const size_t component, std::index_sequence<I...> /*meta*/) const {
    return (*tensor)[component](indices[I]...);
  }
};

/// Get a lightweight view of a tensor at a specific grid point index.
template <typename TensorType, typename... Is>
KOKKOS_FUNCTION AtIndexView<TensorType, sizeof...(Is)> make_at_index_view(
    const TensorType& tensor, const Is&... i) {
  return {&tensor, {static_cast<size_t>(i)...}};
}

template <typename TensorType, size_t NumIndices>
KOKKOS_FUNCTION decltype(auto) get(
    const AtIndexView<TensorType, NumIndices>& tensor_at_index_view) {
  static_assert(TensorType::rank() == 0,
                "This overload is only for rank-0 tensors (scalars).");
  return tensor_at_index_view[0];
}

template <typename TensorType>
struct AtIndexView1D {
  const TensorType* tensor;
  int idx;  // 32-bit

  KOKKOS_FUNCTION static constexpr size_t size() { return TensorType::size(); }

  template <typename... Ns>
  KOKKOS_FUNCTION decltype(auto) get(const Ns&... n) const {
    return (*tensor)[TensorType::get_storage_index(n...)](idx);
  }

  KOKKOS_FUNCTION decltype(auto) operator[](const size_t component) const {
    return (*tensor)[component](idx);
  }
};

template <typename TensorType>
KOKKOS_FUNCTION AtIndexView1D<TensorType> make_at_index_view_1d(
    const TensorType& tensor, int i) {
  return {&tensor, i};
}

/*!
 * \brief Set the `Tensor` at a specific grid point index.
 *
 * Works with a `Tensor` of `Kokkos::View`s. Can also be made to work with
 * `DataVector` if needed (use `operator[]` instead of `operator()`).
 */
template <typename TensorType, std::integral... Is,
          typename ValueType = typename TensorType::type::value_type>
KOKKOS_FUNCTION void set_at_index(
    const gsl::not_null<TensorType*> tensor,
    const TensorMetafunctions::swap_type<ValueType, TensorType>& value,
    const Is&... i) {
  for (size_t component = 0; component < value.size(); ++component) {
    (*tensor)[component](i...) = value[component];
  }
}

namespace Tags {

/*!
 * \brief Tag representing a specific grid point index of a tensor.
 *
 * This tag replaces the `Tensor` data type with `double` so that it can be
 * used with pointwise operations, e.g. in a `Kokkos::parallel_for` kernel.
 */
template <typename Tag>
struct AtIndex : db::PrefixTag {
 private:
  using TensorType = typename Tag::type;
  using ValueType = typename TensorType::type::value_type;

 public:
  using tag = Tag;
  using type = TensorMetafunctions::swap_type<ValueType, TensorType>;
};

}  // namespace Tags
