// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

#include "DataStructures/Tensor/Expressions/SpatialSpacetimeIndex.hpp"
#include "DataStructures/Tensor/Expressions/TensorExpression.hpp"
#include "DataStructures/Tensor/Expressions/TensorIndex.hpp"
#include "DataStructures/Tensor/Expressions/TimeIndex.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Symmetry.hpp"
#include "Utilities/ForceInline.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"

namespace tenex {
namespace detail {
template <typename IndexType1, typename IndexType2>
constexpr bool kronecker_delta_index_types_compatible() {
  return std::is_same_v<typename IndexType1::Frame,
                        typename IndexType2::Frame> and
         ((IndexType1::index_type == IndexType2::index_type and
           IndexType1::dim == IndexType2::dim) or
          (IndexType1::index_type == IndexType::Spacetime and
           IndexType1::dim == IndexType2::dim + 1) or
          (IndexType2::index_type == IndexType::Spacetime and
           IndexType1::dim + 1 == IndexType2::dim));
}

}  // namespace detail

template <typename IndexType1, typename IndexType2, typename ArgsList>
struct KroneckerDeltaAsExpression;

template <typename IndexType1, typename IndexType2,
          template <typename...> class ArgsList, typename Arg1, typename Arg2>
struct KroneckerDeltaAsExpression<IndexType1, IndexType2, ArgsList<Arg1, Arg2>>
    : public TensorExpression<KroneckerDeltaAsExpression<IndexType1, IndexType2,
                                                         ArgsList<Arg1, Arg2>>,
                              double, Symmetry<2, 1>,
                              tmpl::list<IndexType1, IndexType2>,
                              ArgsList<Arg1, Arg2>> {
  static_assert(tt::is_tensor_index<Arg1>::value and
                    tt::is_tensor_index<Arg2>::value,
                "Kronecker delta arguments must be TensorIndex types.");
  static_assert(Arg1::valence != Arg2::valence,
                "Kronecker delta indices must have opposite valence.");
  static_assert(
      Arg1::valence == IndexType1::ul and Arg2::valence == IndexType2::ul,
      "Kronecker delta generic index valences must match the valence of their "
      "associated index types.");
  static_assert(
      detail::kronecker_delta_index_types_compatible<IndexType1, IndexType2>(),
      "Kronecker delta index types must be compatible (same frame and same "
      "index type/dimension, or a valid spatial/spacetime pairing).");
  static_assert(
      not(Arg1::is_spacetime and
          IndexType1::index_type == IndexType::Spatial) and
          not(Arg2::is_spacetime and
              IndexType2::index_type == IndexType::Spatial),
      "Cannot use spacetime generic indices with spatial index types.");

  using type = double;
  using symmetry = Symmetry<2, 1>;
  using index_list = tmpl::list<IndexType1, IndexType2>;
  using args_list = ArgsList<Arg1, Arg2>;
  static constexpr auto num_tensor_indices = 2;

  static constexpr size_t num_ops_left_child = 0;
  static constexpr size_t num_ops_right_child = 0;
  static constexpr size_t num_ops_subtree = 0;
  static constexpr size_t height_relative_to_closest_tensor_leaf_in_subtree =
      std::numeric_limits<size_t>::max();

  static constexpr bool is_primary_end = true;
  static constexpr size_t num_ops_to_evaluate_primary_left_child = 0;
  static constexpr size_t num_ops_to_evaluate_primary_right_child = 0;
  static constexpr size_t num_ops_to_evaluate_primary_subtree = 0;
  static constexpr bool is_primary_start = false;
  static constexpr bool primary_child_subtree_contains_primary_start = false;
  static constexpr bool primary_subtree_contains_primary_start = false;

  KroneckerDeltaAsExpression() = default;
  ~KroneckerDeltaAsExpression() override = default;

  template <typename LhsTensor>
  SPECTRE_ALWAYS_INLINE void assert_lhs_tensor_not_in_rhs_expression(
      const gsl::not_null<LhsTensor*> /*lhs_tensor*/) const {}

  template <typename LhsTensorIndices, typename LhsTensor>
  SPECTRE_ALWAYS_INLINE void assert_lhs_tensorindices_same_in_rhs(
      const gsl::not_null<LhsTensor*> /*lhs_tensor*/) const {}

  size_t get_rhs_tensor_component_size() const = delete;

  SPECTRE_ALWAYS_INLINE type
  get(const std::array<size_t, num_tensor_indices>& multi_index) const {
    return gsl::at(multi_index, 0) == gsl::at(multi_index, 1) ? 1.0 : 0.0;
  }

  template <typename ResultType>
  SPECTRE_ALWAYS_INLINE type
  get_primary(const ResultType& /*result_component*/,
              const std::array<size_t, num_tensor_indices>& multi_index) const {
    return get(multi_index);
  }

  template <typename ResultType>
  void evaluate_primary_subtree(
      ResultType&,
      const std::array<size_t, num_tensor_indices>&) const = delete;
};

template <typename T, size_t Position, typename NewArg, typename NewIndexType,
          typename IndexList = typename T::index_list,
          typename ArgsList = typename T::args_list>
struct TensorIndexSubstitution;

template <typename T, size_t Position, typename NewArg, typename NewIndexType,
          template <typename...> class IndexList, typename... Indices,
          template <typename...> class ArgsList, typename... Args>
struct TensorIndexSubstitution<T, Position, NewArg, NewIndexType,
                               IndexList<Indices...>, ArgsList<Args...>>
    : public TensorExpression<
          TensorIndexSubstitution<T, Position, NewArg, NewIndexType,
                                  IndexList<Indices...>, ArgsList<Args...>>,
          typename T::type, typename T::symmetry,
          tmpl::replace_at<IndexList<Indices...>, tmpl::size_t<Position>,
                           NewIndexType>,
          tmpl::replace_at<ArgsList<Args...>, tmpl::size_t<Position>, NewArg>> {
  static_assert(Position < sizeof...(Args),
                "Index substitution out of bounds.");
  static_assert(tt::is_tensor_index<NewArg>::value,
                "Replacement generic index must be a TensorIndex.");

  using old_arg = tmpl::at_c<ArgsList<Args...>, Position>;
  using old_index_type = tmpl::at_c<IndexList<Indices...>, Position>;

  static_assert(tt::is_tensor_index<old_arg>::value,
                "Index substitution requires TensorIndex arguments.");
  static_assert(old_arg::valence == NewArg::valence,
                "Replacement generic index must have the same valence.");
  static_assert(NewArg::valence == NewIndexType::ul,
                "Replacement generic index valence must match replacement "
                "index type valence.");
  static_assert(detail::kronecker_delta_index_types_compatible<old_index_type,
                                                               NewIndexType>(),
                "Replacement index type must be compatible with the original "
                "index type.");

  using type = typename T::type;
  using symmetry = typename T::symmetry;
  using index_list = tmpl::replace_at<IndexList<Indices...>,
                                      tmpl::size_t<Position>, NewIndexType>;
  using args_list =
      tmpl::replace_at<ArgsList<Args...>, tmpl::size_t<Position>, NewArg>;
  static constexpr auto num_tensor_indices = tmpl::size<index_list>::value;
  static constexpr bool substitution_may_fail =
      NewIndexType::index_type == IndexType::Spacetime and
      (not NewArg::is_spacetime or
       old_index_type::index_type == IndexType::Spatial or
       (old_index_type::index_type == IndexType::Spacetime and
        not old_arg::is_spacetime));

  static constexpr size_t num_ops_left_child = T::num_ops_left_child;
  static constexpr size_t num_ops_right_child = T::num_ops_right_child;
  static constexpr size_t num_ops_subtree = T::num_ops_subtree;
  static constexpr size_t height_relative_to_closest_tensor_leaf_in_subtree =
      T::height_relative_to_closest_tensor_leaf_in_subtree !=
              std::numeric_limits<size_t>::max()
          ? T::height_relative_to_closest_tensor_leaf_in_subtree + 1
          : T::height_relative_to_closest_tensor_leaf_in_subtree;

  static constexpr bool is_primary_end = T::is_primary_end;
  static constexpr size_t num_ops_to_evaluate_primary_left_child =
      T::num_ops_to_evaluate_primary_left_child;
  static constexpr size_t num_ops_to_evaluate_primary_right_child =
      T::num_ops_to_evaluate_primary_right_child;
  static constexpr size_t num_ops_to_evaluate_primary_subtree =
      T::num_ops_to_evaluate_primary_subtree;
  static constexpr bool is_primary_start = T::is_primary_start;
  static constexpr bool primary_child_subtree_contains_primary_start =
      T::primary_child_subtree_contains_primary_start;
  static constexpr bool primary_subtree_contains_primary_start =
      T::primary_subtree_contains_primary_start;

  explicit TensorIndexSubstitution(T t) : t_(std::move(t)) {}
  ~TensorIndexSubstitution() override = default;

  template <typename LhsTensor>
  SPECTRE_ALWAYS_INLINE void assert_lhs_tensor_not_in_rhs_expression(
      const gsl::not_null<LhsTensor*> lhs_tensor) const {
    t_.assert_lhs_tensor_not_in_rhs_expression(lhs_tensor);
  }

  template <typename LhsTensorIndices, typename LhsTensor>
  SPECTRE_ALWAYS_INLINE void assert_lhs_tensorindices_same_in_rhs(
      const gsl::not_null<LhsTensor*> lhs_tensor) const {
    t_.template assert_lhs_tensorindices_same_in_rhs<LhsTensorIndices>(
        lhs_tensor);
  }

  SPECTRE_ALWAYS_INLINE size_t get_rhs_tensor_component_size() const {
    return t_.get_rhs_tensor_component_size();
  }

  SPECTRE_ALWAYS_INLINE decltype(auto) get(
      const std::array<size_t, num_tensor_indices>& multi_index) const {
    auto operand_multi_index = multi_index;
    if constexpr (substitution_may_fail) {
      if (not map_substituted_index(
              make_not_null(&gsl::at(operand_multi_index, Position)),
              gsl::at(multi_index, Position))) {
        auto reference_multi_index = operand_multi_index;
        gsl::at(reference_multi_index, Position) = 0;
        return make_with_value<type>(t_.get(reference_multi_index), 0.0);
      }
      return type{t_.get(operand_multi_index)};
    } else {
      (void)map_substituted_index(
          make_not_null(&gsl::at(operand_multi_index, Position)),
          gsl::at(multi_index, Position));
      return t_.get(operand_multi_index);
    }
  }

  template <typename ResultType>
  SPECTRE_ALWAYS_INLINE decltype(auto) get_primary(
      const ResultType& result_component,
      const std::array<size_t, num_tensor_indices>& multi_index) const {
    auto operand_multi_index = multi_index;
    if constexpr (substitution_may_fail) {
      if (not map_substituted_index(
              make_not_null(&gsl::at(operand_multi_index, Position)),
              gsl::at(multi_index, Position))) {
        auto reference_multi_index = operand_multi_index;
        gsl::at(reference_multi_index, Position) = 0;
        return make_with_value<type>(t_.get(reference_multi_index), 0.0);
      }
      return type{t_.get_primary(result_component, operand_multi_index)};
    } else {
      (void)map_substituted_index(
          make_not_null(&gsl::at(operand_multi_index, Position)),
          gsl::at(multi_index, Position));
      return t_.get_primary(result_component, operand_multi_index);
    }
  }

  template <typename ResultType>
  SPECTRE_ALWAYS_INLINE void evaluate_primary_subtree(
      ResultType& result_component,
      const std::array<size_t, num_tensor_indices>& multi_index) const {
    auto operand_multi_index = multi_index;
    if (map_substituted_index(
            make_not_null(&gsl::at(operand_multi_index, Position)),
            gsl::at(multi_index, Position))) {
      t_.evaluate_primary_subtree(result_component, operand_multi_index);
    }
  }

  SPECTRE_ALWAYS_INLINE const T& operand_expression() const { return t_; }

 private:
  SPECTRE_ALWAYS_INLINE static bool map_substituted_index(
      const gsl::not_null<size_t*> old_concrete_index,
      const size_t new_concrete_index) {
    size_t spacetime_component = 0;
    if constexpr (NewIndexType::index_type == IndexType::Spacetime) {
      if constexpr (not NewArg::is_spacetime) {
        if (new_concrete_index == 0) {
          return false;
        }
      }
      spacetime_component = new_concrete_index;
    } else {
      spacetime_component = new_concrete_index + 1;
    }

    if constexpr (old_index_type::index_type == IndexType::Spacetime) {
      if constexpr (not old_arg::is_spacetime) {
        if (spacetime_component == 0) {
          return false;
        }
      }
      *old_concrete_index = spacetime_component;
    } else {
      if (spacetime_component == 0) {
        return false;
      }
      *old_concrete_index = spacetime_component - 1;
    }
    return true;
  }

  T t_;
};

template <typename OldArg, typename NewArg, typename NewIndexType, typename T>
SPECTRE_ALWAYS_INLINE constexpr auto substitute_tensor_index(T t) {
  static_assert(tmpl::list_contains<typename T::args_list, OldArg>::value,
                "Cannot substitute a generic index that is not present in the "
                "expression.");
  return TensorIndexSubstitution<
      T, tmpl::index_of<typename T::args_list, OldArg>::value, NewArg,
      NewIndexType>{std::move(t)};
}

template <typename OldArg, typename NewArg, typename NewIndexType,
          typename IndexType1, typename IndexType2,
          template <typename...> class ArgsList, typename Arg1, typename Arg2>
SPECTRE_ALWAYS_INLINE constexpr auto substitute_tensor_index(
    KroneckerDeltaAsExpression<IndexType1, IndexType2, ArgsList<Arg1, Arg2>>
    /*delta*/) {
  static_assert(tmpl::list_contains<tmpl::list<Arg1, Arg2>, OldArg>::value,
                "Cannot substitute a generic index that is not present in the "
                "expression.");
  if constexpr (std::is_same_v<OldArg, Arg1>) {
    return KroneckerDeltaAsExpression<NewIndexType, IndexType2,
                                      tmpl::list<NewArg, Arg2>>{};
  } else {
    static_assert(std::is_same_v<OldArg, Arg2>,
                  "Cannot substitute a generic index that is not present in "
                  "the expression.");
    return KroneckerDeltaAsExpression<IndexType1, NewIndexType,
                                      tmpl::list<Arg1, NewArg>>{};
  }
}

template <typename IndexType1, typename IndexType2, typename Arg1,
          typename Arg2,
          Requires<tt::is_tensor_index<Arg1>::value and
                   tt::is_tensor_index<Arg2>::value> = nullptr>
SPECTRE_ALWAYS_INLINE constexpr auto kronecker_delta(Arg1 /*arg1*/,
                                                     Arg2 /*arg2*/) {
  return KroneckerDeltaAsExpression<IndexType1, IndexType2,
                                    tmpl::list<Arg1, Arg2>>{};
}

template <size_t Dim, typename Frame, typename Arg1, typename Arg2,
          Requires<tt::is_tensor_index<Arg1>::value and
                   tt::is_tensor_index<Arg2>::value> = nullptr>
SPECTRE_ALWAYS_INLINE constexpr auto spatial_kronecker_delta(Arg1 arg1,
                                                             Arg2 arg2) {
  using index_type1 = SpatialIndex<Dim, Arg1::valence, Frame>;
  using index_type2 = SpatialIndex<Dim, Arg2::valence, Frame>;
  return kronecker_delta<index_type1, index_type2>(arg1, arg2);
}

template <size_t Dim, typename Frame, typename Arg1, typename Arg2,
          Requires<tt::is_tensor_index<Arg1>::value and
                   tt::is_tensor_index<Arg2>::value> = nullptr>
SPECTRE_ALWAYS_INLINE constexpr auto spacetime_kronecker_delta(Arg1 arg1,
                                                               Arg2 arg2) {
  using index_type1 = SpacetimeIndex<Dim, Arg1::valence, Frame>;
  using index_type2 = SpacetimeIndex<Dim, Arg2::valence, Frame>;
  return kronecker_delta<index_type1, index_type2>(arg1, arg2);
}
}  // namespace tenex
