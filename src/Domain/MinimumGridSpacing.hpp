// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"  // IWYU pragma: keep
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Index;
template <size_t Dim>
class Mesh;
/// \endcond

/// \ingroup ComputationalDomainGroup
/// Finds the minimum coordinate distance between grid points.
template <size_t Dim, typename Frame>
double minimum_grid_spacing(
    const Index<Dim>& extents,
    const tnsr::I<DataVector, Dim, Frame>& coords) noexcept;

namespace domain {
namespace Tags {
// @{
/// \ingroup ComputationalDomainGroup
/// \ingroup DataBoxTagsGroup
/// The minimum coordinate distance between grid points.
template <size_t Dim, typename Frame>
struct MinimumGridSpacing : db::SimpleTag {
  using type = double;
};

template <size_t Dim, typename Frame>
struct MinimumGridSpacingCompute : MinimumGridSpacing<Dim, Frame>,
                                   db::ComputeTag {
  using base = MinimumGridSpacing<Dim, Frame>;
  using return_type = double;
  static void function(
      const gsl::not_null<double*> result, const ::Mesh<Dim>& mesh,
      const tnsr::I<DataVector, Dim, Frame>& coordinates) noexcept {
    static double call_once = [&]() {
      return minimum_grid_spacing(mesh.extents(), coordinates);
    }();
    *result = call_once;
  }
  using argument_tags = tmpl::list<Mesh<Dim>, Coordinates<Dim, Frame>>;
};
// @}
}  // namespace Tags
}  // namespace domain
