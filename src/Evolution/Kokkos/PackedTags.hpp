// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"

namespace evolution::Kokkos::Tags {

template <typename System>
struct PackedTopology : db::SimpleTag {
  using type = evolution::Kokkos::PackedTopology<System>;
};

template <typename System>
struct PackedGeometry : db::SimpleTag {
  using type = evolution::Kokkos::PackedGeometry<System>;
};

template <typename System>
struct PackedBoundaryMetadata : db::SimpleTag {
  using type = evolution::Kokkos::PackedBoundaryMetadata<System>;
};

template <typename System>
struct PackedEvolutionState : db::SimpleTag {
  using type = evolution::Kokkos::PackedEvolutionState<System>;
};

template <typename System>
struct PackedBoundaryScratch : db::SimpleTag {
  using type = evolution::Kokkos::PackedBoundaryScratch<System>;
};

}  // namespace evolution::Kokkos::Tags
