// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Executables/ScalarWave/Batched/PackedDataBundles.hpp"

namespace ScalarWave::Batched::Tags {

struct PackedTopology : db::SimpleTag {
  using type = ScalarWave::Batched::PackedTopology;
};

struct PackedGeometry : db::SimpleTag {
  using type = ScalarWave::Batched::PackedGeometry;
};

struct PackedBoundaryMetadata : db::SimpleTag {
  using type = ScalarWave::Batched::PackedBoundaryMetadata;
};

struct PackedEvolutionState : db::SimpleTag {
  using type = ScalarWave::Batched::PackedEvolutionState;
};

struct PackedBoundaryScratch : db::SimpleTag {
  using type = ScalarWave::Batched::PackedBoundaryScratch;
};

}  // namespace ScalarWave::Batched::Tags
