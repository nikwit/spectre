// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ComputeInternalBoundaryTermsBatched {
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedBoundaryScratch>;
  using argument_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedTopology,
                 ScalarWave::Batched::Tags::PackedBoundaryMetadata>;

  static void apply(
      const gsl::not_null<
          ScalarWave::Batched::Tags::PackedBoundaryScratch::type*>
          packed_boundary_scratch,
      const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
      const ScalarWave::Batched::Tags::PackedBoundaryMetadata::type&
          packed_boundary_metadata);
};

}  // namespace ScalarWave::Actions
