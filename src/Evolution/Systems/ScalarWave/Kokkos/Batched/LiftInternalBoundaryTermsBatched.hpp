// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct LiftInternalBoundaryTermsBatched {
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedTopology,
                 ScalarWave::Batched::Tags::PackedGeometry,
                 ScalarWave::Batched::Tags::PackedBoundaryScratch>;

  static void apply(
      const gsl::not_null<
          ScalarWave::Batched::Tags::PackedEvolutionState::type*>
          packed_evolution_state,
      const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
      const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry,
      const ScalarWave::Batched::Tags::PackedBoundaryScratch::type&
          packed_boundary_scratch);
};

}  // namespace ScalarWave::Actions
