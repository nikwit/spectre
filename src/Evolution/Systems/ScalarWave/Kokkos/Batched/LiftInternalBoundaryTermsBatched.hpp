// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct LiftInternalBoundaryTermsBatched {
  using return_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedEvolutionState<ScalarWave::System<3>>>;
  using argument_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>,
      evolution::Kokkos::Tags::PackedGeometry<ScalarWave::System<3>>,
      evolution::Kokkos::Tags::PackedBoundaryScratch<ScalarWave::System<3>>>;

  static void apply(
      const gsl::not_null<evolution::Kokkos::Tags::PackedEvolutionState<
          ScalarWave::System<3>>::type*>
          packed_evolution_state,
      const evolution::Kokkos::Tags::PackedTopology<
          ScalarWave::System<3>>::type& packed_topology,
      const evolution::Kokkos::Tags::PackedGeometry<
          ScalarWave::System<3>>::type& packed_geometry,
      const evolution::Kokkos::Tags::PackedBoundaryScratch<
          ScalarWave::System<3>>::type& packed_boundary_scratch);
};

}  // namespace ScalarWave::Actions
