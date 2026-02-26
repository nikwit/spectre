// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

struct ComputeInternalBoundaryTermsBatched {
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_boundary_scratch_tag =
      evolution::Kokkos::Tags::PackedBoundaryScratch<system>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;
  using packed_boundary_metadata_tag =
      evolution::Kokkos::Tags::PackedBoundaryMetadata<system>;

  using return_tags = tmpl::list<packed_boundary_scratch_tag>;
  using argument_tags =
      tmpl::list<packed_topology_tag, packed_boundary_metadata_tag>;

  static void apply(
      const gsl::not_null<typename packed_boundary_scratch_tag::type*>
          packed_boundary_scratch,
      const typename packed_topology_tag::type& packed_topology,
      const typename packed_boundary_metadata_tag::type&
          packed_boundary_metadata);
};

}  // namespace gh::Actions
