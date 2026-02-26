// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

struct PackageLocalFacesBatched {
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_boundary_scratch_tag =
      evolution::Kokkos::Tags::PackedBoundaryScratch<system>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;
  using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<system>;
  using device_constraint_gamma0_tag = gh::KokkosTags::DeviceConstraintGamma0;
  using device_constraint_gamma1_tag = gh::KokkosTags::DeviceConstraintGamma1;
  using device_constraint_gamma2_tag = gh::KokkosTags::DeviceConstraintGamma2;

  using return_tags = tmpl::list<packed_boundary_scratch_tag>;
  using argument_tags =
      tmpl::list<packed_topology_tag, packed_geometry_tag,
                 packed_evolution_state_tag, device_constraint_gamma0_tag,
                 device_constraint_gamma1_tag, device_constraint_gamma2_tag>;

  static void apply(
      const gsl::not_null<typename packed_boundary_scratch_tag::type*>
          packed_boundary_scratch,
      const typename packed_topology_tag::type& packed_topology,
      const typename packed_geometry_tag::type& packed_geometry,
      const typename packed_evolution_state_tag::type& packed_evolution_state,
      const typename device_constraint_gamma0_tag::type& device_gamma0,
      const typename device_constraint_gamma1_tag::type& device_gamma1,
      const typename device_constraint_gamma2_tag::type& device_gamma2);
};

}  // namespace gh::Actions
