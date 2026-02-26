// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct PackageLocalFacesBatched {
  using system = ScalarWave::System<3>;
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;
  using return_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedBoundaryScratch<system>>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<system>,
                 evolution::Kokkos::Tags::PackedGeometry<system>,
                 device_constraint_gamma2_tag,
                 evolution::Kokkos::Tags::PackedEvolutionState<system>>;

  static void apply(
      const gsl::not_null<
          evolution::Kokkos::Tags::PackedBoundaryScratch<system>::type*>
          packed_boundary_scratch,
      const evolution::Kokkos::Tags::PackedTopology<system>::type&
          packed_topology,
      const evolution::Kokkos::Tags::PackedGeometry<system>::type&
          packed_geometry,
      const typename device_constraint_gamma2_tag::type&
          device_constraint_gamma2,
      const evolution::Kokkos::Tags::PackedEvolutionState<system>::type&
          packed_evolution_state);
};

}  // namespace ScalarWave::Actions
