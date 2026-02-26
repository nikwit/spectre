// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched {
  static constexpr size_t volume_dim = 3;
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;
  using return_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedEvolutionState<ScalarWave::System<3>>>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>,
                 evolution::Kokkos::Tags::PackedGeometry<ScalarWave::System<3>>,
                 device_constraint_gamma2_tag,
                 ::Tags::Time,
                 domain::Tags::ExternalBoundaryConditions<volume_dim>>;

  static void apply(
      const gsl::not_null<evolution::Kokkos::Tags::PackedEvolutionState<
          ScalarWave::System<3>>::type*>
          packed_evolution_state,
      const evolution::Kokkos::Tags::PackedTopology<
          ScalarWave::System<3>>::type& packed_topology,
      const evolution::Kokkos::Tags::PackedGeometry<
          ScalarWave::System<3>>::type& packed_geometry,
      const typename device_constraint_gamma2_tag::type&
          device_constraint_gamma2,
      const double time,
      const typename domain::Tags::ExternalBoundaryConditions<volume_dim>::type&
          external_boundary_conditions_by_block);
};

}  // namespace ScalarWave::Actions
