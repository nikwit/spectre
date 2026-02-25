// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched {
  static constexpr size_t volume_dim = 3;
  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedTopology,
                 ScalarWave::Batched::Tags::PackedGeometry, ::Tags::Time,
                 domain::Tags::ExternalBoundaryConditions<volume_dim>>;

  static void apply(
      const gsl::not_null<
          ScalarWave::Batched::Tags::PackedEvolutionState::type*>
          packed_evolution_state,
      const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
      const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry,
      const double time,
      const typename domain::Tags::ExternalBoundaryConditions<volume_dim>::type&
          external_boundary_conditions_by_block);
};

}  // namespace ScalarWave::Actions
