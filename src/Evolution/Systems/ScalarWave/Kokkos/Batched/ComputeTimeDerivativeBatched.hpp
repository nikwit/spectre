// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ComputeTimeDerivativeBatched {
 private:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using volume_time_derivative_terms =
      typename system::compute_volume_time_derivative_terms;
  using device_gradient_tags =
      db::wrap_tags_in<::Tags::MirrorView, typename system::gradient_variables>;
  using device_derivative_tags =
      db::wrap_tags_in<::Tags::deriv, device_gradient_tags,
                       tmpl::size_t<volume_dim>, Frame::Inertial>;

 public:
  static void compute_time_derivative_batched_volume_impl(
      const gsl::not_null<
          ScalarWave::Batched::Tags::PackedEvolutionState::type*>
          packed_evolution_state,
      const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
      const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry);

  using return_tags =
      tmpl::list<ScalarWave::Batched::Tags::PackedEvolutionState>;
  using argument_tags = tmpl::list<ScalarWave::Batched::Tags::PackedTopology,
                                   ScalarWave::Batched::Tags::PackedGeometry>;

  static void apply(
      const gsl::not_null<
          ScalarWave::Batched::Tags::PackedEvolutionState::type*>
          packed_evolution_state,
      const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
      const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry);
};

}  // namespace ScalarWave::Actions
