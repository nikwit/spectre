// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
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
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;

 public:
  static void compute_time_derivative_batched_volume_impl(
      const gsl::not_null<evolution::Kokkos::Tags::PackedEvolutionState<
          ScalarWave::System<3>>::type*>
          packed_evolution_state,
      const typename device_constraint_gamma2_tag::type&
          device_constraint_gamma2,
      const evolution::Kokkos::Tags::PackedTopology<
          ScalarWave::System<3>>::type& packed_topology,
      const evolution::Kokkos::Tags::PackedGeometry<
          ScalarWave::System<3>>::type& packed_geometry);

  using return_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedEvolutionState<ScalarWave::System<3>>>;
  using argument_tags = tmpl::list<
      device_constraint_gamma2_tag,
      evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>,
      evolution::Kokkos::Tags::PackedGeometry<ScalarWave::System<3>>>;

  static void apply(
      const gsl::not_null<evolution::Kokkos::Tags::PackedEvolutionState<
          ScalarWave::System<3>>::type*>
          packed_evolution_state,
      const typename device_constraint_gamma2_tag::type&
          device_constraint_gamma2,
      const evolution::Kokkos::Tags::PackedTopology<
          ScalarWave::System<3>>::type& packed_topology,
      const evolution::Kokkos::Tags::PackedGeometry<
          ScalarWave::System<3>>::type& packed_geometry);
};

}  // namespace ScalarWave::Actions
