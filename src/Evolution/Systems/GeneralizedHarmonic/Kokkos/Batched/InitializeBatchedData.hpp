// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Time/History.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename Metavariables>
struct InitializeBatchedData {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;
  using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<system>;
  using device_constraint_gamma0_tag = gh::KokkosTags::DeviceConstraintGamma0;
  using device_constraint_gamma1_tag = gh::KokkosTags::DeviceConstraintGamma1;
  using device_constraint_gamma2_tag = gh::KokkosTags::DeviceConstraintGamma2;

 public:
  using return_tags =
      tmpl::list<packed_geometry_tag, packed_evolution_state_tag,
                 device_constraint_gamma0_tag, device_constraint_gamma1_tag,
                 device_constraint_gamma2_tag>;
  using argument_tags = tmpl::list<packed_topology_tag>;
  using const_global_cache_tags = tmpl::list<>;

  static void apply(
      const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const gsl::not_null<typename packed_evolution_state_tag::type*>
          packed_evolution_state,
      const gsl::not_null<typename device_constraint_gamma0_tag::type*>
          device_constraint_gamma0,
      const gsl::not_null<typename device_constraint_gamma1_tag::type*>
          device_constraint_gamma1,
      const gsl::not_null<typename device_constraint_gamma2_tag::type*>
          device_constraint_gamma2,
      const typename packed_topology_tag::type& packed_topology) {
    using device_variables_type =
        typename packed_evolution_state_tag::type::device_variables_type;
    using device_dt_variables_type =
        typename packed_evolution_state_tag::type::device_dt_variables_type;
    using device_step_start_type =
        typename packed_evolution_state_tag::type::device_step_start_type;
    using device_derivative_history_type = typename packed_evolution_state_tag::
        type::device_derivative_history_type;

    packed_evolution_state->device_variables =
        device_variables_type(packed_topology.total_points);
    packed_evolution_state->device_dt_variables =
        device_dt_variables_type(packed_topology.total_points);
    packed_evolution_state->device_step_start =
        device_step_start_type(packed_topology.total_points);
    packed_evolution_state->device_derivative_history =
        device_derivative_history_type(
            "GhBatchedDerivativeHistory", TimeSteppers::history_max_substeps,
            packed_topology.total_points,
            device_dt_variables_type::number_of_independent_components);

    *device_constraint_gamma0 = typename device_constraint_gamma0_tag::type(
        "GhBatchedConstraintGamma0", packed_topology.total_points);
    *device_constraint_gamma1 = typename device_constraint_gamma1_tag::type(
        "GhBatchedConstraintGamma1", packed_topology.total_points);
    *device_constraint_gamma2 = typename device_constraint_gamma2_tag::type(
        "GhBatchedConstraintGamma2", packed_topology.total_points);

    packed_geometry->element_inverse_jacobian_device =
        ::Kokkos::View<double***>("GhBatchedElementInverseJacobian",
                                  packed_topology.local_element_ids.size(),
                                  packed_topology.points_per_element, 9);
    packed_geometry->inertial_coordinates_host = {
        {DataVector{packed_topology.total_points},
         DataVector{packed_topology.total_points},
         DataVector{packed_topology.total_points}}};

    if (packed_topology.total_points > 0) {
      ::Kokkos::deep_copy(packed_evolution_state->device_variables.view(), 0.0);
      ::Kokkos::deep_copy(packed_evolution_state->device_dt_variables.view(),
                          0.0);
      ::Kokkos::deep_copy(packed_evolution_state->device_step_start.view(),
                          0.0);
      ::Kokkos::deep_copy(packed_evolution_state->device_derivative_history,
                          0.0);
      ::Kokkos::deep_copy(get(*device_constraint_gamma0), 0.0);
      ::Kokkos::deep_copy(get(*device_constraint_gamma1), 0.0);
      ::Kokkos::deep_copy(get(*device_constraint_gamma2), 0.0);
      ::Kokkos::deep_copy(packed_geometry->element_inverse_jacobian_device,
                          0.0);
    }
  }
};

}  // namespace gh::Actions
