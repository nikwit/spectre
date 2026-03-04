// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/InitializeBatchedData.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <type_traits>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/TypeTraits.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "Time/History.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#endif

namespace gh::Actions {

void InitializeBatchedData::apply(
    const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const gsl::not_null<typename device_constraint_gamma0_tag::type*>
        device_constraint_gamma0,
    const gsl::not_null<typename device_constraint_gamma1_tag::type*>
        device_constraint_gamma1,
    const gsl::not_null<typename device_constraint_gamma2_tag::type*>
        device_constraint_gamma2,
    const typename packed_topology_tag::type& packed_topology,
    const Domain<volume_dim>& domain,
    const evolution::initial_data::InitialData& initial_data,
    const double initial_time,
    const typename damping_function_gamma0_tag::DampingFunctionType&
        damping_function_gamma0,
    const typename damping_function_gamma1_tag::DampingFunctionType&
        damping_function_gamma1,
    const typename damping_function_gamma2_tag::DampingFunctionType&
        damping_function_gamma2) {
#ifdef SPECTRE_KOKKOS
  const auto execution_space_name = ::Kokkos::DefaultExecutionSpace::name();
  const int execution_space_concurrency =
      ::Kokkos::DefaultExecutionSpace{}.concurrency();
  bool found_gpu = false;
  int cuda_device_count = 0;
#if defined(KOKKOS_ENABLE_CUDA)
  found_gpu = ::cudaGetDeviceCount(&cuda_device_count) == ::cudaSuccess and
              cuda_device_count > 0;
#endif
  Parallel::printf(
      "Kokkos startup (GH batched): initialized=%s, execution_space=%s, "
      "concurrency=%d, cuda_device_count=%d, found_gpu=%s\n",
      ::Kokkos::is_initialized() ? "true" : "false", execution_space_name,
      execution_space_concurrency, cuda_device_count,
      found_gpu ? "true" : "false");
#endif

  using host_variables_type = typename system::variables_tag::type;
  using host_dt_tag =
      db::add_tag_prefix<::Tags::dt, typename system::variables_tag>;
  using host_dt_variables_type = typename host_dt_tag::type;
  using device_variables_type =
      typename packed_evolution_state_tag::type::device_variables_type;
  using device_dt_variables_type =
      typename packed_evolution_state_tag::type::device_dt_variables_type;
  using device_step_start_type =
      typename packed_evolution_state_tag::type::device_step_start_type;
  using device_derivative_history_type =
      typename packed_evolution_state_tag::type::device_derivative_history_type;
  using device_constraint_gamma0_space =
      typename device_constraint_gamma0_tag::type::value_type::memory_space;
  using device_constraint_gamma1_space =
      typename device_constraint_gamma1_tag::type::value_type::memory_space;
  using device_constraint_gamma2_space =
      typename device_constraint_gamma2_tag::type::value_type::memory_space;
  using initial_data_classes =
      tmpl::append<gh::Solutions::all_solutions<volume_dim>,
                   tmpl::list<gh::NumericInitialData>>;

  packed_evolution_state->device_variables =
      device_variables_type(packed_topology.total_points);
  packed_evolution_state->device_dt_variables =
      device_dt_variables_type(packed_topology.total_points);
  packed_evolution_state->device_step_start =
      device_step_start_type(packed_topology.total_points);
  const size_t num_substeps = TimeSteppers::history_max_substeps;
  const size_t num_points = packed_topology.total_points;
  const size_t num_components =
      device_dt_variables_type::number_of_independent_components;
  const size_t point_stride = 1;
  const size_t component_stride = num_points > 0 ? num_points : 1;
  const size_t substep_stride =
      component_stride * (num_components > 0 ? num_components : 1);
  const ::Kokkos::LayoutStride derivative_history_layout(
      num_substeps, substep_stride, num_points, point_stride, num_components,
      component_stride);
  packed_evolution_state->device_derivative_history =
      device_derivative_history_type("GhBatchedDerivativeHistory",
                                     derivative_history_layout);

  *device_constraint_gamma0 = typename device_constraint_gamma0_tag::type(
      "GhBatchedConstraintGamma0", packed_topology.total_points);
  *device_constraint_gamma1 = typename device_constraint_gamma1_tag::type(
      "GhBatchedConstraintGamma1", packed_topology.total_points);
  *device_constraint_gamma2 = typename device_constraint_gamma2_tag::type(
      "GhBatchedConstraintGamma2", packed_topology.total_points);

  packed_geometry->element_inverse_jacobian_device =
      typename packed_geometry_tag::type::device_inverse_jacobian_type(
          "GhBatchedElementInverseJacobian",
          packed_topology.local_element_ids.size(),
          packed_topology.points_per_element, 9);
  packed_geometry->inertial_coordinates_host = {
      {DataVector{packed_topology.total_points},
       DataVector{packed_topology.total_points},
       DataVector{packed_topology.total_points}}};

  host_variables_type host_vars{packed_topology.total_points};
  host_dt_variables_type host_dt_vars{packed_topology.total_points};
  std::array<DataVector, volume_dim> host_inertial_coordinates{
      {DataVector{packed_topology.total_points},
       DataVector{packed_topology.total_points},
       DataVector{packed_topology.total_points}}};
  tnsr::I<DataVector, volume_dim, Frame::Grid> host_grid_coordinates{
      packed_topology.total_points};
  auto host_inverse_jacobian = ::Kokkos::create_mirror_view(
      packed_geometry->element_inverse_jacobian_device);
  if (host_dt_vars.size() > 0) {
    std::fill(host_dt_vars.data(), host_dt_vars.data() + host_dt_vars.size(),
              0.0);
  }

  call_with_dynamic_type<void, initial_data_classes>(
      &initial_data, [&domain, initial_time, &host_vars, &host_inverse_jacobian,
                      &host_inertial_coordinates, &host_grid_coordinates,
                      &packed_topology](const auto* const data_or_solution) {
        using initial_data_subclass = std::decay_t<decltype(*data_or_solution)>;
        if constexpr (is_analytic_data_v<initial_data_subclass> or
                      is_analytic_solution_v<initial_data_subclass>) {
          const size_t total_points = packed_topology.total_points;
          const auto& element_ids = packed_topology.local_element_ids;
          for (size_t e = 0; e < element_ids.size(); ++e) {
            const auto& element_id = element_ids[e];
            const auto& block = domain.blocks()[element_id.block_id()];
            ASSERT(not block.is_time_dependent(),
                   "EvolveGhNoBlackHoleKokkosBatched initial-data packing "
                   "currently supports only time-independent blocks.");

            const Mesh<volume_dim> mesh{
                gsl::at(packed_topology.element_extents_host, e),
                packed_topology.uniform_basis_host,
                packed_topology.uniform_quadrature_host};
            const auto logical_coords = logical_coordinates(mesh);
            const auto element_map =
                ElementMap<volume_dim, Frame::Grid>{element_id, block};
            const auto grid_coords = element_map(logical_coords);
            const auto inverse_jacobian =
                element_map.inv_jacobian(logical_coords);
            tnsr::I<DataVector, volume_dim, Frame::Inertial> inertial_coords{
                mesh.number_of_grid_points()};
            for (size_t d = 0; d < volume_dim; ++d) {
              inertial_coords.get(d) = grid_coords.get(d);
            }

            host_variables_type element_vars{mesh.number_of_grid_points()};
            element_vars.assign_subset(evolution::Initialization::initial_data(
                *data_or_solution, inertial_coords, initial_time,
                typename host_variables_type::tags_list{}));

            ASSERT(e + 1 < packed_topology.element_point_offsets_host.size(),
                   "Element point-offset metadata out of range for packed "
                   "element "
                       << e << ".");
            const size_t point_begin =
                packed_topology.element_point_offsets_host[e];
            const size_t point_end =
                packed_topology.element_point_offsets_host[e + 1];
            const size_t element_points = point_end - point_begin;
            ASSERT(element_points == mesh.number_of_grid_points(),
                   "Packed element point range mismatch for "
                       << element_id << ": packed points=" << element_points
                       << " mesh points=" << mesh.number_of_grid_points()
                       << ".");
            ASSERT(element_points == packed_topology.points_per_element,
                   "EvolveGhNoBlackHoleKokkosBatched currently requires "
                   "uniform points per element.");

            for (size_t component = 0;
                 component <
                 host_variables_type::number_of_independent_components;
                 ++component) {
              for (size_t p = 0; p < element_points; ++p) {
                host_vars.data()[component * total_points + point_begin + p] =
                    element_vars.data()[component * element_points + p];
              }
            }
            for (size_t d = 0; d < volume_dim; ++d) {
              for (size_t p = 0; p < element_points; ++p) {
                host_inertial_coordinates[d][point_begin + p] =
                    inertial_coords.get(d)[p];
                host_grid_coordinates.get(d)[point_begin + p] =
                    grid_coords.get(d)[p];
              }
            }
            for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
              for (size_t inertial_d = 0; inertial_d < volume_dim;
                   ++inertial_d) {
                for (size_t p = 0; p < element_points; ++p) {
                  host_inverse_jacobian(e, p,
                                        logical_d * volume_dim + inertial_d) =
                      inverse_jacobian.get(logical_d, inertial_d)[p];
                }
              }
            }
          }
        } else {
          ERROR(
              "Trying to use EvolveGhNoBlackHoleKokkosBatched with an "
              "initial data class that is not marked as analytic solution "
              "or analytic data.");
        }
      });

  Scalar<DataVector> host_gamma0{packed_topology.total_points};
  Scalar<DataVector> host_gamma1{packed_topology.total_points};
  Scalar<DataVector> host_gamma2{packed_topology.total_points};
  const domain::FunctionsOfTimeMap functions_of_time{};
  damping_function_gamma0(make_not_null(&host_gamma0), host_grid_coordinates,
                          initial_time, functions_of_time);
  damping_function_gamma1(make_not_null(&host_gamma1), host_grid_coordinates,
                          initial_time, functions_of_time);
  damping_function_gamma2(make_not_null(&host_gamma2), host_grid_coordinates,
                          initial_time, functions_of_time);

  for (size_t d = 0; d < volume_dim; ++d) {
    packed_geometry->inertial_coordinates_device.get(d) =
        ::Kokkos::View<double*>("GhBatchedInertialCoordinates",
                                packed_topology.total_points);
    auto host_coords_component = ::Kokkos::create_mirror_view(
        packed_geometry->inertial_coordinates_device.get(d));
    for (size_t s = 0; s < packed_topology.total_points; ++s) {
      host_coords_component(s) = host_inertial_coordinates[d][s];
    }
    ::Kokkos::deep_copy(packed_geometry->inertial_coordinates_device.get(d),
                        host_coords_component);
  }

  packed_geometry->inertial_coordinates_host =
      std::move(host_inertial_coordinates);
  packed_evolution_state->device_variables = copy_to_device(host_vars);
  packed_evolution_state->device_step_start = copy_to_device(host_vars);
  packed_evolution_state->device_dt_variables = copy_to_device(host_dt_vars);
  ::Kokkos::deep_copy(packed_geometry->element_inverse_jacobian_device,
                      host_inverse_jacobian);
  *device_constraint_gamma0 = copy_to_device(
      host_gamma0, tmpl::type_<device_constraint_gamma0_space>{});
  *device_constraint_gamma1 = copy_to_device(
      host_gamma1, tmpl::type_<device_constraint_gamma1_space>{});
  *device_constraint_gamma2 = copy_to_device(
      host_gamma2, tmpl::type_<device_constraint_gamma2_space>{});

  if (packed_topology.total_points > 0) {
    ::Kokkos::deep_copy(packed_evolution_state->device_derivative_history, 0.0);
  }
}

}  // namespace gh::Actions
