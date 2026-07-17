// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/InitializeBatchedData.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <variant>

#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/SetPiAndPhiFromConstraints.tpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateConstraintDampingBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/TypeTraits.hpp"
#include "IO/Importers/Actions/ReadVolumeData.hpp"
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

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache>
void InitializeBatchedDataWithTag<FunctionsOfTimeTag,
                                  AddFunctionsOfTimeToConstGlobalCache>::
    apply(const gsl::not_null<typename packed_geometry_tag::type*>
              packed_geometry,
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
          const domain::FunctionsOfTimeMap& functions_of_time,
          const evolution::initial_data::InitialData& initial_data,
          const double initial_time,
          const gh::gauges::GaugeCondition& gauge_condition,
          const typename damping_function_gamma0_tag::DampingFunctionType&
              damping_function_gamma0,
          const typename damping_function_gamma1_tag::DampingFunctionType&
              damping_function_gamma1,
          const typename damping_function_gamma2_tag::DampingFunctionType&
              damping_function_gamma2) {
  if (not gauge_condition.is_harmonic()) {
    ERROR("GH batched evolution currently supports only the Harmonic gauge.");
  }
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

  const size_t num_elements = packed_topology.local_element_ids.size();
  const size_t points_per_element = packed_topology.points_per_element;
  ASSERT(points_per_element > 0,
         "GH batched initialization requires positive points per element.");
  ASSERT(num_elements * points_per_element == packed_topology.total_points,
         "Packed topology metadata is inconsistent with total packed points.");
  packed_geometry->logical_to_grid_inverse_jacobian_host =
      typename packed_geometry_tag::type::host_inverse_jacobian_cache_type(
          "GhBatchedLogicalToGridInverseJacobian", num_elements,
          points_per_element, 9);
  packed_geometry->grid_coordinates_host = {
      {DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0}}};
  packed_geometry->inertial_coordinates_host = {
      {DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0}}};
  packed_geometry->mesh_velocity_host = {
      {DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0}}};

  host_variables_type host_vars{packed_topology.total_points};
  host_dt_variables_type host_dt_vars{packed_topology.total_points};
  std::array<DataVector, volume_dim> host_inertial_coordinates{
      {DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0},
       DataVector{packed_topology.total_points, 0.0}}};
  if (host_dt_vars.size() > 0) {
    std::fill(host_dt_vars.data(), host_dt_vars.data() + host_dt_vars.size(),
              0.0);
  }

  std::vector<Mesh<volume_dim>> element_meshes{};
  element_meshes.reserve(num_elements);
  std::vector<tnsr::I<DataVector, volume_dim, Frame::ElementLogical>>
      element_logical_coordinates{};
  element_logical_coordinates.reserve(num_elements);
  std::vector<InverseJacobian<DataVector, volume_dim, Frame::ElementLogical,
                              Frame::Inertial>>
      element_logical_to_inertial_inverse_jacobians{};
  element_logical_to_inertial_inverse_jacobians.reserve(num_elements);
  std::unordered_map<ElementId<volume_dim>,
                     std::pair<tnsr::I<DataVector, volume_dim, Frame::Inertial>,
                               Mesh<volume_dim>>>
      numeric_import_targets{};

  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& block = domain.blocks()[element_id.block_id()];
    const Mesh<volume_dim> mesh{
        gsl::at(packed_topology.element_extents_host, e),
        packed_topology.uniform_basis_host,
        packed_topology.uniform_quadrature_host};
    const auto logical_coords = logical_coordinates(mesh);
    const auto element_map =
        ElementMap<volume_dim, Frame::Grid>{element_id, block};
    const auto grid_coords = element_map(logical_coords);
    const auto logical_to_grid_inverse_jacobian =
        element_map.inv_jacobian(logical_coords);
    tnsr::I<DataVector, volume_dim, Frame::Inertial> inertial_coords{
        mesh.number_of_grid_points()};
    InverseJacobian<DataVector, volume_dim, Frame::ElementLogical,
                    Frame::Inertial>
        logical_to_inertial_inverse_jacobian{mesh.number_of_grid_points()};
    if (block.is_time_dependent()) {
      const auto map_data =
          block.moving_mesh_grid_to_inertial_map()
              .coords_frame_velocity_jacobians(grid_coords, initial_time,
                                               functions_of_time);
      inertial_coords = std::get<0>(map_data);
      const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
      for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
        for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
          auto& composed_inverse_jacobian =
              logical_to_inertial_inverse_jacobian.get(logical_d, inertial_d);
          composed_inverse_jacobian =
              logical_to_grid_inverse_jacobian.get(logical_d, 0) *
              grid_to_inertial_inverse_jacobian.get(0, inertial_d);
          for (size_t grid_d = 1; grid_d < volume_dim; ++grid_d) {
            composed_inverse_jacobian +=
                logical_to_grid_inverse_jacobian.get(logical_d, grid_d) *
                grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d);
          }
        }
      }
    } else {
      for (size_t d = 0; d < volume_dim; ++d) {
        inertial_coords.get(d) = grid_coords.get(d);
      }
      for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
        for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
          logical_to_inertial_inverse_jacobian.get(logical_d, inertial_d) =
              logical_to_grid_inverse_jacobian.get(logical_d, inertial_d);
        }
      }
    }

    ASSERT(e + 1 < packed_topology.element_point_offsets_host.size(),
           "Element point-offset metadata out of range for packed element "
               << e << ".");
    const size_t point_begin = packed_topology.element_point_offsets_host[e];
    const size_t point_end = packed_topology.element_point_offsets_host[e + 1];
    const size_t element_points = point_end - point_begin;
    ASSERT(element_points == mesh.number_of_grid_points(),
           "Packed element point range mismatch for "
               << element_id << ": packed points=" << element_points
               << " mesh points=" << mesh.number_of_grid_points() << ".");
    ASSERT(element_points == packed_topology.points_per_element,
           "EvolveGhNoBlackHoleKokkosBatched currently requires uniform points "
           "per element.");

    for (size_t d = 0; d < volume_dim; ++d) {
      for (size_t p = 0; p < element_points; ++p) {
        host_inertial_coordinates[d][point_begin + p] =
            inertial_coords.get(d)[p];
        packed_geometry->grid_coordinates_host[d][point_begin + p] =
            grid_coords.get(d)[p];
      }
    }
    for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
      for (size_t grid_d = 0; grid_d < volume_dim; ++grid_d) {
        for (size_t p = 0; p < element_points; ++p) {
          packed_geometry->logical_to_grid_inverse_jacobian_host(
              e, p, logical_d * volume_dim + grid_d) =
              logical_to_grid_inverse_jacobian.get(logical_d, grid_d)[p];
        }
      }
    }

    element_meshes.push_back(mesh);
    element_logical_coordinates.push_back(logical_coords);
    element_logical_to_inertial_inverse_jacobians.push_back(
        logical_to_inertial_inverse_jacobian);
    numeric_import_targets.insert_or_assign(
        element_id, std::make_pair(inertial_coords, mesh));
  }

  const auto pack_element_vars = [&host_vars, &packed_topology](
                                     const host_variables_type& element_vars,
                                     const size_t element_index) {
    const size_t total_points = packed_topology.total_points;
    const size_t point_begin =
        packed_topology.element_point_offsets_host[element_index];
    const size_t point_end =
        packed_topology.element_point_offsets_host[element_index + 1];
    const size_t element_points = point_end - point_begin;
    for (size_t component = 0;
         component < host_variables_type::number_of_independent_components;
         ++component) {
      for (size_t p = 0; p < element_points; ++p) {
        host_vars.data()[component * total_points + point_begin + p] =
            element_vars.data()[component * element_points + p];
      }
    }
  };

  call_with_dynamic_type<void, initial_data_classes>(
      &initial_data,
      [&domain, &functions_of_time, initial_time, &gauge_condition, &host_vars,
       &packed_topology, &element_meshes, &element_logical_coordinates,
       &element_logical_to_inertial_inverse_jacobians, &numeric_import_targets,
       &pack_element_vars](const auto* const data_or_solution) {
        using initial_data_subclass = std::decay_t<decltype(*data_or_solution)>;
        if constexpr (is_analytic_data_v<initial_data_subclass> or
                      is_analytic_solution_v<initial_data_subclass>) {
          const auto& element_ids = packed_topology.local_element_ids;
          for (size_t e = 0; e < element_ids.size(); ++e) {
            const auto& element_id = element_ids[e];
            const auto& mesh = gsl::at(element_meshes, e);
            const auto& inertial_coords =
                numeric_import_targets.at(element_id).first;

            host_variables_type element_vars{mesh.number_of_grid_points()};
            element_vars.assign_subset(evolution::Initialization::initial_data(
                *data_or_solution, inertial_coords, initial_time,
                typename host_variables_type::tags_list{}));
            pack_element_vars(element_vars, e);
          }
        } else if constexpr (std::is_same_v<initial_data_subclass,
                                            gh::NumericInitialData>) {
          tuples::tagged_tuple_from_typelist<db::wrap_tags_in<
              importers::Tags::Selected, gh::NumericInitialData::all_vars>>
              selected_fields{};
          data_or_solution->select_for_import(make_not_null(&selected_fields));
          auto numeric_data_by_element =
              importers::detail::read_volume_data_for_elements<
                  volume_dim, gh::NumericInitialData::all_vars>(
                  data_or_solution->importer_options(), numeric_import_targets,
                  selected_fields);

          const bool set_pi_and_phi_from_constraints =
              not std::holds_alternative<gh::NumericInitialData::GhVars>(
                  data_or_solution->selected_variables());
          const auto& element_ids = packed_topology.local_element_ids;
          for (size_t e = 0; e < element_ids.size(); ++e) {
            const auto& element_id = element_ids[e];
            auto numeric_data =
                std::move(numeric_data_by_element.at(element_id));
            const auto& mesh = gsl::at(element_meshes, e);
            const auto& inverse_jacobian =
                gsl::at(element_logical_to_inertial_inverse_jacobians, e);

            tnsr::aa<DataVector, volume_dim, Frame::Inertial> spacetime_metric{
                mesh.number_of_grid_points()};
            tnsr::aa<DataVector, volume_dim, Frame::Inertial> pi{
                mesh.number_of_grid_points()};
            tnsr::iaa<DataVector, volume_dim, Frame::Inertial> phi{
                mesh.number_of_grid_points()};
            data_or_solution->set_initial_data(
                make_not_null(&spacetime_metric), make_not_null(&pi),
                make_not_null(&phi), make_not_null(&numeric_data), mesh,
                inverse_jacobian);

            const auto& block = domain.blocks()[element_id.block_id()];
            const auto element_map =
                ElementMap<volume_dim, Frame::Grid>{element_id, block};
            const auto grid_to_inertial_map =
                block.is_time_dependent()
                    ? block.moving_mesh_grid_to_inertial_map().get_clone()
                    : domain::make_coordinate_map_base<Frame::Grid,
                                                       Frame::Inertial>(
                          domain::CoordinateMaps::Identity<volume_dim>{});
            gh::gauges::SetPiAndPhiFromConstraints<
                gh::Solutions::all_solutions<volume_dim>,
                volume_dim>::impl(make_not_null(&pi), make_not_null(&phi),
                                  initial_time, mesh, element_map,
                                  *grid_to_inertial_map, functions_of_time,
                                  gsl::at(element_logical_coordinates, e),
                                  spacetime_metric, gauge_condition,
                                  set_pi_and_phi_from_constraints);

            host_variables_type element_vars{mesh.number_of_grid_points()};
            get<gr::Tags::SpacetimeMetric<DataVector, volume_dim>>(
                element_vars) = std::move(spacetime_metric);
            get<gh::Tags::Pi<DataVector, volume_dim>>(element_vars) =
                std::move(pi);
            get<gh::Tags::Phi<DataVector, volume_dim>>(element_vars) =
                std::move(phi);
            pack_element_vars(element_vars, e);
          }
        } else {
          ERROR(
              "Trying to use EvolveGhNoBlackHoleKokkosBatched with an "
              "initial data class that is not marked as analytic solution "
              "or analytic data.");
        }
      });

  packed_geometry->inertial_coordinates_host =
      std::move(host_inertial_coordinates);
  packed_evolution_state->device_variables = copy_to_device(host_vars);
  packed_evolution_state->device_step_start = copy_to_device(host_vars);
  packed_evolution_state->device_dt_variables = copy_to_device(host_dt_vars);
  UpdateFrameFieldsBatchedWithTag<
      FunctionsOfTimeTag,
      AddFunctionsOfTimeToConstGlobalCache>::apply(packed_geometry,
                                                   packed_topology, domain,
                                                   functions_of_time,
                                                   initial_time);
  UpdateConstraintDampingBatchedWithTag<FunctionsOfTimeTag,
                                        AddFunctionsOfTimeToConstGlobalCache>::
      apply(device_constraint_gamma0, device_constraint_gamma1,
            device_constraint_gamma2, *packed_geometry, functions_of_time,
            initial_time, damping_function_gamma0, damping_function_gamma1,
            damping_function_gamma2);

  if (packed_topology.total_points > 0) {
    ::Kokkos::deep_copy(packed_evolution_state->device_derivative_history, 0.0);
  }
}

template struct InitializeBatchedDataWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize>;
template struct InitializeBatchedDataWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize, false>;
template struct InitializeBatchedDataWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize>;
template struct InitializeBatchedDataWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize, false>;

}  // namespace gh::Actions
