// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/InitialExtents.hpp"
#include "Domain/Creators/Tags/InitialRefinementLevels.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/CreateInitialMesh.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/TypeTraits.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

template <typename Metavariables>
struct InitializeBatchedData {
 private:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using variables_tag = typename system::variables_tag;
  using host_variables_type = typename variables_tag::type;
  using host_dt_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
  using host_dt_variables_type = typename host_dt_tag::type;

 public:
  using return_tags = tmpl::list<ScalarWave::Batched::Tags::DeviceData>;
  using argument_tags =
      tmpl::list<::domain::Tags::Domain<volume_dim>,
                 evolution::initial_data::Tags::InitialData,
                 domain::Tags::InitialRefinementLevels<volume_dim>,
                 domain::Tags::InitialExtents<volume_dim>,
                 evolution::dg::Tags::Quadrature, ::Tags::Time>;
  using const_global_cache_tags =
      tmpl::list<::domain::Tags::Domain<volume_dim>,
                 evolution::initial_data::Tags::InitialData>;

  static void apply(
      const gsl::not_null<ScalarWave::Batched::Tags::DeviceData::type*>
          device_data,
      const Domain<volume_dim>& domain,
      const evolution::initial_data::InitialData& initial_data,
      const std::vector<std::array<size_t, volume_dim>>&
          initial_refinement_levels,
      const std::vector<std::array<size_t, volume_dim>>& initial_extents,
      const Spectral::Quadrature& quadrature, const double initial_time) {
    host_variables_type host_vars{device_data->total_points()};
    host_dt_variables_type host_dt_vars{device_data->total_points()};
    std::array<DataVector, volume_dim> host_inertial_coordinates{
        {DataVector{device_data->total_points()},
         DataVector{device_data->total_points()},
         DataVector{device_data->total_points()}}};
    auto host_inverse_jacobian = ::Kokkos::create_mirror_view(
        device_data->element_inverse_jacobian_device());
    if (host_dt_vars.size() > 0) {
      std::fill(host_dt_vars.data(), host_dt_vars.data() + host_dt_vars.size(),
                0.0);
    }
    std::array<Spectral::Basis, volume_dim> uniform_basis{};
    std::array<Spectral::Quadrature, volume_dim> uniform_quadrature{};
    bool have_uniform_mesh_metadata = false;

    using derived_classes =
        tmpl::at<typename Metavariables::factory_creation::factory_classes,
                 evolution::initial_data::InitialData>;
    call_with_dynamic_type<void, derived_classes>(
        &initial_data,
        [&domain, &initial_refinement_levels, &initial_extents, quadrature,
         initial_time, &host_vars, &host_inverse_jacobian, &uniform_basis,
         &uniform_quadrature, &have_uniform_mesh_metadata,
         &host_inertial_coordinates,
         total_points = device_data->total_points(),
         &element_ids = device_data->local_element_ids(),
         &device_data](const auto* const data_or_solution) {
          using initial_data_subclass =
              std::decay_t<decltype(*data_or_solution)>;
          if constexpr (is_analytic_data_v<initial_data_subclass> or
                        is_analytic_solution_v<initial_data_subclass>) {
            for (size_t e = 0; e < element_ids.size(); ++e) {
              const auto& element_id = element_ids[e];
              const auto& block = domain.blocks()[element_id.block_id()];
              ASSERT(not block.is_time_dependent(),
                     "EvolveScalarWaveKokkosBatched initial-data packing "
                     "currently supports only time-independent blocks.");

              const auto element =
                  ::domain::Initialization::create_initial_element(
                      element_id, domain.blocks(), initial_refinement_levels);
              const auto mesh = ::domain::Initialization::create_initial_mesh(
                  initial_extents, element, quadrature);
              if (not have_uniform_mesh_metadata) {
                for (size_t d = 0; d < volume_dim; ++d) {
                  uniform_basis[d] = mesh.basis(d);
                  uniform_quadrature[d] = mesh.quadrature(d);
                }
                have_uniform_mesh_metadata = true;
              } else {
                for (size_t d = 0; d < volume_dim; ++d) {
                  ASSERT(uniform_basis[d] == mesh.basis(d),
                         "EvolveScalarWaveKokkosBatched currently requires "
                         "uniform basis across all packed elements.");
                  ASSERT(uniform_quadrature[d] == mesh.quadrature(d),
                         "EvolveScalarWaveKokkosBatched currently requires "
                         "uniform quadrature across all packed elements.");
                }
              }
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
              element_vars.assign_subset(
                  evolution::Initialization::initial_data(
                      *data_or_solution, inertial_coords, initial_time,
                      typename host_variables_type::tags_list{}));

              const auto [point_begin, point_end] =
                  device_data->element_point_range(e);
              const size_t element_points = point_end - point_begin;
              ASSERT(element_points == mesh.number_of_grid_points(),
                     "Packed element point range mismatch for "
                         << element_id << ": packed points=" << element_points
                         << " mesh points=" << mesh.number_of_grid_points()
                         << ".");
              ASSERT(element_points == device_data->points_per_element(),
                     "EvolveScalarWaveKokkosBatched currently requires "
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
                "Trying to use EvolveScalarWaveKokkosBatched with an initial "
                "data class that is not marked as analytic solution or "
                "analytic data.");
          }
        });

    ASSERT(have_uniform_mesh_metadata,
           "Failed to gather uniform mesh metadata for packed elements.");
    device_data->set_uniform_mesh_metadata(device_data->uniform_extents_host(),
                                           uniform_basis, uniform_quadrature);
    device_data->set_inertial_coordinates_host(
        std::move(host_inertial_coordinates));
    device_data->device_variables() = copy_to_device(host_vars);
    device_data->device_step_start() = copy_to_device(host_vars);
    device_data->device_dt_variables() = copy_to_device(host_dt_vars);
    ::Kokkos::deep_copy(device_data->element_inverse_jacobian_device(),
                        host_inverse_jacobian);
    if (device_data->total_points() > 0) {
      ::Kokkos::deep_copy(device_data->device_derivative_history(), 0.0);
    }
  }
};

}  // namespace ScalarWave::Actions
