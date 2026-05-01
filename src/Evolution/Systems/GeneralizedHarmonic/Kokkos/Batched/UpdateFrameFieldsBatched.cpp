// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"

#include <array>
#include <cstddef>
#include <tuple>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

template <size_t VolumeDim>
void resize_host_array_if_needed(
    const gsl::not_null<std::array<DataVector, VolumeDim>*> host_data,
    const size_t total_points) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    if (gsl::at(*host_data, d).size() != total_points) {
      gsl::at(*host_data, d) = DataVector{total_points, 0.0};
    }
  }
}

template <size_t VolumeDim>
void ensure_frame_field_storage(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t num_elements, const size_t points_per_element,
    const size_t total_points) {
  using packed_geometry_type =
      evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>;

  if (packed_geometry->element_inverse_jacobian_device.extent(0) !=
          num_elements or
      packed_geometry->element_inverse_jacobian_device.extent(1) !=
          points_per_element or
      packed_geometry->element_inverse_jacobian_device.extent(2) !=
          VolumeDim * VolumeDim) {
    packed_geometry->element_inverse_jacobian_device =
        typename packed_geometry_type::device_inverse_jacobian_type{
            "GhBatchedElementInverseJacobian", num_elements, points_per_element,
            VolumeDim * VolumeDim};
  }

  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->inertial_coordinates_host), total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->mesh_velocity_host), total_points);

  for (size_t d = 0; d < VolumeDim; ++d) {
    if (packed_geometry->inertial_coordinates_device.get(d).extent(0) !=
        total_points) {
      packed_geometry->inertial_coordinates_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedInertialCoordinates", total_points);
    }
    if (packed_geometry->mesh_velocity_device.get(d).extent(0) !=
        total_points) {
      packed_geometry->mesh_velocity_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedMeshVelocity", total_points);
    }
  }
}

template <size_t VolumeDim>
void copy_host_array_to_device(
    const std::array<DataVector, VolumeDim>& host_data,
    const gsl::not_null<typename evolution::Kokkos::PackedGeometry<
        gh::System<VolumeDim>>::device_inertial_coordinates_type*>
        device_data) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    auto host_component = ::Kokkos::create_mirror_view(device_data->get(d));
    for (size_t s = 0; s < host_data[d].size(); ++s) {
      host_component(s) = host_data[d][s];
    }
    ::Kokkos::deep_copy(device_data->get(d), host_component);
  }
}

}  // namespace

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache>
void UpdateFrameFieldsBatchedWithTag<FunctionsOfTimeTag,
                                     AddFunctionsOfTimeToConstGlobalCache>::
    apply(const gsl::not_null<typename packed_geometry_tag::type*>
              packed_geometry,
          const typename packed_topology_tag::type& packed_topology,
          const Domain<volume_dim>& domain,
          const domain::FunctionsOfTimeMap& functions_of_time,
          const double time) {
  const size_t total_points = packed_topology.total_points;
  if (total_points == 0) {
    return;
  }

  const size_t num_elements = packed_topology.local_element_ids.size();
  const size_t points_per_element = packed_topology.points_per_element;
  ASSERT(points_per_element > 0,
         "Packed points-per-element must be positive for GH batched updates.");
  ASSERT(num_elements * points_per_element == total_points,
         "Packed topology metadata is inconsistent with total packed points.");

  ASSERT(packed_geometry->logical_to_grid_inverse_jacobian_host.extent(0) ==
                 num_elements and
             packed_geometry->logical_to_grid_inverse_jacobian_host.extent(1) ==
                 points_per_element and
             packed_geometry->logical_to_grid_inverse_jacobian_host.extent(2) ==
                 volume_dim * volume_dim,
         "Packed logical-to-grid inverse Jacobian cache must be initialized by "
         "domain-static initialization before UpdateFrameFieldsBatched.");

  ensure_frame_field_storage<volume_dim>(packed_geometry, num_elements,
                                         points_per_element, total_points);

  for (size_t d = 0; d < volume_dim; ++d) {
    ASSERT(packed_geometry->grid_coordinates_host[d].size() == total_points,
           "Packed grid coordinates are not initialized. Expected "
               << total_points << " points in component " << d << ", got "
               << packed_geometry->grid_coordinates_host[d].size() << ".");
  }

  auto host_inverse_jacobian = ::Kokkos::create_mirror_view(
      packed_geometry->element_inverse_jacobian_device);

  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& block = domain.blocks()[element_id.block_id()];

    const size_t point_begin = packed_topology.element_point_offsets_host[e];
    const size_t point_end = packed_topology.element_point_offsets_host[e + 1];
    const size_t element_points = point_end - point_begin;
    ASSERT(element_points == points_per_element,
           "GH batched frame update assumes uniform points per element.");

    tnsr::I<DataVector, volume_dim, Frame::Grid> grid_coords{element_points};
    for (size_t d = 0; d < volume_dim; ++d) {
      for (size_t p = 0; p < element_points; ++p) {
        grid_coords.get(d)[p] =
            packed_geometry->grid_coordinates_host[d][point_begin + p];
      }
    }

    if (block.is_time_dependent()) {
      const auto map_data = block.moving_mesh_grid_to_inertial_map()
                                .coords_frame_velocity_jacobians(
                                    grid_coords, time, functions_of_time);
      const auto& inertial_coords = std::get<0>(map_data);
      const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
      const auto& mesh_velocity = std::get<3>(map_data);

      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              inertial_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] =
              mesh_velocity.get(d)[p];
        }

        for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
          for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
            double composed_inverse_jacobian = 0.0;
            for (size_t grid_d = 0; grid_d < volume_dim; ++grid_d) {
              composed_inverse_jacobian +=
                  packed_geometry->logical_to_grid_inverse_jacobian_host(
                      e, p, logical_d * volume_dim + grid_d) *
                  grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p];
            }
            host_inverse_jacobian(e, p, logical_d * volume_dim + inertial_d) =
                composed_inverse_jacobian;
          }
        }
      }
    } else {
      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              grid_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] = 0.0;
        }
        for (size_t component = 0; component < volume_dim * volume_dim;
             ++component) {
          host_inverse_jacobian(e, p, component) =
              packed_geometry->logical_to_grid_inverse_jacobian_host(e, p,
                                                                     component);
        }
      }
    }
  }

  copy_host_array_to_device<volume_dim>(
      packed_geometry->inertial_coordinates_host,
      make_not_null(&packed_geometry->inertial_coordinates_device));
  copy_host_array_to_device<volume_dim>(
      packed_geometry->mesh_velocity_host,
      make_not_null(&packed_geometry->mesh_velocity_device));
  ::Kokkos::deep_copy(packed_geometry->element_inverse_jacobian_device,
                      host_inverse_jacobian);
}

template struct UpdateFrameFieldsBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize>;
template struct UpdateFrameFieldsBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize, false>;

}  // namespace gh::Actions
