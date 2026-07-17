// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BatchedMapPlan.hpp"
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
  const auto ensure_host_jacobian =
      [num_elements, points_per_element](
          const gsl::not_null<
              typename packed_geometry_type::host_inverse_jacobian_cache_type*>
              jacobian,
          const std::string& label) {
        if (jacobian->extent(0) != num_elements or
            jacobian->extent(1) != points_per_element or
            jacobian->extent(2) != VolumeDim * VolumeDim) {
          *jacobian =
              typename packed_geometry_type::host_inverse_jacobian_cache_type{
                  label, num_elements, points_per_element,
                  VolumeDim * VolumeDim};
        }
      };
  ensure_host_jacobian(
      make_not_null(
          &packed_geometry->logical_to_distorted_inverse_jacobian_host),
      "GhBatchedLogicalToDistortedInverseJacobian");
  ensure_host_jacobian(
      make_not_null(&packed_geometry->grid_to_distorted_inverse_jacobian_host),
      "GhBatchedGridToDistortedInverseJacobian");
  ensure_host_jacobian(
      make_not_null(&packed_geometry->distorted_to_inertial_jacobian_host),
      "GhBatchedDistortedToInertialJacobian");
  ensure_host_jacobian(
      make_not_null(
          &packed_geometry->distorted_to_inertial_inverse_jacobian_host),
      "GhBatchedDistortedToInertialInverseJacobian");

  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->distorted_coordinates_host),
      total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->inertial_coordinates_host), total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->mesh_velocity_host), total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->grid_to_distorted_frame_velocity_host),
      total_points);

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

// Upload the domain-static geometry (grid coordinates and the logical-to-grid
// inverse Jacobian) to the device once. These never change during an
// evolution, so this is a one-time cost guarded by the current device extents.
template <size_t VolumeDim>
void ensure_static_device_geometry(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t num_elements, const size_t points_per_element,
    const size_t total_points) {
  using packed_geometry_type =
      evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>;

  bool coordinates_ready = true;
  for (size_t d = 0; d < VolumeDim; ++d) {
    if (packed_geometry->grid_coordinates_device.get(d).extent(0) !=
        total_points) {
      coordinates_ready = false;
    }
  }
  if (not coordinates_ready) {
    for (size_t d = 0; d < VolumeDim; ++d) {
      packed_geometry->grid_coordinates_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedGridCoordinates", total_points);
    }
    copy_host_array_to_device<VolumeDim>(
        packed_geometry->grid_coordinates_host,
        make_not_null(&packed_geometry->grid_coordinates_device));
  }

  if (packed_geometry->logical_to_grid_inverse_jacobian_device.extent(0) !=
          num_elements or
      packed_geometry->logical_to_grid_inverse_jacobian_device.extent(1) !=
          points_per_element or
      packed_geometry->logical_to_grid_inverse_jacobian_device.extent(2) !=
          VolumeDim * VolumeDim) {
    packed_geometry->logical_to_grid_inverse_jacobian_device =
        typename packed_geometry_type::device_inverse_jacobian_type{
            "GhBatchedLogicalToGridInverseJacobianDevice", num_elements,
            points_per_element, VolumeDim * VolumeDim};
    ::Kokkos::deep_copy(
        packed_geometry->logical_to_grid_inverse_jacobian_device,
        packed_geometry->logical_to_grid_inverse_jacobian_host);
  }
}

// Device frame-field update for a domain whose every block is evaluated on the
// device (a mix of BatchedMapKind::Static and BatchedMapKind::RigidRotationZ).
// Writes the device-resident geometry consumed by the RHS/boundary kernels
// directly on the device, dispatching per element on the map kind. `theta` and
// `theta_dot` are the shared rigid-rotation angle and angular velocity at the
// current time (unused for static-only domains). This is the enum-dispatch form
// of the device map path; the per-element kind is also the natural key for a
// later group-by-kind launch optimization.
template <size_t VolumeDim>
void update_frame_fields_device(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t points_per_element, const size_t total_points,
    const double theta, const double theta_dot) {
  const auto grid_coordinates = packed_geometry->grid_coordinates_device;
  const auto inertial_coordinates =
      packed_geometry->inertial_coordinates_device;
  const auto mesh_velocity = packed_geometry->mesh_velocity_device;
  const auto logical_to_grid_inverse_jacobian =
      packed_geometry->logical_to_grid_inverse_jacobian_device;
  const auto element_inverse_jacobian =
      packed_geometry->element_inverse_jacobian_device;
  const auto element_map_kind = packed_geometry->element_map_kind_device;
  const size_t ppe = points_per_element;
  constexpr int static_kind = static_cast<int>(BatchedMapKind::Static);
  ::Kokkos::parallel_for(
      "GhUpdateFrameFields", total_points,
      KOKKOS_LAMBDA(const int linear_index_int) {
        const size_t s = static_cast<size_t>(linear_index_int);
        const size_t e = s / ppe;
        const size_t p = s % ppe;
        if (element_map_kind(e) == static_kind) {
          for (size_t d = 0; d < VolumeDim; ++d) {
            inertial_coordinates.get(d)(s) = grid_coordinates.get(d)(s);
            mesh_velocity.get(d)(s) = 0.0;
          }
          for (size_t component = 0; component < VolumeDim * VolumeDim;
               ++component) {
            element_inverse_jacobian(e, p, component) =
                logical_to_grid_inverse_jacobian(e, p, component);
          }
        } else {
          // BatchedMapKind::RigidRotationZ
          double inertial[3];
          double velocity[3];
          double grid_to_inertial_inv_jac[9];
          rigid_rotation_z_frame_fields(
              grid_coordinates.get(0)(s), grid_coordinates.get(1)(s),
              grid_coordinates.get(2)(s), theta, theta_dot, inertial, velocity,
              grid_to_inertial_inv_jac);
          for (size_t d = 0; d < VolumeDim; ++d) {
            inertial_coordinates.get(d)(s) = inertial[d];
            mesh_velocity.get(d)(s) = velocity[d];
          }
          // Compose logical->grid (static) with grid->inertial (rotation):
          // invjac(logical, inertial) = sum_grid L2G(logical, grid) *
          //                             G2I(grid, inertial).
          for (size_t logical_d = 0; logical_d < VolumeDim; ++logical_d) {
            for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
              double composed = 0.0;
              for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
                composed += logical_to_grid_inverse_jacobian(
                                e, p, logical_d * VolumeDim + grid_d) *
                            grid_to_inertial_inv_jac[grid_d * VolumeDim +
                                                     inertial_d];
              }
              element_inverse_jacobian(e, p,
                                       logical_d * VolumeDim + inertial_d) =
                  composed;
            }
          }
        }
      });
}

// Copy a device coordinate/velocity array back into its host mirror. Used to
// keep the host geometry arrays current for the host-side observation and
// apparent-horizon consumers when the device path evaluates a moving map.
template <size_t VolumeDim>
void copy_device_array_to_host(
    const typename evolution::Kokkos::PackedGeometry<
        gh::System<VolumeDim>>::device_inertial_coordinates_type& device_data,
    const gsl::not_null<std::array<DataVector, VolumeDim>*> host_data) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    const auto mirror = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, device_data.get(d));
    const size_t num_points = device_data.get(d).extent(0);
    if (gsl::at(*host_data, d).size() != num_points) {
      gsl::at(*host_data, d) = DataVector{num_points, 0.0};
    }
    for (size_t s = 0; s < num_points; ++s) {
      gsl::at(*host_data, d)[s] = mirror(s);
    }
  }
}

// Verify that a block's grid-to-inertial map is, at the given time, the rigid
// z-rotation closed form used by the device path (rigid_rotation_z_frame_fields
// with the angle and angular velocity from `fot_name`). Compares the host map's
// coordinates, mesh velocity, and grid-to-inertial inverse Jacobian at the
// element's grid points. This is the safety gate for the name-based
// recognizer: only if the host map provably matches the device closed form is
// the block evaluated on the device; otherwise it falls back to the host path.
template <size_t VolumeDim>
bool verify_rigid_rotation_z(
    const Block<VolumeDim>& block,
    const tnsr::I<DataVector, VolumeDim, Frame::Grid>& grid_coords,
    const std::string& fot_name, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  const auto map_data =
      block.moving_mesh_grid_to_inertial_map().coords_frame_velocity_jacobians(
          grid_coords, time, functions_of_time);
  const auto& inertial_coords = std::get<0>(map_data);
  const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
  const auto& mesh_velocity = std::get<3>(map_data);

  const auto angle_and_deriv =
      functions_of_time.at(fot_name)->func_and_deriv(time);
  if (angle_and_deriv[0].size() != 1 or angle_and_deriv[1].size() != 1) {
    return false;
  }
  const double theta = angle_and_deriv[0][0];
  const double theta_dot = angle_and_deriv[1][0];

  const double tolerance = 1.0e-11;
  const auto within_tolerance = [tolerance](const double a, const double b) {
    return std::abs(a - b) <= tolerance * (1.0 + std::abs(b));
  };

  const size_t num_points = grid_coords.get(0).size();
  for (size_t p = 0; p < num_points; ++p) {
    double inertial[3];
    double velocity[3];
    double grid_to_inertial_inv_jac[9];
    rigid_rotation_z_frame_fields(
        grid_coords.get(0)[p], grid_coords.get(1)[p], grid_coords.get(2)[p],
        theta, theta_dot, inertial, velocity, grid_to_inertial_inv_jac);
    for (size_t d = 0; d < VolumeDim; ++d) {
      if (not within_tolerance(inertial_coords.get(d)[p], inertial[d])) {
        return false;
      }
      if (not within_tolerance(mesh_velocity.get(d)[p], velocity[d])) {
        return false;
      }
    }
    for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
      for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
        if (not within_tolerance(
                grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p],
                grid_to_inertial_inv_jac[grid_d * VolumeDim + inertial_d])) {
          return false;
        }
      }
    }
  }
  return true;
}

// Build the per-element device map plan once: classify every element's block
// as Static, RigidRotationZ (verified against the host map), or HostEvaluated,
// record the shared rotation function-of-time name, and upload the per-element
// kind array to the device.
template <size_t VolumeDim>
void build_map_plan(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const evolution::Kokkos::PackedTopology<gh::System<VolumeDim>>&
        packed_topology,
    const Domain<VolumeDim>& domain,
    const domain::FunctionsOfTimeMap& functions_of_time, const double time) {
  const size_t num_elements = packed_topology.local_element_ids.size();
  const size_t points_per_element = packed_topology.points_per_element;
  packed_geometry->element_map_kind_host.assign(
      num_elements, static_cast<int>(BatchedMapKind::HostEvaluated));

  std::map<size_t, BatchedMapKind> kind_by_block{};
  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& block = domain.blocks()[element_id.block_id()];
    const size_t block_id = element_id.block_id();

    BatchedMapKind kind = BatchedMapKind::HostEvaluated;
    const auto cached = kind_by_block.find(block_id);
    if (cached != kind_by_block.end()) {
      kind = cached->second;
    } else {
      kind = classify_block_map(block);
      if (kind != BatchedMapKind::Static) {
        // Attempt to recognize a device-supported time-dependent map. A rigid
        // z-rotation has a single scalar function of time and no distorted
        // frame; verify it against the host map before enabling the device
        // path.
        const auto& names =
            block.moving_mesh_grid_to_inertial_map().function_of_time_names();
        if (names.size() == 1 and not block.has_distorted_frame()) {
          const std::string& name = *names.begin();
          const size_t point_begin =
              packed_topology.element_point_offsets_host[e];
          tnsr::I<DataVector, VolumeDim, Frame::Grid> grid_coords{
              points_per_element};
          for (size_t d = 0; d < VolumeDim; ++d) {
            for (size_t p = 0; p < points_per_element; ++p) {
              grid_coords.get(d)[p] =
                  packed_geometry->grid_coordinates_host[d][point_begin + p];
            }
          }
          if (verify_rigid_rotation_z<VolumeDim>(block, grid_coords, name, time,
                                                 functions_of_time)) {
            kind = BatchedMapKind::RigidRotationZ;
            if (packed_geometry->rotation_function_of_time_name.empty()) {
              packed_geometry->rotation_function_of_time_name = name;
            } else {
              ASSERT(packed_geometry->rotation_function_of_time_name == name,
                     "Batched device rigid rotation currently requires a "
                     "single shared rotation function of time, but found both '"
                         << packed_geometry->rotation_function_of_time_name
                         << "' and '" << name << "'.");
            }
          }
        }
      }
      kind_by_block.emplace(block_id, kind);
    }
    packed_geometry->element_map_kind_host[e] = static_cast<int>(kind);
  }

  packed_geometry->element_map_kind_device =
      ::Kokkos::View<int*>("GhBatchedElementMapKind", num_elements);
  auto host_kind = ::Kokkos::create_mirror_view(
      packed_geometry->element_map_kind_device);
  for (size_t e = 0; e < num_elements; ++e) {
    host_kind(e) = packed_geometry->element_map_kind_host[e];
  }
  ::Kokkos::deep_copy(packed_geometry->element_map_kind_device, host_kind);
  packed_geometry->map_plan_initialized = true;
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

  // Classify each block's map once. Static and rigid-z-rotation blocks are
  // evaluated on the device; any other time-dependent map falls back to the
  // host path below.
  if (not packed_geometry->map_plan_initialized) {
    ensure_static_device_geometry<volume_dim>(
        packed_geometry, num_elements, points_per_element, total_points);
    build_map_plan<volume_dim>(packed_geometry, packed_topology, domain,
                               functions_of_time, time);
  }

  bool has_host_evaluated = false;
  bool all_static = true;
  for (size_t e = 0; e < num_elements; ++e) {
    const int kind = packed_geometry->element_map_kind_host[e];
    if (kind == static_cast<int>(BatchedMapKind::HostEvaluated)) {
      has_host_evaluated = true;
      break;
    }
    if (kind != static_cast<int>(BatchedMapKind::Static)) {
      all_static = false;
    }
  }

  if (not has_host_evaluated) {
    ensure_static_device_geometry<volume_dim>(
        packed_geometry, num_elements, points_per_element, total_points);
    double theta = 0.0;
    double theta_dot = 0.0;
    if (not packed_geometry->rotation_function_of_time_name.empty()) {
      const auto angle_and_deriv =
          functions_of_time
              .at(packed_geometry->rotation_function_of_time_name)
              ->func_and_deriv(time);
      theta = angle_and_deriv[0][0];
      theta_dot = angle_and_deriv[1][0];
    }
    update_frame_fields_device<volume_dim>(
        packed_geometry, points_per_element, total_points, theta, theta_dot);
    // For a static domain the host geometry arrays are constant and were set
    // during initialization, so leave them untouched. When a moving map is
    // evaluated on the device, refresh the host coordinate and mesh-velocity
    // arrays so the host-side observation and apparent-horizon consumers see
    // the current geometry. (Narrowing this to observation/AH times is a later
    // milestone.)
    if (not all_static) {
      copy_device_array_to_host<volume_dim>(
          packed_geometry->inertial_coordinates_device,
          make_not_null(&packed_geometry->inertial_coordinates_host));
      copy_device_array_to_host<volume_dim>(
          packed_geometry->mesh_velocity_device,
          make_not_null(&packed_geometry->mesh_velocity_host));
    }
    return;
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
      const auto& grid_to_inertial_jacobian = std::get<2>(map_data);
      const auto& mesh_velocity = std::get<3>(map_data);

      std::optional<std::tuple<
          tnsr::I<DataVector, volume_dim, Frame::Distorted>,
          InverseJacobian<DataVector, volume_dim, Frame::Grid,
                          Frame::Distorted>,
          Jacobian<DataVector, volume_dim, Frame::Grid, Frame::Distorted>,
          tnsr::I<DataVector, volume_dim, Frame::Distorted>>>
          grid_to_distorted_map_data{};
      std::optional<std::tuple<
          tnsr::I<DataVector, volume_dim, Frame::Inertial>,
          InverseJacobian<DataVector, volume_dim, Frame::Distorted,
                          Frame::Inertial>,
          Jacobian<DataVector, volume_dim, Frame::Distorted, Frame::Inertial>,
          tnsr::I<DataVector, volume_dim, Frame::Inertial>>>
          distorted_to_inertial_map_data{};
      if (block.has_distorted_frame()) {
        grid_to_distorted_map_data =
            block.moving_mesh_grid_to_distorted_map()
                .coords_frame_velocity_jacobians(grid_coords, time,
                                                 functions_of_time);
        distorted_to_inertial_map_data =
            block.moving_mesh_distorted_to_inertial_map()
                .coords_frame_velocity_jacobians(
                    std::get<0>(*grid_to_distorted_map_data), time,
                    functions_of_time);
      }

      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              inertial_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] =
              mesh_velocity.get(d)[p];
          packed_geometry->distorted_coordinates_host[d][point_begin + p] =
              block.has_distorted_frame()
                  ? std::get<0>(*grid_to_distorted_map_data).get(d)[p]
                  : grid_coords.get(d)[p];
          packed_geometry
              ->grid_to_distorted_frame_velocity_host[d][point_begin + p] =
              block.has_distorted_frame()
                  ? std::get<3>(*grid_to_distorted_map_data).get(d)[p]
                  : 0.0;
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
        for (size_t source_d = 0; source_d < volume_dim; ++source_d) {
          for (size_t target_d = 0; target_d < volume_dim; ++target_d) {
            const double grid_to_distorted_inverse_jacobian =
                block.has_distorted_frame()
                    ? std::get<1>(*grid_to_distorted_map_data)
                          .get(source_d, target_d)[p]
                    : (source_d == target_d ? 1.0 : 0.0);
            packed_geometry->grid_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                grid_to_distorted_inverse_jacobian;
            packed_geometry->distorted_to_inertial_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                block.has_distorted_frame()
                    ? std::get<2>(*distorted_to_inertial_map_data)
                          .get(source_d, target_d)[p]
                    : grid_to_inertial_jacobian.get(source_d, target_d)[p];
            packed_geometry->distorted_to_inertial_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                block.has_distorted_frame()
                    ? std::get<1>(*distorted_to_inertial_map_data)
                          .get(source_d, target_d)[p]
                    : grid_to_inertial_inverse_jacobian.get(source_d,
                                                            target_d)[p];

            double composed_logical_to_distorted_inverse_jacobian = 0.0;
            for (size_t grid_d = 0; grid_d < volume_dim; ++grid_d) {
              composed_logical_to_distorted_inverse_jacobian +=
                  packed_geometry->logical_to_grid_inverse_jacobian_host(
                      e, p, source_d * volume_dim + grid_d) *
                  (block.has_distorted_frame()
                       ? std::get<1>(*grid_to_distorted_map_data)
                             .get(grid_d, target_d)[p]
                       : (grid_d == target_d ? 1.0 : 0.0));
            }
            packed_geometry->logical_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                composed_logical_to_distorted_inverse_jacobian;
          }
        }
      }
    } else {
      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->distorted_coordinates_host[d][point_begin + p] =
              grid_coords.get(d)[p];
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              grid_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] = 0.0;
          packed_geometry
              ->grid_to_distorted_frame_velocity_host[d][point_begin + p] = 0.0;
        }
        for (size_t component = 0; component < volume_dim * volume_dim;
             ++component) {
          host_inverse_jacobian(e, p, component) =
              packed_geometry->logical_to_grid_inverse_jacobian_host(e, p,
                                                                     component);
          packed_geometry->logical_to_distorted_inverse_jacobian_host(
              e, p, component) =
              packed_geometry->logical_to_grid_inverse_jacobian_host(e, p,
                                                                     component);
        }
        for (size_t source_d = 0; source_d < volume_dim; ++source_d) {
          for (size_t target_d = 0; target_d < volume_dim; ++target_d) {
            packed_geometry->grid_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
            packed_geometry->distorted_to_inertial_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
            packed_geometry->distorted_to_inertial_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
          }
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
template struct UpdateFrameFieldsBatchedWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize>;
template struct UpdateFrameFieldsBatchedWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize, false>;

}  // namespace gh::Actions
