// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CoordinateMaps/TimeDependent/MapDescriptor.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Kokkos/SystemTraits.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace evolution::Kokkos {

template <typename System>
using packed_boundary_metadata_storage_t =
    typename PackedSystemTraits<System>::boundary_metadata_storage;

template <typename System>
using packed_boundary_scratch_storage_t =
    typename PackedSystemTraits<System>::boundary_scratch_storage;

template <typename System>
struct PackedTopology {
  static constexpr size_t volume_dim = System::volume_dim;
  using device_face_to_volume_index_map_type =
      typename evolution::Kokkos::Tags::DeviceFaceToVolumeIndexMap<
          volume_dim>::type;

  std::vector<ElementId<volume_dim>> local_element_ids{};
  std::vector<Element<volume_dim>> local_elements{};
  std::map<ElementId<volume_dim>, size_t> element_index_by_id{};
  std::vector<std::array<size_t, volume_dim>> element_extents_host{};
  std::vector<size_t> element_point_offsets_host{};
  size_t total_points{0};
  size_t points_per_element{0};
  std::array<size_t, volume_dim> uniform_extents_host{};
  std::array<Spectral::Basis, volume_dim> uniform_basis_host{};
  std::array<Spectral::Quadrature, volume_dim> uniform_quadrature_host{};
  ::Kokkos::View<size_t*> element_point_offsets_device{};
  ::Kokkos::View<size_t* [volume_dim]> element_extents_device {};
  device_face_to_volume_index_map_type device_face_to_volume_index_map{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedTopology, but pup "
        "is not supported with kokkos");
  }
};

template <typename System>
struct PackedGeometry {
  static constexpr size_t volume_dim = System::volume_dim;
  using device_inertial_coordinates_type =
      typename evolution::Kokkos::Tags::DeviceInertialCoordinates<
          volume_dim>::type;
  using device_mesh_velocity_type =
      typename evolution::Kokkos::Tags::DeviceInertialCoordinates<
          volume_dim>::type;
  using device_inverse_jacobian_type =
      ::Kokkos::View<double***, ::Kokkos::LayoutRight>;
  using host_inverse_jacobian_cache_type =
      ::Kokkos::View<double***, ::Kokkos::LayoutRight, ::Kokkos::HostSpace>;
  using device_affine_map_data_type =
      ::Kokkos::View<double**, ::Kokkos::LayoutRight>;
  using device_shape_coefficients_type =
      ::Kokkos::View<double***, ::Kokkos::LayoutRight>;
  using device_shape_point_data_type =
      ::Kokkos::View<double***, ::Kokkos::LayoutRight>;
  std::array<DataVector, volume_dim> grid_coordinates_host{};
  std::array<DataVector, volume_dim> distorted_coordinates_host{};
  std::array<DataVector, volume_dim> inertial_coordinates_host{};
  std::array<DataVector, volume_dim> mesh_velocity_host{};
  std::array<DataVector, volume_dim> grid_to_distorted_frame_velocity_host{};
  host_inverse_jacobian_cache_type logical_to_grid_inverse_jacobian_host{};
  host_inverse_jacobian_cache_type logical_to_distorted_inverse_jacobian_host{};
  host_inverse_jacobian_cache_type grid_to_distorted_inverse_jacobian_host{};
  host_inverse_jacobian_cache_type distorted_to_inertial_jacobian_host{};
  host_inverse_jacobian_cache_type
      distorted_to_inertial_inverse_jacobian_host{};
  device_inertial_coordinates_type inertial_coordinates_device{};
  device_mesh_velocity_type mesh_velocity_device{};
  device_inverse_jacobian_type element_inverse_jacobian_device{};
  // Intermediate-frame geometry is device-resident for map compositions with
  // a distorted frame, then mirrored to the existing host caches for
  // observation and apparent-horizon consumers.
  device_inertial_coordinates_type distorted_coordinates_device{};
  device_mesh_velocity_type grid_to_distorted_frame_velocity_device{};
  device_inverse_jacobian_type logical_to_distorted_inverse_jacobian_device{};
  device_inverse_jacobian_type grid_to_distorted_inverse_jacobian_device{};
  device_inverse_jacobian_type distorted_to_inertial_jacobian_device{};
  device_inverse_jacobian_type distorted_to_inertial_inverse_jacobian_device{};
  // Device-resident copies of the domain-static geometry. These are uploaded
  // once (they never change during an evolution) and are consumed by the
  // on-device frame-field update for blocks whose coordinate map is evaluated
  // on the device (see gh::Actions::BatchedMapKind).
  device_inertial_coordinates_type grid_coordinates_device{};
  device_inverse_jacobian_type logical_to_grid_inverse_jacobian_device{};
  // Per-element device coordinate-map plan: the gh::Actions::BatchedMapKind of
  // each element's block, fixed at initialization. Drives the on-device
  // frame-field kernel's per-element dispatch (and is the natural key for a
  // later group-by-map-kind launch optimization). Empty until built.
  std::vector<int> element_map_kind_host{};
  ::Kokkos::View<int*> element_map_kind_device{};
  ::Kokkos::View<int*> element_has_distorted_frame_device{};
  // Exact host-side descriptions of each block's grid-to-inertial map. An
  // empty optional marks a composition containing a map that has not yet
  // implemented the descriptor interface. These are captured once while the
  // map plan is built and are used to stage compact device parameters without
  // probing the type-erased coordinate map at every time step.
  std::vector<std::optional<
      std::vector<domain::CoordinateMaps::TimeDependent::MapDescriptor>>>
      block_map_descriptors_host{};
  // Per-element coefficients for maps classified as ComposedAffine. The
  // packed layout contains the coordinate matrix and translation, the affine
  // frame-velocity matrix and translation, and the grid-to-inertial inverse
  // Jacobian. Coefficients are refreshed once per stage and consumed directly
  // by the device frame-field kernel.
  device_affine_map_data_type affine_map_data_device{};
  // Shape-map data. The coefficient view stores, per distinct Shape
  // function-of-time group, the spectral expansions of lambda, dt(lambda),
  // and its three Cartesian angular derivatives. Elements select their shared
  // coefficient set through shape_group_index_device. Point data stores the
  // static transition function, its radius-divided value, and its Cartesian
  // gradient.
  device_shape_coefficients_type shape_coefficients_device{};
  device_shape_point_data_type shape_point_data_device{};
  ::Kokkos::View<size_t*> shape_group_index_device{};
  std::vector<size_t> shape_group_index_host{};
  std::vector<size_t> shape_group_representative_element_host{};
  ::Kokkos::View<double**, ::Kokkos::LayoutRight>
      shape_recurrence_alpha_device{};
  ::Kokkos::View<double**, ::Kokkos::LayoutRight>
      shape_recurrence_beta_device{};
  ::Kokkos::View<double**, ::Kokkos::LayoutRight> shape_recurrence_pmm_device{};
  ::Kokkos::View<size_t**, ::Kokkos::LayoutRight>
      shape_recurrence_indices_device{};
  ::Kokkos::View<double* [volume_dim]> shape_centers_device {};
  ::Kokkos::View<size_t*> shape_lmax_device{};
  std::vector<size_t> shape_lmax_host{};
  // Name of the single scalar function of time driving the on-device
  // time-dependent map blocks (empty if there are none): the rotation angle for
  // RigidRotationZ or the scale factor for UniformRadialScaling. Its value and
  // first time derivative are staged to the device each stage.
  std::string device_map_function_of_time_name{};
  bool map_plan_initialized{false};
  // Device map evaluation invalidates the host coordinate-map caches. Events
  // synchronize them lazily, avoiding device-to-host transfers on ordinary
  // RHS stages.
  bool host_frame_fields_dirty{false};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedGeometry, but pup "
        "is not supported with kokkos");
  }
};

template <typename System>
struct PackedBoundaryMetadata : packed_boundary_metadata_storage_t<System> {
 private:
  using boundary_metadata_storage = packed_boundary_metadata_storage_t<System>;

 public:
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedBoundaryMetadata, "
        "but pup is not supported with kokkos");
  }
};

template <typename System>
struct PackedEvolutionState {
  using device_variables_type =
      typename evolution::Kokkos::Tags::DeviceVariables<System>::type;
  using device_dt_variables_type =
      typename evolution::Kokkos::Tags::DeviceDtVariables<System>::type;
  using device_step_start_type =
      typename evolution::Kokkos::Tags::DeviceStepStart<System>::type;
  using device_derivative_history_type =
      typename evolution::Kokkos::Tags::DeviceDerivativeHistory<System>::type;

  device_variables_type device_variables{};
  device_dt_variables_type device_dt_variables{};
  device_step_start_type device_step_start{};
  device_derivative_history_type device_derivative_history{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedEvolutionState, but "
        "pup is not supported with kokkos");
  }
};

template <typename System>
struct PackedBoundaryScratch : packed_boundary_scratch_storage_t<System> {
 private:
  using boundary_scratch_storage = packed_boundary_scratch_storage_t<System>;

 public:
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedBoundaryScratch, "
        "but pup is not supported with kokkos");
  }
};

}  // namespace evolution::Kokkos
