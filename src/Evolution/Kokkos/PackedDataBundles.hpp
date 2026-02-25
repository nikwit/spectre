// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
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
using packed_evolution_state_extras_t =
    typename PackedSystemTraits<System>::evolution_state_extras;

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
  ::Kokkos::View<size_t* [volume_dim]> element_extents_device{};
  device_face_to_volume_index_map_type device_face_to_volume_index_map{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | local_element_ids;
    p | local_elements;
    p | element_index_by_id;
    p | element_extents_host;
    p | element_point_offsets_host;
    p | total_points;
    p | points_per_element;
    p | uniform_extents_host;
    p | uniform_basis_host;
    p | uniform_quadrature_host;

    size_t point_offsets_size = element_point_offsets_device.extent(0);
    size_t extents_size = element_extents_device.extent(0);
    size_t face_to_volume_points = 0;
    for (size_t d = 0; d < volume_dim; ++d) {
      face_to_volume_points +=
          gsl::at(device_face_to_volume_index_map, d).first.extent(0);
      face_to_volume_points +=
          gsl::at(device_face_to_volume_index_map, d).second.extent(0);
    }
    p | point_offsets_size;
    p | extents_size;
    p | face_to_volume_points;
    if (point_offsets_size != 0 or extents_size != 0 or
        face_to_volume_points != 0) {
      ERROR(
          "PUP for non-empty evolution::Kokkos::PackedTopology device "
          "metadata is currently unsupported.");
    }
  }
};

template <typename System>
struct PackedGeometry {
  static constexpr size_t volume_dim = System::volume_dim;
  std::array<DataVector, volume_dim> inertial_coordinates_host{};
  ::Kokkos::View<double***> element_inverse_jacobian_device{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | inertial_coordinates_host;
    size_t inverse_jacobian_num_elements =
        element_inverse_jacobian_device.extent(0);
    p | inverse_jacobian_num_elements;
    if (inverse_jacobian_num_elements != 0) {
      ERROR(
          "PUP for non-empty evolution::Kokkos::PackedGeometry device "
          "metadata is currently unsupported.");
    }
  }
};

template <typename System>
struct PackedBoundaryMetadata : packed_boundary_metadata_storage_t<System> {
 private:
  using boundary_metadata_storage = packed_boundary_metadata_storage_t<System>;

 public:
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) { boundary_metadata_storage::pup(p); }
};

template <typename System>
struct PackedEvolutionState : packed_evolution_state_extras_t<System> {
 private:
  using evolution_state_extras = packed_evolution_state_extras_t<System>;

 public:
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
  void pup(PUP::er& p) {
    size_t device_points = device_variables.number_of_grid_points();
    p | device_points;
    if (device_points != 0) {
      ERROR(
          "PUP for non-empty evolution::Kokkos::PackedEvolutionState is "
          "currently unsupported.");
    }
    evolution_state_extras::pup(p);
  }
};

template <typename System>
struct PackedBoundaryScratch : packed_boundary_scratch_storage_t<System> {
 private:
  using boundary_scratch_storage = packed_boundary_scratch_storage_t<System>;

 public:
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) { boundary_scratch_storage::pup(p); }
};

}  // namespace evolution::Kokkos
