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
  using device_inverse_jacobian_type =
      ::Kokkos::View<double***, ::Kokkos::LayoutRight>;
  std::array<DataVector, volume_dim> inertial_coordinates_host{};
  device_inertial_coordinates_type inertial_coordinates_device{};
  device_inverse_jacobian_type element_inverse_jacobian_device{};

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
