// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include <pup.h>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Executables/ScalarWave/Batched/BoundaryBatchMetadata.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Batched {

struct PackedTopology {
  static constexpr size_t volume_dim = 3;
  using device_face_to_volume_index_map_type =
      typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<
          volume_dim>::type;

  std::vector<ElementId<volume_dim>> local_element_ids{};
  std::vector<Element<volume_dim>> local_elements{};
  std::map<ElementId<volume_dim>, size_t> element_index_by_id{};
  std::vector<std::array<size_t, volume_dim>> element_extents_host{};
  std::vector<size_t> element_point_offsets_host{};
  size_t total_points{0};
  size_t points_per_element{0};
  std::array<size_t, volume_dim> uniform_extents_host{{0, 0, 0}};
  std::array<Spectral::Basis, volume_dim> uniform_basis_host{
      {Spectral::Basis::Legendre, Spectral::Basis::Legendre,
       Spectral::Basis::Legendre}};
  std::array<Spectral::Quadrature, volume_dim> uniform_quadrature_host{
      {Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::GaussLobatto,
       Spectral::Quadrature::GaussLobatto}};
  ::Kokkos::View<size_t*> element_point_offsets_device{};
  ::Kokkos::View<size_t* [volume_dim]> element_extents_device {};
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
          "PUP for non-empty ScalarWave::Batched::PackedTopology device "
          "metadata is currently unsupported.");
    }
  }
};

struct PackedGeometry {
  static constexpr size_t volume_dim = 3;
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
          "PUP for non-empty ScalarWave::Batched::PackedGeometry device "
          "metadata is currently unsupported.");
    }
  }
};

struct PackedBoundaryMetadata {
  OrientedRemoteFaceIndexMaps oriented_remote_face_index_for_local{};
  MortarMetadataMaps mortar_metadata{};
  FaceBoundaryMetadataArray boundary_correction_face_metadata{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | mortar_metadata;
    p | boundary_correction_face_metadata;
    size_t oriented_map_count = oriented_remote_face_index_for_local.size();
    p | oriented_map_count;
    if (oriented_map_count != 0) {
      ERROR(
          "PUP for non-empty "
          "ScalarWave::Batched::PackedBoundaryMetadata oriented maps is "
          "currently unsupported.");
    }
  }
};

struct PackedEvolutionState {
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using device_variables_type =
      typename ScalarWave::KokkosTags::DeviceVariables<system>::type;
  using device_dt_variables_type =
      typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type;
  using device_step_start_type =
      typename ScalarWave::KokkosTags::DeviceStepStart<system>::type;
  using device_derivative_history_type =
      typename ScalarWave::KokkosTags::DeviceDerivativeHistory<system>::type;
  using device_constraint_gamma2_type =
      typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type;

  device_variables_type device_variables{};
  device_dt_variables_type device_dt_variables{};
  device_step_start_type device_step_start{};
  device_derivative_history_type device_derivative_history{};
  device_constraint_gamma2_type device_constraint_gamma2{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    size_t device_points = device_variables.number_of_grid_points();
    p | device_points;
    if (device_points != 0) {
      ERROR(
          "PUP for non-empty ScalarWave::Batched::PackedEvolutionState is "
          "currently unsupported.");
    }
  }
};

struct PackedBoundaryScratch {
  static constexpr size_t volume_dim = 3;
  static constexpr size_t number_of_faces = boundary_number_of_faces;
  using package_field_tags =
      typename ScalarWave::BoundaryCorrections::UpwindPenalty<
          volume_dim>::dg_package_field_tags;
  using device_package_field_tags =
      db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;
  using packaged_face_data_storage_type =
      std::array<Variables<device_package_field_tags>, number_of_faces>;
  using dt_boundary_tags =
      tmpl::list<::Tags::dt<ScalarWave::Tags::Psi>,
                 ::Tags::dt<ScalarWave::Tags::Pi>,
                 ::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>;
  using device_dt_boundary_tags =
      db::wrap_tags_in<::Tags::MirrorView, dt_boundary_tags>;
  using internal_boundary_terms_storage_type =
      std::array<Variables<device_dt_boundary_tags>, number_of_faces>;

  packaged_face_data_storage_type packaged_face_data_for_all_elements{};
  internal_boundary_terms_storage_type
      internal_boundary_terms_for_all_elements{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    size_t packaged_points = 0;
    for (const auto& packaged_face_data : packaged_face_data_for_all_elements) {
      packaged_points += packaged_face_data.number_of_grid_points();
    }
    size_t internal_points = 0;
    for (const auto& internal_boundary_terms :
         internal_boundary_terms_for_all_elements) {
      internal_points += internal_boundary_terms.number_of_grid_points();
    }
    p | packaged_points;
    p | internal_points;
    if (packaged_points != 0 or internal_points != 0) {
      ERROR(
          "PUP for non-empty ScalarWave::Batched::PackedBoundaryScratch is "
          "currently unsupported.");
    }
  }
};

}  // namespace ScalarWave::Batched
