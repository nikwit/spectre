// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include <pup.h>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tags/MirrorView.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Executables/ScalarWave/Batched/BoundaryBatchMetadata.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution::Kokkos {

template <typename System>
struct PackedEvolutionStateExtras {
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

template <typename System>
struct PackedBoundaryMetadataStorage {
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

template <typename System>
struct PackedBoundaryScratchStorage {
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

template <typename System>
struct PackedSystemTraits {
  using evolution_state_extras = PackedEvolutionStateExtras<System>;
  using boundary_metadata_storage = PackedBoundaryMetadataStorage<System>;
  using boundary_scratch_storage = PackedBoundaryScratchStorage<System>;
};

template <size_t Dim>
struct PackedEvolutionStateExtras<ScalarWave::System<Dim>> {
  using device_constraint_gamma2_type =
      typename ::Tags::MirrorView<ScalarWave::Tags::ConstraintGamma2>::type;
  device_constraint_gamma2_type device_constraint_gamma2{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    size_t device_points = get(device_constraint_gamma2).extent(0);
    p | device_points;
    if (device_points != 0) {
      ERROR(
          "PUP for non-empty ScalarWave packed constraint-gamma2 state is "
          "currently unsupported.");
    }
  }
};

template <size_t Dim>
struct PackedBoundaryMetadataStorage<ScalarWave::System<Dim>> {
  static constexpr size_t volume_dim = Dim;
  static_assert(
      volume_dim == ScalarWave::Batched::boundary_volume_dim,
      "PackedBoundaryMetadataStorage currently supports only ScalarWave 3D "
      "boundary metadata.");

  ScalarWave::Batched::OrientedRemoteFaceIndexMaps
      oriented_remote_face_index_for_local{};
  ScalarWave::Batched::MortarMetadataMaps mortar_metadata{};
  ScalarWave::Batched::FaceBoundaryMetadataArray
      boundary_correction_face_metadata{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | mortar_metadata;
    p | boundary_correction_face_metadata;
    size_t oriented_map_count = oriented_remote_face_index_for_local.size();
    p | oriented_map_count;
    if (oriented_map_count != 0) {
      ERROR(
          "PUP for non-empty ScalarWave packed boundary metadata oriented "
          "maps is currently unsupported.");
    }
  }
};

template <size_t Dim>
struct PackedBoundaryScratchStorage<ScalarWave::System<Dim>> {
  static constexpr size_t volume_dim = Dim;
  static_assert(
      volume_dim == ScalarWave::Batched::boundary_volume_dim,
      "PackedBoundaryScratchStorage currently supports only ScalarWave 3D "
      "boundary metadata.");
  static constexpr size_t number_of_faces =
      ScalarWave::Batched::boundary_number_of_faces;
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
          "PUP for non-empty ScalarWave packed boundary scratch is currently "
          "unsupported.");
    }
  }
};

}  // namespace evolution::Kokkos
