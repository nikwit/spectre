// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include <pup.h>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tags/MirrorView.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Kokkos/BoundaryBatchMetadata.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution::Kokkos {

template <typename System>
struct PackedBoundaryMetadataStorage {
  static constexpr size_t volume_dim = System::volume_dim;
  static_assert(
      volume_dim == evolution::Kokkos::Batched::boundary_volume_dim,
      "PackedBoundaryMetadataStorage currently supports only 3D boundary "
      "metadata.");

  evolution::Kokkos::Batched::OrientedRemoteFaceIndexMaps
      oriented_remote_face_index_for_local{};
  evolution::Kokkos::Batched::MortarMetadataMaps mortar_metadata{};
  evolution::Kokkos::Batched::FaceBoundaryMetadataArray
      boundary_correction_face_metadata{};
  std::array<::Kokkos::View<int*>,
             evolution::Kokkos::Batched::boundary_number_of_faces>
      external_face_mask_for_all_elements{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for "
        "evolution::Kokkos::PackedBoundaryMetadataStorage,"
        " but pup is not supported with kokkos");
  }
};

template <typename System, typename BoundaryCorrection>
struct PackedBoundaryScratchStorage {
  static constexpr size_t volume_dim = System::volume_dim;
  static_assert(
      volume_dim == evolution::Kokkos::Batched::boundary_volume_dim,
      "PackedBoundaryScratchStorage currently supports only 3D boundary "
      "metadata.");
  static constexpr size_t number_of_faces =
      evolution::Kokkos::Batched::boundary_number_of_faces;
  using package_field_tags = typename BoundaryCorrection::dg_package_field_tags;
  using device_package_field_tags =
      db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;
  using packaged_face_data_storage_type =
      std::array<Variables<device_package_field_tags>, number_of_faces>;
  using dt_boundary_tags =
      db::wrap_tags_in<::Tags::dt, typename System::variables_tag::tags_list>;
  using device_dt_boundary_tags =
      db::wrap_tags_in<::Tags::MirrorView, dt_boundary_tags>;
  using internal_boundary_terms_storage_type =
      std::array<Variables<device_dt_boundary_tags>, number_of_faces>;

  packaged_face_data_storage_type packaged_face_data_for_all_elements{};
  internal_boundary_terms_storage_type
      internal_boundary_terms_for_all_elements{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedBoundaryScratchStorage,"
        " but pup is not supported with kokkos");
  }
};

template <typename System>
struct PackedBoundaryScratchStorage<System, void> {
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::PackedBoundaryScratchStorage,"
        " but pup is not supported with kokkos");
  }
};

template <typename System>
struct PackedSystemTraits {
  using boundary_metadata_storage = PackedBoundaryMetadataStorage<System>;
  using boundary_scratch_storage = PackedBoundaryScratchStorage<System, void>;
};

template <size_t Dim>
struct PackedSystemTraits<ScalarWave::System<Dim>> {
  using boundary_metadata_storage =
      PackedBoundaryMetadataStorage<ScalarWave::System<Dim>>;
  using boundary_scratch_storage = PackedBoundaryScratchStorage<
      ScalarWave::System<Dim>,
      ScalarWave::BoundaryCorrections::UpwindPenalty<Dim>>;
};

template <size_t Dim>
struct PackedSystemTraits<gh::System<Dim>> {
  using boundary_metadata_storage =
      PackedBoundaryMetadataStorage<gh::System<Dim>>;
  using boundary_scratch_storage =
      PackedBoundaryScratchStorage<gh::System<Dim>,
                                   gh::BoundaryCorrections::UpwindPenalty<Dim>>;
};

}  // namespace evolution::Kokkos
