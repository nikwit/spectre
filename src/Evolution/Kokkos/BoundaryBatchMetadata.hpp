// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <pup.h>

#include "Domain/Structure/DirectionalIdMap.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/SegmentSize.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace evolution::Kokkos::Batched {

static constexpr size_t boundary_volume_dim = 3;
static constexpr size_t boundary_number_of_faces = 2 * boundary_volume_dim;
static constexpr size_t boundary_projection_group_count = 16;

struct BoundaryCorrectionWorkItem {
  size_t local_element_index{0};
  size_t remote_element_index{0};
  size_t remote_face_id{0};
  size_t oriented_remote_face_index_offset{0};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::Batched::"
        "BoundaryCorrectionWorkItem, but pup is not supported with kokkos");
  }
};

struct MortarMetadata {
  Mesh<boundary_volume_dim - 1> mortar_mesh{};
  std::array<Spectral::SegmentSize, boundary_volume_dim - 1> mortar_size{};
  bool needs_projection{false};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR("Tried to call pup for evolution::Kokkos::Batched::MortarMetadata, "
          "but pup is not supported with kokkos");
  }
};

struct ProjectionGroupMetadata {
  std::array<Spectral::SegmentSize, boundary_volume_dim - 1> mortar_size{
      {Spectral::SegmentSize::Uninitialized,
       Spectral::SegmentSize::Uninitialized}};
  Mesh<boundary_volume_dim - 1> mortar_mesh{};
  ::Kokkos::View<BoundaryCorrectionWorkItem*> work_items{};
  ::Kokkos::View<size_t*> oriented_remote_face_indices{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR(
        "Tried to call pup for evolution::Kokkos::Batched::"
        "ProjectionGroupMetadata, but pup is not supported with kokkos");
  }
};

struct FaceBoundaryMetadata {
  ::Kokkos::View<BoundaryCorrectionWorkItem*> no_projection_work_items{};
  ::Kokkos::View<size_t*> no_projection_oriented_remote_face_indices{};
  std::array<ProjectionGroupMetadata, boundary_projection_group_count>
      projection_groups{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {
    ERROR("Tried to call pup for evolution::Kokkos::Batched::"
          "FaceBoundaryMetadata, but pup is not supported with kokkos");
  }
};

using OrientedRemoteFaceIndexMap =
    DirectionalIdMap<boundary_volume_dim, ::Kokkos::View<size_t*>>;
using MortarMetadataMap = DirectionalIdMap<boundary_volume_dim, MortarMetadata>;
using OrientedRemoteFaceIndexMaps = std::vector<OrientedRemoteFaceIndexMap>;
using MortarMetadataMaps = std::vector<MortarMetadataMap>;
using FaceBoundaryMetadataArray =
    std::array<FaceBoundaryMetadata, boundary_number_of_faces>;

}  // namespace evolution::Kokkos::Batched
