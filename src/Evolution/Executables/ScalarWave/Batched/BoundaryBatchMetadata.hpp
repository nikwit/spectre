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

namespace ScalarWave::Batched {

static constexpr size_t boundary_volume_dim = 3;
static constexpr size_t boundary_number_of_faces = 2 * boundary_volume_dim;
static constexpr size_t boundary_projection_group_count = 16;

struct BoundaryCorrectionWorkItem {
  size_t local_element_index{0};
  size_t remote_element_index{0};
  size_t remote_face_id{0};
  size_t oriented_remote_face_index_offset{0};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | local_element_index;
    p | remote_element_index;
    p | remote_face_id;
    p | oriented_remote_face_index_offset;
  }
};

struct MortarMetadata {
  Mesh<boundary_volume_dim - 1> mortar_mesh{};
  std::array<Spectral::SegmentSize, boundary_volume_dim - 1> mortar_size{};
  bool needs_projection{false};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | mortar_mesh;
    p | mortar_size;
    p | needs_projection;
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
  void pup(PUP::er& p) {
    p | mortar_size;
    p | mortar_mesh;
    size_t work_items_size = work_items.extent(0);
    size_t oriented_indices_size = oriented_remote_face_indices.extent(0);
    p | work_items_size;
    p | oriented_indices_size;
    if (work_items_size != 0 or oriented_indices_size != 0) {
      ERROR(
          "PUP for non-empty ScalarWave::Batched::ProjectionGroupMetadata is "
          "currently unsupported.");
    }
  }
};

struct FaceBoundaryMetadata {
  ::Kokkos::View<BoundaryCorrectionWorkItem*> no_projection_work_items{};
  ::Kokkos::View<size_t*> no_projection_oriented_remote_face_indices{};
  std::array<ProjectionGroupMetadata, boundary_projection_group_count>
      projection_groups{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    size_t no_projection_work_items_size = no_projection_work_items.extent(0);
    size_t no_projection_oriented_indices_size =
        no_projection_oriented_remote_face_indices.extent(0);
    p | no_projection_work_items_size;
    p | no_projection_oriented_indices_size;
    p | projection_groups;
    if (no_projection_work_items_size != 0 or
        no_projection_oriented_indices_size != 0) {
      ERROR(
          "PUP for non-empty ScalarWave::Batched::FaceBoundaryMetadata is "
          "currently unsupported.");
    }
  }
};

using OrientedRemoteFaceIndexMap =
    DirectionalIdMap<boundary_volume_dim, ::Kokkos::View<size_t*>>;
using MortarMetadataMap = DirectionalIdMap<boundary_volume_dim, MortarMetadata>;
using OrientedRemoteFaceIndexMaps = std::vector<OrientedRemoteFaceIndexMap>;
using MortarMetadataMaps = std::vector<MortarMetadataMap>;
using FaceBoundaryMetadataArray =
    std::array<FaceBoundaryMetadata, boundary_number_of_faces>;

}  // namespace ScalarWave::Batched
