// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

#include "DataStructures/Index.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/SegmentSize.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct InitializeBoundaryBatchMetadata {
  static constexpr size_t volume_dim = 3;
  static constexpr size_t number_of_faces =
      ScalarWave::Batched::boundary_number_of_faces;
  static constexpr size_t projection_group_count =
      ScalarWave::Batched::boundary_projection_group_count;

  using return_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedBoundaryMetadata<ScalarWave::System<3>>,
      evolution::Kokkos::Tags::PackedBoundaryScratch<ScalarWave::System<3>>>;
  using argument_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>>;

  static constexpr size_t side_index(const Side side) {
    return side == Side::Upper ? static_cast<size_t>(1)
                               : static_cast<size_t>(0);
  }
  static constexpr size_t face_index(const size_t sliced_dim,
                                     const size_t side_i) {
    return 2 * sliced_dim + side_i;
  }
  static size_t face_index(const Direction<volume_dim>& direction) {
    return face_index(direction.dimension(), side_index(direction.side()));
  }

  static void apply(
      const gsl::not_null<evolution::Kokkos::Tags::PackedBoundaryMetadata<
          ScalarWave::System<3>>::type*>
          packed_boundary_metadata,
      const gsl::not_null<evolution::Kokkos::Tags::PackedBoundaryScratch<
          ScalarWave::System<3>>::type*>
          packed_boundary_scratch,
      const evolution::Kokkos::Tags::PackedTopology<
          ScalarWave::System<3>>::type& packed_topology) {
    using BoundaryCorrectionWorkItem =
        ScalarWave::Batched::BoundaryCorrectionWorkItem;
    using MortarMetadata = ScalarWave::Batched::MortarMetadata;
    using ProjectionGroupMetadata =
        ScalarWave::Batched::ProjectionGroupMetadata;
    using FaceBoundaryMetadata = ScalarWave::Batched::FaceBoundaryMetadata;
    using packed_boundary_scratch_type =
        evolution::Kokkos::PackedBoundaryScratch<ScalarWave::System<3>>;
    using device_package_field_tags =
        typename packed_boundary_scratch_type::device_package_field_tags;
    using device_dt_boundary_tags =
        typename packed_boundary_scratch_type::device_dt_boundary_tags;

    struct ProjectionGroupHostData {
      std::array<Spectral::SegmentSize, volume_dim - 1> mortar_size{
          {Spectral::SegmentSize::Uninitialized,
           Spectral::SegmentSize::Uninitialized}};
      Mesh<volume_dim - 1> mortar_mesh{};
      std::vector<BoundaryCorrectionWorkItem> work_items{};
      std::vector<size_t> oriented_remote_face_indices{};
    };
    struct FaceBoundaryHostData {
      std::vector<BoundaryCorrectionWorkItem> no_projection_work_items{};
      std::vector<size_t> no_projection_oriented_remote_face_indices{};
      std::array<ProjectionGroupHostData, projection_group_count>
          projection_groups{};
    };

    packed_boundary_metadata->oriented_remote_face_index_for_local.clear();
    packed_boundary_metadata->oriented_remote_face_index_for_local.resize(
        packed_topology.local_elements.size());
    packed_boundary_metadata->mortar_metadata.clear();
    packed_boundary_metadata->mortar_metadata.resize(
        packed_topology.local_elements.size());

    ASSERT(packed_topology.uniform_extents_host[0] > 0 and
               packed_topology.uniform_extents_host[1] > 0 and
               packed_topology.uniform_extents_host[2] > 0,
           "Cannot initialize boundary-batch metadata before uniform extents "
           "are known.");
    ASSERT(packed_topology.uniform_extents_host[0] ==
                   packed_topology.uniform_extents_host[1] and
               packed_topology.uniform_extents_host[1] ==
                   packed_topology.uniform_extents_host[2],
           "Boundary-batch metadata currently assumes uniform p across "
           "dimensions.");

    const Mesh<volume_dim> mesh{packed_topology.uniform_extents_host,
                                packed_topology.uniform_basis_host,
                                packed_topology.uniform_quadrature_host};
    const Index<volume_dim> extents_index{
        packed_topology.uniform_extents_host[0],
        packed_topology.uniform_extents_host[1],
        packed_topology.uniform_extents_host[2]};
    const size_t num_face_points = mesh.slice_away(0).number_of_grid_points();

    for (size_t e = 0; e < packed_topology.local_elements.size(); ++e) {
      const auto& element = packed_topology.local_elements[e];
      auto& oriented_remote_face_index_for_local =
          packed_boundary_metadata->oriented_remote_face_index_for_local[e];
      auto& mortar_metadata = packed_boundary_metadata->mortar_metadata[e];
      for (const auto& [direction, neighbors] : element.neighbors()) {
        const size_t sliced_dim = direction.dimension();
        const auto slice_extents = extents_index.slice_away(sliced_dim);
        const Mesh<volume_dim - 1> local_face_mesh =
            mesh.slice_away(sliced_dim);
        for (const auto& neighbor : neighbors) {
          const auto& orientation = neighbors.orientation(neighbor);
          std::vector<double> local_grid_point_indices(slice_extents.product());
          std::iota(local_grid_point_indices.begin(),
                    local_grid_point_indices.end(), 0.0);
          const std::vector<double> local_index_in_remote_order =
              orient_variables_on_slice(local_grid_point_indices, slice_extents,
                                        sliced_dim, orientation);

          std::vector<size_t> remote_index_for_local(slice_extents.product(),
                                                     0);
          std::vector<bool> seen_local_index(slice_extents.product(), false);
          for (size_t remote_index = 0;
               remote_index < local_index_in_remote_order.size();
               ++remote_index) {
            const size_t local_index =
                static_cast<size_t>(local_index_in_remote_order[remote_index]);
            ASSERT(
                local_index < remote_index_for_local.size(),
                "Oriented face-index map produced out-of-range local index.");
            ASSERT(not seen_local_index[local_index],
                   "Oriented face-index map produced duplicate local index.");
            remote_index_for_local[local_index] = remote_index;
            seen_local_index[local_index] = true;
          }
          for (size_t local_index = 0; local_index < seen_local_index.size();
               ++local_index) {
            ASSERT(seen_local_index[local_index],
                   "Oriented face-index map is missing local index "
                       << local_index << ".");
          }

          ::Kokkos::View<size_t*> device_map(
              "BatchedRemoteFaceIndexForLocalFaceIndex",
              remote_index_for_local.size());
          auto host_map = ::Kokkos::create_mirror_view(device_map);
          for (size_t i = 0; i < remote_index_for_local.size(); ++i) {
            host_map(i) = remote_index_for_local[i];
          }
          ::Kokkos::deep_copy(device_map, host_map);
          const DirectionalId<volume_dim> mortar_id{direction, neighbor};
          oriented_remote_face_index_for_local.insert_or_assign(
              mortar_id, std::move(device_map));

          const Direction<volume_dim> remote_direction =
              orientation(direction.opposite());
          const Mesh<volume_dim - 1> remote_face_mesh =
              mesh.slice_away(remote_direction.dimension());
          const Mesh<volume_dim - 1> oriented_remote_face_mesh =
              orient_mesh_on_slice(remote_face_mesh,
                                   remote_direction.dimension(), orientation);
          MortarMetadata mortar{};
          mortar.mortar_mesh =
              ::dg::mortar_mesh(local_face_mesh, oriented_remote_face_mesh);
          mortar.mortar_size = ::dg::mortar_size(element.id(), neighbor,
                                                 sliced_dim, orientation);
          mortar.needs_projection = Spectral::needs_projection(
              local_face_mesh, mortar.mortar_mesh, mortar.mortar_size);
          mortar_metadata.insert_or_assign(mortar_id, std::move(mortar));
        }
      }
    }

    std::array<FaceBoundaryHostData, number_of_faces> host_metadata{};
    for (size_t e = 0; e < packed_topology.local_elements.size(); ++e) {
      const auto& element = packed_topology.local_elements[e];
      const auto& oriented_remote_face_index_for_local =
          packed_boundary_metadata->oriented_remote_face_index_for_local[e];
      const auto& mortar_metadata =
          packed_boundary_metadata->mortar_metadata[e];
      for (const auto& [direction, neighbors] : element.neighbors()) {
        const size_t local_face_id = face_index(direction);
        auto& local_face_metadata = host_metadata[local_face_id];
        for (const auto& neighbor_id : neighbors) {
          const DirectionalId<volume_dim> mortar_id{direction, neighbor_id};
          ASSERT(oriented_remote_face_index_for_local.count(mortar_id) == 1,
                 "Missing oriented remote-face index map for " << mortar_id
                                                               << ".");
          ASSERT(mortar_metadata.count(mortar_id) == 1,
                 "Missing mortar metadata for " << mortar_id << ".");
          const auto& mortar = mortar_metadata.at(mortar_id);

          const auto& orientation = neighbors.orientation(neighbor_id);
          const Direction<volume_dim> remote_direction =
              orientation(direction.opposite());
          const auto remote_face_index_for_local_face_index =
              oriented_remote_face_index_for_local.at(mortar_id);
          ASSERT(remote_face_index_for_local_face_index.extent(0) ==
                     num_face_points,
                 "Remote face-orientation map size mismatch for mortar "
                     << mortar_id << ".");
          const auto host_remote_face_index_for_local_face_index =
              ::Kokkos::create_mirror_view_and_copy(
                  ::Kokkos::HostSpace{},
                  remote_face_index_for_local_face_index);

          BoundaryCorrectionWorkItem work_item{
              e, packed_topology.element_index_by_id.at(neighbor_id),
              face_index(remote_direction), 0};
          if (mortar.needs_projection) {
            const size_t projection_group_key =
                static_cast<size_t>(gsl::at(mortar.mortar_size, 0)) +
                4 * static_cast<size_t>(gsl::at(mortar.mortar_size, 1));
            ASSERT(projection_group_key < projection_group_count,
                   "Projection-group key out of range.");
            auto& projection_group =
                local_face_metadata.projection_groups[projection_group_key];
            if (projection_group.work_items.empty()) {
              projection_group.mortar_size = mortar.mortar_size;
              projection_group.mortar_mesh = mortar.mortar_mesh;
            } else {
              ASSERT(projection_group.mortar_size == mortar.mortar_size,
                     "Mixed mortar sizes in projection group.");
              ASSERT(projection_group.mortar_mesh.extents() ==
                         mortar.mortar_mesh.extents(),
                     "Mixed mortar extents in projection group.");
            }
            work_item.oriented_remote_face_index_offset =
                projection_group.oriented_remote_face_indices.size();
            projection_group.oriented_remote_face_indices.reserve(
                projection_group.oriented_remote_face_indices.size() +
                num_face_points);
            for (size_t i = 0; i < num_face_points; ++i) {
              projection_group.oriented_remote_face_indices.push_back(
                  host_remote_face_index_for_local_face_index(i));
            }
            projection_group.work_items.push_back(work_item);
          } else {
            work_item.oriented_remote_face_index_offset =
                local_face_metadata.no_projection_oriented_remote_face_indices
                    .size();
            local_face_metadata.no_projection_oriented_remote_face_indices
                .reserve(
                    local_face_metadata
                        .no_projection_oriented_remote_face_indices.size() +
                    num_face_points);
            for (size_t i = 0; i < num_face_points; ++i) {
              local_face_metadata.no_projection_oriented_remote_face_indices
                  .push_back(host_remote_face_index_for_local_face_index(i));
            }
            local_face_metadata.no_projection_work_items.push_back(work_item);
          }
        }
      }
    }

    for (size_t face_id = 0; face_id < number_of_faces; ++face_id) {
      auto& face_metadata =
          packed_boundary_metadata->boundary_correction_face_metadata[face_id];
      const auto& host_face_metadata = host_metadata[face_id];
      face_metadata.no_projection_work_items =
          copy_work_items_to_device(host_face_metadata.no_projection_work_items,
                                    "BatchedNoProjectionMortarWorkItems");
      face_metadata.no_projection_oriented_remote_face_indices =
          copy_indices_to_device(
              host_face_metadata.no_projection_oriented_remote_face_indices,
              "BatchedNoProjectionOrientedRemoteFaceIndices");
      for (size_t projection_group_key = 0;
           projection_group_key < projection_group_count;
           ++projection_group_key) {
        auto& projection_group =
            face_metadata.projection_groups[projection_group_key];
        const auto& host_projection_group =
            host_face_metadata.projection_groups[projection_group_key];
        projection_group.mortar_size = host_projection_group.mortar_size;
        projection_group.mortar_mesh = host_projection_group.mortar_mesh;
        projection_group.work_items = copy_work_items_to_device(
            host_projection_group.work_items, "BatchedProjectionWorkItems");
        projection_group.oriented_remote_face_indices = copy_indices_to_device(
            host_projection_group.oriented_remote_face_indices,
            "BatchedProjectionOrientedRemoteFaceIndices");
      }
    }

    for (auto& packaged_face_data :
         packed_boundary_scratch->packaged_face_data_for_all_elements) {
      packaged_face_data = Variables<device_package_field_tags>{
          packed_topology.local_elements.size() * num_face_points};
      ::Kokkos::deep_copy(packaged_face_data.view(), 0.0);
    }
    for (auto& internal_boundary_terms :
         packed_boundary_scratch->internal_boundary_terms_for_all_elements) {
      internal_boundary_terms = Variables<device_dt_boundary_tags>{
          packed_topology.local_elements.size() * num_face_points};
      ::Kokkos::deep_copy(internal_boundary_terms.view(), 0.0);
    }
  }

 private:
  static ::Kokkos::View<ScalarWave::Batched::BoundaryCorrectionWorkItem*>
  copy_work_items_to_device(
      const std::vector<ScalarWave::Batched::BoundaryCorrectionWorkItem>&
          host_data,
      const char* label) {
    ::Kokkos::View<ScalarWave::Batched::BoundaryCorrectionWorkItem*>
        device_data{label, host_data.size()};
    auto host_view = ::Kokkos::create_mirror_view(device_data);
    for (size_t i = 0; i < host_data.size(); ++i) {
      host_view(i) = host_data[i];
    }
    ::Kokkos::deep_copy(device_data, host_view);
    return device_data;
  }

  static ::Kokkos::View<size_t*> copy_indices_to_device(
      const std::vector<size_t>& host_data, const char* label) {
    ::Kokkos::View<size_t*> device_data{label, host_data.size()};
    auto host_view = ::Kokkos::create_mirror_view(device_data);
    for (size_t i = 0; i < host_data.size(); ++i) {
      host_view(i) = host_data[i];
    }
    ::Kokkos::deep_copy(device_data, host_view);
    return device_data;
  }
};

}  // namespace ScalarWave::Actions
