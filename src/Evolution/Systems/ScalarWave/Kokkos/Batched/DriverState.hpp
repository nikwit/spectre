// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "DataStructures/Index.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/InitialExtents.hpp"
#include "Domain/Creators/Tags/InitialRefinementLevels.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Batched::Initialization {

struct DriverState {
  using packed_system = ScalarWave::System<3>;
  using packed_topology_tag =
      evolution::Kokkos::Tags::PackedTopology<packed_system>;
  using packed_geometry_tag =
      evolution::Kokkos::Tags::PackedGeometry<packed_system>;
  using packed_boundary_metadata_tag =
      evolution::Kokkos::Tags::PackedBoundaryMetadata<packed_system>;
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<packed_system>;
  using packed_boundary_scratch_tag =
      evolution::Kokkos::Tags::PackedBoundaryScratch<packed_system>;
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;

  using const_global_cache_tags = tmpl::list<::domain::Tags::Domain<3>>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options =
      tmpl::list<domain::Tags::InitialRefinementLevels<3>,
                 domain::Tags::InitialExtents<3>,
                 evolution::dg::Tags::Quadrature>;
  using simple_tags =
      tmpl::list<packed_topology_tag, packed_geometry_tag,
                 packed_boundary_metadata_tag, packed_evolution_state_tag,
                 packed_boundary_scratch_tag, device_constraint_gamma2_tag>;
  using compute_tags = tmpl::list<>;

  using return_tags = tmpl::list<packed_topology_tag>;
  using argument_tags = tmpl::list<
      ::domain::Tags::Domain<3>, ::domain::Tags::InitialRefinementLevels<3>,
      ::domain::Tags::InitialExtents<3>, evolution::dg::Tags::Quadrature>;

  static void apply(
      const gsl::not_null<typename packed_topology_tag::type*> packed_topology,
      const ::Domain<3>& domain,
      const std::vector<std::array<size_t, 3>>& initial_refinement_levels,
      const std::vector<std::array<size_t, 3>>& initial_extents,
      const Spectral::Quadrature& quadrature) {
    constexpr size_t volume_dim = 3;
    ASSERT(initial_refinement_levels.size() == domain.blocks().size(),
           "Initial refinement levels must contain one entry per block. Have "
               << initial_refinement_levels.size() << " levels for "
               << domain.blocks().size() << " blocks.");
    ASSERT(initial_extents.size() == domain.blocks().size(),
           "Initial extents must contain one entry per block. Have "
               << initial_extents.size() << " extents for "
               << domain.blocks().size() << " blocks.");

    packed_topology->local_element_ids.clear();
    packed_topology->local_elements.clear();
    packed_topology->element_index_by_id.clear();
    packed_topology->element_extents_host.clear();
    packed_topology->element_point_offsets_host.clear();
    packed_topology->element_point_offsets_host.push_back(0);
    packed_topology->points_per_element = 0;
    std::array<size_t, volume_dim> first_extents{};
    bool have_first_extents = false;

    for (const auto& block : domain.blocks()) {
      const size_t block_id = block.id();
      const auto block_element_ids = initial_element_ids(
          block_id, gsl::at(initial_refinement_levels, block_id));
      const auto& block_extents = gsl::at(initial_extents, block_id);
      if (not have_first_extents) {
        first_extents = block_extents;
        have_first_extents = true;
      } else {
        ASSERT(
            block_extents == first_extents,
            "EvolveScalarWaveKokkosBatched currently requires uniform extents "
            "across all blocks/elements. Found extents mismatch: first="
                << first_extents[0] << "," << first_extents[1] << ","
                << first_extents[2] << " current=" << block_extents[0] << ","
                << block_extents[1] << "," << block_extents[2] << ".");
      }
      const size_t block_points_per_element =
          block_extents[0] * block_extents[1] * block_extents[2];
      if (packed_topology->points_per_element == 0) {
        packed_topology->points_per_element = block_points_per_element;
      } else {
        ASSERT(packed_topology->points_per_element == block_points_per_element,
               "EvolveScalarWaveKokkosBatched currently requires uniform "
               "points per element.");
      }

      for (const auto& element_id : block_element_ids) {
        const auto element = ::domain::Initialization::create_initial_element(
            element_id, domain.blocks(), initial_refinement_levels);
        packed_topology->element_index_by_id.insert_or_assign(
            element_id, packed_topology->local_element_ids.size());
        packed_topology->local_element_ids.push_back(element_id);
        packed_topology->local_elements.push_back(element);
        packed_topology->element_extents_host.push_back(block_extents);
        packed_topology->element_point_offsets_host.push_back(
            packed_topology->element_point_offsets_host.back() +
            packed_topology->points_per_element);
      }
    }
    if (have_first_extents) {
      packed_topology->uniform_extents_host = first_extents;
    }
    packed_topology->uniform_basis_host = {{Spectral::Basis::Legendre,
                                            Spectral::Basis::Legendre,
                                            Spectral::Basis::Legendre}};
    packed_topology->uniform_quadrature_host = {quadrature, quadrature,
                                                quadrature};

    packed_topology->total_points =
        packed_topology->element_point_offsets_host.empty()
            ? 0
            : packed_topology->element_point_offsets_host.back();

    packed_topology->element_point_offsets_device = ::Kokkos::View<size_t*>(
        "BatchedElementPointOffsets",
        packed_topology->element_point_offsets_host.size());
    auto host_point_offsets = ::Kokkos::create_mirror_view(
        packed_topology->element_point_offsets_device);
    for (size_t i = 0; i < packed_topology->element_point_offsets_host.size();
         ++i) {
      host_point_offsets(i) = packed_topology->element_point_offsets_host[i];
    }
    ::Kokkos::deep_copy(packed_topology->element_point_offsets_device,
                        host_point_offsets);

    packed_topology->element_extents_device = ::Kokkos::View<size_t* [3]>(
        "BatchedElementExtents", packed_topology->element_extents_host.size());
    auto host_extents =
        ::Kokkos::create_mirror_view(packed_topology->element_extents_device);
    for (size_t e = 0; e < packed_topology->element_extents_host.size(); ++e) {
      for (size_t d = 0; d < volume_dim; ++d) {
        host_extents(e, d) = packed_topology->element_extents_host[e][d];
      }
    }
    ::Kokkos::deep_copy(packed_topology->element_extents_device, host_extents);

    ASSERT(packed_topology->uniform_extents_host[0] > 0 and
               packed_topology->uniform_extents_host[1] > 0 and
               packed_topology->uniform_extents_host[2] > 0,
           "Cannot initialize batched face-to-volume index maps before uniform "
           "extents are known.");
    const Index<3> extents_index{packed_topology->uniform_extents_host[0],
                                 packed_topology->uniform_extents_host[1],
                                 packed_topology->uniform_extents_host[2]};
    const auto [indices_buffer, volume_and_slice_index_map] =
        volume_and_slice_indices(extents_index);
    (void)indices_buffer;
    for (size_t d = 0; d < volume_dim; ++d) {
      const size_t num_face_points = extents_index.slice_away(d).product();
      auto& lower_face_to_volume_index =
          gsl::at(packed_topology->device_face_to_volume_index_map, d).first;
      auto& upper_face_to_volume_index =
          gsl::at(packed_topology->device_face_to_volume_index_map, d).second;
      lower_face_to_volume_index = ::Kokkos::View<size_t*>(
          "BatchedLowerFaceToVolumeIndexMap", num_face_points);
      upper_face_to_volume_index = ::Kokkos::View<size_t*>(
          "BatchedUpperFaceToVolumeIndexMap", num_face_points);
      auto host_lower_face_to_volume_index =
          ::Kokkos::create_mirror_view(lower_face_to_volume_index);
      auto host_upper_face_to_volume_index =
          ::Kokkos::create_mirror_view(upper_face_to_volume_index);
      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).first) {
        host_lower_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }
      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).second) {
        host_upper_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }
      ::Kokkos::deep_copy(lower_face_to_volume_index,
                          host_lower_face_to_volume_index);
      ::Kokkos::deep_copy(upper_face_to_volume_index,
                          host_upper_face_to_volume_index);
    }
  }
};

}  // namespace ScalarWave::Batched::Initialization
