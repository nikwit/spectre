// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Amr/NeighborsOfChild.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Helpers.hpp"
#include "Domain/Amr/Info.hpp"
#include "Domain/Amr/NewNeighborIds.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/FaceType.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace amr {
template <size_t VolumeDim>
std::pair<DirectionMap<VolumeDim, Neighbors<VolumeDim>>,
          DirectionalIdMap<VolumeDim, Mesh<VolumeDim>>>
neighbors_of_child(
    const Element<VolumeDim>& parent, const Info<VolumeDim>& parent_info,
    const std::unordered_map<ElementId<VolumeDim>, Info<VolumeDim>>&
        parent_neighbor_info,
    const ElementId<VolumeDim>& child_id) {
  std::pair<DirectionMap<VolumeDim, Neighbors<VolumeDim>>,
            DirectionalIdMap<VolumeDim, Mesh<VolumeDim>>>
      result;

  const auto sibling_id = [&child_id](const size_t dim) {
    auto sibling_segment_ids = child_id.segment_ids();
    auto& segment_id_to_change = gsl::at(sibling_segment_ids, dim);
    segment_id_to_change = segment_id_to_change.id_of_sibling();
    return ElementId<VolumeDim>{child_id.block_id(), sibling_segment_ids,
                                child_id.grid_index()};
  };

  for (const auto& [direction, old_neighbors] : parent.neighbors()) {
    if (not old_neighbors.are_conforming()) {
      // Nonconforming block coordinates are not related by an orientation
      // map. For p changes keep the connectivity, and store a single
      // neighbor's mesh in its own frame, as in DgDomain initialization.
      // The many-neighbor side intentionally stores no neighbor meshes.
      const auto changes_h = [](const auto& flags) {
        return std::any_of(flags.begin(), flags.end(), [](const auto flag) {
          return flag == Flag::Split or flag == Flag::Join;
        });
      };
      if (child_id != parent.id() or changes_h(parent_info.flags)) {
        ERROR(
            "AMR across nonconforming block interfaces currently supports "
            "p refinement only");
      }
      for (const auto& neighbor : old_neighbors) {
        const auto& info = parent_neighbor_info.at(neighbor);
        if (changes_h(info.flags)) {
          ERROR(
              "AMR across nonconforming block interfaces currently supports "
              "p refinement only");
        }
        if (parent.face_types().at(direction) ==
            domain::FaceType::SingleNonconforming) {
          result.second.emplace(DirectionalId{direction, neighbor},
                                info.new_mesh);
        }
      }
      result.first.emplace(direction, old_neighbors);
      continue;
    }
    const auto dim = direction.dimension();
    if (gsl::at(parent_info.flags, dim) == Flag::Split and
        has_potential_sibling(child_id, direction)) {
      const auto id = sibling_id(dim);
      result.first.emplace(
          direction, Neighbors<VolumeDim>{
                         {id}, OrientationMap<VolumeDim>::create_aligned()});
      result.second.insert({{direction, id}, parent_info.new_mesh});
    } else {
      const auto new_neighbor_ids_and_meshes = amr::new_neighbor_ids(
          child_id, direction, old_neighbors, parent_neighbor_info);
      std::unordered_set<ElementId<VolumeDim>> new_neighbor_ids;
      for (const auto& [id, mesh] : new_neighbor_ids_and_meshes) {
        result.second.insert_or_assign({direction, id}, mesh);
        new_neighbor_ids.insert(id);
      }

      result.first.emplace(
          direction,
          Neighbors<VolumeDim>{new_neighbor_ids, old_neighbors.orientations(),
                               old_neighbors.are_conforming()});
    }
  }

  for (const auto& direction : parent.external_boundaries()) {
    const auto dim = direction.dimension();
    if (gsl::at(parent_info.flags, dim) == Flag::Split and
        has_potential_sibling(child_id, direction)) {
      const auto id = sibling_id(dim);
      result.first.emplace(
          direction, Neighbors<VolumeDim>{
                         {id}, OrientationMap<VolumeDim>::create_aligned()});
      result.second.insert({{direction, id}, parent_info.new_mesh});
    }
  }

  return result;
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                \
  template std::pair<DirectionMap<DIM(data), Neighbors<DIM(data)>>,         \
                     DirectionalIdMap<DIM(data), Mesh<DIM(data)>>>          \
  neighbors_of_child(                                                       \
      const Element<DIM(data)>& parent, const Info<DIM(data)>& parent_info, \
      const std::unordered_map<ElementId<DIM(data)>, Info<DIM(data)>>&      \
          parent_neighbor_info,                                             \
      const ElementId<DIM(data)>& child_id);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef DIM
#undef INSTANTIATE
}  // namespace amr
