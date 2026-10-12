// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace amr {
using CriterionTensorData = std::map<std::string, std::vector<DataVector>>;

/// A read-only snapshot for a prospective-parent criterion. The round number
/// distinguishes data that arrived before the local preparation action.
template <size_t Dim>
struct CriterionElementData {
  size_t round{0};
  Mesh<Dim> mesh{};
  CriterionTensorData tensors{};
  void pup(PUP::er& p) {
    p | round;
    p | mesh;
    p | tensors;
  }
};

/// Ephemeral sibling data, retained until AdjustDomain. Only the round is
/// inherited by new elements, so all elements rotate the join axis together.
template <size_t Dim>
struct CriterionNeighborData {
  size_t round{0};
  bool waiting{false};
  std::optional<std::pair<ElementId<Dim>, CriterionElementData<Dim>>> sibling{};
  void pup(PUP::er& p) {
    p | round;
    p | waiting;
    p | sibling;
  }
  void clear() {
    waiting = false;
    sibling.reset();
  }
};

/// Return a face neighbor that is an exact sibling along the selected axis.
/// Nonconforming siblings and neighbors across blocks are not candidates.
template <size_t Dim>
std::optional<ElementId<Dim>> criterion_sibling(const Element<Dim>& element,
                                                const size_t axis) {
  auto segments = element.id().segment_ids();
  if (segments[axis].refinement_level() == 0) {
    return std::nullopt;
  }
  const Direction<Dim> direction{axis, segments[axis].side_of_sibling()};
  segments[axis] = segments[axis].id_of_sibling();
  const ElementId<Dim> sibling{element.id().block_id(), segments,
                               element.id().grid_index()};
  const auto& neighbors = element.neighbors();
  const auto it = neighbors.find(direction);
  if (it != neighbors.end() and it->second.ids().size() == 1 and
      it->second.ids().count(sibling) == 1) {
    return sibling;
  }
  return std::nullopt;
}

namespace Tags {
template <size_t Dim>
struct CriterionNeighborData : db::SimpleTag {
  using type = amr::CriterionNeighborData<Dim>;
};
}  // namespace Tags
}  // namespace amr
