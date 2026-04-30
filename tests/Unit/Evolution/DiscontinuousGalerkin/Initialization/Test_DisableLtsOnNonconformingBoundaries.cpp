// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/DiscontinuousGalerkin/Initialization/DisableLtsOnNonconformingBoundaries.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Block.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Structure/BlockNeighbors.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/Structure/Topology.hpp"
#include "Domain/Tags.hpp"
#include "Framework/TestingFramework.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/FixedLtsRatio.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"

namespace {
std::vector<Block<2>> make_blocks() {
  const auto aligned = OrientationMap<2>::create_aligned();
  std::vector<Block<2>> blocks{};
  blocks.emplace_back(
      nullptr, 0,
      DirectionMap<2, BlockNeighbors<2>>{
          {Direction<2>::upper_xi(),
           BlockNeighbors<2>{{1, 2}, {{1, aligned}, {2, aligned}}, false}}},
      "Annulus", std::array{domain::Topology::I1, domain::Topology::I1});
  blocks.emplace_back(
      nullptr, 1,
      DirectionMap<2, BlockNeighbors<2>>{
          {Direction<2>::lower_xi(),
           BlockNeighbors<2>{{0}, {{0, aligned}}, false}},
          {Direction<2>::upper_eta(), BlockNeighbors<2>{2, aligned}}},
      "Wedge0", std::array{domain::Topology::I1, domain::Topology::I1});
  blocks.emplace_back(
      nullptr, 2,
      DirectionMap<2, BlockNeighbors<2>>{
          {Direction<2>::lower_xi(),
           BlockNeighbors<2>{{0}, {{0, aligned}}, false}},
          {Direction<2>::lower_eta(), BlockNeighbors<2>{1, aligned}}},
      "Wedge1", std::array{domain::Topology::I1, domain::Topology::I1});
  blocks.emplace_back(nullptr, 3, DirectionMap<2, BlockNeighbors<2>>{},
                      "Isolated",
                      std::array{domain::Topology::I1, domain::Topology::I1});
  return blocks;
}

void check_element(const Element<2>& element,
                   const std::optional<size_t>& expected_ratio) {
  const Slab slab{2.0, 10.0};
  auto box = db::create<db::AddSimpleTags<
      Tags::FixedLtsRatio, domain::Tags::Element<2>, Tags::TimeStep>>(
      std::optional<size_t>{}, element, slab.duration() / 8);

  db::mutate_apply<
      evolution::dg::Initialization::DisableLtsOnNonconformingBoundaries<2>>(
      make_not_null(&box));

  CHECK(db::get<Tags::FixedLtsRatio>(box) == expected_ratio);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.DG.Initialization.DisableLtsOnNonconformingBoundaries",
    "[Evolution][Unit]") {
  const auto blocks = make_blocks();
  const std::vector<std::array<size_t, 2>> initial_refinement_levels{
      std::array{1_st, 0_st}, std::array{0_st, 1_st}, std::array{0_st, 1_st},
      std::array{0_st, 0_st}};

  check_element(
      domain::create_initial_element(
          ElementId<2>{0, std::array{SegmentId{1, 1}, SegmentId{0, 0}}},
          blocks, initial_refinement_levels),
      8_st);
  check_element(
      domain::create_initial_element(
          ElementId<2>{1, std::array{SegmentId{0, 0}, SegmentId{1, 0}}},
          blocks, initial_refinement_levels),
      8_st);
  check_element(domain::create_initial_element(ElementId<2>{3}, blocks,
                                              initial_refinement_levels),
                std::nullopt);
}
