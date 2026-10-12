// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/ArrayCollection/IsDgElementCollection.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Amr/Actions/EvaluateRefinementCriteria.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Tags/Criteria.hpp"
#include "ParallelAlgorithms/Amr/CriterionData.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace amr::Actions {
/// Receive a sibling snapshot. Preparation and receipt may occur in either
/// order. Evaluate only after both are complete; no element mutates its fields
/// until the subsequent AdjustDomain phase after Charm++ quiescence.
struct ReceiveCriterionData {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ElementId<Metavariables::volume_dim>& element_id,
                    const ElementId<Metavariables::volume_dim>& sender,
                    CriterionElementData<Metavariables::volume_dim> data) {
    constexpr size_t Dim = Metavariables::volume_dim;
    auto& state = db::get_mutable_reference<Tags::CriterionNeighborData<Dim>>(
        make_not_null(&box));
    ASSERT(not state.sibling.has_value(), "Duplicate sibling snapshot");
    ASSERT(data.round == state.round or data.round == state.round + 1,
           "Stale sibling snapshot");
    const auto sibling = criterion_sibling(
        db::get<domain::Tags::Element<Dim>>(box), (data.round - 1) % Dim);
    ASSERT(sibling.has_value() and *sibling == sender,
           "Snapshot is not from the selected face sibling");
    state.sibling.emplace(sender, std::move(data));
    if (state.waiting) {
      ASSERT(state.sibling->second.round == state.round,
             "Sibling AMR rounds differ");
      state.waiting = false;
      EvaluateRefinementCriteria::apply<ParallelComponent>(box, cache,
                                                           element_id);
    }
  }
};

/// Local criteria retain their original path. Criteria that request sibling
/// data exchange nodal tensors along one rotating axis per AMR round. Only
/// exact face siblings participate, so there are at most two children in a
/// certified join. Data are cleared by AdjustDomain and the round counter is
/// inherited through both splitting and joining.
struct PrepareRefinementCriteria {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ElementId<Metavariables::volume_dim>& element_id) {
    constexpr size_t Dim = Metavariables::volume_dim;
    const auto& criteria = db::get<amr::Criteria::Tags::Criteria>(box);
    const bool needs_data = alg::any_of(criteria, [](const auto& criterion) {
      return criterion->needs_sibling_data();
    });
    if (not needs_data) {
      EvaluateRefinementCriteria::apply<ParallelComponent>(box, cache,
                                                           element_id);
      return;
    }
    if constexpr (not db::tag_is_retrievable_v<Tags::CriterionNeighborData<Dim>,
                                               db::DataBox<DbTags>> or
                  Parallel::is_dg_element_collection_v<ParallelComponent> or
                  Metavariables::amr::keep_coarse_grids) {
      ERROR(
          "Sibling criteria require an AMR element array without retained "
          "coarse grids, initialized by amr::Initialization::Initialize.");
    } else {
      auto& state = db::get_mutable_reference<Tags::CriterionNeighborData<Dim>>(
          make_not_null(&box));
      ASSERT(not state.waiting, "AMR preparation repeated before evaluation");
      ++state.round;
      const auto sibling = criterion_sibling(
          db::get<domain::Tags::Element<Dim>>(box), (state.round - 1) % Dim);
      if (not sibling.has_value()) {
        ASSERT(not state.sibling.has_value(), "Unexpected sibling snapshot");
        EvaluateRefinementCriteria::apply<ParallelComponent>(box, cache,
                                                             element_id);
        return;
      }
      CriterionElementData<Dim> data{};
      data.round = state.round;
      data.mesh = db::get<domain::Tags::Mesh<Dim>>(box);
      for (const auto& criterion : criteria) {
        if (criterion->needs_sibling_data()) {
          auto tensors = criterion->template sibling_data<Metavariables>(box);
          data.tensors.merge(tensors);
        }
      }
      auto& array = Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::simple_action<ReceiveCriterionData>(array[*sibling], element_id,
                                                    std::move(data));
      if (state.sibling.has_value()) {
        ASSERT(state.sibling->second.round == state.round,
               "Sibling AMR rounds differ");
        EvaluateRefinementCriteria::apply<ParallelComponent>(box, cache,
                                                             element_id);
      } else {
        state.waiting = true;
      }
    }
  }
};
}  // namespace amr::Actions
