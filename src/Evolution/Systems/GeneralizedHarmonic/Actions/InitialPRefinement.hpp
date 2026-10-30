// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Tags/Flags.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Bbh/PhaseControl/InitialPRefinement.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/PhaseControl/ContributeToPhaseChangeReduction.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Tags/Criteria.hpp"

namespace gh::Actions {
namespace detail {
struct MeshBeforeInitialPRefinement : db::SimpleTag {
  using type = Mesh<3>;
};
}  // namespace detail
/// Enter the startup-only decision mode after initial data and gauge
/// constraints have been set, including after each resampling iteration.
struct PrepareInitialPRefinement {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename Component>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, const tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex&,
      ActionList, const Component*) {
    if (bbh::phase_control::initial_p_refinement_enabled(cache)) {
      const auto& criteria = db::get<amr::Criteria::Tags::Criteria>(box);
      if (criteria.size() != 1 or dynamic_cast<const SpecSpectral*>(
                                      criteria.front().get()) == nullptr) {
        ERROR(
            "InitialPRefinement requires SpecSpectral as the sole AMR "
            "criterion.");
      }
      db::mutate<spec_amr::HistoryTag>(
          [](const auto history) { history->initializing = true; },
          make_not_null(&box));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

/// Run only after refinement criteria and neighbor decisions have quiesced.
/// Reject h changes, report local policy limits, and save the pre-adjustment
/// mesh. CheckInitialPRefinementResult counts the finalized mesh changes.
struct CheckInitialPRefinement {
  using simple_tags = tmpl::list<detail::MeshBeforeInitialPRefinement>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename Component>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, const tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& id,
      ActionList, const Component*) {
    const auto& history = db::get<spec_amr::HistoryTag>(box);
    const auto& flags = db::get<amr::Tags::Info<3>>(box).flags;
    bool changed = false;
    size_t limited = 0;
    for (size_t d = 0; d < 3; ++d) {
      if (flags[d] == amr::Flag::Split or flags[d] == amr::Flag::Join) {
        ERROR("InitialPRefinement cannot change the h mesh.");
      }
      changed |= flags[d] == amr::Flag::IncreaseResolution or
                 flags[d] == amr::Flag::DecreaseResolution;
      limited += history.decision[d] == amr::Flag::IncreaseResolution and
                 flags[d] == amr::Flag::DoNothing;
    }
    const auto& current_mesh = db::get<domain::Tags::Mesh<3>>(box);
    db::mutate<detail::MeshBeforeInitialPRefinement>(
        [&current_mesh](const auto mesh) { *mesh = current_mesh; },
        make_not_null(&box));
    using namespace bbh::phase_control::initial_p_refinement;
    Parallel::contribute_to_phase_change_reduction<Component>(
        tuples::TaggedTuple<ChangedElements, LimitedDirections,
                            CyclingElements>{
            size_t{changed}, limited, size_t{history.initial_cycle_detected}},
        cache, id);
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

/// Count actual mesh changes after central policies and projection have
/// completed. Local proposals may be cancelled or enlarged by an object floor.
struct CheckInitialPRefinementResult {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename Component>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, const tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& id,
      ActionList, const Component*) {
    if (db::get<spec_amr::HistoryTag>(box).initializing) {
      const bool changed = db::get<domain::Tags::Mesh<3>>(box) !=
                           db::get<detail::MeshBeforeInitialPRefinement>(box);
      using bbh::phase_control::initial_p_refinement::ChangedElements;
      Parallel::contribute_to_phase_change_reduction<Component>(
          tuples::TaggedTuple<ChangedElements>{size_t{changed}}, cache, id);
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

/// Leave startup mode before registering observers and initializing the time
/// stepper. No evolution call count, prediction or cooldown was accumulated.
struct FinishInitialPRefinement {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename Component>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, const tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Metavariables>&, const ArrayIndex&, ActionList,
      const Component*) {
    db::mutate<spec_amr::HistoryTag>(
        [](const auto history) {
          history->initializing = false;
          history->initial_meshes.clear();
          history->initial_cycle_detected = false;
        },
        make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace gh::Actions
