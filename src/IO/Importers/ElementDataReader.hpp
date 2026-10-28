// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Importers/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Algorithms/AlgorithmNodegroup.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/Gsl.hpp"

namespace importers {

namespace detail {
template <size_t Dim>
struct InitializeElementDataReader;
struct ResetCompletedImports;
}  // namespace detail

/*!
 * \brief A nodegroup parallel component that reads in a volume data file and
 * distributes its data to elements of an array parallel component.
 *
 * Each element of the array parallel component must register itself before
 * data can be sent to it. To do so, invoke
 * `importers::Actions::RegisterWithElementDataReader` on each element. In a
 * subsequent phase you can then invoke
 * `importers::ThreadedActions::ReadVolumeData` on the `ElementDataReader`
 * component to read in the file and distribute its data to the registered
 * elements.
 *
 * \see Dev guide on \ref dev_guide_importing
 */
template <typename Metavariables>
struct ElementDataReader {
  static constexpr size_t Dim = Metavariables::volume_dim;

  using chare_type = Parallel::Algorithms::Nodegroup;
  static constexpr bool checkpoint_data = false;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<detail::InitializeElementDataReader<Dim>>>,
      Parallel::PhaseActions<Parallel::Phase::RegisterWithElementDataReader,
                             tmpl::list<detail::ResetCompletedImports>>>;
  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void initialize(
      Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<ElementDataReader>(local_cache)
        .start_phase(next_phase);
  }
};

namespace detail {
/// Registration may update coordinates after initial p refinement. Permit the
/// next import to resample even when its file/observation identifier is
/// unchanged. The phase boundary ensures all registration updates finish before
/// reading.
struct ResetCompletedImports {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box, const tuples::TaggedTuple<InboxTags...>&,
      const Parallel::GlobalCache<Metavariables>&, const ArrayIndex&,
      ActionList, const ParallelComponent*) {
    db::mutate<Tags::ElementDataAlreadyRead>(
        [](const auto already_read) { already_read->clear(); },
        make_not_null(&box));
    return {Parallel::AlgorithmExecution::Pause, 0};
  }
};

template <size_t Dim>
struct InitializeElementDataReader {
  using simple_tags =
      tmpl::list<Tags::RegisteredElements<Dim>, Tags::ElementDataAlreadyRead>;
  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& /*box*/,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    return {Parallel::AlgorithmExecution::Pause, std::nullopt};
  }
};
}  // namespace detail

}  // namespace importers
