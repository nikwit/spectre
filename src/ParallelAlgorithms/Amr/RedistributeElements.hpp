// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementDistribution.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Domain/Tags/ElementDistribution.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Algorithms/AlgorithmSingletonDeclarations.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Parallel/Tags/Section.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace amr::redistribution {
template <typename Metavariables>
struct Component;

namespace Tags {
struct TargetProc : db::SimpleTag {
  using type = std::optional<int>;
};
template <size_t Dim>
struct Costs : db::SimpleTag {
  // Current cost and processor. All records arrive in
  // CollectElementDistribution before any are used in PlanElementDistribution.
  using type = std::unordered_map<ElementId<Dim>, std::pair<double, int>>;
  using option_tags = tmpl::list<>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options() { return {}; }
};
}  // namespace Tags

/// Reuse the initial distribution algorithm on a complete fixed h mesh with
/// arbitrary per-element p. Reject nonuniform h and multigrid layouts before
/// sending any migration commands. Ignored processors remain unavailable.
template <size_t Dim>
std::map<ElementId<Dim>, int> make_plan(
    const std::vector<Block<Dim>>& blocks,
    const typename Tags::Costs<Dim>::type& records, const size_t num_procs,
    const std::unordered_set<size_t>& ignored, const bool weighted) {
  if (blocks.empty() or records.empty()) {
    ERROR("Redistribution requires a nonempty domain and element set.");
  }
  if (ignored.size() >= num_procs) {
    ERROR("Redistribution requires at least one available processor.");
  }
  for (const auto proc : ignored) {
    if (proc >= num_procs) {
      ERROR("Invalid ignored processor in redistribution.");
    }
  }
  std::vector<std::optional<std::array<size_t, Dim>>> levels(blocks.size());
  std::unordered_map<ElementId<Dim>, double> costs{};
  for (const auto& [id, record] : records) {
    if (id.block_id() >= blocks.size() or id.grid_index() != 0 or
        not std::isfinite(record.first) or record.first <= 0. or
        record.second < 0 or static_cast<size_t>(record.second) >= num_procs) {
      ERROR("Invalid element or cost in redistribution: " << id);
    }
    auto& block_levels = levels[id.block_id()];
    if (block_levels and *block_levels != id.refinement_levels()) {
      ERROR("Redistribution currently requires uniform h within each block.");
    }
    block_levels = id.refinement_levels();
    costs.emplace(id, record.first);
  }
  std::vector<std::array<size_t, Dim>> refinement{};
  std::vector<ElementId<Dim>> ordered_ids{};
  for (size_t b = 0; b < blocks.size(); ++b) {
    if (not levels[b]) {
      ERROR("Missing block in redistribution: " << b);
    }
    refinement.push_back(*levels[b]);
    for (const auto& id : initial_element_ids(b, *levels[b])) {
      if (records.count(id) != 1) {
        ERROR("Missing element in redistribution: " << id);
      }
      ordered_ids.push_back(id);
    }
  }
  if (ordered_ids.size() != records.size()) {
    ERROR("Unexpected elements in redistribution.");
  }
  std::map<ElementId<Dim>, int> result{};
  if (weighted) {
    const domain::BlockZCurveProcDistribution<Dim> distribution{
        costs, num_procs - ignored.size(), blocks, refinement, ignored};
    for (const auto& id : ordered_ids) {
      result.emplace(id, distribution.get_proc_for_element(id));
    }
  } else {
    size_t proc = 0;
    for (const auto& id : ordered_ids) {
      while (ignored.count(proc) != 0) {
        proc = (proc + 1) % num_procs;
      }
      result.emplace(id, proc);
      proc = (proc + 1) % num_procs;
    }
  }
  return result;
}

namespace Actions {
struct SetTargetProc {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*id*/, const int target) {
    db::mutate<Tags::TargetProc>(
        [target](const auto destination) {
          if (destination->has_value()) {
            ERROR("Duplicate redistribution destination.");
          }
          *destination = target;
        },
        make_not_null(&box));
  }
};

struct RecordCost {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*index*/,
                    const ElementId<Metavariables::volume_dim>& id,
                    const double cost, const int proc) {
    db::mutate<Tags::Costs<Metavariables::volume_dim>>(
        [&](const auto records) {
          if (not records->emplace(id, std::pair{cost, proc}).second) {
            ERROR("Duplicate redistribution cost for " << id);
          }
        },
        make_not_null(&box));
  }
};

struct Collect {
  using simple_tags = tmpl::list<Tags::TargetProc>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Metavariables::volume_dim>& id, ActionList /*meta*/,
      const ParallelComponent* /*meta*/) {
    constexpr size_t dim = Metavariables::volume_dim;
    // Section reduction cookies cannot currently be migrated safely. Reject
    // active sections during collection, before any destinations are applied.
    tmpl::for_each<DbTags>([&box](auto tag_v) {
      using Tag = typename decltype(tag_v)::type;
      if constexpr (tt::is_a_v<Parallel::Tags::Section, Tag>) {
        if (db::get<Tag>(box).has_value()) {
          ERROR("Redistribution does not support active Charm array sections.");
        }
      }
    });
    const auto weight = Parallel::get<domain::Tags::ElementDistribution>(cache);
    const double cost = domain::element_cost(
        db::get<domain::Tags::Mesh<dim>>(box),
        db::get<domain::Tags::Coordinates<dim, Frame::Grid>>(box),
        weight.value_or(domain::ElementWeight::Uniform));
    Parallel::simple_action<RecordCost>(
        Parallel::get_parallel_component<Component<Metavariables>>(cache), id,
        cost, Parallel::my_proc<int>(cache));
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};

struct Plan {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& /*index*/,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    constexpr size_t dim = Metavariables::volume_dim;
    const auto& records = db::get<Tags::Costs<dim>>(box);
    const size_t num_procs = Parallel::number_of_procs<size_t>(cache);
    const auto& ignored = cache.get_resource_info().procs_to_ignore();
    const auto plan = make_plan(
        Parallel::get<domain::Tags::Domain<dim>>(cache).blocks(), records,
        num_procs, ignored,
        Parallel::get<domain::Tags::ElementDistribution>(cache).has_value());
    std::vector<double> before(num_procs, 0.), after(num_procs, 0.);
    size_t moved = 0;
    auto& elements = Parallel::get_parallel_component<
        typename Metavariables::amr::element_array>(cache);
    for (const auto& [id, target] : plan) {
      const auto& [cost, source] = records.at(id);
      before[static_cast<size_t>(source)] += cost;
      after[static_cast<size_t>(target)] += cost;
      moved += source != target;
      Parallel::simple_action<SetTargetProc>(elements[id], target);
    }
    double total = 0.;
    for (const double cost : after) {
      total += cost;
    }
    const double mean = total / (num_procs - ignored.size());
    Parallel::printf(
        "RedistributeElements: %zu of %zu elements moving; predicted "
        "max/mean processor cost %.6f -> %.6f.\n",
        moved, plan.size(),
        *std::max_element(before.begin(), before.end()) / mean,
        *std::max_element(after.begin(), after.end()) / mean);
    db::mutate<Tags::Costs<dim>>([](const auto costs) { costs->clear(); },
                                 make_not_null(&box));
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};

struct Migrate {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& id,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    const auto destination = db::get<Tags::TargetProc>(box);
    if (not destination) {
      ERROR("Missing redistribution destination for " << id);
    }
    db::mutate<Tags::TargetProc>([](const auto target) { target->reset(); },
                                 make_not_null(&box));
    if (*destination != Parallel::my_proc<int>(cache)) {
      // A non-inline entry method runs only after this action has returned.
      Parallel::get_parallel_component<ParallelComponent>(cache)[id].migrate(
          *destination);
    }
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};

}  // namespace Actions

/// Coordinator for deterministic redistribution of an ordinary DG array.
/// Phase boundaries provide quiescence between collection, plan delivery and
/// migration. No restart, topology change, or measured Charm LB is involved.
template <typename Metavariables>
struct Component {
  static std::string name() { return "ElementRedistribution"; }
  using metavariables = Metavariables;
  using chare_type = Parallel::Algorithms::Singleton;
  static constexpr bool checkpoint_data = true;
  using simple_tags_from_options =
      tmpl::list<Tags::Costs<Metavariables::volume_dim>>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Metavariables::volume_dim>,
                 domain::Tags::ElementDistribution>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::PlanElementDistribution, tmpl::list<Actions::Plan>>>;
  static void execute_next_phase(
      const Parallel::Phase phase,
      Parallel::CProxy_GlobalCache<Metavariables>& cache) {
    Parallel::get_parallel_component<Component>(*Parallel::local_branch(cache))
        .start_phase(phase);
  }
};
}  // namespace amr::redistribution
