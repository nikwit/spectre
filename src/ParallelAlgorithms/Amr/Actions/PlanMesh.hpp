// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "Domain/Amr/Tags/Flags.hpp"
#include "Domain/Amr/Tags/NeighborFlags.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Tag.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Amr/Actions/AdjustDomain.hpp"
#include "ParallelAlgorithms/Amr/MeshPlan.hpp"
#include "ParallelAlgorithms/Amr/Policies/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Literals.hpp"

namespace amr {
template <typename Metavariables>
struct Component;
}

namespace amr::planning {
/// Coordinator state only. AMR phase quiescence separates planning rounds.
/// No element is created/destroyed until both reductions have completed.
struct State {
  size_t round = 0;
  enum class Stage { Idle, Collecting, Installing } stage = Stage::Idle;
  size_t expected_installations = 0;
  void pup(PUP::er& p) {
    p | round;
    p | stage;
    p | expected_installations;
  }
};
namespace Tags {
struct State : db::SimpleTag {
  using type = planning::State;
  using option_tags = tmpl::list<>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options() { return {}; }
};
}  // namespace Tags

/// Tree reduction of compact metadata, with duplicate-element detection.
struct CombineRequests {
  Requests<3> operator()(Requests<3> a, const Requests<3>& b) const {
    for (const auto& [id, request] : b) {
      if (not a.emplace(id, request).second) {
        ERROR("Duplicate element in central AMR collection: " << id);
      }
    }
    return a;
  }
};
using RequestReduction = Parallel::ReductionData<
    Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
    Parallel::ReductionDatum<Requests<3>, CombineRequests>>;
using InstalledReduction = Parallel::ReductionData<
    Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
    Parallel::ReductionDatum<size_t, funcl::Plus<>>>;

/// Determine the hard angular cutoff from the actual element volume filter.
/// Smooth attenuation below the cutoff is deliberately not counted as killing
/// modes. Unknown filters remain unknown instead of assuming a top-four filter.
template <typename Metavariables, typename DbTags>
std::optional<size_t> filtered_modes(const db::DataBox<DbTags>& box,
                                     const Mesh<3>& mesh) {
  if (mesh.basis(1) != Spectral::Basis::SphericalHarmonic) {
    return std::nullopt;
  }
  if constexpr (requires { typename Metavariables::system::variables_tag; }) {
    using VariablesTag = typename Metavariables::system::variables_tag;
    using FilterTag =
        Filters::Tags::SpectralFilter<3, typename VariablesTag::tags_list>;
    if constexpr (db::tag_is_retrievable_v<FilterTag, db::DataBox<DbTags>>) {
      const auto& filter = db::get<FilterTag>(box);
      const auto response = filter.volume_filter_response(mesh);
      if (not response or (*response)[1].size() != mesh.extents(1) or
          (*response)[2] != (*response)[1]) {
        return std::nullopt;
      }
      const auto& angular = (*response)[1];
      size_t retained = 0;
      bool zero_seen = false;
      for (const double value : angular) {
        if (not std::isfinite(value) or value < 0.0 or value > 1.0 or
            (zero_seen and value != 0.0)) {
          return std::nullopt;
        }
        zero_seen |= value == 0.0;
        retained += value != 0.0 ? 1 : 0;
      }
      if (retained == 0) {
        return std::nullopt;
      }
      return angular.size() - retained;
    }
  }
  return std::nullopt;
}

namespace Actions {
/// All recipients have installed their own and their neighbors' final meshes.
/// Only now may the existing distributed projection actions mutate fields.
template <typename ApplyAction = amr::Actions::AdjustDomain>
struct Apply {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const int /*index*/, const size_t round,
                    const size_t installed) {
    auto& state = db::get_mutable_reference<Tags::State>(make_not_null(&box));
    if (state.stage != State::Stage::Installing or round != state.round or
        installed != state.expected_installations) {
      ERROR("Incomplete or stale central AMR plan acknowledgement");
    }
    state.stage = State::Stage::Idle;
    state.expected_installations = 0;
    Parallel::simple_action<ApplyAction>(
        Parallel::get_parallel_component<
            typename Metavariables::amr::element_array>(cache));
  }
};

template <typename ApplyAction = amr::Actions::AdjustDomain>
struct Install {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ElementId<3>& id, const size_t round,
                    const Mesh<3>& expected_current_mesh,
                    const Decision<3>& decision,
                    const std::unordered_map<ElementId<3>, Info<3>>& neighbors,
                    const std::optional<size_t> expected_filtered_modes) {
    if (db::get<domain::Tags::Mesh<3>>(box) != expected_current_mesh or
        filtered_modes<Metavariables>(box, decision.info.new_mesh) !=
            expected_filtered_modes) {
      ERROR("Mesh or filter changed while installing central AMR plan for "
            << id);
    }
    db::mutate<amr::Tags::Info<3>, amr::Tags::NeighborInfo<3>>(
        [&decision, &neighbors](const auto info, const auto neighbor_info) {
          *info = decision.info;
          *neighbor_info = neighbors;
        },
        make_not_null(&box));
    Parallel::contribute_to_reduction<Apply<ApplyAction>>(
        InstalledReduction{round, 1_st},
        Parallel::get_parallel_component<ParallelComponent>(cache)[id],
        Parallel::get_parallel_component<amr::Component<Metavariables>>(cache));
  }
};

template <typename ApplyAction = amr::Actions::AdjustDomain>
struct Plan {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const int /*index*/, const size_t round,
                    const Requests<3>& requests) {
    auto& state = db::get_mutable_reference<Tags::State>(make_not_null(&box));
    if (state.stage != State::Stage::Collecting or round != state.round) {
      ERROR("Unexpected central AMR proposal round");
    }
    const auto plan =
        make_plan(requests, Parallel::get<domain::Tags::Domain<3>>(cache),
                  Parallel::get<amr::Tags::Policies>(cache));
    state.stage = State::Stage::Installing;
    state.expected_installations = plan.size();
    auto& array = Parallel::get_parallel_component<
        typename Metavariables::amr::element_array>(cache);
    for (const auto& [id, decision] : plan) {
      std::unordered_map<ElementId<3>, Info<3>> neighbors;
      for (const auto& [direction, ns] : requests.at(id).element.neighbors()) {
        (void)direction;
        for (const auto& neighbor : ns.ids()) {
          neighbors.emplace(neighbor, plan.at(neighbor).info);
        }
      }
      Parallel::simple_action<Install<ApplyAction>>(
          array[id], round, requests.at(id).current_mesh, decision,
          std::move(neighbors), requests.at(id).filtered_angular_modes);
    }
  }
};

template <typename ApplyAction = amr::Actions::AdjustDomain>
struct Collect {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ElementId<3>& id, const size_t round) {
    const auto& proposal = db::get<amr::Tags::Info<3>>(box);
    Requests<3> requests;
    requests.emplace(
        id, Request<3>{db::get<domain::Tags::Element<3>>(box),
                       db::get<domain::Tags::Mesh<3>>(box), proposal,
                       filtered_modes<Metavariables>(box, proposal.new_mesh)});
    Parallel::contribute_to_reduction<Plan<ApplyAction>>(
        RequestReduction{round, std::move(requests)},
        Parallel::get_parallel_component<ParallelComponent>(cache)[id],
        Parallel::get_parallel_component<amr::Component<Metavariables>>(cache));
  }
};

template <typename ApplyAction = amr::Actions::AdjustDomain>
struct Begin {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const int /*index*/) {
    auto& state = db::get_mutable_reference<Tags::State>(make_not_null(&box));
    if (state.stage != State::Stage::Idle) {
      ERROR("Overlapping central AMR planning rounds");
    }
    ++state.round;
    state.stage = State::Stage::Collecting;
    Parallel::simple_action<Collect<ApplyAction>>(
        Parallel::get_parallel_component<
            typename Metavariables::amr::element_array>(cache),
        state.round);
  }
};
}  // namespace Actions
}  // namespace amr::planning
