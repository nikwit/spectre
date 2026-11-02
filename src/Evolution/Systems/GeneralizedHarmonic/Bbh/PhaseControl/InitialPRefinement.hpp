// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <utility>

#include "DataStructures/TaggedTuple.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/PhaseChange.hpp"
#include "Parallel/PhaseControl/PhaseControlTags.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::bbh::phase_control {
namespace initial_p_refinement {
struct State {
  size_t iteration{0};
  bool finished{false};
  bool changed{false};
  void pup(PUP::er& p) {
    p | iteration;
    p | finished;
    p | changed;
  }
};
struct StateTag {
  using type = State;
  struct combine_method {
    State operator()(const State&, const State&) const {
      ERROR("Initial p-refinement state is owned by Main.");
    }
  };
  using main_combine_method = combine_method;
};
struct ChangedElements {
  using type = size_t;
  using combine_method = funcl::Plus<>;
  using main_combine_method = combine_method;
};
struct LimitedDirections : ChangedElements {};
struct CyclingElements : ChangedElements {};
}  // namespace initial_p_refinement

/// Iterate p adjustment before time-stepper self-start, resampling initial data
/// after every change. Configure once in PhaseChangeAndTriggers with
/// `Trigger: {Not: Always}`: startup is automatic, with no evolution trigger.
/// Requires the BBH startup actions and SpecSpectral as the sole criterion.
/// The initial h mesh and element ownership are retained. Reaching the
/// iteration limit or a nonfatal p cap is reported, not claimed to satisfy the
/// targets.
class InitialPRefinement : public PhaseChange {
 public:
  struct MaxIterations {
    using type = size_t;
    static size_t lower_bound() { return 1; }
    static constexpr Options::String help =
        "Maximum initial p-adjustment iterations (SpEC BBH uses 40).";
  };
  using options = tmpl::list<MaxIterations>;
  static constexpr Options::String help =
      "Adjust the initial p mesh before evolution and reload the initial data. "
      "Use Trigger: {Not: Always}; this phase change runs automatically at "
      "startup.";
  InitialPRefinement() = default;
  explicit InitialPRefinement(size_t max_iterations);
  explicit InitialPRefinement(CkMigrateMessage* msg) : PhaseChange(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(InitialPRefinement);  // NOLINT
  void pup(PUP::er& p) override;
  using argument_tags = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using phase_change_tags_and_combines =
      tmpl::list<initial_p_refinement::StateTag,
                 initial_p_refinement::ChangedElements,
                 initial_p_refinement::LimitedDirections,
                 initial_p_refinement::CyclingElements>;
  template <typename Metavariables>
  using participating_components = tmpl::list<>;

  template <typename ParallelComponent, typename ArrayIndex,
            typename Metavariables>
  void contribute_phase_data_impl(Parallel::GlobalCache<Metavariables>&,
                                  const ArrayIndex&) const {}

  template <typename... Tags>
  void initialize_phase_data_impl(
      const gsl::not_null<tuples::TaggedTuple<Tags...>*> data) const {
    using namespace initial_p_refinement;
    tuples::get<StateTag>(*data) = State{};
    tuples::get<ChangedElements>(*data) = 0;
    tuples::get<LimitedDirections>(*data) = 0;
    tuples::get<CyclingElements>(*data) = 0;
  }

  template <typename... Tags, typename Metavariables>
  std::optional<std::pair<Parallel::Phase, PhaseControl::ArbitrationStrategy>>
  arbitrate_phase_change_impl(
      const gsl::not_null<tuples::TaggedTuple<Tags...>*> data,
      const Parallel::Phase current_phase,
      const Parallel::GlobalCache<Metavariables>& /*cache*/) const {
    using namespace initial_p_refinement;
    using Parallel::Phase;
    auto& state = tuples::get<StateTag>(*data);
    if (state.finished) {
      return std::nullopt;
    }
    const auto next = next_phase(make_not_null(&state), current_phase,
                                 tuples::get<ChangedElements>(*data),
                                 tuples::get<LimitedDirections>(*data),
                                 tuples::get<CyclingElements>(*data));
    if (current_phase == Phase::CheckInitialPRefinement) {
      tuples::get<ChangedElements>(*data) = 0;
      tuples::get<LimitedDirections>(*data) = 0;
      tuples::get<CyclingElements>(*data) = 0;
    }
    if (next.has_value()) {
      return std::pair{*next, PhaseControl::RunPhaseImmediately};
    }
    return std::nullopt;
  }

 private:
  std::optional<Parallel::Phase> next_phase(
      gsl::not_null<initial_p_refinement::State*> state,
      Parallel::Phase current_phase, size_t changed, size_t limited,
      size_t cycling) const;
  size_t max_iterations_{40};
};

/// Check configuration, rejecting duplicate controllers that would share state.
template <typename Metavariables>
bool initial_p_refinement_enabled(
    const Parallel::GlobalCache<Metavariables>& cache) {
  size_t count = 0;
  for (const auto& entry :
       Parallel::get<PhaseControl::Tags::PhaseChangeAndTriggers>(cache)) {
    for (const auto& change : entry.phase_changes) {
      count += dynamic_cast<const InitialPRefinement*>(change.get()) != nullptr;
    }
  }
  if (count > 1) {
    ERROR("Configure InitialPRefinement only once.");
  }
  return count == 1;
}
}  // namespace gh::bbh::phase_control
