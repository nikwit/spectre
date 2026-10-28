// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <utility>

#include "DataStructures/TaggedTuple.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/ContributeToPhaseChangeReduction.hpp"
#include "Parallel/PhaseControl/PhaseChange.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::bbh::phase_control {
namespace redistribution {
struct State {
  std::optional<Parallel::Phase> return_phase{};
  void pup(PUP::er& p) { p | return_phase; }
};
struct StateTag {
  using type = State;
  struct combine_method {
    State operator()(const State&, const State&) const {
      ERROR("Redistribution state is owned by Main.");
    }
  };
  using main_combine_method = combine_method;
};
struct Requested {
  using type = bool;
  using combine_method = funcl::Or<>;
  using main_combine_method = combine_method;
};
}  // namespace redistribution

/// Recompute placement from current mesh costs at synchronized slab boundaries.
/// Configure once in PhaseChangeAndTriggers, after AMR if both can trigger
/// together. Redistribution has no junk-traversal timing restriction.
class RedistributeElements : public PhaseChange {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help =
      "Redistribute the fixed h mesh using current p costs and the initial "
      "partitioner. No Charm load balancer or junk-traversal delay.";
  RedistributeElements() = default;
  explicit RedistributeElements(CkMigrateMessage* msg) : PhaseChange(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(RedistributeElements);  // NOLINT
  void pup(PUP::er& p) override;
  using argument_tags = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using phase_change_tags_and_combines =
      tmpl::list<redistribution::StateTag, redistribution::Requested>;
  template <typename Metavariables>
  using participating_components =
      tmpl::list<typename Metavariables::amr::element_array>;

  template <typename ParallelComponent, typename ArrayIndex,
            typename Metavariables>
  void contribute_phase_data_impl(Parallel::GlobalCache<Metavariables>& cache,
                                  const ArrayIndex& id) const {
    Parallel::contribute_to_phase_change_reduction<ParallelComponent>(
        tuples::TaggedTuple<redistribution::Requested>{true}, cache, id);
  }

  template <typename... Tags>
  void initialize_phase_data_impl(
      const gsl::not_null<tuples::TaggedTuple<Tags...>*> data) const {
    tuples::get<redistribution::StateTag>(*data) = redistribution::State{};
    tuples::get<redistribution::Requested>(*data) = false;
  }

  template <typename... Tags, typename Metavariables>
  std::optional<std::pair<Parallel::Phase, PhaseControl::ArbitrationStrategy>>
  arbitrate_phase_change_impl(
      const gsl::not_null<tuples::TaggedTuple<Tags...>*> data,
      const Parallel::Phase phase,
      const Parallel::GlobalCache<Metavariables>& /*cache*/) const {
    return next_phase(
        make_not_null(&tuples::get<redistribution::StateTag>(*data)),
        make_not_null(&tuples::get<redistribution::Requested>(*data)), phase);
  }

 private:
  std::optional<std::pair<Parallel::Phase, PhaseControl::ArbitrationStrategy>>
  next_phase(gsl::not_null<redistribution::State*> state,
             gsl::not_null<bool*> requested, Parallel::Phase phase) const;
};
}  // namespace gh::bbh::phase_control
