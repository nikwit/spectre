// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Bbh/PhaseControl/RedistributeElements.hpp"

namespace gh::bbh::phase_control {
void RedistributeElements::pup(PUP::er& p) { PhaseChange::pup(p); }

std::optional<std::pair<Parallel::Phase, PhaseControl::ArbitrationStrategy>>
RedistributeElements::next_phase(
    const gsl::not_null<redistribution::State*> state,
    const gsl::not_null<bool*> requested, const Parallel::Phase phase) const {
  using Parallel::Phase;
  using PhaseControl::RunPhaseImmediately;
  if (state->return_phase) {
    if (phase == Phase::CollectElementDistribution) {
      return std::pair{Phase::PlanElementDistribution, RunPhaseImmediately};
    }
    if (phase == Phase::PlanElementDistribution) {
      return std::pair{Phase::RedistributeElements, RunPhaseImmediately};
    }
    if (phase == Phase::RedistributeElements) {
      return std::pair{Phase::UpdateSections, RunPhaseImmediately};
    }
    if (phase == Phase::UpdateSections) {
      const auto result = *state->return_phase;
      state->return_phase.reset();
      return std::pair{result, PhaseControl::PermitAdditionalJumps};
    }
    ERROR("Unexpected phase during element redistribution: " << phase);
  }
  if (*requested) {
    *requested = false;
    state->return_phase = phase;
    return std::pair{Phase::CollectElementDistribution, RunPhaseImmediately};
  }
  return std::nullopt;
}

PUP::able::PUP_ID RedistributeElements::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::bbh::phase_control
