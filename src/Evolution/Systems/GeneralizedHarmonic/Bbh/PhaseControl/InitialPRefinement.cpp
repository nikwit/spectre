// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Bbh/PhaseControl/InitialPRefinement.hpp"

namespace gh::bbh::phase_control {
InitialPRefinement::InitialPRefinement(const size_t max_iterations)
    : max_iterations_(max_iterations) {
  if (max_iterations == 0) {
    ERROR("InitialPRefinement requires MaxIterations > 0.");
  }
}

void InitialPRefinement::pup(PUP::er& p) {
  PhaseChange::pup(p);
  p | max_iterations_;
}

std::optional<Parallel::Phase> InitialPRefinement::next_phase(
    const gsl::not_null<initial_p_refinement::State*> state,
    const Parallel::Phase phase, const size_t changed, const size_t limited,
    const size_t cycling) const {
  using Parallel::Phase;
  if (phase == Phase::InitializeInitialDataDependentQuantities) {
    if (state->iteration == max_iterations_) {
      Parallel::printf(
          "InitialPRefinement: reached %zu iterations; final initial "
          "data have been resampled. Targets may remain unmet.\n",
          max_iterations_);
      state->finished = true;
      return Phase::Register;
    }
    ++state->iteration;
    return Phase::EvaluateAmrCriteria;
  }
  if (state->iteration == 0) {
    return std::nullopt;
  }
  if (phase == Phase::EvaluateAmrCriteria) {
    return Phase::CheckInitialPRefinement;
  }
  if (phase == Phase::CheckInitialPRefinement) {
    Parallel::printf(
        "InitialPRefinement iteration=%zu changed_elements=%zu "
        "policy_limited_directions=%zu cycling_elements=%zu\n",
        state->iteration, changed, limited, cycling);
    // Even when no element changes, finish neighbor bookkeeping and clear
    // flags.
    return Phase::AdjustDomain;
  }
  if (phase == Phase::AdjustDomain) {
    return Phase::CheckDomain;
  }
  if (phase == Phase::CheckDomain) {
    state->changed = changed > 0;
    Parallel::printf("InitialPRefinement realized_changed_elements=%zu\n",
                     changed);
    if (state->changed) {
      return Phase::RegisterWithElementDataReader;
    }
    Parallel::printf("InitialPRefinement: no further permitted p changes.\n");
    state->finished = true;
    return Phase::Register;
  }
  return std::nullopt;
}

PUP::able::PUP_ID InitialPRefinement::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::bbh::phase_control
