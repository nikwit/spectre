// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ObserveTimeStepBatchedEvent.hpp"

namespace gh::Events {
PUP::able::PUP_ID ObserveNormsBatched::my_PUP_ID = 0;  // NOLINT
PUP::able::PUP_ID CompletionBatched::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::Events
