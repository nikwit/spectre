// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/Actions.hpp"

namespace ScalarWave::Batched::Actions {

template <typename Metavariables>
using InitializeDriver =
    ScalarWave::Actions::InitializeBatchedData<Metavariables>;
using ComputeTimeDerivativeBatched =
    ScalarWave::Actions::ComputeTimeDerivativeBatched;
using ApplyBoundaryCorrectionsToTimeDerivativeBatched =
    ScalarWave::Actions::ApplyBoundaryCorrectionsToTimeDerivativeBatched;
using RecordTimeStepperDataBatched =
    ScalarWave::Actions::RecordTimeStepperDataBatched;
using UpdateUBatched = ScalarWave::Actions::UpdateUBatched;
using CleanHistoryBatched = ScalarWave::Actions::CleanHistoryBatched;

}  // namespace ScalarWave::Batched::Actions
