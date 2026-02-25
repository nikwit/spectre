// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/Actions.hpp"

namespace ScalarWave::Batched::Actions {

template <typename Metavariables>
using InitializeDriver =
    ScalarWave::Actions::InitializeBatchedData<Metavariables>;
using InitializeBoundaryBatchMetadata =
    ScalarWave::Actions::InitializeBoundaryBatchMetadata;
using ComputeTimeDerivativeBatched =
    ScalarWave::Actions::ComputeTimeDerivativeBatched;
using PackageLocalFacesBatched = ScalarWave::Actions::PackageLocalFacesBatched;
using ComputeInternalBoundaryTermsBatched =
    ScalarWave::Actions::ComputeInternalBoundaryTermsBatched;
using LiftInternalBoundaryTermsBatched =
    ScalarWave::Actions::LiftInternalBoundaryTermsBatched;
using ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched = ScalarWave::
    Actions::ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched;
using RecordTimeStepperDataBatched =
    ScalarWave::Actions::RecordTimeStepperDataBatched;
using UpdateUBatched = ScalarWave::Actions::UpdateUBatched;
using CleanHistoryBatched = ScalarWave::Actions::CleanHistoryBatched;

}  // namespace ScalarWave::Batched::Actions
