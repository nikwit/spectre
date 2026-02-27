// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/Actions.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"

namespace gh::Batched::Actions {

template <typename Metavariables>
using InitializeDriver = gh::Actions::InitializeBatchedData;
using InitializeBoundaryBatchMetadata =
    gh::Actions::InitializeBoundaryBatchMetadata;
using ComputeTimeDerivativeBatched = gh::Actions::ComputeTimeDerivativeBatched;
using PackageLocalFacesBatched = gh::Actions::PackageLocalFacesBatched;
using ComputeInternalBoundaryTermsBatched =
    gh::Actions::ComputeInternalBoundaryTermsBatched;
using LiftInternalBoundaryTermsBatched =
    gh::Actions::LiftInternalBoundaryTermsBatched;
using ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched =
    gh::Actions::ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched;
template <typename FilterType>
using FilterKokkosBatched = gh::Actions::FilterKokkosBatched<FilterType>;

}  // namespace gh::Batched::Actions
