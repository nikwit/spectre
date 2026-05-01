// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/Actions.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Batched::Actions {

template <typename Metavariables>
using functions_of_time_compute_tag =
    Parallel::Tags::FromGlobalCache<domain::Tags::FunctionsOfTimeInitialize,
                                    Metavariables>;

template <typename Metavariables>
struct InitializeDriver : gh::Actions::InitializeBatchedDataWithTag<
                              domain::Tags::FunctionsOfTimeInitialize, false> {
  using compute_tags = tmpl::list<functions_of_time_compute_tag<Metavariables>>;
};

template <typename Metavariables>
struct UpdateFrameFieldsBatched
    : gh::Actions::UpdateFrameFieldsBatchedWithTag<
          domain::Tags::FunctionsOfTimeInitialize, false> {
  using compute_tags = tmpl::list<functions_of_time_compute_tag<Metavariables>>;
};

template <typename Metavariables>
struct UpdateConstraintDampingBatched
    : gh::Actions::UpdateConstraintDampingBatchedWithTag<
          domain::Tags::FunctionsOfTimeInitialize, false> {
  using compute_tags = tmpl::list<functions_of_time_compute_tag<Metavariables>>;
};

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
