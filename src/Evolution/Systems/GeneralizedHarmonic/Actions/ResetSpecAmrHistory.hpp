// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::Actions {
/// Forget decisions, spatial targets, target shifts and cooldowns computed
/// using old AMR criteria. All evolution and controller state is retained.
struct ResetSpecAmrHistory {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*id*/) {
    db::mutate<spec_amr::HistoryTag>(
        [](const auto history) { *history = spec_amr::History{}; },
        make_not_null(&box));
  }
};
}  // namespace gh::Actions
