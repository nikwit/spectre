// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/SpectralFilters.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Tag.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"

namespace evolution::dg::Actions {
/// Refresh an element's filter after a restart overlay changes its options.
/// Receive the new filters explicitly: another node's cache may not yet have
/// been updated when this simple action runs.
template <size_t Dim, typename TagList>
struct ResetSpectralFilter {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(
      db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*id*/,
      const std::vector<std::unique_ptr<Filters::Filter<Dim, TagList>>>&
          filters) {
    const auto& element = db::get<domain::Tags::Element<Dim>>(box);
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    db::mutate<Filters::Tags::SpectralFilter<Dim, TagList>>(
        [&filters, &element, &mesh](const auto filter) {
          Initialization::SpectralFilters<Dim, TagList>::apply(filter, filters,
                                                               element, mesh);
        },
        make_not_null(&box));
    if constexpr (db::tag_is_retrievable_v<Tags::MortarDataHistory<Dim>,
                                           db::DataBox<DbTags>>) {
      db::mutate<Tags::MortarDataHistory<Dim>>(
          [](const auto histories) {
            for (auto& [id, history] : *histories) {
              (void)id;
              history.clear_coupling_cache();
            }
          },
          make_not_null(&box));
    }
  }
};
}  // namespace evolution::dg::Actions
