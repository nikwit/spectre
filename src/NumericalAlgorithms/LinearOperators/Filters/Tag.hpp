// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/OptionTags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags/OptionsGroup.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Filter.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Utilities/CloneUniquePtrs.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/TMPL.hpp"

namespace Filters::OptionTags {
/// \brief Option tag for the list of spectral filters read from the input file.
///
/// Parsed under the `Filtering` key in the input file, this tag holds a
/// `std::vector` of heap-allocated `Filters::Filter<Dim, TagList>` objects.
///
/// \tparam Dim Spatial dimension of the element mesh.
/// \tparam TagList `tmpl::list` of tensor tags in the `Variables` to filter.
template <size_t Dim, typename TagList>
struct SpectralFilters {
  using group = dg::OptionTags::DiscontinuousGalerkinGroup;
  static std::string name() { return "Filtering"; }
  static constexpr Options::String help =
      "A vector/list of the different filters that are applied in different "
      "elements.";
  using type = std::vector<std::unique_ptr<Filters::Filter<Dim, TagList>>>;
};
}  // namespace Filters::OptionTags

namespace Filters::Tags {
/// \brief DataBox tag for the spectral filters applied during DG time
/// integration in a single element.
///
/// Holds `Filters::Filter<Dim, TagList>`s that are applied in the element.
/// `TagList` is the `tmpl::list` of tensor tags in the `Variables` to be
/// filtered.
///
/// \tparam Dim Spatial dimension of the element mesh.
/// \tparam TagList `tmpl::list` of tensor tags in the `Variables` to filter.
template <size_t Dim, typename TagList>
struct SpectralFilter : db::SimpleTag {
  using type = std::unique_ptr<Filters::Filter<Dim, TagList>>;
};

/// \brief DataBox tag for the ordered list of all spectral filters for
/// different topologies applied during DG time integration.
///
/// Each element holds a local copy of the filters that should be applied in it.
///
/// \tparam Dim Spatial dimension of the element mesh.
/// \tparam TagList `tmpl::list` of tensor tags in the `Variables` to filter.
template <size_t Dim, typename TagList>
struct SpectralFilters : db::SimpleTag {
  using type = std::vector<std::unique_ptr<Filters::Filter<Dim, TagList>>>;
  static constexpr bool pass_metavariables = false;
  static constexpr bool is_overlayable = true;

  using option_tags =
      tmpl::list<Filters::OptionTags::SpectralFilters<Dim, TagList>,
                 domain::OptionTags::DomainCreator<Dim>>;

  static type create_from_options(
      const type& spectral_filters,
      const std::unique_ptr<DomainCreator<Dim>>& domain_creator) {
    auto filters = clone_unique_ptrs(spectral_filters);
    const auto& block_names = domain_creator->block_names();
    const auto& block_groups = domain_creator->block_groups();
    for (auto& filter : filters) {
      filter->set_blocks_to_filter(block_names, block_groups);
    }
    return filters;
  }

  /// A checkpoint also contains a local clone on each element. Updating only
  /// the global list would silently leave those clones using the old cadence.
  template <typename Metavariables>
  static void on_options_changed(Parallel::GlobalCache<Metavariables>& cache) {
    if constexpr (requires {
                    typename Metavariables::
                        spectral_filter_options_changed_action;
                    typename Metavariables::amr::element_array;
                  }) {
      if (Parallel::my_node<size_t>(cache) == 0) {
        Parallel::printf(
            "Spectral filter options changed; refreshing per-element "
            "filters.\n");
        Parallel::simple_action<
            typename Metavariables::spectral_filter_options_changed_action>(
            Parallel::get_parallel_component<
                typename Metavariables::amr::element_array>(cache),
            clone_unique_ptrs(Parallel::get<SpectralFilters>(cache)));
      }
    } else {
      ERROR(
          "This executable does not support changing spectral filters on "
          "a checkpoint restart: it must provide an action to refresh "
          "the per-element filters.");
    }
  }
};
}  // namespace Filters::Tags
