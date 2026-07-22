// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/OptionTags.hpp"
#include "Options/String.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/StdHelpers.hpp"

namespace OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * \brief Groups the filtering configurations in the input file.
 */
struct FilteringGroup {
  static std::string name() { return "Filtering"; }
  static constexpr Options::String help = "Options for filtering";
};

/*!
 * \ingroup OptionTagsGroup
 * \brief The option tag that retrieves the parameters for the filter
 * from the input file
 */
template <typename FilterType>
struct Filter {
  static std::string name() { return pretty_type::name<FilterType>(); }
  static constexpr Options::String help = "Options for the filter";
  using type = FilterType;
  using group = FilteringGroup;
};
}  // namespace OptionTags

namespace Filters::Tags {
/// Compatibility global-cache tag for the concrete filters used by the
/// Kokkos drivers. The production DG path uses `SpectralFilters` instead.
template <typename FilterType>
struct Filter : db::SimpleTag {
  using type = FilterType;
  template <typename Metavariables>
  using option_tags =
      tmpl::list<::OptionTags::Filter<FilterType>,
                 domain::OptionTags::DomainCreator<Metavariables::volume_dim>>;

  static constexpr bool pass_metavariables = true;
  template <typename Metavariables>
  static FilterType create_from_options(
      const FilterType& filter,
      const std::unique_ptr<DomainCreator<Metavariables::volume_dim>>&
          domain_creator) {
    const auto& blocks_to_filter = filter.blocks_to_filter();
    if (blocks_to_filter.has_value()) {
      const auto& block_names = domain_creator->block_names();
      const auto& block_groups = domain_creator->block_groups();
      if (block_names.empty()) {
        ERROR(
            "The domain chosen doesn't use block names, but the Filter tag "
            "has specified block names to use.");
      }
      for (const std::string& block_to_filter : blocks_to_filter.value()) {
        if (not alg::found(block_names, block_to_filter) and
            block_groups.count(block_to_filter) == 0) {
          ERROR("Specified block (group) name '"
                << block_to_filter
                << "' is not a block name or block group. Existing blocks "
                   "are:\n"
                << block_names << "\nExisting block groups are:\n"
                << keys_of(block_groups));
        }
      }
    }
    return filter;
  }
};
}  // namespace Filters::Tags
