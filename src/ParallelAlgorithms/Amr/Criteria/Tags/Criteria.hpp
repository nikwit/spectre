// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "ParallelAlgorithms/Amr/Tags.hpp"
#include "Utilities/Serialization/Serialize.hpp"
#include "Utilities/TMPL.hpp"

namespace amr::Criteria {
/// Option tags for AMR criteria
namespace OptionTags {
/// \ingroup OptionTagsGroup
/// Options for AMR criteria
struct Criteria {
  static constexpr Options::String help = "Options for AMR criteria";
  using type = std::vector<std::unique_ptr<amr::Criterion>>;
  using group = amr::OptionTags::AmrGroup;
};
}  // namespace OptionTags

namespace Tags {
/// The set of adaptive mesh refinement criteria
struct Criteria : db::SimpleTag {
  using type = std::vector<std::unique_ptr<amr::Criterion>>;
  using option_tags = tmpl::list<amr::Criteria::OptionTags::Criteria>;

  static constexpr bool pass_metavariables = false;
  static constexpr bool is_overlayable = true;
  static type create_from_options(const type& value) {
    return {serialize_and_deserialize<type>(value)};
  }

  /// Invalidate criterion-dependent element state when restart options change.
  /// The action must not read the new criteria: other nodes may still be
  /// updating their caches. Restart waits for both operations to quiesce.
  template <typename Metavariables>
  static void on_options_changed(Parallel::GlobalCache<Metavariables>& cache) {
    if constexpr (requires {
                    typename Metavariables::amr::criteria_change_action;
                  }) {
      if (Parallel::my_node<size_t>(cache) == 0) {
        Parallel::printf(
            "AMR criteria changed; resetting per-element criterion history.\n");
        Parallel::simple_action<
            typename Metavariables::amr::criteria_change_action>(
            Parallel::get_parallel_component<
                typename Metavariables::amr::element_array>(cache));
      }
    }
  }
};
}  // namespace Tags
}  // namespace amr::Criteria
