// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "Time/Tags/FixedLtsRatio.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t VolumeDim>
class Element;
class TimeDelta;
namespace domain::Tags {
template <size_t VolumeDim>
struct Element;
}  // namespace domain::Tags
namespace gsl {
template <class T>
class not_null;
}  // namespace gsl
namespace Tags {
struct TimeStep;
}  // namespace Tags
/// \endcond

namespace evolution::dg::Initialization {
/*!
 * \brief Disable independent local time-stepping on elements with
 * nonconforming boundaries.
 *
 * Nonconforming DG mortars use `TimeSteppingPolicy::EqualRate`, so both sides
 * of the interface must use the same sequence of times. This mutator sets
 * `Tags::FixedLtsRatio` on elements touching a nonconforming face, pinning
 * their LTS step to the initial fraction of a slab.
 */
template <size_t Dim>
struct DisableLtsOnNonconformingBoundaries {
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<::Tags::FixedLtsRatio>;
  using compute_tags = tmpl::list<>;

  using return_tags = tmpl::list<::Tags::FixedLtsRatio>;
  using argument_tags =
      tmpl::list<domain::Tags::Element<Dim>, ::Tags::TimeStep>;

  static void apply(gsl::not_null<std::optional<size_t>*> fixed_lts_ratio,
                    const Element<Dim>& element,
                    const TimeDelta& time_step);
};
}  // namespace evolution::dg::Initialization
