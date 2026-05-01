// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache = true>
struct UpdateFrameFieldsBatchedWithTag {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;
  using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;

 public:
  using return_tags = tmpl::list<packed_geometry_tag>;
  using argument_tags =
      tmpl::list<packed_topology_tag, ::domain::Tags::Domain<volume_dim>,
                 FunctionsOfTimeTag, ::Tags::Time>;
  using const_global_cache_tags = tmpl::conditional_t<
      AddFunctionsOfTimeToConstGlobalCache,
      tmpl::list<::domain::Tags::Domain<volume_dim>, FunctionsOfTimeTag>,
      tmpl::list<::domain::Tags::Domain<volume_dim>>>;

  static void apply(
      gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const typename packed_topology_tag::type& packed_topology,
      const Domain<volume_dim>& domain,
      const domain::FunctionsOfTimeMap& functions_of_time, double time);
};

using UpdateFrameFieldsBatched =
    UpdateFrameFieldsBatchedWithTag<::domain::Tags::FunctionsOfTimeInitialize>;

}  // namespace gh::Actions
