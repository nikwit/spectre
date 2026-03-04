// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Tags/Filter.hpp"
#include "NumericalAlgorithms/LinearOperators/ExponentialFilter.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename FilterType>
struct FilterKokkosBatched {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<system>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;

 public:
  using return_tags = tmpl::list<packed_evolution_state_tag>;
  using argument_tags =
      tmpl::list<packed_topology_tag, ::Filters::Tags::Filter<FilterType>>;
  using const_global_cache_tags =
      tmpl::list<::Filters::Tags::Filter<FilterType>>;

  static void apply(
      const gsl::not_null<typename packed_evolution_state_tag::type*>
          packed_evolution_state,
      const typename packed_topology_tag::type& packed_topology,
      const FilterType& filter_helper);
};

extern template struct FilterKokkosBatched<Filters::Exponential<0>>;

}  // namespace gh::Actions
