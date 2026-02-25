// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/InitialExtents.hpp"
#include "Domain/Creators/Tags/InitialRefinementLevels.hpp"
#include "Domain/Domain.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Batched::Initialization {

struct DriverState {
  using const_global_cache_tags = tmpl::list<::domain::Tags::Domain<3>>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options =
      tmpl::list<domain::Tags::InitialRefinementLevels<3>,
                 domain::Tags::InitialExtents<3>,
                 evolution::dg::Tags::Quadrature>;
  using simple_tags = tmpl::list<Tags::DeviceData>;
  using compute_tags = tmpl::list<>;

  using return_tags = simple_tags;
  using argument_tags = tmpl::list<
      ::domain::Tags::Domain<3>, ::domain::Tags::InitialRefinementLevels<3>,
      ::domain::Tags::InitialExtents<3>, evolution::dg::Tags::Quadrature>;

  static void apply(
      const gsl::not_null<typename Tags::DeviceData::type*> device_data,
      const ::Domain<3>& domain,
      const std::vector<std::array<size_t, 3>>& initial_refinement_levels,
      const std::vector<std::array<size_t, 3>>& initial_extents,
      const Spectral::Quadrature& /*quadrature*/) {
    device_data->initialize_from_domain(domain, initial_refinement_levels,
                                        initial_extents);
  }
};

}  // namespace ScalarWave::Batched::Initialization
