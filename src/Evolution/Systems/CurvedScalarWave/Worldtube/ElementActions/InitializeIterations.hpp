// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"

namespace CurvedScalarWave::Worldtube::Initialization {
struct InitializeIterations {
  using return_tags =
      tmpl::list<CurvedScalarWave::Worldtube::Tags::CurrentIteration>;
  using argument_tags = tmpl::list<>;
  using simple_tags = return_tags;
  using compute_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  static void apply(const gsl::not_null<size_t*> current_iteration) {
    *current_iteration = 0;
  }
};
}  // namespace CurvedScalarWave::Worldtube::Initialization
