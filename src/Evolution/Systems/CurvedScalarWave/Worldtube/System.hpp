// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
template <size_t Dim>
struct System {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool has_primitive_and_conservative_vars = false;
  // Psi0 and dtPsi0 hold the monopole of the regular field inside the
  // worldtube. At expansion order 2 it is degenerate with the trace of the
  // second-order coefficient and is evolved with an ODE derived from the
  // Klein-Gordon equation. At lower expansion orders the corresponding time
  // derivatives are set to zero and the values are unused.
  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::EvolvedPosition<Dim>, Tags::EvolvedVelocity<Dim>,
                 Tags::Psi0, Tags::dtPsi0>>;
};
}  // namespace CurvedScalarWave::Worldtube
