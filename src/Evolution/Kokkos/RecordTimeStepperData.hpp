// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Time/Kokkos/RecordTimeStepperData.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution::Actions::Kokkos {

template <typename System>
struct RecordTimeStepperData {
 private:
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<System>;

 public:
  using return_tags = tmpl::list<packed_evolution_state_tag>;
  using argument_tags = tmpl::list<::Tags::TimeStepId>;

  static void apply(
      const gsl::not_null<typename packed_evolution_state_tag::type*>
          packed_evolution_state,
      const TimeStepId& time_step_id);
};

}  // namespace evolution::Actions::Kokkos
