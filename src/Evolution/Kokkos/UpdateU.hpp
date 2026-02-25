// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "Evolution/Kokkos/PackedTags.hpp"
#include "Time/Kokkos/UpdateU.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution::Actions::Kokkos {

template <typename System>
struct UpdateU {
 private:
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<System>;

 public:
  using return_tags = tmpl::list<packed_evolution_state_tag>;
  using argument_tags = tmpl::list<::Tags::TimeStepper<TimeStepper>,
                                   ::Tags::TimeStepId, ::Tags::TimeStep>;

  static void apply(
      const gsl::not_null<typename packed_evolution_state_tag::type*>
          packed_evolution_state,
      const TimeStepper& time_stepper, const TimeStepId& time_step_id,
      const TimeDelta& time_step);
};

}  // namespace evolution::Actions::Kokkos
