// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct UpdateUKokkos3D {
 private:
  using system = ScalarWave::System<3>;
  using device_variables_tag = KokkosTags::DeviceVariables<system>;
  using device_step_start_tag = KokkosTags::DeviceStepStart<system>;
  using device_derivative_history_tag =
      KokkosTags::DeviceDerivativeHistory<system>;

 public:
  using return_tags = tmpl::list<device_variables_tag>;
  using argument_tags =
      tmpl::list<::Tags::TimeStepper<TimeStepper>, ::Tags::TimeStepId,
                 ::Tags::TimeStep, device_step_start_tag,
                 device_derivative_history_tag>;

  static void apply(
      gsl::not_null<device_variables_tag::type*> device_vars,
      const TimeStepper& time_stepper, const TimeStepId& time_step_id,
      const TimeDelta& time_step,
      const device_step_start_tag::type& device_step_start,
      const device_derivative_history_tag::type&
          device_derivative_history);
};

}  // namespace ScalarWave::Actions
