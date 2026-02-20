// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace Actions::Kokkos {

template <typename System, template <typename> class DeviceVariablesTag,
          template <typename> class DeviceStepStartTag>
struct CleanHistory {
 private:
  using device_variables_tag = DeviceVariablesTag<System>;
  using device_step_start_tag = DeviceStepStartTag<System>;

 public:
  using return_tags = tmpl::list<device_step_start_tag, device_variables_tag>;
  using argument_tags =
      tmpl::list<::Tags::TimeStepper<TimeStepper>, ::Tags::TimeStepId>;

  static void apply(
      const gsl::not_null<typename device_step_start_tag::type*>
          device_step_start,
      const gsl::not_null<typename device_variables_tag::type*> device_vars,
      const TimeStepper& time_stepper, const TimeStepId& time_step_id) {
    const auto* runge_kutta =
        dynamic_cast<const TimeSteppers::RungeKutta*>(&time_stepper);
    ASSERT(runge_kutta != nullptr,
           "Kokkos CleanHistory currently supports only Runge-Kutta "
           "steppers.");

    if (time_step_id.substep() == runge_kutta->number_of_substeps() - 1) {
      // Preserve the updated solution in device_vars and advance the
      // step-start state for the next full step.
      ::Kokkos::deep_copy(device_step_start->view(), device_vars->view());
    }
  }
};

}  // namespace Actions::Kokkos
