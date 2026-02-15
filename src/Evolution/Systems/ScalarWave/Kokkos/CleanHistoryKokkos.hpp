// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <utility>

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

template <typename System>
struct CleanHistoryKokkos {
 private:
  using device_variables_tag = KokkosTags::DeviceVariables<System>;
  using device_step_start_tag = KokkosTags::DeviceStepStart<System>;

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
           "CleanHistoryKokkos currently supports only Runge-Kutta steppers.");

    if (time_step_id.substep() == runge_kutta->number_of_substeps() - 1) {
      using std::swap;
      swap(*device_step_start, *device_vars);
    }
  }
};

}  // namespace ScalarWave::Actions
