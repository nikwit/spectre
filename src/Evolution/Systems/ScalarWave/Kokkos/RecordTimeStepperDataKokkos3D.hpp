// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct RecordTimeStepperDataKokkos3D {
 private:
  using system = ScalarWave::System<3>;
  using device_dt_variables_tag = KokkosTags::DeviceDtVariables<system>;
  using device_derivative_history_tag =
      KokkosTags::DeviceDerivativeHistory<system>;

 public:
  using return_tags = tmpl::list<device_derivative_history_tag>;
  using argument_tags = tmpl::list<::Tags::TimeStepId, device_dt_variables_tag>;

  static void apply(
      gsl::not_null<device_derivative_history_tag::type*>
          device_derivative_history,
      const TimeStepId& time_step_id,
      const device_dt_variables_tag::type& device_dt);
};

}  // namespace ScalarWave::Actions
