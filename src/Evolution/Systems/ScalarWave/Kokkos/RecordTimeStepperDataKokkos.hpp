// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

template <typename System>
struct RecordTimeStepperDataKokkos {
 private:
  using device_dt_variables_tag = KokkosTags::DeviceDtVariables<System>;
  using device_derivative_history_tag =
      KokkosTags::DeviceDerivativeHistory<System>;

 public:
  using return_tags = tmpl::list<device_derivative_history_tag>;
  using argument_tags = tmpl::list<::Tags::TimeStepId, device_dt_variables_tag>;

  static void apply(
      const gsl::not_null<typename device_derivative_history_tag::type*>
          device_derivative_history,
      const TimeStepId& time_step_id,
      const typename device_dt_variables_tag::type& device_dt) {
    const size_t substep = time_step_id.substep();
    ASSERT(substep < static_cast<size_t>(device_derivative_history->extent(0)),
           "Substep " << substep << " exceeds derivative history size "
                      << device_derivative_history->extent(0));

    constexpr size_t number_of_components =
        device_dt_variables_tag::type::number_of_independent_components;
    const size_t num_points = device_dt.number_of_grid_points();
    const auto dt_view = device_dt.view();
    const auto deriv_history = *device_derivative_history;
    Kokkos::parallel_for(
        "RecordTimeStepperDataKokkos", num_points, KOKKOS_LAMBDA(const int i) {
          for (size_t c = 0; c < number_of_components; ++c) {
            deriv_history(substep, i, c) = dt_view(i, c);
          }
        });
  }
};

}  // namespace ScalarWave::Actions
