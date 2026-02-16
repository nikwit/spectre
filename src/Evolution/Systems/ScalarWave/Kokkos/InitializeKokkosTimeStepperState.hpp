// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/History.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

template <typename System>
struct InitializeKokkosTimeStepperState {
 private:
  using host_dt_variables_tag =
      db::add_tag_prefix<::Tags::dt, typename System::variables_tag>;
  using device_variables_tag = KokkosTags::DeviceVariables<System>;
  using device_dt_variables_tag = KokkosTags::DeviceDtVariables<System>;
  using device_step_start_tag = KokkosTags::DeviceStepStart<System>;
  using device_derivative_history_tag =
      KokkosTags::DeviceDerivativeHistory<System>;

 public:
  using simple_tags =
      tmpl::list<device_variables_tag, device_dt_variables_tag,
                 device_step_start_tag, device_derivative_history_tag>;
  using return_tags =
      tmpl::list<device_variables_tag, device_dt_variables_tag,
                 device_step_start_tag, device_derivative_history_tag>;
  using argument_tags =
      tmpl::list<typename System::variables_tag, host_dt_variables_tag>;

  static void apply(
      const gsl::not_null<typename device_variables_tag::type*> device_vars,
      const gsl::not_null<typename device_dt_variables_tag::type*> device_dt,
      const gsl::not_null<typename device_step_start_tag::type*>
          device_step_start,
      const gsl::not_null<typename device_derivative_history_tag::type*>
          device_derivative_history,
      const typename System::variables_tag::type& host_vars,
      const typename host_dt_variables_tag::type& host_dt_vars) {
    *device_vars = copy_to_device(host_vars);
    *device_dt = copy_to_device(host_dt_vars);
    *device_step_start = copy_to_device(host_vars);

    *device_derivative_history = typename device_derivative_history_tag::type(
        "KokkosDerivativeHistory", TimeSteppers::history_max_substeps,
        device_vars->number_of_grid_points(),
        device_dt_variables_tag::type::number_of_independent_components);
    if (device_vars->number_of_grid_points() > 0) {
      ::Kokkos::deep_copy(*device_derivative_history, 0.0);
    }
  }
};

}  // namespace ScalarWave::Actions
