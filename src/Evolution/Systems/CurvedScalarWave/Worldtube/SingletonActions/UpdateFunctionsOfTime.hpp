// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube::Actions {

/*!
 * \brief Waits for the data from all neighboring elements and changes the slab
 * size if a change in the global time step is detected.
 * \details We check the slab size of the time step id sent by the elements. If
 * this is different from the slab size currently used by the worldtube
 * singleton, we assume a global slab size change has occurred in the elements
 * and adjust the worldtube slab size accordingly.
 */
struct UpdateFunctionsOfTime {
  static constexpr size_t Dim = 3;
  using inbox_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& time = db::get<::Tags::Time>(box);
    const auto& time_step = db::get<::Tags::TimeStep>(box);
    const auto& functions_of_time =
        Parallel::get<::domain::Tags::FunctionsOfTime>(cache);
    const std::string function_of_time_name = "Expansion";
    const auto& function_of_time = functions_of_time.at(function_of_time_name);
    const double current_fot_expiration_time =
        function_of_time->time_bounds()[1];
    const double new_fot_expiration_time = time + time_step.value() * 0.01;
    const double period = 20.;
    const double amplitude = 0.1;
    const double value = amplitude * sin(time / (2. * M_PI * period));
    const DataVector new_derivative(1, value);
    Parallel::printf(MakeString{} << get_output(time) << ", step "
                                  << get_output(time_step.value()) << "\n");

    if (time > current_fot_expiration_time) {
      Parallel::printf(MakeString{} << "Mutating Time from "
                                    << current_fot_expiration_time << " to "
                                    << new_fot_expiration_time << "\n");
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, function_of_time_name, current_fot_expiration_time,
          new_derivative, new_fot_expiration_time);
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
