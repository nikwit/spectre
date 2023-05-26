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
struct CheckFunctionsOfTimeAreReady {
  static constexpr size_t Dim = 3;
  using inbox_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& proxy = ::Parallel::get_parallel_component<ParallelComponent>(
        cache)[array_index];
    const std::string function_of_time_name = "Expansion";
    const auto& time = db::get<::Tags::Time>(box);
    bool is_ready =
        Parallel::mutable_cache_item_is_ready<::domain::Tags::FunctionsOfTime>(
            cache,
            [&proxy, &time, &function_of_time_name](
                const std::unordered_map<
                    std::string,
                    std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                    functions_of_time) {
              const auto& f_of_t = functions_of_time.at(function_of_time_name);
              const double expiration_time = f_of_t->time_bounds()[1];
              if (time > expiration_time) {
                return std::unique_ptr<Parallel::Callback>(
                    new Parallel::PerformAlgorithmCallback(proxy));
              }

              return std::unique_ptr<Parallel::Callback>{};
            });
    Parallel::printf(MakeString{}
                     << "Worldtube "
                     << (is_ready ? " functions of time are ready"
                                  : " functions of time are NOT ready")
                     << "\n");

    return {is_ready ? Parallel::AlgorithmExecution::Continue
                     : Parallel::AlgorithmExecution::Retry,
            std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
