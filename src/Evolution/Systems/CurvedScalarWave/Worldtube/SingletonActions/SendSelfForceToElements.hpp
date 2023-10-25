// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube::Actions {
/*!
 * \brief Sends the regular field to each element abutting the worldtube,
 * evaluated at the grid coordinates of each face
 *
 * \details h-refinement could be accounted for by sending the
 * coefficients of the internal solution directly and have each element evaluate
 * it for themselves.
 */
template <typename Metavariables>
struct SendSelfForceToElements {
  static constexpr size_t Dim = Metavariables::volume_dim;
  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    auto& element_proxies = Parallel::get_parallel_component<
        typename Metavariables::dg_element_array>(cache);
    const auto& faces_grid_coords =
        db::get<Tags::ElementFacesGridCoordinates<Dim>>(box);

    for (const auto& [element_id, _] : faces_grid_coords) {
      Scalar<DataVector> self_force = db::get<Tags::SelfForce>(box);
      Parallel::receive_data<Tags::SelfForceInbox<Dim>>(
          element_proxies[element_id], db::get<::Tags::TimeStepId>(box),
          std::move(self_force));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
