// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <pup.h>
#include <string>
#include <utility>

#include "DataStructures/LinkedMessageId.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/HostElementView.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Component.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Events {

/*!
 * \brief Starts an apparent-horizon find from packed GH Kokkos batched data.
 *
 * This event provides the same volume payload that
 * `ah::Events::FindApparentHorizon` sends from normal DG elements, but obtains
 * that payload by copying the packed batched device state to host and slicing
 * it into element-local views.
 */
template <typename HorizonMetavars>
class FindApparentHorizonBatched : public Event {
 public:
  explicit FindApparentHorizonBatched(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(FindApparentHorizonBatched);  // NOLINT

  using options = tmpl::list<>;
  static constexpr Options::String help =
      "Start an apparent-horizon find from packed GH Kokkos batched data.";

  static std::string name() {
    return "FindApparentHorizonBatched" + HorizonMetavars::name();
  }

  FindApparentHorizonBatched() = default;

  explicit FindApparentHorizonBatched(std::optional<std::string> dependency)
      : dependency_(std::move(dependency)) {}

  using compute_tags_for_observation_box =
      typename HorizonMetavars::compute_tags_on_element;
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<
      evolution::Kokkos::Tags::PackedTopology<gh::System<3>>,
      evolution::Kokkos::Tags::PackedGeometry<gh::System<3>>,
      evolution::Kokkos::Tags::PackedEvolutionState<gh::System<3>>,
      typename HorizonMetavars::time_tag>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const evolution::Kokkos::PackedTopology<gh::System<3>>& packed_topology,
      const evolution::Kokkos::PackedGeometry<gh::System<3>>& packed_geometry,
      const evolution::Kokkos::PackedEvolutionState<gh::System<3>>&
          packed_evolution_state,
      const typename HorizonMetavars::time_tag::type& horizon_time,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& /*observation_value*/) const {
    const auto& blocks_to_interpolate =
        Parallel::get<ah::Tags::BlocksForHorizonFind>(cache);
    ASSERT(blocks_to_interpolate.contains(HorizonMetavars::name()),
           "Blocks to interpolate does not contain target "
               << HorizonMetavars::name());
    const auto& blocks_to_interpolate_for_this_target =
        blocks_to_interpolate.at(HorizonMetavars::name());
    const auto& blocks = Parallel::get<domain::Tags::Domain<3>>(cache).blocks();
    auto& horizon_finder_proxy = Parallel::get_parallel_component<
        ah::Component<Metavariables, HorizonMetavars>>(cache);

    gh::Batched::for_each_host_element_view(
        packed_topology, packed_geometry, packed_evolution_state,
        [&blocks, &blocks_to_interpolate_for_this_target, &horizon_finder_proxy,
         &horizon_time,
         this](const gh::Batched::HostElementView& host_element_view) {
          const auto& block_name =
              blocks[host_element_view.element_id.block_id()].name();
          if (not blocks_to_interpolate_for_this_target.contains(block_name)) {
            return;
          }

          Parallel::simple_action<ah::FindApparentHorizon<HorizonMetavars>>(
              horizon_finder_proxy, horizon_time, host_element_view.element_id,
              host_element_view.mesh,
              gh::Batched::make_ah_source_vars(host_element_view), dependency_);
        });
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*component*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  void pup(PUP::er& p) override {
    Event::pup(p);
    p | dependency_;
  }

 private:
  std::optional<std::string> dependency_;
};

template <typename HorizonMetavars>
PUP::able::PUP_ID FindApparentHorizonBatched<HorizonMetavars>::my_PUP_ID =
    0;  // NOLINT

}  // namespace gh::Events
