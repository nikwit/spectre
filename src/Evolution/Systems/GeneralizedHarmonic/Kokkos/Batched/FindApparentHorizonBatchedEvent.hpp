// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ControlSystem/Measurements/CharSpeed.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/HostElementView.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "NumericalAlgorithms/Strahlkorper/Strahlkorper.hpp"
#include "NumericalAlgorithms/Strahlkorper/StrahlkorperFunctions.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Component.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeVarsToInterpolateToTarget.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/Interpolation/ComputeExcisionBoundaryVolumeQuantities.hpp"
#include "ParallelAlgorithms/Interpolation/ComputeExcisionBoundaryVolumeQuantities.tpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
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
  using packed_geometry_tag =
      evolution::Kokkos::Tags::PackedGeometry<gh::System<3>>;
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
  using return_tags = tmpl::list<packed_geometry_tag>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<gh::System<3>>,
                 evolution::Kokkos::Tags::PackedEvolutionState<gh::System<3>>,
                 typename HorizonMetavars::time_tag>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const evolution::Kokkos::PackedTopology<gh::System<3>>& packed_topology,
      const evolution::Kokkos::PackedEvolutionState<gh::System<3>>&
          packed_evolution_state,
      const typename HorizonMetavars::time_tag::type& horizon_time,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& /*observation_value*/) const {
    gh::Actions::synchronize_frame_fields_to_host<3>(packed_geometry);
    const auto& blocks_to_interpolate =
        Parallel::get<ah::Tags::BlocksForHorizonFind>(cache);
    ASSERT(blocks_to_interpolate.contains(HorizonMetavars::name()),
           "Blocks to interpolate does not contain target "
               << HorizonMetavars::name());
    const auto& blocks_to_interpolate_for_this_target =
        blocks_to_interpolate.at(HorizonMetavars::name());
    const auto& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
    const auto& blocks = domain.blocks();
    const auto& functions_of_time =
        Parallel::get<control_system::Tags::FunctionsOfTimeInitialize>(cache);
    auto& horizon_finder_proxy = Parallel::get_parallel_component<
        ah::Component<Metavariables, HorizonMetavars>>(cache);

    gh::Batched::for_each_host_element_view(
        packed_topology, *packed_geometry, packed_evolution_state,
        [&blocks, &blocks_to_interpolate_for_this_target, &domain,
         &functions_of_time, &horizon_finder_proxy, &horizon_time,
         this](const gh::Batched::HostElementView& host_element_view) {
          const auto& block_name =
              blocks[host_element_view.element_id.block_id()].name();
          if (not blocks_to_interpolate_for_this_target.contains(block_name)) {
            return;
          }

          auto source_vars =
              gh::Batched::make_ah_source_vars(host_element_view);
          using horizon_frame = typename HorizonMetavars::frame;
          if constexpr (std::is_same_v<horizon_frame, Frame::Distorted>) {
            const auto& block = blocks[host_element_view.element_id.block_id()];
            if (not block.has_distorted_frame()) {
              return;
            }
          }
          Variables<ah::vars_to_interpolate_to_target<3, horizon_frame>>
              vars_to_interpolate_to_target{
                  host_element_view.mesh.number_of_grid_points()};
          ah::compute_vars_to_interpolate_to_target(
              make_not_null(&vars_to_interpolate_to_target),
              get<gr::Tags::SpacetimeMetric<DataVector, 3>>(source_vars),
              get<gh::Tags::Pi<DataVector, 3>>(source_vars),
              get<gh::Tags::Phi<DataVector, 3>>(source_vars),
              get<::Tags::deriv<gh::Tags::Phi<DataVector, 3>, tmpl::size_t<3>,
                                Frame::Inertial>>(source_vars),
              horizon_time, domain, host_element_view.mesh,
              host_element_view.element_id, functions_of_time);

          Parallel::simple_action<ah::FindApparentHorizon<HorizonMetavars>>(
              horizon_finder_proxy, horizon_time, host_element_view.element_id,
              host_element_view.mesh, std::move(vars_to_interpolate_to_target),
              dependency_);
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

/*!
 * \brief Sends batched GH Kokkos horizon data to a control-system AH component.
 *
 * Analogous to `FindApparentHorizonBatched` but intended for control-system
 * horizon finders that expect a `LinkedMessageId<double>` time stamp.  The
 * linked message id is assembled from the current observation time and
 * `::Tags::PreviousTriggerTime`, which `EventsAndDenseTriggers` populates
 * before events run.
 */
template <typename HorizonMetavars>
class SendHorizonDataBatchedForControlSystem : public Event {
 public:
  using packed_geometry_tag =
      evolution::Kokkos::Tags::PackedGeometry<gh::System<3>>;
  explicit SendHorizonDataBatchedForControlSystem(
      CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(  // NOLINT
      SendHorizonDataBatchedForControlSystem);

  using options = tmpl::list<>;
  static constexpr Options::String help =
      "Send batched GH Kokkos horizon data to a control-system AH component.";

  static std::string name() {
    return "SendHorizonDataBatchedForControlSystem" + HorizonMetavars::name();
  }

  SendHorizonDataBatchedForControlSystem() = default;

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<packed_geometry_tag>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<gh::System<3>>,
                 evolution::Kokkos::Tags::PackedEvolutionState<gh::System<3>>,
                 ::Tags::PreviousTriggerTime>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const evolution::Kokkos::PackedTopology<gh::System<3>>& packed_topology,
      const evolution::Kokkos::PackedEvolutionState<gh::System<3>>&
          packed_evolution_state,
      const std::optional<double>& previous_trigger_time,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value) const {
    gh::Actions::synchronize_frame_fields_to_host<3>(packed_geometry);
    const LinkedMessageId<double> measurement_id{observation_value.value,
                                                 previous_trigger_time};

    const auto& blocks_to_interpolate =
        Parallel::get<ah::Tags::BlocksForHorizonFind>(cache);
    ASSERT(blocks_to_interpolate.contains(HorizonMetavars::name()),
           "Blocks to interpolate does not contain target "
               << HorizonMetavars::name());
    const auto& blocks_to_interpolate_for_this_target =
        blocks_to_interpolate.at(HorizonMetavars::name());
    const auto& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
    const auto& blocks = domain.blocks();
    const auto& functions_of_time =
        Parallel::get<control_system::Tags::FunctionsOfTimeInitialize>(cache);
    auto& horizon_finder_proxy = Parallel::get_parallel_component<
        ah::Component<Metavariables, HorizonMetavars>>(cache);

    gh::Batched::for_each_host_element_view(
        packed_topology, *packed_geometry, packed_evolution_state,
        [&blocks, &blocks_to_interpolate_for_this_target, &domain,
         &functions_of_time, &horizon_finder_proxy, &measurement_id](
            const gh::Batched::HostElementView& host_element_view) {
          const auto& block_name =
              blocks[host_element_view.element_id.block_id()].name();
          if (not blocks_to_interpolate_for_this_target.contains(block_name)) {
            return;
          }

          auto source_vars =
              gh::Batched::make_ah_source_vars(host_element_view);
          using horizon_frame = typename HorizonMetavars::frame;
          if constexpr (std::is_same_v<horizon_frame, Frame::Distorted>) {
            const auto& block = blocks[host_element_view.element_id.block_id()];
            if (not block.has_distorted_frame()) {
              return;
            }
          }
          Variables<ah::vars_to_interpolate_to_target<3, horizon_frame>>
              vars_to_interpolate_to_target{
                  host_element_view.mesh.number_of_grid_points()};
          ah::compute_vars_to_interpolate_to_target(
              make_not_null(&vars_to_interpolate_to_target),
              get<gr::Tags::SpacetimeMetric<DataVector, 3>>(source_vars),
              get<gh::Tags::Pi<DataVector, 3>>(source_vars),
              get<gh::Tags::Phi<DataVector, 3>>(source_vars),
              get<::Tags::deriv<gh::Tags::Phi<DataVector, 3>, tmpl::size_t<3>,
                                Frame::Inertial>>(source_vars),
              measurement_id, domain, host_element_view.mesh,
              host_element_view.element_id, functions_of_time);

          Parallel::simple_action<ah::FindApparentHorizon<HorizonMetavars>>(
              horizon_finder_proxy, measurement_id,
              host_element_view.element_id, host_element_view.mesh,
              std::move(vars_to_interpolate_to_target),
              std::optional<std::string>{});
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

  void pup(PUP::er& p) override { Event::pup(p); }
};

template <typename HorizonMetavars>
PUP::able::PUP_ID
    SendHorizonDataBatchedForControlSystem<HorizonMetavars>::my_PUP_ID =
        0;  // NOLINT

template <::domain::ObjectLabel Object, typename ControlSystems>
class SendExcisionDataBatchedForSizeControl : public Event {
 private:
  using char_speed_measurement =
      control_system::measurements::CharSpeed<Object>;
  using submeasurement = typename char_speed_measurement::Excision;
  using interpolation_target_tag =
      typename submeasurement::template interpolation_target_tag<
          ControlSystems>;
  using vars_to_interpolate =
      typename interpolation_target_tag::vars_to_interpolate_to_target;

  static ylm::Strahlkorper<Frame::Grid> make_excision_surface(
      const intrp::OptionHolders::Sphere& sphere) {
    ASSERT(sphere.radii.size() == 1,
           "Size-control excision measurements require exactly one radius for "
               << pretty_type::name<interpolation_target_tag>() << ".");
    ASSERT(sphere.angular_ordering == ylm::AngularOrdering::Strahlkorper,
           "Size-control excision measurements require Strahlkorper angular "
           "ordering for "
               << pretty_type::name<interpolation_target_tag>() << ".");
    return ylm::Strahlkorper<Frame::Grid>{sphere.l_max, *sphere.radii.begin(),
                                          sphere.center};
  }

  static Variables<vars_to_interpolate> compute_volume_excision_vars(
      const gh::Batched::HostElementView& host_element_view) {
    Variables<vars_to_interpolate> result{
        host_element_view.mesh.number_of_grid_points()};
    auto source_vars = gh::Batched::make_ah_source_vars(host_element_view);
    Jacobian<DataVector, 3, Frame::ElementLogical, Frame::Distorted>
        unused_logical_to_distorted_jacobian{
            host_element_view.mesh.number_of_grid_points()};
    intrp::ComputeExcisionBoundaryVolumeQuantities::apply(
        make_not_null(&result), source_vars, host_element_view.mesh,
        host_element_view.distorted_to_inertial_jacobian,
        host_element_view.distorted_to_inertial_inverse_jacobian,
        unused_logical_to_distorted_jacobian,
        host_element_view.logical_to_distorted_inverse_jacobian,
        host_element_view.grid_to_distorted_inverse_jacobian,
        host_element_view.mesh_velocity,
        host_element_view.grid_to_distorted_frame_velocity);
    return result;
  }

  static std::vector<BlockLogicalCoords<3>> block_logical_coords_for_element(
      const tnsr::I<DataVector, 3, Frame::Grid>& target_points,
      const Domain<3>& domain,
      const domain::FunctionsOfTimeMap& functions_of_time,
      const ElementId<3>& element_id, const double time) {
    std::vector<BlockLogicalCoords<3>> result{get<0>(target_points).size()};
    const auto& block = domain.blocks()[element_id.block_id()];
    if (not block.has_distorted_frame()) {
      return result;
    }
    for (size_t i = 0; i < get<0>(target_points).size(); ++i) {
      tnsr::I<double, 3, Frame::Grid> target_point{};
      for (size_t d = 0; d < 3; ++d) {
        target_point.get(d) = target_points.get(d)[i];
      }
      const auto block_logical_coordinate =
          block_logical_coordinates_single_point(target_point, block, time,
                                                 functions_of_time);
      if (block_logical_coordinate.has_value()) {
        result[i] = make_id_pair(domain::BlockId(element_id.block_id()),
                                 block_logical_coordinate.value());
      }
    }
    return result;
  }

  static void copy_interpolated_points(
      const gsl::not_null<Variables<vars_to_interpolate>*>
          interpolated_on_target,
      const Variables<vars_to_interpolate>& interpolated_from_element,
      const std::vector<size_t>& offsets,
      const gsl::not_null<std::vector<bool>*> filled_points) {
    for (size_t component = 0;
         component <
         Variables<vars_to_interpolate>::number_of_independent_components;
         ++component) {
      for (size_t i = 0; i < offsets.size(); ++i) {
        interpolated_on_target
            ->data()[component *
                         interpolated_on_target->number_of_grid_points() +
                     offsets[i]] =
            interpolated_from_element.data()[component * offsets.size() + i];
      }
    }
    for (const size_t offset : offsets) {
      (*filled_points)[offset] = true;
    }
  }

 public:
  using packed_geometry_tag =
      evolution::Kokkos::Tags::PackedGeometry<gh::System<3>>;
  explicit SendExcisionDataBatchedForSizeControl(CkMigrateMessage* /*unused*/) {
  }
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(  // NOLINT
      SendExcisionDataBatchedForSizeControl);

  using options = tmpl::list<>;
  static constexpr Options::String help =
      "Send batched GH Kokkos excision-surface data to a Size control system.";

  static std::string name() {
    return "SendExcisionDataBatchedForSizeControl" +
           pretty_type::name<interpolation_target_tag>();
  }

  SendExcisionDataBatchedForSizeControl() = default;

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<packed_geometry_tag>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<gh::System<3>>,
                 evolution::Kokkos::Tags::PackedEvolutionState<gh::System<3>>,
                 ::Tags::PreviousTriggerTime>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const evolution::Kokkos::PackedTopology<gh::System<3>>& packed_topology,
      const evolution::Kokkos::PackedEvolutionState<gh::System<3>>&
          packed_evolution_state,
      const std::optional<double>& previous_trigger_time,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value) const {
    gh::Actions::synchronize_frame_fields_to_host<3>(packed_geometry);
    const LinkedMessageId<double> measurement_id{observation_value.value,
                                                 previous_trigger_time};
    const auto& sphere =
        Parallel::get<intrp::Tags::Sphere<interpolation_target_tag>>(cache);
    const auto excision_surface = make_excision_surface(sphere);
    const auto target_points = ylm::cartesian_coords(excision_surface);
    Variables<vars_to_interpolate> interpolated_on_target{
        get<0>(target_points).size()};
    std::fill(interpolated_on_target.data(),
              interpolated_on_target.data() + interpolated_on_target.size(),
              std::numeric_limits<double>::signaling_NaN());
    std::vector<bool> filled_points(get<0>(target_points).size(), false);

    const auto& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
    const auto& functions_of_time =
        Parallel::get<domain::Tags::FunctionsOfTime>(cache);

    gh::Batched::for_each_host_element_view(
        packed_topology, *packed_geometry, packed_evolution_state,
        [&domain, &filled_points, &functions_of_time, &interpolated_on_target,
         &target_points, time = measurement_id.id](
            const gh::Batched::HostElementView& host_element_view) {
          const auto block_logical_coords = block_logical_coords_for_element(
              target_points, domain, functions_of_time,
              host_element_view.element_id, time);
          const std::vector<ElementId<3>> element_ids{
              {host_element_view.element_id}};
          const auto element_coord_holders =
              element_logical_coordinates(element_ids, block_logical_coords);
          if (not element_coord_holders.contains(
                  host_element_view.element_id)) {
            return;
          }

          const auto& element_coord_holder =
              element_coord_holders.at(host_element_view.element_id);
          auto volume_excision_vars =
              compute_volume_excision_vars(host_element_view);
          const intrp::Irregular<3> interpolator{
              host_element_view.mesh,
              element_coord_holder.element_logical_coords};
          copy_interpolated_points(
              make_not_null(&interpolated_on_target),
              interpolator.interpolate(volume_excision_vars),
              element_coord_holder.offsets, make_not_null(&filled_points));
        });

    ASSERT(std::all_of(filled_points.begin(), filled_points.end(),
                       [](const bool filled) { return filled; }),
           "Failed to interpolate all Size-control excision points for "
               << pretty_type::name<interpolation_target_tag>() << ".");

    auto callback_box = db::create<
        db::AddSimpleTags<ylm::Tags::Strahlkorper<Frame::Grid>,
                          ::Tags::Variables<vars_to_interpolate>>,
        db::AddComputeTags<gr::Tags::DetAndInverseSpatialMetricCompute<
            DataVector, 3, Frame::Distorted>>>(
        excision_surface, std::move(interpolated_on_target));
    control_system::RunCallbacks<submeasurement, ControlSystems>::apply(
        callback_box, cache, measurement_id);
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*component*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  void pup(PUP::er& p) override { Event::pup(p); }
};

template <::domain::ObjectLabel Object, typename ControlSystems>
PUP::able::PUP_ID
    SendExcisionDataBatchedForSizeControl<Object, ControlSystems>::my_PUP_ID =
        0;  // NOLINT

}  // namespace gh::Events
