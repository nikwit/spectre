// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "ControlSystem/Actions/LimitTimeStep.hpp"
#include "ControlSystem/CombinedName.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Size/Factory.hpp"
#include "ControlSystem/ControlErrors/Size/State.hpp"
#include "ControlSystem/FutureMeasurements.hpp"
#include "ControlSystem/Measurements/BothHorizons.hpp"
#include "ControlSystem/Measurements/CharSpeed.hpp"
#include "ControlSystem/Metafunctions.hpp"
#include "ControlSystem/Systems/Expansion.hpp"
#include "ControlSystem/Systems/Rotation.hpp"
#include "ControlSystem/Systems/Shape.hpp"
#include "ControlSystem/Systems/Size.hpp"
#include "ControlSystem/Systems/Skew.hpp"
#include "ControlSystem/Systems/Translation.hpp"
#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "ControlSystem/Tags/FutureMeasurements.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "ControlSystem/Trigger.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Creators/OptionTags.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/Executables/GeneralizedHarmonic/Batched/Actions.hpp"
#include "Evolution/Executables/GeneralizedHarmonic/GeneralizedHarmonicBase.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Kokkos/CleanHistory.hpp"
#include "Evolution/Kokkos/RecordTimeStepperData.hpp"
#include "Evolution/Kokkos/UpdateU.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/FindApparentHorizonBatchedEvent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ObserveTimeStepBatchedEvent.hpp"
#include "Evolution/Tags/Filter.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/LinearOperators/ExponentialFilter.hpp"
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FailedHorizonFind.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ObserveFieldsOnHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ObserveTimeSeriesOnHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Component.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Destination.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Protocols/HorizonMetavars.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/WhenToCheck.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gh::Batched::Initialization {
template <typename Metavariables>
struct ControlSystemMeasurements {
 private:
  using control_system_group =
      control_system::metafunctions::control_systems_with_measurement_t<
          typename Metavariables::control_systems,
          control_system::measurements::BothHorizons>;
  using horizon_metavars =
      control_system::metafunctions::horizon_metavars<control_system_group>;
  using future_measurements_tag =
      control_system::Tags::FutureMeasurements<control_system_group>;

 public:
  using simple_tags =
      tmpl::list<future_measurements_tag, ::Tags::EventsAndDenseTriggers,
                 ::Tags::PreviousTriggerTime>;
  using simple_tags_from_options = tmpl::list<>;
  using compute_tags = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementsPerUpdate>;
  using mutable_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementTimescales>;
  using return_tags = tmpl::push_back<simple_tags, ::Tags::TimeStep>;
  using argument_tags = tmpl::list<::Tags::Time, ::Tags::TimeStep,
                                   Parallel::Tags::GlobalCache<Metavariables>>;

  static void apply(
      const gsl::not_null<control_system::FutureMeasurements*>
          future_measurements,
      const gsl::not_null<EventsAndDenseTriggers*> events_and_dense_triggers,
      const gsl::not_null<std::optional<double>*> previous_trigger_time,
      const gsl::not_null<TimeDelta*> time_step, const double initial_time,
      const TimeDelta& initial_time_step,
      const Parallel::GlobalCache<Metavariables>* const& cache) {
    *time_step = initial_time_step;
    *previous_trigger_time = std::nullopt;

    const auto& measurement_timescales =
        Parallel::get<control_system::Tags::MeasurementTimescales>(*cache);
    const auto& measurement_timescale = *measurement_timescales.at(
        control_system::combined_name<control_system_group>());
    const double initial_measurement_timescale =
        measurement_timescale.func(initial_time)[0][0];
    const bool active = initial_measurement_timescale !=
                        std::numeric_limits<double>::infinity();
    *future_measurements = control_system::FutureMeasurements{
        static_cast<size_t>(
            active ? Parallel::get<control_system::Tags::MeasurementsPerUpdate>(
                         *cache)
                   : 1),
        active ? initial_time : std::numeric_limits<double>::infinity()};
    if (active) {
      future_measurements->update(measurement_timescale);
    }

    // Wire up the dense trigger with one batched event per control-system
    // horizon metavar.
    std::vector<std::unique_ptr<Event>> events;
    tmpl::for_each<horizon_metavars>([&events](auto hmv) {
      using HM = tmpl::type_from<decltype(hmv)>;
      events.push_back(
          std::make_unique<
              gh::Events::SendHorizonDataBatchedForControlSystem<HM>>());
    });
    events_and_dense_triggers->add_trigger_and_events(
        std::make_unique<control_system::Trigger<control_system_group>>(),
        std::move(events));

    double earliest_expiration = std::numeric_limits<double>::infinity();
    for (const auto& function_of_time :
         Parallel::get<control_system::Tags::FunctionsOfTimeInitialize>(
             *cache)) {
      earliest_expiration = std::min(earliest_expiration,
                                     function_of_time.second->time_bounds()[1]);
    }
    const auto start_time = time_step->slab().start();
    if ((start_time + *time_step).value() > earliest_expiration) {
      *time_step = Slab(start_time.value(), earliest_expiration).duration();
    }
  }
};

template <typename Metavariables, ::domain::ObjectLabel Object>
struct SizeControlMeasurements {
 private:
  using measurement = control_system::measurements::CharSpeed<Object>;
  using control_system_group =
      control_system::metafunctions::control_systems_with_measurement_t<
          typename Metavariables::control_systems, measurement>;
  using future_measurements_tag =
      control_system::Tags::FutureMeasurements<control_system_group>;
  using horizon_metavars =
      control_system::metafunctions::horizon_metavars<control_system_group>;

 public:
  using simple_tags = tmpl::list<future_measurements_tag>;
  using simple_tags_from_options = tmpl::list<>;
  using compute_tags = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementsPerUpdate>;
  using mutable_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementTimescales>;
  using return_tags =
      tmpl::push_back<simple_tags, ::Tags::EventsAndDenseTriggers,
                      ::Tags::TimeStep>;
  using argument_tags = tmpl::list<::Tags::Time, ::Tags::TimeStep,
                                   Parallel::Tags::GlobalCache<Metavariables>>;

  static void apply(
      const gsl::not_null<control_system::FutureMeasurements*>
          future_measurements,
      const gsl::not_null<EventsAndDenseTriggers*> events_and_dense_triggers,
      const gsl::not_null<TimeDelta*> time_step, const double initial_time,
      const TimeDelta& initial_time_step,
      const Parallel::GlobalCache<Metavariables>* const& cache) {
    *time_step = initial_time_step;

    const auto& measurement_timescales =
        Parallel::get<control_system::Tags::MeasurementTimescales>(*cache);
    const auto& measurement_timescale = *measurement_timescales.at(
        control_system::combined_name<control_system_group>());
    const double initial_measurement_timescale =
        measurement_timescale.func(initial_time)[0][0];
    const bool active = initial_measurement_timescale !=
                        std::numeric_limits<double>::infinity();
    *future_measurements = control_system::FutureMeasurements{
        static_cast<size_t>(
            active ? Parallel::get<control_system::Tags::MeasurementsPerUpdate>(
                         *cache)
                   : 1),
        active ? initial_time : std::numeric_limits<double>::infinity()};
    if (active) {
      future_measurements->update(measurement_timescale);
    }

    std::vector<std::unique_ptr<Event>> events;
    tmpl::for_each<horizon_metavars>([&events](auto hmv) {
      using HM = tmpl::type_from<decltype(hmv)>;
      events.push_back(
          std::make_unique<
              gh::Events::SendHorizonDataBatchedForControlSystem<HM>>());
    });
    events.push_back(
        std::make_unique<gh::Events::SendExcisionDataBatchedForSizeControl<
            Object, control_system_group>>());
    events_and_dense_triggers->add_trigger_and_events(
        std::make_unique<control_system::Trigger<control_system_group>>(),
        std::move(events));

    double earliest_expiration = std::numeric_limits<double>::infinity();
    for (const auto& function_of_time :
         Parallel::get<control_system::Tags::FunctionsOfTimeInitialize>(
             *cache)) {
      earliest_expiration = std::min(earliest_expiration,
                                     function_of_time.second->time_bounds()[1]);
    }
    const auto start_time = time_step->slab().start();
    if ((start_time + *time_step).value() > earliest_expiration) {
      *time_step = Slab(start_time.value(), earliest_expiration).duration();
    }
  }
};
}  // namespace gh::Batched::Initialization

template <class Metavariables>
struct GhBinaryBlackHoleKokkosDriver {
  using chare_type = Parallel::Algorithms::Singleton;
  static constexpr bool checkpoint_data = true;
  using metavariables = Metavariables;
  using functions_of_time_tag = control_system::Tags::FunctionsOfTimeInitialize;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<
              Initialization::Actions::InitializeItems<
                  Initialization::TimeStepping<
                      Metavariables, typename Metavariables::TimeStepperBase>,
                  gh::Batched::Initialization::DriverState,
                  gh::Batched::Initialization::ControlSystemMeasurements<
                      Metavariables>,
                  gh::Batched::Initialization::SizeControlMeasurements<
                      Metavariables, ::domain::ObjectLabel::A>,
                  gh::Batched::Initialization::SizeControlMeasurements<
                      Metavariables, ::domain::ObjectLabel::B>>,
              Actions::MutateApply<gh::Batched::Actions::InitializeDriver<
                  Metavariables, functions_of_time_tag>>,
              Actions::MutateApply<
                  gh::Batched::Actions::InitializeBoundaryBatchMetadata>,
              // Initialization::TimeStepping leaves Tags::TimeStepId
              // undefined; production executables advance it once in the
              // self-start phase before evolving. Runge-Kutta steppers need
              // no self-start history, so a single AdvanceTime supplies the
              // initial TimeStepId here.
              Actions::AdvanceTime, Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<Parallel::Phase::Register,
                             tmpl::list<Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Evolve,
          tmpl::list<
              Actions::MutateApply<
                  gh::Batched::Actions::UpdateFrameFieldsBatched<
                      Metavariables, functions_of_time_tag>>,
              Actions::MutateApply<
                  gh::Batched::Actions::UpdateConstraintDampingBatched<
                      Metavariables, functions_of_time_tag>>,
              evolution::Actions::RunEventsAndTriggers<
                  Metavariables::local_time_stepping>,
              Actions::MutateApply<
                  gh::Batched::Actions::ComputeTimeDerivativeBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::PackageLocalFacesBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::ComputeInternalBoundaryTermsBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::LiftInternalBoundaryTermsBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::
                      ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched>,
              Actions::MutateApply<
                  evolution::Actions::Kokkos::RecordTimeStepperData<
                      typename Metavariables::system>>,
              gh::Batched::Actions::RunEventsAndDenseTriggersBatched<
                  functions_of_time_tag>,
              control_system::Actions::LimitTimeStep<
                  typename Metavariables::control_systems>,
              Actions::MutateApply<evolution::Actions::Kokkos::UpdateU<
                  typename Metavariables::system>>,
              Actions::MutateApply<gh::Batched::Actions::FilterKokkosBatched<
                  Filters::Exponential<0>>>,
              Actions::MutateApply<evolution::Actions::Kokkos::CleanHistory<
                  typename Metavariables::system>>,
              Actions::AdvanceTime>>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<
        GhBinaryBlackHoleKokkosDriver<Metavariables>>(local_cache)
        .start_phase(next_phase);
  }
};

struct EvolutionMetavarsBinaryBlackHoleKokkos {
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using temporal_id = Tags::TimeStepId;
  using TimeStepperBase = TimeStepper;
  static constexpr bool local_time_stepping = false;
  static constexpr bool use_dg_element_collection = false;

  template <::domain::ObjectLabel Horizon>
  struct Ah : tt::ConformsTo<ah::protocols::HorizonMetavars> {
    static constexpr size_t index = 10 + static_cast<size_t>(Horizon);

    using time_tag = ah::Tags::ObservationTime<index>;
    using frame = ::Frame::Inertial;
    using horizon_find_callbacks =
        tmpl::list<ah::callbacks::ObserveFieldsOnHorizon<
                       ::ah::surface_tags_for_observing, Ah>,
                   ah::callbacks::ObserveTimeSeriesOnHorizon<
                       ::ah::tags_for_observing<frame>, Ah>>;
    using horizon_find_failure_callbacks =
        tmpl::list<ah::callbacks::FailedHorizonFind<Ah, false>>;
    using compute_tags_on_element =
        tmpl::list<ah::Tags::ObservationTimeCompute<index>>;

    static constexpr ah::Destination destination = ah::Destination::Observation;

    static std::string name() {
      return "BatchedObservationAh" + ::domain::name(Horizon);
    }
  };

  using AhA = Ah<::domain::ObjectLabel::A>;
  using AhB = Ah<::domain::ObjectLabel::B>;
  using both_horizons = control_system::measurements::BothHorizons;
  using control_systems =
      tmpl::list<control_system::Systems::Rotation<3, both_horizons>,
                 control_system::Systems::Expansion<2, both_horizons>,
                 control_system::Systems::Translation<2, both_horizons, 2>,
                 control_system::Systems::Skew<2, both_horizons>,
                 control_system::Systems::Shape<::domain::ObjectLabel::A, 2,
                                                both_horizons>,
                 control_system::Systems::Shape<::domain::ObjectLabel::B, 2,
                                                both_horizons>,
                 control_system::Systems::Size<::domain::ObjectLabel::A, 2>,
                 control_system::Systems::Size<::domain::ObjectLabel::B, 2>>;
  static constexpr bool use_control_systems =
      tmpl::size<control_systems>::value > 0;
  using SizeA = control_system::Systems::Size<::domain::ObjectLabel::A, 2>;
  using SizeB = control_system::Systems::Size<::domain::ObjectLabel::B, 2>;
  using size_control_systems = tmpl::list<SizeA, SizeB>;
  using size_excision_target_tags =
      control_system::metafunctions::interpolation_target_tags<
          size_control_systems>;
  using control_system_horizon_metavars =
      control_system::metafunctions::horizon_metavars<control_systems>;
  using char_speed_a_control_systems =
      control_system::metafunctions::control_systems_with_measurement_t<
          control_systems,
          control_system::measurements::CharSpeed<::domain::ObjectLabel::A>>;
  using char_speed_b_control_systems =
      control_system::metafunctions::control_systems_with_measurement_t<
          control_systems,
          control_system::measurements::CharSpeed<::domain::ObjectLabel::B>>;
  using control_components =
      control_system::control_components<EvolutionMetavarsBinaryBlackHoleKokkos,
                                         control_systems>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<DenseTrigger,
                   control_system::control_system_triggers<control_systems>>,
        tmpl::pair<
            Event,
            tmpl::flatten<tmpl::list<
                gh::Events::CompletionBatched,
                gh::Events::FindApparentHorizonBatched<AhA>,
                gh::Events::FindApparentHorizonBatched<AhB>,
                gh::Events::ObserveNormsBatched,
                gh::Events::SendExcisionDataBatchedForSizeControl<
                    ::domain::ObjectLabel::A, char_speed_a_control_systems>,
                gh::Events::SendExcisionDataBatchedForSizeControl<
                    ::domain::ObjectLabel::B, char_speed_b_control_systems>,
                tmpl::transform<
                    control_system_horizon_metavars,
                    tmpl::bind<
                        gh::Events::SendHorizonDataBatchedForControlSystem,
                        tmpl::_1>>>>>,
        tmpl::pair<
            gh::BoundaryConditions::BoundaryCondition<volume_dim>,
            gh::BoundaryConditions::standard_boundary_conditions<volume_dim>>,
        tmpl::pair<gh::gauges::GaugeCondition,
                   tmpl::list<gh::gauges::Harmonic>>,
        tmpl::pair<control_system::size::State,
                   control_system::size::States::factory_creatable_states>,
        tmpl::pair<evolution::initial_data::InitialData,
                   tmpl::append<gh::Solutions::all_solutions<volume_dim>,
                                tmpl::list<gh::NumericInitialData>>>,
        tmpl::pair<MathFunction<1, Frame::Inertial>,
                   MathFunctions::all_math_functions<1, Frame::Inertial>>,
        tmpl::pair<TimeSequence<double>,
                   TimeSequences::all_time_sequences<double>>,
        tmpl::pair<TimeSequence<std::uint64_t>,
                   TimeSequences::all_time_sequences<std::uint64_t>>,
        tmpl::pair<TimeStepper, TimeSteppers::time_steppers>,
        tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                         Triggers::time_triggers>>>;
  };

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::push_back<
          tmpl::at<typename factory_creation::factory_classes, Event>>>;

  using const_global_cache_tags = tmpl::flatten<tmpl::list<
      tmpl::list<gh::gauges::Tags::GaugeCondition,
                 evolution::initial_data::Tags::InitialData,
                 gh::Tags::DampingFunctionGamma0<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma1<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma2<volume_dim, Frame::Grid>,
                 Filters::Tags::Filter<Filters::Exponential<0>>,
                 ::Tags::EventsAndTriggers<Triggers::WhenToCheck::AtSlabs>>,
      tmpl::transform<size_excision_target_tags,
                      tmpl::bind<intrp::Tags::Sphere, tmpl::_1>>>>;
  using mutable_global_cache_tags =
      tmpl::list<control_system::Tags::FunctionsOfTimeInitialize>;

  using batched_driver_component =
      GhBinaryBlackHoleKokkosDriver<EvolutionMetavarsBinaryBlackHoleKokkos>;

  using component_list = tmpl::flatten<tmpl::list<
      observers::ObserverWriter<EvolutionMetavarsBinaryBlackHoleKokkos>,
      ah::Component<EvolutionMetavarsBinaryBlackHoleKokkos, AhA>,
      ah::Component<EvolutionMetavarsBinaryBlackHoleKokkos, AhB>,
      tmpl::transform<
          control_system_horizon_metavars,
          tmpl::bind<ah::Component,
                     tmpl::pin<EvolutionMetavarsBinaryBlackHoleKokkos>,
                     tmpl::_1>>,
      control_components, batched_driver_component>>;

  static constexpr Options::String help{
      "Generalized Harmonic Kokkos batched binary-black-hole executable.\n"
      "Evolves packed elements with batched Kokkos actions and drives the "
      "two-horizon control systems that only require apparent-horizon data."};

  static constexpr auto default_phase_order = std::array<Parallel::Phase, 4>{
      Parallel::Phase::Initialization, Parallel::Phase::Register,
      Parallel::Phase::Evolve, Parallel::Phase::Exit};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
