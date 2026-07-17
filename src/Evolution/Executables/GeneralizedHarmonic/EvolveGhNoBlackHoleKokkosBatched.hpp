// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

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
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Protocols/HorizonMetavars.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/WhenToCheck.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

template <class Metavariables>
struct GhKokkosBatchedDriver {
  using chare_type = Parallel::Algorithms::Singleton;
  static constexpr bool checkpoint_data = true;
  using metavariables = Metavariables;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<
              Initialization::Actions::InitializeItems<
                  Initialization::TimeStepping<
                      Metavariables, typename Metavariables::TimeStepperBase>,
                  gh::Batched::Initialization::DriverState>,
              Actions::MutateApply<
                  gh::Batched::Actions::InitializeDriver<Metavariables>>,
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
              Actions::MutateApply<gh::Batched::Actions::
                                       UpdateFrameFieldsBatched<Metavariables>>,
              Actions::MutateApply<
                  gh::Batched::Actions::UpdateConstraintDampingBatched<
                      Metavariables>>,
              evolution::Actions::RunEventsAndTriggers<
                  Metavariables::local_time_stepping>,
              Actions::MutateApply<
                  gh::Batched::Actions::ComputeTimeDerivativeBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::PackageLocalFacesBatched>,
              // External TimeDerivative-type boundary conditions consume the
              // volume RHS before internal-face lifting, matching the
              // production executables where they are applied within
              // ComputeTimeDerivative.
              Actions::MutateApply<
                  gh::Batched::Actions::
                      ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::ComputeInternalBoundaryTermsBatched>,
              Actions::MutateApply<
                  gh::Batched::Actions::LiftInternalBoundaryTermsBatched>,
              Actions::MutateApply<
                  evolution::Actions::Kokkos::RecordTimeStepperData<
                      typename Metavariables::system>>,
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
    Parallel::get_parallel_component<GhKokkosBatchedDriver<Metavariables>>(
        local_cache)
        .start_phase(next_phase);
  }
};

struct EvolutionMetavarsKokkosBatched {
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

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<Event, tmpl::flatten<tmpl::list<
                              gh::Events::CompletionBatched,
                              gh::Events::FindApparentHorizonBatched<AhA>,
                              gh::Events::FindApparentHorizonBatched<AhB>,
                              gh::Events::ObserveNormsBatched>>>,
        tmpl::pair<
            gh::BoundaryConditions::BoundaryCondition<volume_dim>,
            gh::BoundaryConditions::standard_boundary_conditions<volume_dim>>,
        tmpl::pair<gh::gauges::GaugeCondition,
                   tmpl::list<gh::gauges::Harmonic>>,
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

  using const_global_cache_tags =
      tmpl::list<gh::gauges::Tags::GaugeCondition,
                 evolution::initial_data::Tags::InitialData,
                 gh::Tags::DampingFunctionGamma0<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma1<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma2<volume_dim, Frame::Grid>,
                 Filters::Tags::Filter<Filters::Exponential<0>>,
                 ::Tags::EventsAndTriggers<Triggers::WhenToCheck::AtSlabs>>;
  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize>;

  using batched_driver_component =
      GhKokkosBatchedDriver<EvolutionMetavarsKokkosBatched>;

  using component_list = tmpl::flatten<
      tmpl::list<observers::ObserverWriter<EvolutionMetavarsKokkosBatched>,
                 ah::Component<EvolutionMetavarsKokkosBatched, AhA>,
                 ah::Component<EvolutionMetavarsKokkosBatched, AhB>,
                 batched_driver_component>>;

  static constexpr Options::String help{
      "Generalized Harmonic Kokkos batched scaffold executable.\n"
      "Initializes packed element metadata and evolves packed elements with "
      "batched Kokkos actions."};

  static constexpr auto default_phase_order = std::array<Parallel::Phase, 4>{
      Parallel::Phase::Initialization, Parallel::Phase::Register,
      Parallel::Phase::Evolve, Parallel::Phase::Exit};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
