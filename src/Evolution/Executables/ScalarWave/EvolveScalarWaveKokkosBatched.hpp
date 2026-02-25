// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Actions.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/Batched/ObserveTimeStepBatchedEvent.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/WhenToCheck.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/Factory.hpp"
#include "PointwiseFunctions/InitialDataUtilities/NumericData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "PointwiseFunctions/MathFunctions/Factory.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Inertial;
}  // namespace Frame
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

template <class Metavariables>
struct ScalarWaveKokkosBatchedDriver {
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
                  ScalarWave::Batched::Initialization::DriverState>,
              Actions::MutateApply<ScalarWave::Batched::Actions::
                                       InitializeDriver<Metavariables>>,
              Actions::MutateApply<ScalarWave::Batched::Actions::
                                       InitializeBoundaryBatchMetadata>,
              Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<Parallel::Phase::Register,
                             tmpl::list<Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Evolve,
          tmpl::list<
              evolution::Actions::RunEventsAndTriggers<
                  Metavariables::local_time_stepping>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::ComputeTimeDerivativeBatched>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::PackageLocalFacesBatched>,
              Actions::MutateApply<ScalarWave::Batched::Actions::
                                       ComputeInternalBoundaryTermsBatched>,
              Actions::MutateApply<ScalarWave::Batched::Actions::
                                       LiftInternalBoundaryTermsBatched>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::
                      ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::RecordTimeStepperDataBatched>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::UpdateUBatched>,
              Actions::MutateApply<
                  ScalarWave::Batched::Actions::CleanHistoryBatched>,
              Actions::AdvanceTime>>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<
        ScalarWaveKokkosBatchedDriver<Metavariables>>(local_cache)
        .start_phase(next_phase);
  }
};

struct EvolutionMetavarsKokkosBatched {
  static constexpr size_t volume_dim = 3;

  using system = ScalarWave::System<volume_dim>;
  using initial_data_list = ScalarWave::Solutions::all_solutions<volume_dim>;
  using temporal_id = Tags::TimeStepId;
  using TimeStepperBase = TimeStepper;
  static constexpr bool local_time_stepping = false;
  static constexpr bool use_dg_element_collection = false;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<Event, tmpl::flatten<tmpl::list<
                              ScalarWave::Events::CompletionBatched,
                              ScalarWave::Events::ObserveNormsBatched>>>,
        tmpl::pair<evolution::initial_data::InitialData,
                   tmpl::push_back<initial_data_list,
                                   evolution::initial_data::NumericData>>,
        tmpl::pair<
            ScalarWave::BoundaryConditions::BoundaryCondition<volume_dim>,
            ScalarWave::BoundaryConditions::standard_boundary_conditions<
                volume_dim>>,
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
      observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
          tmpl::at<typename factory_creation::factory_classes, Event>>>>;
  using const_global_cache_tags =
      tmpl::list<evolution::initial_data::Tags::InitialData,
                 domain::Tags::ExternalBoundaryConditions<volume_dim>,
                 ::Tags::EventsAndTriggers<Triggers::WhenToCheck::AtSlabs>>;

  using batched_driver_component =
      ScalarWaveKokkosBatchedDriver<EvolutionMetavarsKokkosBatched>;

  using component_list =
      tmpl::list<observers::ObserverWriter<EvolutionMetavarsKokkosBatched>,
                 batched_driver_component>;

  static constexpr Options::String help{
      "ScalarWave Kokkos batched single-node scaffold executable.\n"
      "Initializes packed element metadata and evolves packed elements with "
      "batched Kokkos actions."};

  static constexpr auto default_phase_order = std::array<Parallel::Phase, 4>{
      Parallel::Phase::Initialization, Parallel::Phase::Register,
      Parallel::Phase::Evolve, Parallel::Phase::Exit};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
