// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "DataStructures/Tensor/IndexType.hpp"
#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/ComputeTags.hpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Initialization/NonconservativeSystem.hpp"
#include "Evolution/Initialization/SetVariables.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/ScalarWave/EnergyDensity.hpp"
#include "Evolution/Systems/ScalarWave/Equations.hpp"
#include "Evolution/Systems/ScalarWave/Initialize.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/CleanHistoryKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/CudaDiagnostics.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosTimeStepperState.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/RecordTimeStepperDataKokkos3D.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/UpdateUKokkos3D.hpp"
#include "Evolution/Systems/ScalarWave/MomentumDensity.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/PhaseControl/Factory.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/Events/Factory.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/Factory.hpp"
#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/PlaneWave.hpp"
#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/RegularSphericalWave.hpp"
#include "PointwiseFunctions/InitialDataUtilities/NumericData.hpp"
#include "PointwiseFunctions/MathFunctions/Factory.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/ChangeSlabSize/Action.hpp"
#include "Time/ChangeSlabSize/Tags.hpp"
#include "Time/StepChoosers/ByBlock.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/Functional.hpp"
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

struct EvolutionMetavarsKokkos {
  static constexpr size_t volume_dim = 3;

  using system = ScalarWave::System<volume_dim>;
  using initial_data_list = ScalarWave::Solutions::all_solutions<volume_dim>;
  using temporal_id = Tags::TimeStepId;
  using TimeStepperBase = TimeStepper;
  static constexpr bool local_time_stepping = false;
  static constexpr bool use_dg_element_collection = false;

  using analytic_solution_fields = typename system::variables_tag::tags_list;
  using deriv_compute = ::Tags::DerivCompute<
      typename system::variables_tag, domain::Tags::Mesh<volume_dim>,
      domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                    Frame::Inertial>,
      typename system::gradient_variables>;
  using analytic_compute = evolution::Tags::AnalyticSolutionsCompute<
      volume_dim, analytic_solution_fields, true, initial_data_list>;
  using error_compute = Tags::ErrorsCompute<analytic_solution_fields>;
  using error_tags = db::wrap_tags_in<Tags::Error, analytic_solution_fields>;

  using observe_fields = tmpl::push_back<
      tmpl::append<typename system::variables_tag::tags_list,
                   typename deriv_compute::type::tags_list, error_tags>,
      ScalarWave::Tags::EnergyDensityCompute<volume_dim>,
      domain::Tags::Coordinates<volume_dim, Frame::Grid>,
      domain::Tags::Coordinates<volume_dim, Frame::Inertial>>;
  using non_tensor_compute_tags = tmpl::list<
      ::Events::Tags::ObserverMeshCompute<volume_dim>,
      ::Events::Tags::ObserverCoordinatesCompute<volume_dim, Frame::Inertial>,
      ::Events::Tags::ObserverDetInvJacobianCompute<Frame::ElementLogical,
                                                    Frame::Inertial>,
      deriv_compute, analytic_compute, error_compute>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<Event,
                   tmpl::flatten<tmpl::list<
                       Events::Completion,
                       ScalarWave::Events::SyncKokkosToHost<system>,
                       dg::Events::field_observations<
                           volume_dim, observe_fields, non_tensor_compute_tags>,
                       Events::time_events<system>>>>,
        tmpl::pair<evolution::initial_data::InitialData,
                   tmpl::push_back<initial_data_list,
                                   evolution::initial_data::NumericData>>,
        tmpl::pair<MathFunction<1, Frame::Inertial>,
                   MathFunctions::all_math_functions<1, Frame::Inertial>>,
        tmpl::pair<PhaseChange, PhaseControl::factory_creatable_classes>,
        tmpl::pair<
            ScalarWave::BoundaryConditions::BoundaryCondition<volume_dim>,
            ScalarWave::BoundaryConditions::standard_boundary_conditions<
                volume_dim>>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                   tmpl::push_back<StepChoosers::standard_step_choosers<system>,
                                   StepChoosers::ByBlock<volume_dim>>>,
        tmpl::pair<StepChooser<StepChooserUse::Slab>,
                   StepChoosers::standard_slab_choosers<system, false>>,
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

  using step_actions = tmpl::flatten<tmpl::list<
      ScalarWave::Actions::ComputeTimeDerivativeKokkos,
      ScalarWave::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos,
      Actions::MutateApply<ScalarWave::Actions::RecordTimeStepperDataKokkos3D>,
      Actions::MutateApply<ScalarWave::Actions::UpdateUKokkos3D>,
      Actions::MutateApply<
          ScalarWave::Actions::CleanHistoryKokkos<ScalarWave::System<3>>>>>;

#ifdef SPECTRE_DEBUG
  using evolve_entry_cuda_diagnostics =
      tmpl::list<ScalarWave::Actions::CheckCudaOnEvolveEntry>;
#else
  using evolve_entry_cuda_diagnostics = tmpl::list<>;
#endif  // SPECTRE_DEBUG

  using const_global_cache_tags =
      tmpl::list<evolution::initial_data::Tags::InitialData>;

  using dg_registration_list =
      tmpl::list<observers::Actions::RegisterEventsWithObservers>;

  using initialization_actions = tmpl::list<
      Initialization::Actions::InitializeItems<
          Initialization::TimeStepping<EvolutionMetavarsKokkos,
                                       TimeStepperBase>,
          evolution::dg::Initialization::Domain<EvolutionMetavarsKokkos>,
          Initialization::TimeStepperHistory<EvolutionMetavarsKokkos>>,
      Initialization::Actions::NonconservativeSystem<system>,
      evolution::Initialization::Actions::SetVariables<
          domain::Tags::Coordinates<volume_dim, Frame::ElementLogical>>,
      ScalarWave::Actions::InitializeConstraints<volume_dim>,
      Initialization::Actions::AddComputeTags<
          StepChoosers::step_chooser_compute_tags<EvolutionMetavarsKokkos,
                                                  local_time_stepping>>,
      Actions::MutateApply<
          ScalarWave::Actions::InitializeKokkosTimeStepperState<system>>,
      Actions::MutateApply<ScalarWave::Actions::InitializeKokkosTags<system>>,
      Actions::MutateApply<
          ScalarWave::Actions::InitializeKokkosBoundaryCommunication<
              volume_dim>>,
      ScalarWave::Actions::CheckCudaAfterInitialization,
      Parallel::Actions::TerminatePhase>;

  using dg_element_array = DgElementArray<
      EvolutionMetavarsKokkos,
      tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              tmpl::append<
                  tmpl::list<ScalarWave::Actions::
                                 CheckCudaOnInitializeTimeStepperHistoryEntry>,
                  SelfStart::self_start_procedure<step_actions, system>,
                  tmpl::list<ScalarWave::Actions::
                                 CheckCudaOnInitializeTimeStepperHistoryExit>>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<Parallel::Phase::Restart,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<
                  // ScalarWave::Actions::CheckCudaOnEvolveEntry,
                  evolution::Actions::RunEventsAndTriggers<local_time_stepping>,
                  Actions::ChangeSlabSize, step_actions, Actions::AdvanceTime,
                  PhaseControl::Actions::ExecutePhaseChange>>>>;

  using component_list =
      tmpl::list<observers::Observer<EvolutionMetavarsKokkos>,
                 observers::ObserverWriter<EvolutionMetavarsKokkos>,
                 dg_element_array>;

  static constexpr Options::String help{
      "Minimal 3D ScalarWave Kokkos executable.\n"
      "Assumes global time stepping and periodic/internal-face communication."};

  static constexpr auto default_phase_order = std::array<Parallel::Phase, 5>{
      Parallel::Phase::Initialization,
      Parallel::Phase::InitializeTimeStepperHistory, Parallel::Phase::Register,
      Parallel::Phase::Evolve, Parallel::Phase::Exit};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
