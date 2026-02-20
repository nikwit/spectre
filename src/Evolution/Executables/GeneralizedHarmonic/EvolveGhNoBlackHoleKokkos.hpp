// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/DiscontinuousGalerkin/InboxTags.hpp"
#include "Evolution/Executables/GeneralizedHarmonic/GeneralizedHarmonicBase.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/CleanHistoryKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosTimeStepperState.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/RecordTimeStepperDataKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/UpdateUKokkos.hpp"
#include "Options/String.hpp"
#include "Parallel/ArrayCollection/DgElementCollection.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/PhaseControl/PhaseControlTags.hpp"
#include "Parallel/Protocols/RegistrationMetavariables.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/GaugeWave.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/ChangeSlabSize/Action.hpp"
#include "Time/ChangeSlabSize/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ProtocolHelpers.hpp"

struct EvolutionMetavarsKokkos {
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using TimeStepperBase = TimeStepper;
  static constexpr bool local_time_stepping =
      TimeStepperBase::local_time_stepping;
  static constexpr bool use_dg_element_collection = false;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}

  using factory_creation =
      detail::FactoryCreation<volume_dim, local_time_stepping>;

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::push_back<
          tmpl::at<typename factory_creation::factory_classes, Event>>>;

  using initialize_initial_data_dependent_quantities_actions =
      tmpl::list<gh::gauges::SetPiAndPhiFromConstraints<
                     gh::Solutions::all_solutions<volume_dim>, volume_dim>,
                 Parallel::Actions::TerminatePhase>;

  using const_global_cache_tags =
      tmpl::list<gh::gauges::Tags::GaugeCondition,
                 evolution::initial_data::Tags::InitialData,
                 gh::Tags::DampingFunctionGamma0<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma1<volume_dim, Frame::Grid>,
                 gh::Tags::DampingFunctionGamma2<volume_dim, Frame::Grid>>;

  using dg_registration_list =
      tmpl::list<observers::Actions::RegisterEventsWithObservers>;

  static constexpr std::array<Parallel::Phase, 8> default_phase_order{
      Parallel::Phase::Initialization,
      Parallel::Phase::RegisterWithElementDataReader,
      Parallel::Phase::ImportInitialData,
      Parallel::Phase::InitializeInitialDataDependentQuantities,
      Parallel::Phase::Register,
      Parallel::Phase::InitializeTimeStepperHistory,
      Parallel::Phase::Evolve,
      Parallel::Phase::Exit};

  using step_actions = tmpl::list<
      gh::Actions::ComputeTimeDerivativeKokkos,
      tmpl::list<
          gh::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos,
          Actions::MutateApply<gh::Actions::RecordTimeStepperDataKokkos<system>>,
          evolution::Actions::RunEventsAndDenseTriggers<tmpl::list<>>,
          control_system::Actions::LimitTimeStep<tmpl::list<>>,
          Actions::MutateApply<gh::Actions::UpdateUKokkos<system>>>,
      Actions::MutateApply<gh::Actions::CleanHistoryKokkos<system>>,
      dg::Actions::Filter<
          Filters::Exponential<0>,
          tmpl::list<gr::Tags::SpacetimeMetric<DataVector, volume_dim>,
                     gh::Tags::Pi<DataVector, volume_dim>,
                     gh::Tags::Phi<DataVector, volume_dim>>>>;

  using initialization_actions = tmpl::list<
      Initialization::Actions::InitializeItems<
          Initialization::TimeStepping<EvolutionMetavarsKokkos,
                                       TimeStepperBase>,
          evolution::dg::Initialization::Domain<EvolutionMetavarsKokkos,
                                                false>,
          Initialization::TimeStepperHistory<EvolutionMetavarsKokkos>>,
      Initialization::Actions::NonconservativeSystem<system>,
      Initialization::Actions::AddComputeTags<::Tags::DerivCompute<
          typename system::variables_tag, domain::Tags::Mesh<volume_dim>,
          domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                        Frame::Inertial>,
          typename system::gradient_variables>>,
      gh::Actions::InitializeGhAnd3Plus1Variables<volume_dim>,
      Initialization::Actions::AddComputeTags<
          tmpl::push_back<StepChoosers::step_chooser_compute_tags<
              EvolutionMetavarsKokkos, local_time_stepping>>>,
      ::evolution::dg::Initialization::Mortars<volume_dim, system>,
      evolution::Actions::InitializeRunEventsAndDenseTriggers,
      Actions::MutateApply<gh::Actions::InitializeKokkosTimeStepperState<
          system>>,
      Actions::MutateApply<gh::Actions::InitializeKokkosTags<system>>,
      Actions::MutateApply<
          gh::Actions::InitializeKokkosBoundaryCommunication<volume_dim>>,
      Parallel::Actions::TerminatePhase>;

  using gh_dg_element_array = DgElementArray<
      EvolutionMetavarsKokkos,
      tmpl::flatten<tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::RegisterWithElementDataReader,
              tmpl::list<importers::Actions::RegisterWithElementDataReader,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::ImportInitialData,
              tmpl::list<gh::Actions::SetInitialData,
                         gh::Actions::ReceiveNumericInitialData,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions, system>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<Parallel::Phase::Restart,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<::evolution::Actions::RunEventsAndTriggers<
                             local_time_stepping>,
                         Actions::ChangeSlabSize, step_actions,
                         Actions::AdvanceTime,
                         PhaseControl::Actions::ExecutePhaseChange>>>>>;

  struct registration
      : tt::ConformsTo<Parallel::protocols::RegistrationMetavariables> {
    using element_registrars =
        tmpl::map<tmpl::pair<gh_dg_element_array, dg_registration_list>>;
  };

  using component_list =
      tmpl::flatten<tmpl::list<observers::Observer<EvolutionMetavarsKokkos>,
                               observers::ObserverWriter<
                                   EvolutionMetavarsKokkos>,
                               mem_monitor::MemoryMonitor<
                                   EvolutionMetavarsKokkos>,
                               importers::ElementDataReader<
                                   EvolutionMetavarsKokkos>,
                               gh_dg_element_array>>;

  static constexpr Options::String help{
      "Generalized Harmonic Kokkos development executable (single-node, "
      "global stepping scaffold)."};
};
