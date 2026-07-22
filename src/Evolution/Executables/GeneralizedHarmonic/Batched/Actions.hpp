// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/Actions.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/Tags.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Time/Kokkos/UpdateU.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeSteppers/RungeKutta.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Batched::Actions {

template <typename Metavariables,
          typename FunctionsOfTimeTag = domain::Tags::FunctionsOfTimeInitialize>
using functions_of_time_compute_tag =
    Parallel::Tags::FromGlobalCache<FunctionsOfTimeTag, Metavariables>;

template <typename Metavariables,
          typename FunctionsOfTimeTag = domain::Tags::FunctionsOfTimeInitialize>
struct InitializeDriver
    : gh::Actions::InitializeBatchedDataWithTag<FunctionsOfTimeTag, false> {
  using compute_tags = tmpl::list<
      functions_of_time_compute_tag<Metavariables, FunctionsOfTimeTag>>;
};

template <typename Metavariables,
          typename FunctionsOfTimeTag = domain::Tags::FunctionsOfTimeInitialize>
struct UpdateFrameFieldsBatched
    : gh::Actions::UpdateFrameFieldsBatchedWithTag<FunctionsOfTimeTag, false> {
  using compute_tags = tmpl::list<
      functions_of_time_compute_tag<Metavariables, FunctionsOfTimeTag>>;
};

template <typename Metavariables,
          typename FunctionsOfTimeTag = domain::Tags::FunctionsOfTimeInitialize>
struct UpdateConstraintDampingBatched
    : gh::Actions::UpdateConstraintDampingBatchedWithTag<FunctionsOfTimeTag,
                                                         false> {
  using compute_tags = tmpl::list<
      functions_of_time_compute_tag<Metavariables, FunctionsOfTimeTag>>;
};

using InitializeBoundaryBatchMetadata =
    gh::Actions::InitializeBoundaryBatchMetadata;
using ComputeTimeDerivativeBatched = gh::Actions::ComputeTimeDerivativeBatched;
using PackageLocalFacesBatched = gh::Actions::PackageLocalFacesBatched;
using ComputeInternalBoundaryTermsBatched =
    gh::Actions::ComputeInternalBoundaryTermsBatched;
using LiftInternalBoundaryTermsBatched =
    gh::Actions::LiftInternalBoundaryTermsBatched;
using ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched =
    gh::Actions::ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched;
template <typename FilterType>
using FilterKokkosBatched = gh::Actions::FilterKokkosBatched<FilterType>;

/// \brief Run dense triggers and events for the batched GH Kokkos Singleton.
///
/// Mirrors the logic of `evolution::Actions::RunEventsAndDenseTriggers` but
/// without the `HistoryEvolvedVariables` dependency, since the batched driver
/// stores history in `PackedEvolutionState` rather than the standard tag.
/// Dense output is implemented directly on the packed Kokkos time-stepper
/// storage.
template <typename FunctionsOfTimeTag = domain::Tags::FunctionsOfTimeInitialize>
struct RunEventsAndDenseTriggersBatched {
 private:
  static double evaluate_polynomial(const std::vector<double>& coefficients,
                                    const double x) {
    double result = 0.0;
    for (auto coefficient = coefficients.rbegin();
         coefficient != coefficients.rend(); ++coefficient) {
      result = result * x + *coefficient;
    }
    return result;
  }

  template <typename PackedEvolutionState>
  static bool dense_update_u(
      const gsl::not_null<PackedEvolutionState*> packed_evolution_state,
      const TimeStepper& time_stepper, const TimeStepId& time_step_id,
      const TimeDelta& time_step, const double output_time) {
    const auto* const runge_kutta =
        dynamic_cast<const TimeSteppers::RungeKutta*>(&time_stepper);
    ASSERT(runge_kutta != nullptr,
           "Batched dense output currently supports only Runge-Kutta time "
           "steppers.");

    const double step_start = time_step_id.step_time().value();
    auto& device_vars = packed_evolution_state->device_variables;
    if (output_time == step_start) {
      ::Kokkos::deep_copy(device_vars.view(),
                          packed_evolution_state->device_step_start.view());
      return true;
    }

    const auto& tableau = runge_kutta->butcher_tableau();
    const size_t number_of_dense_coefficients =
        tableau.dense_coefficients.size();
    if (number_of_dense_coefficients == 0 or
        time_step_id.substep() + 1 < number_of_dense_coefficients) {
      return false;
    }

    constexpr size_t max_supported_coefficients = 8;
    ASSERT(number_of_dense_coefficients <= max_supported_coefficients,
           "Batched dense output supports up to "
               << max_supported_coefficients << " dense coefficients, not "
               << number_of_dense_coefficients);
    ASSERT(number_of_dense_coefficients <=
               static_cast<size_t>(
                   packed_evolution_state->device_derivative_history.extent(0)),
           "Not enough derivative history entries for dense output.");

    const double output_fraction =
        (output_time - step_start) / time_step.value();
    ASSERT(output_fraction >= 0.0 and output_fraction <= 1.0,
           "Dense output requested outside the current step at time "
               << output_time << " with step start " << step_start
               << " and step size " << time_step.value());

    ::Kokkos::Array<double, max_supported_coefficients> coefficients_array{};
    for (size_t i = 0; i < number_of_dense_coefficients; ++i) {
      coefficients_array[i] =
          evaluate_polynomial(tableau.dense_coefficients[i], output_fraction);
    }

    const auto u_view = device_vars.view();
    const auto u0_view = packed_evolution_state->device_step_start.view();
    const auto derivative_history =
        packed_evolution_state->device_derivative_history;
    std::array<size_t, 8> u_strides{};
    u_view.stride(u_strides.data());
    std::array<size_t, 8> u0_strides{};
    u0_view.stride(u0_strides.data());
    std::array<size_t, 8> derivative_history_strides{};
    derivative_history.stride(derivative_history_strides.data());

    ::Actions::Kokkos::detail::update_u_impl(
        u_view.data(), u_strides[0], u_strides[1], u0_view.data(),
        u0_strides[0], u0_strides[1], derivative_history.data(),
        derivative_history_strides[0], derivative_history_strides[1],
        derivative_history_strides[2], coefficients_array.data(),
        number_of_dense_coefficients, time_step.value(), u_view.extent(0),
        u_view.extent(1));
    return true;
  }

 public:
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const component) {
    const auto& time_step_id = db::get<::Tags::TimeStepId>(box);
    if (time_step_id.slab_number() < 0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    using system = typename Metavariables::system;
    using packed_evolution_state_tag =
        evolution::Kokkos::Tags::PackedEvolutionState<system>;
    using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;
    using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;

    auto& events_and_dense_triggers =
        db::get_mutable_reference<::Tags::EventsAndDenseTriggers>(
            make_not_null(&box));

    const auto step_end =
        time_step_id.step_time() + db::get<::Tags::TimeStep>(box);
    const evolution_less_equal<double> before_equal{
        time_step_id.time_runs_forward()};

    // Save the original time so we can restore it after processing triggers.
    const double original_time = db::get<::Tags::Time>(box);
    bool time_was_modified = false;
    auto& packed_evolution_state =
        db::get_mutable_reference<packed_evolution_state_tag>(
            make_not_null(&box));
    typename packed_evolution_state_tag::type::device_variables_type
        original_device_variables{};
    bool variables_were_modified = false;

    const auto restore_state = [&]() {
      if (variables_were_modified) {
        ::Kokkos::deep_copy(packed_evolution_state.device_variables.view(),
                            original_device_variables.view());
        variables_were_modified = false;
      }
      if (time_was_modified) {
        db::mutate<::Tags::Time>(
            [&original_time](const gsl::not_null<double*> time) {
              *time = original_time;
            },
            make_not_null(&box));
        time_was_modified = false;
      }
    };
    const auto save_variables_if_needed = [&]() {
      if (not variables_were_modified) {
        original_device_variables =
            typename packed_evolution_state_tag::type::device_variables_type{
                packed_evolution_state.device_variables.view().extent(0)};
        ::Kokkos::deep_copy(original_device_variables.view(),
                            packed_evolution_state.device_variables.view());
        variables_were_modified = true;
      }
    };
    const auto update_frame_fields = [&](const double time) {
      db::mutate<packed_geometry_tag>(
          [&cache, &time](
              const gsl::not_null<typename packed_geometry_tag::type*>
                  packed_geometry,
              const typename packed_topology_tag::type& packed_topology) {
            gh::Actions::UpdateFrameFieldsBatchedWithTag<FunctionsOfTimeTag,
                                                         false>::
                apply(packed_geometry, packed_topology,
                      Parallel::get<domain::Tags::Domain<system::volume_dim>>(
                          cache),
                      Parallel::get<FunctionsOfTimeTag>(cache), time);
          },
          make_not_null(&box), db::get<packed_topology_tag>(box));
    };

    for (;;) {
      const double next_trigger = events_and_dense_triggers.next_trigger(box);
      if (before_equal(step_end.value(), next_trigger)) {
        restore_state();
        return {Parallel::AlgorithmExecution::Continue, std::nullopt};
      }

      if (db::get<::Tags::Time>(box) != next_trigger) {
        db::mutate<::Tags::Time>(
            [&next_trigger](const gsl::not_null<double*> time) {
              *time = next_trigger;
            },
            make_not_null(&box));
        time_was_modified = true;
      }

      const auto triggered = events_and_dense_triggers.is_ready(
          make_not_null(&box), cache, array_index, component);
      using TriggeringState = EventsAndDenseTriggers::TriggeringState;
      switch (triggered) {
        case TriggeringState::NotReady:
          restore_state();
          return {Parallel::AlgorithmExecution::Retry, std::nullopt};
        case TriggeringState::NeedsEvolvedVariables:
          save_variables_if_needed();
          if (not dense_update_u(make_not_null(&packed_evolution_state),
                                 db::get<::Tags::TimeStepper<TimeStepper>>(box),
                                 time_step_id, db::get<::Tags::TimeStep>(box),
                                 next_trigger)) {
            restore_state();
            return {Parallel::AlgorithmExecution::Continue, std::nullopt};
          }
          update_frame_fields(next_trigger);
          break;
        default:
          break;
      }

      events_and_dense_triggers.run_events(box, cache, array_index, component);
      if (not events_and_dense_triggers.reschedule(make_not_null(&box), cache,
                                                   array_index, component)) {
        restore_state();
        return {Parallel::AlgorithmExecution::Retry, std::nullopt};
      }
    }
  }
};

}  // namespace gh::Batched::Actions
