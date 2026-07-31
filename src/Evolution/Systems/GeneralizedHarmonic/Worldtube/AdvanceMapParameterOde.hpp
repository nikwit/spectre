// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Time/History.hpp"
#include "Time/SelfStart.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Worldtube::Actions {
/*!
 * \brief Advance the map-parameter ODE through the element's own time
 * stepper (the `StepperOde` matcher mode).
 *
 * Runs directly after `ComputeTimeDerivative`, so the acceleration is
 * measured from *this* substep's right-hand sides (zero staleness) and the
 * 26-component state (p, pdot) is recorded and updated on exactly the same
 * (sub)step schedule as the evolved fields, by the same `TimeStepper`.
 * The integration order is slaved to the system history's current order;
 * the ODE is dormant during self-start (exact for Schwarzschild initial
 * data); a repeated or regressed `TimeStepId` (step rejection under local
 * time stepping) rewinds the history before recording.
 */
struct AdvanceMapParameterOde {
  using const_global_cache_tags = tmpl::list<Tags::Matcher>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    static constexpr size_t Dim = 3;
    const auto& config_opt = db::get<Tags::Matcher>(box);
    if (not config_opt.has_value() or not config_opt->stepper_ode) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    if (config_opt->rate_ode or config_opt->second_order_ode or
        config_opt->fit_center_offset) {
      ERROR(
          "WorldtubeMatcher: StepperOde cannot be combined with RateOde, "
          "SecondOrderOde or FitCenterOffset.");
    }
    const auto& element = db::get<domain::Tags::Element<Dim>>(box);
    if (element.id().block_id() != 0 or
        element.external_boundaries().count(Direction<Dim>::lower_xi()) ==
            0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    const TimeStepId& time_step_id = db::get<::Tags::TimeStepId>(box);
    if (SelfStart::is_self_starting(time_step_id)) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    const size_t n_radial = mesh.extents(0);
    const size_t n_theta = mesh.extents(1);
    const size_t n_phi = mesh.extents(2);
    const size_t l_max = n_theta - 1;
    const size_t n_face = n_theta * n_phi;
    const size_t radial_index = config_opt->fit_radial_index;

    const auto face_slice = [n_radial, n_face,
                             radial_index](const DataVector& volume) {
      DataVector face(n_face);
      for (size_t k = 0; k < n_face; ++k) {
        face[k] = volume[k * n_radial + radial_index];
      }
      return face;
    };
    const auto slice_tensor = [&face_slice](auto face_tensor,
                                            const auto& volume_tensor) {
      for (size_t storage = 0; storage < volume_tensor.size(); ++storage) {
        (*face_tensor)[storage] = face_slice(volume_tensor[storage]);
      }
    };

    tnsr::aa<DataVector, Dim> metric_face{};
    tnsr::aa<DataVector, Dim> pi_face{};
    tnsr::iaa<DataVector, Dim> phi_face{};
    tnsr::I<DataVector, Dim> coords_face{};
    tnsr::aa<DataVector, Dim> dt_metric_face{};
    tnsr::aa<DataVector, Dim> dt_pi_face{};
    tnsr::iaa<DataVector, Dim> dt_phi_face{};
    slice_tensor(make_not_null(&metric_face),
                 db::get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&pi_face),
                 db::get<gh::Tags::Pi<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&phi_face),
                 db::get<gh::Tags::Phi<DataVector, Dim>>(box));
    slice_tensor(
        make_not_null(&coords_face),
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box));
    slice_tensor(
        make_not_null(&dt_metric_face),
        db::get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(box));
    slice_tensor(make_not_null(&dt_pi_face),
                 db::get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(box));
    slice_tensor(make_not_null(&dt_phi_face),
                 db::get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(box));

    const ylm::Spherepack& ylm_transform = ylm::get_spherepack_cache(l_max);
    // current pdot state for the exact (pdot-quadratic) Hessian term of the
    // model's second time derivative
    std::array<double, num_map_parameters> pdot_state{};
    {
      const auto& current = db::get<Tags::MapParameters>(box);
      if (current.ode_state.size() == 2 * num_map_parameters) {
        for (size_t a = 0; a < num_map_parameters; ++a) {
          gsl::at(pdot_state, a) =
              current.ode_state[num_map_parameters + a];
        }
      }
    }
    const RateFitResult accel = fit_map_parameter_accelerations(
        metric_face, pi_face, phi_face, dt_metric_face, dt_pi_face,
        dt_phi_face, pdot_state, coords_face, ylm_transform, *config_opt);
    bool finite = true;
    for (size_t a = 0; a < num_map_parameters; ++a) {
      finite = finite and std::isfinite(gsl::at(accel.pdot, a));
    }
    if (not finite) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& stepper = db::get<::Tags::TimeStepper<TimeStepper>>(box);
    const auto& time_step = db::get<::Tags::TimeStep>(box);
    const size_t system_order =
        db::get<::Tags::HistoryEvolvedVariables<
            typename gh::System<Dim>::variables_tag>>(box)
            .integration_order();

    db::mutate<Tags::MapParameters>(
        [&](const gsl::not_null<MapParameterData*> data) {
          auto& history = data->ode_history;
          auto& y = data->ode_state;
          if (y.size() != 2 * num_map_parameters) {
            y = DataVector(2 * num_map_parameters, 0.);
          }
          // transactional rewind for repeated/regressed TimeStepIds (step
          // rejection under local time stepping)
          while (history.size() > 0 or not history.substeps().empty()) {
            const TimeStepId& latest = history.at_step_start()
                                           ? history.back().time_step_id
                                           : history.substeps().back()
                                                 .time_step_id;
            if (latest < time_step_id) {
              break;
            }
            history.undo_latest();
          }
          // one-step methods (RK) require the history order to equal their
          // fixed order and keep no cross-step records; multistep (Adams)
          // ramps with the ODE's own history depth
          const bool multistep = stepper.number_of_past_steps() > 0;
          if (time_step_id.substep() == 0) {
            if (data->ode_step_id == time_step_id and
                data->ode_step_start.size() == y.size()) {
              // retry of a rejected step: y holds the rejected
              // end-of-step value, not the value at the step start
              y = data->ode_step_start;
            } else {
              data->ode_step_start = y;
              data->ode_step_id = time_step_id;
            }
            // Multistep: the order the current history depth supports (the
            // update requires size >= order - 1 after this step's record is
            // inserted), capped by the system's order. Ramps 2, 3, ...
            // over the first steps, exactly as self-start would.
            history.integration_order(
                multistep
                    ? std::clamp(history.size() + 2, size_t{2}, system_order)
                    : system_order);
          }
          const size_t order = history.integration_order();
          DataVector dt_y(2 * num_map_parameters);
          for (size_t a = 0; a < num_map_parameters; ++a) {
            dt_y[a] = y[num_map_parameters + a];
            dt_y[num_map_parameters + a] = gsl::at(accel.pdot, a);
          }
          if (const double gamma = config_opt->gauge_damping; gamma > 0.) {
            // damped-oscillator gauge fixing of the free parameters,
            // unfolded through the pins so the pinned relations stay exact
            static constexpr std::array<size_t, 9> free_indices{
                {0, 1, 2, 3, 7, 8, 9, 10, 11}};
            std::array<double, num_map_parameters> damp{};
            for (const size_t a : free_indices) {
              gsl::at(damp, a) = -2. * gamma * y[num_map_parameters + a] -
                                 gamma * gamma * y[a];
            }
            for (size_t i = 0; i < 3; ++i) {
              gsl::at(damp, 4 + i) =
                  damp[0] * gsl::at(config_opt->center_velocity, i);
            }
            damp[12] = -damp[7] - damp[10];
            for (size_t a = 0; a < num_map_parameters; ++a) {
              dt_y[num_map_parameters + a] += gsl::at(damp, a);
            }
          }
          history.insert(time_step_id, y, dt_y);
          stepper.update_u(make_not_null(&y), history, time_step);
          if (multistep and
              time_step_id.substep() + 1 == stepper.number_of_substeps()) {
            // pre-arm order growth: cleaning keeps order - 2 records, so
            // clean at next step's intended order
            history.integration_order(std::min(order + 1, system_order));
          }
          // mirror the harness: clean after every update; the stepper's
          // implementation no-ops when not applicable
          stepper.clean_history(make_not_null(&history));
          for (size_t a = 0; a < num_map_parameters; ++a) {
            gsl::at(data->p, a) = y[a];
            gsl::at(data->pdot, a) = y[num_map_parameters + a];
            gsl::at(data->pddot, a) = gsl::at(accel.pdot, a);
          }
          data->last_fit_time = time_step_id.substep_time();
          data->valid = true;
        },
        make_not_null(&box));

    // occasional diagnostic row (same subfile/format as FitMapParameters)
    const auto& state = db::get<Tags::MapParameters>(box);
    const double time = time_step_id.substep_time();
    if (time_step_id.substep() == 0 and
        time >= state.previous_fit_time + 0.5 - 1.0e-12) {
      db::mutate<Tags::MapParameters>(
          [&time](const gsl::not_null<MapParameterData*> data) {
            data->previous_fit_time = time;
          },
          make_not_null(&box));
      auto& writer = Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache);
      std::vector<std::string> legend{"Time"};
      static const std::array<std::string, num_map_parameters> names{
          {"qdot0", "b_x", "b_y", "b_z", "v_x", "v_y", "v_z", "s_xx", "s_xy",
           "s_xz", "s_yy", "s_yz", "s_zz"}};
      for (const auto& name : names) {
        legend.push_back(name);
      }
      for (const auto& name : names) {
        legend.push_back("dt_" + name);
      }
      legend.emplace_back("q_x");
      legend.emplace_back("q_y");
      legend.emplace_back("q_z");
      legend.emplace_back("ResidualInitial");
      legend.emplace_back("ResidualFinal");
      legend.emplace_back("Iterations");
      std::vector<double> row;
      row.reserve(3 * num_map_parameters);
      row.push_back(time);
      for (size_t a = 0; a < num_map_parameters; ++a) {
        row.push_back(gsl::at(state.p, a));
      }
      for (size_t a = 0; a < num_map_parameters; ++a) {
        row.push_back(gsl::at(state.pdot, a));
      }
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(0.);
      }
      row.push_back(accel.residual_initial);
      row.push_back(accel.residual_final);
      row.push_back(1.);
      Parallel::threaded_action<
          observers::ThreadedActions::WriteReductionDataRow>(
          writer[0], std::string{"/WorldtubeMatcher"}, std::move(legend),
          std::make_tuple(std::move(row)));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace gh::Worldtube::Actions
