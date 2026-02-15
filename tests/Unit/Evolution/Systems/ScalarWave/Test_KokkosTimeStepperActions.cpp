// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/CleanHistoryKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosTimeStepperState.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/RecordTimeStepperDataKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/UpdateUKokkos.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/OdeIntegration/OdeIntegration.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/ClassicalRungeKutta4.hpp"
#include "Time/TimeSteppers/Rk5Tsitouras.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename Vars>
double get_component(const Vars& vars, const size_t component, const size_t i) {
  if (component == 0) {
    return get(get<ScalarWave::Tags::Psi>(vars))[i];
  }
  if (component == 1) {
    return get(get<ScalarWave::Tags::Pi>(vars))[i];
  }
  return get<ScalarWave::Tags::Phi<1>>(vars).get(0)[i];
}

struct SyncEventMetavars {
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<tmpl::pair<Event, tmpl::list<>>>;
  };
  using component_list = tmpl::list<>;
};
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarWave.KokkosTimeStepperActions.Base",
    "[Unit][Evolution]") {
  using system = ScalarWave::System<1>;
  using initialize_action =
      ScalarWave::Actions::InitializeKokkosTimeStepperState<system>;
  using record_action =
      ScalarWave::Actions::RecordTimeStepperDataKokkos<system>;
  using update_action = ScalarWave::Actions::UpdateUKokkos<system>;
  using clean_action = ScalarWave::Actions::CleanHistoryKokkos<system>;

  using device_vars_tag = ScalarWave::KokkosTags::DeviceVariables<system>;
  using device_dt_tag = ScalarWave::KokkosTags::DeviceDtVariables<system>;
  using device_step_start_tag = ScalarWave::KokkosTags::DeviceStepStart<system>;
  using device_history_tag =
      ScalarWave::KokkosTags::DeviceDerivativeHistory<system>;
  using host_dt_vars_tag =
      db::add_tag_prefix<::Tags::dt, typename system::variables_tag>;

  constexpr size_t num_points = 3;
  constexpr size_t num_components =
      device_vars_tag::type::number_of_independent_components;

  typename system::variables_tag::type host_vars{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_vars)) = DataVector{{1.0, 2.0, 3.0}};
  get(get<ScalarWave::Tags::Pi>(host_vars)) = DataVector{{10.0, 20.0, 30.0}};
  get<ScalarWave::Tags::Phi<1>>(host_vars).get(0) =
      DataVector{{100.0, 200.0, 300.0}};

  typename host_dt_vars_tag::type host_dt_vars{num_points, 0.0};

  device_vars_tag::type device_vars{};
  device_dt_tag::type device_dt{};
  device_step_start_tag::type device_step_start{};
  device_history_tag::type device_derivative_history{};

  initialize_action::apply(
      make_not_null(&device_vars), make_not_null(&device_dt),
      make_not_null(&device_step_start),
      make_not_null(&device_derivative_history), host_vars, host_dt_vars);

  CHECK(device_derivative_history.extent(0) ==
        TimeSteppers::history_max_substeps);
  CHECK(device_derivative_history.extent(1) == num_points);
  CHECK(device_derivative_history.extent(2) == num_components);

  auto history_host = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, device_derivative_history);
  CHECK(history_host(0, 0, 0) == 0.0);
  CHECK(history_host(0, 1, 1) == 0.0);
  CHECK(history_host(0, 2, 2) == 0.0);

  const Slab slab(0.0, 1.0);
  const TimeDelta step = slab.duration();
  const Time step_time = slab.start();
  TimeSteppers::ClassicalRungeKutta4 rk4{};

  for (size_t substep = 0; substep < 4; ++substep) {
    auto dt_view = device_dt.view();
    Kokkos::parallel_for(
        "FillKokkosDt", num_points, KOKKOS_LAMBDA(const int i) {
          for (size_t c = 0; c < num_components; ++c) {
            dt_view(i, c) = static_cast<double>(100 * (substep + 1) + 10 * c +
                                                static_cast<size_t>(i));
          }
        });

    const double substep_fraction =
        substep == 0 ? 0.0 : rk4.butcher_tableau().substep_times[substep - 1];
    const TimeStepId substep_id(
        true, 0, step_time, substep, step,
        step_time.value() + substep_fraction * step.value());
    record_action::apply(make_not_null(&device_derivative_history), substep_id,
                         device_dt);
  }

  history_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                                     device_derivative_history);
  CHECK(history_host(0, 1, 2) == 121.0);
  CHECK(history_host(1, 1, 2) == 221.0);
  CHECK(history_host(2, 1, 2) == 321.0);
  CHECK(history_host(3, 1, 2) == 421.0);

  const TimeStepId final_substep_id(true, 0, step_time, 3, step,
                                    step_time.value());
  update_action::apply(make_not_null(&device_vars), rk4, final_substep_id, step,
                       device_step_start, device_derivative_history);

  typename system::variables_tag::type updated_host_vars{num_points, 0.0};
  copy_to_host(make_not_null(&updated_host_vars), device_vars);

  // RK4 result coefficients
  constexpr std::array<double, 4> b{
      {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0}};
  for (size_t i = 0; i < num_points; ++i) {
    for (size_t c = 0; c < num_components; ++c) {
      const double u0 = get_component(host_vars, c, i);
      const double k0 = static_cast<double>(100 * 1 + 10 * c + i);
      const double k1 = static_cast<double>(100 * 2 + 10 * c + i);
      const double k2 = static_cast<double>(100 * 3 + 10 * c + i);
      const double k3 = static_cast<double>(100 * 4 + 10 * c + i);
      const double expected =
          u0 + step.value() * (b[0] * k0 + b[1] * k1 + b[2] * k2 + b[3] * k3);
      CHECK(get_component(updated_host_vars, c, i) == approx(expected));
    }
  }

  clean_action::apply(make_not_null(&device_step_start),
                      make_not_null(&device_vars), rk4, final_substep_id);

  typename system::variables_tag::type cleaned_step_start{num_points, 0.0};
  typename system::variables_tag::type cleaned_device_vars{num_points, 0.0};
  copy_to_host(make_not_null(&cleaned_step_start), device_step_start);
  copy_to_host(make_not_null(&cleaned_device_vars), device_vars);

  for (size_t i = 0; i < num_points; ++i) {
    for (size_t c = 0; c < num_components; ++c) {
      CHECK(get_component(cleaned_step_start, c, i) ==
            approx(get_component(updated_host_vars, c, i)));
      CHECK(get_component(cleaned_device_vars, c, i) ==
            approx(get_component(host_vars, c, i)));
    }
  }
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarWave.KokkosTimeStepperActions.LinearOde",
    "[Unit][Evolution]") {
  using system = ScalarWave::System<1>;
  using initialize_action =
      ScalarWave::Actions::InitializeKokkosTimeStepperState<system>;
  using record_action =
      ScalarWave::Actions::RecordTimeStepperDataKokkos<system>;
  using update_action = ScalarWave::Actions::UpdateUKokkos<system>;
  using clean_action = ScalarWave::Actions::CleanHistoryKokkos<system>;

  using device_vars_tag = ScalarWave::KokkosTags::DeviceVariables<system>;
  using device_dt_tag = ScalarWave::KokkosTags::DeviceDtVariables<system>;
  using device_step_start_tag = ScalarWave::KokkosTags::DeviceStepStart<system>;
  using device_history_tag =
      ScalarWave::KokkosTags::DeviceDerivativeHistory<system>;
  using host_dt_vars_tag =
      db::add_tag_prefix<::Tags::dt, typename system::variables_tag>;

  constexpr size_t num_points = 2;
  constexpr size_t num_components =
      device_vars_tag::type::number_of_independent_components;
  constexpr size_t num_dofs = num_points * num_components;
  constexpr size_t num_steps = 50;

  auto rhs_coefficient = [](const size_t c, const size_t i) {
    return -0.2 - 0.05 * static_cast<double>(c) - 0.03 * static_cast<double>(i);
  };
  auto index = [](const size_t c, const size_t i) {
    return i * num_components + c;
  };

  typename system::variables_tag::type host_vars{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_vars)) = DataVector{{1.0, -0.5}};
  get(get<ScalarWave::Tags::Pi>(host_vars)) = DataVector{{0.3, 1.2}};
  get<ScalarWave::Tags::Phi<1>>(host_vars).get(0) = DataVector{{-0.7, 2.0}};
  typename host_dt_vars_tag::type host_dt_vars{num_points, 0.0};

  device_vars_tag::type device_vars{};
  device_dt_tag::type device_dt{};
  device_step_start_tag::type device_step_start{};
  device_history_tag::type device_derivative_history{};
  initialize_action::apply(
      make_not_null(&device_vars), make_not_null(&device_dt),
      make_not_null(&device_step_start),
      make_not_null(&device_derivative_history), host_vars, host_dt_vars);

  using state_type = std::array<double, num_dofs>;
  state_type reference_state{};
  for (size_t i = 0; i < num_points; ++i) {
    for (size_t c = 0; c < num_components; ++c) {
      reference_state[index(c, i)] = get_component(host_vars, c, i);
    }
  }

  const auto reference_rhs = [&rhs_coefficient, &index](const state_type& y,
                                                        state_type& dydt,
                                                        const double /*t*/) {
    for (size_t i = 0; i < num_points; ++i) {
      for (size_t c = 0; c < num_components; ++c) {
        dydt[index(c, i)] = rhs_coefficient(c, i) * y[index(c, i)];
      }
    }
  };

  const Slab slab(0.0, 1.0);
  const TimeDelta step = slab.duration() / 50;
  const double dt = step.value();
  const double reference_dt = dt / 200.0;
  Time step_start = slab.start();
  TimeSteppers::Rk5Tsitouras rk4{};
  boost::numeric::odeint::runge_kutta_dopri5<state_type> reference_stepper{};

  for (size_t step_index = 0; step_index < num_steps; ++step_index) {
    for (size_t substep = 0; substep < rk4.number_of_substeps(); ++substep) {
      const auto source_view =
          substep == 0 ? device_step_start.view() : device_vars.view();
      auto dt_view = device_dt.view();
      Kokkos::parallel_for(
          "ComputeLinearRhsOnDevice", num_points, KOKKOS_LAMBDA(const int i) {
            for (size_t c = 0; c < num_components; ++c) {
              const double coeff = -0.2 - 0.05 * static_cast<double>(c) -
                                   0.03 * static_cast<double>(i);
              dt_view(i, c) = coeff * source_view(i, c);
            }
          });

      const double substep_fraction =
          substep == 0 ? 0.0 : rk4.butcher_tableau().substep_times[substep - 1];
      const TimeStepId substep_id(true, 0, step_start, substep, step,
                                  step_start.value() + substep_fraction * dt);
      record_action::apply(make_not_null(&device_derivative_history),
                           substep_id, device_dt);
      update_action::apply(make_not_null(&device_vars), rk4, substep_id, step,
                           device_step_start, device_derivative_history);
    }

    const TimeStepId final_substep_id(true, 0, step_start,
                                      rk4.number_of_substeps() - 1, step,
                                      step_start.value() + dt);
    clean_action::apply(make_not_null(&device_step_start),
                        make_not_null(&device_vars), rk4, final_substep_id);

    boost::numeric::odeint::integrate_const(
        reference_stepper, reference_rhs, reference_state, step_start.value(),
        step_start.value() + dt, reference_dt);

    typename system::variables_tag::type device_step_start_host{num_points,
                                                                0.0};
    copy_to_host(make_not_null(&device_step_start_host), device_step_start);

    for (size_t i = 0; i < num_points; ++i) {
      for (size_t c = 0; c < num_components; ++c) {
        CHECK(get_component(device_step_start_host, c, i) ==
              approx(reference_state[index(c, i)]).epsilon(1.0e-12));
      }
    }

    step_start = step_start + step;
  }
}

SPECTRE_TEST_CASE("Unit.Evolution.Systems.ScalarWave.SyncKokkosToHostEvent",
                  "[Unit][Evolution]") {
  using system = ScalarWave::System<1>;
  using event = ScalarWave::Events::SyncKokkosToHost<system>;

  constexpr size_t num_points = 4;
  typename system::variables_tag::type host_vars{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_vars)) = DataVector{{1.0, 2.0, 3.0, 4.0}};
  get(get<ScalarWave::Tags::Pi>(host_vars)) = DataVector{{5.0, 6.0, 7.0, 8.0}};
  get<ScalarWave::Tags::Phi<1>>(host_vars).get(0) =
      DataVector{{9.0, 10.0, 11.0, 12.0}};

  auto device_vars = copy_to_device(host_vars);

  typename system::variables_tag::type host_copy_target{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_copy_target)) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};
  get(get<ScalarWave::Tags::Pi>(host_copy_target)) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};
  get<ScalarWave::Tags::Phi<1>>(host_copy_target).get(0) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};

  event sync_event{};
  Parallel::GlobalCache<SyncEventMetavars> cache{
      typename Parallel::GlobalCache<SyncEventMetavars>::ConstTagsTuple{}};
  const Event::ObservationValue observation_value{"Time", 0.0};
  sync_event(make_not_null(&host_copy_target), device_vars, cache, size_t{0},
             static_cast<const void*>(nullptr), observation_value);

  CHECK_ITERABLE_APPROX(get(get<ScalarWave::Tags::Psi>(host_copy_target)),
                        get(get<ScalarWave::Tags::Psi>(host_vars)));
  CHECK_ITERABLE_APPROX(get(get<ScalarWave::Tags::Pi>(host_copy_target)),
                        get(get<ScalarWave::Tags::Pi>(host_vars)));
  CHECK_ITERABLE_APPROX(get<ScalarWave::Tags::Phi<1>>(host_copy_target),
                        get<ScalarWave::Tags::Phi<1>>(host_vars));
}
