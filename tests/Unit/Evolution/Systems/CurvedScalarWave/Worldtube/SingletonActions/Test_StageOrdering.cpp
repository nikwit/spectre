// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonChare.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/ClassicalRungeKutta4.hpp"
#include "Time/TimeSteppers/Rk5Tsitouras.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {
using system = System<3>;
using variables_tag = system::variables_tag;
using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
using history_tag = ::Tags::HistoryEvolvedVariables<variables_tag>;
using psi_monopole_tag =
    Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, 3, Frame::Inertial>;

// Only the inbox is needed: check the received polynomial and its derivative
// directly against the analytic field, independently of ReceiveWorldtubeData.
struct ReceiveData {
  using inbox_tags = tmpl::list<Tags::RegularFieldInbox<3>>;
  template <typename... Args>
  static Parallel::iterable_action_return_t apply(Args&&... /*args*/) {
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

template <typename Metavariables>
struct MockElement {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<3>;
  using array_allocation_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 db::AddSimpleTags<>, db::AddComputeTags<>>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             tmpl::list<ReceiveData>>>;
};

template <typename Metavariables>
struct MockWorldtube {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using simple_tags_from_options = tmpl::list<>;
  using component_being_mocked = WorldtubeSingleton<Metavariables>;
  using simple_tags = tmpl::list<
      variables_tag, dt_variables_tag, history_tag, ::Tags::TimeStepId,
      ::Tags::Next<::Tags::TimeStepId>, ::Tags::TimeStep, ::Tags::Time,
      Tags::ElementFacesGridCoordinates<3>, Tags::WorldtubeRadius,
      Tags::GeodesicAcceleration<3>, Tags::BackgroundQuantities<3>,
      psi_monopole_tag,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, 3,
                           Frame::Inertial>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, 3, Frame::Inertial>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, 3,
                           Frame::Inertial>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, 3, Frame::Inertial>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 2, 3,
                           Frame::Inertial>>;
  using compute_tags =
      tmpl::append<time_stepper_ref_tags<TimeStepper>,
                   tmpl::list<Tags::EvolvedParticlePositionVelocityCompute<3>>>;
  using tested_actions =
      tmpl::list<::Actions::MutateApply<UpdateAcceleration>,
                 ::Actions::MutateApply<RecordTimeStepperData<system>>,
                 ::Actions::MutateApply<UpdateU<system, false>>,
                 ::Actions::MutateApply<CleanHistory<system>>,
                 Actions::SendToElements<Metavariables>>;
  // Preserve the actual production order, so this test catches a send moved
  // past UpdateU. Boundary projections are supplied analytically below.
  using step_actions = tmpl::filter<
      typename component_being_mocked::step_actions,
      tmpl::bind<tmpl::list_contains, tmpl::pin<tested_actions>, tmpl::_1>>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 simple_tags, compute_tags>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing, step_actions>>;
};

struct Metavariables {
  static constexpr size_t volume_dim = 3;
  struct factory_creation {
    using factory_classes =
        tmpl::map<tmpl::pair<StepChooser<StepChooserUse::Slab>, tmpl::list<>>>;
  };
  using dg_element_array = MockElement<Metavariables>;
  using component_list =
      tmpl::list<MockWorldtube<Metavariables>, dg_element_array>;
  using const_global_cache_tags =
      tmpl::list<::Tags::ConcreteTimeStepper<TimeStepper>,
                 ::Tags::EventsAndTriggers<Triggers::WhenToCheck::AtSlabs>,
                 Tags::ExpansionOrder, Tags::Charge, Tags::Mass,
                 Tags::MaxIterations, Tags::SelfForceTurnOnTime,
                 Tags::SelfForceTurnOnInterval>;
};

template <typename Stepper>
void test_stages(const size_t order, const bool self_force) {
  CAPTURE(order, self_force);
  using worldtube = MockWorldtube<Metavariables>;
  using element = MockElement<Metavariables>;
  const ElementId<3> element_id{0};
  const Slab slab{1., 2.};
  const TimeDelta time_step = slab.duration() / 10;
  const TimeStepId initial_id{true, 0, slab.start()};
  const Stepper time_stepper{};
  const size_t num_substeps = time_stepper.number_of_substeps();
  ActionTesting::MockRuntimeSystem<Metavariables> runner{
      {std::make_unique<Stepper>(), EventsAndTriggers{}, size_t{2}, 0.1,
       std::optional<double>{2.}, self_force ? size_t{3} : size_t{0},
       std::optional<double>{0.}, std::optional<double>{1.}}};

  // Psi(t,x) = t + p_i x^i is an exact flat-space wave solution. Along an
  // accelerating center, Psi0 = t + p_i x_p^i and dot(Psi0) = 1 + p_i v^i.
  const tnsr::i<double, 3> gradient{{0.2, -0.1, 0.05}};
  const tnsr::I<double, 3> position{{2., -0.4, 0.}};
  const tnsr::I<double, 3> velocity{{0.1, 0.2, 0.}};
  const double radius = 0.5;
  // Reconstructing traces amplifies accumulated RK roundoff by 1/radius^2.
  const Approx coefficient_approx = Approx::custom().epsilon(1.e-12).scale(1.);
  tuples::tagged_tuple_from_typelist<typename worldtube::simple_tags> initial{};
  auto& vars = get<variables_tag>(initial);
  vars.initialize(1, 0.);
  get(get<Tags::Psi0>(vars))[0] = initial_id.substep_time();
  get(get<Tags::dtPsi0>(vars))[0] = 1.;
  for (size_t i = 0; i < 3; ++i) {
    get<Tags::EvolvedPosition<3>>(vars).get(i)[0] = position.get(i);
    get<Tags::EvolvedVelocity<3>>(vars).get(i)[0] = velocity.get(i);
    get(get<Tags::Psi0>(vars))[0] += gradient.get(i) * position.get(i);
    get(get<Tags::dtPsi0>(vars))[0] += gradient.get(i) * velocity.get(i);
  }
  get<dt_variables_tag>(initial).initialize(1, 0.);
  get<history_tag>(initial) = history_tag::type{order};
  get<::Tags::TimeStepId>(initial) = initial_id;
  get<::Tags::Next<::Tags::TimeStepId>>(initial) =
      time_stepper.next_time_id(initial_id, time_step);
  get<::Tags::TimeStep>(initial) = time_step;
  get<::Tags::Time>(initial) = initial_id.substep_time();
  get<Tags::ElementFacesGridCoordinates<3>>(initial).emplace(
      element_id, tnsr::I<DataVector, 3, Frame::Grid>{size_t{1}, 0.});
  get<Tags::WorldtubeRadius>(initial) = radius;
  get<Tags::GeodesicAcceleration<3>>(initial) = tnsr::I<double, 3>{0.};
  auto& background = get<Tags::BackgroundQuantities<3>>(initial);
  auto& metric = get<gr::Tags::SpacetimeMetric<double, 3>>(background);
  auto& inverse_metric =
      get<gr::Tags::InverseSpacetimeMetric<double, 3>>(background);
  metric = tnsr::aa<double, 3>{0.};
  inverse_metric = tnsr::AA<double, 3>{0.};
  for (size_t a = 0; a < 4; ++a) {
    metric.get(a, a) = a == 0 ? -1. : 1.;
    inverse_metric.get(a, a) = a == 0 ? -1. : 1.;
  }
  get<gr::Tags::SpacetimeChristoffelSecondKind<double, 3>>(background) =
      tnsr::Abb<double, 3>{0.};
  get<gr::Tags::TraceSpacetimeChristoffelSecondKind<double, 3>>(background) =
      tnsr::A<double, 3>{0.};
  get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, 3,
                           Frame::Inertial>>(initial) = Scalar<double>{1.};
  get<Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, 3, Frame::Inertial>>(
      initial) = gradient;
  get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, 3,
                           Frame::Inertial>>(initial) = tnsr::i<double, 3>{0.};
  get<Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, 3, Frame::Inertial>>(
      initial) = tnsr::ii<double, 3>{0.};
  get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 2, 3,
                           Frame::Inertial>>(initial) = tnsr::ii<double, 3>{0.};
  ActionTesting::emplace_singleton_component_and_initialize<worldtube>(
      &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
      std::move(initial));
  ActionTesting::emplace_array_component_and_initialize<element>(
      &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
      element_id, {});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  auto& box = ActionTesting::get_databox<worldtube>(make_not_null(&runner), 0);

  for (size_t stage = 0; stage < 2 * num_substeps; ++stage) {
    CAPTURE(stage);
    const auto current_id = db::get<::Tags::TimeStepId>(box);
    const auto next_id = db::get<::Tags::Next<::Tags::TimeStepId>>(box);
    const auto pos_vel = db::get<Tags::ParticlePositionVelocity<3>>(box);
    double central_psi = current_id.substep_time();
    double speed_squared = 0.;
    for (size_t i = 0; i < 3; ++i) {
      central_psi += gradient.get(i) * pos_vel[0].get(i);
      speed_squared += pos_vel[1].get(i) * pos_vel[1].get(i);
    }
    db::mutate<psi_monopole_tag, Tags::BackgroundQuantities<3>>(
        [central_psi, speed_squared](const auto psi, const auto bg) {
          get(*psi) = central_psi;
          get(get<Tags::TimeDilationFactor>(*bg)) =
              1. / sqrt(1. - speed_squared);
        },
        make_not_null(&box));
    for (size_t action = 0;
         action < tmpl::size<typename worldtube::step_actions>::value;
         ++action) {
      ActionTesting::next_action<worldtube>(make_not_null(&runner), 0);
    }
    auto& inbox =
        ActionTesting::get_inbox_tag<element, Tags::RegularFieldInbox<3>>(
            make_not_null(&runner), element_id);
    REQUIRE(inbox.size() == 1);
    REQUIRE(inbox.count(current_id) == 1);
    const auto& sent = inbox.at(current_id);
    const auto& psi = get(get<CurvedScalarWave::Tags::Psi>(sent));
    const auto& dt_psi =
        get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(sent));
    REQUIRE(psi.size() == 10);
    REQUIRE(dt_psi.size() == 10);
    CHECK(psi[0] == approx(central_psi));
    CHECK(dt_psi[0] == approx(1.));
    for (size_t i = 0; i < 3; ++i) {
      CHECK(psi[i + 1] == approx(gradient.get(i)));
      CHECK(dt_psi[i + 1] == approx(0.));
    }
    for (size_t i = 4; i < 10; ++i) {
      CHECK(psi[i] == coefficient_approx(0.));
      CHECK(dt_psi[i] == coefficient_approx(0.));
    }
    // Matching the boundary value alone misses this bug: an erroneous Psi0
    // cancels against the reconstructed trace there, but not in the normal
    // derivative at x - x_p = (radius, 0, 0).
    CHECK(psi[0] + radius * psi[1] + radius * radius * psi[4] ==
          approx(central_psi + radius * gradient.get(0)));
    CHECK(psi[1] + 2. * radius * psi[4] == approx(gradient.get(0)));
    inbox.clear();

    // The ODE really advanced, but its new values must not enter the message
    // for current_id. The linear-field identities are preserved by RK stages.
    const auto next_pos_vel = db::get<Tags::ParticlePositionVelocity<3>>(box);
    double next_psi = next_id.substep_time();
    double next_dot_psi = 1.;
    for (size_t i = 0; i < 3; ++i) {
      next_psi += gradient.get(i) * next_pos_vel[0].get(i);
      next_dot_psi += gradient.get(i) * next_pos_vel[1].get(i);
    }
    CHECK(get(db::get<Tags::Psi0>(box))[0] == approx(next_psi));
    CHECK(get(db::get<Tags::dtPsi0>(box))[0] == approx(next_dot_psi));
    if (self_force) {
      CHECK(next_pos_vel[1].get(0) != velocity.get(0));
    }
    db::mutate<::Tags::TimeStepId, ::Tags::Next<::Tags::TimeStepId>,
               ::Tags::Time>(
        [&](const auto id, const auto next, const auto time) {
          *id = next_id;
          *next = time_stepper.next_time_id(next_id, time_step);
          *time = next_id.substep_time();
        },
        make_not_null(&box));
  }
}

SPECTRE_TEST_CASE("Unit.CurvedScalarWave.Worldtube.StageOrdering", "[Unit]") {
  register_classes_with_charm<TimeSteppers::ClassicalRungeKutta4,
                              TimeSteppers::Rk5Tsitouras>();
  for (const bool self_force : {false, true}) {
    test_stages<TimeSteppers::ClassicalRungeKutta4>(4, self_force);
    test_stages<TimeSteppers::Rk5Tsitouras>(5, self_force);
  }
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
