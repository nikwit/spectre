// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/ObservationBox.hpp"
#include "Domain/Block.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/Domain.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/ResetSpecAmrHistory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Bbh/PhaseControl/InitialPRefinement.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Hypercube.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Tags/Criteria.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using GridJacobian = domain::Tags::Jacobian<3, Frame::Grid, Frame::Inertial>;
using FilterTag = Filters::Tags::SpectralFilter<
    3, ylm::TensorYlm::filter_detail::gh_spacetime_vars_list>;
GridJacobian::type identity_grid_jacobian(const size_t n) {
  GridJacobian::type jac(n, 0.0);
  for (size_t d = 0; d < 3; ++d) {
    jac.get(d, d) = 1.0;
  }
  return jac;
}

using Controller = gh::bbh::phase_control::InitialPRefinement;
struct Metavariables {
  static constexpr size_t volume_dim = 3;
  using component_list = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<PhaseChange, tmpl::list<Controller>>,
                  tmpl::pair<amr::Criterion, tmpl::list<gh::SpecSpectral>>>;
  };
};

Domain<3> make_domain() {
  std::vector<Block<3>> blocks{};
  for (size_t i = 0; i < 2; ++i) {
    blocks.emplace_back(
        domain::make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
            domain::CoordinateMaps::Identity<3>{}),
        i, DirectionMap<3, BlockNeighbors<3>>{}, i == 0 ? "Inner" : "Outer");
  }
  return Domain<3>{std::move(blocks), {}, {{"WaveZone", {"Outer"}}}};
}

struct RetainedState : db::SimpleTag {
  using type = std::vector<double>;
};
template <typename Metavars>
struct RestartElement {
  using metavariables = Metavars;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using const_global_cache_tags = tmpl::list<amr::Criteria::Tags::Criteria>;
  using simple_tags = tmpl::list<gh::spec_amr::HistoryTag, RetainedState>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>>;
};
struct RestartMetavariables : Metavariables {
  using component_list = tmpl::list<RestartElement<RestartMetavariables>>;
  struct amr {
    using element_array = RestartElement<RestartMetavariables>;
    using criteria_change_action = gh::Actions::ResetSpecAmrHistory;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.RestartOverlay",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<RestartMetavariables>();
  using Component = RestartElement<RestartMetavariables>;
  using Criteria = amr::Criteria::Tags::Criteria;
  using History = gh::spec_amr::HistoryTag;
  const auto criteria = [](const double target) {
    Criteria::type result{};
    result.push_back(std::make_unique<gh::SpecSpectral>(
        target, std::array<double, 2>{{1.e-3, 1.e-3}}, true, size_t{7}, 2, -1.,
        0., false, 1.));
    return result;
  };
  ActionTesting::MockRuntimeSystem<RestartMetavariables> runner{
      {criteria(1.e-4)}, {}, {1, 1}};
  gh::spec_amr::History history{};
  history.calls = 17;
  history.time = 1000.;
  history.spectral_target = 1.e-4;
  history.radius_release_time = 900.;
  history.log10_target_shift.fill(-0.5);
  history.last_change.fill(15);
  history.last_change_time.fill(999.);
  history.extents.fill(12);
  history.constraints[0].fill(1.e-5);
  history.decision.fill(amr::Flag::DecreaseResolution);
  const std::vector<double> retained{1., 2., 3.};
  for (size_t node = 0; node < 2; ++node) {
    ActionTesting::emplace_array_component_and_initialize<Component>(
        &runner, ActionTesting::NodeId{node}, ActionTesting::LocalCoreId{0},
        static_cast<int>(node), {history, retained});
    ActionTesting::cache<Component>(runner, static_cast<int>(node))
        .overlay_cache_data({criteria(1.e-4)});
  }
  // An identical configuration preserves all AMR history.
  for (int id = 0; id < 2; ++id) {
    CHECK(ActionTesting::is_simple_action_queue_empty<Component>(runner, id));
    CHECK(serialize(ActionTesting::get_databox_tag<Component, History>(
              runner, id)) == serialize(history));
  }
  // Updating node 1 must not broadcast. Updating node 0 then broadcasts once
  // to every element, including those on another node.
  for (const int node : {1, 0}) {
    ActionTesting::cache<Component>(runner, node)
        .overlay_cache_data({criteria(1.e-5)});
    for (int id = 0; id < 2; ++id) {
      CHECK(ActionTesting::number_of_queued_simple_actions<Component>(
                runner, id) == (node == 0 ? 1 : 0));
    }
  }
  for (int id = 0; id < 2; ++id) {
    ActionTesting::invoke_queued_simple_action<Component>(
        make_not_null(&runner), id);
    CHECK(serialize(ActionTesting::get_databox_tag<Component, History>(
              runner, id)) == serialize(gh::spec_amr::History{}));
    CHECK(ActionTesting::get_databox_tag<Component, RetainedState>(
              runner, id) == retained);
    ActionTesting::cache<Component>(runner, id)
        .overlay_cache_data({criteria(1.e-5)});
    CHECK(ActionTesting::is_simple_action_queue_empty<Component>(runner, id));
  }
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.InitialDecisions",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  Settings settings{1.e-4, {{1.e-3, 1.e-3}}, true, 7, 2, -1., 10.};
  History history{};
  history.initializing = true;
  Diagnostics diagnostics{};
  for (auto& group : diagnostics.spectra) {
    group.fill(Spectrum{1.e-7, 2.e-7, 0.3, 0., false});
  }
  auto extents = std::array<size_t, 3>{{10, 10, 10}};
  const auto down = make_array<3>(amr::Flag::DecreaseResolution);
  const auto none = make_array<3>(amr::Flag::DoNothing);
  // Startup can coarsen on the first call and repeatedly at the same time.
  CHECK(decide_p(make_not_null(&history), diagnostics, extents, 0., settings) ==
        down);
  CHECK(history.calls == 0);
  CHECK(decide_p(make_not_null(&history), diagnostics, extents, 0., settings) ==
        down);
  CHECK(history.initial_meshes.size() == 1);
  extents.fill(9);
  CHECK(decide_p(make_not_null(&history), diagnostics, extents, 0., settings) ==
        down);
  auto restored = serialize_and_deserialize(history);
  CHECK(restored.initial_meshes == history.initial_meshes);
  // Reaching an earlier mesh ends a cycle, rather than oscillating
  // indefinitely.
  extents.fill(10);
  CHECK(decide_p(make_not_null(&restored), diagnostics, extents, 0.,
                 settings) == none);
  CHECK(restored.initial_cycle_detected);
  CHECK(restored.calls == 0);
  CHECK(restored.constraints == History{}.constraints);
  // Ordinary evolution retains the original first-call and time cooldown rules.
  history.initializing = false;
  CHECK(decide_p(make_not_null(&history), diagnostics, extents, 0., settings) ==
        none);
  CHECK(history.calls == 1);
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.JunkGate",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  using History = gh::spec_amr::HistoryTag;
  using Metric = gr::Tags::SpacetimeMetric<DataVector, 3>;
  using Pi = gh::Tags::Pi<DataVector, 3>;
  using Phi = gh::Tags::Phi<DataVector, 3>;
  using GridCoords = domain::Tags::Coordinates<3, Frame::Grid>;
  using Inverse =
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>;
  const Mesh<3> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const size_t n = mesh.number_of_grid_points();
  Metric::type metric(n, 0.0);
  get<0, 0>(metric) = -1.;
  get<1, 1>(metric) = get<2, 2>(metric) = get<3, 3>(metric) = 1.;
  Inverse::type inv_jac(n, 0.0);
  get<0, 0>(inv_jac) = get<1, 1>(inv_jac) = get<2, 2>(inv_jac) = 1.;
  GridCoords::type grid_coords(n, 0.0);
  get<0>(grid_coords) = 10.;
  get<0>(grid_coords)[0] =
      20.;  // The farthest point, not center or mean radius.
  auto box = db::create<tmpl::list<
      domain::Tags::Domain<3>, domain::Tags::Mesh<3>, Metric, Pi, Phi, Inverse,
      GridCoords, domain::Tags::Coordinates<3, Frame::Inertial>, ::Tags::Time,
      gh::gauges::Tags::GaugeCondition, gh::Tags::ConstraintGamma2, History,
      GridJacobian, FilterTag>>(
      make_domain(), mesh, metric, Pi::type(n, 0.), Phi::type(n, 0.), inv_jac,
      grid_coords, tnsr::I<DataVector, 3>(n, 0.), 219.99,
      std::unique_ptr<gh::gauges::GaugeCondition>(new gh::gauges::Harmonic{}),
      Scalar<DataVector>(n, 1.), gh::spec_amr::History{},
      identity_grid_jacobian(n),
      FilterTag::type(new Filters::Hypercube<
                      3, ylm::TensorYlm::filter_detail::gh_spacetime_vars_list>(
          36, true, std::nullopt, false, false, 1, std::nullopt)));
  auto observation = make_observation_box<tmpl::list<>>(make_not_null(&box));
  Parallel::GlobalCache<Metavariables> cache{};
  const auto criterion =
      TestHelpers::test_factory_creation<amr::Criterion, gh::SpecSpectral>(
          "SpecSpectral:\n"
          "  RelativeTruncationTarget: 1.e-4\n"
          "  ConstraintTargets: [1.e-3, 1.e-3]\n"
          "  ForceDropPileUpModes: true\n"
          "  MinimumModes: 7\n"
          "  CallsBetweenChanges: 2\n"
          "  TimeBetweenChanges: 1\n"
          "  TruncationErrorAdjustmentLimit: -1\n"
          "  StartTime:\n"
          "    Time: 100\n"
          "    Blocks: [WaveZone]\n"
          "    RadiusPlusTime: 200\n"
          "  Verbose: false\n");
  db::mutate<History>(
      [&mesh](const auto h) {
        h->calls = 3;
        h->extents = mesh.extents().indices();
      },
      make_not_null(&box));
  const auto none = make_array<3>(amr::Flag::DoNothing);
  const auto down = make_array<3>(amr::Flag::DecreaseResolution);
  // Inner region remains eligible even while the outer region is frozen.
  CHECK(criterion->evaluate(observation, cache, ElementId<3>{0}) == down);
  CHECK(criterion->evaluate_and_update(make_not_null(&observation), cache,
                                       ElementId<3>{1}) == none);
  CHECK(db::get<History>(box).calls == 4);
  CHECK(db::get<History>(box).radius_release_time == 220.);
  CHECK(criterion->evaluate_and_update(make_not_null(&observation), cache,
                                       ElementId<3>{1}) == none);
  CHECK(db::get<History>(box).calls == 4);
  db::mutate<::Tags::Time>([](const auto t) { *t = 220.; },
                           make_not_null(&box));
  observation.reset();
  CHECK(serialize_and_deserialize(criterion)->evaluate(
            observation, cache, ElementId<3>{1}) == down);
  // Changing the sampled coordinates later does not move the cached threshold.
  db::mutate<GridCoords>([](const auto coords) { get<0>(*coords) = 30.; },
                         make_not_null(&box));
  observation.reset();
  CHECK(criterion->evaluate(observation, cache, ElementId<3>{1}) == down);
  // Initialization bypasses both the global and radial gate.
  db::mutate<History, ::Tags::Time, GridCoords>(
      [&grid_coords](const auto h, const auto t, const auto coords) {
        *h = gh::spec_amr::History{};
        h->initializing = true;
        *t = 0.;
        *coords = grid_coords;
      },
      make_not_null(&box));
  observation.reset();
  CHECK(criterion->evaluate_and_update(make_not_null(&observation), cache,
                                       ElementId<3>{1}) == down);
  CHECK(db::get<History>(box).calls == 0);
  CHECK(db::get<History>(box).radius_release_time == 220.);
  // Bad group names must not silently disable the gate.
  CHECK_THROWS_WITH(
      gh::SpecSpectral::RadiusDependentStart(0., {"Typo"}, 200.)
          .release_time(db::get<domain::Tags::Domain<3>>(box), 0, grid_coords),
      Catch::Matchers::ContainsSubstring("The block or group 'Typo'"));
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.StartupPhases",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  using namespace gh::bbh::phase_control::initial_p_refinement;
  using Data = tuples::TaggedTuple<StateTag, ChangedElements, LimitedDirections,
                                   CyclingElements>;
  const Parallel::GlobalCache<Metavariables> cache{};
  const auto controller =
      TestHelpers::test_creation<Controller>("MaxIterations: 2");
  Data data{};
  controller.initialize_phase_data_impl(make_not_null(&data));
  using Parallel::Phase;
  const auto jump = [&](const Phase phase) -> std::optional<Phase> {
    const auto result = controller.arbitrate_phase_change_impl(
        make_not_null(&data), phase, cache);
    return result ? std::optional{result->first} : std::nullopt;
  };
  CHECK(jump(Phase::Initialization) == std::nullopt);
  CHECK(jump(Phase::InitializeInitialDataDependentQuantities) ==
        Phase::EvaluateAmrCriteria);
  CHECK(jump(Phase::EvaluateAmrCriteria) == Phase::CheckInitialPRefinement);
  tuples::get<ChangedElements>(data) = 3;
  CHECK(jump(Phase::CheckInitialPRefinement) == Phase::AdjustDomain);
  CHECK(jump(Phase::AdjustDomain) == Phase::CheckDomain);
  tuples::get<ChangedElements>(data) = 3;
  CHECK(jump(Phase::CheckDomain) == Phase::RegisterWithElementDataReader);
  CHECK(jump(Phase::RegisterWithElementDataReader) == std::nullopt);
  CHECK(jump(Phase::ImportInitialData) == std::nullopt);
  CHECK(jump(Phase::InitializeInitialDataDependentQuantities) ==
        Phase::EvaluateAmrCriteria);
  CHECK(jump(Phase::EvaluateAmrCriteria) == Phase::CheckInitialPRefinement);
  tuples::get<ChangedElements>(data) = 1;
  CHECK(jump(Phase::CheckInitialPRefinement) == Phase::AdjustDomain);
  // Even the final permitted change must be followed by reimport, before exit.
  CHECK(jump(Phase::AdjustDomain) == Phase::CheckDomain);
  tuples::get<ChangedElements>(data) = 1;
  CHECK(jump(Phase::CheckDomain) == Phase::RegisterWithElementDataReader);
  CHECK(jump(Phase::InitializeInitialDataDependentQuantities) ==
        Phase::Register);
  CHECK(tuples::get<StateTag>(data).finished);
  CHECK(jump(Phase::EvaluateAmrCriteria) == std::nullopt);
  // A central floor cancels local coarsening proposals. The loop stops on
  // realized meshes, even though the pre-planning reduction reported changes.
  controller.initialize_phase_data_impl(make_not_null(&data));
  CHECK(jump(Phase::InitializeInitialDataDependentQuantities) ==
        Phase::EvaluateAmrCriteria);
  tuples::get<ChangedElements>(data) = 2;
  tuples::get<LimitedDirections>(data) = 5;
  CHECK(jump(Phase::CheckInitialPRefinement) == Phase::AdjustDomain);
  CHECK(jump(Phase::AdjustDomain) == Phase::CheckDomain);
  CHECK(tuples::get<ChangedElements>(data) == 0);
  CHECK(jump(Phase::CheckDomain) == Phase::Register);
  CHECK(tuples::get<StateTag>(data).iteration == 1);
  CHECK(tuples::get<LimitedDirections>(data) == 0);
  CHECK(serialize_and_deserialize(controller)
            .arbitrate_phase_change_impl(make_not_null(&data), Phase::Restart,
                                         cache) == std::nullopt);
}
