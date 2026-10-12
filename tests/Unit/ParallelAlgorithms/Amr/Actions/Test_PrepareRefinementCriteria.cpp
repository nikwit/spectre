// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Amr/Info.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "ParallelAlgorithms/Amr/Actions/PrepareRefinementCriteria.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Chalmers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
struct Field : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct OtherField : db::SimpleTag {
  using type = Scalar<DataVector>;
};
using Criterion = amr::Criteria::Chalmers<1, tmpl::list<Field, OtherField>>;
template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<1>;
  using const_global_cache_tags =
      tmpl::list<amr::Criteria::Tags::Criteria, amr::Tags::AmrBlocks<1>,
                 amr::Tags::Policies>;
  using simple_tags =
      tmpl::list<domain::Tags::Element<1>, domain::Tags::Mesh<1>,
                 amr::Tags::Info<1>, amr::Tags::NeighborInfo<1>,
                 amr::Tags::CriterionNeighborData<1>, Field, OtherField>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>>;
};
struct Metavariables {
  static constexpr size_t volume_dim = 1;
  using component_list = tmpl::list<Component<Metavariables>>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<amr::Criterion, tmpl::list<Criterion>>>;
  };
  struct amr {
    static constexpr bool keep_coarse_grids = false;
    static constexpr bool p_refine_only_in_event = true;
  };
};
using MyComponent = Component<Metavariables>;

void check_sibling_selection() {
  const ElementId<3> id{1, {{{2, 1}, {3, 6}, {1, 0}}}};
  for (size_t axis = 0; axis < 3; ++axis) {
    CAPTURE(axis);
    auto segments = id.segment_ids();
    const Direction<3> direction{axis, segments[axis].side_of_sibling()};
    segments[axis] = segments[axis].id_of_sibling();
    const ElementId<3> sibling{1, segments};
    const auto aligned = OrientationMap<3>::create_aligned();
    CHECK(amr::criterion_sibling(
              Element<3>{id, {{direction, Neighbors<3>{{sibling}, aligned}}}},
              axis) == sibling);
    CHECK_FALSE(amr::criterion_sibling(Element<3>{id, {}}, axis).has_value());
    // A nonconforming face is not a two-child parent candidate.
    const size_t transverse = (axis + 1) % 3;
    CHECK_FALSE(
        amr::criterion_sibling(
            Element<3>{
                id,
                {{direction,
                  Neighbors<3>{{sibling.id_of_child(transverse, Side::Lower),
                                sibling.id_of_child(transverse, Side::Upper)},
                               aligned}}}},
            axis)
            .has_value());
    // Neither a neighboring block nor another grid can be a sibling.
    for (const ElementId<3> other :
         {ElementId<3>{2, segments}, ElementId<3>{1, segments, 1}}) {
      CHECK_FALSE(
          amr::criterion_sibling(
              Element<3>{id, {{direction, Neighbors<3>{{other}, aligned}}}},
              axis)
              .has_value());
    }
    CHECK_FALSE(amr::criterion_sibling(Element<3>{ElementId<3>{1}, {}}, axis)
                    .has_value());
  }
}

void check_exchange(const bool receive_before_prepare, const bool mismatch,
                    const size_t minimum_level) {
  CAPTURE(receive_before_prepare, mismatch, minimum_level);
  auto criterion = TestHelpers::test_creation<Criterion>(
      "VariablesToMonitor: [Field, OtherField]\n"
      "AbsoluteTarget: 1.e-6\n"
      "RelativeTarget: 0.0\n"
      "SmoothnessThreshold: 0.5\n"
      "PRefineOnly: false\n"
      "CoarseningFactor: 0.1\n"
      "HCoarsening: true\n");
  CHECK(criterion.needs_sibling_data());
  criterion = serialize_and_deserialize(criterion);
  CHECK(criterion.needs_sibling_data());
  std::vector<std::unique_ptr<amr::Criterion>> criteria;
  criteria.emplace_back(std::make_unique<Criterion>(criterion));
  const amr::Limits limits{std::array{minimum_level, size_t{6}},
                           std::array{size_t{2}, size_t{20}}, std::nullopt,
                           std::nullopt, true};
  ActionTesting::MockRuntimeSystem<Metavariables> runner{
      {std::move(criteria), std::nullopt,
       amr::Policies{amr::Isotropy::Anisotropic, limits, true, true}}};
  const ElementId<1> left{0, {{{1, 0}}}};
  const ElementId<1> right{0, {{{1, 1}}}};
  const Mesh<1> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  for (const auto id : {left, right}) {
    const Element<1> element{
        id,
        {{id == left ? Direction<1>::upper_xi() : Direction<1>::lower_xi(),
          Neighbors<1>{{id == left ? right : left},
                       OrientationMap<1>::create_aligned()}}}};
    ActionTesting::emplace_component_and_initialize<MyComponent>(
        &runner, id,
        {element, mesh, amr::Info<1>{},
         std::unordered_map<ElementId<1>, amr::Info<1>>{},
         amr::CriterionNeighborData<1>{},
         Scalar<DataVector>{DataVector(8, 1.0)},
         Scalar<DataVector>{
             DataVector(8, mismatch and id == right ? 2.1 : 2.0)}});
  }
  for (size_t round = 1; round <= 2; ++round) {
    // Simulate the reset performed by AdjustDomain without changing the
    // elements, to exercise a second exchange with a new round number.
    if (round == 2) {
      for (const auto id : {left, right}) {
        auto& box =
            ActionTesting::get_databox<MyComponent>(make_not_null(&runner), id);
        db::mutate<amr::Tags::Info<1>, amr::Tags::NeighborInfo<1>,
                   amr::Tags::CriterionNeighborData<1>>(
            [](const auto info, const auto neighbors, const auto state) {
              *info = amr::Info<1>{};
              neighbors->clear();
              state->clear();
            },
            make_not_null(&box));
      }
    }
    ActionTesting::simple_action<MyComponent,
                                 amr::Actions::PrepareRefinementCriteria>(
        make_not_null(&runner), left);
    CHECK(ActionTesting::get_databox_tag<MyComponent, amr::Tags::Info<1>>(
              runner, left)
              .flags[0] == amr::Flag::Undefined);
    if (receive_before_prepare) {
      ActionTesting::invoke_queued_simple_action<MyComponent>(
          make_not_null(&runner), right);
      CHECK(ActionTesting::get_databox_tag<MyComponent, amr::Tags::Info<1>>(
                runner, right)
                .flags[0] == amr::Flag::Undefined);
    }
    ActionTesting::simple_action<MyComponent,
                                 amr::Actions::PrepareRefinementCriteria>(
        make_not_null(&runner), right);
    for (size_t iteration = 0; iteration < 12; ++iteration) {
      bool work = false;
      for (const auto id : {left, right}) {
        if (not ActionTesting::is_simple_action_queue_empty<MyComponent>(runner,
                                                                         id)) {
          ActionTesting::invoke_queued_simple_action<MyComponent>(
              make_not_null(&runner), id);
          work = true;
        }
      }
      if (not work) {
        break;
      }
    }
    for (const auto id : {left, right}) {
      CHECK(
          ActionTesting::is_simple_action_queue_empty<MyComponent>(runner, id));
      CHECK(ActionTesting::get_databox_tag<MyComponent, amr::Tags::Info<1>>(
                runner, id)
                .flags[0] == (mismatch or minimum_level == 1
                                  ? amr::Flag::DecreaseResolution
                                  : amr::Flag::Join));
      const auto& state =
          ActionTesting::get_databox_tag<MyComponent,
                                         amr::Tags::CriterionNeighborData<1>>(
              runner, id);
      CHECK(state.round == round);
      CHECK_FALSE(state.waiting);
      CHECK(state.sibling.has_value());
      const auto copy = serialize_and_deserialize(state);
      CHECK(copy.round == round);
      CHECK(copy.sibling->second.round == round);
      CHECK(copy.sibling->second.mesh == mesh);
      CHECK_ITERABLE_APPROX(copy.sibling->second.tensors.at("Field")[0],
                            DataVector(8, 1.0));
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Amr.Actions.PrepareRefinementCriteria",
                  "[Unit][ParallelAlgorithms]") {
  register_factory_classes_with_charm<Metavariables>();
  check_sibling_selection();
  for (const bool early : {false, true}) {
    for (const bool mismatch : {false, true}) {
      for (const size_t minimum : {0_st, 1_st}) {
        check_exchange(early, mismatch, minimum);
      }
    }
  }
  CHECK_THROWS_WITH((Criterion{{"Field"}, 1.e-6, 0.0, 0.5, false, 0.0, true}),
                    Catch::Matchers::ContainsSubstring("HCoarsening requires"));
  CHECK_THROWS_WITH((Criterion{{"Field"}, 1.e-6, 0.0, 0.5, true, 0.1, true}),
                    Catch::Matchers::ContainsSubstring("HCoarsening requires"));
}
