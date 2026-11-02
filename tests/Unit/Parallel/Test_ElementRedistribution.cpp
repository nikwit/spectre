// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Tags/InitialRefinementLevels.hpp"
#include "Domain/ElementMap.hpp"
#include "Evolution/DiscontinuousGalerkin/EqualRateLts/EqualRateRegions.tpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/SetupEqualRateRegions.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmGroup.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/CharmMain.tpp"
#include "Parallel/ElementRegistration.hpp"
#include "Parallel/Main.hpp"
#include "Parallel/Protocols/ElementRegistrar.hpp"
#include "Parallel/Protocols/RegistrationMetavariables.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/Amr/RedistributeElements.hpp"
#include "Time/Tags/FixedLtsRatio.hpp"
#include "Time/Tags/LtsMode.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Time.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
struct Metavariables;
template <typename Metavariables>
struct Elements;
template <typename Metavariables>
struct Registrations;
struct Members : db::SimpleTag {
  using type = std::set<ElementId<1>>;
};
struct Payload : db::SimpleTag {
  using type = std::vector<double>;
};
struct ExpectedProc : db::SimpleTag {
  using type = int;
};
struct PendingData {
  using type = std::map<size_t, int>;
  using temporal_id = size_t;
};
struct TestRegions {
  using creation_tags = tmpl::list<>;
  static std::unordered_map<std::string, size_t> regions() {
    return {{"Left", 0}, {"Right", 1}};
  }
  static bool is_in_region(const size_t region, const ElementId<1>& id) {
    const auto i = id.segment_ids()[0].index();
    return region == 0 ? i >= 1 and i <= 3 : i >= 5 and i <= 6;
  }
  void pup(PUP::er&) {}
};
using SectionTag =
    Parallel::Tags::Section<Elements<Metavariables>,
                            evolution::dg::Tags::EqualRateRegionId>;
struct ReductionResults : db::SimpleTag {
  using type = size_t;
};

struct ReceiveSectionSum {
  template <typename ParallelComponent, typename DbTags, typename Meta>
  static void apply(db::DataBox<DbTags>& box, Parallel::GlobalCache<Meta>&,
                    const ElementId<1>&, const size_t label, const size_t count,
                    const size_t sum) {
    SPECTRE_PARALLEL_REQUIRE(count == (label == 0 ? 3 : 2));
    SPECTRE_PARALLEL_REQUIRE(sum == (label == 0 ? 9 : 13));
    db::mutate<ReductionResults>([](const auto n) { ++*n; },
                                 make_not_null(&box));
  }
};
struct ReduceSection {
  template <typename DbTags, typename... InboxTags, typename Meta,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Meta>& cache, const ElementId<1>& id, ActionList,
      const ParallelComponent*) {
    auto& section = db::get_mutable_reference<SectionTag>(make_not_null(&box));
    if (section) {
      auto& array = Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::ReductionData<
          Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
          Parallel::ReductionDatum<size_t, funcl::Plus<>>,
          Parallel::ReductionDatum<size_t, funcl::Plus<>>>
          data{section->id().label, 1, id.segment_ids()[0].index() + 1};
      Parallel::contribute_to_reduction<ReceiveSectionSum>(
          std::move(data), array[id], array, make_not_null(&*section));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct RecordMembership {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const int /*index*/, const ElementId<1>& id,
                    const bool insert) {
    db::mutate<Members>(
        [&](const auto members) {
          if (insert) {
            SPECTRE_PARALLEL_REQUIRE(members->insert(id).second);
          } else {
            SPECTRE_PARALLEL_REQUIRE(members->erase(id) == 1);
          }
        },
        make_not_null(&box));
  }
};
struct Registrar : tt::ConformsTo<Parallel::protocols::ElementRegistrar> {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void perform_registration(const db::DataBox<DbTags>& /*box*/,
                                   Parallel::GlobalCache<Metavariables>& cache,
                                   const ElementId<1>& id) {
    Parallel::simple_action<RecordMembership>(
        Parallel::get_parallel_component<Registrations<Metavariables>>(
            cache)[Parallel::my_proc<int>(cache)],
        id, true);
  }
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void perform_deregistration(
      const db::DataBox<DbTags>& /*box*/,
      Parallel::GlobalCache<Metavariables>& cache, const ElementId<1>& id) {
    Parallel::simple_action<RecordMembership>(
        Parallel::get_parallel_component<Registrations<Metavariables>>(
            cache)[Parallel::my_proc<int>(cache)],
        id, false);
  }
};
struct Initialize {
  using simple_tags = tmpl::list<
      Payload, ExpectedProc, domain::Tags::Mesh<1>,
      domain::Tags::Coordinates<1, Frame::Grid>, domain::Tags::Element<1>,
      domain::Tags::InitialRefinementLevels<1>, Tags::LtsMode,
      Tags::FixedLtsRatio, Tags::TimeStep, SectionTag, ReductionResults,
      evolution::dg::Tags::ConcreteEqualRateRegions<1,
                                                    tmpl::list<TestRegions>>>;
  using compute_tags = tmpl::list<
      evolution::dg::Tags::EqualRateRegionsRef<1, tmpl::list<TestRegions>>>;
  using inbox_tags = tmpl::list<PendingData>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache, const ElementId<1>& id,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    const size_t index = id.segment_ids()[0].index();
    const Mesh<1> mesh{index + 3, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
    db::mutate<Payload, ExpectedProc, domain::Tags::Mesh<1>,
               domain::Tags::Coordinates<1, Frame::Grid>>(
        [&](const auto payload, const auto proc, const auto stored_mesh,
            const auto coords) {
          *payload = std::vector<double>(index + 3, index + 0.25);
          *proc = 0;
          *stored_mesh = mesh;
          *coords =
              ElementMap<1, Frame::Grid>{
                  id,
                  Parallel::get<domain::Tags::Domain<1>>(cache).blocks()[0]}(
                  logical_coordinates(mesh));
        },
        make_not_null(&box));
    tuples::get<PendingData>(inboxes).emplace(99, static_cast<int>(index));
    db::mutate<domain::Tags::Element<1>,
               domain::Tags::InitialRefinementLevels<1>, Tags::LtsMode,
               Tags::FixedLtsRatio, Tags::TimeStep>(
        [&](const auto element, const auto levels, const auto mode,
            const auto ratio, const auto step) {
          *levels = {{3}};
          *element = domain::create_initial_element(
              id, Parallel::get<domain::Tags::Domain<1>>(cache).blocks(),
              *levels);
          *mode = LtsMode::Conservative;
          *ratio = 32;
          *step = Slab{0.0, 1.0}.duration() / 32;
        },
        make_not_null(&box));
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};
struct Register {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache, const ElementId<1>& id,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    Parallel::register_element<ParallelComponent>(box, cache, id);
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
struct RestoreSections {
  template <typename DbTags, typename... InboxTags, typename Meta,
            typename Index, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Meta>& cache, const Index& id, ActionList actions,
      const ParallelComponent* component) {
    // The fresh test also visits Restart. Only a real checkpoint restart
    // discards saved cookies: in the fresh run their trees are still alive.
    if (std::getenv("SPECTRE_RESTART_TEST") != nullptr) {
      return evolution::dg::Initialization::Actions::CreateEqualRateSections<
          1, true>::apply(box, inboxes, cache, id, actions, component);
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
struct SaveTarget {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<1>& /*id*/, ActionList /*meta*/,
      const ParallelComponent* /*meta*/) {
    const int target =
        db::get<amr::redistribution::Tags::TargetProc>(box).value();
    db::mutate<ExpectedProc>([target](const auto proc) { *proc = target; },
                             make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
template <size_t ExpectedReductions>
struct Check {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache, const ElementId<1>& id,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    const auto index = id.segment_ids()[0].index();
    SPECTRE_PARALLEL_REQUIRE(Parallel::number_of_procs<size_t>(cache) >= 2);
    SPECTRE_PARALLEL_REQUIRE(Parallel::my_proc<int>(cache) ==
                             db::get<ExpectedProc>(box));
    SPECTRE_PARALLEL_REQUIRE(db::get<Payload>(box) ==
                             std::vector<double>(index + 3, index + 0.25));
    SPECTRE_PARALLEL_REQUIRE(tuples::get<PendingData>(inboxes).at(99) ==
                             static_cast<int>(index));
    SPECTRE_PARALLEL_REQUIRE(
        not db::get<amr::redistribution::Tags::TargetProc>(box));
    SPECTRE_PARALLEL_REQUIRE(db::get<ReductionResults>(box) ==
                             ExpectedReductions);
    SPECTRE_PARALLEL_REQUIRE(db::get<Tags::FixedLtsRatio>(box) == 32);
    SPECTRE_PARALLEL_REQUIRE(db::get<Tags::TimeStep>(box).fraction() ==
                             Rational(1, 32));
    const auto& section = db::get<SectionTag>(box);
    SPECTRE_PARALLEL_REQUIRE(
        section.has_value() ==
        (TestRegions::is_in_region(0, id) or TestRegions::is_in_region(1, id)));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
struct InitializeMembers {
  using simple_tags = tmpl::list<Members>;
  template <typename... Args>
  static Parallel::iterable_action_return_t apply(Args&&... /*unused*/) {
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};
struct CheckMembers {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache, const int index,
      ActionList /*meta*/, const ParallelComponent* /*meta*/) {
    amr::redistribution::Tags::Costs<1>::type records{};
    for (const auto& id : initial_element_ids(0, std::array<size_t, 1>{{3}})) {
      records.emplace(
          id,
          std::pair{static_cast<double>(id.segment_ids()[0].index() + 3), 0});
    }
    const auto plan = amr::redistribution::make_plan(
        Parallel::get<domain::Tags::Domain<1>>(cache).blocks(), records,
        Parallel::number_of_procs<size_t>(cache),
        cache.get_resource_info().procs_to_ignore(), true);
    std::set<ElementId<1>> expected{};
    for (const auto& [id, target] : plan) {
      if (target == index) {
        expected.insert(id);
      }
    }
    SPECTRE_PARALLEL_REQUIRE(db::get<Members>(box) == expected);
    // This test deliberately starts every element on proc 0. Ensure movement.
    SPECTRE_PARALLEL_REQUIRE(plan.at(records.begin()->first) >= 0);
    SPECTRE_PARALLEL_REQUIRE(
        std::any_of(plan.begin(), plan.end(),
                    [](const auto& entry) { return entry.second != 0; }));
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};

template <typename Metavariables>
struct Elements {
  using metavariables = Metavariables;
  using chare_type = Parallel::Algorithms::Array;
  static constexpr bool checkpoint_data = true;
  using array_index = ElementId<1>;
  using array_allocation_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<Initialize>>,
      Parallel::PhaseActions<Parallel::Phase::Register,
                             tmpl::list<Register,
                                        evolution::dg::Initialization::Actions::
                                            CreateEqualRateSections<1>,
                                        Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Execute,
          tmpl::list<Check<0>, ReduceSection, ReduceSection, ReduceSection,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<Parallel::Phase::CollectElementDistribution,
                             tmpl::list<evolution::dg::Initialization::Actions::
                                            ResetEqualRateSections<1>,
                                        amr::redistribution::Actions::Collect>>,
      Parallel::PhaseActions<
          Parallel::Phase::RedistributeElements,
          tmpl::list<SaveTarget, amr::redistribution::Actions::Migrate>>,
      Parallel::PhaseActions<Parallel::Phase::UpdateSections,
                             tmpl::list<evolution::dg::Initialization::Actions::
                                            CreateEqualRateSections<1>,
                                        Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Restart,
          tmpl::list<RestoreSections, Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<Check<6>, ReduceSection, ReduceSection, ReduceSection,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Solve,
          tmpl::list<Check<12>, ReduceSection, ReduceSection, ReduceSection,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<Parallel::Phase::Cleanup,
                             tmpl::list<Check<18>,
                                        evolution::dg::Initialization::Actions::
                                            ResetEqualRateSections<1>,
                                        Parallel::Actions::TerminatePhase>>>;
  static void allocate_array(Parallel::CProxy_GlobalCache<Metavariables>& cache,
                             const tuples::TaggedTuple<>& items,
                             const tuples::TaggedTuple<>& /*unused*/,
                             const std::unordered_set<size_t>& /*ignored*/) {
    auto& proxy = Parallel::get_parallel_component<Elements>(
        *Parallel::local_branch(cache));
    for (const auto& id : initial_element_ids(0, std::array<size_t, 1>{{3}})) {
      proxy[id].insert(cache, items, 0);
    }
    proxy.doneInserting();
  }
  static void execute_next_phase(
      const Parallel::Phase phase,
      Parallel::CProxy_GlobalCache<Metavariables>& cache) {
    Parallel::get_parallel_component<Elements>(*Parallel::local_branch(cache))
        .start_phase(phase);
  }
};
template <typename Metavariables>
struct Registrations {
  using metavariables = Metavariables;
  using chare_type = Parallel::Algorithms::Group;
  static constexpr bool checkpoint_data = true;
  using simple_tags_from_options = tmpl::list<>;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<Parallel::Phase::Initialization,
                                        tmpl::list<InitializeMembers>>,
                 Parallel::PhaseActions<Parallel::Phase::Testing,
                                        tmpl::list<CheckMembers>>>;
  static void execute_next_phase(
      const Parallel::Phase phase,
      Parallel::CProxy_GlobalCache<Metavariables>& cache) {
    Parallel::get_parallel_component<Registrations>(
        *Parallel::local_branch(cache))
        .start_phase(phase);
  }
};
struct Metavariables {
  void pup(PUP::er& /*p*/) {}
  static constexpr size_t volume_dim = 1;
  struct system {};
  struct amr {
    using element_array = Elements<Metavariables>;
  };
  struct registration
      : tt::ConformsTo<Parallel::protocols::RegistrationMetavariables> {
    using element_registrars =
        tmpl::map<tmpl::pair<Elements<Metavariables>, tmpl::list<Registrar>>>;
  };
  using component_list =
      tmpl::list<Elements<Metavariables>, Registrations<Metavariables>,
                 ::amr::redistribution::Component<Metavariables>>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<DomainCreator<1>,
                             tmpl::list<domain::creators::Rectilinear<1>>>>;
  };
  static constexpr Options::String help =
      "Test actual redistribution, inboxes, registration and checkpoint "
      "restart.";
  static constexpr std::array default_phase_order{
      Parallel::Phase::Initialization,
      Parallel::Phase::Register,
      Parallel::Phase::Execute,
      Parallel::Phase::CollectElementDistribution,
      Parallel::Phase::PlanElementDistribution,
      Parallel::Phase::RedistributeElements,
      Parallel::Phase::UpdateSections,
      Parallel::Phase::Testing,
      Parallel::Phase::WriteCheckpoint,
      Parallel::Phase::Restart,
      Parallel::Phase::Solve,
      Parallel::Phase::Cleanup,
      Parallel::Phase::Exit};
};
void register_types() {
  domain::creators::register_derived_with_charm();
  register_factory_classes_with_charm<Metavariables>();
}
}  // namespace
extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<Metavariables>();
  Parallel::charmxx::register_init_node_and_proc({&register_types}, {});
}
