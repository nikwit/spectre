// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Block.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/Topology.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ApplyTensorYlmFilter.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.tpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/CharmMain.tpp"
#include "Parallel/Main.hpp"
#include "ParallelAlgorithms/Amr/Actions/PlanMesh.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/System/ParallelInfo.hpp"

namespace mesh_plan_test {
struct Meta;
template <typename Metavariables>
struct Elements;
struct VariablesTag {
  using tags_list = ylm::TensorYlm::filter_detail::gh_spacetime_vars_list;
};
using Filter = Filters::SphericalShell<VariablesTag::tags_list>;
using FilterTag = Filters::Tags::SpectralFilter<3, VariablesTag::tags_list>;
struct Applied : db::SimpleTag {
  using type = size_t;
};

ElementId<3> element_id(const size_t object, const size_t i) {
  return ElementId<3>{
      object, std::array{SegmentId{2, i}, SegmentId{0, 0}, SegmentId{0, 0}}};
}
Mesh<3> shell_mesh(const size_t l) {
  return Mesh<3>{
      std::array{5_st, l + 1, 2 * l + 1},
      std::array{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
                 Spectral::Basis::SphericalHarmonic},
      std::array{Spectral::Quadrature::GaussLobatto,
                 Spectral::Quadrature::Gauss,
                 Spectral::Quadrature::Equiangular}};
}
struct DomainTag : domain::Tags::Domain<3> {
  using option_tags = tmpl::list<>;
  static constexpr bool pass_metavariables = false;
  static ::Domain<3> create_from_options() {
    std::unordered_map<std::string, ExcisionSphere<3>> spheres;
    for (size_t b = 0; b < 2; ++b) {
      spheres.emplace(b == 0 ? "ExcisionSphereA" : "ExcisionSphereB",
                      ExcisionSphere<3>{1.,
                                        tnsr::I<double, 3, Frame::Grid>{0.},
                                        {{b, Direction<3>::lower_xi()}}});
    }
    return ::Domain<3>{std::vector<Block<3>>{}, std::move(spheres)};
  }
};
struct PoliciesTag : amr::Tags::Policies {
  using option_tags = tmpl::list<>;
  static constexpr bool pass_metavariables = false;
  static amr::Policies create_from_options() {
    return {amr::Isotropy::Anisotropic,
            amr::Limits{},
            false,
            true,
            {"ExcisionSphereA", "ExcisionSphereB"}};
  }
};
struct Initialize {
  using simple_tags =
      tmpl::list<domain::Tags::Element<3>, domain::Tags::Mesh<3>,
                 amr::Tags::Info<3>, amr::Tags::NeighborInfo<3>, FilterTag,
                 Applied>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>&,
      Parallel::GlobalCache<Metavariables>&, const ElementId<3>& id, ActionList,
      const ParallelComponent*) {
    const auto i = id.segment_id(0).index();
    const auto b = id.block_id();
    DirectionMap<3, Neighbors<3>> ns;
    if (i > 0)
      ns.emplace(Direction<3>::lower_xi(),
                 Neighbors<3>{std::unordered_set{element_id(b, i - 1)},
                              OrientationMap<3>::create_aligned()});
    if (i < 3)
      ns.emplace(Direction<3>::upper_xi(),
                 Neighbors<3>{std::unordered_set{element_id(b, i + 1)},
                              OrientationMap<3>::create_aligned()});
    auto mesh = shell_mesh(5);
    auto proposed = shell_mesh(i == 0 ? (b == 0 ? 8 : 12) : 5);
    ::Initialization::mutate_assign<simple_tags>(
        make_not_null(&box),
        Element<3>{
            id, std::move(ns),
            std::array{domain::Topology::I1, domain::Topology::S2Colatitude,
                       domain::Topology::S2Longitude}},
        mesh, amr::planning::p_decision(mesh, proposed),
        std::unordered_map<ElementId<3>, amr::Info<3>>{},
        std::unique_ptr<Filters::Filter<3, VariablesTag::tags_list>>{
            std::make_unique<Filter>(2, std::nullopt, std::nullopt, true,
                                     std::nullopt, true, false, std::nullopt,
                                     std::nullopt)},
        size_t{0});
    return {Parallel::AlgorithmExecution::Halt, std::nullopt};
  }
};
// The planning protocol is independent of the field projection driver. Check
// exact targets and all neighbor decisions here; AdjustDomain's multi-mode
// projection/application path is tested separately in Test_AdjustDomain.cpp.
struct CheckAndApply {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>&,
                    const ElementId<3>& id) {
    const size_t applied = db::get<Applied>(box);
    const auto& info = db::get<amr::Tags::Info<3>>(box);
    const size_t l = (id.block_id() == 0 ? 8 : 12) - applied;
    SPECTRE_PARALLEL_REQUIRE(info.new_mesh == shell_mesh(l));
    // The first round catches up by 3 or 7 modes. No fields were modified
    // during collection or delivery.
    SPECTRE_PARALLEL_REQUIRE(db::get<domain::Tags::Mesh<3>>(box) ==
                             shell_mesh(applied == 0 ? 5 : l + 1));
    const auto& ns = db::get<amr::Tags::NeighborInfo<3>>(box);
    for (const auto& [direction, neighbors] :
         db::get<domain::Tags::Element<3>>(box).neighbors()) {
      (void)direction;
      for (const auto& neighbor : neighbors.ids()) {
        SPECTRE_PARALLEL_REQUIRE(ns.at(neighbor).new_mesh == shell_mesh(l));
      }
    }
    ::Initialization::mutate_assign<tmpl::list<domain::Tags::Mesh<3>, Applied>>(
        make_not_null(&box), info.new_mesh, applied + 1);
  }
};
struct PrepareSecondRound {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>&,
                    const ElementId<3>& id) {
    SPECTRE_PARALLEL_REQUIRE(db::get<Applied>(box) == 1);
    const auto& mesh = db::get<domain::Tags::Mesh<3>>(box);
    const auto target = shell_mesh((id.block_id() == 0 ? 8 : 12) - 1);
    db::mutate<amr::Tags::Info<3>>(
        [&](const auto info) {
          *info = amr::planning::p_decision(mesh, target);
        },
        make_not_null(&box));
  }
};
struct CheckFinished {
  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>&,
                    const ElementId<3>&) {
    SPECTRE_PARALLEL_REQUIRE(db::get<Applied>(box) == 2);
  }
};
template <typename Metavariables>
struct Elements {
  using metavariables = Metavariables;
  using chare_type = Parallel::Algorithms::Array;
  using array_index = ElementId<3>;
  static constexpr bool checkpoint_data = true;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<Parallel::Phase::Initialization,
                                        tmpl::list<Initialize>>>;
  using simple_tags_from_options = tmpl::list<>;
  using array_allocation_tags = tmpl::list<>;
  static void allocate_array(Parallel::CProxy_GlobalCache<Metavariables>& cache,
                             const tuples::TaggedTuple<>&,
                             const tuples::TaggedTuple<>& = {},
                             const std::unordered_set<size_t>& ignored = {}) {
    SPECTRE_PARALLEL_REQUIRE(ignored.empty());
    auto& array = Parallel::get_parallel_component<Elements>(
        *Parallel::local_branch(cache));
    for (size_t b = 0; b < 2; ++b)
      for (size_t i = 0; i < 4; ++i)
        array[element_id(b, i)].insert(
            cache, tuples::TaggedTuple<>{},
            (b * 4 + i) % static_cast<size_t>(sys::number_of_procs()));
    array.doneInserting();
  }
  static void execute_next_phase(
      const Parallel::Phase phase,
      Parallel::CProxy_GlobalCache<Metavariables>& cache) {
    auto& array = Parallel::get_parallel_component<Elements>(
        *Parallel::local_branch(cache));
    array.start_phase(phase);
    if (phase == Parallel::Phase::Testing)
      Parallel::simple_action<PrepareSecondRound>(array);
    if (phase == Parallel::Phase::CheckDomain)
      Parallel::simple_action<CheckFinished>(array);
  }
};
}  // namespace mesh_plan_test

namespace amr {
// Use the real planner actions with a small coordinator, avoiding dependencies
// on evolution criteria and initial-data readers in this protocol test.
template <>
struct Component<mesh_plan_test::Meta> {
  using metavariables = mesh_plan_test::Meta;
  using chare_type = Parallel::Algorithms::Singleton;
  static constexpr bool checkpoint_data = true;
  using const_global_cache_tags =
      tmpl::list<mesh_plan_test::DomainTag, mesh_plan_test::PoliciesTag>;
  using simple_tags_from_options = tmpl::list<planning::Tags::State>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization, tmpl::list<>>>;
  static void execute_next_phase(
      Parallel::Phase phase,
      Parallel::CProxy_GlobalCache<metavariables>& cache);
};
}  // namespace amr
namespace mesh_plan_test {
struct Meta {
  static constexpr size_t volume_dim = 3;
  using component_list = tmpl::list<amr::Component<Meta>, Elements<Meta>>;
  struct system {
    using variables_tag = VariablesTag;
  };
  struct amr {
    using element_array = Elements<Meta>;
  };
  static constexpr const char* help = "Central AMR mesh-plan protocol test";
  static constexpr bool ignore_unrecognized_command_line_options = false;
  static constexpr std::array default_phase_order{
      Parallel::Phase::Initialization, Parallel::Phase::AdjustDomain,
      Parallel::Phase::Testing,        Parallel::Phase::Execute,
      Parallel::Phase::CheckDomain,    Parallel::Phase::Exit};
  void pup(PUP::er&) {}
};
}  // namespace mesh_plan_test
void amr::Component<mesh_plan_test::Meta>::execute_next_phase(
    const Parallel::Phase phase,
    Parallel::CProxy_GlobalCache<metavariables>& cache) {
  auto& local = *Parallel::local_branch(cache);
  auto& component = Parallel::get_parallel_component<Component>(local);
  component.start_phase(phase);
  if (phase == Parallel::Phase::AdjustDomain or
      phase == Parallel::Phase::Execute)
    Parallel::simple_action<
        planning::Actions::Begin<mesh_plan_test::CheckAndApply>>(component);
}
void register_mesh_plan_filter() {
  register_classes_with_charm<mesh_plan_test::Filter>();
}

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<mesh_plan_test::Meta>();
  Parallel::charmxx::register_init_node_and_proc({&register_mesh_plan_filter},
                                                 {});
}
