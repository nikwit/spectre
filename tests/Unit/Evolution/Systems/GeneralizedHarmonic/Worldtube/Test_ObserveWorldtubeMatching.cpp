// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/MetavariablesTag.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Events/ObserveWorldtubeMatching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace observers::Actions {
struct ContributeReductionData;
}  // namespace observers::Actions

namespace {
namespace helpers = TestHelpers::gh_worldtube;
using Observe = gh::worldtube::Events::ObserveWorldtubeMatching;
using ReductionData = gh::worldtube::Events::MatchingReductionData;

template <typename Metavariables>
struct MockContributeReductionData {
  struct Results {
    observers::ObservationId observation_id;
    std::string subfile_name;
    std::vector<std::string> reduction_names;
    ReductionData reduction_data;
    size_t contributions{0};
  };

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static std::optional<Results> results;

  template <typename ParallelComponent, typename... DbTags, typename ArrayIndex>
  static void apply(db::DataBox<tmpl::list<DbTags...>>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const observers::ObservationId& observation_id,
                    Parallel::ArrayComponentId /*sender_array_id*/,
                    const std::string& subfile_name,
                    const std::vector<std::string>& reduction_names,
                    ReductionData&& reduction_data) {
    if (results) {
      CHECK(results->observation_id == observation_id);
      CHECK(results->subfile_name == subfile_name);
      CHECK(results->reduction_names == reduction_names);
      results->reduction_data.combine(std::move(reduction_data));
      ++results->contributions;
    } else {
      results.emplace();
      *results = {observation_id, subfile_name, reduction_names,
                  std::move(reduction_data), 1};
    }
  }
};

template <typename Metavariables>
std::optional<typename MockContributeReductionData<Metavariables>::Results>
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    MockContributeReductionData<Metavariables>::results{};

template <typename Metavariables>
struct ElementComponent {
  using component_being_mocked = void;
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization, tmpl::list<>>>;
};

template <typename Metavariables>
struct MockObserverComponent {
  using component_being_mocked = observers::Observer<Metavariables>;
  using replace_these_simple_actions =
      tmpl::list<observers::Actions::ContributeReductionData>;
  using with_these_simple_actions =
      tmpl::list<MockContributeReductionData<Metavariables>>;
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockGroupChare;
  using array_index = int;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization, tmpl::list<>>>;
};

struct Metavariables {
  using system = gh::System<3>;
  using component_list = tmpl::list<ElementComponent<Metavariables>,
                                    MockObserverComponent<Metavariables>>;
  using const_global_cache_tags = tmpl::list<>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<tmpl::pair<Event, tmpl::list<Observe>>>;
  };
};

using tag_list = tmpl::list<
    Parallel::Tags::MetavariablesImpl<Metavariables>, ::Tags::Time,
    ::Tags::Variables<typename helpers::EvolvedVariables::tags_list>,
    domain::Tags::Mesh<3>,
    domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>,
    domain::Tags::Element<3>, domain::Tags::Domain<3>,
    domain::Tags::MeshVelocity<3>, gh::worldtube::Tags::KretschmannFaceData<3>>;

void test_options() {
  const auto created =
      TestHelpers::test_creation<std::unique_ptr<Event>, Metavariables>(
          "ObserveWorldtubeMatching:\n"
          "  SubfileName: WorldtubeMatching\n"
          "  Mass: 1.0");
  const auto* const observe = dynamic_cast<const Observe*>(created.get());
  REQUIRE(observe != nullptr);
  CHECK(observe->subfile_path() == "/WorldtubeMatching");
  CHECK(observe->mass() == std::optional<double>{1.0});
  CHECK(observe->needs_evolved_variables());
  const auto deserialized = serialize_and_deserialize(*observe);
  CHECK(deserialized.subfile_path() == "/WorldtubeMatching");
  CHECK(deserialized.mass() == std::optional<double>{1.0});

  const auto without_mass =
      TestHelpers::test_creation<std::unique_ptr<Event>, Metavariables>(
          "ObserveWorldtubeMatching:\n"
          "  SubfileName: WorldtubeMatching\n"
          "  Mass: None");
  CHECK_FALSE(
      dynamic_cast<const Observe*>(without_mass.get())->mass().has_value());
}

// Runs the event on the given elements of the excised Kerr-Schild domain and
// returns the finalized reduction data
template <typename... Setup>
std::optional<typename MockContributeReductionData<Metavariables>::Results>
run_event(
    const Observe& observe,
    const std::vector<std::pair<size_t, std::array<SegmentId, 3>>>& elements,
    const std::vector<bool>& expect_registered) {
  using element_component = ElementComponent<Metavariables>;
  using observer_component = MockObserverComponent<Metavariables>;
  auto& results = MockContributeReductionData<Metavariables>::results;
  results.reset();

  ActionTesting::MockRuntimeSystem<Metavariables> runner{{}};
  ActionTesting::emplace_group_component<observer_component>(&runner);

  const std::array<double, 3> velocity{{0.2, -0.1, 0.15}};
  std::vector<db::compute_databox_type<tag_list>> boxes;
  size_t expected_contributions = 0;
  for (size_t e = 0; e < elements.size(); ++e) {
    auto setup = helpers::wedge_element(12, velocity, 2.5, 4.0,
                                        elements[e].first, elements[e].second);
    // Face data of a previous evaluation, so that the event's backward
    // difference for dt K is meaningful for the moving hole
    const double previous_time = -1.e-3;
    gh::worldtube::KretschmannFaceData<3> face_data{};
    {
      const auto previous_vars = setup.evolved_variables(previous_time);
      gh::worldtube::update_kretschmann_face_data(
          make_not_null(&face_data),
          get<gr::Tags::SpacetimeMetric<DataVector, 3>>(previous_vars),
          get<gh::Tags::Pi<DataVector, 3>>(previous_vars),
          get<gh::Tags::Phi<DataVector, 3>>(previous_vars), setup.mesh,
          setup.inverse_jacobian, setup.element,
          setup.domain.excision_spheres(),
          std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>{},
          previous_time);
      face_data.filtered_moments = gr::np::TidalMoments{};
      face_data.filtered_moments_time = previous_time;
    }
    auto box = db::create<tag_list>(
        Metavariables{}, 0.0, setup.evolved_variables(0.), setup.mesh,
        setup.inverse_jacobian, setup.element, std::move(setup.domain),
        std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>{},
        std::move(face_data));
    const auto ids_to_register =
        observers::get_registration_observation_type_and_key(observe, box);
    CHECK(ids_to_register.has_value() == expect_registered[e]);
    if (ids_to_register.has_value()) {
      CHECK(ids_to_register->first == observers::TypeOfObservation::Reduction);
      CHECK(ids_to_register->second ==
            observers::ObservationKey("/WorldtubeMatching.dat"));
      ++expected_contributions;
    }
    boxes.push_back(std::move(box));
    ActionTesting::emplace_component<element_component>(&runner,
                                                        boxes.size() - 1);
  }

  for (size_t index = 0; index < boxes.size(); ++index) {
    CHECK(static_cast<const Event&>(observe).is_ready(
        boxes[index], ActionTesting::cache<element_component>(runner, index),
        static_cast<element_component::array_index>(index),
        std::add_pointer_t<element_component>{}));
    auto obs_box = make_observation_box<db::AddComputeTags<>>(
        make_not_null(&boxes[index]));
    observe.run(make_not_null(&obs_box),
                ActionTesting::cache<element_component>(runner, index),
                static_cast<element_component::array_index>(index),
                std::add_pointer_t<element_component>{}, {"Time", 0.0});
  }
  for (size_t i = 0; i < expected_contributions; ++i) {
    REQUIRE(
        not runner.template is_simple_action_queue_empty<observer_component>(
            0));
    runner.template invoke_queued_simple_action<observer_component>(0);
  }
  CHECK(runner.template is_simple_action_queue_empty<observer_component>(0));
  if (results.has_value()) {
    results->reduction_data.finalize();
  }
  return results;
}

void test_observe() {
  // Two wedges on the excision sphere and the outer half of a third
  const std::vector<std::pair<size_t, std::array<SegmentId, 3>>> elements{
      {0, {{SegmentId{0, 0}, SegmentId{0, 0}, SegmentId{0, 0}}}},
      {1, {{SegmentId{0, 0}, SegmentId{0, 0}, SegmentId{0, 0}}}},
      {2, {{SegmentId{0, 0}, SegmentId{0, 0}, SegmentId{1, 1}}}}};
  const auto results = run_event(Observe{"WorldtubeMatching", 1.0}, elements,
                                 {true, true, false});
  REQUIRE(results.has_value());
  CHECK(results->contributions == 2);
  CHECK(results->subfile_name == "/WorldtubeMatching");
  CHECK(results->reduction_names == Observe::legend());
  CHECK(results->observation_id ==
        observers::ObservationId(0.0, "/WorldtubeMatching.dat"));
  const auto& data = results->reduction_data.data();
  CHECK(std::get<0>(data) == 0.0);
  CHECK(std::get<1>(data) == 2);
  // A boosted Kerr-Schild hole is exactly type D: the type-D solve and the
  // order-zero target reproduce the face curvature to truncation error,
  // which at 12 points per dimension is far above roundoff; the order-two
  // target fits the truncation error of the face curvature as a tide and is
  // the less accurate of the two
  const double max_abs_psi0 = std::get<2>(data);
  CHECK(max_abs_psi0 > 1.e-4);
  CHECK(std::get<3>(data) > 0.);
  CHECK(std::get<4>(data) > 0.);
  // Coulomb scalar -M/r'^3 with r' between the coordinate radius 2.5 and
  // its Lorentz contraction
  const double lorentz_factor =
      1. / sqrt(1. - square(0.2) - square(0.1) - square(0.15));
  CHECK(std::get<5>(data) > -1. / cube(2.5) * 1.01);
  CHECK(std::get<6>(data) < -1. / cube(2.5 * lorentz_factor) * 0.99);
  CHECK(std::get<7>(data) < 1.e-2 / cube(2.5));
  CHECK(std::get<8>(data) < 1.e-3);
  CHECK(std::get<9>(data) == approx(max_abs_psi0).epsilon(1.e-2));
  CHECK(std::get<10>(data) < 1.e-3 * max_abs_psi0);
  CHECK(std::get<11>(data) == approx(max_abs_psi0).epsilon(1.e-2));
  CAPTURE(std::get<10>(data));
  CAPTURE(std::get<12>(data));
  CAPTURE(max_abs_psi0);
  CHECK(std::get<12>(data) < 1.e-2 * max_abs_psi0);
  CHECK(std::isfinite(std::get<13>(data)));
  CHECK(std::get<14>(data) < 1.);
  CHECK(std::isfinite(std::get<15>(data)));
  CHECK(std::isfinite(std::get<16>(data)));
  CHECK(std::get<15>(data) <= std::get<16>(data));
  CHECK(std::get<17>(data) > 2.5 * 0.99);
  CHECK(std::get<18>(data) < 2.5 * lorentz_factor * 1.01);
  // Zero relaxed moments were stored in the face data: the imposed target is
  // the type-D one
  CHECK(std::get<19>(data) == approx(std::get<9>(data)).epsilon(1.e-6));
  // Coulomb decode: exact type D, so the decoded tide vanishes to
  // truncation error, the target is the type-D one and the areal radius
  // lies in the Lorentz-contraction band
  CHECK(std::get<20>(data) == approx(max_abs_psi0).epsilon(1.e-2));
  CHECK(std::isfinite(std::get<21>(data)));
  CHECK(std::get<22>(data) > 2.5 * 0.99);
  CHECK(std::get<23>(data) < 2.5 * lorentz_factor * 1.01);

  // Without a mass the order-two columns are not evaluated
  const auto without_mass = run_event(
      Observe{"WorldtubeMatching", std::nullopt}, {elements[0]}, {true});
  REQUIRE(without_mass.has_value());
  CHECK(without_mass->contributions == 1);
  const auto& data_without_mass = without_mass->reduction_data.data();
  CHECK(std::get<1>(data_without_mass) == 1);
  CHECK(std::isnan(std::get<11>(data_without_mass)));
  CHECK(std::isnan(std::get<13>(data_without_mass)));
  CHECK(std::isnan(std::get<14>(data_without_mass)));
  CHECK(std::isnan(std::get<19>(data_without_mass)));
  CHECK(std::isnan(std::get<20>(data_without_mass)));
  CHECK(std::isnan(std::get<23>(data_without_mass)));
  CHECK(std::isnan(std::get<18>(data_without_mass)));

  // An element off the excision sphere contributes nothing
  const auto none =
      run_event(Observe{"WorldtubeMatching", 1.0}, {elements[2]}, {false});
  CHECK_FALSE(none.has_value());
}

// A hole at rest on the spherical-harmonic excision face of a shell element:
// in Kerr-Schild coordinates the face curvature is exactly type D with a
// vanishing entering mode, and the invariant boost is the one from the
// Kerr-Schild normal observers to the static observers, tanh(eta) =
// beta^r sqrt(gamma_rr) / alpha = 2M/r.
void test_static_hole_on_spherical_shell() {
  const double radius = 2.5;
  const auto setup =
      helpers::shell_element(12, 10, {{0., 0., 0.}}, radius, 5.0);
  REQUIRE(setup.domain.excision_spheres().size() == 1);
  const auto vars = setup.evolved_variables(0.);
  const Observe observe{"WorldtubeMatching", 1.0};
  const auto reduction = observe.compute_reduction_data(
      0., get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
      get<gh::Tags::Pi<DataVector, 3>>(vars),
      get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
      setup.inverse_jacobian, setup.element, setup.domain,
      std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>{},
      gh::worldtube::KretschmannFaceData<3>{});
  REQUIRE(reduction.has_value());
  const auto data = reduction->data();
  const double coulomb = -1. / cube(radius);
  CAPTURE(data);
  CHECK(std::get<1>(data) == 1);
  // |E| = sqrt(6) M / r^3, B = 0, no radiation, r measured exactly
  CHECK(std::get<2>(data) < 1.e-8 * std::abs(coulomb));
  CHECK(std::get<3>(data) == approx(sqrt(6.) / cube(radius)).epsilon(1.e-6));
  CHECK(std::get<4>(data) < 1.e-8 * std::abs(coulomb));
  CHECK(std::get<5>(data) == approx(coulomb).epsilon(1.e-6));
  CHECK(std::get<6>(data) == approx(coulomb).epsilon(1.e-6));
  CHECK(std::get<7>(data) < 1.e-8 * std::abs(coulomb));
  CHECK(std::get<8>(data) < 1.e-6);
  CHECK(std::get<9>(data) < 1.e-8 * std::abs(coulomb));
  CHECK(std::get<10>(data) < 1.e-8 * std::abs(coulomb));
  CHECK(std::get<14>(data) == approx(2. / radius).epsilon(1.e-3));
  CHECK(std::get<15>(data) == approx(-std::atanh(2. / radius)).epsilon(1.e-3));
  CHECK(std::get<16>(data) == approx(-std::atanh(2. / radius)).epsilon(1.e-3));
  CHECK(std::get<17>(data) == approx(radius).epsilon(1.e-6));
  CHECK(std::get<18>(data) == approx(radius).epsilon(1.e-6));
  CHECK(std::isnan(std::get<19>(data)));
  CHECK(std::get<20>(data) < 1.e-6 * std::abs(coulomb));
  CHECK(std::get<22>(data) == approx(radius).epsilon(1.e-4));
  CHECK(std::get<23>(data) == approx(radius).epsilon(1.e-4));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.ObserveMatching",
    "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  test_options();
  test_observe();
  test_static_hole_on_spherical_shell();
}
