// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/Variables.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/RotScaleTrans.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeRecordingData.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "IO/Observer/Initialize.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveMovingSphere.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/FileSystem.hpp"

namespace {
struct OutputMetavariables;
struct Target {
  using vars_to_interpolate_to_target = gh::worldtube::recording_tags;
  using compute_target_points =
      intrp::TargetPoints::Sphere<Target, Frame::Inertial>;
  template <typename Meta>
  static std::array<double, 3> inertial_center(
      const std::array<double, 3>& grid_center,
      const Parallel::GlobalCache<Meta>& cache, const double time) {
    if constexpr (std::is_same_v<Meta, OutputMetavariables>) {
      return intrp::sphere_detail::excision_center(
          grid_center, Parallel::get<domain::Tags::Domain<3>>(cache), time, {});
    } else {
      return intrp::sphere_detail::excision_center(
          grid_center, Parallel::get<domain::Tags::Domain<3>>(cache), time,
          Parallel::get<domain::Tags::FunctionsOfTime>(cache));
    }
  }
};
using SphereOptions = intrp::OptionHolders::Sphere;
using Points = intrp::TargetPoints::Sphere<Target, Frame::Inertial>;
struct Metavariables {
  using component_list = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::append<Points::const_global_cache_tags,
                   tmpl::list<domain::Tags::Domain<3>>>;
  using mutable_global_cache_tags = tmpl::list<domain::Tags::FunctionsOfTime>;
};

template <typename Metavariables>
struct MockObserverWriter {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockNodeGroupChare;
  using array_index = size_t;
  using const_global_cache_tags = tmpl::list<observers::Tags::ReductionFileName,
                                             observers::Tags::SurfaceFileName>;
  using simple_tags =
      typename observers::Actions::InitializeWriter<Metavariables>::simple_tags;
  using compute_tags = typename observers::Actions::InitializeWriter<
      Metavariables>::compute_tags;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<observers::Actions::InitializeWriter<Metavariables>>>>;
  using component_being_mocked = observers::ObserverWriter<Metavariables>;
};

struct OutputMetavariables {
  using observed_reduction_data_tags = tmpl::list<>;
  static constexpr size_t volume_dim = 3;
  using component_list = tmpl::list<MockObserverWriter<OutputMetavariables>>;
  using const_global_cache_tags =
      tmpl::append<Points::const_global_cache_tags,
                   tmpl::list<domain::Tags::Domain<3>>>;
};

void test_callback() {
  const std::string prefix = "TestMovingSphereCallback";
  if (file_system::check_if_file_exists(prefix + ".h5")) {
    file_system::rm(prefix + ".h5", true);
  }
  using Meta = OutputMetavariables;
  using Writer = MockObserverWriter<Meta>;
  const SphereOptions options{
      3, {{2., 3., -1.}}, 1.2, ylm::AngularOrdering::Strahlkorper};
  Parallel::GlobalCache<Meta>::ConstTagsTuple init{};
  tuples::get<intrp::Tags::Sphere<Target>>(init) = options;
  tuples::get<domain::Tags::Domain<3>>(init) = Domain<3>{
      std::vector<Block<3>>{},
      {{"ExcisionSphereB",
        ExcisionSphere<3>{
            0.5, tnsr::I<double, 3, Frame::Grid>{{2., 3., -1.}}, {}}}}};
  tuples::get<observers::Tags::SurfaceFileName>(init) = prefix;
  tuples::get<observers::Tags::ReductionFileName>(init) = prefix;
  ActionTesting::MockRuntimeSystem<Meta> runner{std::move(init)};
  ActionTesting::emplace_nodegroup_component<Writer>(make_not_null(&runner));
  ActionTesting::set_phase(make_not_null(&runner),
                           Parallel::Phase::Initialization);
  ActionTesting::next_action<Writer>(make_not_null(&runner), 0);
  auto box =
      db::create<tmpl::append<db::AddSimpleTags<::Tags::Variables<
                                  Target::vars_to_interpolate_to_target>>,
                              Points::simple_tags>,
                 db::AddComputeTags<Points::compute_tags>>(
          Variables<Target::vars_to_interpolate_to_target>(28, 1.25),
          ylm::Strahlkorper<Frame::Inertial>{},
          tnsr::I<DataVector, 3, Frame::Inertial>{});
  Points::initialize(make_not_null(&box),
                     ActionTesting::cache<Writer>(runner, size_t{0}));
  intrp::callbacks::ObserveMovingSphere<Target>::apply(
      box, ActionTesting::cache<Writer>(runner, size_t{0}), 2.5);
  REQUIRE(ActionTesting::number_of_queued_threaded_actions<Writer>(runner, 0) ==
          1);
  ActionTesting::invoke_queued_threaded_action<Writer>(make_not_null(&runner),
                                                       0);
  {
    h5::H5File<h5::AccessType::ReadOnly> file(prefix + ".h5");
    const auto& volume = file.get<h5::VolumeData>("Target");
    const auto ids = volume.list_observation_ids();
    REQUIRE(ids.size() == 1);
    CHECK(volume.get_observation_value(ids[0]) == 2.5);
    CHECK(std::get<DataVector>(
              volume.get_tensor_component(ids[0], "Pi_tt").data) ==
          DataVector(28, 1.25));
  }
  file_system::rm(prefix + ".h5", true);
}

void test_geometry_and_output(const ylm::AngularOrdering ordering) {
  auto options = TestHelpers::test_creation<SphereOptions>(
      "LMax: 4\nCenter: [2, 3, -1]\nRadius: [2.3, 1.2]\n"
      "AngularOrdering: Strahlkorper");
  options.angular_ordering = ordering;
  CHECK(serialize_and_deserialize(options) == options);

  const tnsr::I<double, 3, Frame::Grid> grid_center{{2., 3., -1.}};
  Domain<3> domain{
      std::vector<Block<3>>{},
      {{"ExcisionSphereB", ExcisionSphere<3>{0.5, grid_center, {}}}}};
  CHECK_THROWS_WITH(
      intrp::sphere_detail::excision_center({{0., 0., 0.}}, domain, 0., {}),
      Catch::Matchers::ContainsSubstring(
          "must equal a domain excision center"));
  CHECK(intrp::sphere_detail::excision_center(options.center, domain, 0., {}) ==
        std::array<double, 3>{{2., 3., -1.}});
  using Map = domain::CoordinateMaps::TimeDependent::RotScaleTrans<3>;
  domain.inject_time_dependent_map_for_excision_sphere(
      "ExcisionSphereB",
      domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(Map{
          std::pair<std::string, std::string>{"Expansion", "OuterExpansion"},
          std::nullopt, "Translation", 100., 200., Map::BlockRegion::Inner}));
  using Polynomial = domain::FunctionsOfTime::PiecewisePolynomial<2>;
  domain::FunctionsOfTimeMap fots;
  fots["Expansion"] = std::make_unique<Polynomial>(
      0., std::array<DataVector, 3>{{{2.}, {0.2}, {0.}}}, 10.);
  fots["OuterExpansion"] = std::make_unique<Polynomial>(
      0., std::array<DataVector, 3>{{{1.}, {0.}, {0.}}}, 10.);
  fots["Translation"] = std::make_unique<Polynomial>(
      0.,
      std::array<DataVector, 3>{
          {{1., -2., 3.}, {0.1, 0.2, -0.3}, {0., 0., 0.}}},
      10.);
  Parallel::GlobalCache<Metavariables> cache{{options, std::move(domain)},
                                             {std::move(fots)}};
  auto box = db::create<db::AddSimpleTags<Points::simple_tags>,
                        db::AddComputeTags<Points::compute_tags>>();
  Points::initialize(make_not_null(&box), cache);
  const auto offsets = Points::points(box, tmpl::type_<Metavariables>{});
  const size_t n = (options.l_max + 1) * (2 * options.l_max + 1);
  CHECK(get<0>(offsets).size() == 2 * n);

  using Fields = gh::worldtube::recording_tags;
  Variables<Fields> vars(2 * n);
  for (size_t i = 0; i < vars.size(); ++i) {
    vars.data()[i] = 0.125 + static_cast<double>(i);
  }
  const std::string filename = "TestMovingSphere.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  {
    h5::H5File<h5::AccessType::ReadWrite> file(filename);
    auto& volume = file.insert<h5::VolumeData>("GhWorldtube", 0);
    size_t observation = 0;
    // Include out-of-order observations: the cached offsets must never change.
    for (const double time : {0., 1.5, 0.5}) {
      const auto center = Points::center(cache, time);
      CHECK_ITERABLE_APPROX(
          center,
          (std::array<double, 3>{{(2. + 0.2 * time) * 2. + 1. + 0.1 * time,
                                  (2. + 0.2 * time) * 3. - 2. + 0.2 * time,
                                  -(2. + 0.2 * time) + 3. - 0.3 * time}}));
      const auto points = Points::points_at_time(offsets, cache, time);
      size_t shell = 0;
      for (const double radius : options.radii) {
        for (size_t p = 0; p < n; ++p) {
          double r_squared = 0.;
          for (size_t i = 0; i < 3; ++i) {
            r_squared += square(points->get(i)[shell * n + p] - center[i]);
          }
          CHECK(sqrt(r_squared) == approx(radius));
        }
        ++shell;
      }
      CHECK(Points::points(box, tmpl::type_<Metavariables>{}) == offsets);
      const auto data =
          intrp::callbacks::observe_moving_sphere_detail::volume_data(
              options, center, *points, vars);
      REQUIRE(data.size() == 2);
      CHECK(data[0].tensor_components.size() == 211);
      volume.write_volume_data(observation, time, data);
      CHECK(volume.get_observation_value(observation) == time);
      CHECK(volume.get_extents(observation) ==
            std::vector<std::vector<size_t>>(
                2, ordering == ylm::AngularOrdering::Strahlkorper
                       ? std::vector<size_t>{5, 9}
                       : std::vector<size_t>{9, 5}));
      // Read every component back from HDF5, including the GH fields,
      // derivatives, K and dK.
      for (size_t c = 0; c < data[0].tensor_components.size(); ++c) {
        const auto written = volume.get_tensor_component(
            observation, data[0].tensor_components[c].name);
        const auto& values = std::get<DataVector>(written.data);
        REQUIRE(values.size() == 2 * n);
        for (size_t s = 0; s < 2; ++s) {
          const auto& expected =
              std::get<DataVector>(data[s].tensor_components[c].data);
          for (size_t p = 0; p < n; ++p) {
            CHECK(values[s * n + p] == expected[p]);
          }
        }
      }
      ++observation;
    }
  }
  file_system::rm(filename, true);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.NumericalAlgorithms.InterpolationTarget.MovingSphere",
                  "[Unit]") {
  test_geometry_and_output(ylm::AngularOrdering::Strahlkorper);
  test_geometry_and_output(ylm::AngularOrdering::Cce);
  test_callback();
}
