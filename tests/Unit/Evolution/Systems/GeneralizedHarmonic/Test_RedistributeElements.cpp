// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/MinimumGridSpacing.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Bbh/PhaseControl/RedistributeElements.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "ParallelAlgorithms/Amr/RedistributeElements.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Controller = gh::bbh::phase_control::RedistributeElements;
struct Metavariables {
  using component_list = tmpl::list<>;
  struct amr {
    using element_array = int;
  };
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<PhaseChange, tmpl::list<Controller>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.RedistributeElements.Phases",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  using namespace gh::bbh::phase_control::redistribution;
  using Parallel::Phase;
  using Data = tuples::TaggedTuple<StateTag, Requested>;
  const Parallel::GlobalCache<Metavariables> cache{};
  auto controller = TestHelpers::test_creation<Controller>("");
  Data data{};
  controller.initialize_phase_data_impl(make_not_null(&data));
  const auto jump = [&](const Phase phase) -> std::optional<Phase> {
    const auto result = controller.arbitrate_phase_change_impl(
        make_not_null(&data), phase, cache);
    return result ? std::optional{result->first} : std::nullopt;
  };
  CHECK(jump(Phase::Initialization) == std::nullopt);
  CHECK(jump(Phase::InitializeInitialDataDependentQuantities) == std::nullopt);
  CHECK(jump(Phase::AdjustDomain) == std::nullopt);
  CHECK(jump(Phase::Register) == std::nullopt);
  CHECK(jump(Phase::Restart) == std::nullopt);
  CHECK(jump(Phase::Evolve) == std::nullopt);
  for (size_t i = 0; i < 2; ++i) {
    tuples::get<Requested>(data) = true;
    CHECK(jump(Phase::Evolve) == Phase::CollectElementDistribution);
    CHECK_FALSE(tuples::get<Requested>(data));
    CHECK(jump(Phase::CollectElementDistribution) ==
          Phase::PlanElementDistribution);
    // Simulate checkpointed state during a redistribution.
    data = serialize_and_deserialize(data);
    controller = serialize_and_deserialize(controller);
    CHECK(jump(Phase::PlanElementDistribution) == Phase::RedistributeElements);
    CHECK(jump(Phase::RedistributeElements) == Phase::UpdateSections);
    data = serialize_and_deserialize(data);
    controller = serialize_and_deserialize(controller);
    const auto result = controller.arbitrate_phase_change_impl(
        make_not_null(&data), Phase::UpdateSections, cache);
    REQUIRE(result.has_value());
    CHECK(result->first == Phase::Evolve);
    CHECK(result->second == PhaseControl::PermitAdditionalJumps);
  }
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.RedistributeElements.Plan",
                  "[Unit][Evolution]") {
  std::vector<Block<1>> blocks{};
  blocks.emplace_back(
      domain::make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
          domain::CoordinateMaps::Identity<1>{}),
      0, DirectionMap<1, BlockNeighbors<1>>{});
  using namespace amr::redistribution;
  amr::redistribution::Tags::Costs<1>::type records{};
  const auto ids = initial_element_ids(0, std::array<size_t, 1>{{3}});
  for (size_t i = 0; i < ids.size(); ++i) {
    // Same initial p mesh would give identical costs; current p is unequal.
    records.emplace(ids[i], std::pair{i < 4 ? 10. : 1., 0});
  }
  const auto plan = make_plan(blocks, records, 3, {1}, true);
  CHECK(plan.size() == ids.size());
  double costs[3]{0., 0., 0.};
  for (const auto& [id, target] : plan) {
    CHECK(target != 1);
    costs[target] += records.at(id).first;
  }
  CHECK(std::max(costs[0], costs[2]) < 30.);
  CHECK(costs[1] == 0.);
  CHECK(plan ==
        make_plan(blocks, serialize_and_deserialize(records), 3, {1}, true));
  const auto round_robin = make_plan(blocks, records, 3, {1}, false);
  for (size_t i = 0; i < ids.size(); ++i) {
    CHECK(round_robin.at(ids[i]) == (i % 2 == 0 ? 0 : 2));
  }
  auto invalid = records;
  invalid.erase(ids.back());
  CHECK_THROWS(make_plan(blocks, invalid, 3, {1}, true));
  invalid = records;
  invalid.emplace(ElementId<1>{0}, std::pair{1., 0});
  CHECK_THROWS(make_plan(blocks, invalid, 3, {}, true));
  invalid = records;
  invalid.at(ids[0]).first = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS(make_plan(blocks, invalid, 3, {}, true));
  CHECK_THROWS(make_plan(blocks, records, 3, {0, 1, 2}, true));
  CHECK_THROWS(make_plan(blocks, records, 3, {3}, true));

  // The spacing model must use the updated mesh, not the input-file extents.
  double previous = 0.;
  for (const size_t n : {3_st, 8_st}) {
    const Mesh<1> mesh{n, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
    const auto coords = ElementMap<1, Frame::Grid>{
        ids[0], blocks[0]}(logical_coordinates(mesh));
    CHECK(domain::element_cost(mesh, coords, domain::ElementWeight::Uniform) ==
          1.);
    CHECK(domain::element_cost(mesh, coords,
                               domain::ElementWeight::NumGridPoints) == n);
    const double cost = domain::element_cost(
        mesh, coords, domain::ElementWeight::NumGridPointsAndGridSpacing);
    CHECK(cost ==
          approx(n / sqrt(minimum_grid_spacing(mesh.extents(), coords))));
    CHECK(cost > previous);
    previous = cost;
  }
}
