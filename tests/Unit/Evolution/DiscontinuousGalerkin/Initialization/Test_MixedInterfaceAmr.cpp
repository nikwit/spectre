// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Info.hpp"
#include "Domain/Amr/NeighborsOfChild.hpp"
#include "Domain/Amr/UpdateAmrDecision.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/NonconformingSphericalShells.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/CreateInitialMesh.hpp"
#include "Domain/Structure/FaceType.hpp"
#include "Evolution/DiscontinuousGalerkin/EqualRateLts/EqualRateRegions.hpp"
#include "Evolution/DiscontinuousGalerkin/EqualRateLts/EqualRateRegions.tpp"
#include "Evolution/DiscontinuousGalerkin/EqualRateLts/NonconformingEqualRateRegions.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/SetupEqualRateRegions.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "Time/History.hpp"
#include "Time/TimeSteppers/AdamsMoultonPc.hpp"
#include "Utilities/MakeArray.hpp"

namespace {
using Shells = domain::creators::NonconformingSphericalShells;
using Flag = amr::Flag;

std::unique_ptr<DomainCreator<3>> make_creator(const size_t angular_h) {
  return std::make_unique<Shells>(
      2., 3., 4., Shells::Excision{}, std::array{angular_h, size_t{0}},
      size_t{1}, std::array<size_t, 2>{9, 5}, std::array<size_t, 2>{6, 5});
}

Mesh<3> change_mesh(const Mesh<3>& mesh, const int radial, const int angular) {
  auto extents = mesh.extents();
  const bool shell = mesh.basis(1) == Spectral::Basis::SphericalHarmonic;
  for (size_t d = 0; d < 3; ++d) {
    const int change =
        shell ? (d == 0 ? radial : (d == 1 ? angular : 2 * angular))
              : (d == 2 ? radial : angular);
    extents[d] = static_cast<size_t>(static_cast<int>(extents[d]) + change);
  }
  return {extents.indices(), mesh.basis(), mesh.quadrature()};
}

// In particular, the S2 side has 96 wedge neighbors at angular_h=2, more
// than can fit in NeighborMesh. The two sides also have unrelated h levels.
void test_neighbor_bookkeeping() {
  const auto creator = make_creator(2);
  const auto domain = creator->create_domain();
  const auto levels = creator->initial_refinement_levels();
  const auto extents = creator->initial_extents();
  for (const auto& id : {ElementId<3>{6, {{{1, 0}, {0, 0}, {0, 0}}}},
                         ElementId<3>{2, {{{2, 0}, {2, 1}, {0, 0}}}}}) {
    CAPTURE(id);
    const auto element =
        domain::create_initial_element(id, domain.blocks(), levels);
    const auto mesh =
        domain::create_initial_mesh(extents, element, Spectral::Basis::Legendre,
                                    Spectral::Quadrature::GaussLobatto);
    const amr::Info<3> info{make_array<3>(Flag::IncreaseResolution),
                            change_mesh(mesh, 1, 1)};
    std::unordered_map<ElementId<3>, amr::Info<3>> neighbor_info{};
    for (const auto& [direction, neighbors] : element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        const auto neighbor_element =
            domain::create_initial_element(neighbor, domain.blocks(), levels);
        neighbor_info.emplace(
            neighbor,
            amr::Info<3>{make_array<3>(Flag::DecreaseResolution),
                         change_mesh(domain::create_initial_mesh(
                                         extents, neighbor_element,
                                         Spectral::Basis::Legendre,
                                         Spectral::Quadrature::GaussLobatto),
                                     -1, -1)});
        auto flags = info.flags;
        CHECK_FALSE(
            amr::update_amr_decision(make_not_null(&flags), element, neighbor,
                                     neighbor_info.at(neighbor).flags, true));
        CHECK(flags == info.flags);
      }
    }
    const auto [new_neighbors, meshes] =
        amr::neighbors_of_child(element, info, neighbor_info, id);
    CHECK(new_neighbors == element.neighbors());
    size_t expected_meshes = 0;
    bool checked_h_restriction = false;
    for (const auto& [direction, neighbors] : element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        const DirectionalId<3> key{direction, neighbor};
        const auto& raw_mesh = neighbor_info.at(neighbor).new_mesh;
        if (element.face_types().at(direction) ==
            domain::FaceType::MultipleNonconforming) {
          CHECK_FALSE(meshes.contains(key));
        } else {
          ++expected_meshes;
          CHECK(meshes.at(key) ==
                (neighbors.are_conforming()
                     ? neighbors.orientation(neighbor).inverse_map()(raw_mesh)
                     : raw_mesh));
        }
        if (not neighbors.are_conforming() and not checked_h_restriction) {
          checked_h_restriction = true;
          auto flags = info.flags;
          flags[0] = Flag::Split;
          CHECK_THROWS_WITH(
              amr::update_amr_decision(make_not_null(&flags), element, neighbor,
                                       neighbor_info.at(neighbor).flags, true),
              Catch::Matchers::ContainsSubstring("p refinement only"));
          auto h_info = info;
          h_info.flags = flags;
          CHECK_THROWS_WITH(
              amr::neighbors_of_child(element, h_info, neighbor_info, id),
              Catch::Matchers::ContainsSubstring("p refinement only"));
        }
      }
    }
    CHECK(meshes.size() == expected_meshes);
    CHECK(checked_h_restriction);
  }
}

// Low angular modes and a radial polynomial, exactly representable before
// and after each projection, with nonconstant geometric and evolved data.
template <size_t Dim>
DataVector profile(const Mesh<Dim>& mesh) {
  const auto x = logical_coordinates(mesh);
  if constexpr (Dim == 2) {
    if (mesh.basis(0) == Spectral::Basis::SphericalHarmonic) {
      return 2. + cos(get<0>(x)) + sin(get<0>(x)) * cos(get<1>(x));
    }
    return 2. + get<0>(x) + square(get<1>(x));
  } else {
    if (mesh.basis(1) == Spectral::Basis::SphericalHarmonic) {
      return 2. + square(get<0>(x)) + cos(get<1>(x)) +
             sin(get<1>(x)) * cos(get<2>(x));
    }
    return 2. + get<0>(x) + square(get<1>(x)) + get<2>(x);
  }
}

struct Field : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct Metavariables {
  static constexpr size_t volume_dim = 3;
  struct system {
    using variables_tag = ::Tags::Variables<tmpl::list<Field>>;
  };
};
using Vars = Metavariables::system::variables_tag::type;
using DtVars = TimeSteppers::History<Vars>::DerivVars;

Vars fields(const Mesh<3>& mesh, const double time) {
  Vars result(mesh.number_of_grid_points());
  get(get<Field>(result)) = (1. + time) * profile(mesh);
  return result;
}

DtVars derivatives(const Mesh<3>& mesh, const double time) {
  return DtVars{fields(mesh, 2. * time)};
}

void check_volume_history(const TimeSteppers::History<Vars>& history,
                          const Mesh<3>& mesh) {
  CHECK(history.integration_order() == 4);
  CHECK(history.size() == 3);
  CHECK(history.substeps().size() == 1);
  const auto check_entry = [&mesh](const auto& entry) {
    const auto t = entry.time_step_id.substep_time();
    REQUIRE(entry.value.has_value());
    CHECK_VARIABLES_APPROX(*entry.value, fields(mesh, t));
    CHECK_VARIABLES_APPROX(entry.derivative, derivatives(mesh, t));
  };
  for (const auto& entry : history) {
    check_entry(entry);
  }
  check_entry(history.substeps()[0]);
}

void test_projectors(const bool lts, const bool shell_host) {
  CAPTURE(lts, shell_host);
  const auto creator = make_creator(0);
  const auto domain = creator->create_domain();
  const auto levels = creator->initial_refinement_levels();
  const auto extents = creator->initial_extents();
  const ElementId<3> id = shell_host
                              ? ElementId<3>{6, {{{1, 0}, {0, 0}, {0, 0}}}}
                              : ElementId<3>{2};
  const auto element =
      domain::create_initial_element(id, domain.blocks(), levels);
  auto mesh =
      domain::create_initial_mesh(extents, element, Spectral::Basis::Legendre,
                                  Spectral::Quadrature::GaussLobatto);
  const auto neighbor_meshes = [&](const int radial, const int angular) {
    DirectionalIdMap<3, Mesh<3>> result{};
    for (const auto& [direction, neighbors] : element.neighbors()) {
      if (element.face_types().at(direction) ==
          domain::FaceType::MultipleNonconforming) {
        continue;
      }
      for (const auto& neighbor : neighbors) {
        const auto neighbor_element =
            domain::create_initial_element(neighbor, domain.blocks(), levels);
        auto neighbor_mesh = change_mesh(
            domain::create_initial_mesh(extents, neighbor_element,
                                        Spectral::Basis::Legendre,
                                        Spectral::Quadrature::GaussLobatto),
            radial, angular);
        if (neighbors.are_conforming()) {
          neighbor_mesh =
              neighbors.orientation(neighbor).inverse_map()(neighbor_mesh);
        }
        result.emplace(DirectionalId<3>{direction, neighbor}, neighbor_mesh);
      }
    }
    return result;
  };
  const auto mode = lts ? LtsMode::Conservative : LtsMode::Off;
  const Slab slab{0., 1.};
  auto time_step = slab.duration() / 8;
  const TimeStepId current{true, 0, slab.start() + slab.duration() / 2};
  auto next = TimeStepId{true, 0, current.step_time() + time_step};
  const auto old_next = next;
  const auto old_step = time_step;
  const evolution::dg::EqualRateRegions<
      3, tmpl::list<evolution::dg::NonconformingEqualRateRegions<3>>>
      regions{creator};
  auto data = evolution::dg::Initialization::detail::empty_mortar_data(element);
  auto infos = evolution::dg::Initialization::detail::mortar_infos(
      domain, element, mesh, neighbor_meshes(0, 0), mode);
  auto [mortars, clocks, normals] =
      evolution::dg::Initialization::detail::mortars_apply_impl(
          element, next, mesh, neighbor_meshes(0, 0));
  const auto set_regions = [&](auto* target) {
    if (lts) {
      std::optional<size_t> fixed_ratio{};
      evolution::dg::Initialization::SetupLocalEqualRateRegion<3>::apply(
          make_not_null(&fixed_ratio), make_not_null(target), element, regions,
          time_step);
      CHECK(fixed_ratio == 8);
    }
  };
  set_regions(&infos);
  const auto initial_infos = infos;
  const auto initial_clocks = clocks;
  evolution::dg::Tags::MortarDataHistory<3>::type histories{};
  // Equal-rate faces deliberately have no BoundaryHistory, while the outer
  // radial S2 face has a populated conservative LTS history.
  for (const auto& [key, info] : infos) {
    if (info.time_stepping_policy() !=
        evolution::dg::TimeSteppingPolicy::Conservative) {
      continue;
    }
    auto& h = histories[key];
    for (int i = 0; i < 3; ++i) {
      const TimeStepId past{true, 0, slab.start() + time_step * i};
      const auto face = mesh.on_interface(key.direction().dimension());
      evolution::dg::MortarData<3> local{};
      local.mortar_mesh = mortars.at(key);
      local.mortar_data = profile(*local.mortar_mesh);
      local.face_mesh = face;
      local.face_normal_magnitude = Scalar<DataVector>{profile(face)};
      local.face_det_jacobian = Scalar<DataVector>{2. * profile(face)};
      local.volume_mesh = mesh;
      local.volume_det_inv_jacobian = Scalar<DataVector>{profile(mesh)};
      h.local().insert(past, 4, local);
      h.remote().insert(past, 4, local);
    }
  }
  CHECK(histories.size() == (lts and shell_host ? 1 : 0));
  const auto initial_histories = histories;
  TimeSteppers::History<Vars> history{4};
  for (int i = 0; i < 3; ++i) {
    const TimeStepId past{true, 0, slab.start() + time_step * i};
    const double t = past.substep_time();
    history.insert(past, fields(mesh, t), derivatives(mesh, t));
  }
  const TimeStepId substep{
      true, 0, history.back().time_step_id.step_time(), 1, time_step, 0.375};
  history.insert(substep, fields(mesh, 0.375), derivatives(mesh, 0.375));
  const auto initial_history = history;
  DtVars dt_vars(mesh.number_of_grid_points(), 0.);
  const TimeSteppers::AdamsMoultonPc<true> time_stepper{4};
  // Increase, decrease, radial-only, angular-only, neighbor-only, no change.
  for (const auto& resolution :
       std::array<std::array<int, 4>, 6>{{{1, 1, 1, 1},
                                          {-1, -1, -1, -1},
                                          {1, 0, 0, 0},
                                          {0, 1, 0, 1},
                                          {0, 0, 1, 1},
                                          {0, 0, 1, 1}}}) {
    CAPTURE(resolution);
    const auto new_mesh = change_mesh(mesh, resolution[0], resolution[1]);
    const auto neighbors = neighbor_meshes(resolution[2], resolution[3]);
    for (auto& [direction, normal] : normals) {
      normal.emplace(
          mesh.on_interface(direction.dimension()).number_of_grid_points(), 1.);
    }
    evolution::dg::Initialization::ProjectMortars<3>::apply(
        make_not_null(&data), make_not_null(&mortars), make_not_null(&infos),
        make_not_null(&clocks), make_not_null(&normals),
        make_not_null(&histories), domain, new_mesh, element, neighbors,
        current, mode, std::pair{mesh, element});
    auto expected_infos = evolution::dg::Initialization::detail::mortar_infos(
        domain, element, new_mesh, neighbors, mode);
    set_regions(&expected_infos);
    CHECK(infos == expected_infos);
    CHECK(clocks == initial_clocks);
    CHECK(data ==
          evolution::dg::Initialization::detail::empty_mortar_data(element));
    for (const auto& [direction, normal] : normals) {
      CHECK(normal.has_value() ==
            (mesh.on_interface(direction.dimension()) ==
             new_mesh.on_interface(direction.dimension())));
    }
    for (const auto& [key, info] : infos) {
      CHECK(info.time_stepping_policy() ==
            initial_infos.at(key).time_stepping_policy());
      if (not element.neighbors().at(key.direction()).are_conforming()) {
        CHECK(mortars.at(key) ==
              new_mesh.on_interface(key.direction().dimension()));
      }
      if (info.interpolator().has_value()) {
        const auto& interpolator = *info.interpolator();
        const auto face = new_mesh.on_interface(key.direction().dimension());
        const DataVector source(
            interpolator.neighbor_mortar_mesh().number_of_grid_points(), 2.5);
        CHECK_ITERABLE_APPROX(interpolator.interpolate_to_host(source),
                              DataVector(face.number_of_grid_points(), 2.5));
        CHECK_ITERABLE_APPROX(
            interpolator.interpolate_to_neighbor(
                DataVector(face.number_of_grid_points(), 3.5)),
            DataVector(interpolator.interpolated_neighbor_data_offsets().size(),
                       3.5));
        test_serialization(interpolator);
      }
    }
    CHECK(histories.size() == initial_histories.size());
    for (const auto& [key, h] : histories) {
      const auto face = new_mesh.on_interface(key.direction().dimension());
      CHECK(h.local().size() == initial_histories.at(key).local().size());
      CHECK(h.remote().size() == initial_histories.at(key).remote().size());
      h.local().for_each([&](const TimeStepId& t, const auto& entry) {
        CHECK(h.local().integration_order(t) == 4);
        CHECK(entry.face_mesh == face);
        CHECK(entry.volume_mesh == new_mesh);
        CHECK_ITERABLE_APPROX(get(*entry.face_normal_magnitude), profile(face));
        CHECK_ITERABLE_APPROX(get(*entry.face_det_jacobian),
                              2. * profile(face));
        CHECK_ITERABLE_APPROX(get(*entry.volume_det_inv_jacobian),
                              profile(new_mesh));
        const auto& old = initial_histories.at(key).local().data(t);
        CHECK(entry.mortar_mesh == old.mortar_mesh);
        CHECK(entry.mortar_data == old.mortar_data);
      });
      h.remote().for_each([&](const TimeStepId& t, const auto& entry) {
        CHECK(entry == initial_histories.at(key).remote().data(t));
        CHECK(h.remote().integration_order(t) == 4);
      });
    }
    ::Initialization::ProjectTimeStepperHistory<Metavariables>::apply(
        make_not_null(&dt_vars), make_not_null(&history), make_not_null(&next),
        make_not_null(&time_step), new_mesh, id, time_stepper, infos,
        std::pair{mesh, element});
    check_volume_history(history, new_mesh);
    CHECK(next == old_next);
    CHECK(time_step == old_step);
    CHECK(dt_vars.number_of_grid_points() == new_mesh.number_of_grid_points());
    for (size_t i = 0; i < history.size(); ++i) {
      CHECK(history[i].time_step_id == initial_history[i].time_step_id);
    }
    CHECK(history.substeps()[0].time_step_id ==
          initial_history.substeps()[0].time_step_id);
    // Exercise restart/redistribution serialization between successive AMR
    // calls.
    infos = serialize_and_deserialize(infos);
    history = serialize_and_deserialize(history);
    histories = serialize_and_deserialize(histories);
    mesh = new_mesh;
  }
}
}  // namespace

// [[TimeOut, 60]]
SPECTRE_TEST_CASE("Unit.Evolution.DG.Initialization.MixedInterfaceAmr",
                  "[Unit][Evolution]") {
  domain::creators::register_derived_with_charm();
  test_neighbor_bookkeeping();
  for (const bool lts : {false, true}) {
    test_projectors(lts, false);
    test_projectors(lts, true);
  }
}
