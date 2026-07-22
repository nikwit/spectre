// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BatchedMapPlan.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/DriverState.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/InitializeBatchedData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateConstraintDampingBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/ConstraintDamping/Constant.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace {

constexpr size_t Dim = 3;
using system = gh::System<Dim>;
using hardcoded_solution = gh::Solutions::WrappedGr<gr::Solutions::KerrSchild>;
using packed_topology_type = evolution::Kokkos::PackedTopology<system>;
using packed_geometry_type = evolution::Kokkos::PackedGeometry<system>;
using packed_evolution_state_type =
    evolution::Kokkos::PackedEvolutionState<system>;
using device_gamma0_type = gh::KokkosTags::DeviceConstraintGamma0::type;
using device_gamma1_type = gh::KokkosTags::DeviceConstraintGamma1::type;
using device_gamma2_type = gh::KokkosTags::DeviceConstraintGamma2::type;

struct DomainStaticCacheSnapshot {
  std::array<DataVector, Dim> grid_coordinates{};
  std::vector<double> logical_to_grid_inverse_jacobian{};
  size_t num_elements = 0;
  size_t points_per_element = 0;
};

DomainStaticCacheSnapshot snapshot_domain_static_caches(
    const packed_geometry_type& packed_geometry) {
  DomainStaticCacheSnapshot snapshot{};
  snapshot.grid_coordinates = packed_geometry.grid_coordinates_host;
  snapshot.num_elements =
      packed_geometry.logical_to_grid_inverse_jacobian_host.extent(0);
  snapshot.points_per_element =
      packed_geometry.logical_to_grid_inverse_jacobian_host.extent(1);
  snapshot.logical_to_grid_inverse_jacobian.resize(
      snapshot.num_elements * snapshot.points_per_element * Dim * Dim);
  size_t offset = 0;
  for (size_t e = 0; e < snapshot.num_elements; ++e) {
    for (size_t p = 0; p < snapshot.points_per_element; ++p) {
      for (size_t component = 0; component < Dim * Dim; ++component) {
        snapshot.logical_to_grid_inverse_jacobian[offset++] =
            packed_geometry.logical_to_grid_inverse_jacobian_host(e, p,
                                                                  component);
      }
    }
  }
  return snapshot;
}

void check_domain_static_caches_unchanged(
    const DomainStaticCacheSnapshot& expected_snapshot,
    const packed_geometry_type& packed_geometry) {
  for (size_t d = 0; d < Dim; ++d) {
    REQUIRE(expected_snapshot.grid_coordinates[d].size() ==
            packed_geometry.grid_coordinates_host[d].size());
    for (size_t s = 0; s < expected_snapshot.grid_coordinates[d].size(); ++s) {
      CHECK(expected_snapshot.grid_coordinates[d][s] ==
            packed_geometry.grid_coordinates_host[d][s]);
    }
  }

  REQUIRE(expected_snapshot.num_elements ==
          packed_geometry.logical_to_grid_inverse_jacobian_host.extent(0));
  REQUIRE(expected_snapshot.points_per_element ==
          packed_geometry.logical_to_grid_inverse_jacobian_host.extent(1));
  REQUIRE(Dim * Dim ==
          packed_geometry.logical_to_grid_inverse_jacobian_host.extent(2));
  size_t offset = 0;
  for (size_t e = 0; e < expected_snapshot.num_elements; ++e) {
    for (size_t p = 0; p < expected_snapshot.points_per_element; ++p) {
      for (size_t component = 0; component < Dim * Dim; ++component) {
        CHECK(expected_snapshot.logical_to_grid_inverse_jacobian[offset++] ==
              packed_geometry.logical_to_grid_inverse_jacobian_host(e, p,
                                                                    component));
      }
    }
  }
}

template <typename DeviceGammaType>
Scalar<DataVector> copy_gamma_to_host(const DeviceGammaType& device_gamma) {
  const size_t num_points = get(device_gamma).extent(0);
  Scalar<DataVector> host_gamma{num_points};
  auto host_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{},
                                                         get(device_gamma));
  for (size_t s = 0; s < num_points; ++s) {
    get(host_gamma)[s] = host_mirror(s);
  }
  return host_gamma;
}

void check_frame_fields_and_damping_at_time(
    const packed_topology_type& packed_topology,
    const packed_geometry_type& packed_geometry,
    const device_gamma0_type& device_gamma0,
    const device_gamma1_type& device_gamma1,
    const device_gamma2_type& device_gamma2,
    const std::array<double, Dim>& translation_velocity, const double time,
    const double expected_gamma0, const double expected_gamma1,
    const double expected_gamma2) {
  const size_t total_points = packed_topology.total_points;
  REQUIRE(total_points > 0);

  const auto host_inverse_jacobian = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, packed_geometry.element_inverse_jacobian_device);
  const auto host_gamma0 = copy_gamma_to_host(device_gamma0);
  const auto host_gamma1 = copy_gamma_to_host(device_gamma1);
  const auto host_gamma2 = copy_gamma_to_host(device_gamma2);

  const Approx approx = Approx::custom().epsilon(1.0e-12).scale(1.0);
  for (size_t d = 0; d < Dim; ++d) {
    const auto inertial_component_device = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{},
        packed_geometry.inertial_coordinates_device.get(d));
    const auto mesh_velocity_component_device =
        Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace{}, packed_geometry.mesh_velocity_device.get(d));
    for (size_t s = 0; s < total_points; ++s) {
      const double grid_coord = packed_geometry.grid_coordinates_host[d][s];
      CHECK(packed_geometry.inertial_coordinates_host[d][s] ==
            approx(grid_coord + translation_velocity[d] * time));
      CHECK(inertial_component_device(s) ==
            approx(grid_coord + translation_velocity[d] * time));
      CHECK(mesh_velocity_component_device(s) ==
            approx(translation_velocity[d]));
    }
  }

  for (size_t e = 0; e < packed_topology.local_element_ids.size(); ++e) {
    const size_t point_begin = packed_topology.element_point_offsets_host[e];
    const size_t point_end = packed_topology.element_point_offsets_host[e + 1];
    const size_t element_points = point_end - point_begin;
    for (size_t p = 0; p < element_points; ++p) {
      for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
        for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
          const size_t component = logical_d * Dim + inertial_d;
          CHECK(host_inverse_jacobian(e, p, component) ==
                approx(packed_geometry.logical_to_grid_inverse_jacobian_host(
                    e, p, component)));
        }
      }
    }
  }

  for (size_t s = 0; s < total_points; ++s) {
    CHECK(get(host_gamma0)[s] == approx(expected_gamma0));
    CHECK(get(host_gamma1)[s] == approx(expected_gamma1));
    CHECK(get(host_gamma2)[s] == approx(expected_gamma2));
  }
}

void test_update_frame_fields_and_constraint_damping_batched() {
  const std::array<double, Dim> velocity{{0.25, -0.15, 0.05}};
  const domain::creators::Rectilinear<Dim> domain_creator{
      {{-1.2, -0.8, -0.6}},
      {{1.1, 0.9, 0.7}},
      {{1, 0, 0}},
      {{4, 4, 4}},
      {{true, true, true}},
      {},
      std::make_unique<
          domain::creators::time_dependence::UniformTranslation<Dim, 0>>(
          0.0, velocity)};

  const auto domain = domain_creator.create_domain();
  const auto initial_refinement_levels =
      domain_creator.initial_refinement_levels();
  const auto initial_extents = domain_creator.initial_extents();
  const auto functions_of_time = domain_creator.functions_of_time();

  packed_topology_type packed_topology{};
  gh::Batched::Initialization::DriverState::apply(
      make_not_null(&packed_topology), domain, initial_refinement_levels,
      initial_extents, Spectral::Quadrature::GaussLobatto);

  packed_geometry_type packed_geometry{};
  packed_evolution_state_type packed_evolution_state{};
  device_gamma0_type device_gamma0{};
  device_gamma1_type device_gamma1{};
  device_gamma2_type device_gamma2{};

  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  const hardcoded_solution solution{1.0, zero_spin, zero_center};
  const gh::gauges::Harmonic gauge_condition{};

  constexpr double gamma0_constant = 0.4;
  constexpr double gamma1_constant = -0.7;
  constexpr double gamma2_constant = 1.3;
  const ConstraintDamping::Constant<Dim, Frame::Grid> damping_gamma0{
      gamma0_constant};
  const ConstraintDamping::Constant<Dim, Frame::Grid> damping_gamma1{
      gamma1_constant};
  const ConstraintDamping::Constant<Dim, Frame::Grid> damping_gamma2{
      gamma2_constant};

  gh::Actions::InitializeBatchedData::apply(
      make_not_null(&packed_geometry), make_not_null(&packed_evolution_state),
      make_not_null(&device_gamma0), make_not_null(&device_gamma1),
      make_not_null(&device_gamma2), packed_topology, domain, functions_of_time,
      solution, 0.0, gauge_condition, damping_gamma0, damping_gamma1,
      damping_gamma2);
  const auto domain_static_cache_snapshot =
      snapshot_domain_static_caches(packed_geometry);

  const std::array<double, 2> test_times{{0.3, 1.6}};
  REQUIRE(packed_topology.local_element_ids.size() >= 2);
  for (size_t time_index = 0; time_index < test_times.size(); ++time_index) {
    const double time = test_times[time_index];
    gh::Actions::UpdateFrameFieldsBatched::apply(
        make_not_null(&packed_geometry), packed_topology, domain,
        functions_of_time, time);
    REQUIRE(packed_geometry.element_map_kind_host.size() ==
            packed_topology.local_element_ids.size());
    for (size_t e = 0; e < packed_geometry.element_map_kind_host.size(); ++e) {
      const auto expected_kind =
          time_index == 1 and e == 0
              ? gh::Actions::BatchedMapKind::HostEvaluated
              : gh::Actions::BatchedMapKind::ComposedAffine;
      CHECK(packed_geometry.element_map_kind_host[e] ==
            static_cast<int>(expected_kind));
    }
    check_domain_static_caches_unchanged(domain_static_cache_snapshot,
                                         packed_geometry);
    gh::Actions::UpdateConstraintDampingBatched::apply(
        make_not_null(&device_gamma0), make_not_null(&device_gamma1),
        make_not_null(&device_gamma2), packed_geometry, functions_of_time, time,
        damping_gamma0, damping_gamma1, damping_gamma2);
    check_domain_static_caches_unchanged(domain_static_cache_snapshot,
                                         packed_geometry);

    check_frame_fields_and_damping_at_time(
        packed_topology, packed_geometry, device_gamma0, device_gamma1,
        device_gamma2, velocity, time, gamma0_constant, gamma1_constant,
        gamma2_constant);

    if (time_index == 0) {
      // Simulate a mixed-map domain at the second time: the first element uses
      // the host representation while the remaining elements retain the
      // verified device representation of the same map. This isolates the
      // mixed dispatch and checks that the combined device/host geometry is
      // coherent without requiring a large multi-block domain fixture.
      packed_geometry.element_map_kind_host[0] =
          static_cast<int>(gh::Actions::BatchedMapKind::HostEvaluated);
      auto host_map_kinds = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{}, packed_geometry.element_map_kind_device);
      host_map_kinds(0) =
          static_cast<int>(gh::Actions::BatchedMapKind::HostEvaluated);
      Kokkos::deep_copy(packed_geometry.element_map_kind_device,
                        host_map_kinds);
    }
  }
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic."
    "UpdateFrameFieldsAndConstraintDampingBatched",
    "[Unit][Evolution]") {
  test_update_frame_fields_and_constraint_damping_batched();
}
