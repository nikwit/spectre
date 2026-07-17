// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <vector>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/HostElementView.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {

constexpr size_t Dim = 3;
using system = gh::System<Dim>;
using packed_topology_type = evolution::Kokkos::PackedTopology<system>;
using packed_geometry_type = evolution::Kokkos::PackedGeometry<system>;
using packed_evolution_state_type =
    evolution::Kokkos::PackedEvolutionState<system>;
using variables_type = typename system::variables_tag::type;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, Dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, Dim, Frame::Inertial>;
using deriv_phi_tag =
    ::Tags::deriv<phi_tag, tmpl::size_t<Dim>, Frame::Inertial>;

struct PackedFixture {
  packed_topology_type topology{};
  packed_geometry_type geometry{};
  packed_evolution_state_type state{};
  variables_type host_vars{};
};

double packed_value(const size_t component, const size_t point) {
  return 1000.0 * static_cast<double>(component + 1) +
         static_cast<double>(point);
}

PackedFixture make_fixture() {
  PackedFixture fixture{};
  const std::array<size_t, Dim> extents{{3, 3, 3}};
  const Mesh<Dim> mesh{extents, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  const size_t points_per_element = mesh.number_of_grid_points();
  const size_t num_elements = 2;
  const size_t total_points = num_elements * points_per_element;

  fixture.topology.local_element_ids = {
      ElementId<Dim>{0, {{SegmentId{1, 0}, SegmentId{0, 0}, SegmentId{0, 0}}}},
      ElementId<Dim>{0, {{SegmentId{1, 1}, SegmentId{0, 0}, SegmentId{0, 0}}}}};
  fixture.topology.local_elements = {
      Element<Dim>{fixture.topology.local_element_ids[0], {}},
      Element<Dim>{fixture.topology.local_element_ids[1], {}}};
  fixture.topology.element_index_by_id = {
      {fixture.topology.local_element_ids[0], 0},
      {fixture.topology.local_element_ids[1], 1}};
  fixture.topology.element_extents_host = {extents, extents};
  fixture.topology.element_point_offsets_host = {0, points_per_element,
                                                 total_points};
  fixture.topology.total_points = total_points;
  fixture.topology.points_per_element = points_per_element;
  fixture.topology.uniform_extents_host = extents;
  fixture.topology.uniform_basis_host = mesh.basis();
  fixture.topology.uniform_quadrature_host = mesh.quadrature();

  fixture.host_vars = variables_type{total_points};
  for (size_t component = 0;
       component < variables_type::number_of_independent_components;
       ++component) {
    for (size_t point = 0; point < total_points; ++point) {
      fixture.host_vars.data()[component * total_points + point] =
          packed_value(component, point);
    }
  }
  fixture.state.device_variables = copy_to_device(fixture.host_vars);

  fixture.geometry.grid_coordinates_host = {{DataVector{total_points},
                                             DataVector{total_points},
                                             DataVector{total_points}}};
  fixture.geometry.distorted_coordinates_host = {{DataVector{total_points},
                                                  DataVector{total_points},
                                                  DataVector{total_points}}};
  fixture.geometry.inertial_coordinates_host = {{DataVector{total_points},
                                                 DataVector{total_points},
                                                 DataVector{total_points}}};
  fixture.geometry.mesh_velocity_host = {{DataVector{total_points},
                                          DataVector{total_points},
                                          DataVector{total_points}}};
  fixture.geometry.grid_to_distorted_frame_velocity_host = {
      {DataVector{total_points}, DataVector{total_points},
       DataVector{total_points}}};
  for (size_t d = 0; d < Dim; ++d) {
    for (size_t point = 0; point < total_points; ++point) {
      fixture.geometry.grid_coordinates_host[d][point] =
          10.0 * static_cast<double>(d + 1) + static_cast<double>(point);
      fixture.geometry.distorted_coordinates_host[d][point] =
          15.0 * static_cast<double>(d + 1) + static_cast<double>(point);
      fixture.geometry.inertial_coordinates_host[d][point] =
          20.0 * static_cast<double>(d + 1) + static_cast<double>(point);
      fixture.geometry.mesh_velocity_host[d][point] =
          30.0 * static_cast<double>(d + 1) + static_cast<double>(point);
      fixture.geometry.grid_to_distorted_frame_velocity_host[d][point] =
          35.0 * static_cast<double>(d + 1) + static_cast<double>(point);
    }
  }

  using host_jacobian_type =
      packed_geometry_type::host_inverse_jacobian_cache_type;
  fixture.geometry.logical_to_distorted_inverse_jacobian_host =
      host_jacobian_type{"TestGhHostElementViewLogicalToDistortedInvJacobian",
                         num_elements, points_per_element, Dim * Dim};
  fixture.geometry.grid_to_distorted_inverse_jacobian_host =
      host_jacobian_type{"TestGhHostElementViewGridToDistortedInvJacobian",
                         num_elements, points_per_element, Dim * Dim};
  fixture.geometry.distorted_to_inertial_jacobian_host =
      host_jacobian_type{"TestGhHostElementViewDistortedToInertialJacobian",
                         num_elements, points_per_element, Dim * Dim};
  fixture.geometry.distorted_to_inertial_inverse_jacobian_host =
      host_jacobian_type{"TestGhHostElementViewDistortedToInertialInvJacobian",
                         num_elements, points_per_element, Dim * Dim};
  for (size_t e = 0; e < num_elements; ++e) {
    for (size_t p = 0; p < points_per_element; ++p) {
      for (size_t source_d = 0; source_d < Dim; ++source_d) {
        for (size_t target_d = 0; target_d < Dim; ++target_d) {
          const size_t component = source_d * Dim + target_d;
          const double offset = static_cast<double>(10 * e + p + component);
          fixture.geometry.logical_to_distorted_inverse_jacobian_host(
              e, p, component) = 100.0 + offset;
          fixture.geometry.grid_to_distorted_inverse_jacobian_host(
              e, p, component) = 200.0 + offset;
          fixture.geometry.distorted_to_inertial_jacobian_host(
              e, p, component) = 300.0 + offset;
          fixture.geometry.distorted_to_inertial_inverse_jacobian_host(
              e, p, component) = 400.0 + offset;
        }
      }
    }
  }

  fixture.geometry.element_inverse_jacobian_device =
      packed_geometry_type::device_inverse_jacobian_type{
          "TestGhHostElementViewInverseJacobian", num_elements,
          points_per_element, Dim * Dim};
  auto host_inverse_jacobian = Kokkos::create_mirror_view(
      fixture.geometry.element_inverse_jacobian_device);
  for (size_t e = 0; e < num_elements; ++e) {
    for (size_t p = 0; p < points_per_element; ++p) {
      for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
        for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
          host_inverse_jacobian(e, p, logical_d * Dim + inertial_d) =
              logical_d == inertial_d
                  ? 1.0 + 0.1 * static_cast<double>(e)
                  : 0.01 *
                        static_cast<double>(1 + e + p + logical_d + inertial_d);
        }
      }
    }
  }
  Kokkos::deep_copy(fixture.geometry.element_inverse_jacobian_device,
                    host_inverse_jacobian);

  return fixture;
}

void check_view_matches_fixture(const gh::Batched::HostElementView& view,
                                const PackedFixture& fixture,
                                const size_t element_index) {
  const size_t point_begin =
      fixture.topology.element_point_offsets_host[element_index];
  const size_t points_per_element = fixture.topology.points_per_element;
  CHECK(view.element_id == fixture.topology.local_element_ids[element_index]);
  CHECK(view.element == fixture.topology.local_elements[element_index]);
  CHECK(view.mesh.extents() ==
        Index<Dim>{fixture.topology.element_extents_host[element_index]});
  CHECK(view.mesh.number_of_grid_points() == points_per_element);

  for (size_t component = 0;
       component < variables_type::number_of_independent_components;
       ++component) {
    for (size_t p = 0; p < points_per_element; ++p) {
      CHECK(view.evolved_vars.data()[component * points_per_element + p] ==
            packed_value(component, point_begin + p));
    }
  }

  for (size_t d = 0; d < Dim; ++d) {
    for (size_t p = 0; p < points_per_element; ++p) {
      CHECK(view.grid_coordinates.get(d)[p] ==
            fixture.geometry.grid_coordinates_host[d][point_begin + p]);
      CHECK(view.distorted_coordinates.get(d)[p] ==
            fixture.geometry.distorted_coordinates_host[d][point_begin + p]);
      CHECK(view.inertial_coordinates.get(d)[p] ==
            fixture.geometry.inertial_coordinates_host[d][point_begin + p]);
      CHECK(view.mesh_velocity.get(d)[p] ==
            fixture.geometry.mesh_velocity_host[d][point_begin + p]);
      CHECK(view.grid_to_distorted_frame_velocity.get(d)[p] ==
            fixture.geometry
                .grid_to_distorted_frame_velocity_host[d][point_begin + p]);
    }
  }

  const auto host_inverse_jacobian = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, fixture.geometry.element_inverse_jacobian_device);
  for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
    for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
      for (size_t p = 0; p < points_per_element; ++p) {
        CHECK(view.logical_to_inertial_inverse_jacobian.get(logical_d,
                                                            inertial_d)[p] ==
              host_inverse_jacobian(element_index, p,
                                    logical_d * Dim + inertial_d));
      }
    }
  }
  for (size_t source_d = 0; source_d < Dim; ++source_d) {
    for (size_t target_d = 0; target_d < Dim; ++target_d) {
      const size_t component = source_d * Dim + target_d;
      for (size_t p = 0; p < points_per_element; ++p) {
        CHECK(view.logical_to_distorted_inverse_jacobian.get(source_d,
                                                             target_d)[p] ==
              fixture.geometry.logical_to_distorted_inverse_jacobian_host(
                  element_index, p, component));
        CHECK(view.grid_to_distorted_inverse_jacobian.get(source_d,
                                                          target_d)[p] ==
              fixture.geometry.grid_to_distorted_inverse_jacobian_host(
                  element_index, p, component));
        CHECK(view.distorted_to_inertial_jacobian.get(source_d, target_d)[p] ==
              fixture.geometry.distorted_to_inertial_jacobian_host(
                  element_index, p, component));
        CHECK(view.distorted_to_inertial_inverse_jacobian.get(source_d,
                                                              target_d)[p] ==
              fixture.geometry.distorted_to_inertial_inverse_jacobian_host(
                  element_index, p, component));
      }
    }
  }
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.HostElementViewBatched",
    "[Unit][Evolution]") {
  const auto fixture = make_fixture();
  size_t element_index = 0;
  gh::Batched::for_each_host_element_view(
      fixture.topology, fixture.geometry, fixture.state,
      [&fixture, &element_index](const gh::Batched::HostElementView& view) {
        check_view_matches_fixture(view, fixture, element_index);

        const auto source_vars = gh::Batched::make_ah_source_vars(view);
        CHECK(source_vars.number_of_grid_points() ==
              view.mesh.number_of_grid_points());
        CHECK_ITERABLE_APPROX(get<spacetime_metric_tag>(source_vars),
                              get<spacetime_metric_tag>(view.evolved_vars));
        CHECK_ITERABLE_APPROX(get<pi_tag>(source_vars),
                              get<pi_tag>(view.evolved_vars));
        CHECK_ITERABLE_APPROX(get<phi_tag>(source_vars),
                              get<phi_tag>(view.evolved_vars));

        const auto expected_deriv_phi =
            partial_derivative(get<phi_tag>(view.evolved_vars), view.mesh,
                               view.logical_to_inertial_inverse_jacobian);
        CHECK_ITERABLE_APPROX(get<deriv_phi_tag>(source_vars),
                              expected_deriv_phi);
        ++element_index;
      });
  CHECK(element_index == fixture.topology.local_element_ids.size());
}
