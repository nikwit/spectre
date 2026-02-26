// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/AnalyticChristoffel.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeTimeDerivativeBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace {

constexpr size_t Dim = 3;
using system = gh::System<Dim>;
using hardcoded_solution = gh::Solutions::WrappedGr<gr::Solutions::KerrSchild>;

using variables_tag = typename system::variables_tag;
using variables_type = typename variables_tag::type;
using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
using dt_variables_type = typename dt_variables_tag::type;
using derivative_variables_tags =
    db::wrap_tags_in<::Tags::deriv, typename system::gradient_variables,
                     tmpl::size_t<Dim>, Frame::Inertial>;

using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, Dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, Dim, Frame::Inertial>;

using d_spacetime_metric_tag =
    ::Tags::deriv<spacetime_metric_tag, tmpl::size_t<Dim>, Frame::Inertial>;
using d_pi_tag = ::Tags::deriv<pi_tag, tmpl::size_t<Dim>, Frame::Inertial>;
using d_phi_tag = ::Tags::deriv<phi_tag, tmpl::size_t<Dim>, Frame::Inertial>;

dt_variables_type compute_host_reference(
    const variables_type& vars, const gh::Tags::ConstraintGamma0::type& gamma0,
    const gh::Tags::ConstraintGamma1::type& gamma1,
    const gh::Tags::ConstraintGamma2::type& gamma2, const Mesh<Dim>& mesh,
    const double time,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian) {
  const size_t num_points = mesh.number_of_grid_points();

  Variables<derivative_variables_tags> host_partial_derivs{num_points};
  partial_derivatives(make_not_null(&host_partial_derivs), vars, mesh,
                      inverse_jacobian);

  using all_solutions = gh::Solutions::all_solutions<Dim>;
  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  gh::gauges::AnalyticChristoffel gauge_condition{
      std::make_unique<hardcoded_solution>(1.0, zero_spin, zero_center)};

  tnsr::aa<DataVector, Dim, Frame::Inertial> host_dt_spacetime_metric{
      num_points};
  tnsr::aa<DataVector, Dim, Frame::Inertial> host_dt_pi{num_points};
  tnsr::iaa<DataVector, Dim, Frame::Inertial> host_dt_phi{num_points};
  Variables<typename gh::TimeDerivative<all_solutions, Dim>::temporary_tags>
      host_buffer{num_points};

  gh::TimeDerivative<all_solutions, Dim>::apply(
      make_not_null(&host_dt_spacetime_metric), make_not_null(&host_dt_pi),
      make_not_null(&host_dt_phi),
      make_not_null(&get<gh::Tags::ConstraintGamma1>(host_buffer)),
      make_not_null(&get<gh::Tags::ConstraintGamma2>(host_buffer)),
      make_not_null(&get<gh::Tags::GaugeH<DataVector, Dim>>(host_buffer)),
      make_not_null(
          &get<gh::Tags::SpacetimeDerivGaugeH<DataVector, Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::Gamma1Gamma2>(host_buffer)),
      make_not_null(&get<gh::Tags::HalfPiTwoNormals>(host_buffer)),
      make_not_null(&get<gh::Tags::NormalDotOneIndexConstraint>(host_buffer)),
      make_not_null(&get<gh::Tags::Gamma1Plus1>(host_buffer)),
      make_not_null(&get<gh::Tags::PiOneNormal<Dim>>(host_buffer)),
      make_not_null(
          &get<gh::Tags::GaugeConstraint<DataVector, Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::HalfPhiTwoNormals<Dim>>(host_buffer)),
      make_not_null(
          &get<gh::Tags::ShiftDotThreeIndexConstraint<Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::MeshVelocityDotThreeIndexConstraint<Dim>>(
          host_buffer)),
      make_not_null(&get<gh::Tags::PhiOneNormal<Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::PiSecondIndexUp<Dim>>(host_buffer)),
      make_not_null(
          &get<gh::Tags::ThreeIndexConstraint<DataVector, Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::PhiFirstIndexUp<Dim>>(host_buffer)),
      make_not_null(&get<gh::Tags::PhiThirdIndexUp<Dim>>(host_buffer)),
      make_not_null(
          &get<gh::Tags::SpacetimeChristoffelFirstKindThirdIndexUp<Dim>>(
              host_buffer)),
      make_not_null(&get<gr::Tags::Lapse<DataVector>>(host_buffer)),
      make_not_null(&get<gr::Tags::Shift<DataVector, Dim>>(host_buffer)),
      make_not_null(
          &get<gr::Tags::InverseSpatialMetric<DataVector, Dim>>(host_buffer)),
      make_not_null(&get<gr::Tags::DetSpatialMetric<DataVector>>(host_buffer)),
      make_not_null(
          &get<gr::Tags::SqrtDetSpatialMetric<DataVector>>(host_buffer)),
      make_not_null(
          &get<gr::Tags::InverseSpacetimeMetric<DataVector, Dim>>(host_buffer)),
      make_not_null(
          &get<gr::Tags::SpacetimeChristoffelFirstKind<DataVector, Dim>>(
              host_buffer)),
      make_not_null(
          &get<gr::Tags::SpacetimeChristoffelSecondKind<DataVector, Dim>>(
              host_buffer)),
      make_not_null(
          &get<gr::Tags::TraceSpacetimeChristoffelFirstKind<DataVector, Dim>>(
              host_buffer)),
      make_not_null(
          &get<gr::Tags::SpacetimeNormalVector<DataVector, Dim>>(host_buffer)),
      get<d_spacetime_metric_tag>(host_partial_derivs),
      get<d_pi_tag>(host_partial_derivs), get<d_phi_tag>(host_partial_derivs),
      get<spacetime_metric_tag>(vars), get<pi_tag>(vars), get<phi_tag>(vars),
      gamma0, gamma1, gamma2, gauge_condition, mesh, time, inertial_coords,
      inverse_jacobian,
      std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>{std::nullopt});

  dt_variables_type dt_reference{num_points, 0.0};
  get<::Tags::dt<spacetime_metric_tag>>(dt_reference) =
      host_dt_spacetime_metric;
  get<::Tags::dt<pi_tag>>(dt_reference) = host_dt_pi;
  get<::Tags::dt<phi_tag>>(dt_reference) = host_dt_phi;
  return dt_reference;
}

void test_compute_time_derivative_batched_matches_host() {
  const Mesh<Dim> mesh{std::array<size_t, Dim>{{4, 3, 5}},
                       Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  const size_t num_points = mesh.number_of_grid_points();

  const auto logical_coords = logical_coordinates(mesh);
  tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords{num_points};
  for (size_t d = 0; d < Dim; ++d) {
    inertial_coords.get(d) = logical_coords.get(d);
  }

  InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian{num_points};
  for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
    for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
      inverse_jacobian.get(logical_d, inertial_d) =
          DataVector{num_points, logical_d == inertial_d ? 1.0 : 0.0};
    }
  }

  const Element<Dim> element{ElementId<Dim>{0}, {}};
  const double time = 0.2;
  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  const hardcoded_solution solution{1.0, zero_spin, zero_center};

  const auto analytic_vars =
      solution.variables(inertial_coords, time,
                         tmpl::list<spacetime_metric_tag, pi_tag, phi_tag>{});
  variables_type host_vars{num_points, 0.0};
  get<spacetime_metric_tag>(host_vars) =
      get<spacetime_metric_tag>(analytic_vars);
  get<pi_tag>(host_vars) = get<pi_tag>(analytic_vars);
  get<phi_tag>(host_vars) = get<phi_tag>(analytic_vars);

  gh::Tags::ConstraintGamma0::type host_gamma0{num_points, 0.0};
  gh::Tags::ConstraintGamma1::type host_gamma1{num_points, 0.0};
  gh::Tags::ConstraintGamma2::type host_gamma2{num_points, 0.0};
  get(host_gamma0) =
      0.9 + 0.1 * inertial_coords.get(0) - 0.05 * inertial_coords.get(1);
  get(host_gamma1) = -0.1 + 0.04 * inertial_coords.get(2);
  get(host_gamma2) =
      0.5 + 0.07 * inertial_coords.get(0) + 0.02 * inertial_coords.get(1);

  evolution::Kokkos::PackedTopology<system> packed_topology{};
  packed_topology.local_element_ids.push_back(element.id());
  packed_topology.local_elements.push_back(element);
  packed_topology.element_index_by_id.insert_or_assign(element.id(), 0);
  packed_topology.element_extents_host.push_back(mesh.extents().indices());
  packed_topology.element_point_offsets_host = {0, num_points};
  packed_topology.total_points = num_points;
  packed_topology.points_per_element = num_points;
  packed_topology.uniform_extents_host = mesh.extents().indices();
  packed_topology.uniform_basis_host = {mesh.basis(0), mesh.basis(1),
                                        mesh.basis(2)};
  packed_topology.uniform_quadrature_host = {
      mesh.quadrature(0), mesh.quadrature(1), mesh.quadrature(2)};

  evolution::Kokkos::PackedGeometry<system> packed_geometry{};
  packed_geometry.inertial_coordinates_host = {
      {inertial_coords.get(0), inertial_coords.get(1), inertial_coords.get(2)}};
  packed_geometry.element_inverse_jacobian_device = ::Kokkos::View<double***>(
      "TestGhBatchedElementInverseJacobian", 1, num_points, Dim * Dim);
  auto host_inverse_jacobian = ::Kokkos::create_mirror_view(
      packed_geometry.element_inverse_jacobian_device);
  for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
    for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
      for (size_t s = 0; s < num_points; ++s) {
        host_inverse_jacobian(0, s, logical_d * Dim + inertial_d) =
            inverse_jacobian.get(logical_d, inertial_d)[s];
      }
    }
  }
  ::Kokkos::deep_copy(packed_geometry.element_inverse_jacobian_device,
                      host_inverse_jacobian);

  evolution::Kokkos::PackedEvolutionState<system> packed_evolution_state{};
  packed_evolution_state.device_variables = copy_to_device(host_vars);
  packed_evolution_state.device_dt_variables =
      typename evolution::Kokkos::PackedEvolutionState<
          system>::device_dt_variables_type{num_points};
  ::Kokkos::deep_copy(packed_evolution_state.device_dt_variables.view(), 0.0);

  typename gh::KokkosTags::DeviceConstraintGamma0::type device_gamma0{
      "DeviceGamma0", num_points};
  typename gh::KokkosTags::DeviceConstraintGamma1::type device_gamma1{
      "DeviceGamma1", num_points};
  typename gh::KokkosTags::DeviceConstraintGamma2::type device_gamma2{
      "DeviceGamma2", num_points};

  auto host_gamma0_mirror = ::Kokkos::create_mirror_view(get(device_gamma0));
  auto host_gamma1_mirror = ::Kokkos::create_mirror_view(get(device_gamma1));
  auto host_gamma2_mirror = ::Kokkos::create_mirror_view(get(device_gamma2));
  std::copy_n(get(host_gamma0).data(), num_points, host_gamma0_mirror.data());
  std::copy_n(get(host_gamma1).data(), num_points, host_gamma1_mirror.data());
  std::copy_n(get(host_gamma2).data(), num_points, host_gamma2_mirror.data());
  ::Kokkos::deep_copy(get(device_gamma0), host_gamma0_mirror);
  ::Kokkos::deep_copy(get(device_gamma1), host_gamma1_mirror);
  ::Kokkos::deep_copy(get(device_gamma2), host_gamma2_mirror);

  gh::Actions::ComputeTimeDerivativeBatched::apply(
      make_not_null(&packed_evolution_state), packed_topology, packed_geometry,
      device_gamma0, device_gamma1, device_gamma2);

  dt_variables_type dt_batched{num_points, 0.0};
  copy_to_host(make_not_null(&dt_batched),
               packed_evolution_state.device_dt_variables);
  const dt_variables_type dt_reference =
      compute_host_reference(host_vars, host_gamma0, host_gamma1, host_gamma2,
                             mesh, time, inertial_coords, inverse_jacobian);

  const Approx approx = Approx::custom().epsilon(1e-12).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get<::Tags::dt<spacetime_metric_tag>>(dt_batched),
      get<::Tags::dt<spacetime_metric_tag>>(dt_reference), approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<pi_tag>>(dt_batched),
                               get<::Tags::dt<pi_tag>>(dt_reference), approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<phi_tag>>(dt_batched),
                               get<::Tags::dt<phi_tag>>(dt_reference), approx);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.ComputeTimeDerivativeBatched",
    "[Unit][Evolution]") {
  test_compute_time_derivative_batched_matches_host();
}
