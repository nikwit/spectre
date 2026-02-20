// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/AnalyticChristoffel.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Time/Slab.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
constexpr size_t Dim = 3;
using system = gh::System<Dim>;

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

using external_boundary_conditions_type =
    typename domain::Tags::ExternalBoundaryConditions<Dim>::type;

void test_volume_terms_match_host_time_derivative() {
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

  const Slab slab{0.0, 1.0};
  const TimeStepId time_step_id{true, 0, slab.start()};
  const double time = 0.2;
  const Element<Dim> element{ElementId<Dim>{0}, {}};

  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  using hardcoded_solution =
      gh::Solutions::WrappedGr<gr::Solutions::KerrSchild>;
  const hardcoded_solution solution{1.0, zero_spin, zero_center};

  const auto analytic_vars =
      solution.variables(inertial_coords, time,
                         tmpl::list<spacetime_metric_tag, pi_tag, phi_tag>{});
  variables_type host_vars{num_points, 0.0};
  get<spacetime_metric_tag>(host_vars) =
      get<spacetime_metric_tag>(analytic_vars);
  get<pi_tag>(host_vars) = get<pi_tag>(analytic_vars);
  get<phi_tag>(host_vars) = get<phi_tag>(analytic_vars);

  gh::Tags::ConstraintGamma0::type host_gamma0{num_points};
  gh::Tags::ConstraintGamma1::type host_gamma1{num_points};
  gh::Tags::ConstraintGamma2::type host_gamma2{num_points};
  get(host_gamma0) =
      0.9 + 0.1 * inertial_coords.get(0) - 0.05 * inertial_coords.get(1);
  get(host_gamma1) = -0.1 + 0.04 * inertial_coords.get(2);
  get(host_gamma2) =
      0.5 + 0.07 * inertial_coords.get(0) + 0.02 * inertial_coords.get(1);

  Variables<derivative_variables_tags> host_partial_derivs{num_points};
  partial_derivatives(make_not_null(&host_partial_derivs), host_vars, mesh,
                      inverse_jacobian);

  using all_solutions = gh::Solutions::all_solutions<Dim>;
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
      get<spacetime_metric_tag>(host_vars), get<pi_tag>(host_vars),
      get<phi_tag>(host_vars), host_gamma0, host_gamma1, host_gamma2,
      gauge_condition, mesh, time, inertial_coords, inverse_jacobian,
      std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>{std::nullopt});

  typename gh::KokkosTags::DeviceVariables<system>::type device_vars =
      copy_to_device(host_vars);
  typename gh::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  typename gh::KokkosTags::DeviceInverseJacobian<Dim>::type
      device_inverse_jacobian{};
  typename gh::KokkosTags::DeviceInertialCoordinates<Dim>::type
      device_inertial_coordinates{};
  typename gh::KokkosTags::DeviceConstraintGamma0::type device_gamma0{};
  typename gh::KokkosTags::DeviceConstraintGamma1::type device_gamma1{};
  typename gh::KokkosTags::DeviceConstraintGamma2::type device_gamma2{};
  typename gh::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type
      device_face_to_volume_index_map{};
  typename gh::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
      device_face_unit_normal_covector{};
  typename gh::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
      device_face_normal_magnitude{};
  typename gh::KokkosTags::DeviceMortarData<Dim>::type device_mortar_data{};
  typename evolution::dg::Tags::MortarMesh<Dim>::type mortar_meshes{};
  typename evolution::dg::Tags::MortarInfo<Dim>::type mortar_infos{};

  gh::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&device_inverse_jacobian),
      make_not_null(&device_inertial_coordinates),
      make_not_null(&device_gamma0), make_not_null(&device_gamma1),
      make_not_null(&device_gamma2),
      make_not_null(&device_face_to_volume_index_map),
      make_not_null(&device_face_unit_normal_covector),
      make_not_null(&device_face_normal_magnitude),
      make_not_null(&device_mortar_data), inverse_jacobian, host_gamma0,
      host_gamma1, host_gamma2, inertial_coords, mesh, element, mortar_meshes,
      mortar_infos);

  typename gh::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
      outgoing_boundary_data{};
  typename gh::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type
      external_boundary_data{};
  external_boundary_conditions_type external_boundary_conditions_by_block{1};

  gh::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      make_not_null(&external_boundary_data), device_vars,
      device_inverse_jacobian, device_inertial_coordinates, device_gamma0,
      device_gamma1, device_gamma2, device_face_to_volume_index_map,
      device_face_unit_normal_covector, device_mortar_data, mortar_meshes,
      host_gamma0, host_gamma1, host_gamma2,
      external_boundary_conditions_by_block, time, mesh, element, time_step_id);

  dt_variables_type dt_kokkos{num_points, 0.0};
  copy_to_host(make_not_null(&dt_kokkos), device_dt);

  const Approx approx = Approx::custom().epsilon(1e-8).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<spacetime_metric_tag>>(dt_kokkos),
                               host_dt_spacetime_metric, approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<pi_tag>>(dt_kokkos), host_dt_pi,
                               approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<phi_tag>>(dt_kokkos), host_dt_phi,
                               approx);

  CHECK(outgoing_boundary_data.empty());
  CHECK(external_boundary_data.empty());
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.ComputeTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_volume_terms_match_host_time_derivative();
}
