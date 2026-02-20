// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/InternalMortarDataImpl.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/UpwindPenalty.hpp"
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
using boundary_correction = gh::BoundaryCorrections::UpwindPenalty<Dim>;
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
using dg_package_field_tags = typename boundary_correction::dg_package_field_tags;
using packaged_data_by_mortar_map = dg::MortarMap<Dim, DataVector>;

using external_boundary_conditions_type =
    typename domain::Tags::ExternalBoundaryConditions<Dim>::type;

struct SetupData {
  Element<Dim> element{};
  Mesh<Dim> mesh{};
  variables_type vars{};
  gh::Tags::ConstraintGamma0::type gamma0{};
  gh::Tags::ConstraintGamma1::type gamma1{};
  gh::Tags::ConstraintGamma2::type gamma2{};
  tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords{};
  domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                Frame::Inertial>::type inverse_jacobian{};
  TimeStepId time_step_id{};
  std::unique_ptr<domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>>
      moving_mesh_map{};
  double time{};
};

struct DomainData {
  Domain<Dim> domain;
  std::vector<std::array<size_t, Dim>> initial_refinement_levels;
  std::vector<ElementId<Dim>> element_ids;
  Mesh<Dim> mesh;
};

DomainData make_internal_boundary_domain_data() {
  const domain::creators::Rectilinear<Dim> domain_creator{{{-1.0, -0.4, 0.2}},
                                                          {{1.2, 0.9, 2.3}},
                                                          {{1, 2, 1}},
                                                          {{4, 5, 4}},
                                                          {{true, true, true}}};
  auto domain = domain_creator.create_domain();
  auto initial_refinement_levels = domain_creator.initial_refinement_levels();
  auto element_ids = initial_element_ids(initial_refinement_levels);
  const Mesh<Dim> mesh{domain_creator.initial_extents()[0],
                       Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  return DomainData{std::move(domain), std::move(initial_refinement_levels),
                    std::move(element_ids), mesh};
}

SetupData make_setup_data(
    const ElementId<Dim>& element_id, const Domain<Dim>& domain_object,
    const std::vector<std::array<size_t, Dim>>& initial_refinement_levels,
    const Mesh<Dim>& mesh) {
  SetupData setup{};
  setup.element = domain::Initialization::create_initial_element(
      element_id, domain_object.blocks(), initial_refinement_levels);
  setup.mesh = mesh;

  const ElementMap<Dim, Frame::Inertial> element_map{
      setup.element.id(),
      domain_object.blocks()[setup.element.id().block_id()]};
  const auto logical_coords = logical_coordinates(setup.mesh);
  setup.inertial_coords = element_map(logical_coords);
  setup.inverse_jacobian = element_map.inv_jacobian(logical_coords);

  const Slab slab{0.0, 1.0};
  setup.time_step_id = TimeStepId{true, 0, slab.start()};
  setup.time = 0.2;

  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  const hardcoded_solution solution{1.0, zero_spin, zero_center};
  const auto analytic_vars =
      solution.variables(setup.inertial_coords, setup.time,
                         tmpl::list<spacetime_metric_tag, pi_tag, phi_tag>{});
  const size_t num_points = setup.mesh.number_of_grid_points();
  setup.vars = variables_type{num_points, 0.0};
  get<spacetime_metric_tag>(setup.vars) = get<spacetime_metric_tag>(analytic_vars);
  get<pi_tag>(setup.vars) = get<pi_tag>(analytic_vars);
  get<phi_tag>(setup.vars) = get<phi_tag>(analytic_vars);

  setup.gamma0 = gh::Tags::ConstraintGamma0::type{num_points, 0.0};
  setup.gamma1 = gh::Tags::ConstraintGamma1::type{num_points, 0.0};
  setup.gamma2 = gh::Tags::ConstraintGamma2::type{num_points, 0.0};
  get(setup.gamma0) =
      0.9 + 0.1 * setup.inertial_coords.get(0) - 0.05 * setup.inertial_coords.get(1);
  get(setup.gamma1) = -0.1 + 0.04 * setup.inertial_coords.get(2);
  get(setup.gamma2) =
      0.5 + 0.07 * setup.inertial_coords.get(0) + 0.02 * setup.inertial_coords.get(1);

  setup.moving_mesh_map =
      domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
          domain::CoordinateMaps::Identity<Dim>{});
  return setup;
}

void compute_host_reference(
    const SetupData& setup,
    const gsl::not_null<dt_variables_type*> dt_reference,
    const gsl::not_null<packaged_data_by_mortar_map*> packaged_data_by_mortar) {
  const size_t num_points = setup.mesh.number_of_grid_points();
  Variables<derivative_variables_tags> host_partial_derivs{num_points};
  partial_derivatives(make_not_null(&host_partial_derivs), setup.vars, setup.mesh,
                      setup.inverse_jacobian);

  using all_solutions = gh::Solutions::all_solutions<Dim>;
  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  gh::gauges::AnalyticChristoffel gauge_condition{
      std::make_unique<hardcoded_solution>(1.0, zero_spin, zero_center)};

  tnsr::aa<DataVector, Dim, Frame::Inertial> host_dt_spacetime_metric{num_points};
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
      get<spacetime_metric_tag>(setup.vars), get<pi_tag>(setup.vars),
      get<phi_tag>(setup.vars), setup.gamma0, setup.gamma1, setup.gamma2,
      gauge_condition, setup.mesh, setup.time, setup.inertial_coords,
      setup.inverse_jacobian,
      std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>{std::nullopt});

  dt_reference->initialize(num_points, 0.0);
  get<::Tags::dt<spacetime_metric_tag>>(*dt_reference) = host_dt_spacetime_metric;
  get<::Tags::dt<pi_tag>>(*dt_reference) = host_dt_pi;
  get<::Tags::dt<phi_tag>>(*dt_reference) = host_dt_phi;

  dg::MortarMap<Dim, Mesh<Dim>> neighbor_mesh{};
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    for (const auto& neighbor : neighbors) {
      neighbor_mesh.insert(
          {DirectionalId<Dim>{direction, neighbor}, setup.mesh});
    }
  }

  auto mortar_infos =
      evolution::dg::Initialization::detail::mortar_infos<Dim>(setup.element);
  auto [host_mortar_meshes, mortar_next_temporal_id,
        normal_covector_and_magnitude] =
      evolution::dg::Initialization::detail::mortars_apply_impl<Dim>(
          setup.element, setup.time_step_id, setup.mesh, neighbor_mesh);
  (void)mortar_next_temporal_id;
  auto host_mortar_data =
      evolution::dg::Initialization::detail::empty_mortar_data<Dim>(
          setup.element);

  using fluxes_tags =
      db::wrap_tags_in<::Tags::Flux, typename system::flux_variables,
                       tmpl::size_t<Dim>, Frame::Inertial>;
  using temporary_tags_for_face =
      typename boundary_correction::dg_package_data_temporary_tags;
  using dg_package_data_projected_tags =
      tmpl::append<typename system::variables_tag::tags_list, fluxes_tags,
                   temporary_tags_for_face, tmpl::list<>>;
  using all_face_temporary_tags = tmpl::remove_duplicates<tmpl::push_back<
      tmpl::append<dg_package_data_projected_tags,
                   evolution::dg::Actions::detail::inverse_spatial_metric_tag<
                       system>>,
      evolution::dg::Actions::detail::OneOverNormalVectorMagnitude,
      evolution::dg::Actions::detail::NormalVector<Dim>>>;
  using vars_face_temporaries = Variables<all_face_temporary_tags>;
  using dg_packaged_data_vars_on_face = Variables<dg_package_field_tags>;

  size_t max_num_face_points = 0;
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    (void)neighbors;
    max_num_face_points = std::max(
        max_num_face_points,
        setup.mesh.slice_away(direction.dimension()).number_of_grid_points());
  }

  std::vector<double> face_temporary_buffer(
      vars_face_temporaries::number_of_independent_components *
      max_num_face_points);
  std::vector<double> packaged_data_buffer(
      dg_packaged_data_vars_on_face::number_of_independent_components *
      max_num_face_points);
  auto face_temporaries_span = gsl::make_span(face_temporary_buffer);
  auto packaged_data_span = gsl::make_span(packaged_data_buffer);
  const boundary_correction correction{};
  Variables<fluxes_tags> volume_fluxes{};
  evolution::dg::Actions::detail::internal_mortar_data_impl<
      system, Dim, boundary_correction,
      typename gh::TimeDerivative<all_solutions, Dim>::temporary_tags>(
      make_not_null(&normal_covector_and_magnitude),
      make_not_null(&host_mortar_data), make_not_null(&face_temporaries_span),
      make_not_null(&packaged_data_span), correction, setup.vars, volume_fluxes,
      host_buffer, static_cast<const Variables<tmpl::list<>>*>(nullptr),
      setup.element, setup.mesh, host_mortar_meshes, mortar_infos,
      *setup.moving_mesh_map, std::nullopt, setup.inverse_jacobian);

  packaged_data_by_mortar->clear();
  for (const auto& [mortar_id, mortar_data_holder] : host_mortar_data) {
    const auto& local_mortar = mortar_data_holder.local();
    REQUIRE(local_mortar.mortar_data.has_value());
    DataVector packaged_data{local_mortar.mortar_data->size(), 0.0};
    std::copy_n(local_mortar.mortar_data->data(),
                local_mortar.mortar_data->size(), packaged_data.data());
    packaged_data_by_mortar->insert_or_assign(mortar_id,
                                              std::move(packaged_data));
  }
}

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
      device_face_unit_normal_covector, device_face_normal_magnitude,
      device_mortar_data, mortar_meshes, host_gamma0, host_gamma1, host_gamma2,
      external_boundary_conditions_by_block, time, mesh, element, time_step_id);

  dt_variables_type dt_kokkos{num_points, 0.0};
  copy_to_host(make_not_null(&dt_kokkos), device_dt);

  const Approx approx = Approx::custom().epsilon(1e-13).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<spacetime_metric_tag>>(dt_kokkos),
                               host_dt_spacetime_metric, approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<pi_tag>>(dt_kokkos), host_dt_pi,
                               approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<phi_tag>>(dt_kokkos), host_dt_phi,
                               approx);

  CHECK(outgoing_boundary_data.empty());
  CHECK(external_boundary_data.empty());
}

void test_compute_and_package_matches_host_internal_reference() {
  const DomainData domain_data = make_internal_boundary_domain_data();

  for (const auto& element_id : domain_data.element_ids) {
    INFO("Checking element " << element_id);
    const SetupData setup = make_setup_data(
        element_id, domain_data.domain, domain_data.initial_refinement_levels,
        domain_data.mesh);

    dg::MortarMap<Dim, Mesh<Dim>> neighbor_mesh{};
    for (const auto& [direction, neighbors] : setup.element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        neighbor_mesh.insert(
            {DirectionalId<Dim>{direction, neighbor}, setup.mesh});
      }
    }
    auto mortar_infos =
        evolution::dg::Initialization::detail::mortar_infos<Dim>(setup.element);
    auto [mortar_meshes, mortar_next_temporal_ids,
          normal_covector_and_magnitude] =
        evolution::dg::Initialization::detail::mortars_apply_impl<Dim>(
            setup.element, setup.time_step_id, setup.mesh, neighbor_mesh);
    (void)mortar_next_temporal_ids;
    (void)normal_covector_and_magnitude;

    typename gh::KokkosTags::DeviceVariables<system>::type device_vars =
        copy_to_device(setup.vars);
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

    gh::Actions::InitializeKokkosTags<system>::apply(
        make_not_null(&device_inverse_jacobian),
        make_not_null(&device_inertial_coordinates),
        make_not_null(&device_gamma0), make_not_null(&device_gamma1),
        make_not_null(&device_gamma2),
        make_not_null(&device_face_to_volume_index_map),
        make_not_null(&device_face_unit_normal_covector),
        make_not_null(&device_face_normal_magnitude),
        make_not_null(&device_mortar_data), setup.inverse_jacobian, setup.gamma0,
        setup.gamma1, setup.gamma2, setup.inertial_coords, setup.mesh,
        setup.element, mortar_meshes, mortar_infos);

    typename gh::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
        outgoing_boundary_data{};
    typename gh::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type
        external_boundary_data{};
    external_boundary_conditions_type external_boundary_conditions_by_block{
        domain_data.domain.blocks().size()};
    gh::Actions::ComputeTimeDerivativeKokkos::apply(
        make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
        make_not_null(&external_boundary_data), device_vars,
        device_inverse_jacobian, device_inertial_coordinates, device_gamma0,
        device_gamma1, device_gamma2, device_face_to_volume_index_map,
        device_face_unit_normal_covector, device_face_normal_magnitude,
        device_mortar_data, mortar_meshes, setup.gamma0, setup.gamma1,
        setup.gamma2,
        external_boundary_conditions_by_block, setup.time, setup.mesh,
        setup.element, setup.time_step_id);

    dt_variables_type dt_kokkos{setup.mesh.number_of_grid_points(), 0.0};
    copy_to_host(make_not_null(&dt_kokkos), device_dt);

    dt_variables_type dt_reference{setup.mesh.number_of_grid_points(), 0.0};
    packaged_data_by_mortar_map host_packaged_data_by_mortar{};
    compute_host_reference(setup, make_not_null(&dt_reference),
                           make_not_null(&host_packaged_data_by_mortar));

    const Approx approx = Approx::custom().epsilon(1e-10).scale(1.0);
    CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<spacetime_metric_tag>>(dt_kokkos),
                                 get<::Tags::dt<spacetime_metric_tag>>(dt_reference),
                                 approx);
    CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<pi_tag>>(dt_kokkos),
                                 get<::Tags::dt<pi_tag>>(dt_reference), approx);
    CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<phi_tag>>(dt_kokkos),
                                 get<::Tags::dt<phi_tag>>(dt_reference), approx);

    CHECK(outgoing_boundary_data.size() == setup.element.number_of_neighbors());
    CHECK(external_boundary_data.empty());

    for (const auto& [direction, neighbors] : setup.element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        const DirectionalId<Dim> mortar_id{direction, neighbor};
        CHECK(outgoing_boundary_data.count(mortar_id) == 1);
        CHECK(host_packaged_data_by_mortar.count(mortar_id) == 1);
        const auto& packaged_host = host_packaged_data_by_mortar.at(mortar_id);

        const auto& kokkos_boundary_data = outgoing_boundary_data.at(mortar_id);
        CHECK(kokkos_boundary_data.validity_range == setup.time_step_id);
        CHECK(kokkos_boundary_data.integration_order == 0);

        Variables<dg_package_field_tags> packaged_kokkos{
            mortar_meshes.at(mortar_id).number_of_grid_points(), 0.0};
        copy_to_host(make_not_null(&packaged_kokkos),
                     kokkos_boundary_data.boundary_correction_data);
        CHECK(packaged_kokkos.size() == packaged_host.size());
        CHECK_ITERABLE_CUSTOM_APPROX(
            gsl::make_span(packaged_kokkos.data(), packaged_kokkos.size()),
            gsl::make_span(packaged_host.data(), packaged_host.size()), approx);
      }
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.ComputeTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_volume_terms_match_host_time_derivative();
  test_compute_and_package_matches_host_internal_reference();
}
