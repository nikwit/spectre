// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/Interval.hpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/Tags.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/AnalyticChristoffel.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BoundaryCorrectionHelpers.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeInternalBoundaryTermsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeTimeDerivativeBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/DriverState.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/InitializeBoundaryBatchMetadata.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/LiftInternalBoundaryTermsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/PackageLocalFacesBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/Phase.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsBashforth.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
constexpr size_t Dim = 3;
using system = gh::System<Dim>;
using hardcoded_solution = gh::Solutions::WrappedGr<gr::Solutions::KerrSchild>;

using variables_type = typename system::variables_tag::type;
using dt_variables_tags =
    db::wrap_tags_in<::Tags::dt, typename system::variables_tag::tags_list>;
using dt_variables_type = Variables<dt_variables_tags>;

using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, Dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, Dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, Dim, Frame::Inertial>;
using d_spacetime_metric_tag =
    ::Tags::deriv<spacetime_metric_tag, tmpl::size_t<Dim>, Frame::Inertial>;
using d_pi_tag = ::Tags::deriv<pi_tag, tmpl::size_t<Dim>, Frame::Inertial>;
using d_phi_tag = ::Tags::deriv<phi_tag, tmpl::size_t<Dim>, Frame::Inertial>;
using derivative_variables_tags =
    db::wrap_tags_in<::Tags::deriv, typename system::gradient_variables,
                     tmpl::size_t<Dim>, Frame::Inertial>;

using outgoing_boundary_data_map =
    typename gh::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type;
using incoming_boundary_data_map =
    typename gh::KokkosTags::IncomingBoundaryCorrectionData<Dim>::type;
using external_boundary_data_map =
    typename gh::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type;
using external_boundary_conditions_type =
    typename domain::Tags::ExternalBoundaryConditions<Dim>::type;
using face_sums_by_direction = DirectionMap<Dim, dt_variables_type>;

using packed_topology_type = evolution::Kokkos::PackedTopology<system>;
using packed_geometry_type = evolution::Kokkos::PackedGeometry<system>;
using packed_boundary_metadata_type =
    evolution::Kokkos::PackedBoundaryMetadata<system>;
using packed_evolution_state_type =
    evolution::Kokkos::PackedEvolutionState<system>;
using packed_boundary_scratch_type =
    evolution::Kokkos::PackedBoundaryScratch<system>;
using device_gamma0_type = gh::KokkosTags::DeviceConstraintGamma0::type;
using device_gamma1_type = gh::KokkosTags::DeviceConstraintGamma1::type;
using device_gamma2_type = gh::KokkosTags::DeviceConstraintGamma2::type;

constexpr size_t num_spacetime_metric_components = (Dim + 1) * (Dim + 2) / 2;
constexpr size_t num_phi_components = Dim * num_spacetime_metric_components;
constexpr size_t num_char_speed_components = Dim + 1;
constexpr size_t offset_v_spacetime_metric = 0;
constexpr size_t offset_v_zero =
    offset_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t offset_v_plus = offset_v_zero + num_phi_components;
constexpr size_t offset_v_minus =
    offset_v_plus + num_spacetime_metric_components;
constexpr size_t offset_normal_times_v_plus =
    offset_v_minus + num_spacetime_metric_components;
constexpr size_t offset_normal_times_v_minus =
    offset_normal_times_v_plus + num_phi_components;
constexpr size_t offset_gamma2_v_spacetime_metric =
    offset_normal_times_v_minus + num_phi_components;
constexpr size_t offset_char_speeds =
    offset_gamma2_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t packaged_boundary_data_components =
    offset_char_speeds + num_char_speed_components;

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
  std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> mesh_velocity{};
  std::optional<Scalar<DataVector>> div_mesh_velocity{};
  TimeStepId time_step_id{};
  double time{};

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
  dg::MortarMap<Dim, Mesh<Dim - 1>> mortar_meshes{};
  dg::MortarMap<Dim, evolution::dg::MortarInfo<Dim>> mortar_infos{};
};

struct DomainData {
  Domain<Dim> domain;
  std::vector<std::array<size_t, Dim>> initial_refinement_levels;
  std::vector<std::array<size_t, Dim>> initial_extents;
  std::vector<ElementId<Dim>> element_ids;
  Mesh<Dim> mesh;
  external_boundary_conditions_type external_boundary_conditions;
  domain::FunctionsOfTimeMap functions_of_time;
};

struct HostBoundaryData {
  dt_variables_type volume_dt{};
  outgoing_boundary_data_map outgoing{};
  external_boundary_data_map external{};
};

struct BatchedData {
  packed_topology_type packed_topology{};
  packed_geometry_type packed_geometry{};
  packed_boundary_metadata_type packed_boundary_metadata{};
  packed_evolution_state_type packed_evolution_state{};
  packed_boundary_scratch_type packed_boundary_scratch{};
  device_gamma0_type device_gamma0{};
  device_gamma1_type device_gamma1{};
  device_gamma2_type device_gamma2{};
};

dt_variables_type pack_element_dt_data(
    const packed_topology_type& packed_topology,
    const std::map<ElementId<Dim>, SetupData>& setup_by_id,
    const std::map<ElementId<Dim>, dt_variables_type>& dt_by_element) {
  dt_variables_type packed_dt{packed_topology.total_points, 0.0};
  const size_t num_elements = packed_topology.local_elements.size();
  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& setup = setup_by_id.at(element_id);
    const auto& element_dt = dt_by_element.at(element_id);
    const size_t offset = packed_topology.element_point_offsets_host[e];
    const size_t num_points = setup.mesh.number_of_grid_points();
    for (size_t component = 0;
         component < dt_variables_type::number_of_independent_components;
         ++component) {
      std::copy_n(
          element_dt.data() + component * num_points, num_points,
          packed_dt.data() + component * packed_topology.total_points + offset);
    }
  }
  return packed_dt;
}

Direction<Dim> direction_from_face_id(const size_t face_id) {
  const size_t sliced_dim = face_id / 2;
  const Side side = face_id % 2 == 0 ? Side::Lower : Side::Upper;
  return Direction<Dim>{sliced_dim, side};
}

DomainData make_periodic_domain_data(
    const bool use_time_dependent_map = false) {
  std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
      time_dependence{};
  if (use_time_dependent_map) {
    time_dependence = std::make_unique<
        domain::creators::time_dependence::UniformTranslation<Dim, 0>>(
        0.0, std::array<double, Dim>{{0.2, -0.1, 0.05}});
  }
  const domain::creators::Rectilinear<Dim> domain_creator{
      {{-1.0, -0.4, 0.2}},
      {{1.2, 0.9, 2.3}},
      {{1, 2, 1}},
      {{4, 4, 4}},
      {{true, true, true}},
      {},
      std::move(time_dependence)};
  auto domain = domain_creator.create_domain();
  auto initial_refinement_levels = domain_creator.initial_refinement_levels();
  auto initial_extents = domain_creator.initial_extents();
  auto element_ids = initial_element_ids(initial_refinement_levels);
  auto external_boundary_conditions =
      domain_creator.external_boundary_conditions();
  auto functions_of_time = domain_creator.functions_of_time();
  const Mesh<Dim> mesh{initial_extents[0], Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  return DomainData{std::move(domain),
                    std::move(initial_refinement_levels),
                    std::move(initial_extents),
                    std::move(element_ids),
                    mesh,
                    std::move(external_boundary_conditions),
                    std::move(functions_of_time)};
}

DomainData make_external_domain_data(
    const bool use_time_dependent_map = false) {
  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  std::array<
      std::array<std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>,
                 2>,
      Dim>
      boundary_conditions{};
  for (size_t d = 0; d < Dim; ++d) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      boundary_conditions[d][side_i] =
          std::make_unique<gh::BoundaryConditions::DirichletAnalytic<Dim>>(
              std::make_unique<hardcoded_solution>(1.0, zero_spin,
                                                   zero_center));
    }
  }

  std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
      time_dependence{};
  if (use_time_dependent_map) {
    time_dependence = std::make_unique<
        domain::creators::time_dependence::UniformTranslation<Dim, 0>>(
        0.0, std::array<double, Dim>{{0.0, 0.0, 0.0}});
  }

  const domain::creators::Rectilinear<Dim> domain_creator{
      {{-0.8, -0.6, -0.4}},
      {{0.8, 0.6, 0.4}},
      {{0, 0, 0}},
      {{4, 4, 4}},
      std::move(boundary_conditions),
      {},
      std::move(time_dependence)};
  auto domain = domain_creator.create_domain();
  auto initial_refinement_levels = domain_creator.initial_refinement_levels();
  auto initial_extents = domain_creator.initial_extents();
  auto element_ids = initial_element_ids(initial_refinement_levels);
  auto external_boundary_conditions =
      domain_creator.external_boundary_conditions();
  auto functions_of_time = domain_creator.functions_of_time();
  const Mesh<Dim> mesh{initial_extents[0], Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  return DomainData{std::move(domain),
                    std::move(initial_refinement_levels),
                    std::move(initial_extents),
                    std::move(element_ids),
                    mesh,
                    std::move(external_boundary_conditions),
                    std::move(functions_of_time)};
}

SetupData make_setup_data(
    const ElementId<Dim>& element_id, const Domain<Dim>& domain_object,
    const std::vector<std::array<size_t, Dim>>& initial_refinement_levels,
    const Mesh<Dim>& mesh,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  SetupData setup{};
  setup.element = domain::Initialization::create_initial_element(
      element_id, domain_object.blocks(), initial_refinement_levels);
  setup.mesh = mesh;

  const auto& block = domain_object.blocks()[setup.element.id().block_id()];
  const ElementMap<Dim, Frame::Grid> element_map{setup.element.id(), block};
  const auto logical_coords = logical_coordinates(setup.mesh);
  const auto grid_coords = element_map(logical_coords);
  const auto logical_to_grid_inverse_jacobian =
      element_map.inv_jacobian(logical_coords);

  const Slab slab{0.0, 1.0};
  setup.time_step_id = TimeStepId{true, 0, slab.start()};
  setup.time = 0.2;

  if (block.is_time_dependent()) {
    const auto map_data = block.moving_mesh_grid_to_inertial_map()
                              .coords_frame_velocity_jacobians(
                                  grid_coords, setup.time, functions_of_time);
    setup.inertial_coords = std::get<0>(map_data);
    const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
    setup.mesh_velocity = std::get<3>(map_data);
    setup.div_mesh_velocity =
        Scalar<DataVector>{setup.mesh.number_of_grid_points(), 0.0};
    setup.inverse_jacobian =
        domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                      Frame::Inertial>::type{
            setup.mesh.number_of_grid_points(), 0.0};
    for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
      for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
        for (size_t p = 0; p < setup.mesh.number_of_grid_points(); ++p) {
          setup.inverse_jacobian.get(logical_d, inertial_d)[p] = 0.0;
          for (size_t grid_d = 0; grid_d < Dim; ++grid_d) {
            setup.inverse_jacobian.get(logical_d, inertial_d)[p] +=
                logical_to_grid_inverse_jacobian.get(logical_d, grid_d)[p] *
                grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p];
          }
        }
      }
    }
  } else {
    setup.inertial_coords = tnsr::I<DataVector, Dim, Frame::Inertial>{
        setup.mesh.number_of_grid_points()};
    setup.inverse_jacobian =
        domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                      Frame::Inertial>::type{
            setup.mesh.number_of_grid_points()};
    for (size_t d = 0; d < Dim; ++d) {
      setup.inertial_coords.get(d) = grid_coords.get(d);
    }
    for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
      for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
        setup.inverse_jacobian.get(logical_d, inertial_d) =
            logical_to_grid_inverse_jacobian.get(logical_d, inertial_d);
      }
    }
    setup.mesh_velocity = std::nullopt;
    setup.div_mesh_velocity = std::nullopt;
  }

  const std::array<double, 3> zero_spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> zero_center{{0.0, 0.0, 0.0}};
  const hardcoded_solution solution{1.0, zero_spin, zero_center};
  const auto analytic_vars =
      solution.variables(setup.inertial_coords, setup.time,
                         tmpl::list<spacetime_metric_tag, pi_tag, phi_tag>{});
  const size_t num_points = setup.mesh.number_of_grid_points();
  setup.vars = variables_type{num_points, 0.0};
  get<spacetime_metric_tag>(setup.vars) =
      get<spacetime_metric_tag>(analytic_vars);
  get<pi_tag>(setup.vars) = get<pi_tag>(analytic_vars);
  get<phi_tag>(setup.vars) = get<phi_tag>(analytic_vars);

  setup.gamma0 = gh::Tags::ConstraintGamma0::type{num_points, 0.0};
  setup.gamma1 = gh::Tags::ConstraintGamma1::type{num_points, 0.0};
  setup.gamma2 = gh::Tags::ConstraintGamma2::type{num_points, 0.0};
  get(setup.gamma0) = 0.9 + 0.1 * setup.inertial_coords.get(0) -
                      0.05 * setup.inertial_coords.get(1);
  get(setup.gamma1) = -0.1 + 0.04 * setup.inertial_coords.get(2);
  get(setup.gamma2) = 0.5 + 0.07 * setup.inertial_coords.get(0) +
                      0.02 * setup.inertial_coords.get(1);

  dg::MortarMap<Dim, Mesh<Dim>> neighbor_mesh{};
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    for (const auto& neighbor : neighbors) {
      neighbor_mesh.insert(
          {DirectionalId<Dim>{direction, neighbor}, setup.mesh});
    }
  }
  setup.mortar_infos =
      evolution::dg::Initialization::detail::mortar_infos<Dim>(setup.element);
  auto [mortar_meshes, mortar_next_temporal_ids,
        normal_covector_and_magnitude] =
      evolution::dg::Initialization::detail::mortars_apply_impl<Dim>(
          setup.element, setup.time_step_id, setup.mesh, neighbor_mesh);
  setup.mortar_meshes = std::move(mortar_meshes);
  (void)mortar_next_temporal_ids;
  (void)normal_covector_and_magnitude;

  gh::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&setup.device_inverse_jacobian),
      make_not_null(&setup.device_inertial_coordinates),
      make_not_null(&setup.device_gamma0), make_not_null(&setup.device_gamma1),
      make_not_null(&setup.device_gamma2),
      make_not_null(&setup.device_face_to_volume_index_map),
      make_not_null(&setup.device_face_unit_normal_covector),
      make_not_null(&setup.device_face_normal_magnitude),
      make_not_null(&setup.device_mortar_data), setup.inverse_jacobian,
      setup.gamma0, setup.gamma1, setup.gamma2, setup.inertial_coords,
      setup.mesh, setup.element, setup.mortar_meshes, setup.mortar_infos);
  return setup;
}

HostBoundaryData compute_host_boundary_data(
    const SetupData& setup, const external_boundary_conditions_type&
                                external_boundary_conditions_by_block) {
  auto device_vars = copy_to_device(setup.vars);
  typename gh::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  HostBoundaryData result{};
  gh::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&result.outgoing),
      make_not_null(&result.external), device_vars,
      setup.device_inverse_jacobian, setup.device_inertial_coordinates,
      setup.device_gamma0, setup.device_gamma1, setup.device_gamma2,
      setup.device_face_to_volume_index_map,
      setup.device_face_unit_normal_covector,
      setup.device_face_normal_magnitude, setup.device_mortar_data,
      setup.mortar_meshes, setup.gamma0, setup.gamma1, setup.gamma2,
      external_boundary_conditions_by_block, setup.time, setup.mesh,
      setup.element, setup.time_step_id);
  result.volume_dt = dt_variables_type{setup.mesh.number_of_grid_points(), 0.0};
  copy_to_host(make_not_null(&result.volume_dt), device_dt);
  return result;
}

dt_variables_type compute_host_volume_dt_reference(const SetupData& setup) {
  const size_t num_points = setup.mesh.number_of_grid_points();
  Variables<derivative_variables_tags> host_partial_derivs{num_points};
  partial_derivatives(make_not_null(&host_partial_derivs), setup.vars,
                      setup.mesh, setup.inverse_jacobian);

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
      get<spacetime_metric_tag>(setup.vars), get<pi_tag>(setup.vars),
      get<phi_tag>(setup.vars), setup.gamma0, setup.gamma1, setup.gamma2,
      gauge_condition, setup.mesh, setup.time, setup.inertial_coords,
      setup.inverse_jacobian, setup.mesh_velocity);

  dt_variables_type result{num_points, 0.0};
  get<::Tags::dt<spacetime_metric_tag>>(result) = host_dt_spacetime_metric;
  get<::Tags::dt<pi_tag>>(result) = host_dt_pi;
  get<::Tags::dt<phi_tag>>(result) = host_dt_phi;
  return result;
}

outgoing_boundary_data_map compute_host_outgoing_boundary_data_moving_map(
    const SetupData& setup) {
  outgoing_boundary_data_map result{};
  const auto& spacetime_metric = get<spacetime_metric_tag>(setup.vars);
  const auto& pi = get<pi_tag>(setup.vars);
  const auto& phi = get<phi_tag>(setup.vars);

  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    const size_t sliced_dim = direction.dimension();
    const size_t num_face_points =
        setup.mesh.slice_away(sliced_dim).number_of_grid_points();
    const auto face_to_volume_index =
        direction.side() == Side::Upper
            ? gsl::at(setup.device_face_to_volume_index_map, sliced_dim).second
            : gsl::at(setup.device_face_to_volume_index_map, sliced_dim).first;
    const auto face_to_volume_index_host = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, face_to_volume_index);
    const double outward = direction.side() == Side::Upper ? 1.0 : -1.0;

    gh::KokkosTags::BoundaryCorrectionDataBuffer<Dim> packaged_face_data{
        num_face_points};
    auto packaged_face_data_host =
        Kokkos::create_mirror_view(packaged_face_data.view());

    for (size_t face_index = 0; face_index < num_face_points; ++face_index) {
      const size_t volume_index = face_to_volume_index_host(face_index);

      tnsr::i<double, Dim, Frame::Inertial> unnormalized_normal_covector{};
      for (size_t d = 0; d < Dim; ++d) {
        unnormalized_normal_covector.get(d) =
            outward * setup.inverse_jacobian.get(sliced_dim, d)[volume_index];
      }

      tnsr::I<double, Dim, Frame::Inertial> mesh_velocity_at_point{};
      if (setup.mesh_velocity.has_value()) {
        for (size_t d = 0; d < Dim; ++d) {
          mesh_velocity_at_point.get(d) =
              setup.mesh_velocity->get(d)[volume_index];
        }
      } else {
        for (size_t d = 0; d < Dim; ++d) {
          mesh_velocity_at_point.get(d) = 0.0;
        }
      }

      tnsr::aa<double, Dim, Frame::Inertial> v_spacetime_metric{};
      tnsr::iaa<double, Dim, Frame::Inertial> v_zero{};
      tnsr::aa<double, Dim, Frame::Inertial> v_plus{};
      tnsr::aa<double, Dim, Frame::Inertial> v_minus{};
      tnsr::iaa<double, Dim, Frame::Inertial> normal_times_v_plus{};
      tnsr::iaa<double, Dim, Frame::Inertial> normal_times_v_minus{};
      tnsr::aa<double, Dim, Frame::Inertial> gamma2_v_spacetime_metric{};
      tnsr::a<double, Dim, Frame::Inertial> char_speeds{};

      tnsr::aa<double, Dim, Frame::Inertial> spacetime_metric_at_point{};
      tnsr::aa<double, Dim, Frame::Inertial> pi_at_point{};
      tnsr::iaa<double, Dim, Frame::Inertial> phi_at_point{};
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          spacetime_metric_at_point.get(a, b) =
              spacetime_metric.get(a, b)[volume_index];
          pi_at_point.get(a, b) = pi.get(a, b)[volume_index];
          for (size_t d = 0; d < Dim; ++d) {
            phi_at_point.get(d, a, b) = phi.get(d, a, b)[volume_index];
          }
        }
      }

      gh::Actions::detail::compute_packaged_boundary_data_at_point<Dim>(
          make_not_null(&v_spacetime_metric), make_not_null(&v_zero),
          make_not_null(&v_plus), make_not_null(&v_minus),
          make_not_null(&normal_times_v_plus),
          make_not_null(&normal_times_v_minus),
          make_not_null(&gamma2_v_spacetime_metric),
          make_not_null(&char_speeds), spacetime_metric_at_point, pi_at_point,
          phi_at_point, get(setup.gamma1)[volume_index],
          get(setup.gamma2)[volume_index], unnormalized_normal_covector,
          mesh_velocity_at_point);

      size_t component_offset = 0;
      for (size_t c = 0; c < v_spacetime_metric.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) =
            v_spacetime_metric[c];
      }
      component_offset += v_spacetime_metric.size();
      for (size_t c = 0; c < v_zero.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) = v_zero[c];
      }
      component_offset += v_zero.size();
      for (size_t c = 0; c < v_plus.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) = v_plus[c];
      }
      component_offset += v_plus.size();
      for (size_t c = 0; c < v_minus.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) = v_minus[c];
      }
      component_offset += v_minus.size();
      for (size_t c = 0; c < normal_times_v_plus.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) =
            normal_times_v_plus[c];
      }
      component_offset += normal_times_v_plus.size();
      for (size_t c = 0; c < normal_times_v_minus.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) =
            normal_times_v_minus[c];
      }
      component_offset += normal_times_v_minus.size();
      for (size_t c = 0; c < gamma2_v_spacetime_metric.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) =
            gamma2_v_spacetime_metric[c];
      }
      component_offset += gamma2_v_spacetime_metric.size();
      for (size_t c = 0; c < char_speeds.size(); ++c) {
        packaged_face_data_host(face_index, component_offset + c) =
            char_speeds[c];
      }
    }

    Kokkos::deep_copy(packaged_face_data.view(), packaged_face_data_host);

    for (const auto& neighbor : neighbors) {
      const DirectionalId<Dim> mortar_id{direction, neighbor};
      gh::KokkosTags::BoundaryCorrectionData<Dim> boundary_data{};
      boundary_data.boundary_correction_data = packaged_face_data;
      result.insert_or_assign(mortar_id, std::move(boundary_data));
    }
  }
  return result;
}

incoming_boundary_data_map build_incoming_boundary_data(
    const SetupData& receiver_setup,
    const std::map<ElementId<Dim>, outgoing_boundary_data_map>&
        outgoing_by_element) {
  incoming_boundary_data_map incoming_boundary_data{};
  for (const auto& [direction, neighbors] :
       receiver_setup.element.neighbors()) {
    REQUIRE(neighbors.size() == 1);
    const auto& sender_id = *neighbors.begin();
    const auto& orientation = neighbors.orientation(sender_id);
    REQUIRE(orientation.is_aligned());

    const auto& sender_outgoing = outgoing_by_element.at(sender_id);
    const DirectionalId<Dim> sender_mortar_id{direction.opposite(),
                                              receiver_setup.element.id()};
    REQUIRE(sender_outgoing.count(sender_mortar_id) == 1);

    auto data_for_receiver = sender_outgoing.at(sender_mortar_id);
    incoming_boundary_data.insert_or_assign(
        DirectionalId<Dim>{direction, sender_id}, std::move(data_for_receiver));
  }
  return incoming_boundary_data;
}

template <typename ViewType>
void unpack_packaged_boundary_data(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        v_spacetime_metric,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> v_zero,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> v_plus,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> v_minus,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        normal_times_v_plus,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        normal_times_v_minus,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        gamma2_v_spacetime_metric,
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*> char_speeds,
    const ViewType& boundary_data_view, const size_t num_points) {
  for (size_t p = 0; p < num_points; ++p) {
    size_t component_index = offset_v_spacetime_metric;
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        v_spacetime_metric->get(a, b)[p] =
            boundary_data_view(p, component_index++);
      }
    }
    component_index = offset_v_zero;
    for (size_t d = 0; d < Dim; ++d) {
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          v_zero->get(d, a, b)[p] = boundary_data_view(p, component_index++);
        }
      }
    }
    component_index = offset_v_plus;
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        v_plus->get(a, b)[p] = boundary_data_view(p, component_index++);
      }
    }
    component_index = offset_v_minus;
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        v_minus->get(a, b)[p] = boundary_data_view(p, component_index++);
      }
    }
    component_index = offset_normal_times_v_plus;
    for (size_t d = 0; d < Dim; ++d) {
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          normal_times_v_plus->get(d, a, b)[p] =
              boundary_data_view(p, component_index++);
        }
      }
    }
    component_index = offset_normal_times_v_minus;
    for (size_t d = 0; d < Dim; ++d) {
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          normal_times_v_minus->get(d, a, b)[p] =
              boundary_data_view(p, component_index++);
        }
      }
    }
    component_index = offset_gamma2_v_spacetime_metric;
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        gamma2_v_spacetime_metric->get(a, b)[p] =
            boundary_data_view(p, component_index++);
      }
    }
    component_index = offset_char_speeds;
    for (size_t a = 0; a < Dim + 1; ++a) {
      char_speeds->get(a)[p] = boundary_data_view(p, component_index++);
    }
  }
}

double metric_normal_magnitude_at_volume_index(
    const variables_type& vars, const size_t volume_index,
    const tnsr::i<double, Dim, Frame::Inertial>& unnormalized_normal_covector) {
  const auto& spacetime_metric = get<spacetime_metric_tag>(vars);
  const double g00 = spacetime_metric.get(1, 1)[volume_index];
  const double g01 = spacetime_metric.get(1, 2)[volume_index];
  const double g02 = spacetime_metric.get(1, 3)[volume_index];
  const double g11 = spacetime_metric.get(2, 2)[volume_index];
  const double g12 = spacetime_metric.get(2, 3)[volume_index];
  const double g22 = spacetime_metric.get(3, 3)[volume_index];
  const double det_spatial_metric = g00 * (g11 * g22 - g12 * g12) -
                                    g01 * (g01 * g22 - g12 * g02) +
                                    g02 * (g01 * g12 - g11 * g02);
  const double inv_det = 1.0 / det_spatial_metric;

  tnsr::II<double, Dim, Frame::Inertial> inverse_spatial_metric{};
  inverse_spatial_metric.get(0, 0) = (g11 * g22 - g12 * g12) * inv_det;
  inverse_spatial_metric.get(0, 1) = (g02 * g12 - g01 * g22) * inv_det;
  inverse_spatial_metric.get(0, 2) = (g01 * g12 - g02 * g11) * inv_det;
  inverse_spatial_metric.get(1, 1) = (g00 * g22 - g02 * g02) * inv_det;
  inverse_spatial_metric.get(1, 2) = (g02 * g01 - g00 * g12) * inv_det;
  inverse_spatial_metric.get(2, 2) = (g00 * g11 - g01 * g01) * inv_det;

  tnsr::I<double, Dim, Frame::Inertial> unnormalized_normal_vector{};
  double normal_magnitude_squared = 0.0;
  for (size_t i = 0; i < Dim; ++i) {
    unnormalized_normal_vector.get(i) = 0.0;
    for (size_t j = 0; j < Dim; ++j) {
      unnormalized_normal_vector.get(i) += inverse_spatial_metric.get(i, j) *
                                           unnormalized_normal_covector.get(j);
    }
    normal_magnitude_squared +=
        unnormalized_normal_vector.get(i) * unnormalized_normal_covector.get(i);
  }
  return sqrt(normal_magnitude_squared);
}

face_sums_by_direction compute_host_internal_face_sums(
    const SetupData& setup,
    const outgoing_boundary_data_map& local_outgoing_boundary_data,
    const incoming_boundary_data_map& incoming_boundary_data) {
  face_sums_by_direction result{};
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    const Mesh<Dim - 1> face_mesh =
        setup.mesh.slice_away(direction.dimension());
    dt_variables_type dt_boundary_correction_on_face_sum{
        face_mesh.number_of_grid_points(), 0.0};

    for (const auto& neighbor : neighbors) {
      const DirectionalId<Dim> mortar_id{direction, neighbor};
      REQUIRE(local_outgoing_boundary_data.count(mortar_id) == 1);
      REQUIRE(incoming_boundary_data.count(mortar_id) == 1);
      const auto& local_boundary_data =
          local_outgoing_boundary_data.at(mortar_id);
      const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);
      const auto& mortar_mesh = setup.mortar_meshes.at(mortar_id);
      const auto& mortar_data = setup.device_mortar_data.at(mortar_id);
      const auto local_packaged_data_view = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{},
          local_boundary_data.boundary_correction_data.view());
      const auto remote_packaged_data_view =
          Kokkos::create_mirror_view_and_copy(
              Kokkos::HostSpace{},
              remote_boundary_data.boundary_correction_data.view());
      REQUIRE(local_packaged_data_view.extent(1) ==
              packaged_boundary_data_components);
      REQUIRE(remote_packaged_data_view.extent(1) ==
              packaged_boundary_data_components);

      dt_variables_type dt_boundary_correction_on_mortar{
          mortar_mesh.number_of_grid_points(), 0.0};
      auto& dt_spacetime_metric_on_mortar =
          get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              dt_boundary_correction_on_mortar);
      auto& dt_pi_on_mortar = get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(
          dt_boundary_correction_on_mortar);
      auto& dt_phi_on_mortar = get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(
          dt_boundary_correction_on_mortar);
      tnsr::aa<DataVector, Dim, Frame::Inertial> local_v_spacetime_metric{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> local_v_zero{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial> local_v_plus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial> local_v_minus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> local_normal_times_v_plus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> local_normal_times_v_minus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial>
          local_gamma2_v_spacetime_metric{mortar_mesh.number_of_grid_points(),
                                          0.0};
      tnsr::a<DataVector, Dim, Frame::Inertial> local_char_speeds{
          mortar_mesh.number_of_grid_points(), 0.0};

      tnsr::aa<DataVector, Dim, Frame::Inertial> remote_v_spacetime_metric{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> remote_v_zero{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial> remote_v_plus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial> remote_v_minus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> remote_normal_times_v_plus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::iaa<DataVector, Dim, Frame::Inertial> remote_normal_times_v_minus{
          mortar_mesh.number_of_grid_points(), 0.0};
      tnsr::aa<DataVector, Dim, Frame::Inertial>
          remote_gamma2_v_spacetime_metric{mortar_mesh.number_of_grid_points(),
                                           0.0};
      tnsr::a<DataVector, Dim, Frame::Inertial> remote_char_speeds{
          mortar_mesh.number_of_grid_points(), 0.0};

      unpack_packaged_boundary_data(
          make_not_null(&local_v_spacetime_metric),
          make_not_null(&local_v_zero), make_not_null(&local_v_plus),
          make_not_null(&local_v_minus),
          make_not_null(&local_normal_times_v_plus),
          make_not_null(&local_normal_times_v_minus),
          make_not_null(&local_gamma2_v_spacetime_metric),
          make_not_null(&local_char_speeds), local_packaged_data_view,
          mortar_mesh.number_of_grid_points());
      unpack_packaged_boundary_data(
          make_not_null(&remote_v_spacetime_metric),
          make_not_null(&remote_v_zero), make_not_null(&remote_v_plus),
          make_not_null(&remote_v_minus),
          make_not_null(&remote_normal_times_v_plus),
          make_not_null(&remote_normal_times_v_minus),
          make_not_null(&remote_gamma2_v_spacetime_metric),
          make_not_null(&remote_char_speeds), remote_packaged_data_view,
          mortar_mesh.number_of_grid_points());

      const gh::BoundaryCorrections::UpwindPenalty<Dim> correction{};
      correction.dg_boundary_terms(
          make_not_null(&dt_spacetime_metric_on_mortar),
          make_not_null(&dt_pi_on_mortar), make_not_null(&dt_phi_on_mortar),
          local_v_spacetime_metric, local_v_zero, local_v_plus, local_v_minus,
          local_normal_times_v_plus, local_normal_times_v_minus,
          local_gamma2_v_spacetime_metric, local_char_speeds,
          remote_v_spacetime_metric, remote_v_zero, remote_v_plus,
          remote_v_minus, remote_normal_times_v_plus,
          remote_normal_times_v_minus, remote_gamma2_v_spacetime_metric,
          remote_char_speeds, ::dg::Formulation::StrongInertial);

      dt_variables_type dt_boundary_correction_on_face{
          face_mesh.number_of_grid_points(), 0.0};
      if (mortar_data.needs_projection) {
        ::dg::project_from_mortar(
            make_not_null(&dt_boundary_correction_on_face),
            dt_boundary_correction_on_mortar, face_mesh, mortar_mesh,
            mortar_data.mortar_size);
      } else {
        dt_boundary_correction_on_face = dt_boundary_correction_on_mortar;
      }

      auto& dt_spacetime_metric_face_sum =
          get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              dt_boundary_correction_on_face_sum);
      auto& dt_pi_face_sum = get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(
          dt_boundary_correction_on_face_sum);
      auto& dt_phi_face_sum = get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(
          dt_boundary_correction_on_face_sum);
      const auto& dt_spacetime_metric_face =
          get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              dt_boundary_correction_on_face);
      const auto& dt_pi_face = get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(
          dt_boundary_correction_on_face);
      const auto& dt_phi_face = get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(
          dt_boundary_correction_on_face);
      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          dt_spacetime_metric_face_sum.get(a, b) +=
              dt_spacetime_metric_face.get(a, b);
          dt_pi_face_sum.get(a, b) += dt_pi_face.get(a, b);
          for (size_t d = 0; d < Dim; ++d) {
            dt_phi_face_sum.get(d, a, b) += dt_phi_face.get(d, a, b);
          }
        }
      }
    }
    result.insert_or_assign(direction,
                            std::move(dt_boundary_correction_on_face_sum));
  }
  CHECK(setup.element.external_boundaries().empty());
  return result;
}

void lift_host_internal_face_sums(
    const gsl::not_null<dt_variables_type*> host_dt_reference,
    const SetupData& setup, const face_sums_by_direction& face_sums) {
  for (const auto& [direction, dt_boundary_correction_on_face_sum] :
       face_sums) {
    const Mesh<Dim - 1> face_mesh =
        setup.mesh.slice_away(direction.dimension());
    const size_t sliced_dim = direction.dimension();
    const size_t extent_perpendicular_to_boundary =
        setup.mesh.extents(sliced_dim);
    const double lift_prefactor =
        -0.5 * static_cast<double>(extent_perpendicular_to_boundary *
                                   (extent_perpendicular_to_boundary - 1));

    const auto host_face_to_volume_index =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_to_volume_index_map, sliced_dim)
                      .second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_to_volume_index_map, sliced_dim)
                      .first);
    const auto host_face_unit_normal_covector =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_unit_normal_covector, sliced_dim)
                      .second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_unit_normal_covector, sliced_dim)
                      .first);
    const auto host_face_normal_magnitude =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_normal_magnitude, sliced_dim)
                      .second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_normal_magnitude, sliced_dim)
                      .first);

    for (size_t face_index = 0; face_index < face_mesh.number_of_grid_points();
         ++face_index) {
      const size_t volume_index = host_face_to_volume_index(face_index);

      tnsr::i<double, Dim, Frame::Inertial> unnormalized_normal_covector{};
      for (size_t d = 0; d < Dim; ++d) {
        unnormalized_normal_covector.get(d) =
            host_face_unit_normal_covector(face_index, d) *
            host_face_normal_magnitude(face_index);
      }
      const double normal_magnitude = metric_normal_magnitude_at_volume_index(
          setup.vars, volume_index, unnormalized_normal_covector);
      const double lifted_factor = lift_prefactor * normal_magnitude;

      for (size_t a = 0; a < Dim + 1; ++a) {
        for (size_t b = a; b < Dim + 1; ++b) {
          get<::Tags::dt<spacetime_metric_tag>>(*host_dt_reference)
              .get(a, b)[volume_index] +=
              lifted_factor * get<::Tags::dt<spacetime_metric_tag>>(
                                  dt_boundary_correction_on_face_sum)
                                  .get(a, b)[face_index];
          get<::Tags::dt<pi_tag>>(*host_dt_reference).get(a, b)[volume_index] +=
              lifted_factor *
              get<::Tags::dt<pi_tag>>(dt_boundary_correction_on_face_sum)
                  .get(a, b)[face_index];
          for (size_t d = 0; d < Dim; ++d) {
            get<::Tags::dt<phi_tag>>(*host_dt_reference)
                .get(d, a, b)[volume_index] +=
                lifted_factor *
                get<::Tags::dt<phi_tag>>(dt_boundary_correction_on_face_sum)
                    .get(d, a, b)[face_index];
          }
        }
      }
    }
  }
}

BatchedData initialize_batched_data(
    const DomainData& domain_data,
    const std::map<ElementId<Dim>, SetupData>& setup_by_id) {
  BatchedData batched{};
  gh::Batched::Initialization::DriverState::apply(
      make_not_null(&batched.packed_topology), domain_data.domain,
      domain_data.initial_refinement_levels, domain_data.initial_extents,
      domain_data.mesh.quadrature(0));

  const size_t total_points = batched.packed_topology.total_points;
  const size_t points_per_element = batched.packed_topology.points_per_element;
  const size_t num_elements = batched.packed_topology.local_elements.size();

  variables_type host_vars{total_points, 0.0};
  dt_variables_type host_dt_vars{total_points, 0.0};
  gh::Tags::ConstraintGamma0::type host_gamma0{total_points};
  gh::Tags::ConstraintGamma1::type host_gamma1{total_points};
  gh::Tags::ConstraintGamma2::type host_gamma2{total_points};
  std::array<DataVector, Dim> host_inertial_coords{
      {DataVector{total_points, 0.0}, DataVector{total_points, 0.0},
       DataVector{total_points, 0.0}}};
  std::array<DataVector, Dim> host_mesh_velocity{
      {DataVector{total_points, 0.0}, DataVector{total_points, 0.0},
       DataVector{total_points, 0.0}}};

  batched.packed_geometry.element_inverse_jacobian_device =
      packed_geometry_type::device_inverse_jacobian_type{
          "TestGhBatchedBoundaryElementInverseJacobian", num_elements,
          points_per_element, Dim * Dim};
  auto host_inverse_jacobian = Kokkos::create_mirror_view(
      batched.packed_geometry.element_inverse_jacobian_device);

  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = batched.packed_topology.local_element_ids[e];
    const auto& setup = setup_by_id.at(element_id);
    const size_t offset = batched.packed_topology.element_point_offsets_host[e];
    const size_t num_points = setup.mesh.number_of_grid_points();
    REQUIRE(num_points == points_per_element);

    for (size_t component = 0;
         component < variables_type::number_of_independent_components;
         ++component) {
      std::copy_n(setup.vars.data() + component * num_points, num_points,
                  host_vars.data() + component * total_points + offset);
    }

    for (size_t p = 0; p < num_points; ++p) {
      get(host_gamma0)[offset + p] = get(setup.gamma0)[p];
      get(host_gamma1)[offset + p] = get(setup.gamma1)[p];
      get(host_gamma2)[offset + p] = get(setup.gamma2)[p];
      for (size_t d = 0; d < Dim; ++d) {
        host_inertial_coords[d][offset + p] = setup.inertial_coords.get(d)[p];
        host_mesh_velocity[d][offset + p] = setup.mesh_velocity.has_value()
                                                ? setup.mesh_velocity->get(d)[p]
                                                : 0.0;
      }
    }

    for (size_t logical_d = 0; logical_d < Dim; ++logical_d) {
      for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
        for (size_t p = 0; p < num_points; ++p) {
          host_inverse_jacobian(e, p, logical_d * Dim + inertial_d) =
              setup.inverse_jacobian.get(logical_d, inertial_d)[p];
        }
      }
    }
  }

  for (size_t d = 0; d < Dim; ++d) {
    batched.packed_geometry.inertial_coordinates_device.get(d) =
        Kokkos::View<double*>("TestGhBatchedBoundaryInertialCoordinates",
                              total_points);
    auto host_coords_component = Kokkos::create_mirror_view(
        batched.packed_geometry.inertial_coordinates_device.get(d));
    for (size_t s = 0; s < total_points; ++s) {
      host_coords_component(s) = host_inertial_coords[d][s];
    }
    Kokkos::deep_copy(
        batched.packed_geometry.inertial_coordinates_device.get(d),
        host_coords_component);
    batched.packed_geometry.mesh_velocity_device.get(d) = Kokkos::View<double*>(
        "TestGhBatchedBoundaryMeshVelocity", total_points);
    auto host_mesh_velocity_component = Kokkos::create_mirror_view(
        batched.packed_geometry.mesh_velocity_device.get(d));
    for (size_t s = 0; s < total_points; ++s) {
      host_mesh_velocity_component(s) = host_mesh_velocity[d][s];
    }
    Kokkos::deep_copy(batched.packed_geometry.mesh_velocity_device.get(d),
                      host_mesh_velocity_component);
  }
  batched.packed_geometry.inertial_coordinates_host = host_inertial_coords;
  batched.packed_geometry.mesh_velocity_host = host_mesh_velocity;
  Kokkos::deep_copy(batched.packed_geometry.element_inverse_jacobian_device,
                    host_inverse_jacobian);

  batched.packed_evolution_state.device_variables = copy_to_device(host_vars);
  batched.packed_evolution_state.device_dt_variables =
      copy_to_device(host_dt_vars);

  using gamma0_space = typename device_gamma0_type::value_type::memory_space;
  using gamma1_space = typename device_gamma1_type::value_type::memory_space;
  using gamma2_space = typename device_gamma2_type::value_type::memory_space;
  batched.device_gamma0 =
      copy_to_device(host_gamma0, tmpl::type_<gamma0_space>{});
  batched.device_gamma1 =
      copy_to_device(host_gamma1, tmpl::type_<gamma1_space>{});
  batched.device_gamma2 =
      copy_to_device(host_gamma2, tmpl::type_<gamma2_space>{});

  gh::Actions::InitializeBoundaryBatchMetadata::apply(
      make_not_null(&batched.packed_boundary_metadata),
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology);

  return batched;
}

void test_periodic_boundary_actions_batched_match_host(
    const bool use_time_dependent_map = false) {
  CAPTURE(use_time_dependent_map);
  const DomainData domain_data =
      make_periodic_domain_data(use_time_dependent_map);

  std::map<ElementId<Dim>, SetupData> setup_by_id{};
  for (const auto& element_id : domain_data.element_ids) {
    setup_by_id.insert_or_assign(
        element_id,
        make_setup_data(element_id, domain_data.domain,
                        domain_data.initial_refinement_levels, domain_data.mesh,
                        domain_data.functions_of_time));
  }

  std::map<ElementId<Dim>, outgoing_boundary_data_map> outgoing_by_element{};
  std::map<ElementId<Dim>, dt_variables_type> host_volume_dt_by_element{};
  for (const auto& [element_id, setup] : setup_by_id) {
    auto host_boundary_data = compute_host_boundary_data(
        setup, domain_data.external_boundary_conditions);
    CHECK(host_boundary_data.external.empty());
    if (use_time_dependent_map) {
      host_volume_dt_by_element.insert_or_assign(
          element_id, compute_host_volume_dt_reference(setup));
      outgoing_by_element.insert_or_assign(
          element_id, compute_host_outgoing_boundary_data_moving_map(setup));
    } else {
      host_volume_dt_by_element.insert_or_assign(element_id,
                                                 host_boundary_data.volume_dt);
      outgoing_by_element.insert_or_assign(element_id,
                                           host_boundary_data.outgoing);
    }
  }

  std::map<ElementId<Dim>, incoming_boundary_data_map> incoming_by_element{};
  for (const auto& [element_id, setup] : setup_by_id) {
    incoming_by_element.insert_or_assign(
        element_id, build_incoming_boundary_data(setup, outgoing_by_element));
  }

  std::map<ElementId<Dim>, face_sums_by_direction> host_face_sums_by_element{};
  for (const auto& [element_id, setup] : setup_by_id) {
    host_face_sums_by_element.insert_or_assign(
        element_id, compute_host_internal_face_sums(
                        setup, outgoing_by_element.at(element_id),
                        incoming_by_element.at(element_id)));
  }

  BatchedData batched = initialize_batched_data(domain_data, setup_by_id);
  if (not use_time_dependent_map) {
    for (size_t d = 0; d < Dim; ++d) {
      const auto mesh_velocity_view = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{},
          batched.packed_geometry.mesh_velocity_device.get(d));
      for (size_t s = 0; s < batched.packed_topology.total_points; ++s) {
        CHECK(std::abs(mesh_velocity_view(s)) < 1.0e-14);
      }
    }
  }
  gh::Actions::ComputeTimeDerivativeBatched::apply(
      make_not_null(&batched.packed_evolution_state),
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology,
      batched.packed_geometry, batched.device_gamma0, batched.device_gamma1,
      batched.device_gamma2);
  dt_variables_type dt_volume_from_batched{batched.packed_topology.total_points,
                                           0.0};
  copy_to_host(make_not_null(&dt_volume_from_batched),
               batched.packed_evolution_state.device_dt_variables);
  gh::Actions::PackageLocalFacesBatched::apply(
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology,
      batched.packed_geometry, batched.packed_evolution_state,
      batched.device_gamma0, batched.device_gamma1, batched.device_gamma2);

  const Approx approx = Approx::custom().epsilon(1e-8).scale(1.0);
  const size_t num_face_points =
      domain_data.mesh.slice_away(0).number_of_grid_points();
  const size_t num_elements = batched.packed_topology.local_elements.size();

  for (size_t face_id = 0; face_id < 2 * Dim; ++face_id) {
    const Direction<Dim> direction = direction_from_face_id(face_id);
    const auto batched_face_data_view = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{}, batched.packed_boundary_scratch
                                 .packaged_face_data_for_all_elements[face_id]
                                 .view());

    for (size_t e = 0; e < num_elements; ++e) {
      const auto& element = batched.packed_topology.local_elements[e];
      REQUIRE(element.neighbors().count(direction) == 1);
      const auto& neighbors = element.neighbors().at(direction);
      REQUIRE(neighbors.size() == 1);
      const auto neighbor = *neighbors.begin();
      const DirectionalId<Dim> mortar_id{direction, neighbor};
      const auto& host_outgoing =
          outgoing_by_element.at(element.id()).at(mortar_id);
      const auto host_outgoing_view = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{}, host_outgoing.boundary_correction_data.view());
      REQUIRE(host_outgoing_view.extent(0) == num_face_points);
      REQUIRE(host_outgoing_view.extent(1) ==
              packaged_boundary_data_components);

      INFO("Comparing packaged face data for element "
           << element.id() << " on face_id " << face_id);
      for (size_t p = 0; p < num_face_points; ++p) {
        for (size_t c = 0; c < packaged_boundary_data_components; ++c) {
          CHECK(batched_face_data_view(e * num_face_points + p, c) ==
                approx(host_outgoing_view(p, c)));
        }
      }
    }
  }

  gh::Actions::ComputeInternalBoundaryTermsBatched::apply(
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology,
      batched.packed_boundary_metadata);

  for (size_t face_id = 0; face_id < 2 * Dim; ++face_id) {
    const Direction<Dim> direction = direction_from_face_id(face_id);
    const auto dt_face_sum_view = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{},
        batched.packed_boundary_scratch
            .internal_boundary_terms_for_all_elements[face_id]
            .view());
    const size_t num_components = dt_face_sum_view.extent(1);
    for (size_t e = 0; e < num_elements; ++e) {
      const auto& element_id = batched.packed_topology.local_element_ids[e];
      const auto& host_face_sum =
          host_face_sums_by_element.at(element_id).at(direction);
      REQUIRE(host_face_sum.number_of_grid_points() == num_face_points);
      REQUIRE(dt_variables_type::number_of_independent_components ==
              num_components);
      INFO("Comparing internal boundary terms for element "
           << element_id << " on face_id " << face_id);
      for (size_t p = 0; p < num_face_points; ++p) {
        for (size_t c = 0; c < num_components; ++c) {
          const double host_value =
              host_face_sum.data()[c * num_face_points + p];
          CHECK(dt_face_sum_view(e * num_face_points + p, c) ==
                approx(host_value));
        }
      }
    }
  }

  gh::Actions::LiftInternalBoundaryTermsBatched::apply(
      make_not_null(&batched.packed_evolution_state), batched.packed_topology,
      batched.packed_geometry, batched.packed_boundary_scratch);

  dt_variables_type dt_from_batched{batched.packed_topology.total_points, 0.0};
  copy_to_host(make_not_null(&dt_from_batched),
               batched.packed_evolution_state.device_dt_variables);
  dt_variables_type dt_boundary_from_batched = dt_from_batched;
  dt_boundary_from_batched -= dt_volume_from_batched;

  dt_variables_type dt_host_reference{batched.packed_topology.total_points,
                                      0.0};
  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = batched.packed_topology.local_element_ids[e];
    const auto& setup = setup_by_id.at(element_id);
    dt_variables_type element_dt{setup.mesh.number_of_grid_points(), 0.0};
    lift_host_internal_face_sums(make_not_null(&element_dt), setup,
                                 host_face_sums_by_element.at(element_id));
    const size_t offset = batched.packed_topology.element_point_offsets_host[e];
    const size_t num_points = setup.mesh.number_of_grid_points();
    for (size_t component = 0;
         component < dt_variables_type::number_of_independent_components;
         ++component) {
      std::copy_n(element_dt.data() + component * num_points, num_points,
                  dt_host_reference.data() +
                      component * batched.packed_topology.total_points +
                      offset);
    }
  }

  CHECK_ITERABLE_CUSTOM_APPROX(
      get<::Tags::dt<spacetime_metric_tag>>(dt_boundary_from_batched),
      get<::Tags::dt<spacetime_metric_tag>>(dt_host_reference), approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get<::Tags::dt<pi_tag>>(dt_boundary_from_batched),
      get<::Tags::dt<pi_tag>>(dt_host_reference), approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get<::Tags::dt<phi_tag>>(dt_boundary_from_batched),
      get<::Tags::dt<phi_tag>>(dt_host_reference), approx);

  const auto dt_host_volume_reference = pack_element_dt_data(
      batched.packed_topology, setup_by_id, host_volume_dt_by_element);
  dt_variables_type dt_host_full_reference = dt_host_reference;
  dt_host_full_reference += dt_host_volume_reference;
  CHECK_ITERABLE_CUSTOM_APPROX(
      get<::Tags::dt<spacetime_metric_tag>>(dt_from_batched),
      get<::Tags::dt<spacetime_metric_tag>>(dt_host_full_reference), approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<pi_tag>>(dt_from_batched),
                               get<::Tags::dt<pi_tag>>(dt_host_full_reference),
                               approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<phi_tag>>(dt_from_batched),
                               get<::Tags::dt<phi_tag>>(dt_host_full_reference),
                               approx);
}

void test_external_boundary_apply_batched_matches_host(
    const bool use_time_dependent_map = false) {
  const DomainData domain_data =
      make_external_domain_data(use_time_dependent_map);
  REQUIRE(domain_data.element_ids.size() == 1);

  const auto element_id = domain_data.element_ids[0];
  const SetupData setup = make_setup_data(
      element_id, domain_data.domain, domain_data.initial_refinement_levels,
      domain_data.mesh, domain_data.functions_of_time);

  std::map<ElementId<Dim>, SetupData> setup_by_id{};
  setup_by_id.insert_or_assign(element_id, setup);
  BatchedData batched = initialize_batched_data(domain_data, setup_by_id);
  gh::Actions::ComputeTimeDerivativeBatched::apply(
      make_not_null(&batched.packed_evolution_state),
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology,
      batched.packed_geometry, batched.device_gamma0, batched.device_gamma1,
      batched.device_gamma2);

  gh::Actions::PackageLocalFacesBatched::apply(
      make_not_null(&batched.packed_boundary_scratch), batched.packed_topology,
      batched.packed_geometry, batched.packed_evolution_state,
      batched.device_gamma0, batched.device_gamma1, batched.device_gamma2);

  dt_variables_type dt_before_external{batched.packed_topology.total_points,
                                       0.0};
  copy_to_host(make_not_null(&dt_before_external),
               batched.packed_evolution_state.device_dt_variables);

  gh::Actions::ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched::apply(
      make_not_null(&batched.packed_evolution_state), batched.packed_topology,
      batched.packed_geometry, batched.packed_boundary_metadata,
      batched.device_gamma0, batched.device_gamma1, batched.device_gamma2,
      batched.packed_boundary_scratch);

  dt_variables_type dt_from_batched{batched.packed_topology.total_points, 0.0};
  copy_to_host(make_not_null(&dt_from_batched),
               batched.packed_evolution_state.device_dt_variables);

  size_t changed_entries = 0;
  for (size_t i = 0; i < dt_from_batched.size(); ++i) {
    if (std::abs(dt_from_batched.data()[i] - dt_before_external.data()[i]) >
        1.0e-12) {
      ++changed_entries;
    }
    CHECK(std::isfinite(dt_from_batched.data()[i]));
  }
  CHECK(changed_entries > 0);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.BoundaryActionsBatched",
    "[Unit][Evolution]") {
  test_periodic_boundary_actions_batched_match_host(false);
  test_periodic_boundary_actions_batched_match_host(true);
  test_external_boundary_apply_batched_matches_host(false);
  test_external_boundary_apply_batched_matches_host(true);
}
