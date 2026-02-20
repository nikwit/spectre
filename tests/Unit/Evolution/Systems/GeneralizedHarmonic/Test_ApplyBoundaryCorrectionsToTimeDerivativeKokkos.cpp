// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
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
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Time/Slab.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
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

using outgoing_boundary_data_map =
    typename gh::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type;
using incoming_boundary_data_map =
    typename gh::KokkosTags::IncomingBoundaryCorrectionData<Dim>::type;
using external_boundary_data_map =
    typename gh::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type;
using external_boundary_conditions_type =
    typename domain::Tags::ExternalBoundaryConditions<Dim>::type;

constexpr size_t num_spacetime_metric_components = (Dim + 1) * (Dim + 2) / 2;
constexpr size_t num_phi_components = Dim * num_spacetime_metric_components;
constexpr size_t num_char_speed_components = Dim + 1;
constexpr size_t offset_v_spacetime_metric = 0;
constexpr size_t offset_v_zero =
    offset_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t offset_v_plus = offset_v_zero + num_phi_components;
constexpr size_t offset_v_minus = offset_v_plus + num_spacetime_metric_components;
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
  std::vector<ElementId<Dim>> element_ids;
  Mesh<Dim> mesh;
};

DomainData make_domain_data() {
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
  const auto analytic_vars = solution.variables(
      setup.inertial_coords, setup.time,
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

  dg::MortarMap<Dim, Mesh<Dim>> neighbor_mesh{};
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    for (const auto& neighbor : neighbors) {
      neighbor_mesh.insert(
          {DirectionalId<Dim>{direction, neighbor}, setup.mesh});
    }
  }
  setup.mortar_infos =
      evolution::dg::Initialization::detail::mortar_infos<Dim>(setup.element);
  auto [mortar_meshes, mortar_next_temporal_ids, normal_covector_and_magnitude] =
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
      setup.gamma0, setup.gamma1, setup.gamma2, setup.inertial_coords, setup.mesh,
      setup.element, setup.mortar_meshes, setup.mortar_infos);
  return setup;
}

outgoing_boundary_data_map compute_outgoing_boundary_data(
    const SetupData& setup, const size_t number_of_blocks) {
  auto device_vars = copy_to_device(setup.vars);
  typename gh::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  outgoing_boundary_data_map outgoing_boundary_data{};
  external_boundary_data_map external_boundary_data{};
  const external_boundary_conditions_type external_boundary_conditions_by_block{
      number_of_blocks};

  gh::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      make_not_null(&external_boundary_data), device_vars,
      setup.device_inverse_jacobian, setup.device_inertial_coordinates,
      setup.device_gamma0, setup.device_gamma1, setup.device_gamma2,
      setup.device_face_to_volume_index_map, setup.device_face_unit_normal_covector,
      setup.device_face_normal_magnitude, setup.device_mortar_data,
      setup.mortar_meshes, setup.gamma0, setup.gamma1, setup.gamma2,
      external_boundary_conditions_by_block, setup.time, setup.mesh, setup.element,
      setup.time_step_id);
  CHECK(external_boundary_data.empty());
  return outgoing_boundary_data;
}

incoming_boundary_data_map build_incoming_boundary_data(
    const SetupData& receiver_setup,
    const std::map<ElementId<Dim>, outgoing_boundary_data_map>&
        outgoing_by_element) {
  incoming_boundary_data_map incoming_boundary_data{};
  for (const auto& [direction, neighbors] : receiver_setup.element.neighbors()) {
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
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
        char_speeds,
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
      unnormalized_normal_vector.get(i) +=
          inverse_spatial_metric.get(i, j) * unnormalized_normal_covector.get(j);
    }
    normal_magnitude_squared +=
        unnormalized_normal_vector.get(i) * unnormalized_normal_covector.get(i);
  }
  return sqrt(normal_magnitude_squared);
}

void add_host_boundary_corrections(
    const gsl::not_null<dt_variables_type*> host_dt_reference,
    const SetupData& setup,
    const outgoing_boundary_data_map& local_outgoing_boundary_data,
    const incoming_boundary_data_map& incoming_boundary_data) {
  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    const Mesh<Dim - 1> face_mesh = setup.mesh.slice_away(direction.dimension());
    dt_variables_type dt_boundary_correction_on_face_sum{
        face_mesh.number_of_grid_points(), 0.0};

    for (const auto& neighbor : neighbors) {
      const DirectionalId<Dim> mortar_id{direction, neighbor};
      REQUIRE(local_outgoing_boundary_data.count(mortar_id) == 1);
      REQUIRE(incoming_boundary_data.count(mortar_id) == 1);
      const auto& local_boundary_data = local_outgoing_boundary_data.at(mortar_id);
      const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);
      const auto& mortar_mesh = setup.mortar_meshes.at(mortar_id);
      const auto& mortar_data = setup.device_mortar_data.at(mortar_id);
      const auto local_packaged_data_view = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{}, local_boundary_data.boundary_correction_data.view());
      const auto remote_packaged_data_view = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace{}, remote_boundary_data.boundary_correction_data.view());
      REQUIRE(local_packaged_data_view.extent(1) ==
              packaged_boundary_data_components);
      REQUIRE(remote_packaged_data_view.extent(1) ==
              packaged_boundary_data_components);

      dt_variables_type dt_boundary_correction_on_mortar{
          mortar_mesh.number_of_grid_points(), 0.0};
      auto& dt_spacetime_metric_on_mortar =
          get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              dt_boundary_correction_on_mortar);
      auto& dt_pi_on_mortar =
          get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(
              dt_boundary_correction_on_mortar);
      auto& dt_phi_on_mortar =
          get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(
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
          make_not_null(&local_v_spacetime_metric), make_not_null(&local_v_zero),
          make_not_null(&local_v_plus), make_not_null(&local_v_minus),
          make_not_null(&local_normal_times_v_plus),
          make_not_null(&local_normal_times_v_minus),
          make_not_null(&local_gamma2_v_spacetime_metric),
          make_not_null(&local_char_speeds), local_packaged_data_view,
          mortar_mesh.number_of_grid_points());
      unpack_packaged_boundary_data(
          make_not_null(&remote_v_spacetime_metric), make_not_null(&remote_v_zero),
          make_not_null(&remote_v_plus), make_not_null(&remote_v_minus),
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

    const size_t sliced_dim = direction.dimension();
    const size_t extent_perpendicular_to_boundary = setup.mesh.extents(sliced_dim);
    const double lift_prefactor =
        -0.5 * static_cast<double>(extent_perpendicular_to_boundary *
                                   (extent_perpendicular_to_boundary - 1));

    const auto host_face_to_volume_index =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_to_volume_index_map, sliced_dim).second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_to_volume_index_map, sliced_dim).first);
    const auto host_face_unit_normal_covector =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_unit_normal_covector, sliced_dim)
                      .second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_unit_normal_covector, sliced_dim).first);
    const auto host_face_normal_magnitude =
        direction.side() == Side::Upper
            ? Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_normal_magnitude, sliced_dim).second)
            : Kokkos::create_mirror_view_and_copy(
                  Kokkos::HostSpace{},
                  gsl::at(setup.device_face_normal_magnitude, sliced_dim).first);

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
          get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              *host_dt_reference)
              .get(a, b)[volume_index] +=
              lifted_factor *
              get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
                  dt_boundary_correction_on_face_sum)
                  .get(a, b)[face_index];
          get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(*host_dt_reference)
              .get(a, b)[volume_index] +=
              lifted_factor * get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(
                                  dt_boundary_correction_on_face_sum)
                                  .get(a, b)[face_index];
          for (size_t d = 0; d < Dim; ++d) {
            get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(*host_dt_reference)
                .get(d, a, b)[volume_index] +=
                lifted_factor *
                get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(
                    dt_boundary_correction_on_face_sum)
                    .get(d, a, b)[face_index];
          }
        }
      }
    }
  }

  CHECK(setup.element.external_boundaries().empty());
}

void test_apply_boundary_corrections_to_time_derivative_kokkos_matches_host() {
  const DomainData domain_data = make_domain_data();

  std::vector<SetupData> setups{};
  setups.reserve(domain_data.element_ids.size());
  for (const auto& element_id : domain_data.element_ids) {
    setups.push_back(make_setup_data(element_id, domain_data.domain,
                                     domain_data.initial_refinement_levels,
                                     domain_data.mesh));
  }

  std::map<ElementId<Dim>, outgoing_boundary_data_map> outgoing_by_element{};
  for (const auto& setup : setups) {
    outgoing_by_element.insert_or_assign(
        setup.element.id(),
        compute_outgoing_boundary_data(setup, domain_data.domain.blocks().size()));
  }

  std::map<ElementId<Dim>, incoming_boundary_data_map> incoming_by_element{};
  for (const auto& setup : setups) {
    incoming_by_element.insert_or_assign(
        setup.element.id(),
        build_incoming_boundary_data(setup, outgoing_by_element));
  }

  for (const auto& setup : setups) {
    INFO("Checking element " << setup.element.id());
    const size_t num_points = setup.mesh.number_of_grid_points();
    dt_variables_type host_dt_initial{num_points, 0.0};
    auto device_dt = copy_to_device(host_dt_initial);
    auto device_vars = copy_to_device(setup.vars);

    gh::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos::apply(
        make_not_null(&device_dt), device_vars,
        outgoing_by_element.at(setup.element.id()),
        incoming_by_element.at(setup.element.id()), external_boundary_data_map{},
        setup.device_face_to_volume_index_map,
        setup.device_face_unit_normal_covector, setup.device_face_normal_magnitude,
        setup.device_mortar_data, setup.mortar_meshes, setup.mesh, setup.element);

    dt_variables_type host_dt_from_kokkos{num_points, 0.0};
    copy_to_host(make_not_null(&host_dt_from_kokkos), device_dt);

    dt_variables_type host_dt_reference{num_points, 0.0};
    add_host_boundary_corrections(make_not_null(&host_dt_reference), setup,
                                  outgoing_by_element.at(setup.element.id()),
                                  incoming_by_element.at(setup.element.id()));

    const Approx approx = Approx::custom().epsilon(1e-8).scale(1.0);
    CHECK_ITERABLE_CUSTOM_APPROX(get<::Tags::dt<spacetime_metric_tag>>(
                                     host_dt_from_kokkos),
                                 get<::Tags::dt<spacetime_metric_tag>>(
                                     host_dt_reference),
                                 approx);
    CHECK_ITERABLE_CUSTOM_APPROX(
        get<::Tags::dt<pi_tag>>(host_dt_from_kokkos),
        get<::Tags::dt<pi_tag>>(host_dt_reference), approx);
    CHECK_ITERABLE_CUSTOM_APPROX(
        get<::Tags::dt<phi_tag>>(host_dt_from_kokkos),
        get<::Tags::dt<phi_tag>>(host_dt_reference), approx);
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic."
    "ApplyBoundaryCorrectionsToTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_apply_boundary_corrections_to_time_derivative_kokkos_matches_host();
}
