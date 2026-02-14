// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
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
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/InternalMortarDataImpl.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/Executables/ScalarWave/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Executables/ScalarWave/InitializeKokkosTags.hpp"
#include "Evolution/Executables/ScalarWave/KokkosBoundaryCommunication.hpp"
#include "Evolution/Executables/ScalarWave/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Time/Slab.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
constexpr size_t Dim = 3;
using system = ScalarWave::System<Dim>;
using boundary_correction = ScalarWave::BoundaryCorrections::UpwindPenalty<Dim>;

using dt_variables_tags =
    db::wrap_tags_in<::Tags::dt, typename system::variables_tag::tags_list>;
using flux_variables_tags =
    db::wrap_tags_in<::Tags::Flux, typename system::flux_variables,
                     tmpl::size_t<Dim>, Frame::Inertial>;
using derivative_variables_tags =
    db::wrap_tags_in<::Tags::deriv, typename system::gradient_variables,
                     tmpl::size_t<Dim>, Frame::Inertial>;
using temporary_tags =
    typename system::compute_volume_time_derivative_terms::temporary_tags;
using dg_package_field_tags =
    typename boundary_correction::dg_package_field_tags;
using packaged_data_by_mortar_map = dg::MortarMap<Dim, DataVector>;

struct SetupData {
  Element<Dim> element{};
  Mesh<Dim> mesh{};
  typename system::variables_tag::type vars{};
  ScalarWave::Tags::ConstraintGamma2::type constraint_gamma2{};
  tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords{};
  domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                Frame::Inertial>::type inverse_jacobian{};
  TimeStepId time_step_id{};
  std::unique_ptr<domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>>
      moving_mesh_map{};
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

  const size_t num_points = setup.mesh.number_of_grid_points();
  setup.vars = typename system::variables_tag::type{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(setup.vars)) =
      0.2 + 0.1 * setup.inertial_coords.get(0) -
      0.03 * setup.inertial_coords.get(1) + 0.07 * setup.inertial_coords.get(2);
  get(get<ScalarWave::Tags::Pi>(setup.vars)) =
      -0.5 + 0.08 * setup.inertial_coords.get(0) +
      0.04 * setup.inertial_coords.get(1) - 0.02 * setup.inertial_coords.get(2);
  for (size_t d = 0; d < Dim; ++d) {
    get<ScalarWave::Tags::Phi<Dim>>(setup.vars).get(d) =
        0.15 * static_cast<double>(d + 1) +
        (0.01 * static_cast<double>(d + 2)) * setup.inertial_coords.get(0) -
        0.02 * setup.inertial_coords.get(1) +
        0.03 * setup.inertial_coords.get(2);
  }

  setup.constraint_gamma2 =
      ScalarWave::Tags::ConstraintGamma2::type{num_points, 0.0};
  get(setup.constraint_gamma2) = 0.6 - 0.05 * setup.inertial_coords.get(2);

  const Slab slab(0.0, 1.0);
  setup.time_step_id = TimeStepId(true, 0, slab.start());
  setup.moving_mesh_map =
      domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
          domain::CoordinateMaps::Identity<Dim>{});
  return setup;
}

void compute_host_reference(
    const SetupData& setup,
    const gsl::not_null<Variables<dt_variables_tags>*> dt_reference,
    const gsl::not_null<packaged_data_by_mortar_map*> packaged_data_by_mortar) {
  const size_t num_points = setup.mesh.number_of_grid_points();
  Variables<flux_variables_tags> volume_fluxes{
      setup.mesh.number_of_grid_points()};
  Variables<derivative_variables_tags> partial_derivs{
      setup.mesh.number_of_grid_points()};
  Variables<temporary_tags> temporaries{setup.mesh.number_of_grid_points()};

  partial_derivatives(make_not_null(&partial_derivs), setup.vars, setup.mesh,
                      setup.inverse_jacobian);

  dt_reference->initialize(num_points, 0.0);
  system::compute_volume_time_derivative_terms::apply(
      make_not_null(&get<::Tags::dt<ScalarWave::Tags::Psi>>(*dt_reference)),
      make_not_null(&get<::Tags::dt<ScalarWave::Tags::Pi>>(*dt_reference)),
      make_not_null(
          &get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(*dt_reference)),
      make_not_null(&get<ScalarWave::Tags::ConstraintGamma2>(temporaries)),
      get<::Tags::deriv<ScalarWave::Tags::Psi, tmpl::size_t<Dim>,
                        Frame::Inertial>>(partial_derivs),
      get<::Tags::deriv<ScalarWave::Tags::Pi, tmpl::size_t<Dim>,
                        Frame::Inertial>>(partial_derivs),
      get<::Tags::deriv<ScalarWave::Tags::Phi<Dim>, tmpl::size_t<Dim>,
                        Frame::Inertial>>(partial_derivs),
      get<ScalarWave::Tags::Pi>(setup.vars),
      get<ScalarWave::Tags::Phi<Dim>>(setup.vars), setup.constraint_gamma2);

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
      dg_package_data_projected_tags,
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
  evolution::dg::Actions::detail::internal_mortar_data_impl<
      system, Dim, boundary_correction, temporary_tags>(
      make_not_null(&normal_covector_and_magnitude),
      make_not_null(&host_mortar_data), make_not_null(&face_temporaries_span),
      make_not_null(&packaged_data_span), correction, setup.vars, volume_fluxes,
      temporaries, static_cast<const Variables<tmpl::list<>>*>(nullptr),
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

void check_compute_time_derivative_kokkos_matches_host_internal_reference_for_element(
    const SetupData& setup) {
  const size_t num_points = setup.mesh.number_of_grid_points();

  // Kokkos implementation under test
  typename ScalarWave::KokkosTags::DeviceVariables<system>::type device_vars =
      copy_to_device(setup.vars);
  typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  typename ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>::type
      device_inverse_jacobian{};
  typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type
      device_constraint_gamma2{};
  typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type
      device_face_to_volume_index_map{};
  typename ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
      outgoing_boundary_data{};

  ScalarWave::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&device_inverse_jacobian),
      make_not_null(&device_constraint_gamma2),
      make_not_null(&device_face_to_volume_index_map), setup.inverse_jacobian,
      setup.constraint_gamma2, setup.mesh);

  ScalarWave::Actions::ComputeTimeDerivativeKokkos<Dim, system>::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      device_vars, device_inverse_jacobian, device_constraint_gamma2,
      device_face_to_volume_index_map, setup.mesh, setup.element,
      setup.time_step_id);

  Variables<dt_variables_tags> dt_kokkos{num_points, 0.0};
  copy_to_host(make_not_null(&dt_kokkos), device_dt);

  // Host reference using direct volume terms and internal_mortar_data_impl.
  Variables<dt_variables_tags> dt_reference{num_points, 0.0};
  packaged_data_by_mortar_map host_packaged_data_by_mortar{};
  compute_host_reference(setup, make_not_null(&dt_reference),
                         make_not_null(&host_packaged_data_by_mortar));

  CHECK_ITERABLE_APPROX(get<::Tags::dt<ScalarWave::Tags::Psi>>(dt_kokkos),
                        get<::Tags::dt<ScalarWave::Tags::Psi>>(dt_reference));
  CHECK_ITERABLE_APPROX(get<::Tags::dt<ScalarWave::Tags::Pi>>(dt_kokkos),
                        get<::Tags::dt<ScalarWave::Tags::Pi>>(dt_reference));
  CHECK_ITERABLE_APPROX(
      get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(dt_kokkos),
      get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(dt_reference));

  CHECK(outgoing_boundary_data.size() == setup.element.number_of_neighbors());

  for (const auto& [direction, neighbors] : setup.element.neighbors()) {
    for (const auto& neighbor : neighbors) {
      const DirectionalId<Dim> mortar_id{direction, neighbor};
      CHECK(outgoing_boundary_data.count(mortar_id) == 1);
      CHECK(host_packaged_data_by_mortar.count(mortar_id) == 1);
      const auto& packaged_host = host_packaged_data_by_mortar.at(mortar_id);

      const auto& kokkos_boundary_data = outgoing_boundary_data.at(mortar_id);
      CHECK(kokkos_boundary_data.volume_mesh == setup.mesh);
      const Mesh<Dim - 1> face_mesh =
          setup.mesh.slice_away(direction.dimension());
      CHECK(kokkos_boundary_data.boundary_correction_mesh == face_mesh);
      CHECK(kokkos_boundary_data.validity_range == setup.time_step_id);
      CHECK(kokkos_boundary_data.integration_order == 0);

      Variables<dg_package_field_tags> packaged_kokkos{
          kokkos_boundary_data.boundary_correction_mesh.number_of_grid_points(),
          0.0};
      copy_to_host(make_not_null(&packaged_kokkos),
                   kokkos_boundary_data.boundary_correction_data);
      CHECK(packaged_kokkos.size() == packaged_host.size());
      CHECK_ITERABLE_APPROX(
          gsl::make_span(packaged_kokkos.data(), packaged_kokkos.size()),
          gsl::make_span(packaged_host.data(), packaged_host.size()));
    }
  }
}

void test_compute_time_derivative_kokkos_matches_host_internal_reference() {
  const DomainData domain_data = make_domain_data();
  for (const auto& element_id : domain_data.element_ids) {
    INFO("Checking element " << element_id);
    const SetupData setup = make_setup_data(
        element_id, domain_data.domain, domain_data.initial_refinement_levels,
        domain_data.mesh);
    check_compute_time_derivative_kokkos_matches_host_internal_reference_for_element(
        setup);
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarWave.ComputeTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_compute_time_derivative_kokkos_matches_host_internal_reference();
}
