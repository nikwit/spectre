// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/SliceVariables.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/InterfaceHelpers.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/InternalMortarDataImpl.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/LiftFlux.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/RegularSphericalWave.hpp"
#include "PointwiseFunctions/MathFunctions/Gaussian.hpp"
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
using external_boundary_conditions_type =
    domain::Tags::ExternalBoundaryConditions<Dim>::type;

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
  double time{};
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

DomainData make_external_boundary_domain_data() {
  const domain::creators::Rectilinear<Dim> domain_creator{
      {{-1.0, -0.8, -0.6}},
      {{1.1, 0.7, 0.9}},
      {{0, 0, 0}},
      {{4, 5, 6}},
      {{false, false, false}}};
  auto domain = domain_creator.create_domain();
  auto initial_refinement_levels = domain_creator.initial_refinement_levels();
  auto element_ids = initial_element_ids(initial_refinement_levels);
  const Mesh<Dim> mesh{domain_creator.initial_extents()[0],
                       Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  return DomainData{std::move(domain), std::move(initial_refinement_levels),
                    std::move(element_ids), mesh};
}

external_boundary_conditions_type make_empty_external_boundary_conditions(
    const size_t number_of_blocks) {
  return external_boundary_conditions_type{number_of_blocks};
}

external_boundary_conditions_type make_dirichlet_external_boundary_conditions(
    const size_t number_of_blocks) {
  external_boundary_conditions_type boundary_conditions{number_of_blocks};

  auto make_bc =
      []() -> std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> {
    return std::make_unique<
        ScalarWave::BoundaryConditions::DirichletAnalytic<Dim>>(
        std::make_unique<ScalarWave::Solutions::RegularSphericalWave>(
            std::make_unique<MathFunctions::Gaussian<1, Frame::Inertial>>(
                1.0, 0.7, 0.0)));
  };

  for (size_t block_id = 0; block_id < number_of_blocks; ++block_id) {
    auto& block_boundary_conditions = boundary_conditions[block_id];
    for (size_t d = 0; d < Dim; ++d) {
      for (const Side side : {Side::Lower, Side::Upper}) {
        block_boundary_conditions.emplace(Direction<Dim>{d, side}, make_bc());
      }
    }
  }

  return boundary_conditions;
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
  setup.time = slab.start().value();
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

void apply_host_external_dirichlet_corrections(
    const SetupData& setup,
    const external_boundary_conditions_type&
        external_boundary_conditions_by_block,
    const gsl::not_null<Variables<dt_variables_tags>*> dt_reference) {
  const auto& external_boundary_conditions =
      external_boundary_conditions_by_block.at(setup.element.id().block_id());
  const boundary_correction correction{};
  const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>
      face_mesh_velocity{std::nullopt};
  const std::optional<Scalar<DataVector>> normal_dot_mesh_velocity{
      std::nullopt};

  for (const Direction<Dim>& direction : setup.element.external_boundaries()) {
    const auto* const dirichlet_analytic = dynamic_cast<
        const ScalarWave::BoundaryConditions::DirichletAnalytic<Dim>*>(
        external_boundary_conditions.at(direction).get());
    REQUIRE(dirichlet_analytic != nullptr);

    const size_t sliced_dim = direction.dimension();
    const Mesh<Dim - 1> face_mesh = setup.mesh.slice_away(sliced_dim);
    const size_t num_face_points = face_mesh.number_of_grid_points();
    const size_t slice_at = index_to_slice_at(setup.mesh.extents(), direction);

    Scalar<DataVector> interior_psi{num_face_points};
    Scalar<DataVector> interior_pi{num_face_points};
    tnsr::i<DataVector, Dim, Frame::Inertial> interior_phi{num_face_points};
    Scalar<DataVector> interior_gamma2{num_face_points};
    tnsr::i<DataVector, Dim, Frame::Inertial> interior_normal_covector{
        num_face_points};
    Scalar<DataVector> face_normal_magnitude{num_face_points};
    tnsr::I<DataVector, Dim, Frame::Inertial> face_coords{num_face_points};

    for (SliceIterator si(setup.mesh.extents(), sliced_dim, slice_at); si;
         ++si) {
      const size_t face_index = si.slice_offset();
      const size_t volume_index = si.volume_offset();

      get(interior_psi)[face_index] =
          get(get<ScalarWave::Tags::Psi>(setup.vars))[volume_index];
      get(interior_pi)[face_index] =
          get(get<ScalarWave::Tags::Pi>(setup.vars))[volume_index];
      for (size_t d = 0; d < Dim; ++d) {
        interior_phi.get(d)[face_index] =
            get<ScalarWave::Tags::Phi<Dim>>(setup.vars).get(d)[volume_index];
        face_coords.get(d)[face_index] =
            setup.inertial_coords.get(d)[volume_index];
      }
      get(interior_gamma2)[face_index] =
          get(setup.constraint_gamma2)[volume_index];

      double normal_magnitude_squared = 0.0;
      for (size_t d = 0; d < Dim; ++d) {
        const double unnormalized_normal_component =
            direction.sign() *
            setup.inverse_jacobian.get(sliced_dim, d)[volume_index];
        interior_normal_covector.get(d)[face_index] =
            unnormalized_normal_component;
        normal_magnitude_squared +=
            unnormalized_normal_component * unnormalized_normal_component;
      }
      const double normal_magnitude = std::sqrt(normal_magnitude_squared);
      get(face_normal_magnitude)[face_index] = normal_magnitude;
      for (size_t d = 0; d < Dim; ++d) {
        interior_normal_covector.get(d)[face_index] /= normal_magnitude;
      }
    }

    Variables<dg_package_field_tags> internal_packaged_data{num_face_points};
    (void)correction.dg_package_data(
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<0>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<1>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<2>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<3>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<4>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<5>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<6>>>(
            internal_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<7>>>(
            internal_packaged_data)),
        interior_psi, interior_pi, interior_phi, interior_gamma2,
        interior_normal_covector, face_mesh_velocity, normal_dot_mesh_velocity);

    Scalar<DataVector> exterior_psi{num_face_points};
    Scalar<DataVector> exterior_pi{num_face_points};
    tnsr::i<DataVector, Dim, Frame::Inertial> exterior_phi{num_face_points};
    Scalar<DataVector> exterior_gamma2{num_face_points};
    const auto error_message = dirichlet_analytic->dg_ghost(
        make_not_null(&exterior_psi), make_not_null(&exterior_pi),
        make_not_null(&exterior_phi), make_not_null(&exterior_gamma2),
        face_mesh_velocity, interior_normal_covector, face_coords,
        interior_gamma2, setup.time);
    REQUIRE_FALSE(error_message.has_value());

    tnsr::i<DataVector, Dim, Frame::Inertial> exterior_normal_covector{
        num_face_points};
    for (size_t d = 0; d < Dim; ++d) {
      exterior_normal_covector.get(d) = -interior_normal_covector.get(d);
    }

    Variables<dg_package_field_tags> external_packaged_data{num_face_points};
    (void)correction.dg_package_data(
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<0>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<1>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<2>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<3>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<4>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<5>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<6>>>(
            external_packaged_data)),
        make_not_null(&get<tmpl::at<dg_package_field_tags, tmpl::size_t<7>>>(
            external_packaged_data)),
        exterior_psi, exterior_pi, exterior_phi, exterior_gamma2,
        exterior_normal_covector, face_mesh_velocity, normal_dot_mesh_velocity);

    Variables<dt_variables_tags> boundary_corrections_on_face{num_face_points};
    correction.dg_boundary_terms(
        make_not_null(&get<::Tags::dt<ScalarWave::Tags::Psi>>(
            boundary_corrections_on_face)),
        make_not_null(&get<::Tags::dt<ScalarWave::Tags::Pi>>(
            boundary_corrections_on_face)),
        make_not_null(&get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(
            boundary_corrections_on_face)),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<0>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<1>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<2>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<3>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<4>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<5>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<6>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<7>>>(
            internal_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<0>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<1>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<2>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<3>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<4>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<5>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<6>>>(
            external_packaged_data),
        get<tmpl::at<dg_package_field_tags, tmpl::size_t<7>>>(
            external_packaged_data),
        ::dg::Formulation::StrongInertial);

    ::dg::lift_flux(make_not_null(&boundary_corrections_on_face),
                    setup.mesh.extents(direction.dimension()),
                    face_normal_magnitude);
    add_slice_to_data(dt_reference, boundary_corrections_on_face,
                      setup.mesh.extents(), sliced_dim, slice_at);
  }
}

void check_compute_and_package_for_element(const SetupData& setup,
                                           const size_t number_of_blocks) {
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
  typename ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
      device_face_unit_normal_covector{};
  typename ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
      device_face_normal_magnitude{};
  typename ScalarWave::KokkosTags::DeviceMortarData<Dim>::type
      device_mortar_data{};
  typename ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
      outgoing_boundary_data{};
  typename ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type
      external_boundary_data{};

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

  ScalarWave::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&device_inverse_jacobian),
      make_not_null(&device_constraint_gamma2),
      make_not_null(&device_face_to_volume_index_map),
      make_not_null(&device_face_unit_normal_covector),
      make_not_null(&device_face_normal_magnitude),
      make_not_null(&device_mortar_data), setup.inverse_jacobian,
      setup.constraint_gamma2, setup.mesh, setup.element, mortar_meshes,
      mortar_infos);

  // The fused iterable action computes, packages, and sends. This unit test
  // exercises the shared compute+package path directly.
  const auto external_boundary_conditions =
      make_empty_external_boundary_conditions(number_of_blocks);
  ScalarWave::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      make_not_null(&external_boundary_data), device_vars,
      device_inverse_jacobian, device_constraint_gamma2,
      device_face_to_volume_index_map, device_face_unit_normal_covector,
      device_mortar_data, mortar_meshes, setup.constraint_gamma2,
      external_boundary_conditions, setup.inertial_coords, setup.time,
      setup.mesh, setup.element, setup.time_step_id);

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
      CHECK_ITERABLE_APPROX(
          gsl::make_span(packaged_kokkos.data(), packaged_kokkos.size()),
          gsl::make_span(packaged_host.data(), packaged_host.size()));
    }
  }
}

void test_compute_and_package_matches_host_internal_reference() {
  const DomainData domain_data = make_domain_data();
  for (const auto& element_id : domain_data.element_ids) {
    INFO("Checking element " << element_id);
    const SetupData setup = make_setup_data(
        element_id, domain_data.domain, domain_data.initial_refinement_levels,
        domain_data.mesh);
    check_compute_and_package_for_element(setup,
                                          domain_data.domain.blocks().size());
  }
}

void test_packaged_data_populates_inboxes_with_oriented_metadata() {
  using inbox_tag = ScalarWave::KokkosTags::BoundaryCorrectionInbox<Dim, false>;
  using inbox_type = typename inbox_tag::type;

  const DomainData domain_data = make_domain_data();

  for (const auto& element_id : domain_data.element_ids) {
    INFO("Checking sender element " << element_id);
    const SetupData setup = make_setup_data(
        element_id, domain_data.domain, domain_data.initial_refinement_levels,
        domain_data.mesh);
    typename ScalarWave::KokkosTags::DeviceVariables<system>::type device_vars =
        copy_to_device(setup.vars);
    typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type
        device_dt{};
    typename ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>::type
        device_inverse_jacobian{};
    typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type
        device_constraint_gamma2{};
    typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type
        device_face_to_volume_index_map{};
    typename ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
        device_face_unit_normal_covector{};
    typename ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
        device_face_normal_magnitude{};
    typename ScalarWave::KokkosTags::DeviceMortarData<Dim>::type
        device_mortar_data{};
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

    ScalarWave::Actions::InitializeKokkosTags<system>::apply(
        make_not_null(&device_inverse_jacobian),
        make_not_null(&device_constraint_gamma2),
        make_not_null(&device_face_to_volume_index_map),
        make_not_null(&device_face_unit_normal_covector),
        make_not_null(&device_face_normal_magnitude),
        make_not_null(&device_mortar_data), setup.inverse_jacobian,
        setup.constraint_gamma2, setup.mesh, setup.element, mortar_meshes,
        mortar_infos);

    typename ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
        outgoing_boundary_data{};
    typename ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type
        external_boundary_data{};
    const auto external_boundary_conditions =
        make_empty_external_boundary_conditions(
            domain_data.domain.blocks().size());
    ScalarWave::Actions::ComputeTimeDerivativeKokkos::apply(
        make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
        make_not_null(&external_boundary_data), device_vars,
        device_inverse_jacobian, device_constraint_gamma2,
        device_face_to_volume_index_map, device_face_unit_normal_covector,
        device_mortar_data, mortar_meshes, setup.constraint_gamma2,
        external_boundary_conditions, setup.inertial_coords, setup.time,
        setup.mesh, setup.element, setup.time_step_id);

    Variables<dt_variables_tags> unused_dt{
        domain_data.mesh.number_of_grid_points(), 0.0};
    packaged_data_by_mortar_map host_packaged_data{};
    compute_host_reference(setup, make_not_null(&unused_dt),
                           make_not_null(&host_packaged_data));

    for (const auto& [direction, neighbors] : setup.element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        const auto& orientation = neighbors.orientation(neighbor);
        const DirectionalId<Dim> sender_mortar_id{direction, neighbor};
        const DirectionalId<Dim> expected_neighbor_inbox_key{
            orientation(direction.opposite()), setup.element.id()};
        REQUIRE(outgoing_boundary_data.count(sender_mortar_id) == 1);

        auto data_for_neighbor = outgoing_boundary_data.at(sender_mortar_id);

        inbox_type neighbor_inbox{};
        const bool inbox_ready = inbox_tag::insert_into_inbox(
            make_not_null(&neighbor_inbox), setup.time_step_id,
            std::make_pair(expected_neighbor_inbox_key,
                           std::move(data_for_neighbor)));
        CHECK_FALSE(inbox_ready);

        const auto inbox_record = neighbor_inbox.find(setup.time_step_id);
        REQUIRE(inbox_record != neighbor_inbox.end());
        REQUIRE(inbox_record->second.count(expected_neighbor_inbox_key) == 1);
        REQUIRE(host_packaged_data.count(sender_mortar_id) == 1);

        const auto& received_data =
            inbox_record->second.at(expected_neighbor_inbox_key);
        CHECK(received_data.validity_range == setup.time_step_id);
        CHECK(received_data.integration_order == 0);
        const auto& expected_host_packaged_data =
            host_packaged_data.at(sender_mortar_id);
        REQUIRE(received_data.boundary_correction_data.size() ==
                expected_host_packaged_data.size());

        Variables<dg_package_field_tags> received_packaged_data{
            received_data.boundary_correction_data.number_of_grid_points(),
            0.0};
        copy_to_host(make_not_null(&received_packaged_data),
                     received_data.boundary_correction_data);
        CHECK(received_packaged_data.size() ==
              expected_host_packaged_data.size());
        CHECK_ITERABLE_APPROX(
            gsl::make_span(received_packaged_data.data(),
                           received_packaged_data.size()),
            gsl::make_span(expected_host_packaged_data.data(),
                           expected_host_packaged_data.size()));
      }
    }
  }
}

void test_initialize_kokkos_tags_sets_device_mortar_data() {
  const DomainData domain_data = make_domain_data();

  for (const auto& element_id : domain_data.element_ids) {
    INFO("Checking initialized mortar metadata for element " << element_id);
    const SetupData setup = make_setup_data(
        element_id, domain_data.domain, domain_data.initial_refinement_levels,
        domain_data.mesh);

    typename ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>::type
        device_inverse_jacobian{};
    typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type
        device_constraint_gamma2{};
    typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type
        device_face_to_volume_index_map{};
    typename ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
        device_face_unit_normal_covector{};
    typename ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
        device_face_normal_magnitude{};
    typename ScalarWave::KokkosTags::DeviceMortarData<Dim>::type
        device_mortar_data{};

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

    ScalarWave::Actions::InitializeKokkosTags<system>::apply(
        make_not_null(&device_inverse_jacobian),
        make_not_null(&device_constraint_gamma2),
        make_not_null(&device_face_to_volume_index_map),
        make_not_null(&device_face_unit_normal_covector),
        make_not_null(&device_face_normal_magnitude),
        make_not_null(&device_mortar_data), setup.inverse_jacobian,
        setup.constraint_gamma2, setup.mesh, setup.element, mortar_meshes,
        mortar_infos);

    CHECK(device_mortar_data.size() == setup.element.number_of_neighbors());
    for (const auto& [direction, neighbors] : setup.element.neighbors()) {
      const size_t sliced_dim = direction.dimension();
      const Mesh<Dim - 1> face_mesh = setup.mesh.slice_away(sliced_dim);
      for (const auto& neighbor : neighbors) {
        const DirectionalId<Dim> mortar_id{direction, neighbor};
        REQUIRE(device_mortar_data.count(mortar_id) == 1);
        const auto& initialized_mortar_data = device_mortar_data.at(mortar_id);

        const auto& mortar_mesh = mortar_meshes.at(mortar_id);
        const auto& mortar_size = mortar_infos.at(mortar_id).mortar_size();
        CHECK(initialized_mortar_data.needs_projection ==
              Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size));
        CHECK(initialized_mortar_data.mortar_size == mortar_size);

        std::vector<double> grid_point_indices(
            mortar_mesh.number_of_grid_points());
        std::iota(grid_point_indices.begin(), grid_point_indices.end(), 0.0);
        const std::vector<double> expected_oriented_indices =
            orient_variables_on_slice(grid_point_indices, mortar_mesh.extents(),
                                      sliced_dim,
                                      neighbors.orientation(neighbor));
        const auto actual_oriented_indices =
            Kokkos::create_mirror_view_and_copy(
                Kokkos::HostSpace{},
                initialized_mortar_data
                    .oriented_mortar_grid_point_source_index);
        REQUIRE(actual_oriented_indices.extent(0) ==
                expected_oriented_indices.size());
        for (size_t i = 0; i < expected_oriented_indices.size(); ++i) {
          CHECK(actual_oriented_indices(i) ==
                static_cast<size_t>(expected_oriented_indices[i]));
        }
      }
    }
  }
}

void test_compute_time_derivative_matches_host_with_external_dirichlet() {
  const DomainData domain_data = make_external_boundary_domain_data();
  REQUIRE(domain_data.element_ids.size() == 1);

  const SetupData setup =
      make_setup_data(domain_data.element_ids[0], domain_data.domain,
                      domain_data.initial_refinement_levels, domain_data.mesh);
  REQUIRE_FALSE(setup.element.external_boundaries().empty());
  REQUIRE(setup.element.number_of_neighbors() == 0);

  typename ScalarWave::KokkosTags::DeviceVariables<system>::type device_vars =
      copy_to_device(setup.vars);
  typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  typename ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>::type
      device_inverse_jacobian{};
  typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type
      device_constraint_gamma2{};
  typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type
      device_face_to_volume_index_map{};
  typename ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
      device_face_unit_normal_covector{};
  typename ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
      device_face_normal_magnitude{};
  typename ScalarWave::KokkosTags::DeviceMortarData<Dim>::type
      device_mortar_data{};
  typename ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type
      outgoing_boundary_data{};
  typename ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type
      external_boundary_data{};

  dg::MortarMap<Dim, Mesh<Dim>> neighbor_mesh{};
  auto mortar_infos =
      evolution::dg::Initialization::detail::mortar_infos<Dim>(setup.element);
  auto [mortar_meshes, mortar_next_temporal_ids,
        normal_covector_and_magnitude] =
      evolution::dg::Initialization::detail::mortars_apply_impl<Dim>(
          setup.element, setup.time_step_id, setup.mesh, neighbor_mesh);
  (void)mortar_next_temporal_ids;
  (void)normal_covector_and_magnitude;

  ScalarWave::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&device_inverse_jacobian),
      make_not_null(&device_constraint_gamma2),
      make_not_null(&device_face_to_volume_index_map),
      make_not_null(&device_face_unit_normal_covector),
      make_not_null(&device_face_normal_magnitude),
      make_not_null(&device_mortar_data), setup.inverse_jacobian,
      setup.constraint_gamma2, setup.mesh, setup.element, mortar_meshes,
      mortar_infos);

  const auto external_boundary_conditions =
      make_dirichlet_external_boundary_conditions(
          domain_data.domain.blocks().size());
  ScalarWave::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      make_not_null(&external_boundary_data), device_vars,
      device_inverse_jacobian, device_constraint_gamma2,
      device_face_to_volume_index_map, device_face_unit_normal_covector,
      device_mortar_data, mortar_meshes, setup.constraint_gamma2,
      external_boundary_conditions, setup.inertial_coords, setup.time,
      setup.mesh, setup.element, setup.time_step_id);

  const typename ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<
      Dim>::type incoming_boundary_data{};
  ScalarWave::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), outgoing_boundary_data, incoming_boundary_data,
      external_boundary_data, device_face_to_volume_index_map,
      device_face_normal_magnitude, device_mortar_data, mortar_meshes,
      setup.mesh, setup.element);

  Variables<dt_variables_tags> dt_kokkos{setup.mesh.number_of_grid_points(),
                                         0.0};
  copy_to_host(make_not_null(&dt_kokkos), device_dt);

  Variables<dt_variables_tags> dt_reference{setup.mesh.number_of_grid_points(),
                                            0.0};
  packaged_data_by_mortar_map host_packaged_data_by_mortar{};
  compute_host_reference(setup, make_not_null(&dt_reference),
                         make_not_null(&host_packaged_data_by_mortar));
  apply_host_external_dirichlet_corrections(setup, external_boundary_conditions,
                                            make_not_null(&dt_reference));

  CHECK_ITERABLE_APPROX(get<::Tags::dt<ScalarWave::Tags::Psi>>(dt_kokkos),
                        get<::Tags::dt<ScalarWave::Tags::Psi>>(dt_reference));
  CHECK_ITERABLE_APPROX(get<::Tags::dt<ScalarWave::Tags::Pi>>(dt_kokkos),
                        get<::Tags::dt<ScalarWave::Tags::Pi>>(dt_reference));
  CHECK_ITERABLE_APPROX(
      get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(dt_kokkos),
      get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(dt_reference));
  CHECK(outgoing_boundary_data.size() ==
        setup.element.external_boundaries().size());
  CHECK(external_boundary_data.size() ==
        setup.element.external_boundaries().size());
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarWave.ComputeTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_compute_and_package_matches_host_internal_reference();
  test_packaged_data_populates_inboxes_with_oriented_metadata();
  test_initialize_kokkos_tags_sets_device_mortar_data();
  test_compute_time_derivative_matches_host_with_external_dirichlet();
}
