// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
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
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/ComputeTimeDerivativeKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/InitializeKokkosTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {
constexpr size_t Dim = 3;
using system = ScalarWave::System<Dim>;
using boundary_correction = ScalarWave::BoundaryCorrections::UpwindPenalty<Dim>;

using dt_variables_tags =
    db::wrap_tags_in<::Tags::dt, typename system::variables_tag::tags_list>;
using dg_package_field_tags =
    typename boundary_correction::dg_package_field_tags;
template <size_t I>
using package_field_tag = tmpl::at<dg_package_field_tags, tmpl::size_t<I>>;

using outgoing_boundary_data_map =
    typename ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>::type;
using incoming_boundary_data_map =
    typename ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<Dim>::type;
using external_boundary_data_map =
    typename ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>::type;
using device_face_to_volume_index_map_type =
    typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>::type;
using apply_action =
    ScalarWave::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos;
using inbox_tag = ScalarWave::KokkosTags::BoundaryCorrectionInbox<Dim, false>;

struct SetupData {
  Element<Dim> element{};
  Mesh<Dim> mesh{};
  typename system::variables_tag::type vars{};
  tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords{};
  ScalarWave::Tags::ConstraintGamma2::type constraint_gamma2{};
  domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                Frame::Inertial>::type inverse_jacobian{};
  typename ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>::type
      device_inverse_jacobian{};
  typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type
      device_constraint_gamma2{};
  TimeStepId time_step_id{};
  device_face_to_volume_index_map_type device_face_to_volume_index_map{};
  typename ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<Dim>::type
      device_face_unit_normal_covector{};
  typename ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>::type
      device_face_normal_magnitude{};
  typename ScalarWave::KokkosTags::DeviceMortarData<Dim>::type
      device_mortar_data{};
  dg::MortarMap<Dim, Mesh<Dim - 1>> mortar_meshes{};
  dg::MortarMap<Dim, evolution::dg::MortarInfo<Dim>> mortar_infos{};
  double time{};
};

struct DomainData {
  Domain<Dim> domain;
  std::vector<std::array<size_t, Dim>> initial_refinement_levels;
  std::vector<ElementId<Dim>> element_ids;
  Mesh<Dim> mesh;
};

struct TestMetavariables {
  using component_list = tmpl::list<>;
};

struct DummyParallelComponent;

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
  const auto inertial_coords = element_map(logical_coords);
  setup.inertial_coords = inertial_coords;
  setup.inverse_jacobian = element_map.inv_jacobian(logical_coords);

  const size_t num_points = setup.mesh.number_of_grid_points();
  setup.vars = typename system::variables_tag::type{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(setup.vars)) =
      0.2 + 0.1 * inertial_coords.get(0) - 0.03 * inertial_coords.get(1) +
      0.07 * inertial_coords.get(2);
  get(get<ScalarWave::Tags::Pi>(setup.vars)) =
      -0.5 + 0.08 * inertial_coords.get(0) + 0.04 * inertial_coords.get(1) -
      0.02 * inertial_coords.get(2);
  for (size_t d = 0; d < Dim; ++d) {
    get<ScalarWave::Tags::Phi<Dim>>(setup.vars).get(d) =
        0.15 * static_cast<double>(d + 1) +
        (0.01 * static_cast<double>(d + 2)) * inertial_coords.get(0) -
        0.02 * inertial_coords.get(1) + 0.03 * inertial_coords.get(2);
  }

  setup.constraint_gamma2 =
      ScalarWave::Tags::ConstraintGamma2::type{num_points, 0.0};
  get(setup.constraint_gamma2) = 0.6 - 0.05 * inertial_coords.get(2);

  const Slab slab(0.0, 1.0);
  setup.time_step_id = TimeStepId(true, 0, slab.start());
  setup.time = slab.start().value();

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

  ScalarWave::Actions::InitializeKokkosTags<system>::apply(
      make_not_null(&setup.device_inverse_jacobian),
      make_not_null(&setup.device_constraint_gamma2),
      make_not_null(&setup.device_face_to_volume_index_map),
      make_not_null(&setup.device_face_unit_normal_covector),
      make_not_null(&setup.device_face_normal_magnitude),
      make_not_null(&setup.device_mortar_data), setup.inverse_jacobian,
      setup.constraint_gamma2, setup.mesh, setup.element, setup.mortar_meshes,
      setup.mortar_infos);
  return setup;
}

outgoing_boundary_data_map compute_outgoing_boundary_data(
    const SetupData& setup) {
  typename ScalarWave::KokkosTags::DeviceVariables<system>::type device_vars =
      copy_to_device(setup.vars);
  typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type device_dt{};
  outgoing_boundary_data_map outgoing_boundary_data{};
  external_boundary_data_map external_boundary_data{};
  const auto external_boundary_conditions =
      domain::Tags::ExternalBoundaryConditions<Dim>::type{
          setup.element.id().block_id() + 1};

  ScalarWave::Actions::ComputeTimeDerivativeKokkos::apply(
      make_not_null(&device_dt), make_not_null(&outgoing_boundary_data),
      make_not_null(&external_boundary_data), device_vars,
      setup.device_inverse_jacobian, setup.device_constraint_gamma2,
      setup.device_face_to_volume_index_map,
      setup.device_face_unit_normal_covector, setup.device_mortar_data,
      setup.mortar_meshes, setup.constraint_gamma2,
      external_boundary_conditions, setup.inertial_coords, setup.time,
      setup.mesh, setup.element, setup.time_step_id);
  CHECK(external_boundary_data.empty());
  return outgoing_boundary_data;
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

void add_host_boundary_corrections(
    const gsl::not_null<Variables<dt_variables_tags>*> host_dt_reference,
    const outgoing_boundary_data_map& local_outgoing_boundary_data,
    const incoming_boundary_data_map& incoming_boundary_data,
    const dg::MortarMap<Dim, Mesh<Dim - 1>>& mortar_meshes,
    const domain::Tags::InverseJacobian<
        Dim, Frame::ElementLogical, Frame::Inertial>::type& inverse_jacobian,
    const device_face_to_volume_index_map_type& device_face_to_volume_index_map,
    const Mesh<Dim>& mesh, const Element<Dim>& element) {
  for (const auto& [direction, neighbors] : element.neighbors()) {
    REQUIRE(neighbors.size() == 1);
    const auto& neighbor = *neighbors.begin();
    const DirectionalId<Dim> mortar_id{direction, neighbor};
    REQUIRE(local_outgoing_boundary_data.count(mortar_id) == 1);
    REQUIRE(incoming_boundary_data.count(mortar_id) == 1);

    const auto& local_boundary_data =
        local_outgoing_boundary_data.at(mortar_id);
    const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);
    const auto& mortar_mesh = mortar_meshes.at(mortar_id);
    REQUIRE(
        local_boundary_data.boundary_correction_data.number_of_grid_points() ==
        mortar_mesh.number_of_grid_points());
    REQUIRE(
        remote_boundary_data.boundary_correction_data.number_of_grid_points() ==
        mortar_mesh.number_of_grid_points());

    const size_t num_face_points = mortar_mesh.number_of_grid_points();
    Variables<dg_package_field_tags> local_packaged_data{num_face_points, 0.0};
    Variables<dg_package_field_tags> remote_packaged_data{num_face_points, 0.0};
    copy_to_host(make_not_null(&local_packaged_data),
                 local_boundary_data.boundary_correction_data);
    copy_to_host(make_not_null(&remote_packaged_data),
                 remote_boundary_data.boundary_correction_data);

    const size_t sliced_dim = direction.dimension();
    const auto lower_face_to_volume_index = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{},
        gsl::at(device_face_to_volume_index_map, sliced_dim).first);
    const auto upper_face_to_volume_index = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace{},
        gsl::at(device_face_to_volume_index_map, sliced_dim).second);
    const size_t extent_perpendicular_to_boundary = mesh.extents(sliced_dim);
    const double lift_prefactor =
        -0.5 * static_cast<double>(extent_perpendicular_to_boundary *
                                   (extent_perpendicular_to_boundary - 1));
    const double outward_sign = direction.side() == Side::Upper ? 1.0 : -1.0;

    Scalar<DataVector> psi_boundary_correction{num_face_points, 0.0};
    Scalar<DataVector> pi_boundary_correction{num_face_points, 0.0};
    tnsr::i<DataVector, Dim, Frame::Inertial> phi_boundary_correction{
        num_face_points, 0.0};
    ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<Dim,
                                                                    DataVector>(
        make_not_null(&psi_boundary_correction),
        make_not_null(&pi_boundary_correction),
        make_not_null(&phi_boundary_correction),
        get<package_field_tag<0>>(local_packaged_data),
        get<package_field_tag<1>>(local_packaged_data),
        get<package_field_tag<2>>(local_packaged_data),
        get<package_field_tag<3>>(local_packaged_data),
        get<package_field_tag<4>>(local_packaged_data),
        get<package_field_tag<5>>(local_packaged_data),
        get<package_field_tag<6>>(local_packaged_data),
        get<package_field_tag<7>>(local_packaged_data),
        get<package_field_tag<0>>(remote_packaged_data),
        get<package_field_tag<1>>(remote_packaged_data),
        get<package_field_tag<2>>(remote_packaged_data),
        get<package_field_tag<3>>(remote_packaged_data),
        get<package_field_tag<4>>(remote_packaged_data),
        get<package_field_tag<5>>(remote_packaged_data),
        get<package_field_tag<6>>(remote_packaged_data),
        get<package_field_tag<7>>(remote_packaged_data));

    for (size_t face_index = 0; face_index < num_face_points; ++face_index) {
      const size_t volume_index = direction.side() == Side::Upper
                                      ? upper_face_to_volume_index(face_index)
                                      : lower_face_to_volume_index(face_index);
      double normal_magnitude = 0.0;
      for (size_t d = 0; d < Dim; ++d) {
        const double unnormalized_normal_component =
            outward_sign * inverse_jacobian.get(sliced_dim, d)[volume_index];
        normal_magnitude +=
            unnormalized_normal_component * unnormalized_normal_component;
      }
      normal_magnitude = std::sqrt(normal_magnitude);
      const double lifted_factor = lift_prefactor * normal_magnitude;

      get(get<::Tags::dt<ScalarWave::Tags::Psi>>(
          *host_dt_reference))[volume_index] +=
          lifted_factor * get(psi_boundary_correction)[face_index];
      get(get<::Tags::dt<ScalarWave::Tags::Pi>>(
          *host_dt_reference))[volume_index] +=
          lifted_factor * get(pi_boundary_correction)[face_index];
      for (size_t d = 0; d < Dim; ++d) {
        get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(*host_dt_reference)
            .get(d)[volume_index] +=
            lifted_factor * phi_boundary_correction.get(d)[face_index];
      }
    }
  }
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
    outgoing_by_element.insert_or_assign(setup.element.id(),
                                         compute_outgoing_boundary_data(setup));
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

    Variables<dt_variables_tags> host_dt_initial{num_points, 0.0};
    auto device_dt = copy_to_device(host_dt_initial);

    ScalarWave::Actions::ApplyBoundaryCorrectionsToTimeDerivativeKokkos::apply(
        make_not_null(&device_dt), outgoing_by_element.at(setup.element.id()),
        incoming_by_element.at(setup.element.id()),
        external_boundary_data_map{}, setup.device_face_to_volume_index_map,
        setup.device_face_normal_magnitude, setup.device_mortar_data,
        setup.mortar_meshes, setup.mesh, setup.element);

    Variables<dt_variables_tags> host_dt_from_kokkos{num_points, 0.0};
    copy_to_host(make_not_null(&host_dt_from_kokkos), device_dt);

    Variables<dt_variables_tags> host_dt_reference{num_points, 0.0};
    add_host_boundary_corrections(make_not_null(&host_dt_reference),
                                  outgoing_by_element.at(setup.element.id()),
                                  incoming_by_element.at(setup.element.id()),
                                  setup.mortar_meshes, setup.inverse_jacobian,
                                  setup.device_face_to_volume_index_map,
                                  setup.mesh, setup.element);

    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Psi>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Psi>>(host_dt_reference));
    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Pi>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Pi>>(host_dt_reference));
    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(host_dt_reference));
  }
}

void test_apply_boundary_corrections_iterable_action_matches_host() {
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
    outgoing_by_element.insert_or_assign(setup.element.id(),
                                         compute_outgoing_boundary_data(setup));
  }

  std::map<ElementId<Dim>, incoming_boundary_data_map> incoming_by_element{};
  for (const auto& setup : setups) {
    incoming_by_element.insert_or_assign(
        setup.element.id(),
        build_incoming_boundary_data(setup, outgoing_by_element));
  }

  for (const auto& setup : setups) {
    INFO("Checking iterable action for element " << setup.element.id());

    const size_t num_points = setup.mesh.number_of_grid_points();
    Variables<dt_variables_tags> host_dt_initial{num_points, 0.0};
    auto device_dt_initial = copy_to_device(host_dt_initial);
    auto box = db::create<db::AddSimpleTags<
        ScalarWave::KokkosTags::DeviceDtVariables<system>,
        ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>,
        ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<Dim>,
        ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>,
        ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>,
        ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<Dim>,
        ScalarWave::KokkosTags::DeviceMortarData<Dim>,
        evolution::dg::Tags::MortarMesh<Dim>, domain::Tags::Mesh<Dim>,
        domain::Tags::Element<Dim>, ::Tags::TimeStepId>>(
        std::move(device_dt_initial),
        outgoing_by_element.at(setup.element.id()),
        incoming_boundary_data_map{}, external_boundary_data_map{},
        setup.device_face_to_volume_index_map,
        setup.device_face_normal_magnitude, setup.device_mortar_data,
        setup.mortar_meshes, setup.mesh, setup.element, setup.time_step_id);

    tuples::TaggedTuple<inbox_tag> inboxes{};
    const Parallel::GlobalCache<TestMetavariables> cache{
        typename Parallel::GlobalCache<TestMetavariables>::ConstTagsTuple{}};

    const auto first_result = apply_action::apply(
        box, inboxes, cache, static_cast<size_t>(0), tmpl::list<>{},
        static_cast<const DummyParallelComponent*>(nullptr));
    CHECK(std::get<0>(first_result) == Parallel::AlgorithmExecution::Retry);

    auto& inbox = tuples::get<inbox_tag>(inboxes);
    for (const auto& [directional_id, data] :
         incoming_by_element.at(setup.element.id())) {
      inbox_tag::insert_into_inbox(make_not_null(&inbox), setup.time_step_id,
                                   std::make_pair(directional_id, data));
    }

    const auto second_result = apply_action::apply(
        box, inboxes, cache, static_cast<size_t>(0), tmpl::list<>{},
        static_cast<const DummyParallelComponent*>(nullptr));
    CHECK(std::get<0>(second_result) == Parallel::AlgorithmExecution::Continue);
    CHECK(tuples::get<inbox_tag>(inboxes).empty());

    Variables<dt_variables_tags> host_dt_from_kokkos{num_points, 0.0};
    copy_to_host(
        make_not_null(&host_dt_from_kokkos),
        db::get<ScalarWave::KokkosTags::DeviceDtVariables<system>>(box));

    const auto& inbox_to_databox_boundary_data =
        db::get<ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<Dim>>(
            box);
    CHECK(inbox_to_databox_boundary_data.size() ==
          setup.element.number_of_neighbors());

    Variables<dt_variables_tags> host_dt_reference{num_points, 0.0};
    add_host_boundary_corrections(make_not_null(&host_dt_reference),
                                  outgoing_by_element.at(setup.element.id()),
                                  incoming_by_element.at(setup.element.id()),
                                  setup.mortar_meshes, setup.inverse_jacobian,
                                  setup.device_face_to_volume_index_map,
                                  setup.mesh, setup.element);

    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Psi>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Psi>>(host_dt_reference));
    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Pi>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Pi>>(host_dt_reference));
    CHECK_ITERABLE_APPROX(
        get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(host_dt_from_kokkos),
        get<::Tags::dt<ScalarWave::Tags::Phi<Dim>>>(host_dt_reference));
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.ScalarWave."
    "ApplyBoundaryCorrectionsToTimeDerivativeKokkos",
    "[Unit][Evolution]") {
  test_apply_boundary_corrections_to_time_derivative_kokkos_matches_host();
  test_apply_boundary_corrections_iterable_action_matches_host();
}
