// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/ComputeTimeDerivativeKokkos.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Actions {

void ComputeTimeDerivativeKokkos::orient_boundary_data_for_send(
    const gsl::not_null<Variables<device_package_field_tags>*>
        oriented_boundary_data,
    const Variables<device_package_field_tags>& boundary_data,
    const ::Kokkos::View<size_t*>& oriented_mortar_grid_point_source_index) {
  ScalarWave::Kokkos::orient_each_component_device(
      oriented_boundary_data, boundary_data,
      oriented_mortar_grid_point_source_index);
}

void ComputeTimeDerivativeKokkos::
    compute_volume_terms_and_package_boundary_data(
        const gsl::not_null<device_dt_type*> device_dt,
        const gsl::not_null<outgoing_boundary_data_type*>
            outgoing_boundary_data,
        const gsl::not_null<external_boundary_data_type*>
            external_boundary_data,
        const device_variables_type& device_vars,
        const device_inverse_jacobian_type& device_inverse_jacobian,
        const device_constraint_gamma2_type& device_constraint_gamma2,
        const device_face_to_volume_index_map_type&
            device_face_to_volume_index_map,
        const device_face_unit_normal_covector_type&
            device_face_unit_normal_covector,
        const device_mortar_data_type& device_mortar_data,
        const typename mortar_mesh_tag::type& mortar_meshes,
        const ScalarWave::Tags::ConstraintGamma2::type& host_constraint_gamma2,
        const typename domain::Tags::ExternalBoundaryConditions<
            volume_dim>::type& external_boundary_conditions_by_block,
        const tnsr::I<DataVector, volume_dim, Frame::Inertial>&
            inertial_coordinates,
        const double time, const Mesh<volume_dim>& mesh,
        const Element<volume_dim>& element, const TimeStepId& time_step_id) {
  using device_gradient_tags =
      db::wrap_tags_in<::Tags::MirrorView, typename System::gradient_variables>;
  using device_derivative_tags =
      db::wrap_tags_in<::Tags::deriv, device_gradient_tags,
                       tmpl::size_t<volume_dim>, Frame::Inertial>;

  Variables<device_derivative_tags> device_partial_derivatives{
      mesh.number_of_grid_points()};

  if (device_dt->number_of_grid_points() != mesh.number_of_grid_points()) {
    device_dt->initialize(mesh.number_of_grid_points());
  }

  partial_derivatives(make_not_null(&device_partial_derivatives), device_vars,
                      mesh, device_inverse_jacobian);

  const auto dt_psi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(*device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(*device_dt);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(*device_dt);
  static constexpr size_t Dim = 3;

  ::Kokkos::parallel_for(
      "ComputeTimeDerivativeKokkosVolumeTerms", mesh.number_of_grid_points(),
      KOKKOS_LAMBDA(const int s) {
        Scalar<double> dt_psi_at_s{};
        Scalar<double> dt_pi_at_s{};
        tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_at_s{};
        Scalar<double> result_gamma2_at_s{};

        const auto vars_at_s = make_at_index(device_vars, s);
        const auto derivs_at_s = make_at_index(device_partial_derivatives, s);
        const auto gamma2_at_s = make_at_index(device_constraint_gamma2, s);

        const auto decisions = volume_time_derivative_terms::apply(
            make_not_null(&dt_psi_at_s), make_not_null(&dt_pi_at_s),
            make_not_null(&dt_phi_at_s), make_not_null(&result_gamma2_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<Tags::Psi>, tmpl::size_t<Dim>,
                              Frame::Inertial>>>(derivs_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<Tags::Pi>, tmpl::size_t<Dim>,
                              Frame::Inertial>>>(derivs_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<Tags::Phi<Dim>>,
                              tmpl::size_t<Dim>, Frame::Inertial>>>(
                derivs_at_s),
            get<::Tags::AtIndex<::Tags::MirrorView<Tags::Pi>>>(vars_at_s),
            get<::Tags::AtIndex<::Tags::MirrorView<Tags::Phi<Dim>>>>(vars_at_s),
            gamma2_at_s);
        (void)decisions;
        (void)result_gamma2_at_s;

        get(dt_psi)[s] = get(dt_psi_at_s);
        get(dt_pi)[s] = get(dt_pi_at_s);
        for (size_t d = 0; d < volume_dim; ++d) {
          dt_phi.get(d)[s] = dt_phi_at_s.get(d);
        }
      });

  ASSERT(mesh.quadrature(0) == Spectral::Quadrature::GaussLobatto,
         "ComputeTimeDerivativeKokkos currently supports Gauss-Lobatto "
         "quadrature only.");

  const auto package_boundary_data_on_face = [&](const Direction<volume_dim>&
                                                     direction) {
    const size_t sliced_dim = direction.dimension();
    const size_t num_face_points =
        mesh.slice_away(sliced_dim).number_of_grid_points();
    auto face_to_volume_index =
        direction.side() == Side::Upper
            ? gsl::at(device_face_to_volume_index_map, sliced_dim).second
            : gsl::at(device_face_to_volume_index_map, sliced_dim).first;
    auto face_unit_normal_covector =
        direction.side() == Side::Upper
            ? gsl::at(device_face_unit_normal_covector, sliced_dim).second
            : gsl::at(device_face_unit_normal_covector, sliced_dim).first;

    Variables<device_package_field_tags> packaged_face_data{num_face_points};
    if (num_face_points > 0) {
      ::Kokkos::parallel_for(
          "PackageBoundaryCorrectionDataOnFace", num_face_points,
          KOKKOS_LAMBDA(const int face_index) {
            const size_t volume_index =
                face_to_volume_index(static_cast<size_t>(face_index));

            tnsr::i<double, volume_dim, Frame::Inertial> normal_covector{};
            for (size_t d = 0; d < volume_dim; ++d) {
              normal_covector.get(d) =
                  face_unit_normal_covector(static_cast<size_t>(face_index), d);
            }

            const auto vars_at_s = make_at_index(device_vars, volume_index);
            const auto gamma2_at_s =
                make_at_index(device_constraint_gamma2, volume_index);

            Scalar<double> char_speed_v_psi{};
            tnsr::i<double, volume_dim, Frame::Inertial> char_speed_v_zero{};
            Scalar<double> char_speed_v_plus{};
            Scalar<double> char_speed_v_minus{};
            tnsr::i<double, volume_dim, Frame::Inertial>
                char_speed_n_times_v_plus{};
            tnsr::i<double, volume_dim, Frame::Inertial>
                char_speed_n_times_v_minus{};
            Scalar<double> char_speed_gamma2_v_psi{};
            tnsr::i<double, volume_dim, Frame::Inertial> char_speeds{};
            (void)ScalarWave::BoundaryCorrections::detail::dg_package_data_impl(
                make_not_null(&char_speed_v_psi),
                make_not_null(&char_speed_v_zero),
                make_not_null(&char_speed_v_plus),
                make_not_null(&char_speed_v_minus),
                make_not_null(&char_speed_n_times_v_plus),
                make_not_null(&char_speed_n_times_v_minus),
                make_not_null(&char_speed_gamma2_v_psi),
                make_not_null(&char_speeds),
                get<::Tags::AtIndex<::Tags::MirrorView<Tags::Psi>>>(vars_at_s),
                get<::Tags::AtIndex<::Tags::MirrorView<Tags::Pi>>>(vars_at_s),
                get<::Tags::AtIndex<::Tags::MirrorView<Tags::Phi<volume_dim>>>>(
                    vars_at_s),
                gamma2_at_s, normal_covector,
                static_cast<const Scalar<double>*>(nullptr));

            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<0>>>(
                    packaged_face_data)),
                char_speed_v_psi, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<1>>>(
                    packaged_face_data)),
                char_speed_v_zero, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<2>>>(
                    packaged_face_data)),
                char_speed_v_plus, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<3>>>(
                    packaged_face_data)),
                char_speed_v_minus, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<4>>>(
                    packaged_face_data)),
                char_speed_n_times_v_plus, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<5>>>(
                    packaged_face_data)),
                char_speed_n_times_v_minus, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<6>>>(
                    packaged_face_data)),
                char_speed_gamma2_v_psi, face_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<7>>>(
                    packaged_face_data)),
                char_speeds, face_index);
          });
    }
    return packaged_face_data;
  };

  outgoing_boundary_data->clear();
  external_boundary_data->clear();

  for (const auto& [direction, neighbors] : element.neighbors()) {
    const size_t sliced_dim = direction.dimension();
    const Mesh<volume_dim - 1> face_mesh = mesh.slice_away(sliced_dim);
    auto packaged_face_data = package_boundary_data_on_face(direction);

    for (const auto& neighbor : neighbors) {
      const DirectionalId<volume_dim> mortar_id{direction, neighbor};
      const auto& mortar_mesh = mortar_meshes.at(mortar_id);
      const auto& mortar_data = device_mortar_data.at(mortar_id);
      Variables<device_package_field_tags> packaged_mortar_data{
          mortar_mesh.number_of_grid_points()};
      if (mortar_data.needs_projection) {
        ScalarWave::Kokkos::project_to_mortar_device(
            make_not_null(&packaged_mortar_data), packaged_face_data, face_mesh,
            mortar_mesh, mortar_data.mortar_size);
      } else {
        ASSERT(
            packaged_mortar_data.number_of_grid_points() ==
                packaged_face_data.number_of_grid_points(),
            "Expected identical face and mortar point counts when projection "
            "is not needed.");
        ::Kokkos::deep_copy(packaged_mortar_data.view(),
                            packaged_face_data.view());
      }

      ScalarWave::KokkosTags::BoundaryCorrectionData<volume_dim>
          boundary_data{};
      boundary_data.boundary_correction_data = std::move(packaged_mortar_data);
      boundary_data.validity_range = time_step_id;
      boundary_data.integration_order = 0;
      outgoing_boundary_data->insert_or_assign(mortar_id,
                                               std::move(boundary_data));
    }
  }

  if (element.external_boundaries().empty()) {
    return;
  }

  const auto& external_boundary_conditions =
      external_boundary_conditions_by_block.at(element.id().block_id());
  const ScalarWave::BoundaryCorrections::UpwindPenalty<volume_dim>
      boundary_correction{};
  const std::optional<tnsr::I<DataVector, volume_dim, Frame::Inertial>>
      face_mesh_velocity{std::nullopt};
  const std::optional<Scalar<DataVector>> normal_dot_mesh_velocity{
      std::nullopt};

  for (const Direction<volume_dim>& direction : element.external_boundaries()) {
    const auto& boundary_condition_base =
        *external_boundary_conditions.at(direction);
    const auto* const dirichlet_analytic = dynamic_cast<
        const ScalarWave::BoundaryConditions::DirichletAnalytic<volume_dim>*>(
        &boundary_condition_base);
    if (dirichlet_analytic == nullptr) {
      ERROR(
          "ComputeTimeDerivativeKokkos currently supports only "
          "DirichletAnalytic for external boundaries. Unsupported boundary "
          "condition on direction "
          << direction << ".");
    }

    const size_t sliced_dim = direction.dimension();
    const Mesh<volume_dim - 1> face_mesh = mesh.slice_away(sliced_dim);
    const size_t num_face_points = face_mesh.number_of_grid_points();

    auto face_to_volume_index =
        direction.side() == Side::Upper
            ? gsl::at(device_face_to_volume_index_map, sliced_dim).second
            : gsl::at(device_face_to_volume_index_map, sliced_dim).first;
    auto face_unit_normal_covector =
        direction.side() == Side::Upper
            ? gsl::at(device_face_unit_normal_covector, sliced_dim).second
            : gsl::at(device_face_unit_normal_covector, sliced_dim).first;
    // Gather host data for the ghost boundary condition.
    const auto host_face_to_volume_index =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
                                              face_to_volume_index);
    const auto host_face_unit_normal_covector =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
                                              face_unit_normal_covector);
    tnsr::I<DataVector, volume_dim, Frame::Inertial> coords_on_face{
        num_face_points};
    Scalar<DataVector> interior_gamma2_on_face{num_face_points};
    tnsr::i<DataVector, volume_dim, Frame::Inertial> interior_normal_covector{
        num_face_points};
    for (size_t face_index = 0; face_index < num_face_points; ++face_index) {
      const size_t volume_index = host_face_to_volume_index(face_index);
      for (size_t d = 0; d < volume_dim; ++d) {
        coords_on_face.get(d)[face_index] =
            inertial_coordinates.get(d)[volume_index];
        interior_normal_covector.get(d)[face_index] =
            host_face_unit_normal_covector(face_index, d);
      }
      get(interior_gamma2_on_face)[face_index] =
          get(host_constraint_gamma2)[volume_index];
    }

    Scalar<DataVector> exterior_psi{num_face_points};
    Scalar<DataVector> exterior_pi{num_face_points};
    tnsr::i<DataVector, volume_dim, Frame::Inertial> exterior_phi{
        num_face_points};
    Scalar<DataVector> exterior_gamma2{num_face_points};
    const auto error_message = dirichlet_analytic->dg_ghost(
        make_not_null(&exterior_psi), make_not_null(&exterior_pi),
        make_not_null(&exterior_phi), make_not_null(&exterior_gamma2),
        face_mesh_velocity, interior_normal_covector, coords_on_face,
        interior_gamma2_on_face, time);
    if (error_message.has_value()) {
      ERROR(*error_message << "\n\nIn element: " << element.id()
                           << "\nIn direction: " << direction);
    }

    tnsr::i<DataVector, volume_dim, Frame::Inertial> exterior_normal_covector{
        num_face_points};
    for (size_t d = 0; d < volume_dim; ++d) {
      exterior_normal_covector.get(d) = -interior_normal_covector.get(d);
    }

    Variables<package_field_tags> packaged_exterior_face_data_on_host{
        num_face_points};
    (void)boundary_correction.dg_package_data(
        make_not_null(
            &get<package_field_tag<0>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<1>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<2>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<3>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<4>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<5>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<6>>(packaged_exterior_face_data_on_host)),
        make_not_null(
            &get<package_field_tag<7>>(packaged_exterior_face_data_on_host)),
        exterior_psi, exterior_pi, exterior_phi, exterior_gamma2,
        exterior_normal_covector, face_mesh_velocity, normal_dot_mesh_velocity);

    auto packaged_interior_face_data = package_boundary_data_on_face(direction);
    Variables<device_package_field_tags> packaged_exterior_face_data =
        copy_to_device(packaged_exterior_face_data_on_host);
    const DirectionalId<volume_dim> external_mortar_id{
        direction, ElementId<volume_dim>::external_boundary_id()};

    ScalarWave::KokkosTags::BoundaryCorrectionData<volume_dim>
        local_boundary_data{};
    local_boundary_data.boundary_correction_data =
        std::move(packaged_interior_face_data);
    local_boundary_data.validity_range = time_step_id;
    local_boundary_data.integration_order = 0;
    outgoing_boundary_data->insert_or_assign(external_mortar_id,
                                             std::move(local_boundary_data));

    ScalarWave::KokkosTags::BoundaryCorrectionData<volume_dim>
        external_face_data{};
    external_face_data.boundary_correction_data =
        std::move(packaged_exterior_face_data);
    external_face_data.validity_range = time_step_id;
    external_face_data.integration_order = 0;
    external_boundary_data->insert_or_assign(external_mortar_id,
                                             std::move(external_face_data));
  }
}

}  // namespace ScalarWave::Actions
