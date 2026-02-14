// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Executables/ScalarWave/KokkosBoundaryCommunication.hpp"
#include "Evolution/Executables/ScalarWave/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

// Skeleton action for a Kokkos-native ComputeTimeDerivative path.
// Assumptions for this first draft:
// - global time stepping
// - static mesh
// - no AMR / no subcell
// - same-node communication only
// - boundary correction implementation (UpwindPenalty) not implemented yet
template <size_t Dim, typename System>
struct ComputeTimeDerivativeKokkos {
 private:
  using device_variables_tag = ScalarWave::KokkosTags::DeviceVariables<System>;
  using device_dt_variables_tag =
      ScalarWave::KokkosTags::DeviceDtVariables<System>;
  using device_inverse_jacobian_tag =
      ScalarWave::KokkosTags::DeviceInverseJacobian<Dim>;
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;
  using device_face_to_volume_index_map_tag =
      ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<Dim>;
  using outgoing_boundary_data_tag =
      ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>;
  using volume_time_derivative_terms =
      typename System::compute_volume_time_derivative_terms;
  using package_field_tags =
      typename ScalarWave::BoundaryCorrections::UpwindPenalty<
          Dim>::dg_package_field_tags;
  template <size_t I>
  using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
  using device_package_field_tags =
      db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;

 public:
  using return_tags =
      tmpl::list<device_dt_variables_tag, outgoing_boundary_data_tag>;
  using argument_tags =
      tmpl::list<device_variables_tag, device_inverse_jacobian_tag,
                 device_constraint_gamma2_tag,
                 device_face_to_volume_index_map_tag, domain::Tags::Mesh<Dim>,
                 domain::Tags::Element<Dim>, ::Tags::TimeStepId>;

  static void apply(
      const gsl::not_null<typename device_dt_variables_tag::type*> device_dt,
      const gsl::not_null<typename outgoing_boundary_data_tag::type*>
          outgoing_boundary_data,
      const typename device_variables_tag::type& device_vars,
      const typename device_inverse_jacobian_tag::type& device_inverse_jacobian,
      const typename device_constraint_gamma2_tag::type&
          device_constraint_gamma2,
      const typename device_face_to_volume_index_map_tag::type&
          device_face_to_volume_index_map,
      const Mesh<Dim>& mesh, const Element<Dim>& element,
      const TimeStepId& time_step_id) {
    using device_gradient_tags =
        db::wrap_tags_in<::Tags::MirrorView,
                         typename System::gradient_variables>;
    using device_derivative_tags =
        db::wrap_tags_in<::Tags::deriv, device_gradient_tags, tmpl::size_t<Dim>,
                         Frame::Inertial>;
    Variables<device_derivative_tags> device_partial_derivatives{
        mesh.number_of_grid_points()};

    if (device_dt->number_of_grid_points() != mesh.number_of_grid_points()) {
      device_dt->initialize(mesh.number_of_grid_points());
    }
    partial_derivatives(make_not_null(&device_partial_derivatives), device_vars,
                        mesh, device_inverse_jacobian);

    const auto dt_psi =
        get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(*device_dt);
    const auto dt_pi =
        get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(*device_dt);
    const auto dt_phi =
        get<::Tags::MirrorView<::Tags::dt<Tags::Phi<Dim>>>>(*device_dt);

    Kokkos::parallel_for(
        "ComputeTimeDerivativeKokkosVolumeTerms", mesh.number_of_grid_points(),
        KOKKOS_LAMBDA(const int s) {
          Scalar<double> dt_psi_at_s{};
          Scalar<double> dt_pi_at_s{};
          tnsr::i<double, Dim, Frame::Inertial> dt_phi_at_s{};
          Scalar<double> result_gamma2_at_s{};

          const auto vars_at_s = make_at_index(device_vars, s);
          const auto derivs_at_s = make_at_index(device_partial_derivatives, s);
          const auto gamma2_at_s = make_at_index(device_constraint_gamma2, s);

          const auto decisions = volume_time_derivative_terms::apply(
              make_not_null(&dt_psi_at_s), make_not_null(&dt_pi_at_s),
              make_not_null(&dt_phi_at_s), make_not_null(&result_gamma2_at_s),
              get<::Tags::AtIndex<
                  ::Tags::deriv<::Tags::MirrorView<Tags::Psi>,
                                tmpl::size_t<Dim>, Frame::Inertial>>>(
                  derivs_at_s),
              get<::Tags::AtIndex<
                  ::Tags::deriv<::Tags::MirrorView<Tags::Pi>, tmpl::size_t<Dim>,
                                Frame::Inertial>>>(derivs_at_s),
              get<::Tags::AtIndex<
                  ::Tags::deriv<::Tags::MirrorView<Tags::Phi<Dim>>,
                                tmpl::size_t<Dim>, Frame::Inertial>>>(
                  derivs_at_s),
              get<::Tags::AtIndex<::Tags::MirrorView<Tags::Pi>>>(vars_at_s),
              get<::Tags::AtIndex<::Tags::MirrorView<Tags::Phi<Dim>>>>(
                  vars_at_s),
              gamma2_at_s);
          (void)decisions;
          (void)result_gamma2_at_s;

          get(dt_psi)[s] = get(dt_psi_at_s);
          get(dt_pi)[s] = get(dt_pi_at_s);
          for (size_t d = 0; d < Dim; ++d) {
            dt_phi.get(d)[s] = dt_phi_at_s.get(d);
          }
        });

    outgoing_boundary_data->clear();
    const auto inverse_jacobian = device_inverse_jacobian;

    for (const auto& [direction, neighbors] : element.neighbors()) {
      const size_t sliced_dim = direction.dimension();
      const double outward_sign = direction.side() == Side::Upper ? 1.0 : -1.0;
      const Mesh<Dim - 1> face_mesh = mesh.slice_away(sliced_dim);
      const size_t num_face_points = face_mesh.number_of_grid_points();
      auto face_to_volume_index =
          direction.side() == Side::Upper
              ? gsl::at(device_face_to_volume_index_map, sliced_dim).second
              : gsl::at(device_face_to_volume_index_map, sliced_dim).first;

      Variables<device_package_field_tags> packaged_face_data{num_face_points};

      if (num_face_points > 0) {
        Kokkos::parallel_for(
            "PackageBoundaryCorrectionData", num_face_points,
            KOKKOS_LAMBDA(const int face_index) {
              const size_t volume_index =
                  face_to_volume_index(static_cast<size_t>(face_index));

              std::array<double, Dim> unit_normal_covector{};
              double normal_magnitude_squared = 0.0;
              for (size_t d = 0; d < Dim; ++d) {
                unit_normal_covector[d] =
                    outward_sign * inverse_jacobian.get(sliced_dim, d)[volume_index];
                normal_magnitude_squared +=
                    unit_normal_covector[d] * unit_normal_covector[d];
              }
              const double inverse_normal_magnitude =
                  1.0 / sqrt(normal_magnitude_squared);
              for (size_t d = 0; d < Dim; ++d) {
                unit_normal_covector[d] *= inverse_normal_magnitude;
              }

              tnsr::i<double, Dim, Frame::Inertial> normal_covector{};
              for (size_t d = 0; d < Dim; ++d) {
                normal_covector.get(d) = unit_normal_covector[d];
              }

              const auto vars_at_s = make_at_index(device_vars, volume_index);
              const auto gamma2_at_s =
                  make_at_index(device_constraint_gamma2, volume_index);

              Scalar<double> char_speed_v_psi{};
              tnsr::i<double, Dim, Frame::Inertial> char_speed_v_zero{};
              Scalar<double> char_speed_v_plus{};
              Scalar<double> char_speed_v_minus{};
              tnsr::i<double, Dim, Frame::Inertial> char_speed_n_times_v_plus{};
              tnsr::i<double, Dim, Frame::Inertial> char_speed_n_times_v_minus{};
              Scalar<double> char_speed_gamma2_v_psi{};
              tnsr::i<double, 3, Frame::Inertial> char_speeds{};
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
                  get<::Tags::AtIndex<::Tags::MirrorView<Tags::Phi<Dim>>>>(
                      vars_at_s),
                  gamma2_at_s, normal_covector,
                  static_cast<const Scalar<double>*>(nullptr));

              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<0>>>(packaged_face_data)),
                           char_speed_v_psi, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<1>>>(packaged_face_data)),
                           char_speed_v_zero, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<2>>>(packaged_face_data)),
                           char_speed_v_plus, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<3>>>(packaged_face_data)),
                           char_speed_v_minus, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<4>>>(packaged_face_data)),
                           char_speed_n_times_v_plus, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<5>>>(packaged_face_data)),
                           char_speed_n_times_v_minus, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<6>>>(packaged_face_data)),
                           char_speed_gamma2_v_psi, face_index);
              set_at_index(make_not_null(&get<::Tags::MirrorView<
                               package_field_tag<7>>>(packaged_face_data)),
                           char_speeds, face_index);
            });
      }

      for (const auto& neighbor : neighbors) {
        ScalarWave::KokkosTags::BoundaryCorrectionData<Dim> boundary_data{};
        boundary_data.volume_mesh = mesh;
        boundary_data.boundary_correction_mesh = face_mesh;
        boundary_data.boundary_correction_data = packaged_face_data;
        boundary_data.validity_range = time_step_id;
        boundary_data.integration_order = 0;
        outgoing_boundary_data->insert_or_assign(
            DirectionalId<Dim>{direction, neighbor}, std::move(boundary_data));
      }
    }
  }
};

}  // namespace ScalarWave::Actions
