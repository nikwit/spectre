// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Kokkos/MortarDataPrimitives.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ApplyBoundaryCorrectionsToTimeDerivativeBatched {
 private:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using package_field_tags =
      typename ScalarWave::BoundaryCorrections::UpwindPenalty<
          volume_dim>::dg_package_field_tags;
  template <size_t I>
  using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
  using device_package_field_tags =
      db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;
  using dt_boundary_tags =
      tmpl::list<::Tags::dt<ScalarWave::Tags::Psi>,
                 ::Tags::dt<ScalarWave::Tags::Pi>,
                 ::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>;
  using device_dt_boundary_tags =
      db::wrap_tags_in<::Tags::MirrorView, dt_boundary_tags>;
  static double outward_sign(const Side side) {
    return side == Side::Upper ? 1.0 : -1.0;
  }

 public:
  using return_tags = tmpl::list<ScalarWave::Batched::Tags::DeviceData>;
  using argument_tags =
      tmpl::list<::Tags::Time,
                 domain::Tags::ExternalBoundaryConditions<volume_dim>>;

  static void apply(
      const gsl::not_null<ScalarWave::Batched::Tags::DeviceData::type*>
          device_data,
      const double time,
      const typename domain::Tags::ExternalBoundaryConditions<
          volume_dim>::type& external_boundary_conditions_by_block) {
    if (device_data->total_points() == 0) {
      return;
    }

    ASSERT(device_data->uniform_quadrature_host()[0] ==
               Spectral::Quadrature::GaussLobatto and
               device_data->uniform_quadrature_host()[1] ==
                   Spectral::Quadrature::GaussLobatto and
               device_data->uniform_quadrature_host()[2] ==
                   Spectral::Quadrature::GaussLobatto,
           "ApplyBoundaryCorrectionsToTimeDerivativeBatched currently supports "
           "Gauss-Lobatto quadrature only.");

    const auto& extents = device_data->uniform_extents_host();
    const Mesh<volume_dim> mesh{
        extents, device_data->uniform_basis_host(),
        device_data->uniform_quadrature_host()};
    const auto& face_to_volume_index_map =
        device_data->device_face_to_volume_index_map();
    const auto& elements = device_data->local_elements();
    auto& dt_vars = device_data->device_dt_variables();
    const auto& vars = device_data->device_variables();
    const auto& gamma2 = device_data->device_constraint_gamma2();
    const auto& inverse_jacobian = device_data->element_inverse_jacobian_device();
    const auto host_gamma2 =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, get(gamma2));
    const auto host_inverse_jacobian = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, inverse_jacobian);
    const auto dt_psi =
        get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(dt_vars);
    const auto dt_pi =
        get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(dt_vars);
    const auto dt_phi = get<
        ::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
        dt_vars);

    for (size_t e = 0; e < elements.size(); ++e) {
      const auto& element = elements[e];
      const auto& oriented_remote_face_index_for_local =
          device_data->oriented_remote_face_index_for_local(e);
      const auto& mortar_metadata = device_data->mortar_metadata(e);
      const size_t local_point_offset =
          device_data->element_point_offsets_host()[e];

      for (const auto& [direction, neighbors] : element.neighbors()) {
        const size_t sliced_dim = direction.dimension();
        const Mesh<volume_dim - 1> face_mesh = mesh.slice_away(sliced_dim);
        const size_t num_face_points = face_mesh.number_of_grid_points();
        Variables<device_dt_boundary_tags> dt_boundary_correction_on_face_sum{
            num_face_points};
        ::Kokkos::deep_copy(dt_boundary_correction_on_face_sum.view(), 0.0);
        const auto local_face_to_volume_index =
            direction.side() == Side::Upper
                ? gsl::at(face_to_volume_index_map, sliced_dim).second
                : gsl::at(face_to_volume_index_map, sliced_dim).first;
        const double local_outward_sign = outward_sign(direction.side());

        Variables<device_package_field_tags> local_packaged_face_data{
            num_face_points};
        ::Kokkos::parallel_for(
            "ApplyBoundaryCorrectionsToTimeDerivativeBatchedPackageLocalFace",
            num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
              const size_t face_index = static_cast<size_t>(face_index_int);
              const size_t local_volume_index_in_element =
                  local_face_to_volume_index(face_index);
              const size_t local_volume_index =
                  local_point_offset + local_volume_index_in_element;

              tnsr::i<double, volume_dim, Frame::Inertial> local_normal_covector{};
              double local_normal_magnitude_squared = 0.0;
              for (size_t d = 0; d < volume_dim; ++d) {
                const double local_component =
                    local_outward_sign *
                    inverse_jacobian(e, local_volume_index_in_element,
                                     sliced_dim * volume_dim + d);
                local_normal_covector.get(d) = local_component;
                local_normal_magnitude_squared += local_component * local_component;
              }
              const double local_normal_magnitude =
                  sqrt(local_normal_magnitude_squared);
              for (size_t d = 0; d < volume_dim; ++d) {
                local_normal_covector.get(d) /= local_normal_magnitude;
              }

              const auto local_vars_at_s = make_at_index(vars, local_volume_index);
              const auto local_gamma2_at_s = make_at_index(gamma2, local_volume_index);
              Scalar<double> local_v_psi{};
              tnsr::i<double, volume_dim, Frame::Inertial> local_v_zero{};
              Scalar<double> local_v_plus{};
              Scalar<double> local_v_minus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus{};
              Scalar<double> local_gamma2_v_psi{};
              tnsr::i<double, volume_dim, Frame::Inertial> local_char_speeds{};
              (void)ScalarWave::BoundaryCorrections::detail::dg_package_data_impl(
                  make_not_null(&local_v_psi), make_not_null(&local_v_zero),
                  make_not_null(&local_v_plus), make_not_null(&local_v_minus),
                  make_not_null(&local_normal_times_v_plus),
                  make_not_null(&local_normal_times_v_minus),
                  make_not_null(&local_gamma2_v_psi),
                  make_not_null(&local_char_speeds),
                  get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Psi>>>(
                      local_vars_at_s),
                  get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Pi>>>(
                      local_vars_at_s),
                  get<::Tags::AtIndex<
                      ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>>(
                      local_vars_at_s),
                  local_gamma2_at_s, local_normal_covector,
                  static_cast<const Scalar<double>*>(nullptr));
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<0>>>(
                      local_packaged_face_data)),
                  local_v_psi, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<1>>>(
                      local_packaged_face_data)),
                  local_v_zero, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<2>>>(
                      local_packaged_face_data)),
                  local_v_plus, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<3>>>(
                      local_packaged_face_data)),
                  local_v_minus, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<4>>>(
                      local_packaged_face_data)),
                  local_normal_times_v_plus, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<5>>>(
                      local_packaged_face_data)),
                  local_normal_times_v_minus, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<6>>>(
                      local_packaged_face_data)),
                  local_gamma2_v_psi, face_index);
              set_at_index(
                  make_not_null(&get<::Tags::MirrorView<package_field_tag<7>>>(
                      local_packaged_face_data)),
                  local_char_speeds, face_index);
            });

        for (const auto& neighbor_id : neighbors) {
          const auto& orientation = neighbors.orientation(neighbor_id);
          const size_t remote_element_index =
              device_data->element_index(neighbor_id);
          const size_t remote_point_offset =
              device_data->element_point_offsets_host()[remote_element_index];
          const Direction<volume_dim> remote_direction =
              orientation(direction.opposite());
          const size_t remote_sliced_dim = remote_direction.dimension();
          const DirectionalId<volume_dim> mortar_id{direction, neighbor_id};
          ASSERT(oriented_remote_face_index_for_local.count(mortar_id) == 1,
                 "Missing oriented remote-face index map for " << mortar_id
                                                               << ".");
          ASSERT(mortar_metadata.count(mortar_id) == 1,
                 "Missing mortar metadata for " << mortar_id << ".");
          const auto& mortar = mortar_metadata.at(mortar_id);
          const Mesh<volume_dim - 1>& mortar_mesh = mortar.mortar_mesh;
          const auto& mortar_size = mortar.mortar_size;
          const bool needs_projection = mortar.needs_projection;
          const auto remote_face_to_volume_index =
              remote_direction.side() == Side::Upper
                  ? gsl::at(face_to_volume_index_map, remote_sliced_dim).second
                  : gsl::at(face_to_volume_index_map, remote_sliced_dim).first;
          const double remote_outward_sign = outward_sign(remote_direction.side());
          const auto remote_face_index_for_local_face_index =
              oriented_remote_face_index_for_local.at(mortar_id);

          Variables<device_package_field_tags> remote_packaged_face_data{
              num_face_points};
          ::Kokkos::parallel_for(
              "ApplyBoundaryCorrectionsToTimeDerivativeBatchedPackageRemoteFace",
              num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
                const size_t face_index = static_cast<size_t>(face_index_int);
                const size_t remote_face_index =
                    remote_face_index_for_local_face_index(face_index);
                const size_t remote_volume_index_in_element =
                    remote_face_to_volume_index(remote_face_index);
                const size_t remote_volume_index =
                    remote_point_offset + remote_volume_index_in_element;

                tnsr::i<double, volume_dim, Frame::Inertial> remote_normal_covector{};
                double remote_normal_magnitude_squared = 0.0;
                for (size_t d = 0; d < volume_dim; ++d) {
                  const double remote_component =
                      remote_outward_sign *
                      inverse_jacobian(remote_element_index,
                                       remote_volume_index_in_element,
                                       remote_sliced_dim * volume_dim + d);
                  remote_normal_covector.get(d) = remote_component;
                  remote_normal_magnitude_squared +=
                      remote_component * remote_component;
                }
                const double remote_normal_magnitude =
                    sqrt(remote_normal_magnitude_squared);
                for (size_t d = 0; d < volume_dim; ++d) {
                  remote_normal_covector.get(d) /= remote_normal_magnitude;
                }

                const auto remote_vars_at_s = make_at_index(vars, remote_volume_index);
                const auto remote_gamma2_at_s =
                    make_at_index(gamma2, remote_volume_index);
                Scalar<double> remote_v_psi{};
                tnsr::i<double, volume_dim, Frame::Inertial> remote_v_zero{};
                Scalar<double> remote_v_plus{};
                Scalar<double> remote_v_minus{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_normal_times_v_plus{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_normal_times_v_minus{};
                Scalar<double> remote_gamma2_v_psi{};
                tnsr::i<double, volume_dim, Frame::Inertial> remote_char_speeds{};
                (void)ScalarWave::BoundaryCorrections::detail::dg_package_data_impl(
                    make_not_null(&remote_v_psi), make_not_null(&remote_v_zero),
                    make_not_null(&remote_v_plus), make_not_null(&remote_v_minus),
                    make_not_null(&remote_normal_times_v_plus),
                    make_not_null(&remote_normal_times_v_minus),
                    make_not_null(&remote_gamma2_v_psi),
                    make_not_null(&remote_char_speeds),
                    get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Psi>>>(
                        remote_vars_at_s),
                    get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Pi>>>(
                        remote_vars_at_s),
                    get<::Tags::AtIndex<
                        ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>>(
                        remote_vars_at_s),
                    remote_gamma2_at_s, remote_normal_covector,
                    static_cast<const Scalar<double>*>(nullptr));
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<0>>>(
                        remote_packaged_face_data)),
                    remote_v_psi, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<1>>>(
                        remote_packaged_face_data)),
                    remote_v_zero, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<2>>>(
                        remote_packaged_face_data)),
                    remote_v_plus, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<3>>>(
                        remote_packaged_face_data)),
                    remote_v_minus, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<4>>>(
                        remote_packaged_face_data)),
                    remote_normal_times_v_plus, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<5>>>(
                        remote_packaged_face_data)),
                    remote_normal_times_v_minus, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<6>>>(
                        remote_packaged_face_data)),
                    remote_gamma2_v_psi, face_index);
                set_at_index(
                    make_not_null(&get<::Tags::MirrorView<package_field_tag<7>>>(
                        remote_packaged_face_data)),
                    remote_char_speeds, face_index);
              });

          Variables<device_package_field_tags> local_packaged_mortar_data{
              mortar_mesh.number_of_grid_points()};
          Variables<device_package_field_tags> remote_packaged_mortar_data{
              mortar_mesh.number_of_grid_points()};
          if (needs_projection) {
            ScalarWave::Kokkos::project_to_mortar_device(
                make_not_null(&local_packaged_mortar_data), local_packaged_face_data,
                face_mesh, mortar_mesh, mortar_size);
            ScalarWave::Kokkos::project_to_mortar_device(
                make_not_null(&remote_packaged_mortar_data),
                remote_packaged_face_data, face_mesh, mortar_mesh, mortar_size);
          } else {
            ::Kokkos::deep_copy(local_packaged_mortar_data.view(),
                                local_packaged_face_data.view());
            ::Kokkos::deep_copy(remote_packaged_mortar_data.view(),
                                remote_packaged_face_data.view());
          }

          Variables<device_dt_boundary_tags> dt_boundary_correction_on_mortar{
              mortar_mesh.number_of_grid_points()};
          const auto local_v_psi =
              get<::Tags::MirrorView<package_field_tag<0>>>(
                  local_packaged_mortar_data);
          const auto local_v_zero =
              get<::Tags::MirrorView<package_field_tag<1>>>(
                  local_packaged_mortar_data);
          const auto local_v_plus =
              get<::Tags::MirrorView<package_field_tag<2>>>(
                  local_packaged_mortar_data);
          const auto local_v_minus =
              get<::Tags::MirrorView<package_field_tag<3>>>(
                  local_packaged_mortar_data);
          const auto local_normal_times_v_plus =
              get<::Tags::MirrorView<package_field_tag<4>>>(
                  local_packaged_mortar_data);
          const auto local_normal_times_v_minus =
              get<::Tags::MirrorView<package_field_tag<5>>>(
                  local_packaged_mortar_data);
          const auto local_gamma2_v_psi =
              get<::Tags::MirrorView<package_field_tag<6>>>(
                  local_packaged_mortar_data);
          const auto local_char_speeds =
              get<::Tags::MirrorView<package_field_tag<7>>>(
                  local_packaged_mortar_data);
          const auto remote_v_psi =
              get<::Tags::MirrorView<package_field_tag<0>>>(
                  remote_packaged_mortar_data);
          const auto remote_v_zero =
              get<::Tags::MirrorView<package_field_tag<1>>>(
                  remote_packaged_mortar_data);
          const auto remote_v_plus =
              get<::Tags::MirrorView<package_field_tag<2>>>(
                  remote_packaged_mortar_data);
          const auto remote_v_minus =
              get<::Tags::MirrorView<package_field_tag<3>>>(
                  remote_packaged_mortar_data);
          const auto remote_normal_times_v_plus =
              get<::Tags::MirrorView<package_field_tag<4>>>(
                  remote_packaged_mortar_data);
          const auto remote_normal_times_v_minus =
              get<::Tags::MirrorView<package_field_tag<5>>>(
                  remote_packaged_mortar_data);
          const auto remote_gamma2_v_psi =
              get<::Tags::MirrorView<package_field_tag<6>>>(
                  remote_packaged_mortar_data);
          const auto remote_char_speeds =
              get<::Tags::MirrorView<package_field_tag<7>>>(
                  remote_packaged_mortar_data);
          const auto dt_psi_on_mortar =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                  dt_boundary_correction_on_mortar);
          const auto dt_pi_on_mortar =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                  dt_boundary_correction_on_mortar);
          const auto dt_phi_on_mortar =
              get<::Tags::MirrorView<
                  ::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
                  dt_boundary_correction_on_mortar);

          ::Kokkos::parallel_for(
              "ApplyBoundaryCorrectionsToTimeDerivativeBatchedOnMortar",
              mortar_mesh.number_of_grid_points(),
              KOKKOS_LAMBDA(const int mortar_index_int) {
                const size_t mortar_index = static_cast<size_t>(mortar_index_int);
                tnsr::i<double, volume_dim, Frame::Inertial>
                    local_v_zero_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    local_normal_times_v_plus_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    local_normal_times_v_minus_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    local_char_speeds_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_v_zero_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_normal_times_v_plus_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_normal_times_v_minus_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial>
                    remote_char_speeds_at_mortar{};
                for (size_t d = 0; d < volume_dim; ++d) {
                  local_v_zero_at_mortar.get(d) = local_v_zero.get(d)[mortar_index];
                  local_normal_times_v_plus_at_mortar.get(d) =
                      local_normal_times_v_plus.get(d)[mortar_index];
                  local_normal_times_v_minus_at_mortar.get(d) =
                      local_normal_times_v_minus.get(d)[mortar_index];
                  local_char_speeds_at_mortar.get(d) =
                      local_char_speeds.get(d)[mortar_index];
                  remote_v_zero_at_mortar.get(d) =
                      remote_v_zero.get(d)[mortar_index];
                  remote_normal_times_v_plus_at_mortar.get(d) =
                      remote_normal_times_v_plus.get(d)[mortar_index];
                  remote_normal_times_v_minus_at_mortar.get(d) =
                      remote_normal_times_v_minus.get(d)[mortar_index];
                  remote_char_speeds_at_mortar.get(d) =
                      remote_char_speeds.get(d)[mortar_index];
                }

                Scalar<double> dt_psi_at_mortar{};
                Scalar<double> dt_pi_at_mortar{};
                tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_at_mortar{};
                Scalar<double> local_v_psi_at_mortar{};
                Scalar<double> local_v_plus_at_mortar{};
                Scalar<double> local_v_minus_at_mortar{};
                Scalar<double> local_gamma2_v_psi_at_mortar{};
                Scalar<double> remote_v_psi_at_mortar{};
                Scalar<double> remote_v_plus_at_mortar{};
                Scalar<double> remote_v_minus_at_mortar{};
                Scalar<double> remote_gamma2_v_psi_at_mortar{};
                get(local_v_psi_at_mortar) = get(local_v_psi)[mortar_index];
                get(local_v_plus_at_mortar) = get(local_v_plus)[mortar_index];
                get(local_v_minus_at_mortar) = get(local_v_minus)[mortar_index];
                get(local_gamma2_v_psi_at_mortar) =
                    get(local_gamma2_v_psi)[mortar_index];
                get(remote_v_psi_at_mortar) = get(remote_v_psi)[mortar_index];
                get(remote_v_plus_at_mortar) = get(remote_v_plus)[mortar_index];
                get(remote_v_minus_at_mortar) = get(remote_v_minus)[mortar_index];
                get(remote_gamma2_v_psi_at_mortar) =
                    get(remote_gamma2_v_psi)[mortar_index];
                ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
                    volume_dim, double>(
                    make_not_null(&dt_psi_at_mortar),
                    make_not_null(&dt_pi_at_mortar),
                    make_not_null(&dt_phi_at_mortar), local_v_psi_at_mortar,
                    local_v_zero_at_mortar, local_v_plus_at_mortar,
                    local_v_minus_at_mortar, local_normal_times_v_plus_at_mortar,
                    local_normal_times_v_minus_at_mortar,
                    local_gamma2_v_psi_at_mortar, local_char_speeds_at_mortar,
                    remote_v_psi_at_mortar, remote_v_zero_at_mortar,
                    remote_v_plus_at_mortar, remote_v_minus_at_mortar,
                    remote_normal_times_v_plus_at_mortar,
                    remote_normal_times_v_minus_at_mortar,
                    remote_gamma2_v_psi_at_mortar, remote_char_speeds_at_mortar);
                get(dt_psi_on_mortar)[mortar_index] = get(dt_psi_at_mortar);
                get(dt_pi_on_mortar)[mortar_index] = get(dt_pi_at_mortar);
                for (size_t d = 0; d < volume_dim; ++d) {
                  dt_phi_on_mortar.get(d)[mortar_index] = dt_phi_at_mortar.get(d);
                }
              });

          Variables<device_dt_boundary_tags> dt_boundary_correction_on_face{
              num_face_points};
          if (needs_projection) {
            ScalarWave::Kokkos::project_from_mortar_device(
                make_not_null(&dt_boundary_correction_on_face),
                dt_boundary_correction_on_mortar, face_mesh, mortar_mesh,
                mortar_size);
          } else {
            ::Kokkos::deep_copy(dt_boundary_correction_on_face.view(),
                                dt_boundary_correction_on_mortar.view());
          }
          const auto dt_psi_on_face =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                  dt_boundary_correction_on_face);
          const auto dt_pi_on_face =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                  dt_boundary_correction_on_face);
          const auto dt_phi_on_face =
              get<::Tags::MirrorView<
                  ::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
                  dt_boundary_correction_on_face);
          const auto dt_psi_face_sum =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                  dt_boundary_correction_on_face_sum);
          const auto dt_pi_face_sum =
              get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                  dt_boundary_correction_on_face_sum);
          const auto dt_phi_face_sum =
              get<::Tags::MirrorView<
                  ::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
                  dt_boundary_correction_on_face_sum);
          ::Kokkos::parallel_for(
              "ApplyBoundaryCorrectionsToTimeDerivativeBatchedAccumulateFace",
              num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
                const size_t face_index = static_cast<size_t>(face_index_int);
                get(dt_psi_face_sum)[face_index] += get(dt_psi_on_face)[face_index];
                get(dt_pi_face_sum)[face_index] += get(dt_pi_on_face)[face_index];
                for (size_t d = 0; d < volume_dim; ++d) {
                  dt_phi_face_sum.get(d)[face_index] +=
                      dt_phi_on_face.get(d)[face_index];
                }
              });
        }

        const auto dt_psi_face_sum =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                dt_boundary_correction_on_face_sum);
        const auto dt_pi_face_sum =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                dt_boundary_correction_on_face_sum);
        const auto dt_phi_face_sum =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
                dt_boundary_correction_on_face_sum);
        const double lift_prefactor =
            -0.5 * static_cast<double>(extents[sliced_dim] *
                                       (extents[sliced_dim] - 1));
        ::Kokkos::parallel_for(
            "ApplyBoundaryCorrectionsToTimeDerivativeBatchedLiftFace",
            num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
              const size_t face_index = static_cast<size_t>(face_index_int);
              const size_t local_volume_index_in_element =
                  local_face_to_volume_index(face_index);
              const size_t local_volume_index =
                  local_point_offset + local_volume_index_in_element;
              double local_normal_magnitude_squared = 0.0;
              for (size_t d = 0; d < volume_dim; ++d) {
                const double local_component =
                    local_outward_sign *
                    inverse_jacobian(e, local_volume_index_in_element,
                                     sliced_dim * volume_dim + d);
                local_normal_magnitude_squared += local_component * local_component;
              }
              const double local_normal_magnitude =
                  sqrt(local_normal_magnitude_squared);
              const double lifted_factor =
                  lift_prefactor * local_normal_magnitude;
              get(dt_psi)[local_volume_index] +=
                  lifted_factor * get(dt_psi_face_sum)[face_index];
              get(dt_pi)[local_volume_index] +=
                  lifted_factor * get(dt_pi_face_sum)[face_index];
              for (size_t d = 0; d < volume_dim; ++d) {
                dt_phi.get(d)[local_volume_index] +=
                    lifted_factor * dt_phi_face_sum.get(d)[face_index];
              }
            });
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
              "ApplyBoundaryCorrectionsToTimeDerivativeBatched currently "
              "supports only DirichletAnalytic for external boundaries. "
              "Unsupported boundary condition on direction "
              << direction << ".");
        }

        const size_t sliced_dim = direction.dimension();
        const Mesh<volume_dim - 1> face_mesh = mesh.slice_away(sliced_dim);
        const size_t num_face_points = face_mesh.number_of_grid_points();
        const double local_outward_sign = outward_sign(direction.side());
        const auto local_face_to_volume_index =
            direction.side() == Side::Upper
                ? gsl::at(face_to_volume_index_map, sliced_dim).second
                : gsl::at(face_to_volume_index_map, sliced_dim).first;
        const auto host_face_to_volume_index =
            ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
                                                  local_face_to_volume_index);

        tnsr::I<DataVector, volume_dim, Frame::Inertial> coords_on_face{
            num_face_points};
        Scalar<DataVector> interior_gamma2_on_face{num_face_points};
        tnsr::i<DataVector, volume_dim, Frame::Inertial> interior_normal_covector{
            num_face_points};
        for (size_t face_index = 0; face_index < num_face_points; ++face_index) {
          const size_t local_volume_index_in_element =
              host_face_to_volume_index(face_index);
          const size_t local_volume_index =
              local_point_offset + local_volume_index_in_element;
          double normal_magnitude_squared = 0.0;
          for (size_t d = 0; d < volume_dim; ++d) {
            coords_on_face.get(d)[face_index] =
                device_data->inertial_coordinates_host()[d][local_volume_index];
            const double normal_component =
                local_outward_sign *
                host_inverse_jacobian(e, local_volume_index_in_element,
                                      sliced_dim * volume_dim + d);
            interior_normal_covector.get(d)[face_index] = normal_component;
            normal_magnitude_squared += normal_component * normal_component;
          }
          const double normal_magnitude = std::sqrt(normal_magnitude_squared);
          for (size_t d = 0; d < volume_dim; ++d) {
            interior_normal_covector.get(d)[face_index] /= normal_magnitude;
          }
          get(interior_gamma2_on_face)[face_index] = host_gamma2(local_volume_index);
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

        const Variables<device_package_field_tags> packaged_exterior_face_data =
            copy_to_device(packaged_exterior_face_data_on_host);
        const auto remote_v_psi =
            get<::Tags::MirrorView<package_field_tag<0>>>(
                packaged_exterior_face_data);
        const auto remote_v_zero =
            get<::Tags::MirrorView<package_field_tag<1>>>(
                packaged_exterior_face_data);
        const auto remote_v_plus =
            get<::Tags::MirrorView<package_field_tag<2>>>(
                packaged_exterior_face_data);
        const auto remote_v_minus =
            get<::Tags::MirrorView<package_field_tag<3>>>(
                packaged_exterior_face_data);
        const auto remote_normal_times_v_plus =
            get<::Tags::MirrorView<package_field_tag<4>>>(
                packaged_exterior_face_data);
        const auto remote_normal_times_v_minus =
            get<::Tags::MirrorView<package_field_tag<5>>>(
                packaged_exterior_face_data);
        const auto remote_gamma2_v_psi =
            get<::Tags::MirrorView<package_field_tag<6>>>(
                packaged_exterior_face_data);
        const auto remote_char_speeds =
            get<::Tags::MirrorView<package_field_tag<7>>>(
                packaged_exterior_face_data);
        const double lift_prefactor =
            -0.5 * static_cast<double>(extents[sliced_dim] *
                                       (extents[sliced_dim] - 1));

        ::Kokkos::parallel_for(
            "ApplyBoundaryCorrectionsToTimeDerivativeBatchedExternalFace",
            num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
              const size_t face_index = static_cast<size_t>(face_index_int);
              const size_t local_volume_index_in_element =
                  local_face_to_volume_index(face_index);
              const size_t local_volume_index =
                  local_point_offset + local_volume_index_in_element;

              tnsr::i<double, volume_dim, Frame::Inertial> local_normal_covector{};
              double local_normal_magnitude_squared = 0.0;
              for (size_t d = 0; d < volume_dim; ++d) {
                const double local_component =
                    local_outward_sign *
                    inverse_jacobian(e, local_volume_index_in_element,
                                     sliced_dim * volume_dim + d);
                local_normal_covector.get(d) = local_component;
                local_normal_magnitude_squared += local_component * local_component;
              }
              const double local_normal_magnitude =
                  sqrt(local_normal_magnitude_squared);
              for (size_t d = 0; d < volume_dim; ++d) {
                local_normal_covector.get(d) /= local_normal_magnitude;
              }

              const auto local_vars_at_s = make_at_index(vars, local_volume_index);
              const auto local_gamma2_at_s = make_at_index(gamma2, local_volume_index);

              Scalar<double> local_v_psi{};
              tnsr::i<double, volume_dim, Frame::Inertial> local_v_zero{};
              Scalar<double> local_v_plus{};
              Scalar<double> local_v_minus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus{};
              Scalar<double> local_gamma2_v_psi{};
              tnsr::i<double, volume_dim, Frame::Inertial> local_char_speeds{};
              (void)ScalarWave::BoundaryCorrections::detail::dg_package_data_impl(
                  make_not_null(&local_v_psi), make_not_null(&local_v_zero),
                  make_not_null(&local_v_plus), make_not_null(&local_v_minus),
                  make_not_null(&local_normal_times_v_plus),
                  make_not_null(&local_normal_times_v_minus),
                  make_not_null(&local_gamma2_v_psi),
                  make_not_null(&local_char_speeds),
                  get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Psi>>>(
                      local_vars_at_s),
                  get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Pi>>>(
                      local_vars_at_s),
                  get<::Tags::AtIndex<
                      ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>>(
                      local_vars_at_s),
                  local_gamma2_at_s, local_normal_covector,
                  static_cast<const Scalar<double>*>(nullptr));

              tnsr::i<double, volume_dim, Frame::Inertial> remote_v_zero_at_face{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_plus_at_face{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_minus_at_face{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_char_speeds_at_face{};
              for (size_t d = 0; d < volume_dim; ++d) {
                remote_v_zero_at_face.get(d) = remote_v_zero.get(d)[face_index];
                remote_normal_times_v_plus_at_face.get(d) =
                    remote_normal_times_v_plus.get(d)[face_index];
                remote_normal_times_v_minus_at_face.get(d) =
                    remote_normal_times_v_minus.get(d)[face_index];
                remote_char_speeds_at_face.get(d) =
                    remote_char_speeds.get(d)[face_index];
              }

              Scalar<double> remote_v_psi_at_face{};
              Scalar<double> remote_v_plus_at_face{};
              Scalar<double> remote_v_minus_at_face{};
              Scalar<double> remote_gamma2_v_psi_at_face{};
              get(remote_v_psi_at_face) = get(remote_v_psi)[face_index];
              get(remote_v_plus_at_face) = get(remote_v_plus)[face_index];
              get(remote_v_minus_at_face) = get(remote_v_minus)[face_index];
              get(remote_gamma2_v_psi_at_face) =
                  get(remote_gamma2_v_psi)[face_index];

              Scalar<double> dt_psi_face{};
              Scalar<double> dt_pi_face{};
              tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_face{};
              ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
                  volume_dim, double>(
                  make_not_null(&dt_psi_face), make_not_null(&dt_pi_face),
                  make_not_null(&dt_phi_face), local_v_psi, local_v_zero,
                  local_v_plus, local_v_minus, local_normal_times_v_plus,
                  local_normal_times_v_minus, local_gamma2_v_psi,
                  local_char_speeds, remote_v_psi_at_face, remote_v_zero_at_face,
                  remote_v_plus_at_face, remote_v_minus_at_face,
                  remote_normal_times_v_plus_at_face,
                  remote_normal_times_v_minus_at_face,
                  remote_gamma2_v_psi_at_face, remote_char_speeds_at_face);

              const double lifted_factor =
                  lift_prefactor * local_normal_magnitude;
              get(dt_psi)[local_volume_index] +=
                  lifted_factor * get(dt_psi_face);
              get(dt_pi)[local_volume_index] +=
                  lifted_factor * get(dt_pi_face);
              for (size_t d = 0; d < volume_dim; ++d) {
                dt_phi.get(d)[local_volume_index] +=
                    lifted_factor * dt_phi_face.get(d);
              }
            });
      }
    }
  }
};

}  // namespace ScalarWave::Actions
