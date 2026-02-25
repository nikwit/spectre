// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched.hpp"

#include <cmath>
#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {
namespace {

static constexpr size_t volume_dim = 3;
using exterior_primitive_tags =
    tmpl::list<ScalarWave::Tags::Psi, ScalarWave::Tags::Pi,
               ScalarWave::Tags::Phi<volume_dim>,
               ScalarWave::Tags::ConstraintGamma2>;
using device_exterior_primitive_tags =
    db::wrap_tags_in<::Tags::MirrorView, exterior_primitive_tags>;

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
}

}  // namespace

void ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched::apply(
    const gsl::not_null<ScalarWave::Batched::Tags::PackedEvolutionState::type*>
        packed_evolution_state,
    const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
    const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry,
    const double time,
    const typename domain::Tags::ExternalBoundaryConditions<volume_dim>::type&
        external_boundary_conditions_by_block) {
  if (packed_topology.total_points == 0) {
    return;
  }

  ASSERT(packed_topology.uniform_quadrature_host[0] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[1] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[2] ==
                 Spectral::Quadrature::GaussLobatto,
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched currently "
         "supports Gauss-Lobatto quadrature only.");

  const auto& extents = packed_topology.uniform_extents_host;
  ASSERT(extents[0] == extents[1] and extents[1] == extents[2],
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched assumes "
         "uniform p across dimensions (equal extents in all dimensions).");
  const Mesh<volume_dim> mesh{extents, packed_topology.uniform_basis_host,
                              packed_topology.uniform_quadrature_host};
  const Mesh<volume_dim - 1> uniform_face_mesh = mesh.slice_away(0);
  const size_t num_face_points = uniform_face_mesh.number_of_grid_points();

  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& elements = packed_topology.local_elements;
  const auto& element_point_offsets_host =
      packed_topology.element_point_offsets_host;
  auto& dt_vars = packed_evolution_state->device_dt_variables;
  const auto& vars = packed_evolution_state->device_variables;
  const auto& gamma2 = packed_evolution_state->device_constraint_gamma2;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  const auto& element_point_offsets_device =
      packed_topology.element_point_offsets_device;
  const auto host_gamma2 =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, get(gamma2));
  const auto host_inverse_jacobian = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, inverse_jacobian);
  const auto dt_psi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(dt_vars);
  const auto dt_pi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(dt_vars);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
          dt_vars);
  const std::optional<tnsr::I<DataVector, volume_dim, Frame::Inertial>>
      face_mesh_velocity{std::nullopt};
  const std::optional<Scalar<DataVector>> normal_dot_mesh_velocity{
      std::nullopt};

  const size_t num_elements = elements.size();
  for (size_t sliced_dim = 0; sliced_dim < volume_dim; ++sliced_dim) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const Side side = side_i == 0 ? Side::Lower : Side::Upper;
      const Direction<volume_dim> direction{sliced_dim, side};
      const double local_outward_sign = outward_sign(side);
      const auto local_face_to_volume_index =
          side == Side::Upper
              ? gsl::at(face_to_volume_index_map, sliced_dim).second
              : gsl::at(face_to_volume_index_map, sliced_dim).first;
      const auto host_face_to_volume_index =
          ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
                                                local_face_to_volume_index);
      const double lift_prefactor =
          -0.5 *
          static_cast<double>(extents[sliced_dim] * (extents[sliced_dim] - 1));

      Variables<exterior_primitive_tags> exterior_data_on_host{num_elements *
                                                               num_face_points};
      auto& exterior_psi_all =
          get<ScalarWave::Tags::Psi>(exterior_data_on_host);
      auto& exterior_pi_all = get<ScalarWave::Tags::Pi>(exterior_data_on_host);
      auto& exterior_phi_all =
          get<ScalarWave::Tags::Phi<volume_dim>>(exterior_data_on_host);
      auto& exterior_gamma2_all =
          get<ScalarWave::Tags::ConstraintGamma2>(exterior_data_on_host);

      ::Kokkos::View<int*> external_face_mask{"ExternalFaceMask", num_elements};
      auto host_external_face_mask =
          ::Kokkos::create_mirror_view(external_face_mask);
      for (size_t e = 0; e < num_elements; ++e) {
        host_external_face_mask(e) = 0;
      }
      size_t num_external_elements = 0;

      for (size_t e = 0; e < num_elements; ++e) {
        const auto& element = elements[e];
        if (element.external_boundaries().count(direction) == 0) {
          continue;
        }
        host_external_face_mask(e) = 1;
        ++num_external_elements;

        const auto& external_boundary_conditions =
            external_boundary_conditions_by_block.at(element.id().block_id());
        const auto& boundary_condition_base =
            *external_boundary_conditions.at(direction);
        const auto* const dirichlet_analytic =
            dynamic_cast<const ScalarWave::BoundaryConditions::
                             DirichletAnalytic<volume_dim>*>(
                &boundary_condition_base);
        if (dirichlet_analytic == nullptr) {
          ERROR(
              "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched "
              "currently "
              "supports only DirichletAnalytic for external boundaries. "
              "Unsupported boundary condition on direction "
              << direction << ".");
        }

        const size_t local_point_offset = element_point_offsets_host[e];
        tnsr::I<DataVector, volume_dim, Frame::Inertial> coords_on_face{
            num_face_points};
        Scalar<DataVector> interior_gamma2_on_face{num_face_points};
        tnsr::i<DataVector, volume_dim, Frame::Inertial>
            interior_normal_covector{num_face_points};
        for (size_t face_index = 0; face_index < num_face_points;
             ++face_index) {
          const size_t local_volume_index_in_element =
              host_face_to_volume_index(face_index);
          const size_t local_volume_index =
              local_point_offset + local_volume_index_in_element;
          double normal_magnitude_squared = 0.0;
          for (size_t d = 0; d < volume_dim; ++d) {
            coords_on_face.get(d)[face_index] =
                packed_geometry
                    .inertial_coordinates_host[d][local_volume_index];
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
          get(interior_gamma2_on_face)[face_index] =
              host_gamma2(local_volume_index);
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

        for (size_t face_index = 0; face_index < num_face_points;
             ++face_index) {
          const size_t linear_face_index = e * num_face_points + face_index;
          get(exterior_psi_all)[linear_face_index] =
              get(exterior_psi)[face_index];
          get(exterior_pi_all)[linear_face_index] =
              get(exterior_pi)[face_index];
          for (size_t d = 0; d < volume_dim; ++d) {
            exterior_phi_all.get(d)[linear_face_index] =
                exterior_phi.get(d)[face_index];
          }
          get(exterior_gamma2_all)[linear_face_index] =
              get(exterior_gamma2)[face_index];
        }
      }

      if (num_external_elements == 0) {
        continue;
      }

      ::Kokkos::deep_copy(external_face_mask, host_external_face_mask);
      const Variables<device_exterior_primitive_tags> exterior_data_on_device =
          copy_to_device(exterior_data_on_host);
      const auto exterior_psi = get<::Tags::MirrorView<ScalarWave::Tags::Psi>>(
          exterior_data_on_device);
      const auto exterior_pi = get<::Tags::MirrorView<ScalarWave::Tags::Pi>>(
          exterior_data_on_device);
      const auto exterior_phi =
          get<::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>(
              exterior_data_on_device);
      const auto exterior_gamma2 =
          get<::Tags::MirrorView<ScalarWave::Tags::ConstraintGamma2>>(
              exterior_data_on_device);

      ::Kokkos::parallel_for(
          "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatchedExternalFace",
          num_elements * num_face_points,
          KOKKOS_LAMBDA(const int linear_index_int) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t face_index = linear_index % num_face_points;
            const size_t element_index = linear_index / num_face_points;
            if (external_face_mask(element_index) == 0) {
              return;
            }

            const size_t local_volume_index_in_element =
                local_face_to_volume_index(face_index);
            const size_t local_volume_index =
                element_point_offsets_device(element_index) +
                local_volume_index_in_element;

            tnsr::i<double, volume_dim, Frame::Inertial>
                local_normal_covector{};
            double local_normal_magnitude_squared = 0.0;
            for (size_t d = 0; d < volume_dim; ++d) {
              const double local_component =
                  local_outward_sign *
                  inverse_jacobian(element_index, local_volume_index_in_element,
                                   sliced_dim * volume_dim + d);
              local_normal_covector.get(d) = local_component;
              local_normal_magnitude_squared +=
                  local_component * local_component;
            }
            const double local_normal_magnitude =
                sqrt(local_normal_magnitude_squared);
            for (size_t d = 0; d < volume_dim; ++d) {
              local_normal_covector.get(d) /= local_normal_magnitude;
            }

            const auto local_vars_at_s =
                make_at_index(vars, local_volume_index);
            const auto local_gamma2_at_s =
                make_at_index(gamma2, local_volume_index);

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

            const size_t remote_linear_face_index =
                element_index * num_face_points + face_index;
            Scalar<double> exterior_psi_at_face{};
            Scalar<double> exterior_pi_at_face{};
            tnsr::i<double, volume_dim, Frame::Inertial> exterior_phi_at_face{};
            Scalar<double> exterior_gamma2_at_face{};
            get(exterior_psi_at_face) =
                get(exterior_psi)[remote_linear_face_index];
            get(exterior_pi_at_face) =
                get(exterior_pi)[remote_linear_face_index];
            for (size_t d = 0; d < volume_dim; ++d) {
              exterior_phi_at_face.get(d) =
                  exterior_phi.get(d)[remote_linear_face_index];
            }
            get(exterior_gamma2_at_face) =
                get(exterior_gamma2)[remote_linear_face_index];

            tnsr::i<double, volume_dim, Frame::Inertial>
                exterior_normal_covector{};
            for (size_t d = 0; d < volume_dim; ++d) {
              exterior_normal_covector.get(d) = -local_normal_covector.get(d);
            }
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
                make_not_null(&remote_char_speeds), exterior_psi_at_face,
                exterior_pi_at_face, exterior_phi_at_face,
                exterior_gamma2_at_face, exterior_normal_covector,
                static_cast<const Scalar<double>*>(nullptr));

            Scalar<double> dt_psi_face{};
            Scalar<double> dt_pi_face{};
            tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_face{};
            ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
                volume_dim, double>(
                make_not_null(&dt_psi_face), make_not_null(&dt_pi_face),
                make_not_null(&dt_phi_face), local_v_psi, local_v_zero,
                local_v_plus, local_v_minus, local_normal_times_v_plus,
                local_normal_times_v_minus, local_gamma2_v_psi,
                local_char_speeds, remote_v_psi, remote_v_zero, remote_v_plus,
                remote_v_minus, remote_normal_times_v_plus,
                remote_normal_times_v_minus, remote_gamma2_v_psi,
                remote_char_speeds);

            const double lifted_factor =
                lift_prefactor * local_normal_magnitude;
            get(dt_psi)[local_volume_index] += lifted_factor * get(dt_psi_face);
            get(dt_pi)[local_volume_index] += lifted_factor * get(dt_pi_face);
            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi.get(d)[local_volume_index] +=
                  lifted_factor * dt_phi_face.get(d);
            }
          });
    }
  }
}

}  // namespace ScalarWave::Actions
