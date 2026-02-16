// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Kokkos/MortarDataPrimitives.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Actions {
namespace {

static constexpr size_t volume_dim = 3;
using system = ScalarWave::System<volume_dim>;
using boundary_correction_data_type =
    ScalarWave::KokkosTags::BoundaryCorrectionData<volume_dim>;
using device_dt_type = ScalarWave::KokkosTags::DeviceDtVariables<system>::type;
using device_face_to_volume_index_map_type =
    ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<volume_dim>::type;
using device_face_normal_magnitude_type =
    ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<volume_dim>::type;
using package_field_tags =
    typename ScalarWave::BoundaryCorrections::UpwindPenalty<
        volume_dim>::dg_package_field_tags;
template <size_t I>
using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
using dt_boundary_tags = tmpl::list<::Tags::dt<Tags::Psi>, ::Tags::dt<Tags::Pi>,
                                    ::Tags::dt<Tags::Phi<volume_dim>>>;
using device_dt_boundary_tags =
    db::wrap_tags_in<::Tags::MirrorView, dt_boundary_tags>;

void accumulate_mortar_pair_to_face_sum(
    const gsl::not_null<Variables<device_dt_boundary_tags>*>
        dt_boundary_correction_on_face_sum,
    const boundary_correction_data_type& local_boundary_data,
    const boundary_correction_data_type& remote_boundary_data,
    const Mesh<2>& face_mesh, const Mesh<2>& mortar_mesh,
    const bool needs_projection,
    const std::array<Spectral::SegmentSize, volume_dim - 1>& mortar_size) {
  ASSERT(local_boundary_data.boundary_correction_data.number_of_grid_points() ==
             mortar_mesh.number_of_grid_points(),
         "Local packaged data size does not match mortar mesh.");
  ASSERT(
      remote_boundary_data.boundary_correction_data.number_of_grid_points() ==
          mortar_mesh.number_of_grid_points(),
      "Remote packaged data size does not match mortar mesh.");

  const auto local_v_psi = get<::Tags::MirrorView<package_field_tag<0>>>(
      local_boundary_data.boundary_correction_data);
  const auto local_v_zero = get<::Tags::MirrorView<package_field_tag<1>>>(
      local_boundary_data.boundary_correction_data);
  const auto local_v_plus = get<::Tags::MirrorView<package_field_tag<2>>>(
      local_boundary_data.boundary_correction_data);
  const auto local_v_minus = get<::Tags::MirrorView<package_field_tag<3>>>(
      local_boundary_data.boundary_correction_data);
  const auto local_normal_times_v_plus =
      get<::Tags::MirrorView<package_field_tag<4>>>(
          local_boundary_data.boundary_correction_data);
  const auto local_normal_times_v_minus =
      get<::Tags::MirrorView<package_field_tag<5>>>(
          local_boundary_data.boundary_correction_data);
  const auto local_gamma2_v_psi = get<::Tags::MirrorView<package_field_tag<6>>>(
      local_boundary_data.boundary_correction_data);
  const auto local_char_speeds = get<::Tags::MirrorView<package_field_tag<7>>>(
      local_boundary_data.boundary_correction_data);

  const auto remote_v_psi = get<::Tags::MirrorView<package_field_tag<0>>>(
      remote_boundary_data.boundary_correction_data);
  const auto remote_v_zero = get<::Tags::MirrorView<package_field_tag<1>>>(
      remote_boundary_data.boundary_correction_data);
  const auto remote_v_plus = get<::Tags::MirrorView<package_field_tag<2>>>(
      remote_boundary_data.boundary_correction_data);
  const auto remote_v_minus = get<::Tags::MirrorView<package_field_tag<3>>>(
      remote_boundary_data.boundary_correction_data);
  const auto remote_normal_times_v_plus =
      get<::Tags::MirrorView<package_field_tag<4>>>(
          remote_boundary_data.boundary_correction_data);
  const auto remote_normal_times_v_minus =
      get<::Tags::MirrorView<package_field_tag<5>>>(
          remote_boundary_data.boundary_correction_data);
  const auto remote_gamma2_v_psi =
      get<::Tags::MirrorView<package_field_tag<6>>>(
          remote_boundary_data.boundary_correction_data);
  const auto remote_char_speeds = get<::Tags::MirrorView<package_field_tag<7>>>(
      remote_boundary_data.boundary_correction_data);

  Variables<device_dt_boundary_tags> dt_boundary_correction_on_mortar{
      mortar_mesh.number_of_grid_points()};
  const auto dt_psi_on_mortar = get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(
      dt_boundary_correction_on_mortar);
  const auto dt_pi_on_mortar = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(
      dt_boundary_correction_on_mortar);
  const auto dt_phi_on_mortar =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(
          dt_boundary_correction_on_mortar);

  ::Kokkos::parallel_for(
      "ABCTDComputeBoundaryCorrectionOnMortar",
      mortar_mesh.number_of_grid_points(),
      KOKKOS_LAMBDA(const int mortar_index_int) {
        const size_t mortar_index = static_cast<size_t>(mortar_index_int);

        tnsr::i<double, volume_dim, Frame::Inertial> local_v_zero_at_mortar{};
        tnsr::i<double, volume_dim, Frame::Inertial>
            local_normal_times_v_plus_at_mortar{};
        tnsr::i<double, volume_dim, Frame::Inertial>
            local_normal_times_v_minus_at_mortar{};
        tnsr::i<double, volume_dim, Frame::Inertial>
            local_char_speeds_at_mortar{};
        tnsr::i<double, volume_dim, Frame::Inertial> remote_v_zero_at_mortar{};
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
          remote_v_zero_at_mortar.get(d) = remote_v_zero.get(d)[mortar_index];
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
            make_not_null(&dt_psi_at_mortar), make_not_null(&dt_pi_at_mortar),
            make_not_null(&dt_phi_at_mortar), local_v_psi_at_mortar,
            local_v_zero_at_mortar, local_v_plus_at_mortar,
            local_v_minus_at_mortar, local_normal_times_v_plus_at_mortar,
            local_normal_times_v_minus_at_mortar, local_gamma2_v_psi_at_mortar,
            local_char_speeds_at_mortar, remote_v_psi_at_mortar,
            remote_v_zero_at_mortar, remote_v_plus_at_mortar,
            remote_v_minus_at_mortar, remote_normal_times_v_plus_at_mortar,
            remote_normal_times_v_minus_at_mortar,
            remote_gamma2_v_psi_at_mortar, remote_char_speeds_at_mortar);

        get(dt_psi_on_mortar)[mortar_index] = get(dt_psi_at_mortar);
        get(dt_pi_on_mortar)[mortar_index] = get(dt_pi_at_mortar);
        for (size_t d = 0; d < volume_dim; ++d) {
          dt_phi_on_mortar.get(d)[mortar_index] = dt_phi_at_mortar.get(d);
        }
      });

  Variables<device_dt_boundary_tags> dt_boundary_correction_on_face{
      face_mesh.number_of_grid_points()};
  if (needs_projection) {
    ScalarWave::Kokkos::project_from_mortar_device(
        make_not_null(&dt_boundary_correction_on_face),
        dt_boundary_correction_on_mortar, face_mesh, mortar_mesh, mortar_size);
  } else {
    ASSERT(dt_boundary_correction_on_face.number_of_grid_points() ==
               dt_boundary_correction_on_mortar.number_of_grid_points(),
           "Expected identical face and mortar point counts when projection "
           "is not needed.");
    ::Kokkos::deep_copy(dt_boundary_correction_on_face.view(),
                        dt_boundary_correction_on_mortar.view());
  }

  const auto dt_psi_on_face = get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(
      dt_boundary_correction_on_face);
  const auto dt_pi_on_face = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(
      dt_boundary_correction_on_face);
  const auto dt_phi_on_face =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(
          dt_boundary_correction_on_face);
  const auto dt_psi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(
          *dt_boundary_correction_on_face_sum);
  const auto dt_pi_boundary_sum = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(
      *dt_boundary_correction_on_face_sum);
  const auto dt_phi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(
          *dt_boundary_correction_on_face_sum);

  ::Kokkos::parallel_for(
      "ABCTDAccumulateBoundaryCorrectionOnFace",
      face_mesh.number_of_grid_points(),
      KOKKOS_LAMBDA(const int face_index_int) {
        const size_t face_index = static_cast<size_t>(face_index_int);
        get(dt_psi_boundary_sum)[face_index] += get(dt_psi_on_face)[face_index];
        get(dt_pi_boundary_sum)[face_index] += get(dt_pi_on_face)[face_index];
        for (size_t d = 0; d < volume_dim; ++d) {
          dt_phi_boundary_sum.get(d)[face_index] +=
              dt_phi_on_face.get(d)[face_index];
        }
      });
}

void lift_face_correction_to_volume(
    const gsl::not_null<device_dt_type*> device_dt,
    const Variables<device_dt_boundary_tags>&
        dt_boundary_correction_on_face_sum,
    const Direction<volume_dim>& direction,
    const device_face_to_volume_index_map_type& device_face_to_volume_index_map,
    const device_face_normal_magnitude_type& device_face_normal_magnitude,
    const Mesh<volume_dim>& mesh) {
  const auto dt_psi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(*device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(*device_dt);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(*device_dt);
  const auto dt_psi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(
          dt_boundary_correction_on_face_sum);
  const auto dt_pi_boundary_sum = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(
      dt_boundary_correction_on_face_sum);
  const auto dt_phi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(
          dt_boundary_correction_on_face_sum);
  const size_t num_face_points =
      dt_boundary_correction_on_face_sum.number_of_grid_points();
  const size_t extent_perpendicular_to_boundary =
      mesh.extents(direction.dimension());
  const double lift_prefactor =
      -0.5 * static_cast<double>(extent_perpendicular_to_boundary *
                                 (extent_perpendicular_to_boundary - 1));
  const size_t sliced_dim = direction.dimension();

  auto face_to_volume_index =
      direction.side() == Side::Upper
          ? gsl::at(device_face_to_volume_index_map, sliced_dim).second
          : gsl::at(device_face_to_volume_index_map, sliced_dim).first;
  auto face_normal_magnitude =
      direction.side() == Side::Upper
          ? gsl::at(device_face_normal_magnitude, sliced_dim).second
          : gsl::at(device_face_normal_magnitude, sliced_dim).first;

  ::Kokkos::parallel_for(
      "ABCTDApplyBoundaryCorrectionsToTimeDerivativeKokkosFace",
      num_face_points, KOKKOS_LAMBDA(const int face_index_int) {
        const size_t face_index = static_cast<size_t>(face_index_int);
        const size_t volume_index = face_to_volume_index(face_index);

        const double normal_magnitude = face_normal_magnitude(face_index);
        const double lifted_factor = lift_prefactor * normal_magnitude;

        get(dt_psi)[volume_index] +=
            lifted_factor * get(dt_psi_boundary_sum)[face_index];
        get(dt_pi)[volume_index] +=
            lifted_factor * get(dt_pi_boundary_sum)[face_index];
        for (size_t d = 0; d < volume_dim; ++d) {
          dt_phi.get(d)[volume_index] +=
              lifted_factor * dt_phi_boundary_sum.get(d)[face_index];
        }
      });
}

}  // namespace

void ApplyBoundaryCorrectionsToTimeDerivativeKokkos::
    apply_boundary_corrections_on_device(
        const gsl::not_null<device_dt_type*> device_dt,
        const outgoing_boundary_data_type& outgoing_boundary_data,
        const incoming_boundary_data_type& incoming_boundary_data,
        const external_boundary_data_type& external_boundary_data,
        const device_face_to_volume_index_map_type&
            device_face_to_volume_index_map,
        const device_face_normal_magnitude_type& device_face_normal_magnitude,
        const device_mortar_data_type& device_mortar_data,
        const typename mortar_mesh_tag::type& mortar_meshes,
        const Mesh<3>& mesh, const Element<3>& element) {
  ASSERT(mesh.quadrature(0) == Spectral::Quadrature::GaussLobatto,
         "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
         "Gauss-Lobatto quadrature only.");

  for (const auto& [direction, neighbors] : element.neighbors()) {
    const Mesh<2> face_mesh = mesh.slice_away(direction.dimension());
    const size_t num_face_points = face_mesh.number_of_grid_points();
    Variables<device_dt_boundary_tags> dt_boundary_correction_on_face_sum{
        num_face_points};
    ::Kokkos::deep_copy(dt_boundary_correction_on_face_sum.view(), 0.0);

    for (const auto& neighbor : neighbors) {
      const DirectionalId<3> mortar_id{direction, neighbor};
      ASSERT(
          outgoing_boundary_data.count(mortar_id) == 1,
          "Missing local Kokkos boundary data for mortar " << mortar_id << ".");
      ASSERT(incoming_boundary_data.count(mortar_id) == 1,
             "Missing received Kokkos boundary data for mortar " << mortar_id
                                                                 << ".");

      const auto& local_boundary_data = outgoing_boundary_data.at(mortar_id);
      const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);
      const Mesh<2>& mortar_mesh = mortar_meshes.at(mortar_id);
      ASSERT(local_boundary_data.boundary_correction_data
                     .number_of_grid_points() ==
                 mortar_mesh.number_of_grid_points(),
             "Local packaged mortar data size does not match mortar mesh for "
                 << mortar_id << ".");
      ASSERT(remote_boundary_data.boundary_correction_data
                     .number_of_grid_points() ==
                 mortar_mesh.number_of_grid_points(),
             "Remote packaged mortar data size does not match mortar mesh for "
                 << mortar_id << ".");
      const auto& mortar_data = device_mortar_data.at(mortar_id);
      accumulate_mortar_pair_to_face_sum(
          make_not_null(&dt_boundary_correction_on_face_sum),
          local_boundary_data, remote_boundary_data, face_mesh, mortar_mesh,
          mortar_data.needs_projection, mortar_data.mortar_size);
    }
    lift_face_correction_to_volume(
        device_dt, dt_boundary_correction_on_face_sum, direction,
        device_face_to_volume_index_map, device_face_normal_magnitude, mesh);
  }

  const std::array<Spectral::SegmentSize, volume_dim - 1> full_mortar_size{
      Spectral::SegmentSize::Full, Spectral::SegmentSize::Full};
  for (const Direction<3>& direction : element.external_boundaries()) {
    const Mesh<2> face_mesh = mesh.slice_away(direction.dimension());
    Variables<device_dt_boundary_tags> dt_boundary_correction_on_face_sum{
        face_mesh.number_of_grid_points()};
    ::Kokkos::deep_copy(dt_boundary_correction_on_face_sum.view(), 0.0);

    const DirectionalId<3> external_mortar_id{
        direction, ElementId<3>::external_boundary_id()};
    ASSERT(outgoing_boundary_data.count(external_mortar_id) == 1,
           "Missing local packaged data for external boundary "
               << external_mortar_id << ".");
    ASSERT(external_boundary_data.count(external_mortar_id) == 1,
           "Missing external ghost packaged data for boundary "
               << external_mortar_id << ".");

    accumulate_mortar_pair_to_face_sum(
        make_not_null(&dt_boundary_correction_on_face_sum),
        outgoing_boundary_data.at(external_mortar_id),
        external_boundary_data.at(external_mortar_id), face_mesh, face_mesh,
        false, full_mortar_size);
    lift_face_correction_to_volume(
        device_dt, dt_boundary_correction_on_face_sum, direction,
        device_face_to_volume_index_map, device_face_normal_magnitude, mesh);
  }
}

}  // namespace ScalarWave::Actions
