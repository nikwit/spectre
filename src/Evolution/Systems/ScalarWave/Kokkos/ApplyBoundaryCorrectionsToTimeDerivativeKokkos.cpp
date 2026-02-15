// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"

#include <cmath>
#include <cstddef>
#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime_api.h>
#endif

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/CudaDiagnostics.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Actions {

void ApplyBoundaryCorrectionsToTimeDerivativeKokkos::
    apply_boundary_corrections_on_device(
        const gsl::not_null<device_dt_type*> device_dt,
        const outgoing_boundary_data_type& outgoing_boundary_data,
        const incoming_boundary_data_type& incoming_boundary_data,
        const device_face_to_volume_index_map_type&
            device_face_to_volume_index_map,
        const device_face_normal_magnitude_type& device_face_normal_magnitude,
        const Mesh<3>& mesh, const Element<3>& element) {
  detail::check_cuda_error_and_clear("ApplyBoundaryCorrectionsKokkosActionEntry");

  ASSERT(mesh.quadrature(0) == Spectral::Quadrature::GaussLobatto,
         "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
         "Gauss-Lobatto quadrature only.");

  const auto dt_psi = get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(*device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(*device_dt);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<3>>>>(*device_dt);

  for (const auto& [direction, neighbors] : element.neighbors()) {
    ASSERT(neighbors.size() == 1,
           "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
           "uniform meshes with one neighbor per direction.");
    const auto& neighbor = *neighbors.begin();
    ASSERT(neighbors.orientation(neighbor).is_aligned(),
           "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
           "aligned neighbors only.");
    const DirectionalId<3> mortar_id{direction, neighbor};

    ASSERT(outgoing_boundary_data.count(mortar_id) == 1,
           "Missing local Kokkos boundary data for mortar " << mortar_id
                                                             << ".");
    ASSERT(incoming_boundary_data.count(mortar_id) == 1,
           "Missing received Kokkos boundary data for mortar " << mortar_id
                                                                << ".");

    const auto& local_boundary_data = outgoing_boundary_data.at(mortar_id);
    const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);

    ASSERT(local_boundary_data.boundary_correction_mesh ==
               remote_boundary_data.boundary_correction_mesh,
           "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently assumes "
           "matching face meshes across each mortar.");

    const auto& local_packaged_data =
        local_boundary_data.boundary_correction_data;
    const auto& remote_packaged_data =
        remote_boundary_data.boundary_correction_data;
    const size_t num_face_points =
        local_boundary_data.boundary_correction_mesh.number_of_grid_points();
    const size_t extent_perpendicular_to_boundary =
        mesh.extents(direction.dimension());
    const double lift_prefactor =
        -0.5 *
        static_cast<double>(extent_perpendicular_to_boundary *
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

#ifdef KOKKOS_ENABLE_CUDA
    Kokkos::fence("ApplyBoundaryCorrectionsKokkosPreFaceKernel");
    const auto pre_face_kernel_err = cudaPeekAtLastError();
    ASSERT(pre_face_kernel_err == cudaSuccess,
           "CUDA error before ApplyBoundaryCorrections face kernel: "
               << cudaGetErrorString(pre_face_kernel_err));

    Kokkos::parallel_for("ApplyBoundaryCorrectionsKokkosProbeKernel", 1,
                         KOKKOS_LAMBDA(const int /*unused*/) {});
    Kokkos::fence("ApplyBoundaryCorrectionsKokkosAfterProbeKernel");
    const auto probe_face_kernel_err = cudaPeekAtLastError();
    ASSERT(probe_face_kernel_err == cudaSuccess,
           "Probe kernel failed before ApplyBoundaryCorrections face kernel: "
               << cudaGetErrorString(probe_face_kernel_err));
#endif

    Kokkos::parallel_for(
        "ApplyBoundaryCorrectionsToTimeDerivativeKokkosFace", num_face_points,
        KOKKOS_LAMBDA(const int face_index_int) {
          const size_t face_index = static_cast<size_t>(face_index_int);
          const size_t volume_index = face_to_volume_index(face_index);

          const auto local_packaged_at_face =
              make_at_index(local_packaged_data, face_index);
          const auto remote_packaged_at_face =
              make_at_index(remote_packaged_data, face_index);

          Scalar<double> psi_boundary_correction{};
          Scalar<double> pi_boundary_correction{};
          tnsr::i<double, 3, Frame::Inertial> phi_boundary_correction{};

          ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
              3, double>(
              make_not_null(&psi_boundary_correction),
              make_not_null(&pi_boundary_correction),
              make_not_null(&phi_boundary_correction),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<0>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<1>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<2>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<3>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<4>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<5>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<6>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<7>>>>(
                  local_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<0>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<1>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<2>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<3>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<4>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<5>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<6>>>>(
                  remote_packaged_at_face),
              get<::Tags::AtIndex<::Tags::MirrorView<package_field_tag<7>>>>(
                  remote_packaged_at_face));

          const double normal_magnitude = face_normal_magnitude(face_index);
          const double lifted_factor = lift_prefactor * normal_magnitude;

          get(dt_psi)[volume_index] += lifted_factor * get(psi_boundary_correction);
          get(dt_pi)[volume_index] += lifted_factor * get(pi_boundary_correction);
          for (size_t d = 0; d < 3; ++d) {
            dt_phi.get(d)[volume_index] +=
                lifted_factor * phi_boundary_correction.get(d);
          }
        });

#ifdef KOKKOS_ENABLE_CUDA
    Kokkos::fence("ApplyBoundaryCorrectionsKokkosPostFaceKernel");
    const auto post_face_kernel_err = cudaPeekAtLastError();
    ASSERT(post_face_kernel_err == cudaSuccess,
           "CUDA error after ApplyBoundaryCorrections face kernel: "
               << cudaGetErrorString(post_face_kernel_err));
#endif
  }
}

}  // namespace ScalarWave::Actions
