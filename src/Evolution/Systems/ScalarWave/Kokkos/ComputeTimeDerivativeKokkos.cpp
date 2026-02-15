// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/ComputeTimeDerivativeKokkos.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime_api.h>
#endif

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Actions {

void ComputeTimeDerivativeKokkos::
    compute_volume_terms_and_package_boundary_data(
        const gsl::not_null<device_dt_type*> device_dt,
        const gsl::not_null<outgoing_boundary_data_type*>
            outgoing_boundary_data,
        const device_variables_type& device_vars,
        const device_inverse_jacobian_type& device_inverse_jacobian,
        const device_constraint_gamma2_type& device_constraint_gamma2,
        const device_face_to_volume_index_map_type&
            device_face_to_volume_index_map,
        const device_face_unit_normal_covector_type&
            device_face_unit_normal_covector,
        const Mesh<volume_dim>& mesh, const Element<volume_dim>& element,
        const TimeStepId& time_step_id) {
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence("ComputeTimeDerivativeKokkosActionEntry");
  const auto action_entry_err = cudaPeekAtLastError();
  ASSERT(action_entry_err == cudaSuccess,
         "CUDA error on entry to ComputeTimeDerivativeKokkos action: "
             << cudaGetErrorString(action_entry_err));
#endif

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
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence("ComputeTimeDerivativeKokkosAfterDtInit");
  const auto after_dt_init_err = cudaPeekAtLastError();
  ASSERT(after_dt_init_err == cudaSuccess,
         "CUDA error after DeviceDtVariables initialization: "
             << cudaGetErrorString(after_dt_init_err));
#endif

#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence("ComputeTimeDerivativeKokkosBeforePartialDerivatives");
  const auto before_partial_derivatives_err = cudaPeekAtLastError();
  ASSERT(before_partial_derivatives_err == cudaSuccess,
         "CUDA error before partial_derivatives: "
             << cudaGetErrorString(before_partial_derivatives_err));

  Kokkos::parallel_for("ComputeTimeDerivativeKokkosProbeBeforePartials", 1,
                       KOKKOS_LAMBDA(const int /*unused*/){});
  Kokkos::fence("ComputeTimeDerivativeKokkosAfterProbeBeforePartials");
  const auto probe_before_partials_err = cudaPeekAtLastError();
  ASSERT(probe_before_partials_err == cudaSuccess,
         "Probe kernel failed right before partial_derivatives: "
             << cudaGetErrorString(probe_before_partials_err));
#endif

  try {
    partial_derivatives(make_not_null(&device_partial_derivatives), device_vars,
                        mesh, device_inverse_jacobian);
  } catch (const std::exception& e) {
    ERROR_NO_TRACE(
        "partial_derivatives failed in "
        "ComputeTimeDerivativeKokkos: "
        << e.what());
  }
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence("ComputeTimeDerivativeKokkosAfterPartialDerivatives");
  const auto after_partial_derivatives_err = cudaPeekAtLastError();
  ASSERT(after_partial_derivatives_err == cudaSuccess,
         "CUDA error after partial_derivatives: "
             << cudaGetErrorString(after_partial_derivatives_err));
#endif

  const auto dt_psi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Psi>>>(*device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<Tags::Pi>>>(*device_dt);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<Tags::Phi<volume_dim>>>>(*device_dt);
  static constexpr size_t Dim = 3;

#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::fence("ComputeTimeDerivativeKokkosPreVolumeTerms");
  const auto pre_volume_terms_err = cudaPeekAtLastError();
  ASSERT(pre_volume_terms_err == cudaSuccess,
         "Pre-kernel CUDA error before ComputeTimeDerivativeKokkosVolumeTerms: "
             << cudaGetErrorString(pre_volume_terms_err));

  Kokkos::parallel_for("ComputeTimeDerivativeKokkosProbeKernel", 1,
                       KOKKOS_LAMBDA(const int /*unused*/){});
  Kokkos::fence("ComputeTimeDerivativeKokkosAfterProbeKernel");
  const auto probe_kernel_err = cudaPeekAtLastError();
  ASSERT(probe_kernel_err == cudaSuccess,
         "Probe kernel failed before ComputeTimeDerivativeKokkosVolumeTerms: "
             << cudaGetErrorString(probe_kernel_err));
#endif

  Kokkos::parallel_for(
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

  outgoing_boundary_data->clear();

  for (const auto& [direction, neighbors] : element.neighbors()) {
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

    Variables<device_package_field_tags> packaged_face_data{num_face_points};

    if (num_face_points > 0) {
      Kokkos::parallel_for(
          "PackageBoundaryCorrectionData", num_face_points,
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

    for (const auto& neighbor : neighbors) {
      ScalarWave::KokkosTags::BoundaryCorrectionData<volume_dim>
          boundary_data{};
      boundary_data.volume_mesh = mesh;
      boundary_data.boundary_correction_mesh = face_mesh;
      boundary_data.boundary_correction_data = packaged_face_data;
      boundary_data.validity_range = time_step_id;
      boundary_data.integration_order = 0;
      outgoing_boundary_data->insert_or_assign(
          DirectionalId<volume_dim>{direction, neighbor},
          std::move(boundary_data));
    }
  }
}

}  // namespace ScalarWave::Actions
