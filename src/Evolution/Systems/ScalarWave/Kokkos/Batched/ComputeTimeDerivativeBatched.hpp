// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

struct ComputeTimeDerivativeBatched {
 private:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using volume_time_derivative_terms =
      typename system::compute_volume_time_derivative_terms;
  using device_gradient_tags =
      db::wrap_tags_in<::Tags::MirrorView, typename system::gradient_variables>;
  using device_derivative_tags =
      db::wrap_tags_in<::Tags::deriv, device_gradient_tags,
                       tmpl::size_t<volume_dim>, Frame::Inertial>;

 public:
  static void compute_time_derivative_batched_volume_impl(
      const gsl::not_null<ScalarWave::Batched::Tags::DeviceData::type*>
          device_data) {
    const size_t total_points = device_data->total_points();
    if (total_points == 0) {
      return;
    }

    const auto extents = device_data->uniform_extents_host();
    const size_t points_per_element = device_data->points_per_element();
    ASSERT(extents[0] * extents[1] * extents[2] == points_per_element,
           "Uniform extents and points-per-element mismatch.");
    ASSERT(device_data->num_elements() * points_per_element == total_points,
           "Packed point count mismatch in batched volume derivative.");
    const Mesh<3> mesh{extents, device_data->uniform_basis_host(),
                       device_data->uniform_quadrature_host()};

    auto& device_vars = device_data->device_variables();
    auto& device_dt = device_data->device_dt_variables();
    const auto& gamma2_full = device_data->device_constraint_gamma2();
    Variables<device_derivative_tags> partial_derivatives_all_elements{
        total_points};

    partial_derivatives_batched(
        make_not_null(&partial_derivatives_all_elements), device_vars, mesh,
        device_data->element_inverse_jacobian_device());

    const auto dt_psi =
        get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(device_dt);
    const auto dt_pi =
        get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(device_dt);
    const auto dt_phi =
        get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
            device_dt);
    ::Kokkos::parallel_for(
        "BatchedScalarWaveVolumeTimeDerivative", total_points,
        KOKKOS_LAMBDA(const int s) {
          Scalar<double> dt_psi_at_s{};
          Scalar<double> dt_pi_at_s{};
          tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_at_s{};
          Scalar<double> result_gamma2_at_s{};

          const auto vars_at_s = make_at_index(device_vars, s);
          const auto derivs_at_s =
              make_at_index(partial_derivatives_all_elements, s);
          const auto gamma2_at_s = make_at_index(gamma2_full, s);

          const auto decisions = volume_time_derivative_terms::apply(
              make_not_null(&dt_psi_at_s), make_not_null(&dt_pi_at_s),
              make_not_null(&dt_phi_at_s), make_not_null(&result_gamma2_at_s),
              get<::Tags::AtIndex<
                  ::Tags::deriv<::Tags::MirrorView<ScalarWave::Tags::Psi>,
                                tmpl::size_t<volume_dim>, Frame::Inertial>>>(
                  derivs_at_s),
              get<::Tags::AtIndex<
                  ::Tags::deriv<::Tags::MirrorView<ScalarWave::Tags::Pi>,
                                tmpl::size_t<volume_dim>, Frame::Inertial>>>(
                  derivs_at_s),
              get<::Tags::AtIndex<::Tags::deriv<
                  ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>,
                  tmpl::size_t<volume_dim>, Frame::Inertial>>>(derivs_at_s),
              get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Pi>>>(
                  vars_at_s),
              get<::Tags::AtIndex<
                  ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>>(
                  vars_at_s),
              gamma2_at_s);
          (void)decisions;
          (void)result_gamma2_at_s;

          get(dt_psi)[s] = get(dt_psi_at_s);
          get(dt_pi)[s] = get(dt_pi_at_s);
          for (size_t d = 0; d < volume_dim; ++d) {
            dt_phi.get(d)[s] = dt_phi_at_s.get(d);
          }
        });
  }

  using return_tags = tmpl::list<ScalarWave::Batched::Tags::DeviceData>;
  using argument_tags = tmpl::list<>;

  static void apply(
      const gsl::not_null<ScalarWave::Batched::Tags::DeviceData::type*>
          device_data) {
    compute_time_derivative_batched_volume_impl(device_data);
  }
};

}  // namespace ScalarWave::Actions
