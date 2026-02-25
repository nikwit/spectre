// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/ComputeTimeDerivativeBatched.hpp"

namespace ScalarWave::Actions {

void ComputeTimeDerivativeBatched::compute_time_derivative_batched_volume_impl(
    const gsl::not_null<ScalarWave::Batched::Tags::PackedEvolutionState::type*>
        packed_evolution_state,
    const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
    const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry) {
  const size_t total_points = packed_topology.total_points;
  if (total_points == 0) {
    return;
  }

  const auto extents = packed_topology.uniform_extents_host;
  const size_t points_per_element = packed_topology.points_per_element;
  ASSERT(extents[0] * extents[1] * extents[2] == points_per_element,
         "Uniform extents and points-per-element mismatch.");
  ASSERT(packed_topology.local_elements.size() * points_per_element ==
             total_points,
         "Packed point count mismatch in batched volume derivative.");
  const Mesh<3> mesh{extents, packed_topology.uniform_basis_host,
                     packed_topology.uniform_quadrature_host};

  auto& device_vars = packed_evolution_state->device_variables;
  auto& device_dt = packed_evolution_state->device_dt_variables;
  const auto& gamma2_full = packed_evolution_state->device_constraint_gamma2;
  Variables<device_derivative_tags> partial_derivatives_all_elements{
      total_points};

  partial_derivatives_batched(make_not_null(&partial_derivatives_all_elements),
                              device_vars, mesh,
                              packed_geometry.element_inverse_jacobian_device);

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

void ComputeTimeDerivativeBatched::apply(
    const gsl::not_null<ScalarWave::Batched::Tags::PackedEvolutionState::type*>
        packed_evolution_state,
    const ScalarWave::Batched::Tags::PackedTopology::type& packed_topology,
    const ScalarWave::Batched::Tags::PackedGeometry::type& packed_geometry) {
  compute_time_derivative_batched_volume_impl(packed_evolution_state,
                                              packed_topology, packed_geometry);
}

}  // namespace ScalarWave::Actions
