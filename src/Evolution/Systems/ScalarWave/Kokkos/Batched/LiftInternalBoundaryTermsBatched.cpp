// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/LiftInternalBoundaryTermsBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Domain/Structure/Direction.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Actions {
namespace {

static constexpr size_t volume_dim = 3;
static constexpr size_t number_of_faces = 2 * volume_dim;
using device_dt_boundary_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::device_dt_boundary_tags;

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
}

constexpr size_t face_index(const size_t sliced_dim, const size_t side_i) {
  return 2 * sliced_dim + side_i;
}

}  // namespace

void LiftInternalBoundaryTermsBatched::apply(
    const gsl::not_null<evolution::Kokkos::Tags::PackedEvolutionState<
        ScalarWave::System<3>>::type*>
        packed_evolution_state,
    const evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>::type&
        packed_topology,
    const evolution::Kokkos::Tags::PackedGeometry<ScalarWave::System<3>>::type&
        packed_geometry,
    const evolution::Kokkos::Tags::PackedBoundaryScratch<
        ScalarWave::System<3>>::type& packed_boundary_scratch) {
  if (packed_topology.total_points == 0) {
    return;
  }

  ASSERT(packed_topology.uniform_quadrature_host[0] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[1] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[2] ==
                 Spectral::Quadrature::GaussLobatto,
         "LiftInternalBoundaryTermsBatched currently supports "
         "Gauss-Lobatto quadrature only.");

  const auto& extents = packed_topology.uniform_extents_host;
  ASSERT(extents[0] == extents[1] and extents[1] == extents[2],
         "LiftInternalBoundaryTermsBatched assumes uniform p across "
         "dimensions (equal extents in all dimensions).");
  const Mesh<volume_dim> mesh{extents, packed_topology.uniform_basis_host,
                              packed_topology.uniform_quadrature_host};
  const Mesh<volume_dim - 1> uniform_face_mesh = mesh.slice_away(0);
  const size_t num_face_points = uniform_face_mesh.number_of_grid_points();

  const size_t num_elements = packed_topology.local_elements.size();
  const auto& dt_face_sum_for_all_elements =
      packed_boundary_scratch.internal_boundary_terms_for_all_elements;
  for (size_t face_id = 0; face_id < number_of_faces; ++face_id) {
    ASSERT(dt_face_sum_for_all_elements[face_id].number_of_grid_points() ==
               num_elements * num_face_points,
           "ComputeInternalBoundaryTermsBatched must run before lifting "
           "internal batched boundary terms.");
  }

  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& element_point_offsets_device =
      packed_topology.element_point_offsets_device;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  auto& dt_vars = packed_evolution_state->device_dt_variables;
  const auto dt_psi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(dt_vars);
  const auto dt_pi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(dt_vars);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
          dt_vars);

  for (size_t sliced_dim = 0; sliced_dim < volume_dim; ++sliced_dim) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const Side side = side_i == 0 ? Side::Lower : Side::Upper;
      const size_t local_face_id = face_index(sliced_dim, side_i);
      const auto& dt_face_sum_all = dt_face_sum_for_all_elements[local_face_id];
      const auto dt_psi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
              dt_face_sum_all);
      const auto dt_pi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
              dt_face_sum_all);
      const auto dt_phi_face_sum_all = get<
          ::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
          dt_face_sum_all);
      const auto local_face_to_volume_index =
          side == Side::Upper
              ? gsl::at(face_to_volume_index_map, sliced_dim).second
              : gsl::at(face_to_volume_index_map, sliced_dim).first;
      const double local_outward_sign = outward_sign(side);
      const double lift_prefactor =
          -0.5 *
          static_cast<double>(extents[sliced_dim] * (extents[sliced_dim] - 1));

      ::Kokkos::parallel_for(
          "LiftInternalBoundaryTermsBatchedLiftFace",
          num_elements * num_face_points,
          KOKKOS_LAMBDA(const int linear_index_int) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t face_index_on_face = linear_index % num_face_points;
            const size_t element_index = linear_index / num_face_points;
            const size_t local_face_linear_index =
                element_index * num_face_points + face_index_on_face;
            const size_t local_volume_index_in_element =
                local_face_to_volume_index(face_index_on_face);
            const size_t local_volume_index =
                element_point_offsets_device(element_index) +
                local_volume_index_in_element;
            double local_normal_magnitude_squared = 0.0;
            for (size_t d = 0; d < volume_dim; ++d) {
              const double local_component =
                  local_outward_sign *
                  inverse_jacobian(element_index, local_volume_index_in_element,
                                   sliced_dim * volume_dim + d);
              local_normal_magnitude_squared +=
                  local_component * local_component;
            }
            const double local_normal_magnitude =
                sqrt(local_normal_magnitude_squared);
            const double lifted_factor =
                lift_prefactor * local_normal_magnitude;
            get(dt_psi)[local_volume_index] +=
                lifted_factor *
                get(dt_psi_face_sum_all)[local_face_linear_index];
            get(dt_pi)[local_volume_index] +=
                lifted_factor *
                get(dt_pi_face_sum_all)[local_face_linear_index];
            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi.get(d)[local_volume_index] +=
                  lifted_factor *
                  dt_phi_face_sum_all.get(d)[local_face_linear_index];
            }
          });
    }
  }
}

}  // namespace ScalarWave::Actions
