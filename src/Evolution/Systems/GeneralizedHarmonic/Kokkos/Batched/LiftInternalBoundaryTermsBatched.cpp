// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/LiftInternalBoundaryTermsBatched.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

static constexpr size_t volume_dim = 3;
static constexpr size_t number_of_faces = 2 * volume_dim;
using system = gh::System<volume_dim>;
using device_dt_boundary_tags =
    evolution::Kokkos::PackedBoundaryScratch<system>::device_dt_boundary_tags;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>;

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
}

constexpr size_t face_index(const size_t sliced_dim, const size_t side_i) {
  return 2 * sliced_dim + side_i;
}

KOKKOS_INLINE_FUNCTION void inverse_spatial_metric_and_det(
    const gsl::not_null<tnsr::II<double, volume_dim, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<double*> det_spatial_metric,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& spacetime_metric) {
  const double g00 = spacetime_metric.get(1, 1);
  const double g01 = spacetime_metric.get(1, 2);
  const double g02 = spacetime_metric.get(1, 3);
  const double g11 = spacetime_metric.get(2, 2);
  const double g12 = spacetime_metric.get(2, 3);
  const double g22 = spacetime_metric.get(3, 3);

  *det_spatial_metric = g00 * (g11 * g22 - g12 * g12) -
                        g01 * (g01 * g22 - g12 * g02) +
                        g02 * (g01 * g12 - g11 * g02);
  const double inv_det = 1.0 / *det_spatial_metric;

  inverse_spatial_metric->get(0, 0) = (g11 * g22 - g12 * g12) * inv_det;
  inverse_spatial_metric->get(0, 1) = (g02 * g12 - g01 * g22) * inv_det;
  inverse_spatial_metric->get(0, 2) = (g01 * g12 - g02 * g11) * inv_det;
  inverse_spatial_metric->get(1, 1) = (g00 * g22 - g02 * g02) * inv_det;
  inverse_spatial_metric->get(1, 2) = (g02 * g01 - g00 * g12) * inv_det;
  inverse_spatial_metric->get(2, 2) = (g00 * g11 - g01 * g01) * inv_det;
}

}  // namespace

void LiftInternalBoundaryTermsBatched::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const typename packed_boundary_scratch_tag::type& packed_boundary_scratch) {
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
  const auto& vars = packed_evolution_state->device_variables;
  const auto spacetime_metric =
      get<::Tags::MirrorView<spacetime_metric_tag>>(vars);
  auto& dt_vars = packed_evolution_state->device_dt_variables;
  const auto dt_spacetime_metric =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(dt_vars);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<pi_tag>>>(dt_vars);
  const auto dt_phi = get<::Tags::MirrorView<::Tags::dt<phi_tag>>>(dt_vars);

  for (size_t sliced_dim = 0; sliced_dim < volume_dim; ++sliced_dim) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const Side side = side_i == 0 ? Side::Lower : Side::Upper;
      const size_t local_face_id = face_index(sliced_dim, side_i);
      const auto& dt_face_sum_all = dt_face_sum_for_all_elements[local_face_id];
      const auto dt_spacetime_metric_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
              dt_face_sum_all);
      const auto dt_pi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<pi_tag>>>(dt_face_sum_all);
      const auto dt_phi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<phi_tag>>>(dt_face_sum_all);
      const auto local_face_to_volume_index =
          side == Side::Upper
              ? gsl::at(face_to_volume_index_map, sliced_dim).second
              : gsl::at(face_to_volume_index_map, sliced_dim).first;
      const double local_outward_sign = outward_sign(side);
      const double lift_prefactor =
          -0.5 *
          static_cast<double>(extents[sliced_dim] * (extents[sliced_dim] - 1));

      ::Kokkos::parallel_for(
          "GhLiftInternalBoundaryTermsBatchedLiftFace",
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

            tnsr::aa<double, volume_dim, Frame::Inertial>
                local_spacetime_metric{};
            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                local_spacetime_metric.get(a, b) =
                    spacetime_metric.get(a, b)[local_volume_index];
              }
            }
            tnsr::II<double, volume_dim, Frame::Inertial>
                inverse_spatial_metric{};
            double det_spatial_metric = 0.0;
            inverse_spatial_metric_and_det(
                make_not_null(&inverse_spatial_metric),
                make_not_null(&det_spatial_metric), local_spacetime_metric);
            (void)det_spatial_metric;

            tnsr::i<double, volume_dim, Frame::Inertial>
                unnormalized_normal_covector{};
            for (size_t d = 0; d < volume_dim; ++d) {
              unnormalized_normal_covector.get(d) =
                  local_outward_sign *
                  inverse_jacobian(element_index, local_volume_index_in_element,
                                   sliced_dim * volume_dim + d);
            }
            tnsr::I<double, volume_dim, Frame::Inertial>
                unnormalized_normal_vector{};
            double normal_magnitude_squared = 0.0;
            for (size_t i = 0; i < volume_dim; ++i) {
              unnormalized_normal_vector.get(i) = 0.0;
              for (size_t j = 0; j < volume_dim; ++j) {
                unnormalized_normal_vector.get(i) +=
                    inverse_spatial_metric.get(i, j) *
                    unnormalized_normal_covector.get(j);
              }
              normal_magnitude_squared += unnormalized_normal_vector.get(i) *
                                          unnormalized_normal_covector.get(i);
            }
            const double lifted_factor =
                lift_prefactor * sqrt(normal_magnitude_squared);

            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                dt_spacetime_metric.get(a, b)[local_volume_index] +=
                    lifted_factor * dt_spacetime_metric_face_sum_all.get(
                                        a, b)[local_face_linear_index];
                dt_pi.get(a, b)[local_volume_index] +=
                    lifted_factor *
                    dt_pi_face_sum_all.get(a, b)[local_face_linear_index];
                for (size_t d = 0; d < volume_dim; ++d) {
                  dt_phi.get(d, a, b)[local_volume_index] +=
                      lifted_factor *
                      dt_phi_face_sum_all.get(d, a, b)[local_face_linear_index];
                }
              }
            }
          });
    }
  }
}

}  // namespace gh::Actions
