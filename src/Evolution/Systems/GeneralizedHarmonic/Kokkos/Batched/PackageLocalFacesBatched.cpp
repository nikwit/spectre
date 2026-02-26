// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/PackageLocalFacesBatched.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Variables.hpp"
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
using system = gh::System<volume_dim>;
using device_package_field_tags =
    evolution::Kokkos::PackedBoundaryScratch<system>::device_package_field_tags;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>;

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
}

constexpr size_t local_face_index(const size_t sliced_dim,
                                  const size_t side_i) {
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

KOKKOS_INLINE_FUNCTION void compute_packaged_boundary_data_at_point(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_spacetime_metric,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_zero,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_plus,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_minus,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_n_times_v_minus,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_gamma2_v_spacetime_metric,
    const gsl::not_null<tnsr::a<double, volume_dim, Frame::Inertial>*>
        char_speeds,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>& phi,
    const double gamma1, const double gamma2,
    const tnsr::i<double, volume_dim, Frame::Inertial>&
        unnormalized_normal_covector) {
  tnsr::II<double, volume_dim, Frame::Inertial> inverse_spatial_metric{};
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det(make_not_null(&inverse_spatial_metric),
                                 make_not_null(&det_spatial_metric),
                                 spacetime_metric);
  (void)det_spatial_metric;

  tnsr::I<double, volume_dim, Frame::Inertial> shift{};
  for (size_t i = 0; i < volume_dim; ++i) {
    shift.get(i) = 0.0;
    for (size_t j = 0; j < volume_dim; ++j) {
      shift.get(i) +=
          inverse_spatial_metric.get(i, j) * spacetime_metric.get(0, j + 1);
    }
  }
  double lapse_squared = -spacetime_metric.get(0, 0);
  for (size_t i = 0; i < volume_dim; ++i) {
    lapse_squared += shift.get(i) * spacetime_metric.get(0, i + 1);
  }
  const double lapse = sqrt(lapse_squared);

  tnsr::I<double, volume_dim, Frame::Inertial> normal_vector{};
  double normal_magnitude_squared = 0.0;
  for (size_t i = 0; i < volume_dim; ++i) {
    normal_vector.get(i) = 0.0;
    for (size_t j = 0; j < volume_dim; ++j) {
      normal_vector.get(i) += inverse_spatial_metric.get(i, j) *
                              unnormalized_normal_covector.get(j);
    }
    normal_magnitude_squared +=
        normal_vector.get(i) * unnormalized_normal_covector.get(i);
  }
  const double one_over_normal_magnitude = 1.0 / sqrt(normal_magnitude_squared);
  tnsr::i<double, volume_dim, Frame::Inertial> normal_covector{};
  for (size_t i = 0; i < volume_dim; ++i) {
    normal_vector.get(i) *= one_over_normal_magnitude;
    normal_covector.get(i) =
        unnormalized_normal_covector.get(i) * one_over_normal_magnitude;
  }

  double shift_dot_normal = 0.0;
  for (size_t i = 0; i < volume_dim; ++i) {
    shift_dot_normal += shift.get(i) * normal_covector.get(i);
  }
  shift_dot_normal *= -1.0;

  char_speeds->get(0) = (1.0 + gamma1) * shift_dot_normal;
  char_speeds->get(1) = shift_dot_normal;
  char_speeds->get(2) = lapse + shift_dot_normal;
  char_speeds->get(3) = -lapse + shift_dot_normal;

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      char_speed_gamma2_v_spacetime_metric->get(a, b) =
          gamma2 * spacetime_metric.get(a, b);
    }
  }

  tnsr::aa<double, volume_dim, Frame::Inertial> normal_dot_phi{};
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      normal_dot_phi.get(a, b) = normal_vector.get(0) * phi.get(0, a, b);
      for (size_t i = 1; i < volume_dim; ++i) {
        normal_dot_phi.get(a, b) += normal_vector.get(i) * phi.get(i, a, b);
      }
    }
  }

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      char_speed_v_plus->get(a, b) =
          char_speeds->get(2) *
          (pi.get(a, b) + normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));
      char_speed_v_minus->get(a, b) =
          char_speeds->get(3) *
          (pi.get(a, b) - normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));

      for (size_t i = 0; i < volume_dim; ++i) {
        char_speed_v_zero->get(i, a, b) =
            char_speeds->get(1) *
            (phi.get(i, a, b) -
             normal_covector.get(i) * normal_dot_phi.get(a, b));
      }
    }
  }

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      for (size_t i = 0; i < volume_dim; ++i) {
        char_speed_n_times_v_plus->get(i, a, b) =
            char_speed_v_plus->get(a, b) * normal_covector.get(i);
        char_speed_n_times_v_minus->get(i, a, b) =
            char_speed_v_minus->get(a, b) * normal_covector.get(i);
      }
      char_speed_v_spacetime_metric->get(a, b) =
          char_speeds->get(0) * spacetime_metric.get(a, b);
      char_speed_gamma2_v_spacetime_metric->get(a, b) *= char_speeds->get(0);
    }
  }
}

}  // namespace

void PackageLocalFacesBatched::apply(
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const typename packed_evolution_state_tag::type& packed_evolution_state,
    const typename device_constraint_gamma0_tag::type& /*device_gamma0*/,
    const typename device_constraint_gamma1_tag::type& device_gamma1,
    const typename device_constraint_gamma2_tag::type& device_gamma2) {
  if (packed_topology.total_points == 0) {
    return;
  }

  ASSERT(packed_topology.uniform_quadrature_host[0] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[1] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[2] ==
                 Spectral::Quadrature::GaussLobatto,
         "PackageLocalFacesBatched currently supports Gauss-Lobatto "
         "quadrature only.");

  const auto& extents = packed_topology.uniform_extents_host;
  ASSERT(extents[0] == extents[1] and extents[1] == extents[2],
         "PackageLocalFacesBatched assumes uniform p across dimensions "
         "(equal extents in all dimensions).");
  const Mesh<volume_dim> mesh{extents, packed_topology.uniform_basis_host,
                              packed_topology.uniform_quadrature_host};
  const Mesh<volume_dim - 1> uniform_face_mesh = mesh.slice_away(0);
  const size_t num_face_points = uniform_face_mesh.number_of_grid_points();
  const size_t expected_face_points =
      packed_topology.local_elements.size() * num_face_points;

  auto& packaged_face_data_for_all_elements =
      packed_boundary_scratch->packaged_face_data_for_all_elements;
  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& vars = packed_evolution_state.device_variables;
  const auto& gamma1 = device_gamma1;
  const auto& gamma2 = device_gamma2;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  const auto& element_point_offsets_device =
      packed_topology.element_point_offsets_device;

  for (size_t d = 0; d < volume_dim; ++d) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const Side side = side_i == 0 ? Side::Lower : Side::Upper;
      auto& packaged_face_data =
          packaged_face_data_for_all_elements[local_face_index(d, side_i)];
      if (packaged_face_data.number_of_grid_points() != expected_face_points) {
        packaged_face_data =
            Variables<device_package_field_tags>{expected_face_points};
      }
      const auto packaged_face_data_view = packaged_face_data.view();
      const auto face_to_volume_index =
          side == Side::Upper ? gsl::at(face_to_volume_index_map, d).second
                              : gsl::at(face_to_volume_index_map, d).first;
      const double outward = outward_sign(side);
      ::Kokkos::parallel_for(
          "GhPackageLocalFacesBatched", expected_face_points,
          KOKKOS_LAMBDA(const int linear_index_int) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t face_index = linear_index % num_face_points;
            const size_t element_index = linear_index / num_face_points;
            const size_t local_volume_index_in_element =
                face_to_volume_index(face_index);
            const size_t local_volume_index =
                element_point_offsets_device(element_index) +
                local_volume_index_in_element;

            tnsr::i<double, volume_dim, Frame::Inertial>
                unnormalized_normal_covector{};
            for (size_t i = 0; i < volume_dim; ++i) {
              unnormalized_normal_covector.get(i) =
                  outward * inverse_jacobian(element_index,
                                             local_volume_index_in_element,
                                             d * volume_dim + i);
            }

            const auto vars_at_s = make_at_index(vars, local_volume_index);
            const auto gamma1_at_s = make_at_index(gamma1, local_volume_index);
            const auto gamma2_at_s = make_at_index(gamma2, local_volume_index);

            tnsr::aa<double, volume_dim, Frame::Inertial> v_spacetime_metric{};
            tnsr::iaa<double, volume_dim, Frame::Inertial> v_zero{};
            tnsr::aa<double, volume_dim, Frame::Inertial> v_plus{};
            tnsr::aa<double, volume_dim, Frame::Inertial> v_minus{};
            tnsr::iaa<double, volume_dim, Frame::Inertial>
                normal_times_v_plus{};
            tnsr::iaa<double, volume_dim, Frame::Inertial>
                normal_times_v_minus{};
            tnsr::aa<double, volume_dim, Frame::Inertial>
                gamma2_v_spacetime_metric{};
            tnsr::a<double, volume_dim, Frame::Inertial> char_speeds{};

            compute_packaged_boundary_data_at_point(
                make_not_null(&v_spacetime_metric), make_not_null(&v_zero),
                make_not_null(&v_plus), make_not_null(&v_minus),
                make_not_null(&normal_times_v_plus),
                make_not_null(&normal_times_v_minus),
                make_not_null(&gamma2_v_spacetime_metric),
                make_not_null(&char_speeds),
                get<::Tags::AtIndex<::Tags::MirrorView<spacetime_metric_tag>>>(
                    vars_at_s),
                get<::Tags::AtIndex<::Tags::MirrorView<pi_tag>>>(vars_at_s),
                get<::Tags::AtIndex<::Tags::MirrorView<phi_tag>>>(vars_at_s),
                get(gamma1_at_s), get(gamma2_at_s),
                unnormalized_normal_covector);

            size_t component_offset = 0;
            for (size_t c = 0; c < v_spacetime_metric.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  v_spacetime_metric[c];
            }
            component_offset += v_spacetime_metric.size();
            for (size_t c = 0; c < v_zero.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  v_zero[c];
            }
            component_offset += v_zero.size();
            for (size_t c = 0; c < v_plus.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  v_plus[c];
            }
            component_offset += v_plus.size();
            for (size_t c = 0; c < v_minus.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  v_minus[c];
            }
            component_offset += v_minus.size();
            for (size_t c = 0; c < normal_times_v_plus.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  normal_times_v_plus[c];
            }
            component_offset += normal_times_v_plus.size();
            for (size_t c = 0; c < normal_times_v_minus.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  normal_times_v_minus[c];
            }
            component_offset += normal_times_v_minus.size();
            for (size_t c = 0; c < gamma2_v_spacetime_metric.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  gamma2_v_spacetime_metric[c];
            }
            component_offset += gamma2_v_spacetime_metric.size();
            for (size_t c = 0; c < char_speeds.size(); ++c) {
              packaged_face_data_view(linear_index, component_offset + c) =
                  char_speeds[c];
            }
          });
    }
  }
}

}  // namespace gh::Actions
