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
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BoundaryCorrectionHelpers.hpp"
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

            detail::compute_packaged_boundary_data_at_point<volume_dim>(
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
