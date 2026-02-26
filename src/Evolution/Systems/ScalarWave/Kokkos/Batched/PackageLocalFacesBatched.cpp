// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/PackageLocalFacesBatched.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {
namespace {

static constexpr size_t volume_dim = 3;
using package_field_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::package_field_tags;
template <size_t I>
using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
using device_package_field_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::device_package_field_tags;

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
}

constexpr size_t local_face_index(const size_t sliced_dim,
                                  const size_t side_i) {
  return 2 * sliced_dim + side_i;
}

}  // namespace

void PackageLocalFacesBatched::apply(
    const gsl::not_null<evolution::Kokkos::Tags::PackedBoundaryScratch<
        ScalarWave::System<3>>::type*>
        packed_boundary_scratch,
    const evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>::type&
        packed_topology,
    const evolution::Kokkos::Tags::PackedGeometry<ScalarWave::System<3>>::type&
        packed_geometry,
    const typename ::Tags::MirrorView<ScalarWave::Tags::ConstraintGamma2>::type&
        device_constraint_gamma2,
    const evolution::Kokkos::Tags::PackedEvolutionState<
        ScalarWave::System<3>>::type& packed_evolution_state) {
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
  const auto& gamma2 = device_constraint_gamma2;
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
      const auto face_to_volume_index =
          side == Side::Upper ? gsl::at(face_to_volume_index_map, d).second
                              : gsl::at(face_to_volume_index_map, d).first;
      const double outward = outward_sign(side);
      ::Kokkos::parallel_for(
          "PackageLocalFacesBatched", expected_face_points,
          KOKKOS_LAMBDA(const int linear_index_int) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t face_index = linear_index % num_face_points;
            const size_t element_index = linear_index / num_face_points;
            const size_t local_volume_index_in_element =
                face_to_volume_index(face_index);
            const size_t local_volume_index =
                element_point_offsets_device(element_index) +
                local_volume_index_in_element;

            tnsr::i<double, volume_dim, Frame::Inertial> normal_covector{};
            double normal_magnitude_squared = 0.0;
            for (size_t i = 0; i < volume_dim; ++i) {
              const double component =
                  outward * inverse_jacobian(element_index,
                                             local_volume_index_in_element,
                                             d * volume_dim + i);
              normal_covector.get(i) = component;
              normal_magnitude_squared += component * component;
            }
            const double normal_magnitude = sqrt(normal_magnitude_squared);
            for (size_t i = 0; i < volume_dim; ++i) {
              normal_covector.get(i) /= normal_magnitude;
            }

            const auto vars_at_s = make_at_index(vars, local_volume_index);
            const auto gamma2_at_s = make_at_index(gamma2, local_volume_index);
            Scalar<double> v_psi{};
            tnsr::i<double, volume_dim, Frame::Inertial> v_zero{};
            Scalar<double> v_plus{};
            Scalar<double> v_minus{};
            tnsr::i<double, volume_dim, Frame::Inertial> normal_times_v_plus{};
            tnsr::i<double, volume_dim, Frame::Inertial> normal_times_v_minus{};
            Scalar<double> gamma2_v_psi{};
            tnsr::i<double, volume_dim, Frame::Inertial> char_speeds{};
            (void)ScalarWave::BoundaryCorrections::detail::dg_package_data_impl(
                make_not_null(&v_psi), make_not_null(&v_zero),
                make_not_null(&v_plus), make_not_null(&v_minus),
                make_not_null(&normal_times_v_plus),
                make_not_null(&normal_times_v_minus),
                make_not_null(&gamma2_v_psi), make_not_null(&char_speeds),
                get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Psi>>>(
                    vars_at_s),
                get<::Tags::AtIndex<::Tags::MirrorView<ScalarWave::Tags::Pi>>>(
                    vars_at_s),
                get<::Tags::AtIndex<
                    ::Tags::MirrorView<ScalarWave::Tags::Phi<volume_dim>>>>(
                    vars_at_s),
                gamma2_at_s, normal_covector,
                static_cast<const Scalar<double>*>(nullptr));

            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<0>>>(
                    packaged_face_data)),
                v_psi, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<1>>>(
                    packaged_face_data)),
                v_zero, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<2>>>(
                    packaged_face_data)),
                v_plus, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<3>>>(
                    packaged_face_data)),
                v_minus, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<4>>>(
                    packaged_face_data)),
                normal_times_v_plus, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<5>>>(
                    packaged_face_data)),
                normal_times_v_minus, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<6>>>(
                    packaged_face_data)),
                gamma2_v_psi, linear_index);
            set_at_index(
                make_not_null(&get<::Tags::MirrorView<package_field_tag<7>>>(
                    packaged_face_data)),
                char_speeds, linear_index);
          });
    }
  }
}

}  // namespace ScalarWave::Actions
