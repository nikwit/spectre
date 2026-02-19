// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarInfo.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename System>
struct InitializeKokkosTags {
 private:
  static constexpr size_t volume_dim = System::volume_dim;
  using host_inverse_jacobian_tag =
      domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                    Frame::Inertial>;
  using device_inverse_jacobian_tag =
      KokkosTags::DeviceInverseJacobian<volume_dim>;
  using device_inverse_jacobian_space =
      typename device_inverse_jacobian_tag::type::value_type::memory_space;

  using device_constraint_gamma0_tag = KokkosTags::DeviceConstraintGamma0;
  using device_constraint_gamma1_tag = KokkosTags::DeviceConstraintGamma1;
  using device_constraint_gamma2_tag = KokkosTags::DeviceConstraintGamma2;
  using device_constraint_gamma0_space =
      typename device_constraint_gamma0_tag::type::value_type::memory_space;
  using device_constraint_gamma1_space =
      typename device_constraint_gamma1_tag::type::value_type::memory_space;
  using device_constraint_gamma2_space =
      typename device_constraint_gamma2_tag::type::value_type::memory_space;

  using device_face_to_volume_index_map_tag =
      KokkosTags::DeviceFaceToVolumeIndexMap<volume_dim>;
  using device_face_unit_normal_covector_tag =
      KokkosTags::DeviceFaceUnitNormalCovector<volume_dim>;
  using device_face_normal_magnitude_tag =
      KokkosTags::DeviceFaceNormalMagnitude<volume_dim>;
  using device_mortar_data_tag = KokkosTags::DeviceMortarData<volume_dim>;
  using mortar_mesh_tag = evolution::dg::Tags::MortarMesh<volume_dim>;
  using mortar_info_tag = evolution::dg::Tags::MortarInfo<volume_dim>;

 public:
  using simple_tags = tmpl::list<
      device_inverse_jacobian_tag, device_constraint_gamma0_tag,
      device_constraint_gamma1_tag, device_constraint_gamma2_tag,
      device_face_to_volume_index_map_tag, device_face_unit_normal_covector_tag,
      device_face_normal_magnitude_tag, device_mortar_data_tag>;
  using return_tags = simple_tags;
  using argument_tags = tmpl::list<
      host_inverse_jacobian_tag, gh::Tags::ConstraintGamma0,
      gh::Tags::ConstraintGamma1, gh::Tags::ConstraintGamma2,
      domain::Tags::Mesh<volume_dim>, domain::Tags::Element<volume_dim>,
      mortar_mesh_tag, mortar_info_tag>;

  static void apply(
      const gsl::not_null<typename device_inverse_jacobian_tag::type*>
          device_inverse_jacobian,
      const gsl::not_null<typename device_constraint_gamma0_tag::type*>
          device_constraint_gamma0,
      const gsl::not_null<typename device_constraint_gamma1_tag::type*>
          device_constraint_gamma1,
      const gsl::not_null<typename device_constraint_gamma2_tag::type*>
          device_constraint_gamma2,
      const gsl::not_null<typename device_face_to_volume_index_map_tag::type*>
          device_face_to_volume_index_map,
      const gsl::not_null<
          typename device_face_unit_normal_covector_tag::type*>
          device_face_unit_normal_covector,
      const gsl::not_null<typename device_face_normal_magnitude_tag::type*>
          device_face_normal_magnitude,
      const gsl::not_null<typename device_mortar_data_tag::type*>
          device_mortar_data,
      const typename host_inverse_jacobian_tag::type& host_inverse_jacobian,
      const gh::Tags::ConstraintGamma0::type& host_constraint_gamma0,
      const gh::Tags::ConstraintGamma1::type& host_constraint_gamma1,
      const gh::Tags::ConstraintGamma2::type& host_constraint_gamma2,
      const Mesh<volume_dim>& mesh, const Element<volume_dim>& element,
      const typename mortar_mesh_tag::type& mortar_meshes,
      const typename mortar_info_tag::type& mortar_infos) {
    *device_inverse_jacobian = copy_to_device(
        host_inverse_jacobian, tmpl::type_<device_inverse_jacobian_space>{});
    *device_constraint_gamma0 = copy_to_device(
        host_constraint_gamma0, tmpl::type_<device_constraint_gamma0_space>{});
    *device_constraint_gamma1 = copy_to_device(
        host_constraint_gamma1, tmpl::type_<device_constraint_gamma1_space>{});
    *device_constraint_gamma2 = copy_to_device(
        host_constraint_gamma2, tmpl::type_<device_constraint_gamma2_space>{});

    const auto [indices_buffer, volume_and_slice_index_map] =
        volume_and_slice_indices(mesh.extents());
    (void)indices_buffer;

    for (size_t d = 0; d < volume_dim; ++d) {
      const size_t num_face_points = mesh.extents().slice_away(d).product();
      auto& lower_face_to_volume_index =
          gsl::at(*device_face_to_volume_index_map, d).first;
      auto& upper_face_to_volume_index =
          gsl::at(*device_face_to_volume_index_map, d).second;
      auto& lower_face_unit_normal_covector =
          gsl::at(*device_face_unit_normal_covector, d).first;
      auto& upper_face_unit_normal_covector =
          gsl::at(*device_face_unit_normal_covector, d).second;
      auto& lower_face_normal_magnitude =
          gsl::at(*device_face_normal_magnitude, d).first;
      auto& upper_face_normal_magnitude =
          gsl::at(*device_face_normal_magnitude, d).second;

      lower_face_to_volume_index =
          ::Kokkos::View<size_t*>("GhKokkosLowerFaceToVolumeIndexMap",
                                  num_face_points);
      upper_face_to_volume_index =
          ::Kokkos::View<size_t*>("GhKokkosUpperFaceToVolumeIndexMap",
                                  num_face_points);
      lower_face_unit_normal_covector =
          ::Kokkos::View<double**>("GhKokkosLowerFaceUnitNormalCovector",
                                   num_face_points, volume_dim);
      upper_face_unit_normal_covector =
          ::Kokkos::View<double**>("GhKokkosUpperFaceUnitNormalCovector",
                                   num_face_points, volume_dim);
      lower_face_normal_magnitude =
          ::Kokkos::View<double*>("GhKokkosLowerFaceNormalMagnitude",
                                  num_face_points);
      upper_face_normal_magnitude =
          ::Kokkos::View<double*>("GhKokkosUpperFaceNormalMagnitude",
                                  num_face_points);

      auto host_lower_face_to_volume_index =
          ::Kokkos::create_mirror_view(lower_face_to_volume_index);
      auto host_upper_face_to_volume_index =
          ::Kokkos::create_mirror_view(upper_face_to_volume_index);
      auto host_lower_face_unit_normal_covector =
          ::Kokkos::create_mirror_view(lower_face_unit_normal_covector);
      auto host_upper_face_unit_normal_covector =
          ::Kokkos::create_mirror_view(upper_face_unit_normal_covector);
      auto host_lower_face_normal_magnitude =
          ::Kokkos::create_mirror_view(lower_face_normal_magnitude);
      auto host_upper_face_normal_magnitude =
          ::Kokkos::create_mirror_view(upper_face_normal_magnitude);

      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).first) {
        host_lower_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }
      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).second) {
        host_upper_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }

      for (size_t face_index = 0; face_index < num_face_points; ++face_index) {
        for (const bool upper_side : {false, true}) {
          const auto& face_to_volume_index =
              upper_side ? host_upper_face_to_volume_index
                         : host_lower_face_to_volume_index;
          auto& face_unit_normal_covector =
              upper_side ? host_upper_face_unit_normal_covector
                         : host_lower_face_unit_normal_covector;
          auto& face_normal_magnitude =
              upper_side ? host_upper_face_normal_magnitude
                         : host_lower_face_normal_magnitude;
          const double outward_sign = upper_side ? 1.0 : -1.0;
          const size_t volume_index = face_to_volume_index(face_index);

          double normal_magnitude_squared = 0.0;
          for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
            const double unnormalized_normal_component =
                outward_sign *
                host_inverse_jacobian.get(d, inertial_d)[volume_index];
            face_unit_normal_covector(face_index, inertial_d) =
                unnormalized_normal_component;
            normal_magnitude_squared +=
                unnormalized_normal_component * unnormalized_normal_component;
          }
          const double normal_magnitude = std::sqrt(normal_magnitude_squared);
          face_normal_magnitude(face_index) = normal_magnitude;
          for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
            face_unit_normal_covector(face_index, inertial_d) /=
                normal_magnitude;
          }
        }
      }

      ::Kokkos::deep_copy(lower_face_to_volume_index,
                          host_lower_face_to_volume_index);
      ::Kokkos::deep_copy(upper_face_to_volume_index,
                          host_upper_face_to_volume_index);
      ::Kokkos::deep_copy(lower_face_unit_normal_covector,
                          host_lower_face_unit_normal_covector);
      ::Kokkos::deep_copy(upper_face_unit_normal_covector,
                          host_upper_face_unit_normal_covector);
      ::Kokkos::deep_copy(lower_face_normal_magnitude,
                          host_lower_face_normal_magnitude);
      ::Kokkos::deep_copy(upper_face_normal_magnitude,
                          host_upper_face_normal_magnitude);
    }

    device_mortar_data->clear();
    for (const auto& [direction, neighbors] : element.neighbors()) {
      const size_t sliced_dim = direction.dimension();
      const Mesh<volume_dim - 1> face_mesh = mesh.slice_away(sliced_dim);
      for (const auto& neighbor : neighbors) {
        const DirectionalId<volume_dim> mortar_id{direction, neighbor};
        const auto& mortar_mesh = mortar_meshes.at(mortar_id);
        const auto& mortar_size = mortar_infos.at(mortar_id).mortar_size();
        const auto& orientation = neighbors.orientation(neighbor);

        std::vector<double> grid_point_indices(
            mortar_mesh.number_of_grid_points());
        std::iota(grid_point_indices.begin(), grid_point_indices.end(), 0.0);
        const std::vector<double> oriented_grid_point_indices =
            orient_variables_on_slice(grid_point_indices, mortar_mesh.extents(),
                                      sliced_dim, orientation);

        KokkosTags::MortarData<volume_dim> mortar_data{};
        mortar_data.needs_projection =
            Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size);
        mortar_data.mortar_size = mortar_size;
        mortar_data.oriented_mortar_grid_point_source_index = ::Kokkos::View<
            size_t*>("GhKokkosOrientedMortarGridPointSourceIndex",
                     mortar_mesh.number_of_grid_points());
        auto host_oriented_mortar_grid_point_source_index =
            ::Kokkos::create_mirror_view(
                mortar_data.oriented_mortar_grid_point_source_index);
        for (size_t i = 0; i < oriented_grid_point_indices.size(); ++i) {
          host_oriented_mortar_grid_point_source_index(i) =
              static_cast<size_t>(oriented_grid_point_indices[i]);
        }
        ::Kokkos::deep_copy(mortar_data.oriented_mortar_grid_point_source_index,
                            host_oriented_mortar_grid_point_source_index);

        device_mortar_data->insert_or_assign(mortar_id, std::move(mortar_data));
      }
    }
  }
};

}  // namespace gh::Actions
