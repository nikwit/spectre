// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Executables/ScalarWave/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

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
  using device_constraint_gamma2_tag = KokkosTags::DeviceConstraintGamma2;
  using device_constraint_gamma2_space =
      typename device_constraint_gamma2_tag::type::value_type::memory_space;
  using device_face_to_volume_index_map_tag =
      KokkosTags::DeviceFaceToVolumeIndexMap<volume_dim>;

 public:
  using simple_tags = tmpl::list<device_inverse_jacobian_tag,
                                 device_constraint_gamma2_tag,
                                 device_face_to_volume_index_map_tag>;
  using return_tags = simple_tags;
  using argument_tags =
      tmpl::list<host_inverse_jacobian_tag, ScalarWave::Tags::ConstraintGamma2,
                 domain::Tags::Mesh<volume_dim>>;

  static void apply(
      const gsl::not_null<typename device_inverse_jacobian_tag::type*>
          device_inverse_jacobian,
      const gsl::not_null<typename device_constraint_gamma2_tag::type*>
          device_constraint_gamma2,
      const gsl::not_null<typename device_face_to_volume_index_map_tag::type*>
          device_face_to_volume_index_map,
      const typename host_inverse_jacobian_tag::type& host_inverse_jacobian,
      const ScalarWave::Tags::ConstraintGamma2::type& host_constraint_gamma2,
      const Mesh<volume_dim>& mesh) {
    *device_inverse_jacobian = copy_to_device(
        host_inverse_jacobian, tmpl::type_<device_inverse_jacobian_space>{});
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
      lower_face_to_volume_index =
          Kokkos::View<size_t*>("KokkosLowerFaceToVolumeIndexMap",
                                num_face_points);
      upper_face_to_volume_index =
          Kokkos::View<size_t*>("KokkosUpperFaceToVolumeIndexMap",
                                num_face_points);

      auto host_lower_face_to_volume_index =
          Kokkos::create_mirror_view(lower_face_to_volume_index);
      auto host_upper_face_to_volume_index =
          Kokkos::create_mirror_view(upper_face_to_volume_index);
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
      Kokkos::deep_copy(lower_face_to_volume_index,
                        host_lower_face_to_volume_index);
      Kokkos::deep_copy(upper_face_to_volume_index,
                        host_upper_face_to_volume_index);
    }
  }
};

}  // namespace ScalarWave::Actions

