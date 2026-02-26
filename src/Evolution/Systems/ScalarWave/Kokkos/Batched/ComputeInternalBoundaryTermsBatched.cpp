// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/Kokkos/Batched/ComputeInternalBoundaryTermsBatched.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {
namespace {

static constexpr size_t volume_dim = 3;
static constexpr size_t number_of_faces = 2 * volume_dim;
using package_field_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::package_field_tags;
template <size_t I>
using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
using device_package_field_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::device_package_field_tags;
using device_dt_boundary_tags = evolution::Kokkos::PackedBoundaryScratch<
    ScalarWave::System<3>>::device_dt_boundary_tags;
using face_boundary_metadata = evolution::Kokkos::Batched::FaceBoundaryMetadata;
using projection_group_metadata =
    evolution::Kokkos::Batched::ProjectionGroupMetadata;

constexpr size_t side_index(const Side side) {
  return side == Side::Upper ? static_cast<size_t>(1) : static_cast<size_t>(0);
}

constexpr size_t face_index(const size_t sliced_dim, const size_t side_i) {
  return 2 * sliced_dim + side_i;
}

size_t face_index(const Direction<volume_dim>& direction) {
  return face_index(direction.dimension(), side_index(direction.side()));
}

template <typename TagsList>
void apply_tensor_product_projection_2d_batched(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& input, const size_t num_work_items,
    const size_t source_points_dim_0, const size_t source_points_dim_1,
    const size_t target_points_dim_0, const size_t target_points_dim_1,
    const MatrixViewRO& matrix_dim_0, const MatrixViewRO& matrix_dim_1) {
  const auto input_view = input.view();
  auto result_view = result->view();
  const size_t num_components = input_view.extent(1);
  const size_t source_points_per_work_item =
      source_points_dim_0 * source_points_dim_1;
  const size_t target_points_per_work_item =
      target_points_dim_0 * target_points_dim_1;

  ASSERT(input.number_of_grid_points() ==
             num_work_items * source_points_per_work_item,
         "Input size mismatch in batched 2D projection.");
  ASSERT(result->number_of_grid_points() ==
             num_work_items * target_points_per_work_item,
         "Result size mismatch in batched 2D projection.");

  ::Kokkos::View<double***> projected_dim_0(
      "ComputeInternalBoundaryTermsBatchedProjectedDim0", num_work_items,
      target_points_dim_0 * source_points_dim_1, num_components);
  ::Kokkos::parallel_for(
      "ComputeInternalBoundaryTermsBatchedProjectDim0Batched",
      ::Kokkos::RangePolicy<size_t>{0, num_work_items * target_points_dim_0 *
                                           source_points_dim_1 *
                                           num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i1 = remaining % source_points_dim_1;
        remaining /= source_points_dim_1;
        const size_t i0 = remaining % target_points_dim_0;
        const size_t work_item_index = remaining / target_points_dim_0;

        double sum = 0.0;
        for (size_t k0 = 0; k0 < source_points_dim_0; ++k0) {
          const size_t source_point =
              k0 + source_points_dim_0 * i1 +
              work_item_index * source_points_per_work_item;
          sum += matrix_dim_0(i0, k0) * input_view(source_point, component);
        }
        projected_dim_0(work_item_index, i0 + target_points_dim_0 * i1,
                        component) = sum;
      });

  ::Kokkos::parallel_for(
      "ComputeInternalBoundaryTermsBatchedProjectDim1Batched",
      ::Kokkos::RangePolicy<size_t>{0, num_work_items * target_points_dim_0 *
                                           target_points_dim_1 *
                                           num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i1 = remaining % target_points_dim_1;
        remaining /= target_points_dim_1;
        const size_t i0 = remaining % target_points_dim_0;
        const size_t work_item_index = remaining / target_points_dim_0;

        double sum = 0.0;
        for (size_t k1 = 0; k1 < source_points_dim_1; ++k1) {
          sum += matrix_dim_1(i1, k1) *
                 projected_dim_0(work_item_index, i0 + target_points_dim_0 * k1,
                                 component);
        }
        const size_t target_point =
            i0 + target_points_dim_0 * i1 +
            work_item_index * target_points_per_work_item;
        result_view(target_point, component) = sum;
      });
}

template <typename TagsList>
void project_to_mortar_device_batched(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& vars, const size_t num_work_items,
    const Mesh<2>& face_mesh, const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size) {
  if (not Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size)) {
    ASSERT(result->number_of_grid_points() == vars.number_of_grid_points(),
           "Cannot skip batched projection for incompatible result size.");
    ::Kokkos::deep_copy(result->view(), vars.view());
    return;
  }

  const auto& matrix_dim_0 =
      Spectral::projection_matrix_parent_to_child_on_device(
          face_mesh.slice_through(0), mortar_mesh.slice_through(0),
          gsl::at(mortar_size, 0));
  const auto& matrix_dim_1 =
      Spectral::projection_matrix_parent_to_child_on_device(
          face_mesh.slice_through(1), mortar_mesh.slice_through(1),
          gsl::at(mortar_size, 1));
  apply_tensor_product_projection_2d_batched(
      result, vars, num_work_items, face_mesh.extents(0), face_mesh.extents(1),
      mortar_mesh.extents(0), mortar_mesh.extents(1), matrix_dim_0,
      matrix_dim_1);
}

template <typename TagsList>
void project_from_mortar_device_batched(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& vars, const size_t num_work_items,
    const Mesh<2>& face_mesh, const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size) {
  ASSERT(Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size),
         "project_from_mortar_device_batched should not be called when no "
         "projection is needed.");

  const auto& matrix_dim_0 =
      Spectral::projection_matrix_child_to_parent_on_device(
          mortar_mesh.slice_through(0), face_mesh.slice_through(0),
          gsl::at(mortar_size, 0));
  const auto& matrix_dim_1 =
      Spectral::projection_matrix_child_to_parent_on_device(
          mortar_mesh.slice_through(1), face_mesh.slice_through(1),
          gsl::at(mortar_size, 1));
  apply_tensor_product_projection_2d_batched(
      result, vars, num_work_items, mortar_mesh.extents(0),
      mortar_mesh.extents(1), face_mesh.extents(0), face_mesh.extents(1),
      matrix_dim_0, matrix_dim_1);
}

}  // namespace

void ComputeInternalBoundaryTermsBatched::apply(
    const gsl::not_null<evolution::Kokkos::Tags::PackedBoundaryScratch<
        ScalarWave::System<3>>::type*>
        packed_boundary_scratch,
    const evolution::Kokkos::Tags::PackedTopology<ScalarWave::System<3>>::type&
        packed_topology,
    const evolution::Kokkos::Tags::PackedBoundaryMetadata<
        ScalarWave::System<3>>::type& packed_boundary_metadata) {
  if (packed_topology.total_points == 0) {
    return;
  }

  ASSERT(packed_topology.uniform_quadrature_host[0] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[1] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[2] ==
                 Spectral::Quadrature::GaussLobatto,
         "ComputeInternalBoundaryTermsBatched currently supports "
         "Gauss-Lobatto quadrature only.");

  const auto& extents = packed_topology.uniform_extents_host;
  ASSERT(extents[0] == extents[1] and extents[1] == extents[2],
         "ComputeInternalBoundaryTermsBatched assumes uniform p across "
         "dimensions (equal extents in all dimensions).");
  const Mesh<volume_dim> mesh{extents, packed_topology.uniform_basis_host,
                              packed_topology.uniform_quadrature_host};
  const Mesh<volume_dim - 1> uniform_face_mesh = mesh.slice_away(0);
  const size_t num_face_points = uniform_face_mesh.number_of_grid_points();
  const auto& elements = packed_topology.local_elements;

  const auto& packaged_face_data_for_all_elements =
      packed_boundary_scratch->packaged_face_data_for_all_elements;
  for (size_t face_id = 0; face_id < number_of_faces; ++face_id) {
    ASSERT(
        packaged_face_data_for_all_elements[face_id].number_of_grid_points() ==
            elements.size() * num_face_points,
        "PackageLocalFacesBatched must run before computing internal batched "
        "boundary terms.");
  }
  const auto packaged_face_data_view_0_0 =
      packaged_face_data_for_all_elements[face_index(0, 0)].view();
  const auto packaged_face_data_view_0_1 =
      packaged_face_data_for_all_elements[face_index(0, 1)].view();
  const auto packaged_face_data_view_1_0 =
      packaged_face_data_for_all_elements[face_index(1, 0)].view();
  const auto packaged_face_data_view_1_1 =
      packaged_face_data_for_all_elements[face_index(1, 1)].view();
  const auto packaged_face_data_view_2_0 =
      packaged_face_data_for_all_elements[face_index(2, 0)].view();
  const auto packaged_face_data_view_2_1 =
      packaged_face_data_for_all_elements[face_index(2, 1)].view();

  auto& dt_face_sum_for_all_elements =
      packed_boundary_scratch->internal_boundary_terms_for_all_elements;
  for (size_t face_id = 0; face_id < number_of_faces; ++face_id) {
    if (dt_face_sum_for_all_elements[face_id].number_of_grid_points() !=
        elements.size() * num_face_points) {
      dt_face_sum_for_all_elements[face_id] =
          Variables<device_dt_boundary_tags>{elements.size() * num_face_points};
    }
  }

  for (size_t d = 0; d < volume_dim; ++d) {
    const Mesh<volume_dim - 1>& face_mesh = uniform_face_mesh;
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      auto& dt_face_sum_all =
          dt_face_sum_for_all_elements[face_index(d, side_i)];
      ::Kokkos::deep_copy(dt_face_sum_all.view(), 0.0);

      const face_boundary_metadata& cached_face_metadata =
          packed_boundary_metadata
              .boundary_correction_face_metadata[face_index(d, side_i)];

      const auto& local_packaged_face_data_all =
          packaged_face_data_for_all_elements[face_index(d, side_i)];
      const auto local_v_psi_all =
          get<::Tags::MirrorView<package_field_tag<0>>>(
              local_packaged_face_data_all);
      const auto local_v_zero_all =
          get<::Tags::MirrorView<package_field_tag<1>>>(
              local_packaged_face_data_all);
      const auto local_v_plus_all =
          get<::Tags::MirrorView<package_field_tag<2>>>(
              local_packaged_face_data_all);
      const auto local_v_minus_all =
          get<::Tags::MirrorView<package_field_tag<3>>>(
              local_packaged_face_data_all);
      const auto local_normal_times_v_plus_all =
          get<::Tags::MirrorView<package_field_tag<4>>>(
              local_packaged_face_data_all);
      const auto local_normal_times_v_minus_all =
          get<::Tags::MirrorView<package_field_tag<5>>>(
              local_packaged_face_data_all);
      const auto local_gamma2_v_psi_all =
          get<::Tags::MirrorView<package_field_tag<6>>>(
              local_packaged_face_data_all);
      const auto local_char_speeds_all =
          get<::Tags::MirrorView<package_field_tag<7>>>(
              local_packaged_face_data_all);
      const auto dt_psi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
              dt_face_sum_all);
      const auto dt_pi_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
              dt_face_sum_all);
      const auto dt_phi_face_sum_all = get<
          ::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
          dt_face_sum_all);
      const auto& mortar_work_items =
          cached_face_metadata.no_projection_work_items;
      const auto& oriented_remote_face_indices =
          cached_face_metadata.no_projection_oriented_remote_face_indices;
      if (mortar_work_items.extent(0) > 0) {
        ::Kokkos::parallel_for(
            "ComputeInternalBoundaryTermsBatchedNoProjectionMortars",
            mortar_work_items.extent(0) * num_face_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              const size_t face_index_on_face = linear_index % num_face_points;
              const size_t mortar_work_item_index =
                  linear_index / num_face_points;
              const auto mortar_work_item =
                  mortar_work_items(mortar_work_item_index);

              const size_t local_linear_index =
                  mortar_work_item.local_element_index * num_face_points +
                  face_index_on_face;
              tnsr::i<double, volume_dim, Frame::Inertial> local_v_zero{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus{};
              tnsr::i<double, volume_dim, Frame::Inertial> local_char_speeds{};
              for (size_t i = 0; i < volume_dim; ++i) {
                local_v_zero.get(i) =
                    local_v_zero_all.get(i)[local_linear_index];
                local_normal_times_v_plus.get(i) =
                    local_normal_times_v_plus_all.get(i)[local_linear_index];
                local_normal_times_v_minus.get(i) =
                    local_normal_times_v_minus_all.get(i)[local_linear_index];
                local_char_speeds.get(i) =
                    local_char_speeds_all.get(i)[local_linear_index];
              }
              Scalar<double> local_v_psi{};
              Scalar<double> local_v_plus{};
              Scalar<double> local_v_minus{};
              Scalar<double> local_gamma2_v_psi{};
              get(local_v_psi) = get(local_v_psi_all)[local_linear_index];
              get(local_v_plus) = get(local_v_plus_all)[local_linear_index];
              get(local_v_minus) = get(local_v_minus_all)[local_linear_index];
              get(local_gamma2_v_psi) =
                  get(local_gamma2_v_psi_all)[local_linear_index];

              const size_t remote_face_index = oriented_remote_face_indices(
                  mortar_work_item.oriented_remote_face_index_offset +
                  face_index_on_face);
              const size_t remote_linear_index =
                  mortar_work_item.remote_element_index * num_face_points +
                  remote_face_index;
              const auto& remote_packaged_face_data_all =
                  packaged_face_data_for_all_elements[mortar_work_item
                                                          .remote_face_id];
              const auto remote_v_psi_all =
                  get<::Tags::MirrorView<package_field_tag<0>>>(
                      remote_packaged_face_data_all);
              const auto remote_v_zero_all =
                  get<::Tags::MirrorView<package_field_tag<1>>>(
                      remote_packaged_face_data_all);
              const auto remote_v_plus_all =
                  get<::Tags::MirrorView<package_field_tag<2>>>(
                      remote_packaged_face_data_all);
              const auto remote_v_minus_all =
                  get<::Tags::MirrorView<package_field_tag<3>>>(
                      remote_packaged_face_data_all);
              const auto remote_normal_times_v_plus_all =
                  get<::Tags::MirrorView<package_field_tag<4>>>(
                      remote_packaged_face_data_all);
              const auto remote_normal_times_v_minus_all =
                  get<::Tags::MirrorView<package_field_tag<5>>>(
                      remote_packaged_face_data_all);
              const auto remote_gamma2_v_psi_all =
                  get<::Tags::MirrorView<package_field_tag<6>>>(
                      remote_packaged_face_data_all);
              const auto remote_char_speeds_all =
                  get<::Tags::MirrorView<package_field_tag<7>>>(
                      remote_packaged_face_data_all);
              tnsr::i<double, volume_dim, Frame::Inertial> remote_v_zero{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_plus{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_minus{};
              tnsr::i<double, volume_dim, Frame::Inertial> remote_char_speeds{};
              for (size_t i = 0; i < volume_dim; ++i) {
                remote_v_zero.get(i) =
                    remote_v_zero_all.get(i)[remote_linear_index];
                remote_normal_times_v_plus.get(i) =
                    remote_normal_times_v_plus_all.get(i)[remote_linear_index];
                remote_normal_times_v_minus.get(i) =
                    remote_normal_times_v_minus_all.get(i)[remote_linear_index];
                remote_char_speeds.get(i) =
                    remote_char_speeds_all.get(i)[remote_linear_index];
              }
              Scalar<double> remote_v_psi{};
              Scalar<double> remote_v_plus{};
              Scalar<double> remote_v_minus{};
              Scalar<double> remote_gamma2_v_psi{};
              get(remote_v_psi) = get(remote_v_psi_all)[remote_linear_index];
              get(remote_v_plus) = get(remote_v_plus_all)[remote_linear_index];
              get(remote_v_minus) =
                  get(remote_v_minus_all)[remote_linear_index];
              get(remote_gamma2_v_psi) =
                  get(remote_gamma2_v_psi_all)[remote_linear_index];

              Scalar<double> dt_psi_at_face{};
              Scalar<double> dt_pi_at_face{};
              tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_at_face{};
              ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
                  volume_dim, double>(
                  make_not_null(&dt_psi_at_face), make_not_null(&dt_pi_at_face),
                  make_not_null(&dt_phi_at_face), local_v_psi, local_v_zero,
                  local_v_plus, local_v_minus, local_normal_times_v_plus,
                  local_normal_times_v_minus, local_gamma2_v_psi,
                  local_char_speeds, remote_v_psi, remote_v_zero, remote_v_plus,
                  remote_v_minus, remote_normal_times_v_plus,
                  remote_normal_times_v_minus, remote_gamma2_v_psi,
                  remote_char_speeds);

              ::Kokkos::atomic_add(
                  &get(dt_psi_face_sum_all)[local_linear_index],
                  get(dt_psi_at_face));
              ::Kokkos::atomic_add(&get(dt_pi_face_sum_all)[local_linear_index],
                                   get(dt_pi_at_face));
              for (size_t i = 0; i < volume_dim; ++i) {
                ::Kokkos::atomic_add(
                    &dt_phi_face_sum_all.get(i)[local_linear_index],
                    dt_phi_at_face.get(i));
              }
            });
      }

      for (const projection_group_metadata& projection_group :
           cached_face_metadata.projection_groups) {
        if (projection_group.work_items.extent(0) == 0) {
          continue;
        }

        const size_t num_projection_work_items =
            projection_group.work_items.extent(0);
        const size_t num_mortar_points =
            projection_group.mortar_mesh.number_of_grid_points();
        const auto& projection_work_items = projection_group.work_items;
        const auto& projection_oriented_remote_face_indices =
            projection_group.oriented_remote_face_indices;

        Variables<device_package_field_tags> local_packaged_face_data_batched{
            num_projection_work_items * num_face_points};
        Variables<device_package_field_tags> remote_packaged_face_data_batched{
            num_projection_work_items * num_face_points};
        const auto local_packaged_face_data_all_view =
            local_packaged_face_data_all.view();
        auto local_packaged_face_data_batched_view =
            local_packaged_face_data_batched.view();
        auto remote_packaged_face_data_batched_view =
            remote_packaged_face_data_batched.view();
        const size_t num_package_components =
            local_packaged_face_data_batched_view.extent(1);
        ::Kokkos::parallel_for(
            "ComputeInternalBoundaryTermsBatchedGatherProjectionLocalFaces",
            num_projection_work_items * num_face_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              const size_t face_index_on_face = linear_index % num_face_points;
              const size_t work_item_index = linear_index / num_face_points;
              const auto work_item = projection_work_items(work_item_index);
              const size_t source_linear_index =
                  work_item.local_element_index * num_face_points +
                  face_index_on_face;
              for (size_t component = 0; component < num_package_components;
                   ++component) {
                local_packaged_face_data_batched_view(linear_index, component) =
                    local_packaged_face_data_all_view(source_linear_index,
                                                      component);
              }
            });
        ::Kokkos::parallel_for(
            "ComputeInternalBoundaryTermsBatchedGatherProjectionRemoteFaces",
            num_projection_work_items * num_face_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              const size_t face_index_on_face = linear_index % num_face_points;
              const size_t work_item_index = linear_index / num_face_points;
              const auto work_item = projection_work_items(work_item_index);
              const size_t remote_face_index =
                  projection_oriented_remote_face_indices(
                      work_item.oriented_remote_face_index_offset +
                      face_index_on_face);
              const size_t source_linear_index =
                  work_item.remote_element_index * num_face_points +
                  remote_face_index;
              auto remote_packaged_face_data_all_view =
                  packaged_face_data_view_2_1;
              if (work_item.remote_face_id == 0) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_0_0;
              } else if (work_item.remote_face_id == 1) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_0_1;
              } else if (work_item.remote_face_id == 2) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_1_0;
              } else if (work_item.remote_face_id == 3) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_1_1;
              } else if (work_item.remote_face_id == 4) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_2_0;
              }
              for (size_t component = 0; component < num_package_components;
                   ++component) {
                remote_packaged_face_data_batched_view(linear_index,
                                                       component) =
                    remote_packaged_face_data_all_view(source_linear_index,
                                                       component);
              }
            });

        Variables<device_package_field_tags> local_packaged_mortar_data_batched{
            num_projection_work_items * num_mortar_points};
        Variables<device_package_field_tags>
            remote_packaged_mortar_data_batched{num_projection_work_items *
                                                num_mortar_points};
        project_to_mortar_device_batched(
            make_not_null(&local_packaged_mortar_data_batched),
            local_packaged_face_data_batched, num_projection_work_items,
            face_mesh, projection_group.mortar_mesh,
            projection_group.mortar_size);
        project_to_mortar_device_batched(
            make_not_null(&remote_packaged_mortar_data_batched),
            remote_packaged_face_data_batched, num_projection_work_items,
            face_mesh, projection_group.mortar_mesh,
            projection_group.mortar_size);

        Variables<device_dt_boundary_tags> dt_boundary_on_mortar_batched{
            num_projection_work_items * num_mortar_points};
        const auto local_v_psi_on_mortar =
            get<::Tags::MirrorView<package_field_tag<0>>>(
                local_packaged_mortar_data_batched);
        const auto local_v_zero_on_mortar =
            get<::Tags::MirrorView<package_field_tag<1>>>(
                local_packaged_mortar_data_batched);
        const auto local_v_plus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<2>>>(
                local_packaged_mortar_data_batched);
        const auto local_v_minus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<3>>>(
                local_packaged_mortar_data_batched);
        const auto local_normal_times_v_plus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<4>>>(
                local_packaged_mortar_data_batched);
        const auto local_normal_times_v_minus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<5>>>(
                local_packaged_mortar_data_batched);
        const auto local_gamma2_v_psi_on_mortar =
            get<::Tags::MirrorView<package_field_tag<6>>>(
                local_packaged_mortar_data_batched);
        const auto local_char_speeds_on_mortar =
            get<::Tags::MirrorView<package_field_tag<7>>>(
                local_packaged_mortar_data_batched);
        const auto remote_v_psi_on_mortar =
            get<::Tags::MirrorView<package_field_tag<0>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_v_zero_on_mortar =
            get<::Tags::MirrorView<package_field_tag<1>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_v_plus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<2>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_v_minus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<3>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_normal_times_v_plus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<4>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_normal_times_v_minus_on_mortar =
            get<::Tags::MirrorView<package_field_tag<5>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_gamma2_v_psi_on_mortar =
            get<::Tags::MirrorView<package_field_tag<6>>>(
                remote_packaged_mortar_data_batched);
        const auto remote_char_speeds_on_mortar =
            get<::Tags::MirrorView<package_field_tag<7>>>(
                remote_packaged_mortar_data_batched);
        const auto dt_psi_on_mortar =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                dt_boundary_on_mortar_batched);
        const auto dt_pi_on_mortar =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                dt_boundary_on_mortar_batched);
        const auto dt_phi_on_mortar = get<
            ::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
            dt_boundary_on_mortar_batched);
        ::Kokkos::parallel_for(
            "ComputeInternalBoundaryTermsBatchedProjectionMortarBoundaryTerms",
            num_projection_work_items * num_mortar_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_v_zero_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  local_char_speeds_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_v_zero_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_plus_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_minus_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial>
                  remote_char_speeds_at_mortar{};
              for (size_t i = 0; i < volume_dim; ++i) {
                local_v_zero_at_mortar.get(i) =
                    local_v_zero_on_mortar.get(i)[linear_index];
                local_normal_times_v_plus_at_mortar.get(i) =
                    local_normal_times_v_plus_on_mortar.get(i)[linear_index];
                local_normal_times_v_minus_at_mortar.get(i) =
                    local_normal_times_v_minus_on_mortar.get(i)[linear_index];
                local_char_speeds_at_mortar.get(i) =
                    local_char_speeds_on_mortar.get(i)[linear_index];
                remote_v_zero_at_mortar.get(i) =
                    remote_v_zero_on_mortar.get(i)[linear_index];
                remote_normal_times_v_plus_at_mortar.get(i) =
                    remote_normal_times_v_plus_on_mortar.get(i)[linear_index];
                remote_normal_times_v_minus_at_mortar.get(i) =
                    remote_normal_times_v_minus_on_mortar.get(i)[linear_index];
                remote_char_speeds_at_mortar.get(i) =
                    remote_char_speeds_on_mortar.get(i)[linear_index];
              }

              Scalar<double> local_v_psi_at_mortar{};
              Scalar<double> local_v_plus_at_mortar{};
              Scalar<double> local_v_minus_at_mortar{};
              Scalar<double> local_gamma2_v_psi_at_mortar{};
              Scalar<double> remote_v_psi_at_mortar{};
              Scalar<double> remote_v_plus_at_mortar{};
              Scalar<double> remote_v_minus_at_mortar{};
              Scalar<double> remote_gamma2_v_psi_at_mortar{};
              get(local_v_psi_at_mortar) =
                  get(local_v_psi_on_mortar)[linear_index];
              get(local_v_plus_at_mortar) =
                  get(local_v_plus_on_mortar)[linear_index];
              get(local_v_minus_at_mortar) =
                  get(local_v_minus_on_mortar)[linear_index];
              get(local_gamma2_v_psi_at_mortar) =
                  get(local_gamma2_v_psi_on_mortar)[linear_index];
              get(remote_v_psi_at_mortar) =
                  get(remote_v_psi_on_mortar)[linear_index];
              get(remote_v_plus_at_mortar) =
                  get(remote_v_plus_on_mortar)[linear_index];
              get(remote_v_minus_at_mortar) =
                  get(remote_v_minus_on_mortar)[linear_index];
              get(remote_gamma2_v_psi_at_mortar) =
                  get(remote_gamma2_v_psi_on_mortar)[linear_index];

              Scalar<double> dt_psi_at_mortar{};
              Scalar<double> dt_pi_at_mortar{};
              tnsr::i<double, volume_dim, Frame::Inertial> dt_phi_at_mortar{};
              ScalarWave::BoundaryCorrections::detail::dg_boundary_terms_impl<
                  volume_dim, double>(
                  make_not_null(&dt_psi_at_mortar),
                  make_not_null(&dt_pi_at_mortar),
                  make_not_null(&dt_phi_at_mortar), local_v_psi_at_mortar,
                  local_v_zero_at_mortar, local_v_plus_at_mortar,
                  local_v_minus_at_mortar, local_normal_times_v_plus_at_mortar,
                  local_normal_times_v_minus_at_mortar,
                  local_gamma2_v_psi_at_mortar, local_char_speeds_at_mortar,
                  remote_v_psi_at_mortar, remote_v_zero_at_mortar,
                  remote_v_plus_at_mortar, remote_v_minus_at_mortar,
                  remote_normal_times_v_plus_at_mortar,
                  remote_normal_times_v_minus_at_mortar,
                  remote_gamma2_v_psi_at_mortar, remote_char_speeds_at_mortar);
              get(dt_psi_on_mortar)[linear_index] = get(dt_psi_at_mortar);
              get(dt_pi_on_mortar)[linear_index] = get(dt_pi_at_mortar);
              for (size_t i = 0; i < volume_dim; ++i) {
                dt_phi_on_mortar.get(i)[linear_index] = dt_phi_at_mortar.get(i);
              }
            });

        Variables<device_dt_boundary_tags> dt_boundary_on_face_batched{
            num_projection_work_items * num_face_points};
        project_from_mortar_device_batched(
            make_not_null(&dt_boundary_on_face_batched),
            dt_boundary_on_mortar_batched, num_projection_work_items, face_mesh,
            projection_group.mortar_mesh, projection_group.mortar_size);
        const auto dt_psi_on_face =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Psi>>>(
                dt_boundary_on_face_batched);
        const auto dt_pi_on_face =
            get<::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Pi>>>(
                dt_boundary_on_face_batched);
        const auto dt_phi_on_face = get<
            ::Tags::MirrorView<::Tags::dt<ScalarWave::Tags::Phi<volume_dim>>>>(
            dt_boundary_on_face_batched);
        ::Kokkos::parallel_for(
            "ComputeInternalBoundaryTermsBatchedAccumulateProjectionFaces",
            num_projection_work_items * num_face_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              const size_t face_index_on_face = linear_index % num_face_points;
              const size_t work_item_index = linear_index / num_face_points;
              const auto work_item = projection_work_items(work_item_index);
              const size_t local_linear_index =
                  work_item.local_element_index * num_face_points +
                  face_index_on_face;
              ::Kokkos::atomic_add(
                  &get(dt_psi_face_sum_all)[local_linear_index],
                  get(dt_psi_on_face)[linear_index]);
              ::Kokkos::atomic_add(&get(dt_pi_face_sum_all)[local_linear_index],
                                   get(dt_pi_on_face)[linear_index]);
              for (size_t i = 0; i < volume_dim; ++i) {
                ::Kokkos::atomic_add(
                    &dt_phi_face_sum_all.get(i)[local_linear_index],
                    dt_phi_on_face.get(i)[linear_index]);
              }
            });
      }
    }
  }
}

}  // namespace ScalarWave::Actions
