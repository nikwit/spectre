// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeInternalBoundaryTermsBatched.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
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
using packed_boundary_scratch_type =
    evolution::Kokkos::PackedBoundaryScratch<system>;
using device_package_field_tags =
    packed_boundary_scratch_type::device_package_field_tags;
using device_dt_boundary_tags =
    packed_boundary_scratch_type::device_dt_boundary_tags;
using projection_workspace_type =
    typename packed_boundary_scratch_type::projection_workspace_type;
using package_storage_type =
    typename Variables<device_package_field_tags>::storage_type;
using face_boundary_metadata = evolution::Kokkos::Batched::FaceBoundaryMetadata;
using projection_group_metadata =
    evolution::Kokkos::Batched::ProjectionGroupMetadata;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;

constexpr size_t num_spacetime_metric_components =
    (volume_dim + 1) * (volume_dim + 2) / 2;
constexpr size_t num_phi_components =
    volume_dim * num_spacetime_metric_components;
constexpr size_t num_char_speed_components = volume_dim + 1;
constexpr size_t offset_v_spacetime_metric = 0;
constexpr size_t offset_v_zero =
    offset_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t offset_v_plus = offset_v_zero + num_phi_components;
constexpr size_t offset_v_minus =
    offset_v_plus + num_spacetime_metric_components;
constexpr size_t offset_normal_times_v_plus =
    offset_v_minus + num_spacetime_metric_components;
constexpr size_t offset_normal_times_v_minus =
    offset_normal_times_v_plus + num_phi_components;
constexpr size_t offset_gamma2_v_spacetime_metric =
    offset_normal_times_v_minus + num_phi_components;
constexpr size_t offset_char_speeds =
    offset_gamma2_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t packaged_boundary_data_components =
    offset_char_speeds + num_char_speed_components;

constexpr size_t side_index(const Side side) {
  return side == Side::Upper ? static_cast<size_t>(1) : static_cast<size_t>(0);
}

constexpr size_t face_index(const size_t sliced_dim, const size_t side_i) {
  return 2 * sliced_dim + side_i;
}

size_t face_index(const Direction<volume_dim>& direction) {
  return face_index(direction.dimension(), side_index(direction.side()));
}

template <typename StorageType>
void apply_tensor_product_projection_2d_batched(
    const StorageType& result_view, const size_t result_num_points,
    const StorageType& input_view, const size_t input_num_points,
    const size_t num_work_items, const size_t source_points_dim_0,
    const size_t source_points_dim_1, const size_t target_points_dim_0,
    const size_t target_points_dim_1, const MatrixViewRO& matrix_dim_0,
    const MatrixViewRO& matrix_dim_1,
    const gsl::not_null<projection_workspace_type*> projected_dim_0_workspace) {
  const size_t num_components = input_view.extent(1);
  const size_t source_points_per_work_item =
      source_points_dim_0 * source_points_dim_1;
  const size_t target_points_per_work_item =
      target_points_dim_0 * target_points_dim_1;

  ASSERT(input_num_points == num_work_items * source_points_per_work_item,
         "Input size mismatch in batched 2D projection.");
  ASSERT(result_num_points == num_work_items * target_points_per_work_item,
         "Result size mismatch in batched 2D projection.");

  if (projected_dim_0_workspace->extent(0) != num_work_items or
      projected_dim_0_workspace->extent(1) !=
          target_points_dim_0 * source_points_dim_1 or
      projected_dim_0_workspace->extent(2) != num_components) {
    *projected_dim_0_workspace = projection_workspace_type(
        "GhComputeInternalBoundaryTermsBatchedProjectedDim0", num_work_items,
        target_points_dim_0 * source_points_dim_1, num_components);
  }
  const auto projected_dim_0 = *projected_dim_0_workspace;
  ::Kokkos::parallel_for(
      "GhComputeInternalBoundaryTermsBatchedProjectDim0",
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
      "GhComputeInternalBoundaryTermsBatchedProjectDim1",
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

void project_to_mortar_package_data_batched(
    const gsl::not_null<Variables<device_package_field_tags>*> result,
    const Variables<device_package_field_tags>& vars,
    const size_t num_work_items, const Mesh<2>& face_mesh,
    const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size,
    const gsl::not_null<projection_workspace_type*> projected_dim_0_workspace) {
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
      result->view(), result->number_of_grid_points(), vars.view(),
      vars.number_of_grid_points(), num_work_items, face_mesh.extents(0),
      face_mesh.extents(1), mortar_mesh.extents(0), mortar_mesh.extents(1),
      matrix_dim_0, matrix_dim_1, projected_dim_0_workspace);
}

void project_from_mortar_dt_data_batched(
    const gsl::not_null<Variables<device_dt_boundary_tags>*> result,
    const Variables<device_dt_boundary_tags>& vars, const size_t num_work_items,
    const Mesh<2>& face_mesh, const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size,
    const gsl::not_null<projection_workspace_type*> projected_dim_0_workspace) {
  ASSERT(Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size),
         "project_from_mortar_dt_data_batched should not be called when no "
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
      result->view(), result->number_of_grid_points(), vars.view(),
      vars.number_of_grid_points(), num_work_items, mortar_mesh.extents(0),
      mortar_mesh.extents(1), face_mesh.extents(0), face_mesh.extents(1),
      matrix_dim_0, matrix_dim_1, projected_dim_0_workspace);
}

KOKKOS_INLINE_FUNCTION double step_function_double(const double value) {
  return value < 0.0 ? 0.0 : 1.0;
}

KOKKOS_INLINE_FUNCTION void load_aa_from_packaged_data(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      tensor->get(a, b) =
          packaged_data_view(point, component_offset + component_index);
      ++component_index;
    }
  }
}

KOKKOS_INLINE_FUNCTION void load_iaa_from_packaged_data(
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  for (size_t d = 0; d < volume_dim; ++d) {
    for (size_t a = 0; a < volume_dim + 1; ++a) {
      for (size_t b = a; b < volume_dim + 1; ++b) {
        tensor->get(d, a, b) =
            packaged_data_view(point, component_offset + component_index);
        ++component_index;
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void load_a_from_packaged_data(
    const gsl::not_null<tnsr::a<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    tensor->get(a) = packaged_data_view(point, component_offset + a);
  }
}

KOKKOS_INLINE_FUNCTION void compute_boundary_terms_at_point(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*> dt_pi,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*> dt_phi,
    const tnsr::aa<double, volume_dim, Frame::Inertial>&
        local_v_spacetime_metric,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>& local_v_zero,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& local_v_plus,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& local_v_minus,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>&
        local_normal_times_v_plus,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>&
        local_normal_times_v_minus,
    const tnsr::aa<double, volume_dim, Frame::Inertial>&
        local_gamma2_v_spacetime_metric,
    const tnsr::a<double, volume_dim, Frame::Inertial>& local_char_speeds,
    const tnsr::aa<double, volume_dim, Frame::Inertial>&
        remote_v_spacetime_metric,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>& remote_v_zero,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& remote_v_plus,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& remote_v_minus,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>&
        remote_normal_times_v_plus,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>&
        remote_normal_times_v_minus,
    const tnsr::aa<double, volume_dim, Frame::Inertial>&
        remote_gamma2_v_spacetime_metric,
    const tnsr::a<double, volume_dim, Frame::Inertial>& remote_char_speeds) {
  const double weighted_lambda_spacetime_metric_int =
      step_function_double(-local_char_speeds.get(0));
  const double weighted_lambda_spacetime_metric_ext =
      -step_function_double(remote_char_speeds.get(0));
  const double weighted_lambda_zero_int =
      step_function_double(-local_char_speeds.get(1));
  const double weighted_lambda_zero_ext =
      -step_function_double(remote_char_speeds.get(1));
  const double weighted_lambda_plus_int =
      step_function_double(-local_char_speeds.get(2));
  const double weighted_lambda_plus_ext =
      -step_function_double(remote_char_speeds.get(2));
  const double weighted_lambda_minus_int =
      step_function_double(-local_char_speeds.get(3));
  const double weighted_lambda_minus_ext =
      -step_function_double(remote_char_speeds.get(3));

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      dt_spacetime_metric->get(a, b) = weighted_lambda_spacetime_metric_ext *
                                           remote_v_spacetime_metric.get(a, b) -
                                       weighted_lambda_spacetime_metric_int *
                                           local_v_spacetime_metric.get(a, b);

      dt_pi->get(a, b) =
          0.5 * (weighted_lambda_plus_ext * remote_v_plus.get(a, b) +
                 weighted_lambda_minus_ext * remote_v_minus.get(a, b)) +
          weighted_lambda_spacetime_metric_ext *
              remote_gamma2_v_spacetime_metric.get(a, b) -
          0.5 * (weighted_lambda_plus_int * local_v_plus.get(a, b) +
                 weighted_lambda_minus_int * local_v_minus.get(a, b)) -
          weighted_lambda_spacetime_metric_int *
              local_gamma2_v_spacetime_metric.get(a, b);

      for (size_t d = 0; d < volume_dim; ++d) {
        dt_phi->get(d, a, b) =
            -0.5 * (weighted_lambda_minus_ext *
                        remote_normal_times_v_minus.get(d, a, b) -
                    weighted_lambda_plus_ext *
                        remote_normal_times_v_plus.get(d, a, b)) +
            weighted_lambda_zero_ext * remote_v_zero.get(d, a, b) -
            0.5 * (weighted_lambda_plus_int *
                       local_normal_times_v_plus.get(d, a, b) -
                   weighted_lambda_minus_int *
                       local_normal_times_v_minus.get(d, a, b)) -
            weighted_lambda_zero_int * local_v_zero.get(d, a, b);
      }
    }
  }
}

}  // namespace

void ComputeInternalBoundaryTermsBatched::apply(
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_boundary_metadata_tag::type&
        packed_boundary_metadata) {
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
  auto& local_packaged_face_data_batched =
      packed_boundary_scratch->projection_local_packaged_face_data;
  auto& remote_packaged_face_data_batched =
      packed_boundary_scratch->projection_remote_packaged_face_data;
  auto& local_packaged_mortar_data_batched =
      packed_boundary_scratch->projection_local_packaged_mortar_data;
  auto& remote_packaged_mortar_data_batched =
      packed_boundary_scratch->projection_remote_packaged_mortar_data;
  auto& dt_boundary_on_mortar_batched =
      packed_boundary_scratch->projection_dt_boundary_on_mortar;
  auto& dt_boundary_on_face_batched =
      packed_boundary_scratch->projection_dt_boundary_on_face;
  auto& projection_workspace = packed_boundary_scratch->projection_workspace;

  for (size_t d = 0; d < volume_dim; ++d) {
    const Mesh<volume_dim - 1>& face_mesh = uniform_face_mesh;
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const size_t local_face_id = face_index(d, side_i);
      auto& dt_face_sum_all = dt_face_sum_for_all_elements[local_face_id];
      const auto dt_face_sum_all_view = dt_face_sum_all.view();
      const size_t dt_face_points = dt_face_sum_all_view.extent(0);
      const size_t dt_face_components = dt_face_sum_all_view.extent(1);
      ::Kokkos::parallel_for(
          "GhComputeInternalBoundaryTermsBatchedZeroFaceSum",
          dt_face_points * dt_face_components,
          KOKKOS_LAMBDA(const size_t linear_index) {
            const size_t point = linear_index / dt_face_components;
            const size_t component = linear_index % dt_face_components;
            dt_face_sum_all_view(point, component) = 0.0;
          });

      const face_boundary_metadata& cached_face_metadata =
          packed_boundary_metadata
              .boundary_correction_face_metadata[local_face_id];
      const auto local_packaged_face_data_all_view =
          packaged_face_data_for_all_elements[local_face_id].view();

      const auto dt_spacetime_metric_face_sum_all =
          get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
              dt_face_sum_all);
      const auto dt_pi_face_sum_all = get<::Tags::MirrorView<
          ::Tags::dt<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>>>(
          dt_face_sum_all);
      const auto dt_phi_face_sum_all = get<::Tags::MirrorView<
          ::Tags::dt<gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>>>(
          dt_face_sum_all);

      const auto& mortar_work_items =
          cached_face_metadata.no_projection_work_items;
      const auto& oriented_remote_face_indices =
          cached_face_metadata.no_projection_oriented_remote_face_indices;
      if (mortar_work_items.extent(0) > 0) {
        ::Kokkos::parallel_for(
            "GhComputeInternalBoundaryTermsBatchedNoProjectionMortars",
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

              tnsr::aa<double, volume_dim, Frame::Inertial>
                  local_v_spacetime_metric{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> local_v_zero{};
              tnsr::aa<double, volume_dim, Frame::Inertial> local_v_plus{};
              tnsr::aa<double, volume_dim, Frame::Inertial> local_v_minus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus{};
              tnsr::aa<double, volume_dim, Frame::Inertial>
                  local_gamma2_v_spacetime_metric{};
              tnsr::a<double, volume_dim, Frame::Inertial> local_char_speeds{};

              load_aa_from_packaged_data(
                  make_not_null(&local_v_spacetime_metric),
                  local_packaged_face_data_all_view, local_linear_index,
                  offset_v_spacetime_metric);
              load_iaa_from_packaged_data(make_not_null(&local_v_zero),
                                          local_packaged_face_data_all_view,
                                          local_linear_index, offset_v_zero);
              load_aa_from_packaged_data(make_not_null(&local_v_plus),
                                         local_packaged_face_data_all_view,
                                         local_linear_index, offset_v_plus);
              load_aa_from_packaged_data(make_not_null(&local_v_minus),
                                         local_packaged_face_data_all_view,
                                         local_linear_index, offset_v_minus);
              load_iaa_from_packaged_data(
                  make_not_null(&local_normal_times_v_plus),
                  local_packaged_face_data_all_view, local_linear_index,
                  offset_normal_times_v_plus);
              load_iaa_from_packaged_data(
                  make_not_null(&local_normal_times_v_minus),
                  local_packaged_face_data_all_view, local_linear_index,
                  offset_normal_times_v_minus);
              load_aa_from_packaged_data(
                  make_not_null(&local_gamma2_v_spacetime_metric),
                  local_packaged_face_data_all_view, local_linear_index,
                  offset_gamma2_v_spacetime_metric);
              load_a_from_packaged_data(make_not_null(&local_char_speeds),
                                        local_packaged_face_data_all_view,
                                        local_linear_index, offset_char_speeds);

              const size_t remote_face_index = oriented_remote_face_indices(
                  mortar_work_item.oriented_remote_face_index_offset +
                  face_index_on_face);
              const size_t remote_linear_index =
                  mortar_work_item.remote_element_index * num_face_points +
                  remote_face_index;

              auto remote_packaged_face_data_all_view =
                  packaged_face_data_view_2_1;
              if (mortar_work_item.remote_face_id == 0) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_0_0;
              } else if (mortar_work_item.remote_face_id == 1) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_0_1;
              } else if (mortar_work_item.remote_face_id == 2) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_1_0;
              } else if (mortar_work_item.remote_face_id == 3) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_1_1;
              } else if (mortar_work_item.remote_face_id == 4) {
                remote_packaged_face_data_all_view =
                    packaged_face_data_view_2_0;
              }

              tnsr::aa<double, volume_dim, Frame::Inertial>
                  remote_v_spacetime_metric{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> remote_v_zero{};
              tnsr::aa<double, volume_dim, Frame::Inertial> remote_v_plus{};
              tnsr::aa<double, volume_dim, Frame::Inertial> remote_v_minus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_plus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_minus{};
              tnsr::aa<double, volume_dim, Frame::Inertial>
                  remote_gamma2_v_spacetime_metric{};
              tnsr::a<double, volume_dim, Frame::Inertial> remote_char_speeds{};
              load_aa_from_packaged_data(
                  make_not_null(&remote_v_spacetime_metric),
                  remote_packaged_face_data_all_view, remote_linear_index,
                  offset_v_spacetime_metric);
              load_iaa_from_packaged_data(make_not_null(&remote_v_zero),
                                          remote_packaged_face_data_all_view,
                                          remote_linear_index, offset_v_zero);
              load_aa_from_packaged_data(make_not_null(&remote_v_plus),
                                         remote_packaged_face_data_all_view,
                                         remote_linear_index, offset_v_plus);
              load_aa_from_packaged_data(make_not_null(&remote_v_minus),
                                         remote_packaged_face_data_all_view,
                                         remote_linear_index, offset_v_minus);
              load_iaa_from_packaged_data(
                  make_not_null(&remote_normal_times_v_plus),
                  remote_packaged_face_data_all_view, remote_linear_index,
                  offset_normal_times_v_plus);
              load_iaa_from_packaged_data(
                  make_not_null(&remote_normal_times_v_minus),
                  remote_packaged_face_data_all_view, remote_linear_index,
                  offset_normal_times_v_minus);
              load_aa_from_packaged_data(
                  make_not_null(&remote_gamma2_v_spacetime_metric),
                  remote_packaged_face_data_all_view, remote_linear_index,
                  offset_gamma2_v_spacetime_metric);
              load_a_from_packaged_data(make_not_null(&remote_char_speeds),
                                        remote_packaged_face_data_all_view,
                                        remote_linear_index,
                                        offset_char_speeds);

              tnsr::aa<double, volume_dim, Frame::Inertial>
                  dt_spacetime_metric_at_face{};
              tnsr::aa<double, volume_dim, Frame::Inertial> dt_pi_at_face{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> dt_phi_at_face{};
              compute_boundary_terms_at_point(
                  make_not_null(&dt_spacetime_metric_at_face),
                  make_not_null(&dt_pi_at_face), make_not_null(&dt_phi_at_face),
                  local_v_spacetime_metric, local_v_zero, local_v_plus,
                  local_v_minus, local_normal_times_v_plus,
                  local_normal_times_v_minus, local_gamma2_v_spacetime_metric,
                  local_char_speeds, remote_v_spacetime_metric, remote_v_zero,
                  remote_v_plus, remote_v_minus, remote_normal_times_v_plus,
                  remote_normal_times_v_minus, remote_gamma2_v_spacetime_metric,
                  remote_char_speeds);

              for (size_t a = 0; a < volume_dim + 1; ++a) {
                for (size_t b = a; b < volume_dim + 1; ++b) {
                  ::Kokkos::atomic_add(&dt_spacetime_metric_face_sum_all.get(
                                           a, b)[local_linear_index],
                                       dt_spacetime_metric_at_face.get(a, b));
                  ::Kokkos::atomic_add(
                      &dt_pi_face_sum_all.get(a, b)[local_linear_index],
                      dt_pi_at_face.get(a, b));
                  for (size_t i = 0; i < volume_dim; ++i) {
                    ::Kokkos::atomic_add(
                        &dt_phi_face_sum_all.get(i, a, b)[local_linear_index],
                        dt_phi_at_face.get(i, a, b));
                  }
                }
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

        const size_t projection_face_points =
            num_projection_work_items * num_face_points;
        if (local_packaged_face_data_batched.number_of_grid_points() !=
            projection_face_points) {
          local_packaged_face_data_batched.initialize(projection_face_points);
        }
        if (remote_packaged_face_data_batched.number_of_grid_points() !=
            projection_face_points) {
          remote_packaged_face_data_batched.initialize(projection_face_points);
        }
        auto local_packaged_face_data_batched_view =
            local_packaged_face_data_batched.view();
        auto remote_packaged_face_data_batched_view =
            remote_packaged_face_data_batched.view();
        const size_t num_package_components =
            local_packaged_face_data_batched_view.extent(1);

        ::Kokkos::parallel_for(
            "GhComputeInternalBoundaryTermsBatchedGatherProjectionLocalFaces",
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
            "GhComputeInternalBoundaryTermsBatchedGatherProjectionRemoteFaces",
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

        const size_t projection_mortar_points =
            num_projection_work_items * num_mortar_points;
        if (local_packaged_mortar_data_batched.number_of_grid_points() !=
            projection_mortar_points) {
          local_packaged_mortar_data_batched.initialize(
              projection_mortar_points);
        }
        if (remote_packaged_mortar_data_batched.number_of_grid_points() !=
            projection_mortar_points) {
          remote_packaged_mortar_data_batched.initialize(
              projection_mortar_points);
        }
        project_to_mortar_package_data_batched(
            make_not_null(&local_packaged_mortar_data_batched),
            local_packaged_face_data_batched, num_projection_work_items,
            face_mesh, projection_group.mortar_mesh,
            projection_group.mortar_size, make_not_null(&projection_workspace));
        project_to_mortar_package_data_batched(
            make_not_null(&remote_packaged_mortar_data_batched),
            remote_packaged_face_data_batched, num_projection_work_items,
            face_mesh, projection_group.mortar_mesh,
            projection_group.mortar_size, make_not_null(&projection_workspace));

        if (dt_boundary_on_mortar_batched.number_of_grid_points() !=
            projection_mortar_points) {
          dt_boundary_on_mortar_batched.initialize(projection_mortar_points);
        }
        const auto dt_spacetime_metric_on_mortar =
            get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
                dt_boundary_on_mortar_batched);
        const auto dt_pi_on_mortar = get<::Tags::MirrorView<
            ::Tags::dt<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>>>(
            dt_boundary_on_mortar_batched);
        const auto dt_phi_on_mortar = get<::Tags::MirrorView<::Tags::dt<
            gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>>>(
            dt_boundary_on_mortar_batched);
        const auto local_packaged_mortar_data_view =
            local_packaged_mortar_data_batched.view();
        const auto remote_packaged_mortar_data_view =
            remote_packaged_mortar_data_batched.view();
        ::Kokkos::parallel_for(
            "GhComputeInternalBoundaryTermsBatchedProjectionMortarBoundaryTerm"
            "s",
            num_projection_work_items * num_mortar_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              tnsr::aa<double, volume_dim, Frame::Inertial>
                  local_v_spacetime_metric{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> local_v_zero{};
              tnsr::aa<double, volume_dim, Frame::Inertial> local_v_plus{};
              tnsr::aa<double, volume_dim, Frame::Inertial> local_v_minus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_plus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  local_normal_times_v_minus{};
              tnsr::aa<double, volume_dim, Frame::Inertial>
                  local_gamma2_v_spacetime_metric{};
              tnsr::a<double, volume_dim, Frame::Inertial> local_char_speeds{};
              load_aa_from_packaged_data(
                  make_not_null(&local_v_spacetime_metric),
                  local_packaged_mortar_data_view, linear_index,
                  offset_v_spacetime_metric);
              load_iaa_from_packaged_data(make_not_null(&local_v_zero),
                                          local_packaged_mortar_data_view,
                                          linear_index, offset_v_zero);
              load_aa_from_packaged_data(make_not_null(&local_v_plus),
                                         local_packaged_mortar_data_view,
                                         linear_index, offset_v_plus);
              load_aa_from_packaged_data(make_not_null(&local_v_minus),
                                         local_packaged_mortar_data_view,
                                         linear_index, offset_v_minus);
              load_iaa_from_packaged_data(
                  make_not_null(&local_normal_times_v_plus),
                  local_packaged_mortar_data_view, linear_index,
                  offset_normal_times_v_plus);
              load_iaa_from_packaged_data(
                  make_not_null(&local_normal_times_v_minus),
                  local_packaged_mortar_data_view, linear_index,
                  offset_normal_times_v_minus);
              load_aa_from_packaged_data(
                  make_not_null(&local_gamma2_v_spacetime_metric),
                  local_packaged_mortar_data_view, linear_index,
                  offset_gamma2_v_spacetime_metric);
              load_a_from_packaged_data(make_not_null(&local_char_speeds),
                                        local_packaged_mortar_data_view,
                                        linear_index, offset_char_speeds);

              tnsr::aa<double, volume_dim, Frame::Inertial>
                  remote_v_spacetime_metric{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> remote_v_zero{};
              tnsr::aa<double, volume_dim, Frame::Inertial> remote_v_plus{};
              tnsr::aa<double, volume_dim, Frame::Inertial> remote_v_minus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_plus{};
              tnsr::iaa<double, volume_dim, Frame::Inertial>
                  remote_normal_times_v_minus{};
              tnsr::aa<double, volume_dim, Frame::Inertial>
                  remote_gamma2_v_spacetime_metric{};
              tnsr::a<double, volume_dim, Frame::Inertial> remote_char_speeds{};
              load_aa_from_packaged_data(
                  make_not_null(&remote_v_spacetime_metric),
                  remote_packaged_mortar_data_view, linear_index,
                  offset_v_spacetime_metric);
              load_iaa_from_packaged_data(make_not_null(&remote_v_zero),
                                          remote_packaged_mortar_data_view,
                                          linear_index, offset_v_zero);
              load_aa_from_packaged_data(make_not_null(&remote_v_plus),
                                         remote_packaged_mortar_data_view,
                                         linear_index, offset_v_plus);
              load_aa_from_packaged_data(make_not_null(&remote_v_minus),
                                         remote_packaged_mortar_data_view,
                                         linear_index, offset_v_minus);
              load_iaa_from_packaged_data(
                  make_not_null(&remote_normal_times_v_plus),
                  remote_packaged_mortar_data_view, linear_index,
                  offset_normal_times_v_plus);
              load_iaa_from_packaged_data(
                  make_not_null(&remote_normal_times_v_minus),
                  remote_packaged_mortar_data_view, linear_index,
                  offset_normal_times_v_minus);
              load_aa_from_packaged_data(
                  make_not_null(&remote_gamma2_v_spacetime_metric),
                  remote_packaged_mortar_data_view, linear_index,
                  offset_gamma2_v_spacetime_metric);
              load_a_from_packaged_data(make_not_null(&remote_char_speeds),
                                        remote_packaged_mortar_data_view,
                                        linear_index, offset_char_speeds);

              tnsr::aa<double, volume_dim, Frame::Inertial>
                  dt_spacetime_metric_at_mortar{};
              tnsr::aa<double, volume_dim, Frame::Inertial> dt_pi_at_mortar{};
              tnsr::iaa<double, volume_dim, Frame::Inertial> dt_phi_at_mortar{};
              compute_boundary_terms_at_point(
                  make_not_null(&dt_spacetime_metric_at_mortar),
                  make_not_null(&dt_pi_at_mortar),
                  make_not_null(&dt_phi_at_mortar), local_v_spacetime_metric,
                  local_v_zero, local_v_plus, local_v_minus,
                  local_normal_times_v_plus, local_normal_times_v_minus,
                  local_gamma2_v_spacetime_metric, local_char_speeds,
                  remote_v_spacetime_metric, remote_v_zero, remote_v_plus,
                  remote_v_minus, remote_normal_times_v_plus,
                  remote_normal_times_v_minus, remote_gamma2_v_spacetime_metric,
                  remote_char_speeds);
              for (size_t a = 0; a < volume_dim + 1; ++a) {
                for (size_t b = a; b < volume_dim + 1; ++b) {
                  dt_spacetime_metric_on_mortar.get(a, b)[linear_index] =
                      dt_spacetime_metric_at_mortar.get(a, b);
                  dt_pi_on_mortar.get(a, b)[linear_index] =
                      dt_pi_at_mortar.get(a, b);
                  for (size_t i = 0; i < volume_dim; ++i) {
                    dt_phi_on_mortar.get(i, a, b)[linear_index] =
                        dt_phi_at_mortar.get(i, a, b);
                  }
                }
              }
            });

        if (dt_boundary_on_face_batched.number_of_grid_points() !=
            projection_face_points) {
          dt_boundary_on_face_batched.initialize(projection_face_points);
        }
        project_from_mortar_dt_data_batched(
            make_not_null(&dt_boundary_on_face_batched),
            dt_boundary_on_mortar_batched, num_projection_work_items, face_mesh,
            projection_group.mortar_mesh, projection_group.mortar_size,
            make_not_null(&projection_workspace));
        const auto dt_spacetime_metric_on_face =
            get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
                dt_boundary_on_face_batched);
        const auto dt_pi_on_face = get<::Tags::MirrorView<
            ::Tags::dt<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>>>(
            dt_boundary_on_face_batched);
        const auto dt_phi_on_face = get<::Tags::MirrorView<::Tags::dt<
            gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>>>(
            dt_boundary_on_face_batched);

        ::Kokkos::parallel_for(
            "GhComputeInternalBoundaryTermsBatchedAccumulateProjectionFaces",
            num_projection_work_items * num_face_points,
            KOKKOS_LAMBDA(const int linear_index_int) {
              const size_t linear_index = static_cast<size_t>(linear_index_int);
              const size_t face_index_on_face = linear_index % num_face_points;
              const size_t work_item_index = linear_index / num_face_points;
              const auto work_item = projection_work_items(work_item_index);
              const size_t local_linear_index =
                  work_item.local_element_index * num_face_points +
                  face_index_on_face;
              for (size_t a = 0; a < volume_dim + 1; ++a) {
                for (size_t b = a; b < volume_dim + 1; ++b) {
                  ::Kokkos::atomic_add(
                      &dt_spacetime_metric_face_sum_all.get(
                          a, b)[local_linear_index],
                      dt_spacetime_metric_on_face.get(a, b)[linear_index]);
                  ::Kokkos::atomic_add(
                      &dt_pi_face_sum_all.get(a, b)[local_linear_index],
                      dt_pi_on_face.get(a, b)[linear_index]);
                  for (size_t i = 0; i < volume_dim; ++i) {
                    ::Kokkos::atomic_add(
                        &dt_phi_face_sum_all.get(i, a, b)[local_linear_index],
                        dt_phi_on_face.get(i, a, b)[linear_index]);
                  }
                }
              }
            });
      }
    }
  }
}

}  // namespace gh::Actions
