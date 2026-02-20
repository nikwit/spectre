// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

constexpr size_t volume_dim = 3;
using system = gh::System<volume_dim>;
using boundary_correction_data_type =
    gh::KokkosTags::BoundaryCorrectionData<volume_dim>;
using device_variables_type = gh::KokkosTags::DeviceVariables<system>::type;
using device_dt_type = gh::KokkosTags::DeviceDtVariables<system>::type;
using device_face_to_volume_index_map_type =
    gh::KokkosTags::DeviceFaceToVolumeIndexMap<volume_dim>::type;
using device_face_unit_normal_covector_type =
    gh::KokkosTags::DeviceFaceUnitNormalCovector<volume_dim>::type;
using device_face_normal_magnitude_type =
    gh::KokkosTags::DeviceFaceNormalMagnitude<volume_dim>::type;
using boundary_data_buffer_type =
    gh::KokkosTags::BoundaryCorrectionDataBuffer<volume_dim>;
using boundary_data_storage_type = typename boundary_data_buffer_type::storage_type;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;
using dt_boundary_tags =
    tmpl::list<::Tags::dt<spacetime_metric_tag>,
               ::Tags::dt<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>,
               ::Tags::dt<gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>>;
using device_dt_boundary_tags =
    db::wrap_tags_in<::Tags::MirrorView, dt_boundary_tags>;
using dt_boundary_storage_type =
    typename Variables<device_dt_boundary_tags>::storage_type;

constexpr size_t num_spacetime_metric_components =
    (volume_dim + 1) * (volume_dim + 2) / 2;
constexpr size_t num_phi_components =
    volume_dim * num_spacetime_metric_components;
constexpr size_t num_char_speed_components = volume_dim + 1;
constexpr size_t offset_v_spacetime_metric = 0;
constexpr size_t offset_v_zero =
    offset_v_spacetime_metric + num_spacetime_metric_components;
constexpr size_t offset_v_plus = offset_v_zero + num_phi_components;
constexpr size_t offset_v_minus = offset_v_plus + num_spacetime_metric_components;
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

namespace detail {

void apply_tensor_product_projection_2d(
    const dt_boundary_storage_type& result_view, const size_t result_num_points,
    const dt_boundary_storage_type& input_view, const size_t input_num_points,
    const MatrixViewRO& matrix_dim_0, const MatrixViewRO& matrix_dim_1) {
  const size_t source_points_dim_0 = matrix_dim_0.extent(1);
  const size_t source_points_dim_1 = matrix_dim_1.extent(1);
  const size_t target_points_dim_0 = matrix_dim_0.extent(0);
  const size_t target_points_dim_1 = matrix_dim_1.extent(0);
  const size_t num_components = input_view.extent(1);

  ASSERT(input_num_points == source_points_dim_0 * source_points_dim_1,
         "Input has " << input_num_points << " points, expected "
                      << source_points_dim_0 << " * " << source_points_dim_1
                      << ".");
  ASSERT(result_num_points == target_points_dim_0 * target_points_dim_1,
         "Result has " << result_num_points
                       << " points, expected " << target_points_dim_0 << " * "
                       << target_points_dim_1 << ".");

  ::Kokkos::View<double**> projected_dim_0(
      "GhKokkosProjectFromMortarDim0",
      target_points_dim_0 * source_points_dim_1, num_components);
  ::Kokkos::parallel_for(
      "GhKokkosProjectFromMortarDim0",
      ::Kokkos::RangePolicy<size_t>{
          0, target_points_dim_0 * source_points_dim_1 * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i0 = remaining % target_points_dim_0;
        const size_t i1 = remaining / target_points_dim_0;

        double sum = 0.0;
        for (size_t k0 = 0; k0 < source_points_dim_0; ++k0) {
          sum += matrix_dim_0(i0, k0) *
                 input_view(k0 + source_points_dim_0 * i1, component);
        }
        projected_dim_0(i0 + target_points_dim_0 * i1, component) = sum;
      });

  ::Kokkos::parallel_for(
      "GhKokkosProjectFromMortarDim1",
      ::Kokkos::RangePolicy<size_t>{
          0, target_points_dim_0 * target_points_dim_1 * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i0 = remaining % target_points_dim_0;
        const size_t i1 = remaining / target_points_dim_0;

        double sum = 0.0;
        for (size_t k1 = 0; k1 < source_points_dim_1; ++k1) {
          sum += matrix_dim_1(i1, k1) *
                 projected_dim_0(i0 + target_points_dim_0 * k1, component);
        }
        result_view(i0 + target_points_dim_0 * i1, component) = sum;
      });
}

void project_from_mortar_device(
    const gsl::not_null<Variables<device_dt_boundary_tags>*> result,
    const Variables<device_dt_boundary_tags>& vars, const Mesh<2>& face_mesh,
    const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size) {
  ASSERT(Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size),
         "project_from_mortar_device should not be called when no projection "
         "is needed.");
  ASSERT(result->number_of_grid_points() == face_mesh.number_of_grid_points(),
         "Result has " << result->number_of_grid_points()
                       << " points, expected "
                       << face_mesh.number_of_grid_points() << ".");
  const auto& matrix_dim_0 =
      Spectral::projection_matrix_child_to_parent_on_device(
          mortar_mesh.slice_through(0), face_mesh.slice_through(0),
          gsl::at(mortar_size, 0));
  const auto& matrix_dim_1 =
      Spectral::projection_matrix_child_to_parent_on_device(
          mortar_mesh.slice_through(1), face_mesh.slice_through(1),
          gsl::at(mortar_size, 1));
  apply_tensor_product_projection_2d(
      result->view(), result->number_of_grid_points(), vars.view(),
      vars.number_of_grid_points(), matrix_dim_0, matrix_dim_1);
}

}  // namespace detail

KOKKOS_INLINE_FUNCTION double step_function_double(const double value) {
  return value < 0.0 ? 0.0 : 1.0;
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

KOKKOS_INLINE_FUNCTION void load_aa_from_boundary_data(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        tensor,
    const boundary_data_storage_type& boundary_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      tensor->get(a, b) =
          boundary_data_view(point, component_offset + component_index);
      ++component_index;
    }
  }
}

KOKKOS_INLINE_FUNCTION void load_iaa_from_boundary_data(
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        tensor,
    const boundary_data_storage_type& boundary_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  for (size_t d = 0; d < volume_dim; ++d) {
    for (size_t a = 0; a < volume_dim + 1; ++a) {
      for (size_t b = a; b < volume_dim + 1; ++b) {
        tensor->get(d, a, b) =
            boundary_data_view(point, component_offset + component_index);
        ++component_index;
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void load_a_from_boundary_data(
    const gsl::not_null<tnsr::a<double, volume_dim, Frame::Inertial>*> tensor,
    const boundary_data_storage_type& boundary_data_view, const size_t point,
    const size_t component_offset) {
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    tensor->get(a) = boundary_data_view(point, component_offset + a);
  }
}

void accumulate_mortar_pair_to_face_sum(
    const gsl::not_null<Variables<device_dt_boundary_tags>*>
        dt_boundary_correction_on_face_sum,
    const boundary_correction_data_type& local_boundary_data,
    const boundary_correction_data_type& remote_boundary_data,
    const Mesh<2>& face_mesh, const Mesh<2>& mortar_mesh,
    const bool needs_projection,
    const std::array<Spectral::SegmentSize, volume_dim - 1>& mortar_size) {
  ASSERT(local_boundary_data.boundary_correction_data.number_of_grid_points() ==
             mortar_mesh.number_of_grid_points(),
         "Local packaged data size does not match mortar mesh.");
  ASSERT(
      remote_boundary_data.boundary_correction_data.number_of_grid_points() ==
          mortar_mesh.number_of_grid_points(),
      "Remote packaged data size does not match mortar mesh.");
  const auto local_boundary_data_view =
      local_boundary_data.boundary_correction_data.view();
  const auto remote_boundary_data_view =
      remote_boundary_data.boundary_correction_data.view();
  ASSERT(
      local_boundary_data_view.extent(1) == packaged_boundary_data_components,
      "Local packaged data has " << local_boundary_data_view.extent(1)
                                 << " components per point, expected "
                                 << packaged_boundary_data_components << ".");
  ASSERT(
      remote_boundary_data_view.extent(1) == packaged_boundary_data_components,
      "Remote packaged data has " << remote_boundary_data_view.extent(1)
                                  << " components per point, expected "
                                  << packaged_boundary_data_components << ".");

  Variables<device_dt_boundary_tags> dt_boundary_correction_on_mortar{
      mortar_mesh.number_of_grid_points()};
  const auto dt_spacetime_metric_on_mortar =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
          dt_boundary_correction_on_mortar);
  const auto dt_pi_on_mortar =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Pi<DataVector, volume_dim,
                                                     Frame::Inertial>>>>(
          dt_boundary_correction_on_mortar);
  const auto dt_phi_on_mortar =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Phi<DataVector, volume_dim,
                                                      Frame::Inertial>>>>(
          dt_boundary_correction_on_mortar);

  ::Kokkos::parallel_for(
      "GhABCTDComputeBoundaryCorrectionOnMortar",
      mortar_mesh.number_of_grid_points(),
      KOKKOS_LAMBDA(const int mortar_index_int) {
        const size_t mortar_index = static_cast<size_t>(mortar_index_int);

        tnsr::aa<double, volume_dim, Frame::Inertial> local_v_spacetime_metric{};
        tnsr::iaa<double, volume_dim, Frame::Inertial> local_v_zero{};
        tnsr::aa<double, volume_dim, Frame::Inertial> local_v_plus{};
        tnsr::aa<double, volume_dim, Frame::Inertial> local_v_minus{};
        tnsr::iaa<double, volume_dim, Frame::Inertial> local_normal_times_v_plus{};
        tnsr::iaa<double, volume_dim, Frame::Inertial> local_normal_times_v_minus{};
        tnsr::aa<double, volume_dim, Frame::Inertial>
            local_gamma2_v_spacetime_metric{};
        tnsr::a<double, volume_dim, Frame::Inertial> local_char_speeds{};
        tnsr::aa<double, volume_dim, Frame::Inertial> remote_v_spacetime_metric{};
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

        load_aa_from_boundary_data(make_not_null(&local_v_spacetime_metric),
                                   local_boundary_data_view, mortar_index,
                                   offset_v_spacetime_metric);
        load_iaa_from_boundary_data(make_not_null(&local_v_zero),
                                    local_boundary_data_view, mortar_index,
                                    offset_v_zero);
        load_aa_from_boundary_data(make_not_null(&local_v_plus),
                                   local_boundary_data_view, mortar_index,
                                   offset_v_plus);
        load_aa_from_boundary_data(make_not_null(&local_v_minus),
                                   local_boundary_data_view, mortar_index,
                                   offset_v_minus);
        load_iaa_from_boundary_data(make_not_null(&local_normal_times_v_plus),
                                    local_boundary_data_view, mortar_index,
                                    offset_normal_times_v_plus);
        load_iaa_from_boundary_data(make_not_null(&local_normal_times_v_minus),
                                    local_boundary_data_view, mortar_index,
                                    offset_normal_times_v_minus);
        load_aa_from_boundary_data(
            make_not_null(&local_gamma2_v_spacetime_metric),
            local_boundary_data_view, mortar_index,
            offset_gamma2_v_spacetime_metric);
        load_a_from_boundary_data(make_not_null(&local_char_speeds),
                                  local_boundary_data_view, mortar_index,
                                  offset_char_speeds);

        load_aa_from_boundary_data(make_not_null(&remote_v_spacetime_metric),
                                   remote_boundary_data_view, mortar_index,
                                   offset_v_spacetime_metric);
        load_iaa_from_boundary_data(make_not_null(&remote_v_zero),
                                    remote_boundary_data_view, mortar_index,
                                    offset_v_zero);
        load_aa_from_boundary_data(make_not_null(&remote_v_plus),
                                   remote_boundary_data_view, mortar_index,
                                   offset_v_plus);
        load_aa_from_boundary_data(make_not_null(&remote_v_minus),
                                   remote_boundary_data_view, mortar_index,
                                   offset_v_minus);
        load_iaa_from_boundary_data(make_not_null(&remote_normal_times_v_plus),
                                    remote_boundary_data_view, mortar_index,
                                    offset_normal_times_v_plus);
        load_iaa_from_boundary_data(
            make_not_null(&remote_normal_times_v_minus),
            remote_boundary_data_view, mortar_index,
            offset_normal_times_v_minus);
        load_aa_from_boundary_data(
            make_not_null(&remote_gamma2_v_spacetime_metric),
            remote_boundary_data_view, mortar_index,
            offset_gamma2_v_spacetime_metric);
        load_a_from_boundary_data(make_not_null(&remote_char_speeds),
                                  remote_boundary_data_view, mortar_index,
                                  offset_char_speeds);

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
            dt_spacetime_metric_on_mortar.get(a, b)[mortar_index] =
                weighted_lambda_spacetime_metric_ext *
                    remote_v_spacetime_metric.get(a, b) -
                weighted_lambda_spacetime_metric_int *
                    local_v_spacetime_metric.get(a, b);

            dt_pi_on_mortar.get(a, b)[mortar_index] =
                0.5 * (weighted_lambda_plus_ext *
                           remote_v_plus.get(a, b) +
                       weighted_lambda_minus_ext *
                           remote_v_minus.get(a, b)) +
                weighted_lambda_spacetime_metric_ext *
                    remote_gamma2_v_spacetime_metric.get(a, b) -
                0.5 * (weighted_lambda_plus_int *
                           local_v_plus.get(a, b) +
                       weighted_lambda_minus_int *
                           local_v_minus.get(a, b)) -
                weighted_lambda_spacetime_metric_int *
                    local_gamma2_v_spacetime_metric.get(a, b);

            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi_on_mortar.get(d, a, b)[mortar_index] =
                  -0.5 * (weighted_lambda_minus_ext *
                              remote_normal_times_v_minus.get(d, a, b) -
                          weighted_lambda_plus_ext *
                              remote_normal_times_v_plus.get(d, a, b)) +
                  weighted_lambda_zero_ext *
                      remote_v_zero.get(d, a, b) -
                  0.5 * (weighted_lambda_plus_int *
                             local_normal_times_v_plus.get(d, a, b) -
                         weighted_lambda_minus_int *
                             local_normal_times_v_minus.get(d, a, b)) -
                  weighted_lambda_zero_int *
                      local_v_zero.get(d, a, b);
            }
          }
        }
      });

  Variables<device_dt_boundary_tags> dt_boundary_correction_on_face{
      face_mesh.number_of_grid_points()};
  if (needs_projection) {
    detail::project_from_mortar_device(
        make_not_null(&dt_boundary_correction_on_face),
        dt_boundary_correction_on_mortar, face_mesh, mortar_mesh, mortar_size);
  } else {
    ASSERT(dt_boundary_correction_on_face.number_of_grid_points() ==
               dt_boundary_correction_on_mortar.number_of_grid_points(),
           "Expected identical face and mortar point counts when projection "
           "is not needed.");
    ::Kokkos::deep_copy(dt_boundary_correction_on_face.view(),
                        dt_boundary_correction_on_mortar.view());
  }

  const auto dt_spacetime_metric_on_face =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
          dt_boundary_correction_on_face);
  const auto dt_pi_on_face =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Pi<DataVector, volume_dim,
                                                     Frame::Inertial>>>>(
          dt_boundary_correction_on_face);
  const auto dt_phi_on_face =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Phi<DataVector, volume_dim,
                                                      Frame::Inertial>>>>(
          dt_boundary_correction_on_face);
  const auto dt_spacetime_metric_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
          *dt_boundary_correction_on_face_sum);
  const auto dt_pi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Pi<DataVector, volume_dim,
                                                     Frame::Inertial>>>>(
          *dt_boundary_correction_on_face_sum);
  const auto dt_phi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Phi<DataVector, volume_dim,
                                                      Frame::Inertial>>>>(
          *dt_boundary_correction_on_face_sum);

  ::Kokkos::parallel_for(
      "GhABCTDAccumulateBoundaryCorrectionOnFace",
      face_mesh.number_of_grid_points(),
      KOKKOS_LAMBDA(const int face_index_int) {
        const size_t face_index = static_cast<size_t>(face_index_int);
        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = a; b < volume_dim + 1; ++b) {
            dt_spacetime_metric_boundary_sum.get(a, b)[face_index] +=
                dt_spacetime_metric_on_face.get(a, b)[face_index];
            dt_pi_boundary_sum.get(a, b)[face_index] +=
                dt_pi_on_face.get(a, b)[face_index];
            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi_boundary_sum.get(d, a, b)[face_index] +=
                  dt_phi_on_face.get(d, a, b)[face_index];
            }
          }
        }
      });
}

void lift_face_correction_to_volume(
    const gsl::not_null<device_dt_type*> device_dt,
    const Variables<device_dt_boundary_tags>&
        dt_boundary_correction_on_face_sum,
    const Direction<volume_dim>& direction,
    const device_variables_type& device_vars,
    const device_face_to_volume_index_map_type& device_face_to_volume_index_map,
    const device_face_unit_normal_covector_type& device_face_unit_normal_covector,
    const device_face_normal_magnitude_type& device_face_normal_magnitude,
    const Mesh<volume_dim>& mesh) {
  const auto dt_spacetime_metric =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(*device_dt);
  const auto dt_pi =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Pi<DataVector, volume_dim,
                                                     Frame::Inertial>>>>(
          *device_dt);
  const auto dt_phi =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Phi<DataVector, volume_dim,
                                                      Frame::Inertial>>>>(
          *device_dt);
  const auto dt_spacetime_metric_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(
          dt_boundary_correction_on_face_sum);
  const auto dt_pi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Pi<DataVector, volume_dim,
                                                     Frame::Inertial>>>>(
          dt_boundary_correction_on_face_sum);
  const auto dt_phi_boundary_sum =
      get<::Tags::MirrorView<::Tags::dt<gh::Tags::Phi<DataVector, volume_dim,
                                                      Frame::Inertial>>>>(
          dt_boundary_correction_on_face_sum);
  const auto spacetime_metric =
      get<::Tags::MirrorView<spacetime_metric_tag>>(device_vars);

  const size_t num_face_points =
      dt_boundary_correction_on_face_sum.number_of_grid_points();
  const size_t extent_perpendicular_to_boundary = mesh.extents(direction.dimension());
  const double lift_prefactor =
      -0.5 * static_cast<double>(extent_perpendicular_to_boundary *
                                 (extent_perpendicular_to_boundary - 1));
  const size_t sliced_dim = direction.dimension();

  const auto face_to_volume_index =
      direction.side() == Side::Upper
          ? gsl::at(device_face_to_volume_index_map, sliced_dim).second
          : gsl::at(device_face_to_volume_index_map, sliced_dim).first;
  const auto face_unit_normal_covector =
      direction.side() == Side::Upper
          ? gsl::at(device_face_unit_normal_covector, sliced_dim).second
          : gsl::at(device_face_unit_normal_covector, sliced_dim).first;
  const auto face_normal_magnitude =
      direction.side() == Side::Upper
          ? gsl::at(device_face_normal_magnitude, sliced_dim).second
          : gsl::at(device_face_normal_magnitude, sliced_dim).first;

  ::Kokkos::parallel_for(
      "GhABCTDLiftFaceCorrectionToVolume", num_face_points,
      KOKKOS_LAMBDA(const int face_index_int) {
        const size_t face_index = static_cast<size_t>(face_index_int);
        const size_t volume_index = face_to_volume_index(face_index);

        tnsr::aa<double, volume_dim, Frame::Inertial> spacetime_metric_at_point{};
        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = a; b < volume_dim + 1; ++b) {
            spacetime_metric_at_point.get(a, b) =
                spacetime_metric.get(a, b)[volume_index];
          }
        }
        tnsr::II<double, volume_dim, Frame::Inertial> inverse_spatial_metric{};
        double det_spatial_metric = 0.0;
        inverse_spatial_metric_and_det(make_not_null(&inverse_spatial_metric),
                                       make_not_null(&det_spatial_metric),
                                       spacetime_metric_at_point);
        (void)det_spatial_metric;

        tnsr::i<double, volume_dim, Frame::Inertial> unnormalized_normal_covector{};
        for (size_t d = 0; d < volume_dim; ++d) {
          unnormalized_normal_covector.get(d) =
              face_unit_normal_covector(face_index, d) *
              face_normal_magnitude(face_index);
        }

        tnsr::I<double, volume_dim, Frame::Inertial> unnormalized_normal_vector{};
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
        const double lifted_factor = lift_prefactor * sqrt(normal_magnitude_squared);

        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = a; b < volume_dim + 1; ++b) {
            dt_spacetime_metric.get(a, b)[volume_index] +=
                lifted_factor *
                dt_spacetime_metric_boundary_sum.get(a, b)[face_index];
            dt_pi.get(a, b)[volume_index] +=
                lifted_factor * dt_pi_boundary_sum.get(a, b)[face_index];
            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi.get(d, a, b)[volume_index] +=
                  lifted_factor * dt_phi_boundary_sum.get(d, a, b)[face_index];
            }
          }
        }
      });
}

}  // namespace

void ApplyBoundaryCorrectionsToTimeDerivativeKokkos::
    apply_boundary_corrections_on_device(
        const gsl::not_null<device_dt_type*> device_dt,
        const device_variables_type& device_vars,
        const outgoing_boundary_data_type& outgoing_boundary_data,
        const incoming_boundary_data_type& incoming_boundary_data,
        const external_boundary_data_type& external_boundary_data,
        const device_face_to_volume_index_map_type& device_face_to_volume_index_map,
        const device_face_unit_normal_covector_type&
            device_face_unit_normal_covector,
        const device_face_normal_magnitude_type& device_face_normal_magnitude,
        const device_mortar_data_type& device_mortar_data,
        const typename mortar_mesh_tag::type& mortar_meshes, const Mesh<3>& mesh,
        const Element<3>& element) {
  ASSERT(mesh.quadrature(0) == Spectral::Quadrature::GaussLobatto,
         "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
         "Gauss-Lobatto quadrature only.");

  for (const auto& [direction, neighbors] : element.neighbors()) {
    const Mesh<2> face_mesh = mesh.slice_away(direction.dimension());
    const size_t num_face_points = face_mesh.number_of_grid_points();
    Variables<device_dt_boundary_tags> dt_boundary_correction_on_face_sum{
        num_face_points};
    ::Kokkos::deep_copy(dt_boundary_correction_on_face_sum.view(), 0.0);

    for (const auto& neighbor : neighbors) {
      const DirectionalId<3> mortar_id{direction, neighbor};
      ASSERT(
          outgoing_boundary_data.count(mortar_id) == 1,
          "Missing local Kokkos boundary data for mortar " << mortar_id << ".");
      ASSERT(incoming_boundary_data.count(mortar_id) == 1,
             "Missing received Kokkos boundary data for mortar " << mortar_id
                                                                 << ".");

      const auto& local_boundary_data = outgoing_boundary_data.at(mortar_id);
      const auto& remote_boundary_data = incoming_boundary_data.at(mortar_id);
      const Mesh<2>& mortar_mesh = mortar_meshes.at(mortar_id);
      ASSERT(local_boundary_data.boundary_correction_data.number_of_grid_points() ==
                 mortar_mesh.number_of_grid_points(),
             "Local packaged mortar data size does not match mortar mesh for "
                 << mortar_id << ".");
      ASSERT(remote_boundary_data.boundary_correction_data
                     .number_of_grid_points() == mortar_mesh.number_of_grid_points(),
             "Remote packaged mortar data size does not match mortar mesh for "
                 << mortar_id << ".");
      const auto& mortar_data = device_mortar_data.at(mortar_id);
      accumulate_mortar_pair_to_face_sum(
          make_not_null(&dt_boundary_correction_on_face_sum), local_boundary_data,
          remote_boundary_data, face_mesh, mortar_mesh, mortar_data.needs_projection,
          mortar_data.mortar_size);
    }
    lift_face_correction_to_volume(
        device_dt, dt_boundary_correction_on_face_sum, direction, device_vars,
        device_face_to_volume_index_map, device_face_unit_normal_covector,
        device_face_normal_magnitude, mesh);
  }

  const std::array<Spectral::SegmentSize, volume_dim - 1> full_mortar_size{
      Spectral::SegmentSize::Full, Spectral::SegmentSize::Full};
  for (const Direction<3>& direction : element.external_boundaries()) {
    const Mesh<2> face_mesh = mesh.slice_away(direction.dimension());
    Variables<device_dt_boundary_tags> dt_boundary_correction_on_face_sum{
        face_mesh.number_of_grid_points()};
    ::Kokkos::deep_copy(dt_boundary_correction_on_face_sum.view(), 0.0);

    const DirectionalId<3> external_mortar_id{
        direction, ElementId<3>::external_boundary_id()};
    ASSERT(outgoing_boundary_data.count(external_mortar_id) == 1,
           "Missing local packaged data for external boundary "
               << external_mortar_id << ".");
    ASSERT(external_boundary_data.count(external_mortar_id) == 1,
           "Missing external ghost packaged data for boundary "
               << external_mortar_id << ".");

    accumulate_mortar_pair_to_face_sum(
        make_not_null(&dt_boundary_correction_on_face_sum),
        outgoing_boundary_data.at(external_mortar_id),
        external_boundary_data.at(external_mortar_id), face_mesh, face_mesh, false,
        full_mortar_size);
    lift_face_correction_to_volume(
        device_dt, dt_boundary_correction_on_face_sum, direction, device_vars,
        device_face_to_volume_index_map, device_face_unit_normal_covector,
        device_face_normal_magnitude, mesh);
  }
}

}  // namespace gh::Actions
