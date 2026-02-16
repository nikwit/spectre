// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/VariablesKokkos.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/SegmentSize.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Kokkos {

namespace detail {

template <typename TagsList>
void apply_tensor_product_projection_2d(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& input, const MatrixViewRO& matrix_dim_0,
    const MatrixViewRO& matrix_dim_1) {
  const size_t source_points_dim_0 = matrix_dim_0.extent(1);
  const size_t source_points_dim_1 = matrix_dim_1.extent(1);
  const size_t target_points_dim_0 = matrix_dim_0.extent(0);
  const size_t target_points_dim_1 = matrix_dim_1.extent(0);
  const auto input_view = input.view();
  auto result_view = result->view();
  const size_t num_components = input_view.extent(1);

  ASSERT(input.number_of_grid_points() ==
             source_points_dim_0 * source_points_dim_1,
         "Input has " << input.number_of_grid_points() << " points, expected "
                      << source_points_dim_0 << " * " << source_points_dim_1
                      << ".");
  ASSERT(result->number_of_grid_points() ==
             target_points_dim_0 * target_points_dim_1,
         "Result has " << result->number_of_grid_points()
                       << " points, expected " << target_points_dim_0 << " * "
                       << target_points_dim_1 << ".");

  ::Kokkos::View<double**> projected_dim_0(
      "ProjectedDim0", target_points_dim_0 * source_points_dim_1,
      num_components);
  ::Kokkos::parallel_for(
      "ProjectDim0",
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
      "ProjectDim1",
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

}  // namespace detail

template <typename TagsList>
void orient_each_component_device(
    const gsl::not_null<Variables<TagsList>*> oriented_variables,
    const Variables<TagsList>& variables,
    const ::Kokkos::View<size_t*>& oriented_offset) {
  const size_t num_points = variables.number_of_grid_points();
  ASSERT(oriented_variables->number_of_grid_points() == num_points,
         "Oriented result has " << oriented_variables->number_of_grid_points()
                                << " points, expected " << num_points << ".");
  ASSERT(oriented_offset.extent(0) == num_points,
         "Orientation map has " << oriented_offset.extent(0)
                                << " points, expected " << num_points << ".");

  auto oriented_view = oriented_variables->view();
  const auto variables_view = variables.view();
  const auto oriented_offset_view = oriented_offset;
  const size_t num_components = variables_view.extent(1);
  ::Kokkos::parallel_for(
      "OrientEachComponentOnDevice",
      ::Kokkos::RangePolicy<size_t>{0, num_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        const size_t component = linear_index % num_components;
        const size_t source_point = linear_index / num_components;
        const size_t target_point = oriented_offset_view(source_point);
        oriented_view(target_point, component) =
            variables_view(source_point, component);
      });
}

template <typename TagsList>
void project_to_mortar_device(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& vars, const Mesh<2>& face_mesh,
    const Mesh<2>& mortar_mesh,
    const std::array<Spectral::SegmentSize, 2>& mortar_size) {
  if (not Spectral::needs_projection(face_mesh, mortar_mesh, mortar_size)) {
    ASSERT(result->number_of_grid_points() == vars.number_of_grid_points(),
           "Cannot skip projection for incompatible result size.");
    ::Kokkos::deep_copy(result->view(), vars.view());
    return;
  }

  ASSERT(result->number_of_grid_points() == mortar_mesh.number_of_grid_points(),
         "Result has " << result->number_of_grid_points()
                       << " points, expected "
                       << mortar_mesh.number_of_grid_points() << ".");
  const auto& matrix_dim_0 =
      Spectral::projection_matrix_parent_to_child_on_device(
          face_mesh.slice_through(0), mortar_mesh.slice_through(0),
          gsl::at(mortar_size, 0));
  const auto& matrix_dim_1 =
      Spectral::projection_matrix_parent_to_child_on_device(
          face_mesh.slice_through(1), mortar_mesh.slice_through(1),
          gsl::at(mortar_size, 1));
  detail::apply_tensor_product_projection_2d(result, vars, matrix_dim_0,
                                             matrix_dim_1);
}

template <typename TagsList>
void project_from_mortar_device(
    const gsl::not_null<Variables<TagsList>*> result,
    const Variables<TagsList>& vars, const Mesh<2>& face_mesh,
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
  detail::apply_tensor_product_projection_2d(result, vars, matrix_dim_0,
                                             matrix_dim_1);
}

}  // namespace ScalarWave::Kokkos
