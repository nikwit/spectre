// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/Matrix.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Tags/Filter.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/MaximumNumberOfPoints.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/StaticCache.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

namespace detail {

using MatrixViewRightLayout =
    Kokkos::View<const double**, Kokkos::LayoutRight,
                 typename Kokkos::DefaultExecutionSpace::memory_space,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
inline MatrixViewRightLayout matrix_on_device(const Matrix& matrix,
                                              const char* const label) {
  Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>
      host_matrix_view{"GhFilterKokkosBatchedHostMatrix", matrix.rows(),
                       matrix.columns()};
  for (size_t i = 0; i < matrix.rows(); ++i) {
    for (size_t j = 0; j < matrix.columns(); ++j) {
      host_matrix_view(i, j) = matrix(i, j);
    }
  }
  Kokkos::View<double**, Kokkos::LayoutRight> device_matrix_view{
      label, matrix.rows(), matrix.columns()};
  Kokkos::deep_copy(device_matrix_view, host_matrix_view);
  return MatrixViewRightLayout{device_matrix_view};
}

template <typename FilterType>
const MatrixViewRightLayout& filter_matrix_on_device_batched(
    const FilterType& filter_helper, const Mesh<1>& mesh) {
  const static FilterType cached_filter = filter_helper;
  ASSERT(cached_filter == filter_helper,
         "FilterKokkosBatched device matrix cache was initialized with "
         "different filter parameters. Use a different FilterIndex for a "
         "different filter configuration.");

  const static auto cache = make_static_cache<
      CacheRange<1_st,
                 Spectral::maximum_number_of_points<Spectral::Basis::Legendre> +
                     1>,
      CacheEnumeration<Spectral::Basis, Spectral::Basis::Legendre,
                       Spectral::Basis::Chebyshev>,
      CacheEnumeration<Spectral::Quadrature, Spectral::Quadrature::Gauss,
                       Spectral::Quadrature::GaussLobatto>>(
      [filter = filter_helper](const size_t extents,
                               const Spectral::Basis basis,
                               const Spectral::Quadrature quadrature) {
        const Mesh<1> matrix_mesh{extents, basis, quadrature};
        return matrix_on_device(filter.filter_matrix(matrix_mesh),
                                "GhFilterKokkosBatchedMatrix");
      });

  return cache(mesh.extents(0), mesh.basis(0), mesh.quadrature(0));
}

inline void apply_filter_matrices_on_device_batched(
    const Kokkos::View<double**>& vars_view, const Mesh<3>& mesh,
    const MatrixViewRightLayout& matrix_dim_0,
    const MatrixViewRightLayout& matrix_dim_1,
    const MatrixViewRightLayout& matrix_dim_2,
    const gsl::not_null<Kokkos::View<double**>*> scratch_0,
    const gsl::not_null<Kokkos::View<double**>*> scratch_1) {
  const auto extents = mesh.extents();
  const size_t n0 = extents[0];
  const size_t n1 = extents[1];
  const size_t n2 = extents[2];
  const size_t points_per_element = extents.product();
  const size_t total_points = vars_view.extent(0);
  const size_t num_components = vars_view.extent(1);
  if (points_per_element == 0 or total_points == 0 or num_components == 0) {
    return;
  }
  ASSERT(total_points % points_per_element == 0,
         "Packed GH filter expected total points (" << total_points
                                                    << ") to be divisible by "
                                                       "points per element ("
                                                    << points_per_element
                                                    << ").");
  const size_t num_elements = total_points / points_per_element;

  if (scratch_0->extent(0) != total_points or
      scratch_0->extent(1) != num_components) {
    *scratch_0 = Kokkos::View<double**>{"GhFilterKokkosBatchedScratch0",
                                        total_points, num_components};
  }
  if (scratch_1->extent(0) != total_points or
      scratch_1->extent(1) != num_components) {
    *scratch_1 = Kokkos::View<double**>{"GhFilterKokkosBatchedScratch1",
                                        total_points, num_components};
  }
  const auto scratch_0_view = *scratch_0;
  const auto scratch_1_view = *scratch_1;

  Kokkos::parallel_for(
      "GhFilterKokkosBatchedDim0_MD",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                             {total_points, num_components}),
      KOKKOS_LAMBDA(const size_t point_index, const size_t component) {
        const int element_index =
            static_cast<int>(point_index / points_per_element);
        int local_point = static_cast<int>(point_index % points_per_element);

        const int i0 = local_point % n0;
        local_point /= n0;
        const int i1 = local_point % n1;
        const int i2 = local_point / n1;

        double sum = 0.0;
        for (int k0 = 0; k0 < n0; ++k0) {
          const int local_source_index = k0 + n0 * (i1 + n1 * i2);
          const size_t source_index =
              static_cast<size_t>(element_index) * points_per_element +
              static_cast<size_t>(local_source_index);

          sum += matrix_dim_0(i0, k0) * vars_view(source_index, component);
        }
        scratch_0_view(point_index, component) = sum;
      });

  Kokkos::parallel_for(
      "GhFilterKokkosBatchedDim1_MD",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                             {total_points, num_components}),
      KOKKOS_LAMBDA(const size_t point_index, const size_t component) {
        const int element_index = point_index / points_per_element;
        int local_point = point_index % points_per_element;

        const int i0 = local_point % n0;
        local_point /= n0;
        const int i1 = local_point % n1;
        const int i2 = local_point / n1;

        double sum = 0.0;
        for (int k1 = 0; k1 < n1; ++k1) {
          const int local_source_index = i0 + n0 * (k1 + n1 * i2);
          const size_t source_index =
              static_cast<size_t>(element_index) * points_per_element +
              static_cast<size_t>(local_source_index);
          sum += matrix_dim_1(i1, k1) * scratch_0_view(source_index, component);
        }
        scratch_1_view(point_index, component) = sum;
      });

  Kokkos::parallel_for(
      "GhFilterKokkosBatchedDim2_MD",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0},
                                             {total_points, num_components}),
      KOKKOS_LAMBDA(const size_t point_index, const size_t component) {
        const int element_index =
            static_cast<int>(point_index / points_per_element);
        int local_point = static_cast<int>(point_index % points_per_element);

        const int i0 = local_point % n0;
        local_point /= n0;
        const int i1 = local_point % n1;
        const int i2 = local_point / n1;

        double sum = 0.0;
        for (int k2 = 0; k2 < n2; ++k2) {
          const int local_source_index = i0 + n0 * (i1 + n1 * k2);
          const size_t source_index =
              static_cast<size_t>(element_index) * points_per_element +
              static_cast<size_t>(local_source_index);

          sum += matrix_dim_2(i2, k2) * scratch_1_view(source_index, component);
        }
        vars_view(point_index, component) = sum;
      });
  (void)num_elements;
}

}  // namespace detail

template <typename FilterType>
struct FilterKokkosBatched {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_evolution_state_tag =
      evolution::Kokkos::Tags::PackedEvolutionState<system>;
  using packed_boundary_scratch_tag =
      evolution::Kokkos::Tags::PackedBoundaryScratch<system>;
  using packed_topology_tag = evolution::Kokkos::Tags::PackedTopology<system>;

 public:
  using return_tags =
      tmpl::list<packed_evolution_state_tag, packed_boundary_scratch_tag>;
  using argument_tags =
      tmpl::list<packed_topology_tag, ::Filters::Tags::Filter<FilterType>>;
  using const_global_cache_tags =
      tmpl::list<::Filters::Tags::Filter<FilterType>>;

  static void apply(
      const gsl::not_null<typename packed_evolution_state_tag::type*>
          packed_evolution_state,
      const gsl::not_null<typename packed_boundary_scratch_tag::type*>
          packed_boundary_scratch,
      const typename packed_topology_tag::type& packed_topology,
      const FilterType& filter_helper) {
    if (packed_topology.total_points == 0 or not filter_helper.enable()) {
      return;
    }
    ASSERT(not filter_helper.blocks_to_filter().has_value(),
           "FilterKokkosBatched does not yet support BlocksToFilter.");

    const Mesh<3> mesh{packed_topology.uniform_extents_host,
                       packed_topology.uniform_basis_host,
                       packed_topology.uniform_quadrature_host};
    const auto& matrix_view_0 = detail::filter_matrix_on_device_batched(
        filter_helper, mesh.slice_through(0));
    const auto& matrix_view_1 = detail::filter_matrix_on_device_batched(
        filter_helper, mesh.slice_through(1));
    const auto& matrix_view_2 = detail::filter_matrix_on_device_batched(
        filter_helper, mesh.slice_through(2));

    detail::apply_filter_matrices_on_device_batched(
        packed_evolution_state->device_variables.view(), mesh, matrix_view_0,
        matrix_view_1, matrix_view_2,
        make_not_null(&packed_boundary_scratch->filter_workspace_0),
        make_not_null(&packed_boundary_scratch->filter_workspace_1));
  }
};

}  // namespace gh::Actions
