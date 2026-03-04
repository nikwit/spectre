// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/FilterKokkosBatched.hpp"

#include <cstddef>

#include "DataStructures/Matrix.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/MaximumNumberOfPoints.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/StaticCache.hpp"

namespace gh::Actions {

namespace detail {

using MatrixViewRightLayout =
    Kokkos::View<const double**, Kokkos::LayoutRight,
                 typename Kokkos::DefaultExecutionSpace::memory_space,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

MatrixViewRightLayout matrix_on_device(const Matrix& matrix,
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

void apply_filter_matrices_on_device_batched(
    const Kokkos::View<double**>& vars_view, const Mesh<3>& mesh,
    const MatrixViewRightLayout& matrix_dim_0,
    const MatrixViewRightLayout& matrix_dim_1,
    const MatrixViewRightLayout& matrix_dim_2) {
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

  using exec_space = typename Kokkos::View<double**>::execution_space;
  using team_policy = Kokkos::TeamPolicy<exec_space>;
  using member_type = team_policy::member_type;

  const int p = static_cast<int>(points_per_element);
  const int n_0 = static_cast<int>(n0);
  const int n_1 = static_cast<int>(n1);
  const int n_2 = static_cast<int>(n2);
  const int num_components_int = static_cast<int>(num_components);
  const size_t bytes_per_team = 2 * static_cast<size_t>(p) * sizeof(double);

  team_policy policy(static_cast<int>(num_elements), Kokkos::AUTO);
  Kokkos::parallel_for(
      "GhFilterKokkosBatchedFusedPerElement",
      policy.set_scratch_size(0, Kokkos::PerTeam(bytes_per_team)),
      KOKKOS_LAMBDA(const member_type& team) {
        const int element_index = team.league_rank();
        const size_t element_base =
            static_cast<size_t>(element_index) * static_cast<size_t>(p);

        double* buffer_0 = static_cast<double*>(team.team_shmem().get_shmem(
            sizeof(double) * static_cast<size_t>(p)));
        double* buffer_1 = static_cast<double*>(team.team_shmem().get_shmem(
            sizeof(double) * static_cast<size_t>(p)));

        for (int component = 0; component < num_components_int; ++component) {
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, p), [&](const int local_point) {
                const size_t point_index =
                    element_base + static_cast<size_t>(local_point);
                buffer_0[local_point] = vars_view(point_index, component);
              });
          team.team_barrier();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, p), [&](const int local_point) {
                int tensor_index = local_point;
                const int i_0 = tensor_index % n_0;
                tensor_index /= n_0;
                const int i_1 = tensor_index % n_1;
                const int i_2 = tensor_index / n_1;

                double sum = 0.0;
                const int base = n_0 * (i_1 + n_1 * i_2);
                for (int k_0 = 0; k_0 < n_0; ++k_0) {
                  sum += matrix_dim_0(i_0, k_0) * buffer_0[k_0 + base];
                }
                buffer_1[local_point] = sum;
              });
          team.team_barrier();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, p), [&](const int local_point) {
                int tensor_index = local_point;
                const int i_0 = tensor_index % n_0;
                tensor_index /= n_0;
                const int i_1 = tensor_index % n_1;
                const int i_2 = tensor_index / n_1;

                double sum = 0.0;
                for (int k_1 = 0; k_1 < n_1; ++k_1) {
                  const int source = i_0 + n_0 * (k_1 + n_1 * i_2);
                  sum += matrix_dim_1(i_1, k_1) * buffer_1[source];
                }
                buffer_0[local_point] = sum;
              });
          team.team_barrier();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, p), [&](const int local_point) {
                int tensor_index = local_point;
                const int i_0 = tensor_index % n_0;
                tensor_index /= n_0;
                const int i_1 = tensor_index % n_1;
                const int i_2 = tensor_index / n_1;

                double sum = 0.0;
                for (int k_2 = 0; k_2 < n_2; ++k_2) {
                  const int source = i_0 + n_0 * (i_1 + n_1 * k_2);
                  sum += matrix_dim_2(i_2, k_2) * buffer_0[source];
                }
                const size_t point_index =
                    element_base + static_cast<size_t>(local_point);
                vars_view(point_index, component) = sum;
              });
          team.team_barrier();
        }
      });
}

}  // namespace detail

template <typename FilterType>
void FilterKokkosBatched<FilterType>::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
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
      matrix_view_1, matrix_view_2);
}

template struct FilterKokkosBatched<Filters::Exponential<0>>;

}  // namespace gh::Actions
