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

  // Fused kernel: one team per element; intermediates in team scratch.
  using exec_space = typename Kokkos::View<double**>::execution_space;
  using team_policy = Kokkos::TeamPolicy<exec_space>;
  using member_type = team_policy::member_type;

  const int P = static_cast<int>(points_per_element);
  const int N0 = static_cast<int>(n0);
  const int N1 = static_cast<int>(n1);
  const int N2 = static_cast<int>(n2);
  const int NC = static_cast<int>(num_components);

  // Two ping-pong buffers of length P (doubles) per team.
  const size_t bytes_per_team = 2ull * static_cast<size_t>(P) * sizeof(double);

  team_policy pol(static_cast<int>(num_elements), Kokkos::AUTO);
  Kokkos::parallel_for(
      "GhFilterKokkosBatchedFusedPerElement",
      pol.set_scratch_size(0, Kokkos::PerTeam(bytes_per_team)),
      KOKKOS_LAMBDA(const member_type& team) {
        const int e = team.league_rank();
        const size_t elem_base =
            static_cast<size_t>(e) * static_cast<size_t>(P);

        double* buf0 =
            (double*)team.team_shmem().get_shmem(sizeof(double) * P);
        double* buf1 =
            (double*)team.team_shmem().get_shmem(sizeof(double) * P);

        for (int component = 0; component < NC; ++component) {
          // Load this element/component into buf0
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, P),
                               [&](const int lp) {
                                 const size_t pidx =
                                     elem_base + static_cast<size_t>(lp);
                                 buf0[lp] = vars_view(pidx, component);
                               });
          team.team_barrier();

          // Dim 0: buf0 -> buf1
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, P), [&](const int lp) {
                int t = lp;
                const int i0 = t % N0;
                t /= N0;
                const int i1 = t % N1;
                const int i2 = t / N1;

                double sum = 0.0;
                const int base = N0 * (i1 + N1 * i2);
                for (int k0 = 0; k0 < N0; ++k0) {
                  sum += matrix_dim_0(i0, k0) * buf0[k0 + base];
                }
                buf1[lp] = sum;
              });
          team.team_barrier();

          // Dim 1: buf1 -> buf0 (ping-pong)
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, P), [&](const int lp) {
                int t = lp;
                const int i0 = t % N0;
                t /= N0;
                const int i1 = t % N1;
                const int i2 = t / N1;

                double sum = 0.0;
                for (int k1 = 0; k1 < N1; ++k1) {
                  const int src = i0 + N0 * (k1 + N1 * i2);
                  sum += matrix_dim_1(i1, k1) * buf1[src];
                }
                buf0[lp] = sum;
              });
          team.team_barrier();

          // Dim 2: buf0 -> vars_view
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, P), [&](const int lp) {
                int t = lp;
                const int i0 = t % N0;
                t /= N0;
                const int i1 = t % N1;
                const int i2 = t / N1;

                double sum = 0.0;
                for (int k2 = 0; k2 < N2; ++k2) {
                  const int src = i0 + N0 * (i1 + N1 * k2);
                  sum += matrix_dim_2(i2, k2) * buf0[src];
                }
                const size_t pidx = elem_base + static_cast<size_t>(lp);
                vars_view(pidx, component) = sum;
              });
          team.team_barrier();
        }
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
