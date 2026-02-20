// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Structure/BlockGroups.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Tags/Filter.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/MaximumNumberOfPoints.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/StaticCache.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace gh::Actions {

namespace detail {

inline MatrixViewRO matrix_on_device(const Matrix& matrix,
                                     const char* const label) {
  Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace>
      host_matrix_view("GhFilterKokkosHostMatrix", matrix.rows(),
                       matrix.columns());
  for (size_t i = 0; i < matrix.rows(); ++i) {
    for (size_t j = 0; j < matrix.columns(); ++j) {
      host_matrix_view(i, j) = matrix(i, j);
    }
  }
  Kokkos::View<double**, Kokkos::LayoutLeft> device_matrix_view(
      label, matrix.rows(), matrix.columns());
  Kokkos::deep_copy(device_matrix_view, host_matrix_view);
  return MatrixViewRO{device_matrix_view};
}

template <typename FilterType>
const MatrixViewRO& filter_matrix_on_device(const FilterType& filter_helper,
                                            const Mesh<1>& mesh) {
  const static FilterType cached_filter = filter_helper;
  ASSERT(cached_filter == filter_helper,
         "FilterKokkos device matrix cache was initialized with different "
         "filter parameters. Use a different FilterIndex for a different "
         "filter configuration.");

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
                                "GhFilterKokkosMatrix");
      });

  return cache(mesh.extents(0), mesh.basis(0), mesh.quadrature(0));
}

void apply_filter_matrices_on_device(const Kokkos::View<double**>& vars_view,
                                     const Mesh<3>& mesh,
                                     const MatrixViewRO& matrix_dim_0,
                                     const MatrixViewRO& matrix_dim_1,
                                     const MatrixViewRO& matrix_dim_2) {
  const auto extents = mesh.extents();
  const size_t num_components = vars_view.extent(1);
  const size_t total_points = extents.product();
  if (total_points == 0 or num_components == 0) {
    return;
  }

  const size_t n0 = extents[0];
  const size_t n1 = extents[1];
  const size_t n2 = extents[2];
  Kokkos::View<double**> scratch_0("GhFilterKokkosScratch0", total_points,
                                   num_components);
  Kokkos::View<double**> scratch_1("GhFilterKokkosScratch1", total_points,
                                   num_components);

  Kokkos::parallel_for(
      "GhFilterKokkosDim0",
      Kokkos::RangePolicy<size_t>{0, total_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i0 = remaining % n0;
        remaining /= n0;
        const size_t i1 = remaining % n1;
        const size_t i2 = remaining / n1;

        double sum = 0.0;
        for (size_t k0 = 0; k0 < n0; ++k0) {
          sum += matrix_dim_0(i0, k0) *
                 vars_view(k0 + n0 * (i1 + n1 * i2), component);
        }
        scratch_0(i0 + n0 * (i1 + n1 * i2), component) = sum;
      });

  Kokkos::parallel_for(
      "GhFilterKokkosDim1",
      Kokkos::RangePolicy<size_t>{0, total_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i0 = remaining % n0;
        remaining /= n0;
        const size_t i1 = remaining % n1;
        const size_t i2 = remaining / n1;

        double sum = 0.0;
        for (size_t k1 = 0; k1 < n1; ++k1) {
          sum += matrix_dim_1(i1, k1) *
                 scratch_0(i0 + n0 * (k1 + n1 * i2), component);
        }
        scratch_1(i0 + n0 * (i1 + n1 * i2), component) = sum;
      });

  Kokkos::parallel_for(
      "GhFilterKokkosDim2",
      Kokkos::RangePolicy<size_t>{0, total_points * num_components},
      KOKKOS_LAMBDA(const size_t linear_index) {
        size_t remaining = linear_index;
        const size_t component = remaining % num_components;
        remaining /= num_components;
        const size_t i0 = remaining % n0;
        remaining /= n0;
        const size_t i1 = remaining % n1;
        const size_t i2 = remaining / n1;

        double sum = 0.0;
        for (size_t k2 = 0; k2 < n2; ++k2) {
          sum += matrix_dim_2(i2, k2) *
                 scratch_1(i0 + n0 * (i1 + n1 * k2), component);
        }
        vars_view(i0 + n0 * (i1 + n1 * i2), component) = sum;
      });
}

}  // namespace detail

template <typename FilterType>
struct FilterKokkos {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using device_variables_tag = gh::KokkosTags::DeviceVariables<system>;
  using device_variables_type = typename device_variables_tag::type;

 public:
  using return_tags = tmpl::list<device_variables_tag>;
  using argument_tags = tmpl::list<domain::Tags::Mesh<volume_dim>,
                                   domain::Tags::Element<volume_dim>>;
  using inbox_tags = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<::domain::Tags::Domain<volume_dim>,
                 ::Filters::Tags::Filter<FilterType>>;

  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent,
            typename Metavariables>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const FilterType& filter_helper =
        Parallel::get<::Filters::Tags::Filter<FilterType>>(cache);
    const auto& element = db::get<domain::Tags::Element<volume_dim>>(box);
    const auto& mesh = db::get<domain::Tags::Mesh<volume_dim>>(box);

    const auto& domain = Parallel::get<domain::Tags::Domain<volume_dim>>(cache);
    const auto& block_groups = domain.block_groups();
    const std::string& block_name =
        domain.blocks()[element.id().block_id()].name();

    bool enable = filter_helper.enable();
    if (enable and filter_helper.blocks_to_filter().has_value()) {
      enable = alg::any_of(
          filter_helper.blocks_to_filter().value(),
          [&block_name, &block_groups](const std::string& block_to_filter) {
            return domain::block_is_in_group(block_name, block_to_filter,
                                             block_groups);
          });
    }
    if (not enable) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& matrix_view_0 =
        detail::filter_matrix_on_device(filter_helper, mesh.slice_through(0));
    const auto& matrix_view_1 =
        detail::filter_matrix_on_device(filter_helper, mesh.slice_through(1));
    const auto& matrix_view_2 =
        detail::filter_matrix_on_device(filter_helper, mesh.slice_through(2));

    db::mutate<device_variables_tag>(
        [&mesh, &matrix_view_0, &matrix_view_1, &matrix_view_2](
            const gsl::not_null<device_variables_type*> device_variables) {
          detail::apply_filter_matrices_on_device(device_variables->view(),
                                                  mesh, matrix_view_0,
                                                  matrix_view_1, matrix_view_2);
        },
        make_not_null(&box));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace gh::Actions
