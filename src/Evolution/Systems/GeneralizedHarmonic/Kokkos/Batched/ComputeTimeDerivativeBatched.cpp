// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeTimeDerivativeBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

constexpr size_t dim = 3;
constexpr size_t spacetime_dim = dim + 1;
constexpr size_t spacetime_symmetric_size =
    spacetime_dim * (spacetime_dim + 1) / 2;
using system = gh::System<dim>;
using packed_evolution_state_type =
    evolution::Kokkos::Tags::PackedEvolutionState<system>::type;
using packed_boundary_scratch_type =
    evolution::Kokkos::Tags::PackedBoundaryScratch<system>::type;
using packed_topology_type =
    evolution::Kokkos::Tags::PackedTopology<system>::type;
using packed_geometry_type =
    evolution::Kokkos::Tags::PackedGeometry<system>::type;
using device_constraint_gamma0_type =
    gh::KokkosTags::DeviceConstraintGamma0::type;
using device_constraint_gamma1_type =
    gh::KokkosTags::DeviceConstraintGamma1::type;
using device_constraint_gamma2_type =
    gh::KokkosTags::DeviceConstraintGamma2::type;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, dim, Frame::Inertial>;
using device_spacetime_metric_tag = ::Tags::MirrorView<spacetime_metric_tag>;
using device_pi_tag = ::Tags::MirrorView<pi_tag>;
using device_phi_tag = ::Tags::MirrorView<phi_tag>;
using device_gauge_h_tag =
    ::Tags::MirrorView<gh::Tags::GaugeH<DataVector, dim>>;
using device_spacetime_deriv_gauge_h_tag =
    ::Tags::MirrorView<gh::Tags::SpacetimeDerivGaugeH<DataVector, dim>>;
using device_gauge_data_type = Variables<
    tmpl::list<device_gauge_h_tag, device_spacetime_deriv_gauge_h_tag>>;
using device_gradient_tags =
    db::wrap_tags_in<::Tags::MirrorView, typename system::gradient_variables>;
using device_derivative_tags =
    db::wrap_tags_in<::Tags::deriv, device_gradient_tags, tmpl::size_t<dim>,
                     Frame::Inertial>;
using device_d_spacetime_metric_tag =
    ::Tags::deriv<device_spacetime_metric_tag, tmpl::size_t<dim>,
                  Frame::Inertial>;
using device_d_pi_tag =
    ::Tags::deriv<device_pi_tag, tmpl::size_t<dim>, Frame::Inertial>;
using device_d_phi_tag =
    ::Tags::deriv<device_phi_tag, tmpl::size_t<dim>, Frame::Inertial>;
using device_inverse_spatial_metric_tag =
    typename packed_boundary_scratch_type::device_inverse_spatial_metric_tag;
using device_shift_tag =
    typename packed_boundary_scratch_type::device_shift_tag;
using device_lapse_tag =
    typename packed_boundary_scratch_type::device_lapse_tag;

KOKKOS_INLINE_FUNCTION void inverse_spatial_metric_and_det(
    const gsl::not_null<tnsr::II<double, dim, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<double*> det_spatial_metric,
    const tnsr::aa<double, dim, Frame::Inertial>& spacetime_metric) {
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

KOKKOS_INLINE_FUNCTION void compute_inverse_spatial_metric_shift_lapse(
    const gsl::not_null<tnsr::II<double, dim, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<tnsr::I<double, dim, Frame::Inertial>*> shift,
    const gsl::not_null<double*> lapse,
    const gsl::not_null<double*> inv_lapse_squared,
    const tnsr::aa<double, dim, Frame::Inertial>& spacetime_metric) {
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det(inverse_spatial_metric,
                                 make_not_null(&det_spatial_metric),
                                 spacetime_metric);
  (void)det_spatial_metric;

  for (size_t i = 0; i < dim; ++i) {
    shift->get(i) = 0.0;
    for (size_t j = 0; j < dim; ++j) {
      shift->get(i) +=
          inverse_spatial_metric->get(i, j) * spacetime_metric.get(0, j + 1);
    }
  }

  double lapse_squared = -spacetime_metric.get(0, 0);
  for (size_t i = 0; i < dim; ++i) {
    lapse_squared += shift->get(i) * spacetime_metric.get(0, i + 1);
  }
  *lapse = sqrt(lapse_squared);
  *inv_lapse_squared = 1.0 / lapse_squared;
}

KOKKOS_INLINE_FUNCTION void compute_inverse_spacetime_metric(
    const gsl::not_null<tnsr::AA<double, dim, Frame::Inertial>*>
        inverse_spacetime_metric,
    const tnsr::II<double, dim, Frame::Inertial>& inverse_spatial_metric,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const double inv_lapse_squared) {
  inverse_spacetime_metric->get(0, 0) = -inv_lapse_squared;
  for (size_t i = 0; i < dim; ++i) {
    inverse_spacetime_metric->get(0, i + 1) = shift.get(i) * inv_lapse_squared;
    for (size_t j = i; j < dim; ++j) {
      inverse_spacetime_metric->get(i + 1, j + 1) =
          inverse_spatial_metric.get(i, j) -
          shift.get(i) * shift.get(j) * inv_lapse_squared;
    }
  }
}

template <std::size_t N, class F, std::size_t... Is>
KOKKOS_INLINE_FUNCTION constexpr void static_for_impl(
    F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t N, class F>
KOKKOS_INLINE_FUNCTION constexpr void static_for(F&& f) {
  static_for_impl<N>(static_cast<F&&>(f), std::make_index_sequence<N>{});
}
KOKKOS_INLINE_FUNCTION size_t symmetric_spacetime_index(const size_t a,
                                                        const size_t b) {
  const size_t mu = a < b ? a : b;
  const size_t nu = a < b ? b : a;
  return mu * (2 * spacetime_dim - mu + 1) / 2 + (nu - mu);
}

KOKKOS_INLINE_FUNCTION std::array<size_t, 2> symmetric_spacetime_indices(
    const size_t symmetric_index) {
  size_t mu = 0;
  size_t remaining = symmetric_index;
  size_t row_width = spacetime_dim;
  while (remaining >= row_width) {
    remaining -= row_width;
    ++mu;
    --row_width;
  }
  return {{mu, mu + remaining}};
}

KOKKOS_INLINE_FUNCTION double get_symmetric_spacetime_component(
    const double* symmetric_components, const size_t a, const size_t b) {
  return symmetric_components[symmetric_spacetime_index(a, b)];
}

KOKKOS_INLINE_FUNCTION double get_phi_component(const double* phi_components,
                                                const size_t i, const size_t a,
                                                const size_t b) {
  return phi_components[i * spacetime_symmetric_size +
                        symmetric_spacetime_index(a, b)];
}

KOKKOS_INLINE_FUNCTION double get_d_phi_component(
    const double* d_phi_components, const size_t m, const size_t i,
    const size_t a, const size_t b) {
  return d_phi_components[(m * dim + i) * spacetime_symmetric_size +
                          symmetric_spacetime_index(a, b)];
}

template <typename TensorType>
KOKKOS_INLINE_FUNCTION constexpr std::array<size_t, dim>
make_spatial_vector_storage_indices() {
  std::array<size_t, dim> result{};
  for (size_t i = 0; i < dim; ++i) {
    result[i] = TensorType::get_storage_index(i);
  }
  return result;
}

template <typename TensorType>
KOKKOS_INLINE_FUNCTION constexpr std::array<size_t, dim * dim>
make_spatial_matrix_storage_indices() {
  std::array<size_t, dim * dim> result{};
  size_t linear_index = 0;
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      result[linear_index] = TensorType::get_storage_index(i, j);
      ++linear_index;
    }
  }
  return result;
}

template <typename TensorType>
KOKKOS_INLINE_FUNCTION constexpr std::array<size_t, spacetime_symmetric_size>
make_spacetime_symmetric_storage_indices() {
  std::array<size_t, spacetime_symmetric_size> result{};
  size_t linear_index = 0;
  for (size_t mu = 0; mu < spacetime_dim; ++mu) {
    for (size_t nu = mu; nu < spacetime_dim; ++nu) {
      result[linear_index] = TensorType::get_storage_index(mu, nu);
      ++linear_index;
    }
  }
  return result;
}

template <typename TensorType>
KOKKOS_INLINE_FUNCTION constexpr std::array<size_t,
                                            dim * spacetime_symmetric_size>
make_spatial_spacetime_symmetric_storage_indices() {
  std::array<size_t, dim * spacetime_symmetric_size> result{};
  size_t linear_index = 0;
  for (size_t i = 0; i < dim; ++i) {
    for (size_t mu = 0; mu < spacetime_dim; ++mu) {
      for (size_t nu = mu; nu < spacetime_dim; ++nu) {
        result[linear_index] = TensorType::get_storage_index(i, mu, nu);
        ++linear_index;
      }
    }
  }
  return result;
}

template <typename TensorType>
KOKKOS_INLINE_FUNCTION constexpr std::array<
    size_t, dim * dim * spacetime_symmetric_size>
make_spatial_spatial_spacetime_symmetric_storage_indices() {
  std::array<size_t, dim * dim * spacetime_symmetric_size> result{};
  size_t linear_index = 0;
  for (size_t m = 0; m < dim; ++m) {
    for (size_t i = 0; i < dim; ++i) {
      for (size_t mu = 0; mu < spacetime_dim; ++mu) {
        for (size_t nu = mu; nu < spacetime_dim; ++nu) {
          result[linear_index] = TensorType::get_storage_index(m, i, mu, nu);
          ++linear_index;
        }
      }
    }
  }
  return result;
}

KOKKOS_INLINE_FUNCTION double spacetime_metric_deriv_component_from_scratch(
    const size_t a, const size_t mu, const size_t nu, const double lapse,
    const double* shift, const double* pi_components,
    const double* phi_components) {
  if (a == 0) {
    double dt_spacetime_metric_mu_nu =
        -lapse * get_symmetric_spacetime_component(pi_components, mu, nu);
    for (size_t m = 0; m < dim; ++m) {
      dt_spacetime_metric_mu_nu +=
          shift[m] * get_phi_component(phi_components, m, mu, nu);
    }
    return dt_spacetime_metric_mu_nu;
  }
  return get_phi_component(phi_components, a - 1, mu, nu);
}

KOKKOS_INLINE_FUNCTION double christoffel_first_kind_component_from_scratch(
    const size_t k, const size_t i, const size_t j, const double lapse,
    const double* shift, const double* pi_components,
    const double* phi_components) {
  return 0.5 * (spacetime_metric_deriv_component_from_scratch(
                    i, j, k, lapse, shift, pi_components, phi_components) +
                spacetime_metric_deriv_component_from_scratch(
                    j, i, k, lapse, shift, pi_components, phi_components) -
                spacetime_metric_deriv_component_from_scratch(
                    k, i, j, lapse, shift, pi_components, phi_components));
}

}  // namespace

void ComputeTimeDerivativeBatched::compute_time_derivative_batched_volume_impl(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename device_constraint_gamma0_tag::type& device_constraint_gamma0,
    const typename device_constraint_gamma1_tag::type& device_constraint_gamma1,
    const typename device_constraint_gamma2_tag::type& device_constraint_gamma2,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry) {
  const size_t total_points = packed_topology.total_points;
  if (total_points == 0) {
    return;
  }

  const auto extents = packed_topology.uniform_extents_host;
  const size_t points_per_element = packed_topology.points_per_element;
  ASSERT(extents[0] * extents[1] * extents[2] == points_per_element,
         "Uniform extents and points-per-element mismatch.");
  ASSERT(packed_topology.local_elements.size() * points_per_element ==
             total_points,
         "Packed point count mismatch in batched volume derivative.");
  const Mesh<dim> mesh{extents, packed_topology.uniform_basis_host,
                       packed_topology.uniform_quadrature_host};

  auto& device_vars = packed_evolution_state->device_variables;
  auto& device_dt = packed_evolution_state->device_dt_variables;

  auto& partial_derivatives_all_elements =
      packed_boundary_scratch->volume_partial_derivatives;
  if (partial_derivatives_all_elements.number_of_grid_points() !=
      total_points) {
    partial_derivatives_all_elements.initialize(total_points);
  }
  partial_derivatives_batched(make_not_null(&partial_derivatives_all_elements),
                              device_vars, mesh,
                              packed_geometry.element_inverse_jacobian_device);

  ASSERT(packed_geometry.inertial_coordinates_device.get(0).extent(0) ==
                 total_points and
             packed_geometry.inertial_coordinates_device.get(1).extent(0) ==
                 total_points and
             packed_geometry.inertial_coordinates_device.get(2).extent(0) ==
                 total_points,
         "Packed inertial device coordinates size mismatch with total points.");
  ASSERT(
      packed_geometry.mesh_velocity_device.get(0).extent(0) == total_points and
          packed_geometry.mesh_velocity_device.get(1).extent(0) ==
              total_points and
          packed_geometry.mesh_velocity_device.get(2).extent(0) == total_points,
      "Packed mesh-velocity device coordinates size mismatch with total "
      "points.");

  auto& device_gauge_data = packed_boundary_scratch->volume_gauge_data;
  if (device_gauge_data.number_of_grid_points() != total_points) {
    device_gauge_data.initialize(total_points);
  }
  ::Kokkos::deep_copy(device_gauge_data.view(), 0.0);

#ifdef KOKKOS_ENABLE_CUDA
  using device_vars_view_type = std::decay_t<decltype(device_vars.view())>;
  using device_dt_view_type = std::decay_t<decltype(device_dt.view())>;
  using device_deriv_view_type =
      std::decay_t<decltype(partial_derivatives_all_elements.view())>;
  static_assert(
      std::is_same_v<typename device_vars_view_type::array_layout,
                     Kokkos::LayoutLeft>,
      "GH batched volume RHS expects device_variables with LayoutLeft so "
      "points are fastest-varying.");
  static_assert(
      std::is_same_v<typename device_dt_view_type::array_layout,
                     Kokkos::LayoutLeft>,
      "GH batched volume RHS expects device_dt_variables with LayoutLeft so "
      "points are fastest-varying.");
  static_assert(std::is_same_v<typename device_deriv_view_type::array_layout,
                               Kokkos::LayoutLeft>,
                "GH batched volume RHS expects volume_partial_derivatives with "
                "LayoutLeft so points are fastest-varying.");
#endif

  const auto assert_points_fastest = [](const auto& view, const char* name) {
    std::array<size_t, 8> strides{};
    view.stride(strides.data());
    ASSERT(strides[0] == 1,
           name << " must have grid points as fastest-varying index "
                << "(stride[0] == 1), but got stride[0] = " << strides[0]
                << " with extents (" << view.extent(0) << ", " << view.extent(1)
                << ") and strides (" << strides[0] << ", " << strides[1]
                << ").");
  };
  assert_points_fastest(device_vars.view(), "device_variables.view()");
  assert_points_fastest(device_dt.view(), "device_dt_variables.view()");
  assert_points_fastest(partial_derivatives_all_elements.view(),
                        "volume_partial_derivatives.view()");

  const auto dt_spacetime_metric =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<pi_tag>>>(device_dt);
  const auto dt_phi = get<::Tags::MirrorView<::Tags::dt<phi_tag>>>(device_dt);
  const auto spacetime_metric = get<device_spacetime_metric_tag>(device_vars);
  const auto pi = get<device_pi_tag>(device_vars);
  const auto phi = get<device_phi_tag>(device_vars);
  const auto d_spacetime_metric =
      get<device_d_spacetime_metric_tag>(partial_derivatives_all_elements);
  const auto d_pi = get<device_d_pi_tag>(partial_derivatives_all_elements);
  const auto d_phi = get<device_d_phi_tag>(partial_derivatives_all_elements);
  const auto mesh_velocity = packed_geometry.mesh_velocity_device;
  const auto gauge_h = get<device_gauge_h_tag>(device_gauge_data);
  const auto spacetime_deriv_gauge_h =
      get<device_spacetime_deriv_gauge_h_tag>(device_gauge_data);
  auto& device_three_plus_one_data =
      packed_boundary_scratch->volume_three_plus_one_data;
  if (device_three_plus_one_data.number_of_grid_points() != total_points) {
    device_three_plus_one_data.initialize(total_points);
  }
  auto cached_inverse_spatial_metric =
      get<device_inverse_spatial_metric_tag>(device_three_plus_one_data);
  auto cached_shift = get<device_shift_tag>(device_three_plus_one_data);
  auto cached_lapse = get<device_lapse_tag>(device_three_plus_one_data);

  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeCacheThreePlusOne", total_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        const auto spacetime_metric_at_s = make_at_index(spacetime_metric, s);

        tnsr::II<double, dim, Frame::Inertial> inverse_spatial_metric{};
        tnsr::I<double, dim, Frame::Inertial> shift{};
        double lapse = 0.0;
        double inv_lapse_squared = 0.0;
        compute_inverse_spatial_metric_shift_lapse(
            make_not_null(&inverse_spatial_metric), make_not_null(&shift),
            make_not_null(&lapse), make_not_null(&inv_lapse_squared),
            spacetime_metric_at_s);
        (void)inv_lapse_squared;

        set_at_index(make_not_null(&cached_inverse_spatial_metric),
                     inverse_spatial_metric, point);
        set_at_index(make_not_null(&cached_shift), shift, point);
        get(cached_lapse)[point] = lapse;
      });

  Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtSpacetimeMetric",
      Kokkos::RangePolicy<Kokkos::Cuda, Kokkos::IndexType<int>>(0,
                                                                total_points),
      KOKKOS_LAMBDA(const int s) {
        const int point = s;

        const double gamma1 = get(make_at_index(device_constraint_gamma1, s));
        const double gamma1_plus_1 = 1.0 + gamma1;

        const auto shift_at_s = make_at_index(cached_shift, s);
        const auto mesh_velocity_at_s = make_at_index(mesh_velocity, s);
        const auto phi_at_s = make_at_index(phi, s);
        const auto pi_at_s = make_at_index(pi, s);
        const double lapse_at_s = get(make_at_index(cached_lapse, s));
        static_for<4>([&](auto mu_c) {
          constexpr int mu = (int)mu_c;
          static_for<4 - mu>([&](auto off_c) {
            constexpr int nu = mu + (int)off_c;
            double dt_piece = -lapse_at_s * get<mu, nu>(pi_at_s);
            double sum = 0.0;
            double mesh_sum = 0.0;
            double mesh_advection = 0.0;
            static_for<dim>([&](auto m_c) {
              constexpr int m = m_c;
              dt_piece += get<m>(shift_at_s) * get<m, mu, nu>(phi_at_s);
              sum += get<m>(shift_at_s) *
                     (get<m, mu, nu>(d_spacetime_metric)[point] -
                      get<m, mu, nu>(phi_at_s));
              mesh_sum += get<m>(mesh_velocity_at_s) *
                          (get<m, mu, nu>(d_spacetime_metric)[point] -
                           get<m, mu, nu>(phi_at_s));
              // Generic moving-mesh advection v^m d_m g_{mu nu}, added for
              // every evolved variable. Production adds this in
              // VolumeTermsImpl, separately from the GH-specific gamma1
              // mesh term (gamma1 * mesh_sum) above.
              mesh_advection += get<m>(mesh_velocity_at_s) *
                                get<m, mu, nu>(d_spacetime_metric)[point];
            });

            get<mu, nu>(dt_spacetime_metric)[point] =
                dt_piece + gamma1_plus_1 * sum + gamma1 * mesh_sum +
                mesh_advection;
          });
        });
      });

  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtPi", total_points,
      KOKKOS_LAMBDA(const int s) {
        constexpr int int_dim = 3;
        constexpr int int_dim_plus_one = 4;

        // -------------------------
        // Inputs
        // -------------------------
        const auto inverse_spatial_metric =
            make_at_index(cached_inverse_spatial_metric, s);
        const double lapse = get(make_at_index(cached_lapse, s));
        const auto shift = make_at_index(cached_shift, s);
        const auto mesh_velocity_at_s = make_at_index(mesh_velocity, s);
        const double gamma0 = get(make_at_index(device_constraint_gamma0, s));
        const double gamma1 = get(make_at_index(device_constraint_gamma1, s));
        const double gamma2 = get(make_at_index(device_constraint_gamma2, s));

        // -------------------------
        // Inverse spacetime metric
        // -------------------------
        double lapse_squared = -get<0, 0>(spacetime_metric)[s];
        static_for<int_dim>([&](auto i_c) {
          constexpr int i = (int)i_c;
          lapse_squared += get<i>(shift) * get<0, i + 1>(spacetime_metric)[s];
        });
        const double inv_lapse_squared = 1.0 / lapse_squared;

        tnsr::AA<double, 3, Frame::Inertial> inverse_spacetime_metric{};
        get<0, 0>(inverse_spacetime_metric) = -inv_lapse_squared;
        static_for<int_dim>([&](auto i_c) {
          constexpr int i = (int)i_c;
          get<0, i + 1>(inverse_spacetime_metric) =
              get<i>(shift) * inv_lapse_squared;
          static_for<int_dim - i>([&](auto off_c) {
            constexpr int j = i + (int)off_c;
            get<i + 1, j + 1>(inverse_spacetime_metric) =
                get<i, j>(inverse_spatial_metric) -
                get<i>(shift) * get<j>(shift) * inv_lapse_squared;
          });
        });

        // -------------------------
        // Normal^A
        // -------------------------
        tnsr::A<double, 3, Frame::Inertial> normal_spacetime_vector{};
        const double inv_lapse = 1.0 / lapse;
        get<0>(normal_spacetime_vector) = inv_lapse;
        static_for<int_dim>([&](auto i_c) {
          constexpr int i = (int)i_c;
          get<i + 1>(normal_spacetime_vector) = -get<i>(shift) * inv_lapse;
        });

        const double gamma1gamma2 = gamma1 * gamma2;
        const double minus_gamma0_lapse = -gamma0 * lapse;

        // =====================================================================
        // Phase B: gauge_constraint + shift_dot_three_index_constraint
        // =====================================================================
        tnsr::a<double, 3, Frame::Inertial> gauge_constraint{};
        tnsr::aa<double, 3, Frame::Inertial> shift_dot_three_index_constraint{};
        tnsr::aa<double, 3, Frame::Inertial>
            mesh_velocity_dot_three_index_constraint{};
        double normal_dot_gauge_constraint = 0.0;

        {
          tnsr::iaa<double, 3, Frame::Inertial> three_index_constraint{};
          static_for<int_dim>([&](auto n_c) {
            constexpr int n = (int)n_c;  // 0..2
            static_for<int_dim_plus_one>([&](auto mu_c) {
              constexpr int mu = (int)mu_c;  // 0..3
              static_for<int_dim_plus_one - mu>([&](auto off_c) {
                constexpr int nu = mu + (int)off_c;  // mu..3
                get<n, mu, nu>(three_index_constraint) =
                    d_spacetime_metric.get(n, mu, nu)[s] -
                    get<n, mu, nu>(phi)[s];
              });
            });
          });

          auto dt_g = [&](auto a_c, auto b_c) -> double {
            constexpr int a = (int)decltype(a_c)::value;  // 0..3
            constexpr int b = (int)decltype(b_c)::value;  // 0..3
            double v = -lapse * get<a, b>(pi)[s];
            static_for<int_dim>([&](auto m_c) {
              constexpr int m = (int)m_c;  // 0..2
              v += get<m>(shift) * get<m, a, b>(phi)[s];
            });
            return v;
          };

          auto d_g = [&](auto i_c, auto a_c, auto b_c) -> double {
            constexpr int i = (int)decltype(i_c)::value;  // 0..2
            constexpr int a = (int)decltype(a_c)::value;  // 0..3
            constexpr int b = (int)decltype(b_c)::value;  // 0..3
            return get<i, a, b>(phi)[s];
          };

          auto dslot_g = [&](auto slot_c, auto a_c, auto b_c) -> double {
            constexpr int slot = (int)decltype(slot_c)::value;  // 0..3
            if constexpr (slot == 0) {
              return dt_g(a_c, b_c);
            } else {
              return d_g(std::integral_constant<int, slot - 1>{}, a_c, b_c);
            }
          };

          tnsr::a<double, 3, Frame::Inertial> trace_christoffel{};
          static_for<int_dim_plus_one>([&](auto a_c) {
            double tr = 0.0;
            static_for<int_dim_plus_one>([&](auto b_c) {
              static_for<int_dim_plus_one>([&](auto c_c) {
                const double Gamma_abc =
                    0.5 * (dslot_g(b_c, c_c, a_c) + dslot_g(c_c, b_c, a_c) -
                           dslot_g(a_c, b_c, c_c));
                tr += Gamma_abc *
                      get<(int)b_c, (int)c_c>(inverse_spacetime_metric);
              });
            });
            get<(int)a_c>(trace_christoffel) = tr;
          });

          static_for<int_dim_plus_one>([&](auto mu_c) {
            constexpr int mu = (int)mu_c;

            get<mu>(gauge_constraint) =
                get<mu>(trace_christoffel) + get<mu>(gauge_h)[s];

            static_for<int_dim_plus_one - mu>([&](auto off_c) {
              constexpr int nu = mu + (int)off_c;
              double sd = 0.0;
              double md = 0.0;
              static_for<int_dim>([&](auto m_c) {
                constexpr int m = (int)m_c;
                sd += get<m>(shift) * get<m, mu, nu>(three_index_constraint);
                md += get<m>(mesh_velocity_at_s) *
                      get<m, mu, nu>(three_index_constraint);
              });
              get<mu, nu>(shift_dot_three_index_constraint) = sd;
              get<mu, nu>(mesh_velocity_dot_three_index_constraint) = md;
            });
          });

          static_for<int_dim_plus_one>([&](auto mu_c) {
            constexpr int mu = (int)mu_c;
            normal_dot_gauge_constraint +=
                get<mu>(normal_spacetime_vector) * get<mu>(gauge_constraint);
          });
        }

        normal_dot_gauge_constraint *= gamma0;

        // =====================================================================
        // Phase C: pi_one_normal and half_pi_two_normals
        // =====================================================================
        tnsr::a<double, 3, Frame::Inertial> pi_one_normal{};
        static_for<int_dim_plus_one>([&](auto mu_c) {
          constexpr int mu = (int)mu_c;
          double v = 0.0;
          static_for<int_dim_plus_one>([&](auto nu_c) {
            constexpr int nu = (int)nu_c;
            v += get<nu>(normal_spacetime_vector) * get<nu, mu>(pi)[s];
          });
          get<mu>(pi_one_normal) = v;
        });

        double half_pi_two_normals = 0.0;
        static_for<int_dim_plus_one>([&](auto mu_c) {
          constexpr int mu = (int)mu_c;
          half_pi_two_normals +=
              get<mu>(normal_spacetime_vector) * get<mu>(pi_one_normal);
        });
        half_pi_two_normals *= 0.5;

        // =====================================================================
        // Phase D: dt_pi prefill
        // =====================================================================
        {
          static_for<int_dim>([&](auto m_c) {
            constexpr int m = (int)m_c;
            constexpr int i = m + 1;
            get<0, i>(dt_pi)[s] =
                minus_gamma0_lapse * get<i>(gauge_constraint) -
                normal_dot_gauge_constraint * get<0, i>(spacetime_metric)[s];
          });
          get<0, 0>(dt_pi)[s] =
              2.0 * minus_gamma0_lapse * get<0>(gauge_constraint) -
              normal_dot_gauge_constraint * get<0, 0>(spacetime_metric)[s];

          static_for<int_dim>([&](auto mu0_c) {
            constexpr int mu = (int)mu0_c + 1;
            static_for<int_dim_plus_one - mu>([&](auto off_c) {
              constexpr int nu = mu + (int)off_c;
              get<mu, nu>(dt_pi)[s] = -normal_dot_gauge_constraint *
                                      get<mu, nu>(spacetime_metric)[s];
            });
          });
        }

        // =====================================================================
        // Phase E: streamed dt_pi accumulation per (mu,nu)
        //   Implements (1): REMOVE pi2_up/phi1_up/phi3_up arrays.
        // =====================================================================
        static_for<int_dim_plus_one>([&](auto mu_c) {
          constexpr int mu = (int)mu_c;

          static_for<int_dim_plus_one - mu>([&](auto off_c) {
            constexpr int nu = mu + (int)off_c;
            const auto nu_c = std::integral_constant<int, nu>{};

            double acc = get<mu, nu>(dt_pi)[s];

            acc -= half_pi_two_normals * get<mu, nu>(pi)[s];
            acc -= get<mu, nu>(spacetime_deriv_gauge_h)[s] +
                   get<nu, mu>(spacetime_deriv_gauge_h)[s];

            // ---- Local metric-derivative accessors (as in your Phase E) ----
            auto dt_g = [&](auto a_c, auto b_c) -> double {
              constexpr int a = (int)decltype(a_c)::value;
              constexpr int b = (int)decltype(b_c)::value;
              double v = -lapse * get<a, b>(pi)[s];
              static_for<int_dim>([&](auto m_c2) {
                constexpr int m2 = (int)m_c2;
                v += get<m2>(shift) * get<m2, a, b>(phi)[s];
              });
              return v;
            };
            auto d_g = [&](auto i_c, auto a_c, auto b_c) -> double {
              constexpr int i = (int)decltype(i_c)::value;
              constexpr int a = (int)decltype(a_c)::value;
              constexpr int b = (int)decltype(b_c)::value;
              return get<i, a, b>(phi)[s];
            };
            auto dslot_g = [&](auto slot_c, auto a_c, auto b_c) -> double {
              constexpr int slot = (int)decltype(slot_c)::value;
              if constexpr (slot == 0) {
                return dt_g(a_c, b_c);
              } else {
                return d_g(std::integral_constant<int, slot - 1>{}, a_c, b_c);
              }
            };

            // ---- On-demand contractions (NO ARRAYS) ----

            // pi2_up(delta) = inv_g^{delta beta} * pi_{nu beta}
            auto pi2_up = [&](auto delta_c2) -> double {
              constexpr int delta = (int)decltype(delta_c2)::value;
              double v = 0.0;
              static_for<int_dim_plus_one>([&](auto beta_c) {
                constexpr int beta = (int)beta_c;
                v += get<delta, beta>(inverse_spacetime_metric) *
                     get<nu, beta>(pi)[s];
              });
              return v;
            };

            // phi1_up(m,delta) = inv_spatial^{m n} * phi_{n mu delta}
            auto phi1_up = [&](auto m_c2, auto delta_c2) -> double {
              constexpr int m = (int)decltype(m_c2)::value;          // 0..2
              constexpr int delta = (int)decltype(delta_c2)::value;  // 0..3
              double v = 0.0;
              static_for<int_dim>([&](auto n_c2) {
                constexpr int n = (int)n_c2;  // 0..2
                v += get<m, n>(inverse_spatial_metric) *
                     get<n, mu, delta>(phi)[s];
              });
              return v;
            };

            // phi3_up(m,alpha) = inv_g^{alpha beta} * phi_{m nu beta}
            auto phi3_up = [&](auto m_c2, auto alpha_c2) -> double {
              constexpr int m = (int)decltype(m_c2)::value;          // 0..2
              constexpr int alpha = (int)decltype(alpha_c2)::value;  // 0..3
              double v = 0.0;
              static_for<int_dim_plus_one>([&](auto beta_c) {
                constexpr int beta = (int)beta_c;
                v += get<alpha, beta>(inverse_spacetime_metric) *
                     get<m, nu, beta>(phi)[s];
              });
              return v;
            };

            // Γ^delta_{mu nu} = inv_g^{delta alpha} * Γ_{alpha mu nu}
            auto Gamma2 = [&](auto delta_c2) -> double {
              constexpr int delta = (int)decltype(delta_c2)::value;
              double v = 0.0;
              static_for<int_dim_plus_one>([&](auto alpha_c) {
                const double Gamma1 = 0.5 * (dslot_g(mu_c, nu_c, alpha_c) +
                                             dslot_g(nu_c, mu_c, alpha_c) -
                                             dslot_g(alpha_c, mu_c, nu_c));
                v +=
                    get<delta, (int)alpha_c>(inverse_spacetime_metric) * Gamma1;
              });
              return v;
            };

            auto Gamma1 = [&](auto k_c, auto i_c, auto j_c) -> double {
              return 0.5 * (dslot_g(i_c, j_c, k_c) + dslot_g(j_c, i_c, k_c) -
                            dslot_g(k_c, i_c, j_c));
            };

            auto Gamma1_3up = [&](auto k_c, auto i_c, auto up_c) -> double {
              constexpr int up = (int)decltype(up_c)::value;
              double v = 0.0;
              static_for<int_dim_plus_one>([&](auto beta_c) {
                v += get<up, (int)beta_c>(inverse_spacetime_metric) *
                     Gamma1(k_c, i_c, beta_c);
              });
              return v;
            };

            // ---- Main delta loop ----
            static_for<int_dim_plus_one>([&](auto delta_c) {
              constexpr int delta = (int)delta_c;

              acc -= 2.0 * get<mu, delta>(pi)[s] * pi2_up(delta_c);
              acc += 2.0 * Gamma2(delta_c) * get<delta>(gauge_h)[s];

              // acc += 2 * sum_n phi1_up(n,delta) * phi3_up(n,delta)
              static_for<int_dim>([&](auto n_c) {
                acc += 2.0 * phi1_up(n_c, delta_c) * phi3_up(n_c, delta_c);
              });

              static_for<int_dim_plus_one>([&](auto alpha_c) {
                acc -= 2.0 * Gamma1_3up(mu_c, alpha_c, delta_c) *
                       Gamma1_3up(nu_c, delta_c, alpha_c);
              });
            });

            // Remaining m,n loops (replace phi1_up[m][nu] with on-demand)
            static_for<int_dim>([&](auto m_c) {
              constexpr int m = (int)m_c;
              constexpr int mp1 = m + 1;

              acc -= get<mp1>(pi_one_normal) * phi1_up(m_c, nu_c);

              static_for<int_dim>([&](auto n_c) {
                acc -= get<m, (int)n_c>(inverse_spatial_metric) *
                       get<m, (int)n_c, mu, nu>(d_phi)[s];
              });
            });

            // finalize
            acc *= lapse;
            acc += gamma1gamma2 * get<mu, nu>(shift_dot_three_index_constraint);
            acc += gamma1gamma2 *
                   get<mu, nu>(mesh_velocity_dot_three_index_constraint);

            static_for<int_dim>([&](auto m_c) {
              constexpr int m = (int)m_c;
              acc += get<m>(shift) * get<m, mu, nu>(d_pi)[s];
              // Generic moving-mesh advection v^m d_m Pi_{mu nu}.
              acc += get<m>(mesh_velocity_at_s) * get<m, mu, nu>(d_pi)[s];
            });

            get<mu, nu>(dt_pi)[s] = acc;
          });
        });
      });

  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtPhi", total_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        constexpr int int_dim = static_cast<int>(dim);
        constexpr int int_dim_plus_one = static_cast<int>(dim + 1);
        const double gamma2 = get(make_at_index(device_constraint_gamma2, s));

        const auto spacetime_metric_at_s = make_at_index(spacetime_metric, s);
        const auto pi_at_s = make_at_index(pi, s);
        const auto phi_at_s = make_at_index(phi, s);
        const auto mesh_velocity_at_s = make_at_index(mesh_velocity, s);

        const auto inverse_spatial_metric =
            make_at_index(cached_inverse_spatial_metric, s);
        const auto lapse = get(make_at_index(cached_lapse, s));
        const auto shift = make_at_index(cached_shift, s);

        tnsr::A<double, dim, Frame::Inertial> normal_spacetime_vector{};
        get<0>(normal_spacetime_vector) = 1.0 / lapse;
        static_for<int_dim>([&](auto i_c) {
          constexpr int i = (int)i_c;
          get<i + 1>(normal_spacetime_vector) = -get<i>(shift) / lapse;
        });

        static_for<int_dim>([&](auto i_c) {
          constexpr int i = (int)i_c;
          double phi_one_normal_i[4];
          static_for<int_dim_plus_one>([&](auto a_c) {
            constexpr int a = (int)a_c;
            phi_one_normal_i[a] = 0.0;
            static_for<int_dim_plus_one>([&](auto b_c) {
              constexpr int b = (int)b_c;
              phi_one_normal_i[a] +=
                  get<b>(normal_spacetime_vector) * get<i, b, a>(phi_at_s);
            });
          });
          double half_phi_two_normals_i = 0.0;

          static_for<int_dim_plus_one>([&](auto a_c) {
            constexpr int a = (int)a_c;
            half_phi_two_normals_i +=
                get<a>(normal_spacetime_vector) * phi_one_normal_i[a];
          });
          half_phi_two_normals_i *= 0.5;
          static_for<int_dim_plus_one>([&](auto mu_c) {
            constexpr int mu = (int)mu_c;
            static_for<int_dim_plus_one - mu>([&](auto offset_c) {
              constexpr int nu = mu + (int)offset_c;

              double dt_phi_i_mu_nu =
                  get<mu, nu>(pi_at_s) * half_phi_two_normals_i -
                  get<i, mu, nu>(d_pi)[point] +
                  gamma2 * (get<i, mu, nu>(d_spacetime_metric)[point] -
                            get<i, mu, nu>(phi_at_s));

              static_for<int_dim>([&](auto n_c) {
                constexpr int n = (int)n_c;
                double phi_1_up_n_mu_nu = 0.0;
                static_for<int_dim>([&](auto m_c) {
                  constexpr int m = (int)m_c;
                  phi_1_up_n_mu_nu += get<n, m>(inverse_spatial_metric) *
                                      get<m, mu, nu>(phi_at_s);
                });
                dt_phi_i_mu_nu += phi_one_normal_i[n + 1] * phi_1_up_n_mu_nu;
              });

              dt_phi_i_mu_nu *= lapse;
              static_for<int_dim>([&](auto m_c) {
                constexpr int m = (int)m_c;
                dt_phi_i_mu_nu +=
                    get<m>(shift) * get<m, i, mu, nu>(d_phi)[point];
                // Generic moving-mesh advection v^m d_m Phi_{i mu nu}.
                dt_phi_i_mu_nu += get<m>(mesh_velocity_at_s) *
                                  get<m, i, mu, nu>(d_phi)[point];
              });
              get<i, mu, nu>(dt_phi)[point] = dt_phi_i_mu_nu;
            });
          });
        });
      });
}

void ComputeTimeDerivativeBatched::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const gh::gauges::GaugeCondition& gauge_condition,
    const typename device_constraint_gamma0_tag::type& device_gamma0,
    const typename device_constraint_gamma1_tag::type& device_gamma1,
    const typename device_constraint_gamma2_tag::type& device_gamma2) {
  if (not gauge_condition.is_harmonic()) {
    ERROR(
        "The GH batched volume derivative currently supports only the "
        "Harmonic gauge.");
  }
  compute_time_derivative_batched_volume_impl(
      packed_evolution_state, packed_boundary_scratch, device_gamma0,
      device_gamma1, device_gamma2, packed_topology, packed_geometry);
}

}  // namespace gh::Actions
