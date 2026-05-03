// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BoundaryCorrectionHelpers.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

static constexpr size_t volume_dim = 3;
using system = gh::System<volume_dim>;
using packed_boundary_scratch_type =
    evolution::Kokkos::Tags::PackedBoundaryScratch<system>::type;
using device_package_field_tags =
    typename packed_boundary_scratch_type::device_package_field_tags;
using package_storage_type =
    typename Variables<device_package_field_tags>::storage_type;
using spacetime_metric_tag =
    gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;
using pi_tag = gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>;
using phi_tag = gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>;
using dt_variables_type =
    db::add_tag_prefix<::Tags::dt, typename system::variables_tag>::type;
using d_spacetime_metric_tag =
    ::Tags::deriv<spacetime_metric_tag, tmpl::size_t<volume_dim>,
                  Frame::Inertial>;
using d_pi_tag =
    ::Tags::deriv<pi_tag, tmpl::size_t<volume_dim>, Frame::Inertial>;
using d_phi_tag =
    ::Tags::deriv<phi_tag, tmpl::size_t<volume_dim>, Frame::Inertial>;
using host_derivative_tags =
    tmpl::list<d_spacetime_metric_tag, d_pi_tag, d_phi_tag>;
using gauge_h_tag = gh::Tags::GaugeH<DataVector, volume_dim>;
using spacetime_deriv_gauge_h_tag =
    gh::Tags::SpacetimeDerivGaugeH<DataVector, volume_dim>;
using host_gauge_tags = tmpl::list<gauge_h_tag, spacetime_deriv_gauge_h_tag>;
using inverse_spatial_metric_tag =
    gr::Tags::InverseSpatialMetric<DataVector, volume_dim, Frame::Inertial>;
using shift_tag = gr::Tags::Shift<DataVector, volume_dim, Frame::Inertial>;
using lapse_tag = gr::Tags::Lapse<DataVector>;
using host_three_plus_one_tags =
    tmpl::list<inverse_spatial_metric_tag, shift_tag, lapse_tag>;

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

double outward_sign(const Side side) {
  return side == Side::Upper ? 1.0 : -1.0;
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

KOKKOS_INLINE_FUNCTION void compute_spacetime_metric_from_3plus1(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        spacetime_metric,
    const Scalar<double>& lapse,
    const tnsr::I<double, volume_dim, Frame::Inertial>& shift,
    const tnsr::ii<double, volume_dim, Frame::Inertial>& spatial_metric) {
  spacetime_metric->get(0, 0) = -get(lapse) * get(lapse);
  for (size_t m = 0; m < volume_dim; ++m) {
    spacetime_metric->get(0, 0) +=
        spatial_metric.get(m, m) * shift.get(m) * shift.get(m);
    for (size_t n = 0; n < m; ++n) {
      spacetime_metric->get(0, 0) +=
          2.0 * spatial_metric.get(m, n) * shift.get(m) * shift.get(n);
    }
  }
  for (size_t i = 0; i < volume_dim; ++i) {
    spacetime_metric->get(0, i + 1) = 0.0;
    for (size_t m = 0; m < volume_dim; ++m) {
      spacetime_metric->get(0, i + 1) +=
          spatial_metric.get(m, i) * shift.get(m);
    }
    for (size_t j = i; j < volume_dim; ++j) {
      spacetime_metric->get(i + 1, j + 1) = spatial_metric.get(i, j);
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_phi_from_3plus1(
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*> phi,
    const Scalar<double>& lapse,
    const tnsr::i<double, volume_dim, Frame::Inertial>& deriv_lapse,
    const tnsr::I<double, volume_dim, Frame::Inertial>& shift,
    const tnsr::iJ<double, volume_dim, Frame::Inertial>& deriv_shift,
    const tnsr::ii<double, volume_dim, Frame::Inertial>& spatial_metric,
    const tnsr::ijj<double, volume_dim, Frame::Inertial>&
        deriv_spatial_metric) {
  for (size_t k = 0; k < volume_dim; ++k) {
    phi->get(k, 0, 0) = -2.0 * get(lapse) * deriv_lapse.get(k);
    for (size_t m = 0; m < volume_dim; ++m) {
      for (size_t n = 0; n < volume_dim; ++n) {
        phi->get(k, 0, 0) +=
            deriv_spatial_metric.get(k, m, n) * shift.get(m) * shift.get(n) +
            2.0 * spatial_metric.get(m, n) * shift.get(m) *
                deriv_shift.get(k, n);
      }
    }
    for (size_t i = 0; i < volume_dim; ++i) {
      phi->get(k, 0, i + 1) = 0.0;
      for (size_t m = 0; m < volume_dim; ++m) {
        phi->get(k, 0, i + 1) +=
            deriv_spatial_metric.get(k, m, i) * shift.get(m) +
            spatial_metric.get(m, i) * deriv_shift.get(k, m);
      }
      for (size_t j = i; j < volume_dim; ++j) {
        phi->get(k, i + 1, j + 1) = deriv_spatial_metric.get(k, i, j);
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_pi_from_3plus1(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*> pi,
    const Scalar<double>& lapse, const Scalar<double>& dt_lapse,
    const tnsr::I<double, volume_dim, Frame::Inertial>& shift,
    const tnsr::I<double, volume_dim, Frame::Inertial>& dt_shift,
    const tnsr::ii<double, volume_dim, Frame::Inertial>& spatial_metric,
    const tnsr::ii<double, volume_dim, Frame::Inertial>& dt_spatial_metric,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>& phi) {
  pi->get(0, 0) = -2.0 * get(lapse) * get(dt_lapse);
  for (size_t m = 0; m < volume_dim; ++m) {
    for (size_t n = 0; n < volume_dim; ++n) {
      pi->get(0, 0) +=
          dt_spatial_metric.get(m, n) * shift.get(m) * shift.get(n) +
          2.0 * spatial_metric.get(m, n) * shift.get(m) * dt_shift.get(n);
    }
  }
  for (size_t i = 0; i < volume_dim; ++i) {
    pi->get(0, i + 1) = 0.0;
    for (size_t m = 0; m < volume_dim; ++m) {
      pi->get(0, i + 1) += dt_spatial_metric.get(m, i) * shift.get(m) +
                           spatial_metric.get(m, i) * dt_shift.get(m);
    }
    for (size_t j = i; j < volume_dim; ++j) {
      pi->get(i + 1, j + 1) = dt_spatial_metric.get(i, j);
    }
  }
  for (size_t mu = 0; mu < volume_dim + 1; ++mu) {
    for (size_t nu = mu; nu < volume_dim + 1; ++nu) {
      for (size_t i = 0; i < volume_dim; ++i) {
        pi->get(mu, nu) -= shift.get(i) * phi.get(i, mu, nu);
      }
      pi->get(mu, nu) /= -get(lapse);
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_hardcoded_schwarzschild_3plus1(
    const gsl::not_null<Scalar<double>*> lapse,
    const gsl::not_null<Scalar<double>*> dt_lapse,
    const gsl::not_null<tnsr::i<double, volume_dim, Frame::Inertial>*>
        deriv_lapse,
    const gsl::not_null<tnsr::I<double, volume_dim, Frame::Inertial>*> shift,
    const gsl::not_null<tnsr::I<double, volume_dim, Frame::Inertial>*> dt_shift,
    const gsl::not_null<tnsr::iJ<double, volume_dim, Frame::Inertial>*>
        deriv_shift,
    const gsl::not_null<tnsr::ii<double, volume_dim, Frame::Inertial>*>
        spatial_metric,
    const gsl::not_null<tnsr::ii<double, volume_dim, Frame::Inertial>*>
        dt_spatial_metric,
    const gsl::not_null<tnsr::ijj<double, volume_dim, Frame::Inertial>*>
        deriv_spatial_metric,
    const tnsr::I<double, volume_dim, Frame::Inertial>& inertial_coordinates) {
  constexpr double mass = 1.0;
  const double x = inertial_coordinates.get(0);
  const double y = inertial_coordinates.get(1);
  const double z = inertial_coordinates.get(2);
  const double r_squared = x * x + y * y + z * z;
  const double r = sqrt(r_squared);
  const double inv_r = 1.0 / r;
  const double inv_r_squared = 1.0 / r_squared;
  const double inv_r_cubed = inv_r * inv_r_squared;
  const double inv_r_fourth = inv_r_squared * inv_r_squared;
  const double inv_r_fifth = inv_r_squared * inv_r_cubed;
  const double two_m_over_r = 2.0 * mass * inv_r;

  const double lapse_squared = 1.0 / (1.0 + two_m_over_r);
  get(*lapse) = sqrt(lapse_squared);
  get(*dt_lapse) = 0.0;

  const std::array<double, volume_dim> coords{{x, y, z}};
  for (size_t i = 0; i < volume_dim; ++i) {
    const double xi = coords[i];
    deriv_lapse->get(i) = mass * xi * inv_r_cubed * get(*lapse) * lapse_squared;
    shift->get(i) = 2.0 * mass * xi * inv_r_squared * lapse_squared;
    dt_shift->get(i) = 0.0;
  }

  for (size_t i = 0; i < volume_dim; ++i) {
    for (size_t j = i; j < volume_dim; ++j) {
      const double delta_ij = i == j ? 1.0 : 0.0;
      spatial_metric->get(i, j) =
          delta_ij + 2.0 * mass * coords[i] * coords[j] * inv_r_cubed;
      dt_spatial_metric->get(i, j) = 0.0;
    }
  }

  for (size_t k = 0; k < volume_dim; ++k) {
    for (size_t i = 0; i < volume_dim; ++i) {
      const double delta_ki = k == i ? 1.0 : 0.0;
      deriv_shift->get(k, i) =
          2.0 * mass * inv_r_squared * lapse_squared * delta_ki -
          4.0 * mass * coords[k] * coords[i] * inv_r_fourth * lapse_squared *
              (1.0 - mass * inv_r * lapse_squared);
      for (size_t j = i; j < volume_dim; ++j) {
        const double delta_kj = k == j ? 1.0 : 0.0;
        deriv_spatial_metric->get(k, i, j) =
            -6.0 * mass * coords[i] * coords[j] * coords[k] * inv_r_fifth +
            2.0 * mass * coords[i] * delta_kj * inv_r_cubed +
            2.0 * mass * coords[j] * delta_ki * inv_r_cubed;
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_hardcoded_schwarzschild_gh_fields(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        spacetime_metric,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*> pi,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*> phi,
    const tnsr::I<double, volume_dim, Frame::Inertial>& inertial_coordinates) {
  Scalar<double> lapse{};
  Scalar<double> dt_lapse{};
  tnsr::i<double, volume_dim, Frame::Inertial> deriv_lapse{};
  tnsr::I<double, volume_dim, Frame::Inertial> shift{};
  tnsr::I<double, volume_dim, Frame::Inertial> dt_shift{};
  tnsr::iJ<double, volume_dim, Frame::Inertial> deriv_shift{};
  tnsr::ii<double, volume_dim, Frame::Inertial> spatial_metric{};
  tnsr::ii<double, volume_dim, Frame::Inertial> dt_spatial_metric{};
  tnsr::ijj<double, volume_dim, Frame::Inertial> deriv_spatial_metric{};

  compute_hardcoded_schwarzschild_3plus1(
      make_not_null(&lapse), make_not_null(&dt_lapse),
      make_not_null(&deriv_lapse), make_not_null(&shift),
      make_not_null(&dt_shift), make_not_null(&deriv_shift),
      make_not_null(&spatial_metric), make_not_null(&dt_spatial_metric),
      make_not_null(&deriv_spatial_metric), inertial_coordinates);

  compute_phi_from_3plus1(phi, lapse, deriv_lapse, shift, deriv_shift,
                          spatial_metric, deriv_spatial_metric);
  compute_pi_from_3plus1(pi, lapse, dt_lapse, shift, dt_shift, spatial_metric,
                         dt_spatial_metric, *phi);
  compute_spacetime_metric_from_3plus1(spacetime_metric, lapse, shift,
                                       spatial_metric);
}

KOKKOS_INLINE_FUNCTION void load_aa_from_packaged_data(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  static_for<volume_dim + 1>([&](auto a_c) {
    constexpr int a = (int)a_c;
    static_for<volume_dim + 1 - a>([&](auto off_c) {
      constexpr int b = a + (int)off_c;
      get<a, b>(*tensor) =
          packaged_data_view(point, component_offset + component_index);
      ++component_index;
    });
  });
}

KOKKOS_INLINE_FUNCTION void load_iaa_from_packaged_data(
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  size_t component_index = 0;
  static_for<volume_dim>([&](auto d_c) {
    constexpr int d = (int)d_c;
    static_for<volume_dim + 1>([&](auto a_c) {
      constexpr int a = (int)a_c;
      static_for<volume_dim + 1 - a>([&](auto off_c) {
        constexpr int b = a + (int)off_c;
        get<d, a, b>(*tensor) =
            packaged_data_view(point, component_offset + component_index);
        ++component_index;
      });
    });
  });
}

KOKKOS_INLINE_FUNCTION void load_a_from_packaged_data(
    const gsl::not_null<tnsr::a<double, volume_dim, Frame::Inertial>*> tensor,
    const package_storage_type& packaged_data_view, const size_t point,
    const size_t component_offset) {
  static_for<volume_dim + 1>([&](auto a_c) {
    constexpr int a = (int)a_c;
    get<a>(*tensor) = packaged_data_view(point, component_offset + a);
  });
}

}  // namespace

void ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const typename packed_boundary_metadata_tag::type& packed_boundary_metadata,
    const typename device_constraint_gamma0_tag::type& /*device_gamma0*/,
    const typename device_constraint_gamma1_tag::type& device_gamma1,
    const typename device_constraint_gamma2_tag::type& device_gamma2,
    const typename packed_boundary_scratch_tag::type& packed_boundary_scratch) {
  if (packed_topology.total_points == 0) {
    return;
  }

  ASSERT(packed_topology.uniform_quadrature_host[0] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[1] ==
                 Spectral::Quadrature::GaussLobatto and
             packed_topology.uniform_quadrature_host[2] ==
                 Spectral::Quadrature::GaussLobatto,
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched currently "
         "supports Gauss-Lobatto quadrature only.");

  const auto& extents = packed_topology.uniform_extents_host;
  ASSERT(extents[0] == extents[1] and extents[1] == extents[2],
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched assumes "
         "uniform p across dimensions (equal extents in all dimensions).");
  const Mesh<volume_dim> mesh{extents, packed_topology.uniform_basis_host,
                              packed_topology.uniform_quadrature_host};
  const Mesh<volume_dim - 1> uniform_face_mesh = mesh.slice_away(0);
  const size_t num_face_points = uniform_face_mesh.number_of_grid_points();
  const size_t total_points = packed_topology.total_points;

  const auto& elements = packed_topology.local_elements;
  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  const auto inverse_jacobian_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, inverse_jacobian);
  const auto& element_point_offsets_device =
      packed_topology.element_point_offsets_device;
  const auto& inertial_coordinates =
      packed_geometry.inertial_coordinates_device;
  const auto& mesh_velocity = packed_geometry.mesh_velocity_device;
  const auto& external_face_mask_for_all_elements =
      packed_boundary_metadata.external_face_mask_for_all_elements;
  ASSERT(mesh_velocity.get(0).extent(0) == packed_topology.total_points and
             mesh_velocity.get(1).extent(0) == packed_topology.total_points and
             mesh_velocity.get(2).extent(0) == packed_topology.total_points,
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched requires "
         "packed mesh velocity with one entry per packed volume point.");

  typename system::variables_tag::type host_vars{total_points};
  dt_variables_type host_dt_vars{total_points};
  Variables<host_derivative_tags> host_derivatives{total_points};
  Variables<host_gauge_tags> host_gauge_data{total_points};
  Variables<host_three_plus_one_tags> host_three_plus_one_data{total_points};
  copy_to_host(make_not_null(&host_vars),
               packed_evolution_state->device_variables);
  copy_to_host(make_not_null(&host_dt_vars),
               packed_evolution_state->device_dt_variables);
  copy_to_host(make_not_null(&host_derivatives),
               packed_boundary_scratch.volume_partial_derivatives);
  copy_to_host(make_not_null(&host_gauge_data),
               packed_boundary_scratch.volume_gauge_data);
  copy_to_host(make_not_null(&host_three_plus_one_data),
               packed_boundary_scratch.volume_three_plus_one_data);

  DataVector gamma1_host{total_points, 0.0};
  DataVector gamma2_host{total_points, 0.0};
  {
    const auto gamma1_device = device_gamma1.get();
    const auto gamma2_device = device_gamma2.get();
    const auto gamma1_mirror = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, gamma1_device);
    const auto gamma2_mirror = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, gamma2_device);
    for (size_t p = 0; p < total_points; ++p) {
      gamma1_host[p] = gamma1_mirror(p);
      gamma2_host[p] = gamma2_mirror(p);
    }
  }

  const auto& spacetime_metric = get<spacetime_metric_tag>(host_vars);
  const auto& pi = get<pi_tag>(host_vars);
  const auto& phi = get<phi_tag>(host_vars);
  auto& dt_spacetime_metric =
      get<::Tags::dt<spacetime_metric_tag>>(host_dt_vars);
  auto& dt_pi = get<::Tags::dt<pi_tag>>(host_dt_vars);
  auto& dt_phi = get<::Tags::dt<phi_tag>>(host_dt_vars);
  const auto& d_spacetime_metric =
      get<d_spacetime_metric_tag>(host_derivatives);
  const auto& d_pi = get<d_pi_tag>(host_derivatives);
  const auto& d_phi = get<d_phi_tag>(host_derivatives);
  const auto& gauge_h = get<gauge_h_tag>(host_gauge_data);
  const auto& spacetime_deriv_gauge_h =
      get<spacetime_deriv_gauge_h_tag>(host_gauge_data);
  const auto& inverse_spatial_metric =
      get<inverse_spatial_metric_tag>(host_three_plus_one_data);
  const auto& shift = get<shift_tag>(host_three_plus_one_data);
  const auto& lapse = get<lapse_tag>(host_three_plus_one_data);

  const auto fill_aa_face =
      [](const gsl::not_null<tnsr::aa<DataVector, volume_dim, Frame::Inertial>*>
             face_tensor,
         const tnsr::aa<DataVector, volume_dim, Frame::Inertial>& volume_tensor,
         const std::vector<size_t>& volume_indices) {
        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = a; b < volume_dim + 1; ++b) {
            for (size_t p = 0; p < volume_indices.size(); ++p) {
              face_tensor->get(a, b)[p] =
                  volume_tensor.get(a, b)[volume_indices[p]];
            }
          }
        }
      };
  const auto fill_iaa_face =
      [](const gsl::not_null<
             tnsr::iaa<DataVector, volume_dim, Frame::Inertial>*>
             face_tensor,
         const tnsr::iaa<DataVector, volume_dim, Frame::Inertial>&
             volume_tensor,
         const std::vector<size_t>& volume_indices) {
        for (size_t i = 0; i < volume_dim; ++i) {
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              for (size_t p = 0; p < volume_indices.size(); ++p) {
                face_tensor->get(i, a, b)[p] =
                    volume_tensor.get(i, a, b)[volume_indices[p]];
              }
            }
          }
        }
      };
  const auto fill_ijaa_face =
      [](const gsl::not_null<
             tnsr::ijaa<DataVector, volume_dim, Frame::Inertial>*>
             face_tensor,
         const tnsr::ijaa<DataVector, volume_dim, Frame::Inertial>&
             volume_tensor,
         const std::vector<size_t>& volume_indices) {
        for (size_t i = 0; i < volume_dim; ++i) {
          for (size_t j = 0; j < volume_dim; ++j) {
            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                for (size_t p = 0; p < volume_indices.size(); ++p) {
                  face_tensor->get(i, j, a, b)[p] =
                      volume_tensor.get(i, j, a, b)[volume_indices[p]];
                }
              }
            }
          }
        }
      };

  const size_t num_elements = elements.size();
  const gh::BoundaryConditions::ConstraintPreservingBjorhus<volume_dim>
      cp_physical_boundary_condition{
          gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType::
              ConstraintPreservingPhysical};

  for (size_t local_face_id = 0; local_face_id < 2 * volume_dim;
       ++local_face_id) {
    const size_t sliced_dim = local_face_id / 2;
    const size_t side_i = local_face_id % 2;
    const Side side = side_i == 0 ? Side::Lower : Side::Upper;
    const auto external_face_mask = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{},
        external_face_mask_for_all_elements[local_face_id]);
    ASSERT(static_cast<size_t>(external_face_mask.extent(0)) == num_elements,
           "External-face mask extent mismatch. "
           "InitializeBoundaryBatchMetadata must run before external "
           "boundary corrections.");

    const auto local_face_to_volume_index_device =
        side == Side::Upper
            ? gsl::at(face_to_volume_index_map, sliced_dim).second
            : gsl::at(face_to_volume_index_map, sliced_dim).first;
    const auto local_face_to_volume_index =
        ::Kokkos::create_mirror_view_and_copy(
            ::Kokkos::HostSpace{}, local_face_to_volume_index_device);
    const double outward = outward_sign(side);
    const double lift_prefactor =
        -0.5 *
        static_cast<double>(extents[sliced_dim] * (extents[sliced_dim] - 1));

    for (size_t element_index = 0; element_index < num_elements;
         ++element_index) {
      if (external_face_mask(element_index) == 0) {
        continue;
      }
      std::vector<size_t> volume_indices(num_face_points);
      std::vector<size_t> local_volume_indices(num_face_points);
      for (size_t p = 0; p < num_face_points; ++p) {
        local_volume_indices[p] = local_face_to_volume_index(p);
        volume_indices[p] =
            packed_topology.element_point_offsets_host[element_index] +
            local_volume_indices[p];
      }

      tnsr::aa<DataVector, volume_dim, Frame::Inertial> face_spacetime_metric{
          num_face_points};
      tnsr::aa<DataVector, volume_dim, Frame::Inertial> face_pi{
          num_face_points};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial> face_phi{
          num_face_points};
      tnsr::aa<DataVector, volume_dim, Frame::Inertial>
          face_dt_spacetime_metric{num_face_points};
      tnsr::aa<DataVector, volume_dim, Frame::Inertial> face_dt_pi{
          num_face_points};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial> face_dt_phi{
          num_face_points};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial>
          face_d_spacetime_metric{num_face_points};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial> face_d_pi{
          num_face_points};
      tnsr::ijaa<DataVector, volume_dim, Frame::Inertial> face_d_phi{
          num_face_points};
      tnsr::I<DataVector, volume_dim, Frame::Inertial> face_coords{
          num_face_points};
      tnsr::I<DataVector, volume_dim, Frame::Inertial> face_mesh_velocity{
          num_face_points};
      tnsr::i<DataVector, volume_dim, Frame::Inertial> normal_covector{
          num_face_points};
      tnsr::I<DataVector, volume_dim, Frame::Inertial> normal_vector{
          num_face_points};
      tnsr::AA<DataVector, volume_dim, Frame::Inertial>
          face_inverse_spacetime_metric{num_face_points};
      tnsr::A<DataVector, volume_dim, Frame::Inertial>
          face_spacetime_unit_normal_vector{num_face_points};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial> three_index_constraint{
          num_face_points};
      tnsr::a<DataVector, volume_dim, Frame::Inertial> face_gauge_h{
          num_face_points};
      tnsr::ab<DataVector, volume_dim, Frame::Inertial>
          face_spacetime_deriv_gauge_h{num_face_points};
      Scalar<DataVector> face_gamma1{num_face_points};
      Scalar<DataVector> face_gamma2{num_face_points};
      Scalar<DataVector> face_lapse{num_face_points};
      tnsr::I<DataVector, volume_dim, Frame::Inertial> face_shift{
          num_face_points};
      tnsr::II<DataVector, volume_dim, Frame::Inertial>
          face_inverse_spatial_metric{num_face_points};
      DataVector normal_magnitude{num_face_points, 0.0};

      fill_aa_face(make_not_null(&face_spacetime_metric), spacetime_metric,
                   volume_indices);
      fill_aa_face(make_not_null(&face_pi), pi, volume_indices);
      fill_iaa_face(make_not_null(&face_phi), phi, volume_indices);
      fill_aa_face(make_not_null(&face_dt_spacetime_metric),
                   dt_spacetime_metric, volume_indices);
      fill_aa_face(make_not_null(&face_dt_pi), dt_pi, volume_indices);
      fill_iaa_face(make_not_null(&face_dt_phi), dt_phi, volume_indices);
      fill_iaa_face(make_not_null(&face_d_spacetime_metric), d_spacetime_metric,
                    volume_indices);
      fill_iaa_face(make_not_null(&face_d_pi), d_pi, volume_indices);
      fill_ijaa_face(make_not_null(&face_d_phi), d_phi, volume_indices);

      for (size_t p = 0; p < num_face_points; ++p) {
        const size_t volume_index = volume_indices[p];
        get(face_gamma1)[p] = gamma1_host[volume_index];
        get(face_gamma2)[p] = gamma2_host[volume_index];
        get(face_lapse)[p] = get(lapse)[volume_index];
        for (size_t d = 0; d < volume_dim; ++d) {
          face_coords.get(d)[p] =
              packed_geometry.inertial_coordinates_host[d][volume_index];
          face_mesh_velocity.get(d)[p] =
              packed_geometry.mesh_velocity_host[d][volume_index];
          face_shift.get(d)[p] = shift.get(d)[volume_index];
          face_gauge_h.get(d + 1)[p] = gauge_h.get(d + 1)[volume_index];
          for (size_t j = d; j < volume_dim; ++j) {
            face_inverse_spatial_metric.get(d, j)[p] =
                inverse_spatial_metric.get(d, j)[volume_index];
          }
        }
        face_gauge_h.get(0)[p] = gauge_h.get(0)[volume_index];
        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = 0; b < volume_dim + 1; ++b) {
            face_spacetime_deriv_gauge_h.get(a, b)[p] =
                spacetime_deriv_gauge_h.get(a, b)[volume_index];
          }
        }

        for (size_t d = 0; d < volume_dim; ++d) {
          normal_covector.get(d)[p] =
              outward * inverse_jacobian_host(element_index,
                                              local_volume_indices[p],
                                              sliced_dim * volume_dim + d);
        }
        double normal_magnitude_squared = 0.0;
        for (size_t i = 0; i < volume_dim; ++i) {
          double component = 0.0;
          for (size_t j = 0; j < volume_dim; ++j) {
            component += face_inverse_spatial_metric.get(i, j)[p] *
                         normal_covector.get(j)[p];
          }
          normal_vector.get(i)[p] = component;
          normal_magnitude_squared += component * normal_covector.get(i)[p];
        }
        normal_magnitude[p] = sqrt(normal_magnitude_squared);
        for (size_t i = 0; i < volume_dim; ++i) {
          normal_covector.get(i)[p] /= normal_magnitude[p];
          normal_vector.get(i)[p] /= normal_magnitude[p];
        }

        const double lapse_at_point = get(face_lapse)[p];
        const double inv_lapse_squared =
            1.0 / (lapse_at_point * lapse_at_point);
        face_inverse_spacetime_metric.get(0, 0)[p] = -inv_lapse_squared;
        face_spacetime_unit_normal_vector.get(0)[p] = 1.0 / lapse_at_point;
        for (size_t i = 0; i < volume_dim; ++i) {
          face_inverse_spacetime_metric.get(0, i + 1)[p] =
              face_shift.get(i)[p] * inv_lapse_squared;
          face_spacetime_unit_normal_vector.get(i + 1)[p] =
              -face_shift.get(i)[p] / lapse_at_point;
          for (size_t j = i; j < volume_dim; ++j) {
            face_inverse_spacetime_metric.get(i + 1, j + 1)[p] =
                face_inverse_spatial_metric.get(i, j)[p] -
                face_shift.get(i)[p] * face_shift.get(j)[p] * inv_lapse_squared;
          }
        }
        for (size_t i = 0; i < volume_dim; ++i) {
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              three_index_constraint.get(i, a, b)[p] =
                  face_d_spacetime_metric.get(i, a, b)[p] -
                  face_phi.get(i, a, b)[p];
            }
          }
        }
      }

      tnsr::aa<DataVector, volume_dim, Frame::Inertial>
          dt_spacetime_metric_correction{num_face_points, 0.0};
      tnsr::aa<DataVector, volume_dim, Frame::Inertial> dt_pi_correction{
          num_face_points, 0.0};
      tnsr::iaa<DataVector, volume_dim, Frame::Inertial> dt_phi_correction{
          num_face_points, 0.0};
      const auto error = cp_physical_boundary_condition.dg_time_derivative(
          make_not_null(&dt_spacetime_metric_correction),
          make_not_null(&dt_pi_correction), make_not_null(&dt_phi_correction),
          std::optional<tnsr::I<DataVector, volume_dim, Frame::Inertial>>{
              face_mesh_velocity},
          normal_covector, normal_vector, face_spacetime_metric, face_pi,
          face_phi, face_coords, face_gamma1, face_gamma2, face_lapse,
          face_shift, face_inverse_spacetime_metric,
          face_spacetime_unit_normal_vector, three_index_constraint,
          face_gauge_h, face_spacetime_deriv_gauge_h, face_dt_spacetime_metric,
          face_dt_pi, face_dt_phi, face_d_spacetime_metric, face_d_pi,
          face_d_phi);
      ASSERT(not error.has_value(),
             "Hard-coded CP Physical GH batched outer boundary failed: "
                 << error.value());

      for (size_t p = 0; p < num_face_points; ++p) {
        const size_t volume_index = volume_indices[p];
        const double lifted_factor = lift_prefactor * normal_magnitude[p];
        for (size_t a = 0; a < volume_dim + 1; ++a) {
          for (size_t b = a; b < volume_dim + 1; ++b) {
            dt_spacetime_metric.get(a, b)[volume_index] +=
                lifted_factor * dt_spacetime_metric_correction.get(a, b)[p];
            dt_pi.get(a, b)[volume_index] +=
                lifted_factor * dt_pi_correction.get(a, b)[p];
            for (size_t d = 0; d < volume_dim; ++d) {
              dt_phi.get(d, a, b)[volume_index] +=
                  lifted_factor * dt_phi_correction.get(d, a, b)[p];
            }
          }
        }
      }
    }
  }
  packed_evolution_state->device_dt_variables = copy_to_device(host_dt_vars);
}

}  // namespace gh::Actions
