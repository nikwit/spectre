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
#include "Utilities/ErrorHandling/Error.hpp"
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

using evolved_variables_type = typename system::variables_tag::type;
using host_derivatives_type = Variables<host_derivative_tags>;
using host_gauge_data_type = Variables<host_gauge_tags>;
using host_three_plus_one_data_type = Variables<host_three_plus_one_tags>;

constexpr size_t num_evolved_components =
    evolved_variables_type::number_of_independent_components;
constexpr size_t num_dt_components =
    dt_variables_type::number_of_independent_components;
constexpr size_t num_derivative_components =
    host_derivatives_type::number_of_independent_components;
constexpr size_t num_gauge_components =
    host_gauge_data_type::number_of_independent_components;
constexpr size_t num_three_plus_one_components =
    host_three_plus_one_data_type::number_of_independent_components;
constexpr size_t bjorhus_offset_evolved = 0;
constexpr size_t bjorhus_offset_dt =
    bjorhus_offset_evolved + num_evolved_components;
constexpr size_t bjorhus_offset_derivatives =
    bjorhus_offset_dt + num_dt_components;
constexpr size_t bjorhus_offset_gauge =
    bjorhus_offset_derivatives + num_derivative_components;
constexpr size_t bjorhus_offset_three_plus_one =
    bjorhus_offset_gauge + num_gauge_components;
constexpr size_t bjorhus_offset_gamma1 =
    bjorhus_offset_three_plus_one + num_three_plus_one_components;
constexpr size_t bjorhus_offset_gamma2 = bjorhus_offset_gamma1 + 1;
constexpr size_t bjorhus_offset_coordinates = bjorhus_offset_gamma2 + 1;
constexpr size_t bjorhus_offset_mesh_velocity =
    bjorhus_offset_coordinates + volume_dim;
constexpr size_t bjorhus_offset_unnormalized_normal =
    bjorhus_offset_mesh_velocity + volume_dim;
constexpr size_t num_bjorhus_input_components =
    bjorhus_offset_unnormalized_normal + volume_dim;

template <typename VariablesType, typename HostView>
void load_face_variables(const gsl::not_null<VariablesType*> face_variables,
                         const HostView& host_input,
                         const size_t face_point_offset,
                         const size_t input_component_offset,
                         const size_t num_face_points) {
  for (size_t component = 0;
       component < VariablesType::number_of_independent_components;
       ++component) {
    for (size_t p = 0; p < num_face_points; ++p) {
      face_variables->data()[component * num_face_points + p] =
          host_input(face_point_offset + p, input_component_offset + component);
    }
  }
}

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
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const typename packed_boundary_metadata_tag::type& packed_boundary_metadata,
    const typename device_constraint_gamma0_tag::type& /*device_gamma0*/,
    const typename device_constraint_gamma1_tag::type& device_gamma1,
    const typename device_constraint_gamma2_tag::type& device_gamma2) {
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

  const auto& elements = packed_topology.local_elements;
  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  const auto& mesh_velocity = packed_geometry.mesh_velocity_device;
  const auto& external_boundary_condition_kind_for_all_elements =
      packed_boundary_metadata
          .external_boundary_condition_kind_for_all_elements;
  ASSERT(mesh_velocity.get(0).extent(0) == packed_topology.total_points and
             mesh_velocity.get(1).extent(0) == packed_topology.total_points and
             mesh_velocity.get(2).extent(0) == packed_topology.total_points,
         "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched requires "
         "packed mesh velocity with one entry per packed volume point.");

  const size_t num_elements = elements.size();
  const auto& external_face_masks_host =
      packed_boundary_metadata.external_face_mask_for_all_elements_host;
  const auto& external_boundary_condition_kinds_host =
      packed_boundary_metadata
          .external_boundary_condition_kind_for_all_elements_host;
  bool has_constraint_preserving_physical_face = false;
  for (size_t local_face_id = 0; local_face_id < 2 * volume_dim;
       ++local_face_id) {
    const auto& host_mask = external_face_masks_host[local_face_id];
    const auto& host_kinds =
        external_boundary_condition_kinds_host[local_face_id];
    ASSERT(host_mask.size() == num_elements,
           "External-face mask extent mismatch. "
           "InitializeBoundaryBatchMetadata must run before external "
           "boundary corrections.");
    ASSERT(host_kinds.size() == num_elements,
           "External-boundary-condition kind extent mismatch. "
           "InitializeBoundaryBatchMetadata must run before external "
           "boundary corrections.");
    bool face_has_demand_outgoing = false;
    for (size_t element_index = 0; element_index < num_elements;
         ++element_index) {
      const auto kind = static_cast<
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind>(
          host_kinds[element_index]);
      if (host_mask[element_index] == 0) {
        if (kind !=
            evolution::Kokkos::Batched::ExternalBoundaryConditionKind::None) {
          ERROR("Nonexternal face "
                << local_face_id << " of element "
                << elements[element_index].id()
                << " has an encoded external boundary-condition kind.");
        }
        continue;
      }
      if (kind == evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
                      DemandOutgoingCharSpeeds) {
        face_has_demand_outgoing = true;
      } else if (kind ==
                 evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
                     ConstraintPreservingPhysical) {
        has_constraint_preserving_physical_face = true;
      } else {
        ERROR("Unsupported encoded external boundary-condition kind "
              << static_cast<int>(kind) << " on face " << local_face_id
              << " of element " << elements[element_index].id() << ".");
      }
    }

    if (face_has_demand_outgoing) {
      const auto& packaged_face_data =
          packed_boundary_scratch
              ->packaged_face_data_for_all_elements[local_face_id];
      const size_t expected_face_points = num_elements * num_face_points;
      ASSERT(packaged_face_data.number_of_grid_points() == expected_face_points,
             "PackageLocalFacesBatched must run before device validation of "
             "DemandOutgoingCharSpeeds. Expected "
                 << expected_face_points << " packaged points on face "
                 << local_face_id << ", but found "
                 << packaged_face_data.number_of_grid_points() << ".");
      const auto packaged_face_data_view = packaged_face_data.view();
      const auto boundary_condition_kind_device =
          external_boundary_condition_kind_for_all_elements[local_face_id];
      constexpr int demand_outgoing_kind = static_cast<int>(
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
              DemandOutgoingCharSpeeds);
      size_t number_of_ingoing_speeds = 0;
      ::Kokkos::parallel_reduce(
          "GhDemandOutgoingCharSpeedsBatched", expected_face_points,
          KOKKOS_LAMBDA(const int linear_index_int,
                        size_t& local_number_of_ingoing_speeds) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t element_index = linear_index / num_face_points;
            if (boundary_condition_kind_device(element_index) !=
                demand_outgoing_kind) {
              return;
            }
            for (size_t speed_index = 0;
                 speed_index < num_char_speed_components; ++speed_index) {
              const double speed = packaged_face_data_view(
                  linear_index, offset_char_speeds + speed_index);
              if (not(speed >= 0.0)) {
                ++local_number_of_ingoing_speeds;
              }
            }
          },
          number_of_ingoing_speeds);
      if (number_of_ingoing_speeds > 0) {
        ERROR("DemandOutgoingCharSpeeds device validation failed on face "
              << local_face_id << " with " << number_of_ingoing_speeds
              << " negative or non-finite characteristic speeds.");
      }
    }
  }

  if (not has_constraint_preserving_physical_face) {
    return;
  }

  const auto evolved_device = packed_evolution_state->device_variables.view();
  const auto dt_device = packed_evolution_state->device_dt_variables.view();
  const auto derivatives_device =
      packed_boundary_scratch->volume_partial_derivatives.view();
  const auto gauge_device = packed_boundary_scratch->volume_gauge_data.view();
  const auto three_plus_one_device =
      packed_boundary_scratch->volume_three_plus_one_data.view();
  const auto gamma1_device = device_gamma1.get();
  const auto gamma2_device = device_gamma2.get();
  const auto inertial_coordinates = packed_geometry.inertial_coordinates_device;
  const auto element_point_offsets_device =
      packed_topology.element_point_offsets_device;
  const gh::BoundaryConditions::ConstraintPreservingBjorhus<volume_dim>
      cp_physical_boundary_condition{
          gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType::
              ConstraintPreservingPhysical};
  constexpr std::array<const char*, 2 * volume_dim> input_labels{
      {"GhBjorhusInput0", "GhBjorhusInput1", "GhBjorhusInput2",
       "GhBjorhusInput3", "GhBjorhusInput4", "GhBjorhusInput5"}};
  constexpr std::array<const char*, 2 * volume_dim> correction_labels{
      {"GhBjorhusCorrection0", "GhBjorhusCorrection1", "GhBjorhusCorrection2",
       "GhBjorhusCorrection3", "GhBjorhusCorrection4", "GhBjorhusCorrection5"}};

  for (size_t local_face_id = 0; local_face_id < 2 * volume_dim;
       ++local_face_id) {
    const size_t sliced_dim = local_face_id / 2;
    const size_t side_i = local_face_id % 2;
    const Side side = side_i == 0 ? Side::Lower : Side::Upper;
    const auto& external_face_mask = external_face_masks_host[local_face_id];
    const auto& external_boundary_condition_kind =
        external_boundary_condition_kinds_host[local_face_id];
    bool face_has_constraint_preserving_physical = false;
    for (const int kind : external_boundary_condition_kind) {
      face_has_constraint_preserving_physical =
          face_has_constraint_preserving_physical or
          kind ==
              static_cast<int>(
                  evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
                      ConstraintPreservingPhysical);
    }
    if (not face_has_constraint_preserving_physical) {
      continue;
    }

    const auto local_face_to_volume_index_device =
        side == Side::Upper
            ? gsl::at(face_to_volume_index_map, sliced_dim).second
            : gsl::at(face_to_volume_index_map, sliced_dim).first;
    const double outward = outward_sign(side);
    const size_t expected_face_points = num_elements * num_face_points;
    auto& input_device =
        packed_boundary_scratch->external_boundary_input_data[local_face_id];
    if (input_device.extent(0) != expected_face_points or
        input_device.extent(1) != num_bjorhus_input_components) {
      input_device = typename packed_boundary_scratch_type::
          external_boundary_workspace_type{input_labels[local_face_id],
                                           expected_face_points,
                                           num_bjorhus_input_components};
    }
    auto& correction_device =
        packed_boundary_scratch->external_boundary_corrections[local_face_id];
    if (correction_device.extent(0) != expected_face_points or
        correction_device.extent(1) != num_dt_components) {
      correction_device = typename packed_boundary_scratch_type::
          external_boundary_workspace_type{correction_labels[local_face_id],
                                           expected_face_points,
                                           num_dt_components};
    }
    const auto boundary_condition_kind_device =
        external_boundary_condition_kind_for_all_elements[local_face_id];
    constexpr int constraint_preserving_physical_kind = static_cast<int>(
        evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
            ConstraintPreservingPhysical);
    ::Kokkos::parallel_for(
        "GhGatherBjorhusBoundaryData", expected_face_points,
        KOKKOS_LAMBDA(const int linear_index_int) {
          const size_t linear_index = static_cast<size_t>(linear_index_int);
          const size_t face_point = linear_index % num_face_points;
          const size_t element_index = linear_index / num_face_points;
          if (boundary_condition_kind_device(element_index) !=
              constraint_preserving_physical_kind) {
            return;
          }
          const size_t local_volume_index =
              local_face_to_volume_index_device(face_point);
          const size_t volume_index =
              element_point_offsets_device(element_index) + local_volume_index;
          for (size_t component = 0; component < num_evolved_components;
               ++component) {
            input_device(linear_index, bjorhus_offset_evolved + component) =
                evolved_device(volume_index, component);
          }
          for (size_t component = 0; component < num_dt_components;
               ++component) {
            input_device(linear_index, bjorhus_offset_dt + component) =
                dt_device(volume_index, component);
          }
          for (size_t component = 0; component < num_derivative_components;
               ++component) {
            input_device(linear_index, bjorhus_offset_derivatives + component) =
                derivatives_device(volume_index, component);
          }
          for (size_t component = 0; component < num_gauge_components;
               ++component) {
            input_device(linear_index, bjorhus_offset_gauge + component) =
                gauge_device(volume_index, component);
          }
          for (size_t component = 0; component < num_three_plus_one_components;
               ++component) {
            input_device(linear_index,
                         bjorhus_offset_three_plus_one + component) =
                three_plus_one_device(volume_index, component);
          }
          input_device(linear_index, bjorhus_offset_gamma1) =
              gamma1_device(volume_index);
          input_device(linear_index, bjorhus_offset_gamma2) =
              gamma2_device(volume_index);
          for (size_t d = 0; d < volume_dim; ++d) {
            input_device(linear_index, bjorhus_offset_coordinates + d) =
                inertial_coordinates.get(d)(volume_index);
            input_device(linear_index, bjorhus_offset_mesh_velocity + d) =
                mesh_velocity.get(d)(volume_index);
            input_device(linear_index, bjorhus_offset_unnormalized_normal + d) =
                outward * inverse_jacobian(element_index, local_volume_index,
                                           sliced_dim * volume_dim + d);
          }
        });
    const auto input_host = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, input_device);
    auto correction_host = ::Kokkos::create_mirror_view(correction_device);
    ::Kokkos::deep_copy(correction_host, 0.0);

    for (size_t element_index = 0; element_index < num_elements;
         ++element_index) {
      if (external_face_mask[element_index] == 0) {
        continue;
      }
      const auto boundary_condition_kind = static_cast<
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind>(
          external_boundary_condition_kind[element_index]);
      if (boundary_condition_kind ==
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind::None) {
        ERROR("External face " << local_face_id << " of element "
                               << elements[element_index].id()
                               << " has no encoded boundary-condition kind.");
      }
      if (boundary_condition_kind ==
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
              DemandOutgoingCharSpeeds) {
        continue;
      }
      if (boundary_condition_kind !=
          evolution::Kokkos::Batched::ExternalBoundaryConditionKind::
              ConstraintPreservingPhysical) {
        ERROR("Unsupported encoded external boundary-condition kind "
              << static_cast<int>(boundary_condition_kind) << " on face "
              << local_face_id << " of element " << elements[element_index].id()
              << ".");
      }
      const size_t face_point_offset = element_index * num_face_points;
      evolved_variables_type face_evolved{num_face_points};
      dt_variables_type face_dt{num_face_points};
      host_derivatives_type face_derivatives{num_face_points};
      host_gauge_data_type face_gauge_data{num_face_points};
      host_three_plus_one_data_type face_three_plus_one{num_face_points};
      load_face_variables(make_not_null(&face_evolved), input_host,
                          face_point_offset, bjorhus_offset_evolved,
                          num_face_points);
      load_face_variables(make_not_null(&face_dt), input_host,
                          face_point_offset, bjorhus_offset_dt,
                          num_face_points);
      load_face_variables(make_not_null(&face_derivatives), input_host,
                          face_point_offset, bjorhus_offset_derivatives,
                          num_face_points);
      load_face_variables(make_not_null(&face_gauge_data), input_host,
                          face_point_offset, bjorhus_offset_gauge,
                          num_face_points);
      load_face_variables(make_not_null(&face_three_plus_one), input_host,
                          face_point_offset, bjorhus_offset_three_plus_one,
                          num_face_points);

      const auto& face_spacetime_metric =
          get<spacetime_metric_tag>(face_evolved);
      const auto& face_pi = get<pi_tag>(face_evolved);
      const auto& face_phi = get<phi_tag>(face_evolved);
      const auto& face_dt_spacetime_metric =
          get<::Tags::dt<spacetime_metric_tag>>(face_dt);
      const auto& face_dt_pi = get<::Tags::dt<pi_tag>>(face_dt);
      const auto& face_dt_phi = get<::Tags::dt<phi_tag>>(face_dt);
      const auto& face_d_spacetime_metric =
          get<d_spacetime_metric_tag>(face_derivatives);
      const auto& face_d_pi = get<d_pi_tag>(face_derivatives);
      const auto& face_d_phi = get<d_phi_tag>(face_derivatives);
      const auto& face_gauge_h = get<gauge_h_tag>(face_gauge_data);
      const auto& face_spacetime_deriv_gauge_h =
          get<spacetime_deriv_gauge_h_tag>(face_gauge_data);
      const auto& face_inverse_spatial_metric =
          get<inverse_spatial_metric_tag>(face_three_plus_one);
      const auto& face_shift = get<shift_tag>(face_three_plus_one);
      const auto& face_lapse = get<lapse_tag>(face_three_plus_one);
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
      Scalar<DataVector> face_gamma1{num_face_points};
      Scalar<DataVector> face_gamma2{num_face_points};
      DataVector normal_magnitude{num_face_points, 0.0};

      for (size_t p = 0; p < num_face_points; ++p) {
        const size_t input_point = face_point_offset + p;
        get(face_gamma1)[p] = input_host(input_point, bjorhus_offset_gamma1);
        get(face_gamma2)[p] = input_host(input_point, bjorhus_offset_gamma2);
        for (size_t d = 0; d < volume_dim; ++d) {
          face_coords.get(d)[p] =
              input_host(input_point, bjorhus_offset_coordinates + d);
          face_mesh_velocity.get(d)[p] =
              input_host(input_point, bjorhus_offset_mesh_velocity + d);
          normal_covector.get(d)[p] =
              input_host(input_point, bjorhus_offset_unnormalized_normal + d);
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

      dt_variables_type face_correction{num_face_points, 0.0};
      auto& dt_spacetime_metric_correction =
          get<::Tags::dt<spacetime_metric_tag>>(face_correction);
      auto& dt_pi_correction = get<::Tags::dt<pi_tag>>(face_correction);
      auto& dt_phi_correction = get<::Tags::dt<phi_tag>>(face_correction);
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
      if (error.has_value()) {
        ERROR("ConstraintPreservingPhysical failed on external face "
              << local_face_id << " of element " << elements[element_index].id()
              << ": " << error.value());
      }

      // TimeDerivative-type boundary corrections are added pointwise at the
      // face collocation points (production uses add_slice_to_data for
      // Gauss-Lobatto meshes); they are not flux corrections, so no mass
      // lifting or normal-magnitude weighting is applied.
      for (size_t component = 0; component < num_dt_components; ++component) {
        for (size_t p = 0; p < num_face_points; ++p) {
          correction_host(face_point_offset + p, component) =
              face_correction.data()[component * num_face_points + p];
        }
      }
    }
    ::Kokkos::deep_copy(correction_device, correction_host);
    ::Kokkos::parallel_for(
        "GhLiftBjorhusBoundaryCorrection", expected_face_points,
        KOKKOS_LAMBDA(const int linear_index_int) {
          const size_t linear_index = static_cast<size_t>(linear_index_int);
          const size_t face_point = linear_index % num_face_points;
          const size_t element_index = linear_index / num_face_points;
          if (boundary_condition_kind_device(element_index) !=
              constraint_preserving_physical_kind) {
            return;
          }
          const size_t volume_index =
              element_point_offsets_device(element_index) +
              local_face_to_volume_index_device(face_point);
          for (size_t component = 0; component < num_dt_components;
               ++component) {
            dt_device(volume_index, component) +=
                correction_device(linear_index, component);
          }
        });
  }
}

}  // namespace gh::Actions
