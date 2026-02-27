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
#include "Domain/Structure/Direction.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
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

constexpr size_t face_index(const size_t sliced_dim, const size_t side_i) {
  return 2 * sliced_dim + side_i;
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

KOKKOS_INLINE_FUNCTION void compute_packaged_boundary_data_at_point(
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_spacetime_metric,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_zero,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_plus,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_v_minus,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::iaa<double, volume_dim, Frame::Inertial>*>
        char_speed_n_times_v_minus,
    const gsl::not_null<tnsr::aa<double, volume_dim, Frame::Inertial>*>
        char_speed_gamma2_v_spacetime_metric,
    const gsl::not_null<tnsr::a<double, volume_dim, Frame::Inertial>*>
        char_speeds,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<double, volume_dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, volume_dim, Frame::Inertial>& phi,
    const double gamma1, const double gamma2,
    const tnsr::i<double, volume_dim, Frame::Inertial>&
        unnormalized_normal_covector) {
  tnsr::II<double, volume_dim, Frame::Inertial> inverse_spatial_metric{};
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det(make_not_null(&inverse_spatial_metric),
                                 make_not_null(&det_spatial_metric),
                                 spacetime_metric);
  (void)det_spatial_metric;

  tnsr::I<double, volume_dim, Frame::Inertial> shift{};
  for (size_t i = 0; i < volume_dim; ++i) {
    shift.get(i) = 0.0;
    for (size_t j = 0; j < volume_dim; ++j) {
      shift.get(i) +=
          inverse_spatial_metric.get(i, j) * spacetime_metric.get(0, j + 1);
    }
  }
  double lapse_squared = -spacetime_metric.get(0, 0);
  for (size_t i = 0; i < volume_dim; ++i) {
    lapse_squared += shift.get(i) * spacetime_metric.get(0, i + 1);
  }
  const double lapse = sqrt(lapse_squared);

  tnsr::I<double, volume_dim, Frame::Inertial> normal_vector{};
  double normal_magnitude_squared = 0.0;
  for (size_t i = 0; i < volume_dim; ++i) {
    normal_vector.get(i) = 0.0;
    for (size_t j = 0; j < volume_dim; ++j) {
      normal_vector.get(i) += inverse_spatial_metric.get(i, j) *
                              unnormalized_normal_covector.get(j);
    }
    normal_magnitude_squared +=
        normal_vector.get(i) * unnormalized_normal_covector.get(i);
  }
  const double one_over_normal_magnitude = 1.0 / sqrt(normal_magnitude_squared);
  tnsr::i<double, volume_dim, Frame::Inertial> normal_covector{};
  for (size_t i = 0; i < volume_dim; ++i) {
    normal_vector.get(i) *= one_over_normal_magnitude;
    normal_covector.get(i) =
        unnormalized_normal_covector.get(i) * one_over_normal_magnitude;
  }

  double shift_dot_normal = 0.0;
  for (size_t i = 0; i < volume_dim; ++i) {
    shift_dot_normal += shift.get(i) * normal_covector.get(i);
  }
  shift_dot_normal *= -1.0;

  char_speeds->get(0) = (1.0 + gamma1) * shift_dot_normal;
  char_speeds->get(1) = shift_dot_normal;
  char_speeds->get(2) = lapse + shift_dot_normal;
  char_speeds->get(3) = -lapse + shift_dot_normal;

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      char_speed_gamma2_v_spacetime_metric->get(a, b) =
          gamma2 * spacetime_metric.get(a, b);
    }
  }

  tnsr::aa<double, volume_dim, Frame::Inertial> normal_dot_phi{};
  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      normal_dot_phi.get(a, b) = normal_vector.get(0) * phi.get(0, a, b);
      for (size_t i = 1; i < volume_dim; ++i) {
        normal_dot_phi.get(a, b) += normal_vector.get(i) * phi.get(i, a, b);
      }
    }
  }

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      char_speed_v_plus->get(a, b) =
          char_speeds->get(2) *
          (pi.get(a, b) + normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));
      char_speed_v_minus->get(a, b) =
          char_speeds->get(3) *
          (pi.get(a, b) - normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));

      for (size_t i = 0; i < volume_dim; ++i) {
        char_speed_v_zero->get(i, a, b) =
            char_speeds->get(1) *
            (phi.get(i, a, b) -
             normal_covector.get(i) * normal_dot_phi.get(a, b));
      }
    }
  }

  for (size_t a = 0; a < volume_dim + 1; ++a) {
    for (size_t b = a; b < volume_dim + 1; ++b) {
      for (size_t i = 0; i < volume_dim; ++i) {
        char_speed_n_times_v_plus->get(i, a, b) =
            char_speed_v_plus->get(a, b) * normal_covector.get(i);
        char_speed_n_times_v_minus->get(i, a, b) =
            char_speed_v_minus->get(a, b) * normal_covector.get(i);
      }
      char_speed_v_spacetime_metric->get(a, b) =
          char_speeds->get(0) * spacetime_metric.get(a, b);
      char_speed_gamma2_v_spacetime_metric->get(a, b) *= char_speeds->get(0);
    }
  }
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

KOKKOS_INLINE_FUNCTION double step_function_double(const double value) {
  return value < 0.0 ? 0.0 : 1.0;
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

  const auto& elements = packed_topology.local_elements;
  const auto& face_to_volume_index_map =
      packed_topology.device_face_to_volume_index_map;
  const auto& inverse_jacobian =
      packed_geometry.element_inverse_jacobian_device;
  const auto& element_point_offsets_device =
      packed_topology.element_point_offsets_device;
  const auto& inertial_coordinates =
      packed_geometry.inertial_coordinates_device;
  const auto& gamma1 = device_gamma1;
  const auto& gamma2 = device_gamma2;
  const auto& packaged_face_data_for_all_elements =
      packed_boundary_scratch.packaged_face_data_for_all_elements;
  const auto& external_face_mask_for_all_elements =
      packed_boundary_metadata.external_face_mask_for_all_elements;

  auto& dt_vars = packed_evolution_state->device_dt_variables;
  const auto dt_spacetime_metric =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(dt_vars);
  const auto dt_pi = get<::Tags::MirrorView<
      ::Tags::dt<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>>>(
      dt_vars);
  const auto dt_phi = get<::Tags::MirrorView<
      ::Tags::dt<gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>>>(
      dt_vars);

  const size_t num_elements = elements.size();
  for (size_t sliced_dim = 0; sliced_dim < volume_dim; ++sliced_dim) {
    for (size_t side_i = 0; side_i < 2; ++side_i) {
      const Side side = side_i == 0 ? Side::Lower : Side::Upper;
      const size_t local_face_id = face_index(sliced_dim, side_i);
      const auto packaged_face_data_view =
          packaged_face_data_for_all_elements[local_face_id].view();
      ASSERT(packaged_face_data_view.extent(1) ==
                 packaged_boundary_data_components,
             "Unexpected number of packaged boundary components for GH batched "
             "external boundary corrections.");
      ASSERT(packaged_face_data_for_all_elements[local_face_id]
                     .number_of_grid_points() == num_elements * num_face_points,
             "PackageLocalFacesBatched must run before external boundary "
             "corrections.");
      const auto external_face_mask =
          external_face_mask_for_all_elements[local_face_id];
      ASSERT(external_face_mask.extent(0) == num_elements,
             "External-face mask extent mismatch. "
             "InitializeBoundaryBatchMetadata must run before external "
             "boundary corrections.");

      const auto local_face_to_volume_index =
          side == Side::Upper
              ? gsl::at(face_to_volume_index_map, sliced_dim).second
              : gsl::at(face_to_volume_index_map, sliced_dim).first;
      const double outward = outward_sign(side);
      const double lift_prefactor =
          -0.5 *
          static_cast<double>(extents[sliced_dim] * (extents[sliced_dim] - 1));

      ::Kokkos::parallel_for(
          "ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched",
          num_elements * num_face_points,
          KOKKOS_LAMBDA(const int linear_index_int) {
            const size_t linear_index = static_cast<size_t>(linear_index_int);
            const size_t face_index_on_face = linear_index % num_face_points;
            const size_t element_index = linear_index / num_face_points;
            if (external_face_mask(element_index) == 0) {
              return;
            }

            const size_t local_volume_index_in_element =
                local_face_to_volume_index(face_index_on_face);
            const size_t local_volume_index =
                element_point_offsets_device(element_index) +
                local_volume_index_in_element;

            tnsr::i<double, volume_dim, Frame::Inertial>
                unnormalized_normal_covector{};
            double normal_magnitude_squared = 0.0;
            for (size_t d = 0; d < volume_dim; ++d) {
              const double component =
                  outward * inverse_jacobian(element_index,
                                             local_volume_index_in_element,
                                             sliced_dim * volume_dim + d);
              unnormalized_normal_covector.get(d) = component;
              normal_magnitude_squared += component * component;
            }
            const double normal_magnitude = sqrt(normal_magnitude_squared);

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
            load_aa_from_packaged_data(make_not_null(&local_v_spacetime_metric),
                                       packaged_face_data_view, linear_index,
                                       offset_v_spacetime_metric);
            load_iaa_from_packaged_data(make_not_null(&local_v_zero),
                                        packaged_face_data_view, linear_index,
                                        offset_v_zero);
            load_aa_from_packaged_data(make_not_null(&local_v_plus),
                                       packaged_face_data_view, linear_index,
                                       offset_v_plus);
            load_aa_from_packaged_data(make_not_null(&local_v_minus),
                                       packaged_face_data_view, linear_index,
                                       offset_v_minus);
            load_iaa_from_packaged_data(
                make_not_null(&local_normal_times_v_plus),
                packaged_face_data_view, linear_index,
                offset_normal_times_v_plus);
            load_iaa_from_packaged_data(
                make_not_null(&local_normal_times_v_minus),
                packaged_face_data_view, linear_index,
                offset_normal_times_v_minus);
            load_aa_from_packaged_data(
                make_not_null(&local_gamma2_v_spacetime_metric),
                packaged_face_data_view, linear_index,
                offset_gamma2_v_spacetime_metric);
            load_a_from_packaged_data(make_not_null(&local_char_speeds),
                                      packaged_face_data_view, linear_index,
                                      offset_char_speeds);

            tnsr::I<double, volume_dim, Frame::Inertial> inertial_coords{};
            for (size_t d = 0; d < volume_dim; ++d) {
              inertial_coords.get(d) =
                  inertial_coordinates.get(d)(local_volume_index);
            }

            tnsr::aa<double, volume_dim, Frame::Inertial>
                exterior_spacetime_metric{};
            tnsr::aa<double, volume_dim, Frame::Inertial> exterior_pi{};
            tnsr::iaa<double, volume_dim, Frame::Inertial> exterior_phi{};
            compute_hardcoded_schwarzschild_gh_fields(
                make_not_null(&exterior_spacetime_metric),
                make_not_null(&exterior_pi), make_not_null(&exterior_phi),
                inertial_coords);

            tnsr::i<double, volume_dim, Frame::Inertial>
                exterior_unnormalized_normal_covector{};
            for (size_t d = 0; d < volume_dim; ++d) {
              exterior_unnormalized_normal_covector.get(d) =
                  -unnormalized_normal_covector.get(d);
            }

            const auto gamma1_at_s = make_at_index(gamma1, local_volume_index);
            const auto gamma2_at_s = make_at_index(gamma2, local_volume_index);
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
            compute_packaged_boundary_data_at_point(
                make_not_null(&remote_v_spacetime_metric),
                make_not_null(&remote_v_zero), make_not_null(&remote_v_plus),
                make_not_null(&remote_v_minus),
                make_not_null(&remote_normal_times_v_plus),
                make_not_null(&remote_normal_times_v_minus),
                make_not_null(&remote_gamma2_v_spacetime_metric),
                make_not_null(&remote_char_speeds), exterior_spacetime_metric,
                exterior_pi, exterior_phi, get(gamma1_at_s), get(gamma2_at_s),
                exterior_unnormalized_normal_covector);

            tnsr::aa<double, volume_dim, Frame::Inertial>
                dt_spacetime_metric_correction{};
            tnsr::aa<double, volume_dim, Frame::Inertial> dt_pi_correction{};
            tnsr::iaa<double, volume_dim, Frame::Inertial> dt_phi_correction{};
            compute_boundary_terms_at_point(
                make_not_null(&dt_spacetime_metric_correction),
                make_not_null(&dt_pi_correction),
                make_not_null(&dt_phi_correction), local_v_spacetime_metric,
                local_v_zero, local_v_plus, local_v_minus,
                local_normal_times_v_plus, local_normal_times_v_minus,
                local_gamma2_v_spacetime_metric, local_char_speeds,
                remote_v_spacetime_metric, remote_v_zero, remote_v_plus,
                remote_v_minus, remote_normal_times_v_plus,
                remote_normal_times_v_minus, remote_gamma2_v_spacetime_metric,
                remote_char_speeds);

            const double lifted_factor = lift_prefactor * normal_magnitude;
            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                dt_spacetime_metric.get(a, b)[local_volume_index] +=
                    lifted_factor * dt_spacetime_metric_correction.get(a, b);
                dt_pi.get(a, b)[local_volume_index] +=
                    lifted_factor * dt_pi_correction.get(a, b);
                for (size_t d = 0; d < volume_dim; ++d) {
                  dt_phi.get(d, a, b)[local_volume_index] +=
                      lifted_factor * dt_phi_correction.get(d, a, b);
                }
              }
            }
          });
    }
  }
}

}  // namespace gh::Actions
