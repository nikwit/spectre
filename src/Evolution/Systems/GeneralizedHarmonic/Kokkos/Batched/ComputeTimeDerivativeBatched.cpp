// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/ComputeTimeDerivativeBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
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
using device_inertial_coordinates_type =
    gh::KokkosTags::DeviceInertialCoordinates<dim>::type;

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
using device_spatial_deriv_gauge_h_tag =
    ::Tags::deriv<device_gauge_h_tag, tmpl::size_t<dim>, Frame::Inertial>;
using device_spatial_deriv_gauge_data_type =
    Variables<tmpl::list<device_spatial_deriv_gauge_h_tag>>;
using device_d_spacetime_metric_tag =
    ::Tags::deriv<device_spacetime_metric_tag, tmpl::size_t<dim>,
                  Frame::Inertial>;
using device_d_pi_tag =
    ::Tags::deriv<device_pi_tag, tmpl::size_t<dim>, Frame::Inertial>;
using device_d_phi_tag =
    ::Tags::deriv<device_phi_tag, tmpl::size_t<dim>, Frame::Inertial>;

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

KOKKOS_INLINE_FUNCTION double spacetime_metric_deriv_component(
    const size_t a, const size_t mu, const size_t nu, const double lapse,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const tnsr::aa<double, dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, dim, Frame::Inertial>& phi) {
  if (a == 0) {
    double dt_spacetime_metric_mu_nu = -lapse * pi.get(mu, nu);
    for (size_t m = 0; m < dim; ++m) {
      dt_spacetime_metric_mu_nu += shift.get(m) * phi.get(m, mu, nu);
    }
    return dt_spacetime_metric_mu_nu;
  }
  return phi.get(a - 1, mu, nu);
}

KOKKOS_INLINE_FUNCTION double christoffel_first_kind_component(
    const size_t k, const size_t i, const size_t j, const double lapse,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const tnsr::aa<double, dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, dim, Frame::Inertial>& phi) {
  return 0.5 *
         (spacetime_metric_deriv_component(i, j, k, lapse, shift, pi, phi) +
          spacetime_metric_deriv_component(j, i, k, lapse, shift, pi, phi) -
          spacetime_metric_deriv_component(k, i, j, lapse, shift, pi, phi));
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

KOKKOS_INLINE_FUNCTION void compute_gauge_from_gh_vars(
    const gsl::not_null<tnsr::a<double, dim, Frame::Inertial>*> gauge_h,
    const tnsr::aa<double, dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<double, dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, dim, Frame::Inertial>& phi) {
  tnsr::II<double, dim, Frame::Inertial> inverse_spatial_metric{};
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det(make_not_null(&inverse_spatial_metric),
                                 make_not_null(&det_spatial_metric),
                                 spacetime_metric);
  (void)det_spatial_metric;

  tnsr::I<double, dim, Frame::Inertial> shift{};
  for (size_t i = 0; i < dim; ++i) {
    shift.get(i) = 0.0;
    for (size_t j = 0; j < dim; ++j) {
      shift.get(i) +=
          inverse_spatial_metric.get(i, j) * spacetime_metric.get(0, j + 1);
    }
  }
  double lapse_squared = -spacetime_metric.get(0, 0);
  for (size_t i = 0; i < dim; ++i) {
    lapse_squared += shift.get(i) * spacetime_metric.get(0, i + 1);
  }
  const double lapse = sqrt(lapse_squared);
  const double inv_lapse_squared = 1.0 / lapse_squared;

  tnsr::AA<double, dim, Frame::Inertial> inverse_spacetime_metric{};
  inverse_spacetime_metric.get(0, 0) = -inv_lapse_squared;
  for (size_t i = 0; i < dim; ++i) {
    inverse_spacetime_metric.get(0, i + 1) = shift.get(i) * inv_lapse_squared;
    for (size_t j = i; j < dim; ++j) {
      inverse_spacetime_metric.get(i + 1, j + 1) =
          inverse_spatial_metric.get(i, j) -
          shift.get(i) * shift.get(j) * inv_lapse_squared;
    }
  }

  tnsr::a<double, dim, Frame::Inertial> trace_christoffel{};
  tnsr::a<double, dim, Frame::Inertial> spacetime_normal_one_form{};
  tnsr::A<double, dim, Frame::Inertial> spacetime_normal_vector{};
  spacetime_normal_one_form.get(0) = -lapse;
  for (size_t a = 1; a < dim + 1; ++a) {
    spacetime_normal_one_form.get(a) = 0.0;
  }
  spacetime_normal_vector.get(0) = 1.0 / lapse;
  for (size_t i = 0; i < dim; ++i) {
    spacetime_normal_vector.get(i + 1) = -shift.get(i) / lapse;
  }

  trace_christoffel.get(0) = 0.0;
  for (size_t b = 0; b < dim + 1; ++b) {
    trace_christoffel.get(0) -=
        0.5 * pi.get(b, b) * inverse_spacetime_metric.get(b, b);
    for (size_t i = 0; i < dim; ++i) {
      trace_christoffel.get(0) -= 0.5 * spacetime_normal_vector.get(i + 1) *
                                  phi.get(i, b, b) *
                                  inverse_spacetime_metric.get(b, b);
    }
    for (size_t c = b + 1; c < dim + 1; ++c) {
      trace_christoffel.get(0) -=
          pi.get(b, c) * inverse_spacetime_metric.get(b, c);
      for (size_t i = 0; i < dim; ++i) {
        trace_christoffel.get(0) -= spacetime_normal_vector.get(i + 1) *
                                    phi.get(i, b, c) *
                                    inverse_spacetime_metric.get(b, c);
      }
    }
  }

  for (size_t a = 1; a < dim + 1; ++a) {
    trace_christoffel.get(a) =
        spacetime_normal_one_form.get(a) * trace_christoffel.get(0);
    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        trace_christoffel.get(a) +=
            inverse_spatial_metric.get(i, j) * phi.get(i, j + 1, a);
      }
    }
    for (size_t b = 0; b < dim + 1; ++b) {
      trace_christoffel.get(a) += spacetime_normal_vector.get(b) * pi.get(b, a);
      trace_christoffel.get(a) -=
          0.5 * phi.get(a - 1, b, b) * inverse_spacetime_metric.get(b, b);
      for (size_t c = b + 1; c < dim + 1; ++c) {
        trace_christoffel.get(a) -=
            phi.get(a - 1, b, c) * inverse_spacetime_metric.get(b, c);
      }
    }
  }

  trace_christoffel.get(0) *= spacetime_normal_one_form.get(0);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      trace_christoffel.get(0) +=
          inverse_spatial_metric.get(i, j) * phi.get(i, j + 1, 0);
    }
  }
  for (size_t b = 0; b < dim + 1; ++b) {
    trace_christoffel.get(0) += spacetime_normal_vector.get(b) * pi.get(b, 0);
  }

  for (size_t a = 0; a < dim + 1; ++a) {
    gauge_h->get(a) = -trace_christoffel.get(a);
  }
}

KOKKOS_INLINE_FUNCTION void compute_spacetime_metric_from_3plus1(
    const gsl::not_null<tnsr::aa<double, dim, Frame::Inertial>*>
        spacetime_metric,
    const Scalar<double>& lapse,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const tnsr::ii<double, dim, Frame::Inertial>& spatial_metric) {
  spacetime_metric->get(0, 0) = -get(lapse) * get(lapse);
  for (size_t m = 0; m < dim; ++m) {
    spacetime_metric->get(0, 0) +=
        spatial_metric.get(m, m) * shift.get(m) * shift.get(m);
    for (size_t n = 0; n < m; ++n) {
      spacetime_metric->get(0, 0) +=
          2.0 * spatial_metric.get(m, n) * shift.get(m) * shift.get(n);
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    spacetime_metric->get(0, i + 1) = 0.0;
    for (size_t m = 0; m < dim; ++m) {
      spacetime_metric->get(0, i + 1) +=
          spatial_metric.get(m, i) * shift.get(m);
    }
    for (size_t j = i; j < dim; ++j) {
      spacetime_metric->get(i + 1, j + 1) = spatial_metric.get(i, j);
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_phi_from_3plus1(
    const gsl::not_null<tnsr::iaa<double, dim, Frame::Inertial>*> phi,
    const Scalar<double>& lapse,
    const tnsr::i<double, dim, Frame::Inertial>& deriv_lapse,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const tnsr::iJ<double, dim, Frame::Inertial>& deriv_shift,
    const tnsr::ii<double, dim, Frame::Inertial>& spatial_metric,
    const tnsr::ijj<double, dim, Frame::Inertial>& deriv_spatial_metric) {
  for (size_t k = 0; k < dim; ++k) {
    phi->get(k, 0, 0) = -2.0 * get(lapse) * deriv_lapse.get(k);
    for (size_t m = 0; m < dim; ++m) {
      for (size_t n = 0; n < dim; ++n) {
        phi->get(k, 0, 0) +=
            deriv_spatial_metric.get(k, m, n) * shift.get(m) * shift.get(n) +
            2.0 * spatial_metric.get(m, n) * shift.get(m) *
                deriv_shift.get(k, n);
      }
    }
    for (size_t i = 0; i < dim; ++i) {
      phi->get(k, 0, i + 1) = 0.0;
      for (size_t m = 0; m < dim; ++m) {
        phi->get(k, 0, i + 1) +=
            deriv_spatial_metric.get(k, m, i) * shift.get(m) +
            spatial_metric.get(m, i) * deriv_shift.get(k, m);
      }
      for (size_t j = i; j < dim; ++j) {
        phi->get(k, i + 1, j + 1) = deriv_spatial_metric.get(k, i, j);
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_pi_from_3plus1(
    const gsl::not_null<tnsr::aa<double, dim, Frame::Inertial>*> pi,
    const Scalar<double>& lapse, const Scalar<double>& dt_lapse,
    const tnsr::I<double, dim, Frame::Inertial>& shift,
    const tnsr::I<double, dim, Frame::Inertial>& dt_shift,
    const tnsr::ii<double, dim, Frame::Inertial>& spatial_metric,
    const tnsr::ii<double, dim, Frame::Inertial>& dt_spatial_metric,
    const tnsr::iaa<double, dim, Frame::Inertial>& phi) {
  pi->get(0, 0) = -2.0 * get(lapse) * get(dt_lapse);
  for (size_t m = 0; m < dim; ++m) {
    for (size_t n = 0; n < dim; ++n) {
      pi->get(0, 0) +=
          dt_spatial_metric.get(m, n) * shift.get(m) * shift.get(n) +
          2.0 * spatial_metric.get(m, n) * shift.get(m) * dt_shift.get(n);
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    pi->get(0, i + 1) = 0.0;
    for (size_t m = 0; m < dim; ++m) {
      pi->get(0, i + 1) += dt_spatial_metric.get(m, i) * shift.get(m) +
                           spatial_metric.get(m, i) * dt_shift.get(m);
    }
    for (size_t j = i; j < dim; ++j) {
      pi->get(i + 1, j + 1) = dt_spatial_metric.get(i, j);
    }
  }
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      for (size_t i = 0; i < dim; ++i) {
        pi->get(mu, nu) -= shift.get(i) * phi.get(i, mu, nu);
      }
      pi->get(mu, nu) /= -get(lapse);
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_hardcoded_schwarzschild_3plus1(
    const gsl::not_null<Scalar<double>*> lapse,
    const gsl::not_null<Scalar<double>*> dt_lapse,
    const gsl::not_null<tnsr::i<double, dim, Frame::Inertial>*> deriv_lapse,
    const gsl::not_null<tnsr::I<double, dim, Frame::Inertial>*> shift,
    const gsl::not_null<tnsr::I<double, dim, Frame::Inertial>*> dt_shift,
    const gsl::not_null<tnsr::iJ<double, dim, Frame::Inertial>*> deriv_shift,
    const gsl::not_null<tnsr::ii<double, dim, Frame::Inertial>*> spatial_metric,
    const gsl::not_null<tnsr::ii<double, dim, Frame::Inertial>*>
        dt_spatial_metric,
    const gsl::not_null<tnsr::ijj<double, dim, Frame::Inertial>*>
        deriv_spatial_metric,
    const tnsr::I<double, dim, Frame::Inertial>& inertial_coordinates) {
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

  const std::array<double, dim> coords{{x, y, z}};
  for (size_t i = 0; i < dim; ++i) {
    const double xi = coords[i];
    deriv_lapse->get(i) = mass * xi * inv_r_cubed * get(*lapse) * lapse_squared;
    shift->get(i) = 2.0 * mass * xi * inv_r_squared * lapse_squared;
    dt_shift->get(i) = 0.0;
  }

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = i; j < dim; ++j) {
      const double delta_ij = i == j ? 1.0 : 0.0;
      spatial_metric->get(i, j) =
          delta_ij + 2.0 * mass * coords[i] * coords[j] * inv_r_cubed;
      dt_spatial_metric->get(i, j) = 0.0;
    }
  }

  for (size_t k = 0; k < dim; ++k) {
    for (size_t i = 0; i < dim; ++i) {
      const double delta_ki = k == i ? 1.0 : 0.0;
      deriv_shift->get(k, i) =
          2.0 * mass * inv_r_squared * lapse_squared * delta_ki -
          4.0 * mass * coords[k] * coords[i] * inv_r_fourth * lapse_squared *
              (1.0 - mass * inv_r * lapse_squared);
      for (size_t j = i; j < dim; ++j) {
        const double delta_kj = k == j ? 1.0 : 0.0;
        deriv_spatial_metric->get(k, i, j) =
            -6.0 * mass * coords[i] * coords[j] * coords[k] * inv_r_fifth +
            2.0 * mass * coords[i] * delta_kj * inv_r_cubed +
            2.0 * mass * coords[j] * delta_ki * inv_r_cubed;
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void compute_hardcoded_analytic_gauge(
    const gsl::not_null<tnsr::a<double, dim, Frame::Inertial>*> gauge_h,
    const tnsr::I<double, dim, Frame::Inertial>& inertial_coordinates) {
  Scalar<double> lapse{};
  Scalar<double> dt_lapse{};
  tnsr::i<double, dim, Frame::Inertial> deriv_lapse{};
  tnsr::I<double, dim, Frame::Inertial> shift{};
  tnsr::I<double, dim, Frame::Inertial> dt_shift{};
  tnsr::iJ<double, dim, Frame::Inertial> deriv_shift{};
  tnsr::ii<double, dim, Frame::Inertial> spatial_metric{};
  tnsr::ii<double, dim, Frame::Inertial> dt_spatial_metric{};
  tnsr::ijj<double, dim, Frame::Inertial> deriv_spatial_metric{};

  compute_hardcoded_schwarzschild_3plus1(
      make_not_null(&lapse), make_not_null(&dt_lapse),
      make_not_null(&deriv_lapse), make_not_null(&shift),
      make_not_null(&dt_shift), make_not_null(&deriv_shift),
      make_not_null(&spatial_metric), make_not_null(&dt_spatial_metric),
      make_not_null(&deriv_spatial_metric), inertial_coordinates);

  tnsr::aa<double, dim, Frame::Inertial> analytic_spacetime_metric{};
  tnsr::iaa<double, dim, Frame::Inertial> analytic_phi{};
  tnsr::aa<double, dim, Frame::Inertial> analytic_pi{};

  compute_phi_from_3plus1(make_not_null(&analytic_phi), lapse, deriv_lapse,
                          shift, deriv_shift, spatial_metric,
                          deriv_spatial_metric);
  compute_pi_from_3plus1(make_not_null(&analytic_pi), lapse, dt_lapse, shift,
                         dt_shift, spatial_metric, dt_spatial_metric,
                         analytic_phi);
  compute_spacetime_metric_from_3plus1(
      make_not_null(&analytic_spacetime_metric), lapse, shift, spatial_metric);
  compute_gauge_from_gh_vars(gauge_h, analytic_spacetime_metric, analytic_pi,
                             analytic_phi);
}

void compute_hardcoded_analytic_gauge_and_spacetime_derivative(
    const gsl::not_null<device_gauge_data_type*> device_gauge_data,
    const gsl::not_null<device_spatial_deriv_gauge_data_type*>
        device_spatial_gauge_deriv,
    const device_inertial_coordinates_type& device_inertial_coordinates,
    const ::Kokkos::View<double***>& element_inverse_jacobian,
    const Mesh<dim>& mesh) {
  const size_t points_per_element = mesh.number_of_grid_points();
  const size_t number_of_elements = element_inverse_jacobian.extent(0);
  const size_t number_of_points = number_of_elements * points_per_element;
  ASSERT(
      device_inertial_coordinates.get(0).extent(0) == number_of_points and
          device_inertial_coordinates.get(1).extent(0) == number_of_points and
          device_inertial_coordinates.get(2).extent(0) == number_of_points,
      "Packed inertial coordinates must have one entry per packed point. "
          << "Expected " << number_of_points << " points from "
          << number_of_elements << " elements x " << points_per_element
          << " points/element, got ["
          << device_inertial_coordinates.get(0).extent(0) << ", "
          << device_inertial_coordinates.get(1).extent(0) << ", "
          << device_inertial_coordinates.get(2).extent(0) << "].");
  if (device_gauge_data->number_of_grid_points() != number_of_points) {
    device_gauge_data->initialize(number_of_points);
  }
  if (device_spatial_gauge_deriv->number_of_grid_points() != number_of_points) {
    device_spatial_gauge_deriv->initialize(number_of_points);
  }

  const auto device_gauge_h = get<device_gauge_h_tag>(*device_gauge_data);
  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeGaugeH", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        tnsr::I<double, dim, Frame::Inertial> inertial_coords_at_s{};
        for (int d = 0; d < static_cast<int>(dim); ++d) {
          inertial_coords_at_s.get(d) =
              device_inertial_coordinates.get(d)(point);
        }
        tnsr::a<double, dim, Frame::Inertial> gauge_h_at_s{};
        compute_hardcoded_analytic_gauge(make_not_null(&gauge_h_at_s),
                                         inertial_coords_at_s);
        for (int a = 0; a < static_cast<int>(dim + 1); ++a) {
          device_gauge_h.get(a)[point] = gauge_h_at_s.get(a);
        }
      });

  partial_derivatives_batched(device_spatial_gauge_deriv, *device_gauge_data,
                              mesh, element_inverse_jacobian);

  const auto device_spatial_deriv_gauge_h =
      get<device_spatial_deriv_gauge_h_tag>(*device_spatial_gauge_deriv);
  const auto device_spacetime_deriv_gauge_h =
      get<device_spacetime_deriv_gauge_h_tag>(*device_gauge_data);
  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDerivGaugeH", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        for (int a = 0; a < static_cast<int>(dim + 1); ++a) {
          device_spacetime_deriv_gauge_h.get(0, a)[point] = 0.0;
          for (int i = 0; i < static_cast<int>(dim); ++i) {
            device_spacetime_deriv_gauge_h.get(i + 1, a)[point] =
                device_spatial_deriv_gauge_h.get(i, a)[point];
          }
        }
      });
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

  auto& device_gauge_data = packed_boundary_scratch->volume_gauge_data;
  if (device_gauge_data.number_of_grid_points() != total_points) {
    device_gauge_data.initialize(total_points);
  }
  auto& device_spatial_gauge_deriv =
      packed_boundary_scratch->volume_spatial_deriv_gauge;
  if (device_spatial_gauge_deriv.number_of_grid_points() != total_points) {
    device_spatial_gauge_deriv.initialize(total_points);
  }
  compute_hardcoded_analytic_gauge_and_spacetime_derivative(
      make_not_null(&device_gauge_data),
      make_not_null(&device_spatial_gauge_deriv),
      packed_geometry.inertial_coordinates_device,
      packed_geometry.element_inverse_jacobian_device, mesh);

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
  const auto gauge_h = get<device_gauge_h_tag>(device_gauge_data);
  const auto spacetime_deriv_gauge_h =
      get<device_spacetime_deriv_gauge_h_tag>(device_gauge_data);

  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtSpacetimeMetric", total_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        const auto gamma1_at_s = make_at_index(device_constraint_gamma1, s);
        const double gamma1 = get(gamma1_at_s);

        const auto spacetime_metric_at_s = make_at_index(spacetime_metric, s);
        const auto pi_at_s = make_at_index(pi, s);
        const auto phi_at_s = make_at_index(phi, s);

        tnsr::II<double, dim, Frame::Inertial> inverse_spatial_metric{};
        tnsr::I<double, dim, Frame::Inertial> shift{};
        double lapse = 0.0;
        double inv_lapse_squared = 0.0;
        compute_inverse_spatial_metric_shift_lapse(
            make_not_null(&inverse_spatial_metric), make_not_null(&shift),
            make_not_null(&lapse), make_not_null(&inv_lapse_squared),
            spacetime_metric_at_s);
        (void)inv_lapse_squared;
        const double gamma1_plus_1 = 1.0 + gamma1;

        for (int mu = 0; mu < static_cast<int>(dim + 1); ++mu) {
          for (int nu = mu; nu < static_cast<int>(dim + 1); ++nu) {
            const double dt_spacetime_metric_piece =
                spacetime_metric_deriv_component(0, mu, nu, lapse, shift,
                                                 pi_at_s, phi_at_s);
            double shift_dot_three_index_constraint_mu_nu = 0.0;
            for (int m = 0; m < static_cast<int>(dim); ++m) {
              shift_dot_three_index_constraint_mu_nu +=
                  shift.get(m) * (d_spacetime_metric.get(m, mu, nu)[point] -
                                  phi_at_s.get(m, mu, nu));
            }
            dt_spacetime_metric.get(mu, nu)[point] =
                dt_spacetime_metric_piece +
                gamma1_plus_1 * shift_dot_three_index_constraint_mu_nu;
          }
        }
      });

  constexpr size_t dt_pi_offset_lapse = 0;
  constexpr size_t dt_pi_offset_half_pi_two_normals = dt_pi_offset_lapse + 1;
  constexpr size_t dt_pi_offset_gamma1gamma2 =
      dt_pi_offset_half_pi_two_normals + 1;
  constexpr size_t dt_pi_offset_minus_gamma0_lapse =
      dt_pi_offset_gamma1gamma2 + 1;
  constexpr size_t dt_pi_offset_normal_dot_gauge_constraint =
      dt_pi_offset_minus_gamma0_lapse + 1;
  constexpr size_t dt_pi_offset_shift =
      dt_pi_offset_normal_dot_gauge_constraint + 1;
  constexpr size_t dt_pi_offset_inverse_spatial_metric =
      dt_pi_offset_shift + dim;
  constexpr size_t dt_pi_offset_inverse_spacetime_metric =
      dt_pi_offset_inverse_spatial_metric + dim * dim;
  constexpr size_t dt_pi_offset_pi_up =
      dt_pi_offset_inverse_spacetime_metric + spacetime_dim * spacetime_dim;
  constexpr size_t dt_pi_offset_pi_one_normal =
      dt_pi_offset_pi_up + spacetime_dim * spacetime_dim;
  constexpr size_t dt_pi_offset_gauge_constraint =
      dt_pi_offset_pi_one_normal + spacetime_dim;
  constexpr size_t dt_pi_offset_gauge_h =
      dt_pi_offset_gauge_constraint + spacetime_dim;
  constexpr size_t dt_pi_offset_spacetime_metric_components =
      dt_pi_offset_gauge_h + spacetime_dim;
  constexpr size_t dt_pi_offset_pi_components =
      dt_pi_offset_spacetime_metric_components + spacetime_symmetric_size;
  constexpr size_t dt_pi_offset_phi_components =
      dt_pi_offset_pi_components + spacetime_symmetric_size;
  constexpr size_t dt_pi_scratch_size_doubles =
      dt_pi_offset_phi_components + dim * spacetime_symmetric_size;
  constexpr size_t dt_pi_scratch_size_bytes =
      dt_pi_scratch_size_doubles * sizeof(double);

  using dt_pi_team_policy = ::Kokkos::TeamPolicy<>;
  dt_pi_team_policy dt_pi_policy(total_points, 32);
  dt_pi_policy = dt_pi_policy.set_scratch_size(
      0, ::Kokkos::PerTeam(dt_pi_scratch_size_bytes));
  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtPi", dt_pi_policy,
      KOKKOS_LAMBDA(const dt_pi_team_policy::member_type& team_member) {
        const size_t point = team_member.league_rank();
        constexpr int int_dim = static_cast<int>(dim);
        constexpr int int_spacetime_dim = static_cast<int>(spacetime_dim);
        const auto gamma0_at_s = make_at_index(device_constraint_gamma0, point);
        const auto gamma1_at_s = make_at_index(device_constraint_gamma1, point);
        const auto gamma2_at_s = make_at_index(device_constraint_gamma2, point);
        const double gamma0 = get(gamma0_at_s);
        const double gamma1 = get(gamma1_at_s);
        const double gamma2 = get(gamma2_at_s);

        double* const shmem = static_cast<double*>(
            team_member.team_shmem().get_shmem(dt_pi_scratch_size_bytes));
        double* const shift = shmem + dt_pi_offset_shift;
        double* const inverse_spatial_metric =
            shmem + dt_pi_offset_inverse_spatial_metric;
        double* const inverse_spacetime_metric =
            shmem + dt_pi_offset_inverse_spacetime_metric;
        double* const pi_up = shmem + dt_pi_offset_pi_up;
        double* const pi_one_normal = shmem + dt_pi_offset_pi_one_normal;
        double* const gauge_constraint = shmem + dt_pi_offset_gauge_constraint;
        double* const gauge_h_components = shmem + dt_pi_offset_gauge_h;
        double* const spacetime_metric_components =
            shmem + dt_pi_offset_spacetime_metric_components;
        double* const pi_components = shmem + dt_pi_offset_pi_components;
        double* const phi_components = shmem + dt_pi_offset_phi_components;

        ::Kokkos::single(::Kokkos::PerTeam(team_member), [&]() {
          for (int mu = 0; mu < int_spacetime_dim; ++mu) {
            gauge_h_components[mu] = gauge_h.get(mu)[point];
            for (int nu = mu; nu < int_spacetime_dim; ++nu) {
              const int symmetric_index =
                  static_cast<int>(symmetric_spacetime_index(mu, nu));
              spacetime_metric_components[symmetric_index] =
                  spacetime_metric.get(mu, nu)[point];
              pi_components[symmetric_index] = pi.get(mu, nu)[point];
              for (int i = 0; i < int_dim; ++i) {
                phi_components[i * spacetime_symmetric_size + symmetric_index] =
                    phi.get(i, mu, nu)[point];
              }
            }
          }

          const double g11 = get_symmetric_spacetime_component(
              spacetime_metric_components, 1, 1);
          const double g12 = get_symmetric_spacetime_component(
              spacetime_metric_components, 1, 2);
          const double g13 = get_symmetric_spacetime_component(
              spacetime_metric_components, 1, 3);
          const double g22 = get_symmetric_spacetime_component(
              spacetime_metric_components, 2, 2);
          const double g23 = get_symmetric_spacetime_component(
              spacetime_metric_components, 2, 3);
          const double g33 = get_symmetric_spacetime_component(
              spacetime_metric_components, 3, 3);
          const double det_spatial_metric = g11 * (g22 * g33 - g23 * g23) -
                                            g12 * (g12 * g33 - g23 * g13) +
                                            g13 * (g12 * g23 - g22 * g13);
          const double inv_det = 1.0 / det_spatial_metric;

          inverse_spatial_metric[0 * dim + 0] =
              (g22 * g33 - g23 * g23) * inv_det;
          inverse_spatial_metric[0 * dim + 1] =
              (g13 * g23 - g12 * g33) * inv_det;
          inverse_spatial_metric[0 * dim + 2] =
              (g12 * g23 - g13 * g22) * inv_det;
          inverse_spatial_metric[1 * dim + 0] =
              inverse_spatial_metric[0 * dim + 1];
          inverse_spatial_metric[1 * dim + 1] =
              (g11 * g33 - g13 * g13) * inv_det;
          inverse_spatial_metric[1 * dim + 2] =
              (g13 * g12 - g11 * g23) * inv_det;
          inverse_spatial_metric[2 * dim + 0] =
              inverse_spatial_metric[0 * dim + 2];
          inverse_spatial_metric[2 * dim + 1] =
              inverse_spatial_metric[1 * dim + 2];
          inverse_spatial_metric[2 * dim + 2] =
              (g11 * g22 - g12 * g12) * inv_det;

          for (int i = 0; i < int_dim; ++i) {
            shift[i] = 0.0;
            for (int j = 0; j < int_dim; ++j) {
              shift[i] += inverse_spatial_metric[i * dim + j] *
                          get_symmetric_spacetime_component(
                              spacetime_metric_components, 0, j + 1);
            }
          }

          double lapse_squared = -get_symmetric_spacetime_component(
              spacetime_metric_components, 0, 0);
          for (int i = 0; i < int_dim; ++i) {
            lapse_squared +=
                shift[i] * get_symmetric_spacetime_component(
                               spacetime_metric_components, 0, i + 1);
          }
          const double lapse = sqrt(lapse_squared);
          const double inv_lapse_squared = 1.0 / lapse_squared;
          shmem[dt_pi_offset_lapse] = lapse;

          for (int a = 0; a < int_spacetime_dim; ++a) {
            for (int b = 0; b < int_spacetime_dim; ++b) {
              inverse_spacetime_metric[a * spacetime_dim + b] = 0.0;
            }
          }
          inverse_spacetime_metric[0 * spacetime_dim + 0] = -inv_lapse_squared;
          for (int i = 0; i < int_dim; ++i) {
            const double zero_i_component = shift[i] * inv_lapse_squared;
            inverse_spacetime_metric[0 * spacetime_dim + (i + 1)] =
                zero_i_component;
            inverse_spacetime_metric[(i + 1) * spacetime_dim + 0] =
                zero_i_component;
            for (int j = i; j < int_dim; ++j) {
              const double ii_component =
                  inverse_spatial_metric[i * dim + j] -
                  shift[i] * shift[j] * inv_lapse_squared;
              inverse_spacetime_metric[(i + 1) * spacetime_dim + (j + 1)] =
                  ii_component;
              inverse_spacetime_metric[(j + 1) * spacetime_dim + (i + 1)] =
                  ii_component;
            }
          }

          for (int nu = 0; nu < int_spacetime_dim; ++nu) {
            for (int delta = 0; delta < int_spacetime_dim; ++delta) {
              pi_up[nu * spacetime_dim + delta] = 0.0;
              for (int beta = 0; beta < int_spacetime_dim; ++beta) {
                pi_up[nu * spacetime_dim + delta] +=
                    inverse_spacetime_metric[delta * spacetime_dim + beta] *
                    get_symmetric_spacetime_component(pi_components, nu, beta);
              }
            }
          }

          double normal_spacetime_vector[spacetime_dim];
          normal_spacetime_vector[0] = 1.0 / lapse;
          for (int i = 0; i < int_dim; ++i) {
            normal_spacetime_vector[i + 1] = -shift[i] / lapse;
          }

          for (int mu = 0; mu < int_spacetime_dim; ++mu) {
            pi_one_normal[mu] = 0.0;
            for (int nu = 0; nu < int_spacetime_dim; ++nu) {
              pi_one_normal[mu] +=
                  normal_spacetime_vector[nu] *
                  get_symmetric_spacetime_component(pi_components, nu, mu);
            }
          }

          double half_pi_two_normals = 0.0;
          for (int mu = 0; mu < int_spacetime_dim; ++mu) {
            half_pi_two_normals +=
                normal_spacetime_vector[mu] * pi_one_normal[mu];
          }
          shmem[dt_pi_offset_half_pi_two_normals] = 0.5 * half_pi_two_normals;

          for (int a = 0; a < int_spacetime_dim; ++a) {
            double trace_christoffel = 0.0;
            for (int b = 0; b < int_spacetime_dim; ++b) {
              for (int c = 0; c < int_spacetime_dim; ++c) {
                trace_christoffel +=
                    christoffel_first_kind_component_from_scratch(
                        a, b, c, lapse, shift, pi_components, phi_components) *
                    inverse_spacetime_metric[b * spacetime_dim + c];
              }
            }
            gauge_constraint[a] = trace_christoffel + gauge_h_components[a];
          }

          double normal_dot_gauge_constraint =
              normal_spacetime_vector[0] * gauge_constraint[0];
          for (int mu = 1; mu < int_spacetime_dim; ++mu) {
            normal_dot_gauge_constraint +=
                normal_spacetime_vector[mu] * gauge_constraint[mu];
          }

          shmem[dt_pi_offset_gamma1gamma2] = gamma1 * gamma2;
          shmem[dt_pi_offset_minus_gamma0_lapse] = -gamma0 * lapse;
          shmem[dt_pi_offset_normal_dot_gauge_constraint] =
              gamma0 * normal_dot_gauge_constraint;
        });

        team_member.team_barrier();

        ::Kokkos::parallel_for(
            ::Kokkos::TeamThreadRange(
                team_member, static_cast<int>(spacetime_symmetric_size)),
            [&](const int component_int) {
              const int component = component_int;
              const auto mu_nu =
                  symmetric_spacetime_indices(static_cast<size_t>(component));
              const int mu = static_cast<int>(mu_nu[0]);
              const int nu = static_cast<int>(mu_nu[1]);

              const double lapse = shmem[dt_pi_offset_lapse];
              const double half_pi_two_normals =
                  shmem[dt_pi_offset_half_pi_two_normals];
              const double gamma1gamma2 = shmem[dt_pi_offset_gamma1gamma2];
              const double minus_gamma0_lapse =
                  shmem[dt_pi_offset_minus_gamma0_lapse];
              const double normal_dot_gauge_constraint =
                  shmem[dt_pi_offset_normal_dot_gauge_constraint];

              double dt_pi_mu_nu = -normal_dot_gauge_constraint *
                                   get_symmetric_spacetime_component(
                                       spacetime_metric_components, mu, nu);
              if (mu == 0) {
                dt_pi_mu_nu = minus_gamma0_lapse * gauge_constraint[nu] -
                              normal_dot_gauge_constraint *
                                  get_symmetric_spacetime_component(
                                      spacetime_metric_components, 0, nu);
                if (nu == 0) {
                  dt_pi_mu_nu = 2.0 * minus_gamma0_lapse * gauge_constraint[0] -
                                normal_dot_gauge_constraint *
                                    get_symmetric_spacetime_component(
                                        spacetime_metric_components, 0, 0);
                }
              }
              dt_pi_mu_nu -=
                  half_pi_two_normals *
                  get_symmetric_spacetime_component(pi_components, mu, nu);
              dt_pi_mu_nu -= spacetime_deriv_gauge_h.get(mu, nu)[point] +
                             spacetime_deriv_gauge_h.get(nu, mu)[point];

              for (int delta = 0; delta < int_spacetime_dim; ++delta) {
                dt_pi_mu_nu -= 2.0 *
                               get_symmetric_spacetime_component(pi_components,
                                                                 mu, delta) *
                               pi_up[nu * spacetime_dim + delta];

                double christoffel_second_kind_delta_mu_nu = 0.0;
                for (int alpha = 0; alpha < int_spacetime_dim; ++alpha) {
                  christoffel_second_kind_delta_mu_nu +=
                      inverse_spacetime_metric[delta * spacetime_dim + alpha] *
                      christoffel_first_kind_component_from_scratch(
                          alpha, mu, nu, lapse, shift, pi_components,
                          phi_components);
                }
                dt_pi_mu_nu += 2.0 * christoffel_second_kind_delta_mu_nu *
                               gauge_h_components[delta];

                for (int n = 0; n < int_dim; ++n) {
                  double phi_1_up_n_mu_delta = 0.0;
                  for (int m = 0; m < int_dim; ++m) {
                    phi_1_up_n_mu_delta +=
                        inverse_spatial_metric[n * dim + m] *
                        get_phi_component(phi_components, m, mu, delta);
                  }

                  double phi_3_up_n_nu_delta = 0.0;
                  for (int beta = 0; beta < int_spacetime_dim; ++beta) {
                    phi_3_up_n_nu_delta +=
                        inverse_spacetime_metric[delta * spacetime_dim + beta] *
                        get_phi_component(phi_components, n, nu, beta);
                  }

                  dt_pi_mu_nu +=
                      2.0 * phi_1_up_n_mu_delta * phi_3_up_n_nu_delta;
                }

                for (int alpha = 0; alpha < int_spacetime_dim; ++alpha) {
                  double christoffel_first_kind_3_up_mu_alpha_delta = 0.0;
                  double christoffel_first_kind_3_up_nu_delta_alpha = 0.0;
                  for (int beta = 0; beta < int_spacetime_dim; ++beta) {
                    christoffel_first_kind_3_up_mu_alpha_delta +=
                        inverse_spacetime_metric[alpha * spacetime_dim + beta] *
                        christoffel_first_kind_component_from_scratch(
                            mu, delta, beta, lapse, shift, pi_components,
                            phi_components);
                    christoffel_first_kind_3_up_nu_delta_alpha +=
                        inverse_spacetime_metric[delta * spacetime_dim + beta] *
                        christoffel_first_kind_component_from_scratch(
                            nu, alpha, beta, lapse, shift, pi_components,
                            phi_components);
                  }
                  dt_pi_mu_nu -= 2.0 *
                                 christoffel_first_kind_3_up_mu_alpha_delta *
                                 christoffel_first_kind_3_up_nu_delta_alpha;
                }
              }

              for (int m = 0; m < int_dim; ++m) {
                double phi_1_up_m_mu_nu = 0.0;
                for (int n = 0; n < int_dim; ++n) {
                  phi_1_up_m_mu_nu +=
                      inverse_spatial_metric[m * dim + n] *
                      get_phi_component(phi_components, n, mu, nu);
                }
                dt_pi_mu_nu -= pi_one_normal[m + 1] * phi_1_up_m_mu_nu;
                for (int n = 0; n < int_dim; ++n) {
                  dt_pi_mu_nu -= inverse_spatial_metric[m * dim + n] *
                                 d_phi.get(m, n, mu, nu)[point];
                }
              }

              dt_pi_mu_nu *= lapse;
              double shift_dot_three_index_constraint_mu_nu = 0.0;
              for (int m = 0; m < int_dim; ++m) {
                shift_dot_three_index_constraint_mu_nu +=
                    shift[m] * (d_spacetime_metric.get(m, mu, nu)[point] -
                                get_phi_component(phi_components, m, mu, nu));
              }
              dt_pi_mu_nu +=
                  gamma1gamma2 * shift_dot_three_index_constraint_mu_nu;
              for (int m = 0; m < int_dim; ++m) {
                dt_pi_mu_nu += shift[m] * d_pi.get(m, mu, nu)[point];
              }
              dt_pi.get(mu, nu)[point] = dt_pi_mu_nu;
            });
      });

  ::Kokkos::parallel_for(
      "GhBatchedComputeTimeDerivativeDtPhi", total_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        constexpr int int_dim = static_cast<int>(dim);
        constexpr int int_dim_plus_one = static_cast<int>(dim + 1);
        const auto gamma2_at_s = make_at_index(device_constraint_gamma2, s);
        const double gamma2 = get(gamma2_at_s);

        const auto spacetime_metric_at_s = make_at_index(spacetime_metric, s);
        const auto pi_at_s = make_at_index(pi, s);
        const auto phi_at_s = make_at_index(phi, s);

        tnsr::II<double, dim, Frame::Inertial> inverse_spatial_metric{};
        tnsr::I<double, dim, Frame::Inertial> shift{};
        double lapse = 0.0;
        double inv_lapse_squared = 0.0;
        compute_inverse_spatial_metric_shift_lapse(
            make_not_null(&inverse_spatial_metric), make_not_null(&shift),
            make_not_null(&lapse), make_not_null(&inv_lapse_squared),
            spacetime_metric_at_s);
        (void)inv_lapse_squared;

        tnsr::A<double, dim, Frame::Inertial> normal_spacetime_vector{};
        normal_spacetime_vector.get(0) = 1.0 / lapse;
        for (int i = 0; i < int_dim; ++i) {
          normal_spacetime_vector.get(i + 1) = -shift.get(i) / lapse;
        }

        for (int i = 0; i < int_dim; ++i) {
          double phi_one_normal_i[4];
          for (int a = 0; a < int_dim_plus_one; ++a) {
            phi_one_normal_i[a] = 0.0;
            for (int b = 0; b < int_dim_plus_one; ++b) {
              phi_one_normal_i[a] +=
                  normal_spacetime_vector.get(b) * phi_at_s.get(i, b, a);
            }
          }
          double half_phi_two_normals_i = 0.0;
          for (int a = 0; a < int_dim_plus_one; ++a) {
            half_phi_two_normals_i +=
                normal_spacetime_vector.get(a) * phi_one_normal_i[a];
          }
          half_phi_two_normals_i *= 0.5;
          for (int mu = 0; mu < int_dim_plus_one; ++mu) {
            for (int nu = mu; nu < int_dim_plus_one; ++nu) {
              double dt_phi_i_mu_nu =
                  pi_at_s.get(mu, nu) * half_phi_two_normals_i -
                  d_pi.get(i, mu, nu)[point] +
                  gamma2 * (d_spacetime_metric.get(i, mu, nu)[point] -
                            phi_at_s.get(i, mu, nu));

              for (int n = 0; n < int_dim; ++n) {
                double phi_1_up_n_mu_nu = 0.0;
                for (int m = 0; m < int_dim; ++m) {
                  phi_1_up_n_mu_nu += inverse_spatial_metric.get(n, m) *
                                      phi_at_s.get(m, mu, nu);
                }
                dt_phi_i_mu_nu += phi_one_normal_i[n + 1] * phi_1_up_n_mu_nu;
              }

              dt_phi_i_mu_nu *= lapse;
              for (int m = 0; m < int_dim; ++m) {
                dt_phi_i_mu_nu += shift.get(m) * d_phi.get(m, i, mu, nu)[point];
              }
              dt_phi.get(i, mu, nu)[point] = dt_phi_i_mu_nu;
            }
          }
        }
      });
}

void ComputeTimeDerivativeBatched::apply(
    const gsl::not_null<typename packed_evolution_state_tag::type*>
        packed_evolution_state,
    const gsl::not_null<typename packed_boundary_scratch_tag::type*>
        packed_boundary_scratch,
    const typename packed_topology_tag::type& packed_topology,
    const typename packed_geometry_tag::type& packed_geometry,
    const typename device_constraint_gamma0_tag::type& device_gamma0,
    const typename device_constraint_gamma1_tag::type& device_gamma1,
    const typename device_constraint_gamma2_tag::type& device_gamma2) {
  compute_time_derivative_batched_volume_impl(
      packed_evolution_state, packed_boundary_scratch, device_gamma0,
      device_gamma1, device_gamma2, packed_topology, packed_geometry);
}

}  // namespace gh::Actions
