// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/ComputeTimeDerivativeKokkos.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {
constexpr size_t dim = 3;
using device_gauge_h_tag =
    ::Tags::MirrorView<gh::Tags::GaugeH<DataVector, dim>>;
using device_spacetime_deriv_gauge_h_tag =
    ::Tags::MirrorView<gh::Tags::SpacetimeDerivGaugeH<DataVector, dim>>;
using device_gauge_data_type = Variables<
    tmpl::list<device_gauge_h_tag, device_spacetime_deriv_gauge_h_tag>>;

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

  // Match gh::trace_christoffel algebra used by host AnalyticChristoffel.
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

KOKKOS_INLINE_FUNCTION void compute_rhs_at_point(
    const gsl::not_null<tnsr::aa<double, dim, Frame::Inertial>*>
        dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<double, dim, Frame::Inertial>*> dt_pi,
    const gsl::not_null<tnsr::iaa<double, dim, Frame::Inertial>*> dt_phi,
    const tnsr::aa<double, dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<double, dim, Frame::Inertial>& pi,
    const tnsr::iaa<double, dim, Frame::Inertial>& phi,
    const tnsr::iaa<double, dim, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<double, dim, Frame::Inertial>& d_pi,
    const tnsr::ijaa<double, dim, Frame::Inertial>& d_phi,
    const tnsr::a<double, dim, Frame::Inertial>& gauge_h,
    const tnsr::ab<double, dim, Frame::Inertial>& spacetime_deriv_gauge_h,
    const double gamma0, const double gamma1, const double gamma2) {
  tnsr::II<double, dim, Frame::Inertial> inverse_spatial_metric{};
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det(make_not_null(&inverse_spatial_metric),
                                 make_not_null(&det_spatial_metric),
                                 spacetime_metric);

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

  // First piece of dt(spacetime_metric), before adding constraints.
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      dt_spacetime_metric->get(mu, nu) = -lapse * pi.get(mu, nu);
      for (size_t m = 0; m < dim; ++m) {
        dt_spacetime_metric->get(mu, nu) += shift.get(m) * phi.get(m, mu, nu);
      }
    }
  }

  tnsr::abb<double, dim, Frame::Inertial> da_spacetime_metric{};
  for (size_t a = 0; a < dim + 1; ++a) {
    for (size_t b = a; b < dim + 1; ++b) {
      da_spacetime_metric.get(0, a, b) = dt_spacetime_metric->get(a, b);
      for (size_t i = 0; i < dim; ++i) {
        da_spacetime_metric.get(i + 1, a, b) = phi.get(i, a, b);
      }
    }
  }

  tnsr::abb<double, dim, Frame::Inertial> christoffel_first_kind{};
  for (size_t k = 0; k < dim + 1; ++k) {
    for (size_t i = 0; i < dim + 1; ++i) {
      for (size_t j = i; j < dim + 1; ++j) {
        christoffel_first_kind.get(k, i, j) =
            0.5 * (da_spacetime_metric.get(i, j, k) +
                   da_spacetime_metric.get(j, i, k) -
                   da_spacetime_metric.get(k, i, j));
      }
    }
  }

  tnsr::a<double, dim, Frame::Inertial> trace_christoffel{};
  for (size_t a = 0; a < dim + 1; ++a) {
    trace_christoffel.get(a) = 0.0;
    for (size_t b = 0; b < dim + 1; ++b) {
      for (size_t c = 0; c < dim + 1; ++c) {
        trace_christoffel.get(a) += christoffel_first_kind.get(a, b, c) *
                                    inverse_spacetime_metric.get(b, c);
      }
    }
  }

  tnsr::A<double, dim, Frame::Inertial> normal_spacetime_vector{};
  normal_spacetime_vector.get(0) = 1.0 / lapse;
  for (size_t i = 0; i < dim; ++i) {
    normal_spacetime_vector.get(i + 1) = -shift.get(i) / lapse;
  }

  const double gamma1gamma2 = gamma1 * gamma2;
  const double gamma1_plus_1 = 1.0 + gamma1;

  tnsr::Iaa<double, dim, Frame::Inertial> phi_1_up{};
  for (size_t m = 0; m < dim; ++m) {
    for (size_t mu = 0; mu < dim + 1; ++mu) {
      for (size_t nu = mu; nu < dim + 1; ++nu) {
        phi_1_up.get(m, mu, nu) = 0.0;
        for (size_t n = 0; n < dim; ++n) {
          phi_1_up.get(m, mu, nu) +=
              inverse_spatial_metric.get(m, n) * phi.get(n, mu, nu);
        }
      }
    }
  }

  tnsr::iaB<double, dim, Frame::Inertial> phi_3_up{};
  for (size_t m = 0; m < dim; ++m) {
    for (size_t nu = 0; nu < dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < dim + 1; ++alpha) {
        phi_3_up.get(m, nu, alpha) = 0.0;
        for (size_t beta = 0; beta < dim + 1; ++beta) {
          phi_3_up.get(m, nu, alpha) +=
              inverse_spacetime_metric.get(alpha, beta) * phi.get(m, nu, beta);
        }
      }
    }
  }

  tnsr::aB<double, dim, Frame::Inertial> pi_2_up{};
  for (size_t nu = 0; nu < dim + 1; ++nu) {
    for (size_t alpha = 0; alpha < dim + 1; ++alpha) {
      pi_2_up.get(nu, alpha) = 0.0;
      for (size_t beta = 0; beta < dim + 1; ++beta) {
        pi_2_up.get(nu, alpha) +=
            inverse_spacetime_metric.get(alpha, beta) * pi.get(nu, beta);
      }
    }
  }

  tnsr::abC<double, dim, Frame::Inertial> christoffel_first_kind_3_up{};
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    for (size_t nu = 0; nu < dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < dim + 1; ++alpha) {
        christoffel_first_kind_3_up.get(mu, nu, alpha) = 0.0;
        for (size_t beta = 0; beta < dim + 1; ++beta) {
          christoffel_first_kind_3_up.get(mu, nu, alpha) +=
              inverse_spacetime_metric.get(alpha, beta) *
              christoffel_first_kind.get(mu, nu, beta);
        }
      }
    }
  }

  tnsr::Abb<double, dim, Frame::Inertial> christoffel_second_kind{};
  for (size_t delta = 0; delta < dim + 1; ++delta) {
    for (size_t mu = 0; mu < dim + 1; ++mu) {
      for (size_t nu = mu; nu < dim + 1; ++nu) {
        christoffel_second_kind.get(delta, mu, nu) = 0.0;
        for (size_t alpha = 0; alpha < dim + 1; ++alpha) {
          christoffel_second_kind.get(delta, mu, nu) +=
              inverse_spacetime_metric.get(delta, alpha) *
              christoffel_first_kind.get(alpha, mu, nu);
        }
      }
    }
  }

  tnsr::a<double, dim, Frame::Inertial> pi_one_normal{};
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    pi_one_normal.get(mu) = 0.0;
    for (size_t nu = 0; nu < dim + 1; ++nu) {
      pi_one_normal.get(mu) += normal_spacetime_vector.get(nu) * pi.get(nu, mu);
    }
  }

  double half_pi_two_normals = 0.0;
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    half_pi_two_normals +=
        normal_spacetime_vector.get(mu) * pi_one_normal.get(mu);
  }
  half_pi_two_normals *= 0.5;

  tnsr::ia<double, dim, Frame::Inertial> phi_one_normal{};
  for (size_t n = 0; n < dim; ++n) {
    for (size_t nu = 0; nu < dim + 1; ++nu) {
      phi_one_normal.get(n, nu) = 0.0;
      for (size_t mu = 0; mu < dim + 1; ++mu) {
        phi_one_normal.get(n, nu) +=
            normal_spacetime_vector.get(mu) * phi.get(n, mu, nu);
      }
    }
  }

  tnsr::i<double, dim, Frame::Inertial> half_phi_two_normals{};
  for (size_t n = 0; n < dim; ++n) {
    half_phi_two_normals.get(n) = 0.0;
    for (size_t mu = 0; mu < dim + 1; ++mu) {
      half_phi_two_normals.get(n) +=
          normal_spacetime_vector.get(mu) * phi_one_normal.get(n, mu);
    }
    half_phi_two_normals.get(n) *= 0.5;
  }

  tnsr::iaa<double, dim, Frame::Inertial> three_index_constraint{};
  for (size_t n = 0; n < dim; ++n) {
    for (size_t mu = 0; mu < dim + 1; ++mu) {
      for (size_t nu = mu; nu < dim + 1; ++nu) {
        three_index_constraint.get(n, mu, nu) =
            d_spacetime_metric.get(n, mu, nu) - phi.get(n, mu, nu);
      }
    }
  }

  tnsr::a<double, dim, Frame::Inertial> gauge_constraint{};
  tnsr::aa<double, dim, Frame::Inertial> shift_dot_three_index_constraint{};
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    gauge_constraint.get(mu) = trace_christoffel.get(mu) + gauge_h.get(mu);
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      shift_dot_three_index_constraint.get(mu, nu) =
          shift.get(0) * three_index_constraint.get(0, mu, nu);
      for (size_t m = 1; m < dim; ++m) {
        shift_dot_three_index_constraint.get(mu, nu) +=
            shift.get(m) * three_index_constraint.get(m, mu, nu);
      }
    }
  }

  double normal_dot_gauge_constraint =
      normal_spacetime_vector.get(0) * gauge_constraint.get(0);
  for (size_t mu = 1; mu < dim + 1; ++mu) {
    normal_dot_gauge_constraint +=
        normal_spacetime_vector.get(mu) * gauge_constraint.get(mu);
  }

  // Full equation for dt(spacetime_metric)
  for (size_t mu = 0; mu < dim + 1; ++mu) {
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      dt_spacetime_metric->get(mu, nu) +=
          gamma1_plus_1 * shift_dot_three_index_constraint.get(mu, nu);
    }
  }

  // Equation for dt(pi)
  normal_dot_gauge_constraint *= gamma0;
  const double minus_gamma0_lapse = -gamma0 * lapse;

  for (size_t i = 1; i < dim + 1; ++i) {
    dt_pi->get(0, i) = minus_gamma0_lapse * gauge_constraint.get(i) -
                       normal_dot_gauge_constraint * spacetime_metric.get(0, i);
  }
  dt_pi->get(0, 0) = 2.0 * minus_gamma0_lapse * gauge_constraint.get(0) -
                     normal_dot_gauge_constraint * spacetime_metric.get(0, 0);
  for (size_t mu = 1; mu < dim + 1; ++mu) {
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      dt_pi->get(mu, nu) =
          -normal_dot_gauge_constraint * spacetime_metric.get(mu, nu);
    }
  }

  for (size_t mu = 0; mu < dim + 1; ++mu) {
    for (size_t nu = mu; nu < dim + 1; ++nu) {
      dt_pi->get(mu, nu) -= half_pi_two_normals * pi.get(mu, nu);
      dt_pi->get(mu, nu) -= spacetime_deriv_gauge_h.get(mu, nu) +
                            spacetime_deriv_gauge_h.get(nu, mu);

      for (size_t delta = 0; delta < dim + 1; ++delta) {
        dt_pi->get(mu, nu) -= 2.0 * pi.get(mu, delta) * pi_2_up.get(nu, delta);
        dt_pi->get(mu, nu) += 2.0 * christoffel_second_kind.get(delta, mu, nu) *
                              gauge_h.get(delta);

        for (size_t n = 0; n < dim; ++n) {
          dt_pi->get(mu, nu) +=
              2.0 * phi_1_up.get(n, mu, delta) * phi_3_up.get(n, nu, delta);
        }

        for (size_t alpha = 0; alpha < dim + 1; ++alpha) {
          dt_pi->get(mu, nu) -=
              2.0 * christoffel_first_kind_3_up.get(mu, alpha, delta) *
              christoffel_first_kind_3_up.get(nu, delta, alpha);
        }
      }

      for (size_t m = 0; m < dim; ++m) {
        dt_pi->get(mu, nu) -=
            pi_one_normal.get(m + 1) * phi_1_up.get(m, mu, nu);
        for (size_t n = 0; n < dim; ++n) {
          dt_pi->get(mu, nu) -=
              inverse_spatial_metric.get(m, n) * d_phi.get(m, n, mu, nu);
        }
      }

      dt_pi->get(mu, nu) *= lapse;
      dt_pi->get(mu, nu) +=
          gamma1gamma2 * shift_dot_three_index_constraint.get(mu, nu);
      for (size_t m = 0; m < dim; ++m) {
        dt_pi->get(mu, nu) += shift.get(m) * d_pi.get(m, mu, nu);
      }
    }
  }

  // Equation for dt(phi)
  for (size_t i = 0; i < dim; ++i) {
    for (size_t mu = 0; mu < dim + 1; ++mu) {
      for (size_t nu = mu; nu < dim + 1; ++nu) {
        dt_phi->get(i, mu, nu) = pi.get(mu, nu) * half_phi_two_normals.get(i) -
                                 d_pi.get(i, mu, nu) +
                                 gamma2 * three_index_constraint.get(i, mu, nu);

        for (size_t n = 0; n < dim; ++n) {
          dt_phi->get(i, mu, nu) +=
              phi_one_normal.get(i, n + 1) * phi_1_up.get(n, mu, nu);
        }

        dt_phi->get(i, mu, nu) *= lapse;
        for (size_t m = 0; m < dim; ++m) {
          dt_phi->get(i, mu, nu) += shift.get(m) * d_phi.get(m, i, mu, nu);
        }
      }
    }
  }
}
}  // namespace

void ComputeTimeDerivativeKokkos::orient_boundary_data_for_send(
    const gsl::not_null<Variables<device_package_field_tags>*>
        oriented_boundary_data,
    const Variables<device_package_field_tags>& boundary_data,
    const ::Kokkos::View<size_t*>& oriented_mortar_grid_point_source_index) {
  (void)oriented_boundary_data;
  (void)boundary_data;
  (void)oriented_mortar_grid_point_source_index;
  ERROR("TODO: Implement GH Kokkos boundary-data orientation for send.");
}

void compute_hardcoded_analytic_gauge_and_spacetime_derivative(
    const gsl::not_null<device_gauge_data_type*> device_gauge_data,
    const gh::KokkosTags::DeviceInertialCoordinates<dim>::type&
        device_inertial_coordinates,
    const gh::KokkosTags::DeviceInverseJacobian<dim>::type&
        device_inverse_jacobian,
    const Mesh<dim>& mesh) {
  const size_t number_of_points = mesh.number_of_grid_points();
  if (device_gauge_data->number_of_grid_points() != number_of_points) {
    device_gauge_data->initialize(number_of_points);
  }

  using device_spatial_deriv_gauge_h_tag =
      ::Tags::deriv<device_gauge_h_tag, tmpl::size_t<dim>, Frame::Inertial>;

  const auto device_gauge_h = get<device_gauge_h_tag>(*device_gauge_data);
  ::Kokkos::parallel_for(
      "GhComputeTimeDerivativeKokkosGaugeH", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        tnsr::I<double, dim, Frame::Inertial> inertial_coords_at_s{};
        for (size_t d = 0; d < dim; ++d) {
          inertial_coords_at_s.get(d) =
              device_inertial_coordinates.get(d)[point];
        }
        tnsr::a<double, dim, Frame::Inertial> gauge_h_at_s{};
        compute_hardcoded_analytic_gauge(make_not_null(&gauge_h_at_s),
                                         inertial_coords_at_s);
        for (size_t a = 0; a < dim + 1; ++a) {
          device_gauge_h.get(a)[point] = gauge_h_at_s.get(a);
        }
      });

  Variables<tmpl::list<device_gauge_h_tag>> device_gauge_h_vars{
      number_of_points};
  const auto gauge_h_src = get<device_gauge_h_tag>(*device_gauge_data);
  const auto gauge_h_dst = get<device_gauge_h_tag>(device_gauge_h_vars);
  ::Kokkos::parallel_for(
      "GhComputeTimeDerivativeKokkosCopyGaugeH", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        for (size_t a = 0; a < dim + 1; ++a) {
          gauge_h_dst.get(a)[point] = gauge_h_src.get(a)[point];
        }
      });

  Variables<tmpl::list<device_spatial_deriv_gauge_h_tag>>
      device_spatial_gauge_deriv{number_of_points};
  partial_derivatives(make_not_null(&device_spatial_gauge_deriv),
                      device_gauge_h_vars, mesh, device_inverse_jacobian);

  const auto device_spatial_deriv_gauge_h =
      get<device_spatial_deriv_gauge_h_tag>(device_spatial_gauge_deriv);
  const auto device_spacetime_deriv_gauge_h =
      get<device_spacetime_deriv_gauge_h_tag>(*device_gauge_data);
  ::Kokkos::parallel_for(
      "GhComputeTimeDerivativeKokkosDerivGaugeH", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        const size_t point = static_cast<size_t>(s);
        for (size_t a = 0; a < dim + 1; ++a) {
          device_spacetime_deriv_gauge_h.get(0, a)[point] = 0.0;
          for (size_t i = 0; i < dim; ++i) {
            device_spacetime_deriv_gauge_h.get(i + 1, a)[point] =
                device_spatial_deriv_gauge_h.get(i, a)[point];
          }
        }
      });
}

void ComputeTimeDerivativeKokkos::
    compute_volume_terms_and_package_boundary_data(
        const gsl::not_null<device_dt_type*> device_dt,
        const gsl::not_null<outgoing_boundary_data_type*>
            outgoing_boundary_data,
        const gsl::not_null<external_boundary_data_type*>
            external_boundary_data,
        const device_variables_type& device_vars,
        const device_inverse_jacobian_type& device_inverse_jacobian,
        const device_inertial_coordinates_type& device_inertial_coordinates,
        const device_constraint_gamma0_type& device_constraint_gamma0,
        const device_constraint_gamma1_type& device_constraint_gamma1,
        const device_constraint_gamma2_type& device_constraint_gamma2,
        const device_face_to_volume_index_map_type&
            device_face_to_volume_index_map,
        const device_face_unit_normal_covector_type&
            device_face_unit_normal_covector,
        const device_mortar_data_type& device_mortar_data,
        const typename mortar_mesh_tag::type& mortar_meshes,
        const gh::Tags::ConstraintGamma0::type& host_constraint_gamma0,
        const gh::Tags::ConstraintGamma1::type& host_constraint_gamma1,
        const gh::Tags::ConstraintGamma2::type& host_constraint_gamma2,
        const typename domain::Tags::ExternalBoundaryConditions<
            volume_dim>::type& external_boundary_conditions_by_block,
        const double time, const Mesh<volume_dim>& mesh,
        const Element<volume_dim>& element, const TimeStepId& time_step_id) {
  (void)device_face_to_volume_index_map;
  (void)device_face_unit_normal_covector;
  (void)device_mortar_data;
  (void)mortar_meshes;
  (void)host_constraint_gamma0;
  (void)host_constraint_gamma1;
  (void)host_constraint_gamma2;
  (void)external_boundary_conditions_by_block;
  (void)time;
  (void)element;
  (void)time_step_id;

  const size_t number_of_points = mesh.number_of_grid_points();
  if (device_dt->number_of_grid_points() != number_of_points) {
    device_dt->initialize(number_of_points);
  }

  using spacetime_metric_tag =
      gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>;
  using pi_tag = gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>;
  using phi_tag = gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>;

  using device_gradient_tags =
      db::wrap_tags_in<::Tags::MirrorView, typename system::gradient_variables>;
  using device_derivative_tags =
      db::wrap_tags_in<::Tags::deriv, device_gradient_tags,
                       tmpl::size_t<volume_dim>, Frame::Inertial>;
  Variables<device_derivative_tags> device_partial_derivatives{
      number_of_points};
  partial_derivatives(make_not_null(&device_partial_derivatives), device_vars,
                      mesh, device_inverse_jacobian);

  device_gauge_data_type device_gauge_data{number_of_points};
  compute_hardcoded_analytic_gauge_and_spacetime_derivative(
      make_not_null(&device_gauge_data), device_inertial_coordinates,
      device_inverse_jacobian, mesh);

  const auto dt_spacetime_metric =
      get<::Tags::MirrorView<::Tags::dt<spacetime_metric_tag>>>(*device_dt);
  const auto dt_pi = get<::Tags::MirrorView<::Tags::dt<pi_tag>>>(*device_dt);
  const auto dt_phi = get<::Tags::MirrorView<::Tags::dt<phi_tag>>>(*device_dt);

  ::Kokkos::parallel_for(
      "GhComputeTimeDerivativeKokkosVolumeTerms", number_of_points,
      KOKKOS_LAMBDA(const int s) {
        tnsr::aa<double, volume_dim, Frame::Inertial>
            dt_spacetime_metric_at_s{};
        tnsr::aa<double, volume_dim, Frame::Inertial> dt_pi_at_s{};
        tnsr::iaa<double, volume_dim, Frame::Inertial> dt_phi_at_s{};

        const auto vars_at_s = make_at_index(device_vars, s);
        const auto derivs_at_s = make_at_index(device_partial_derivatives, s);
        const auto gauge_data_at_s = make_at_index(device_gauge_data, s);
        const auto gamma0_at_s = make_at_index(device_constraint_gamma0, s);
        const auto gamma1_at_s = make_at_index(device_constraint_gamma1, s);
        const auto gamma2_at_s = make_at_index(device_constraint_gamma2, s);

        compute_rhs_at_point(
            make_not_null(&dt_spacetime_metric_at_s),
            make_not_null(&dt_pi_at_s), make_not_null(&dt_phi_at_s),
            get<::Tags::AtIndex<::Tags::MirrorView<spacetime_metric_tag>>>(
                vars_at_s),
            get<::Tags::AtIndex<::Tags::MirrorView<pi_tag>>>(vars_at_s),
            get<::Tags::AtIndex<::Tags::MirrorView<phi_tag>>>(vars_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<spacetime_metric_tag>,
                              tmpl::size_t<volume_dim>, Frame::Inertial>>>(
                derivs_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<pi_tag>,
                              tmpl::size_t<volume_dim>, Frame::Inertial>>>(
                derivs_at_s),
            get<::Tags::AtIndex<
                ::Tags::deriv<::Tags::MirrorView<phi_tag>,
                              tmpl::size_t<volume_dim>, Frame::Inertial>>>(
                derivs_at_s),
            get<::Tags::AtIndex<device_gauge_h_tag>>(gauge_data_at_s),
            get<::Tags::AtIndex<device_spacetime_deriv_gauge_h_tag>>(
                gauge_data_at_s),
            get(gamma0_at_s), get(gamma1_at_s), get(gamma2_at_s));

        for (size_t mu = 0; mu < volume_dim + 1; ++mu) {
          for (size_t nu = mu; nu < volume_dim + 1; ++nu) {
            dt_spacetime_metric.get(mu, nu)[static_cast<size_t>(s)] =
                dt_spacetime_metric_at_s.get(mu, nu);
            dt_pi.get(mu, nu)[static_cast<size_t>(s)] = dt_pi_at_s.get(mu, nu);
            for (size_t i = 0; i < volume_dim; ++i) {
              dt_phi.get(i, mu, nu)[static_cast<size_t>(s)] =
                  dt_phi_at_s.get(i, mu, nu);
            }
          }
        }
      });

  // Boundary packaging is not wired yet for GH Kokkos. Keep boundary maps
  // empty so this action can be exercised in volume-only mode.
  outgoing_boundary_data->clear();
  external_boundary_data->clear();
}

}  // namespace gh::Actions
