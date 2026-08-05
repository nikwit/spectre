// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <pup.h>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/ParseError.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/SpecialRelativity/LorentzBoostMatrix.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace gh::Solutions {
namespace affine_map_model {
namespace {
// The six independent strain components in column order.
constexpr std::array<std::array<size_t, 2>, 6> sigma_pairs{
    {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 2}}}};
}  // namespace

void inverse_metric_combination(
    const gsl::not_null<tnsr::AA<DataVector, 3>*> out,
    const std::array<DataVector, 3>& y, const double mass,
    const double background_weight, const std::array<double, 13>& c) {
  const size_t n_points = y[0].size();
  DataVector rho(n_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    rho += square(gsl::at(y, i));
  }
  rho = sqrt(rho);
  std::array<DataVector, 3> n{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(n, i) = gsl::at(y, i) / rho;
  }

  const DataVector r = rho + mass;
  const DataVector f = 1. - 2. * mass / r;
  const DataVector big_f = 1. + 2. * mass / r;
  const DataVector h = 4. * square(mass) / square(r);
  const DataVector h_prime = -8. * square(mass) / cube(r);
  const DataVector g_perp = square(rho) / square(r);
  const DataVector g_perp_prime = 2. * mass * rho / cube(r);
  const DataVector g_par = -square(mass) / square(r);
  const DataVector g_par_prime = 2. * square(mass) / cube(r);
  const DataVector background_tt = f * square(big_f) - 2. * big_f;
  const DataVector background_tt_prime =
      2. * mass / square(r) * (square(big_f) - 2. * f * big_f + 2.);

  const double c_qdot0 = c[0];
  const std::array<double, 3> c_beta{{c[1], c[2], c[3]}};
  const std::array<double, 3> c_qdot{{c[4], c[5], c[6]}};
  std::array<std::array<double, 3>, 3> c_sigma{};
  for (size_t pair = 0; pair < 6; ++pair) {
    const auto [i, j] = gsl::at(sigma_pairs, pair);
    gsl::at(gsl::at(c_sigma, i), j) = c[7 + pair];
    gsl::at(gsl::at(c_sigma, j), i) = c[7 + pair];
  }

  DataVector beta_dot_n(n_points, 0.);
  DataVector sigma_nn(n_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    beta_dot_n += gsl::at(c_beta, i) * gsl::at(n, i);
    for (size_t j = 0; j < 3; ++j) {
      sigma_nn +=
          gsl::at(gsl::at(c_sigma, i), j) * gsl::at(n, i) * gsl::at(n, j);
    }
  }

  get<0, 0>(*out) = background_weight * background_tt +
                    2. * c_qdot0 * background_tt + 2. * h * beta_dot_n -
                    rho * sigma_nn * background_tt_prime;
  for (size_t i = 0; i < 3; ++i) {
    out->get(0, i + 1) = (background_weight + c_qdot0) * h * gsl::at(n, i) +
                         g_par * beta_dot_n * gsl::at(n, i) +
                         g_perp * gsl::at(c_beta, i) -
                         big_f * (1. + h) * gsl::at(c_qdot, i) +
                         sigma_nn * (h - rho * h_prime) * gsl::at(n, i);
    for (size_t j = i; j < 3; ++j) {
      out->get(i + 1, j + 1) =
          background_weight * (g_par * gsl::at(n, i) * gsl::at(n, j) +
                               (i == j ? 1. : 0.) * g_perp) +
          h * (gsl::at(c_qdot, i) * gsl::at(n, j) +
               gsl::at(c_qdot, j) * gsl::at(n, i)) +
          2. * g_perp * gsl::at(gsl::at(c_sigma, i), j) -
          (i == j ? 1. : 0.) * rho * sigma_nn * g_perp_prime +
          sigma_nn * (2. * g_par - rho * g_par_prime) * gsl::at(n, i) *
              gsl::at(n, j);
    }
  }
}

void spatial_derivative_of_inverse_metric_combination(
    const gsl::not_null<tnsr::iAA<DataVector, 3>*> out,
    const std::array<DataVector, 3>& y, const double mass,
    const double background_weight, const std::array<double, 13>& c) {
  const size_t n_points = y[0].size();
  DataVector rho(n_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    rho += square(gsl::at(y, i));
  }
  rho = sqrt(rho);
  std::array<DataVector, 3> n{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(n, i) = gsl::at(y, i) / rho;
  }

  const DataVector r = rho + mass;
  const DataVector f = 1. - 2. * mass / r;
  const DataVector big_f = 1. + 2. * mass / r;
  const DataVector f_prime = 2. * mass / square(r);
  const DataVector big_f_prime = -2. * mass / square(r);
  const DataVector h = 4. * square(mass) / square(r);
  const DataVector h_prime = -8. * square(mass) / cube(r);
  const DataVector h_pp = 24. * square(mass) / (square(r) * square(r));
  const DataVector g_perp = square(rho) / square(r);
  const DataVector g_perp_prime = 2. * mass * rho / cube(r);
  const DataVector g_perp_pp =
      2. * mass / cube(r) - 6. * mass * rho / (square(r) * square(r));
  const DataVector g_par = -square(mass) / square(r);
  const DataVector g_par_prime = 2. * square(mass) / cube(r);
  const DataVector g_par_pp = -6. * square(mass) / (square(r) * square(r));
  const DataVector background_tt_prime =
      2. * mass / square(r) * (square(big_f) - 2. * f * big_f + 2.);
  const DataVector background_tt_pp =
      -4. * mass / cube(r) * (square(big_f) - 2. * f * big_f + 2.) +
      8. * square(mass) / (square(r) * square(r)) * (f - 2. * big_f);

  const double c_qdot0 = c[0];
  const std::array<double, 3> c_beta{{c[1], c[2], c[3]}};
  const std::array<double, 3> c_qdot{{c[4], c[5], c[6]}};
  std::array<std::array<double, 3>, 3> c_sigma{};
  for (size_t pair = 0; pair < 6; ++pair) {
    const auto [i, j] = gsl::at(sigma_pairs, pair);
    gsl::at(gsl::at(c_sigma, i), j) = c[7 + pair];
    gsl::at(gsl::at(c_sigma, j), i) = c[7 + pair];
  }

  DataVector beta_dot_n(n_points, 0.);
  DataVector sigma_nn(n_points, 0.);
  std::array<DataVector, 3> sigma_n{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(sigma_n, i) = DataVector(n_points, 0.);
  }
  for (size_t i = 0; i < 3; ++i) {
    beta_dot_n += gsl::at(c_beta, i) * gsl::at(n, i);
    for (size_t j = 0; j < 3; ++j) {
      sigma_nn +=
          gsl::at(gsl::at(c_sigma, i), j) * gsl::at(n, i) * gsl::at(n, j);
      gsl::at(sigma_n, i) += gsl::at(gsl::at(c_sigma, i), j) * gsl::at(n, j);
    }
  }

  for (size_t k = 0; k < 3; ++k) {
    const DataVector& nk = gsl::at(n, k);
    // dk of the angular factors
    const DataVector dk_bn = (gsl::at(c_beta, k) - beta_dot_n * nk) / rho;
    const DataVector dk_snn = 2. * (gsl::at(sigma_n, k) - sigma_nn * nk) / rho;
    const auto dk_n = [&n, &rho, &nk, k](const size_t i) -> DataVector {
      return ((i == k ? 1. : 0.) - gsl::at(n, i) * nk) / rho;
    };

    out->get(k, 0, 0) =
        (background_weight + 2. * c_qdot0) * background_tt_prime * nk +
        2. * h_prime * nk * beta_dot_n + 2. * h * dk_bn -
        nk * sigma_nn * background_tt_prime -
        rho * dk_snn * background_tt_prime -
        rho * sigma_nn * background_tt_pp * nk;
    for (size_t i = 0; i < 3; ++i) {
      const DataVector dk_ni = dk_n(i);
      out->get(k, 0, i + 1) =
          (background_weight + c_qdot0) *
              (h_prime * nk * gsl::at(n, i) + h * dk_ni) +
          g_par_prime * nk * beta_dot_n * gsl::at(n, i) +
          g_par * (dk_bn * gsl::at(n, i) + beta_dot_n * dk_ni) +
          g_perp_prime * nk * gsl::at(c_beta, i) -
          (big_f_prime * (1. + h) + big_f * h_prime) * nk * gsl::at(c_qdot, i) +
          (dk_snn * (h - rho * h_prime) - sigma_nn * rho * h_pp * nk) *
              gsl::at(n, i) +
          sigma_nn * (h - rho * h_prime) * dk_ni;
      for (size_t j = i; j < 3; ++j) {
        const DataVector dk_nj = dk_n(j);
        const DataVector dk_ninj =
            dk_ni * gsl::at(n, j) + gsl::at(n, i) * dk_nj;
        out->get(k, i + 1, j + 1) =
            background_weight *
                (g_par_prime * nk * gsl::at(n, i) * gsl::at(n, j) +
                 g_par * dk_ninj + (i == j ? 1. : 0.) * g_perp_prime * nk) +
            h_prime * nk *
                (gsl::at(c_qdot, i) * gsl::at(n, j) +
                 gsl::at(c_qdot, j) * gsl::at(n, i)) +
            h * (gsl::at(c_qdot, i) * dk_nj + gsl::at(c_qdot, j) * dk_ni) +
            2. * g_perp_prime * nk * gsl::at(gsl::at(c_sigma, i), j) -
            (i == j ? 1. : 0.) *
                (nk * sigma_nn * g_perp_prime + rho * dk_snn * g_perp_prime +
                 rho * sigma_nn * g_perp_pp * nk) +
            (dk_snn * (2. * g_par - rho * g_par_prime) +
             sigma_nn * (g_par_prime - rho * g_par_pp) * nk) *
                gsl::at(n, i) * gsl::at(n, j) +
            sigma_nn * (2. * g_par - rho * g_par_prime) * dk_ninj;
      }
    }
  }
}

void first_order_evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, const double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const bool centre_advection) {
  const size_t n_points = get<0>(x).size();
  std::array<DataVector, 3> y{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(y, i) = x.get(i) - gsl::at(center, i);
  }

  set_number_of_grid_points(spacetime_metric, n_points);
  set_number_of_grid_points(pi, n_points);
  set_number_of_grid_points(phi, n_points);

  // G^{ab} = G_0^{ab} + delta G^{ab} is the object derived directly from
  // the first-order coordinate transformation.  Invert only G_0 and form
  // delta g_ab = -(g_0 delta G g_0)_ab instead of exactly inverting the sum,
  // which would silently retain all powers of the first-order coefficients.
  const std::array<double, 13> zero_parameters{};
  tnsr::AA<DataVector, 3> background_inverse(n_points);
  tnsr::AA<DataVector, 3> delta_inverse(n_points);
  inverse_metric_combination(make_not_null(&background_inverse), y, mass, 1.,
                             zero_parameters);
  inverse_metric_combination(make_not_null(&delta_inverse), y, mass, 0., p);
  const tnsr::aa<DataVector, 3> background_metric =
      determinant_and_inverse(background_inverse).second;
  tnsr::aa<DataVector, 3> delta_metric(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          delta_metric.get(a, b) -= background_metric.get(a, c) *
                                    delta_inverse.get(c, d) *
                                    background_metric.get(d, b);
        }
      }
      spacetime_metric->get(a, b) =
          background_metric.get(a, b) + delta_metric.get(a, b);
    }
  }

  // Spatially differentiate the same truncated inverse relation.  The four
  // terms below are precisely the background derivative and the three terms
  // linear in (delta g, delta G); no product of two perturbations is kept.
  tnsr::iAA<DataVector, 3> deriv_background_inverse(n_points);
  tnsr::iAA<DataVector, 3> deriv_delta_inverse(n_points);
  spatial_derivative_of_inverse_metric_combination(
      make_not_null(&deriv_background_inverse), y, mass, 1., zero_parameters);
  spatial_derivative_of_inverse_metric_combination(
      make_not_null(&deriv_delta_inverse), y, mass, 0., p);
  tnsr::iaa<DataVector, 3> background_phi(n_points, 0.);
  tnsr::iaa<DataVector, 3> delta_phi(n_points, 0.);
  for (size_t k = 0; k < 3; ++k) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t c = 0; c < 4; ++c) {
          for (size_t d = 0; d < 4; ++d) {
            background_phi.get(k, a, b) -=
                background_metric.get(a, c) *
                deriv_background_inverse.get(k, c, d) *
                background_metric.get(d, b);
            delta_phi.get(k, a, b) -=
                delta_metric.get(a, c) * deriv_background_inverse.get(k, c, d) *
                    background_metric.get(d, b) +
                background_metric.get(a, c) * deriv_delta_inverse.get(k, c, d) *
                    background_metric.get(d, b) +
                background_metric.get(a, c) *
                    deriv_background_inverse.get(k, c, d) *
                    delta_metric.get(d, b);
          }
        }
        phi->get(k, a, b) =
            background_phi.get(k, a, b) + delta_phi.get(k, a, b);
      }
    }
  }

  // Expand the model lapse and shift directly from the already-linear
  // inverse metric.  Computing them from g_0 + delta g with the full 3+1
  // algebra would reintroduce quadratic and higher powers of p.
  const DataVector background_lapse = 1. / sqrt(-get<0, 0>(background_inverse));
  const DataVector delta_lapse = 0.5 * background_lapse * background_lapse *
                                 background_lapse * get<0, 0>(delta_inverse);
  std::array<DataVector, 3> background_shift{};
  std::array<DataVector, 3> delta_shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(background_shift, i) =
        -background_inverse.get(0, i + 1) / get<0, 0>(background_inverse);
    gsl::at(delta_shift, i) =
        -delta_inverse.get(0, i + 1) / get<0, 0>(background_inverse) +
        background_inverse.get(0, i + 1) * get<0, 0>(delta_inverse) /
            (get<0, 0>(background_inverse) * get<0, 0>(background_inverse));
  }

  // Slow-time ordering: D_t p_A contributes only at O(epsilon^2) and is
  // absent.  The center q^i is a zeroth-order placement, so its physical
  // velocity is O(epsilon) and advects the background metric at first order.
  // Use Phi^(0), since qdot^k Phi^(1)_k is O(epsilon^2).
  tnsr::aa<DataVector, 3> dt_metric_first_order(n_points, 0.);
  if (centre_advection) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t k = 0; k < 3; ++k) {
          dt_metric_first_order.get(a, b) -=
              gsl::at(p, 4 + k) * background_phi.get(k, a, b);
        }
      }
    }
  }

  // Pi = (shift^k Phi_k - partial_t g) / lapse, expanded once about the
  // stationary background.
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      DataVector background_pi(n_points, 0.);
      DataVector delta_pi = -dt_metric_first_order.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        background_pi +=
            gsl::at(background_shift, k) * background_phi.get(k, a, b);
        delta_pi += gsl::at(delta_shift, k) * background_phi.get(k, a, b) +
                    gsl::at(background_shift, k) * delta_phi.get(k, a, b);
      }
      background_pi /= background_lapse;
      delta_pi = delta_pi / background_lapse -
                 background_pi * delta_lapse / background_lapse;
      pi->get(a, b) = background_pi + delta_pi;
    }
  }
}

void first_order_boosted_evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, const double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const std::array<double, 3>& boost_velocity, const bool centre_advection) {
  if (boost_velocity == std::array<double, 3>{{0., 0., 0.}}) {
    first_order_evolved_variables(spacetime_metric, pi, phi, x, mass, center, p,
                                  centre_advection);
    return;
  }
  const size_t n_points = get<0>(x).size();

  tnsr::I<double, 3, Frame::NoFrame> minus_v{};
  for (size_t i = 0; i < 3; ++i) {
    minus_v.get(i) = -gsl::at(boost_velocity, i);
  }
  const auto boost = sr::lorentz_boost_matrix(minus_v);

  // The center is its instantaneous lab position.  Stationarity of the
  // rest-frame fields lets the time-dependent translation cancel explicitly:
  // X^ibar = M^ibar_j (x^j - center^j(t)).
  tnsr::I<DataVector, 3> rest_coords(n_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      rest_coords.get(i) +=
          boost.get(i + 1, j + 1) * (x.get(j) - gsl::at(center, j));
    }
  }

  const std::array<double, 13> zero{};
  tnsr::aa<DataVector, 3> rest_background_metric{};
  tnsr::aa<DataVector, 3> rest_background_pi{};
  tnsr::iaa<DataVector, 3> rest_background_phi{};
  first_order_evolved_variables(make_not_null(&rest_background_metric),
                                make_not_null(&rest_background_pi),
                                make_not_null(&rest_background_phi),
                                rest_coords, mass, {{0., 0., 0.}}, zero, false);
  tnsr::aa<DataVector, 3> rest_full_metric{};
  tnsr::aa<DataVector, 3> rest_full_pi{};
  tnsr::iaa<DataVector, 3> rest_full_phi{};
  first_order_evolved_variables(make_not_null(&rest_full_metric),
                                make_not_null(&rest_full_pi),
                                make_not_null(&rest_full_phi), rest_coords,
                                mass, {{0., 0., 0.}}, p, centre_advection);

  tnsr::aa<DataVector, 3> rest_delta_metric = rest_full_metric;
  tnsr::iaa<DataVector, 3> rest_delta_phi = rest_full_phi;
  for (size_t storage = 0; storage < rest_delta_metric.size(); ++storage) {
    rest_delta_metric[storage] -= rest_background_metric[storage];
  }
  for (size_t storage = 0; storage < rest_delta_phi.size(); ++storage) {
    rest_delta_phi[storage] -= rest_background_phi[storage];
  }

  // In the strict slow-time model the rest-frame background is stationary
  // and the only retained time derivative is center advection of g^(0).
  tnsr::aa<DataVector, 3> rest_delta_dt_metric(n_points, 0.);
  if (centre_advection) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t k = 0; k < 3; ++k) {
          rest_delta_dt_metric.get(a, b) -=
              gsl::at(p, 4 + k) * rest_background_phi.get(k, a, b);
        }
      }
    }
  }

  set_number_of_grid_points(spacetime_metric, n_points);
  set_number_of_grid_points(phi, n_points);
  set_number_of_grid_points(pi, n_points);
  tnsr::aa<DataVector, 3> lab_background_metric(n_points, 0.);
  tnsr::aa<DataVector, 3> lab_delta_metric(n_points, 0.);
  std::array<std::array<std::array<DataVector, 4>, 4>, 4>
      lab_background_deriv{};
  std::array<std::array<std::array<DataVector, 4>, 4>, 4> lab_delta_deriv{};

  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          lab_background_metric.get(a, b) += boost.get(c, a) * boost.get(d, b) *
                                             rest_background_metric.get(c, d);
          lab_delta_metric.get(a, b) +=
              boost.get(c, a) * boost.get(d, b) * rest_delta_metric.get(c, d);
        }
      }
      spacetime_metric->get(a, b) =
          lab_background_metric.get(a, b) + lab_delta_metric.get(a, b);
    }
  }

  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        auto& background_entry =
            gsl::at(gsl::at(gsl::at(lab_background_deriv, mu), a), b);
        auto& delta_entry =
            gsl::at(gsl::at(gsl::at(lab_delta_deriv, mu), a), b);
        background_entry = DataVector(n_points, 0.);
        delta_entry = DataVector(n_points, 0.);
        for (size_t e = 0; e < 4; ++e) {
          for (size_t c = 0; c < 4; ++c) {
            for (size_t d = 0; d < 4; ++d) {
              const double weight =
                  boost.get(e, mu) * boost.get(c, a) * boost.get(d, b);
              if (e > 0) {
                background_entry +=
                    weight * rest_background_phi.get(e - 1, c, d);
                delta_entry += weight * rest_delta_phi.get(e - 1, c, d);
              } else {
                delta_entry += weight * rest_delta_dt_metric.get(c, d);
              }
            }
          }
        }
      }
    }
  }
  for (size_t k = 0; k < 3; ++k) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        phi->get(k, a, b) =
            gsl::at(gsl::at(gsl::at(lab_background_deriv, k + 1), a), b) +
            gsl::at(gsl::at(gsl::at(lab_delta_deriv, k + 1), a), b);
      }
    }
  }

  // Expand the lab lapse, shift, and Pi once about the finite-boost
  // background.  This is the step that prevents the 3+1 reconstruction from
  // reintroducing products of first-order coefficients.
  const auto background_inverse =
      determinant_and_inverse(lab_background_metric).second;
  tnsr::AA<DataVector, 3> delta_inverse(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          delta_inverse.get(a, b) -= background_inverse.get(a, c) *
                                     lab_delta_metric.get(c, d) *
                                     background_inverse.get(d, b);
        }
      }
    }
  }
  const DataVector background_lapse = 1. / sqrt(-get<0, 0>(background_inverse));
  const DataVector delta_lapse = 0.5 * background_lapse * background_lapse *
                                 background_lapse * get<0, 0>(delta_inverse);
  std::array<DataVector, 3> background_shift{};
  std::array<DataVector, 3> delta_shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(background_shift, i) =
        -background_inverse.get(0, i + 1) / get<0, 0>(background_inverse);
    gsl::at(delta_shift, i) =
        -delta_inverse.get(0, i + 1) / get<0, 0>(background_inverse) +
        background_inverse.get(0, i + 1) * get<0, 0>(delta_inverse) /
            square(get<0, 0>(background_inverse));
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      const DataVector& background_dt =
          gsl::at(gsl::at(gsl::at(lab_background_deriv, 0), a), b);
      const DataVector& delta_dt =
          gsl::at(gsl::at(gsl::at(lab_delta_deriv, 0), a), b);
      DataVector background_pi = -background_dt;
      DataVector delta_pi = -delta_dt;
      for (size_t k = 0; k < 3; ++k) {
        const DataVector& background_phi =
            gsl::at(gsl::at(gsl::at(lab_background_deriv, k + 1), a), b);
        const DataVector& delta_phi =
            gsl::at(gsl::at(gsl::at(lab_delta_deriv, k + 1), a), b);
        background_pi += gsl::at(background_shift, k) * background_phi;
        delta_pi += gsl::at(delta_shift, k) * background_phi +
                    gsl::at(background_shift, k) * delta_phi;
      }
      background_pi /= background_lapse;
      delta_pi = delta_pi / background_lapse -
                 background_pi * delta_lapse / background_lapse;
      pi->get(a, b) = background_pi + delta_pi;
    }
  }
}

void evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, const double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const std::array<double, 13>& pdot, const bool centre_advection) {
  const size_t n_points = get<0>(x).size();
  std::array<DataVector, 3> y{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(y, i) = x.get(i) - gsl::at(center, i);
  }

  set_number_of_grid_points(spacetime_metric, n_points);
  set_number_of_grid_points(pi, n_points);
  set_number_of_grid_points(phi, n_points);

  tnsr::AA<DataVector, 3> inverse_metric(n_points);
  inverse_metric_combination(make_not_null(&inverse_metric), y, mass, 1., p);
  *spacetime_metric = determinant_and_inverse(inverse_metric).second;

  // Phi analytically: Phi_kab = d_k g_ab = -(g d_k G^{-1} g)_ab
  tnsr::iAA<DataVector, 3> dk_inverse(n_points);
  spatial_derivative_of_inverse_metric_combination(make_not_null(&dk_inverse),
                                                   y, mass, 1., p);
  for (size_t k = 0; k < 3; ++k) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        phi->get(k, a, b) = 0.;
        for (size_t cc = 0; cc < 4; ++cc) {
          for (size_t d = 0; d < 4; ++d) {
            phi->get(k, a, b) -= spacetime_metric->get(a, cc) *
                                 dk_inverse.get(k, cc, d) *
                                 spacetime_metric->get(d, b);
          }
        }
      }
    }
  }

  // d_t g_ab = -(g S g)_ab with S^{ab} = sum_A pdot_A R_A^{ab}: the
  // coefficient drift of the model at a fixed point
  tnsr::AA<DataVector, 3> rate_direction(n_points);
  inverse_metric_combination(make_not_null(&rate_direction), y, mass, 0., pdot);
  tnsr::aa<DataVector, 3> dt_metric(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          dt_metric.get(a, b) -= spacetime_metric->get(a, c) *
                                 rate_direction.get(c, d) *
                                 spacetime_metric->get(d, b);
        }
      }
    }
  }

  // Centre advection. The map displaces the centre by qdot^i t, a secular
  // piece the instantaneous model drops (a fit at fixed t absorbs the
  // accumulated displacement into the centre offset instead). Its time
  // derivative survives: g(x, t) = g_hat(x - c - qdot t) gives
  // d_t g_ab = -qdot^k d_k g_ab = -qdot^k Phi_kab, evaluated with the Phi
  // just built analytically. This is the only place the velocity reaches
  // Pi; without it qdot^i is determined by g and Phi alone, which is the
  // weak channel (findings 9: freely fitted time jets inflate ~1e3 off the
  // metric while the same quantities come out clean off Pi). It vanishes
  // identically whenever qdot^i = 0, so it cannot perturb any run whose
  // velocity is pinned to a zero CenterVelocity.
  if (centre_advection) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t k = 0; k < 3; ++k) {
          dt_metric.get(a, b) -= gsl::at(p, 4 + k) * phi->get(k, a, b);
        }
      }
    }
  }

  // lapse and shift of the model metric itself;
  // Pi_ab = (beta^k Phi_kab - d_t g_ab) / alpha
  const DataVector lapse = 1. / sqrt(-get<0, 0>(inverse_metric));
  std::array<DataVector, 3> shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(shift, i) =
        -inverse_metric.get(0, i + 1) / get<0, 0>(inverse_metric);
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      pi->get(a, b) = -dt_metric.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        pi->get(a, b) += gsl::at(shift, k) * phi->get(k, a, b);
      }
      pi->get(a, b) /= lapse;
    }
  }
}
void boosted_evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    const gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, const double time, const double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const std::array<double, 13>& pdot,
    const std::array<double, 3>& boost_velocity, const bool centre_advection) {
  if (boost_velocity == std::array<double, 3>{{0., 0., 0.}}) {
    evolved_variables(spacetime_metric, pi, phi, x, mass, center, p, pdot,
                      centre_advection);
    return;
  }
  const size_t n_points = get<0>(x).size();

  // lorentz_boost_matrix(v) is Lambda^a_abar = dx^a_lab / dx^abar_rest, so the
  // matrix taking lab coordinates and lab derivative indices to the rest frame
  // is the one built from -v. Call it M^abar_a = dx^abar / dx^a.
  tnsr::I<double, 3, Frame::NoFrame> minus_v{};
  for (size_t i = 0; i < 3; ++i) {
    minus_v.get(i) = -gsl::at(boost_velocity, i);
  }
  const auto boost = sr::lorentz_boost_matrix(minus_v);

  // rest-frame coordinates of the lab points, measured from the centre:
  // X^abar = M^abar_b (t, x - centre)^b. The hole is at the origin of the
  // rest frame, so pass these to the static solution with a zero centre.
  tnsr::I<DataVector, 3> rest_coords(n_points);
  for (size_t i = 0; i < 3; ++i) {
    rest_coords.get(i) = boost.get(i + 1, 0) * time;
    for (size_t j = 0; j < 3; ++j) {
      rest_coords.get(i) +=
          boost.get(i + 1, j + 1) * (x.get(j) - gsl::at(center, j));
    }
  }

  tnsr::aa<DataVector, 3> rest_metric{};
  tnsr::aa<DataVector, 3> rest_pi{};
  tnsr::iaa<DataVector, 3> rest_phi{};
  // The map parameters are evaluated at the *lab* time. For a pure boost
  // (p = pdot = 0) that is exact; combining a boost with a time-dependent p
  // mixes the two frames' time coordinates and is only first-order
  // consistent, which is why the boost tests run at p = 0.
  evolved_variables(make_not_null(&rest_metric), make_not_null(&rest_pi),
                    make_not_null(&rest_phi), rest_coords, mass, {{0., 0., 0.}},
                    p, pdot, centre_advection);

  // d_nu g_ab in the rest frame: the spatial part is Phi, and the time part
  // comes back out of Pi's definition, dt g = beta^k Phi_k - alpha Pi.
  const auto rest_inverse = determinant_and_inverse(rest_metric).second;
  const DataVector rest_lapse = 1. / sqrt(-get<0, 0>(rest_inverse));
  std::array<std::array<std::array<DataVector, 4>, 4>, 4> rest_deriv{};
  for (size_t nu = 0; nu < 4; ++nu) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        DataVector& entry = gsl::at(gsl::at(gsl::at(rest_deriv, nu), a), b);
        if (nu == 0) {
          entry = -rest_lapse * rest_pi.get(a, b);
          for (size_t k = 0; k < 3; ++k) {
            const DataVector shift_k =
                -rest_inverse.get(0, k + 1) / get<0, 0>(rest_inverse);
            entry += shift_k * rest_phi.get(k, a, b);
          }
        } else {
          entry = rest_phi.get(nu - 1, a, b);
        }
      }
    }
  }

  // g_ab(lab) = M^c_a M^d_b g_cd(rest);
  // d_mu g_ab(lab) = M^e_mu M^c_a M^d_b d_e g_cd(rest).
  set_number_of_grid_points(spacetime_metric, n_points);
  set_number_of_grid_points(phi, n_points);
  set_number_of_grid_points(pi, n_points);
  std::array<std::array<std::array<DataVector, 4>, 4>, 4> lab_deriv{};
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      spacetime_metric->get(a, b) = 0.;
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          spacetime_metric->get(a, b) +=
              boost.get(c, a) * boost.get(d, b) * rest_metric.get(c, d);
        }
      }
    }
  }
  for (size_t mu = 0; mu < 4; ++mu) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        DataVector& entry = gsl::at(gsl::at(gsl::at(lab_deriv, mu), a), b);
        entry = DataVector(n_points, 0.);
        for (size_t e = 0; e < 4; ++e) {
          for (size_t c = 0; c < 4; ++c) {
            for (size_t d = 0; d < 4; ++d) {
              entry += boost.get(e, mu) * boost.get(c, a) * boost.get(d, b) *
                       gsl::at(gsl::at(gsl::at(rest_deriv, e), c), d);
            }
          }
        }
      }
    }
  }
  for (size_t k = 0; k < 3; ++k) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        phi->get(k, a, b) = gsl::at(gsl::at(gsl::at(lab_deriv, k + 1), a), b);
      }
    }
  }
  // Pi from the lab lapse and shift
  const auto lab_inverse = determinant_and_inverse(*spacetime_metric).second;
  const DataVector lab_lapse = 1. / sqrt(-get<0, 0>(lab_inverse));
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      pi->get(a, b) = -gsl::at(gsl::at(gsl::at(lab_deriv, 0), a), b);
      for (size_t k = 0; k < 3; ++k) {
        const DataVector shift_k =
            -lab_inverse.get(0, k + 1) / get<0, 0>(lab_inverse);
        pi->get(a, b) += shift_k * phi->get(k, a, b);
      }
      pi->get(a, b) /= lab_lapse;
    }
  }
}
}  // namespace affine_map_model

AffineMappedHarmonicSchwarzschild::AffineMappedHarmonicSchwarzschild(
    const double mass, const std::array<double, volume_dim>& center,
    std::vector<double> parameter_times,
    std::vector<std::array<double, number_of_parameters>> parameter_values,
    std::vector<std::array<double, number_of_parameters>> parameter_rates,
    const std::array<double, volume_dim>& velocity,
    const Options::Context& context)
    : mass_(mass),
      center_(center),
      parameter_times_(std::move(parameter_times)),
      parameter_values_(std::move(parameter_values)),
      parameter_rates_(std::move(parameter_rates)),
      velocity_(velocity) {
  if (parameter_times_.empty()) {
    PARSE_ERROR(context, "ParameterTimes must not be empty.");
  }
  double v_squared = 0.;
  for (size_t i = 0; i < volume_dim; ++i) {
    v_squared += square(gsl::at(velocity_, i));
  }
  if (v_squared >= 1.) {
    PARSE_ERROR(context, "Velocity must be subluminal, but |v|^2 = "
                             << v_squared << ".");
  }
  if (not std::is_sorted(parameter_times_.begin(), parameter_times_.end()) or
      std::adjacent_find(parameter_times_.begin(), parameter_times_.end()) !=
          parameter_times_.end()) {
    PARSE_ERROR(context, "ParameterTimes must be strictly increasing.");
  }
  if (parameter_values_.size() != parameter_times_.size() or
      parameter_rates_.size() != parameter_times_.size()) {
    PARSE_ERROR(context,
                "ParameterValues and ParameterRates must have one row per "
                "entry of ParameterTimes, but got "
                    << parameter_values_.size() << " and "
                    << parameter_rates_.size() << " rows for "
                    << parameter_times_.size() << " times.");
  }
  build_splines();
}

AffineMappedHarmonicSchwarzschild::AffineMappedHarmonicSchwarzschild(
    CkMigrateMessage* const msg)
    : InitialData(msg) {}

std::unique_ptr<evolution::initial_data::InitialData>
AffineMappedHarmonicSchwarzschild::get_clone() const {
  return std::make_unique<AffineMappedHarmonicSchwarzschild>(*this);
}

void AffineMappedHarmonicSchwarzschild::build_splines() {
  value_splines_.clear();
  rate_splines_.clear();
  if (parameter_times_.size() < 3) {
    return;
  }
  for (size_t a = 0; a < number_of_parameters; ++a) {
    std::vector<double> values(parameter_times_.size());
    std::vector<double> rates(parameter_times_.size());
    for (size_t k = 0; k < parameter_times_.size(); ++k) {
      values[k] = gsl::at(parameter_values_[k], a);
      rates[k] = gsl::at(parameter_rates_[k], a);
    }
    value_splines_.emplace_back(parameter_times_, std::move(values));
    rate_splines_.emplace_back(parameter_times_, std::move(rates));
  }
}

void AffineMappedHarmonicSchwarzschild::map_parameters(
    const gsl::not_null<std::array<double, number_of_parameters>*> values,
    const gsl::not_null<std::array<double, number_of_parameters>*> rates,
    const double time) const {
  if (parameter_times_.size() == 1) {
    *values = parameter_values_[0];
    *rates = parameter_rates_[0];
    return;
  }
  const double clamped =
      std::clamp(time, parameter_times_.front(), parameter_times_.back());
  if (parameter_times_.size() == 2) {
    const double weight = (clamped - parameter_times_[0]) /
                          (parameter_times_[1] - parameter_times_[0]);
    for (size_t a = 0; a < number_of_parameters; ++a) {
      gsl::at(*values, a) = (1. - weight) * gsl::at(parameter_values_[0], a) +
                            weight * gsl::at(parameter_values_[1], a);
      gsl::at(*rates, a) = (1. - weight) * gsl::at(parameter_rates_[0], a) +
                           weight * gsl::at(parameter_rates_[1], a);
    }
    return;
  }
  for (size_t a = 0; a < number_of_parameters; ++a) {
    gsl::at(*values, a) = value_splines_[a](clamped);
    gsl::at(*rates, a) = rate_splines_[a](clamped);
  }
}

AffineMappedHarmonicSchwarzschild::AllVars
AffineMappedHarmonicSchwarzschild::all_variables(
    const tnsr::I<DataVector, volume_dim>& x, const double time) const {
  std::array<double, number_of_parameters> p{};
  std::array<double, number_of_parameters> pdot{};
  map_parameters(make_not_null(&p), make_not_null(&pdot), time);

  const size_t n_points = get<0>(x).size();
  AllVars result{};
  auto& spacetime_metric =
      get<gr::Tags::SpacetimeMetric<DataVector, volume_dim>>(result);
  auto& pi = get<gh::Tags::Pi<DataVector, volume_dim>>(result);
  auto& phi = get<gh::Tags::Phi<DataVector, volume_dim>>(result);
  affine_map_model::boosted_evolved_variables(
      make_not_null(&spacetime_metric), make_not_null(&pi), make_not_null(&phi),
      x, time, mass_, center_, p, pdot, velocity_);

  // lapse and shift from the metric
  auto& lapse = get<gr::Tags::Lapse<DataVector>>(result);
  auto& shift = get<gr::Tags::Shift<DataVector, volume_dim>>(result);
  const auto inverse_metric = determinant_and_inverse(spacetime_metric).second;
  get(lapse) = 1. / sqrt(-get<0, 0>(inverse_metric));
  for (size_t i = 0; i < 3; ++i) {
    shift.get(i) = -inverse_metric.get(0, i + 1) / get<0, 0>(inverse_metric);
  }

  // ADM quantities requested by gh::Actions::SetInitialData
  auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataVector, volume_dim>>(result);
  set_number_of_grid_points(make_not_null(&spatial_metric), n_points);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      spatial_metric.get(i, j) = spacetime_metric.get(i + 1, j + 1);
    }
  }
  tnsr::A<DataVector, 3> spacetime_unit_normal(n_points);
  get<0>(spacetime_unit_normal) = 1. / get(lapse);
  for (size_t i = 0; i < 3; ++i) {
    spacetime_unit_normal.get(i + 1) = -shift.get(i) / get(lapse);
  }
  auto& extrinsic_curvature =
      get<gr::Tags::ExtrinsicCurvature<DataVector, volume_dim>>(result);
  set_number_of_grid_points(make_not_null(&extrinsic_curvature), n_points);
  gh::extrinsic_curvature(make_not_null(&extrinsic_curvature),
                          spacetime_unit_normal, pi, phi);

  return result;
}

void AffineMappedHarmonicSchwarzschild::pup(PUP::er& p) {
  InitialData::pup(p);
  p | mass_;
  p | center_;
  p | parameter_times_;
  p | parameter_values_;
  p | parameter_rates_;
  p | velocity_;
  if (p.isUnpacking()) {
    build_splines();
  }
}

bool operator==(const AffineMappedHarmonicSchwarzschild& lhs,
                const AffineMappedHarmonicSchwarzschild& rhs) {
  return lhs.mass_ == rhs.mass_ and lhs.center_ == rhs.center_ and
         lhs.parameter_times_ == rhs.parameter_times_ and
         lhs.parameter_values_ == rhs.parameter_values_ and
         lhs.parameter_rates_ == rhs.parameter_rates_ and
         lhs.velocity_ == rhs.velocity_;
}

bool operator!=(const AffineMappedHarmonicSchwarzschild& lhs,
                const AffineMappedHarmonicSchwarzschild& rhs) {
  return not(lhs == rhs);
}

// NOLINTNEXTLINE
PUP::able::PUP_ID AffineMappedHarmonicSchwarzschild::my_PUP_ID = 0;
}  // namespace gh::Solutions
