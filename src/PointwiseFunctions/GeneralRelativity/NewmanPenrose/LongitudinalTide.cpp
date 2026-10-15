// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/LongitudinalTide.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <limits>
#include <vector>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
using C = std::complex<double>;
using C4 = std::array<C, 4>;
bool finite(const C z) {
  return std::isfinite(z.real()) and std::isfinite(z.imag());
}
FrameRegistration register_seed(const LongitudinalSeed& seed,
                                const RealMatrix& rotation, const double mass) {
  const size_t n = get(seed.coulomb).size();
  FrameRegistration r{};
  r.coulomb = seed.coulomb;
  r.rotation.a_bar = Scalar<ComplexDataVector>(n);
  r.rotation.b = Scalar<ComplexDataVector>(n);
  r.rotation.x = Scalar<ComplexDataVector>(n);
  r.measured_radius = Scalar<DataVector>(n);
  for (size_t p = 0; p < n; ++p) {
    const C c = get(seed.coulomb)[p];
    const C delta = seed.psi.get(2)[p] / c - 1.;
    // Rationalized root near zero. Compute a and b independently: x/b is
    // singular for an exactly one-sided null rotation (Psi1=0, Psi3!=0).
    const C x = delta / (3. + std::sqrt(9. + 6. * delta));
    const C b = seed.psi.get(1)[p] / (3. * c * (1. + 2. * x));
    const C a = seed.psi.get(3)[p] / (3. * c * (1. + x) * (1. + 2. * x));
    if (not(finite(a) and finite(b) and finite(x) and std::abs(x) < .5))
      ERROR("Longitudinal seed has an invalid type-D branch at " << p);
    get(r.rotation.a_bar)[p] = a;
    get(r.rotation.b)[p] = b;
    get(r.rotation.x)[p] = x;
    get(r.measured_radius)[p] = std::real(std::pow(-mass / c, 1. / 3.));
  }
  r.rotation.predicted_psi =
      type_d_scalars(r.coulomb, r.rotation.a_bar, r.rotation.b);
  r.pulled_back = pull_back(seed.psi, r.rotation.a_bar, r.rotation.b);
  r.member = tangent_boost_member(r.rotation.a_bar, r.rotation.b);
  r.radial_direction = TriadVector(n, 0.);
  r.transverse_velocity = TriadVector(n, 0.);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j) {
      r.radial_direction.get(i) +=
          rotation.get(j, i) * r.member.radial_direction.get(j);
      r.transverse_velocity.get(i) +=
          rotation.get(j, i) * r.member.transverse_velocity.get(j);
    }
  return r;
}
// Registered complex dyad, in the original adapted spacetime tetrad.
C4 rest_dyad(const FrameRegistration& r, const size_t p) {
  const double z = 1. / sqrt(2.);
  C4 l{{z, z, 0., 0.}}, k{{z, -z, 0., 0.}}, m{{0., 0., z, C{0., z}}};
  const C a = -std::conj(get(r.rotation.a_bar)[p]);
  const C b = -get(r.rotation.b)[p];
  C4 nk{}, nm{};
  for (size_t i = 0; i < 4; ++i) {
    nk[i] =
        k[i] + std::conj(a) * m[i] + a * std::conj(m[i]) + std::norm(a) * l[i];
    nm[i] = m[i] + a * l[i];
  }
  for (size_t i = 0; i < 4; ++i)
    nm[i] += b * nk[i];
  return nm;
}
}  // namespace

LongitudinalSeed longitudinal_seed(const WeylScalars& measured) {
  const size_t n = measured.get(1).size();
  LongitudinalSeed result{WeylScalars(n, C{0., 0.}),
                          Scalar<ComplexDataVector>(n),
                          std::numeric_limits<double>::infinity()};
  for (size_t p = 0; p < n; ++p) {
    const C p1 = measured.get(1)[p], p2 = measured.get(2)[p],
            p3 = measured.get(3)[p];
    if (not(finite(p1) and finite(p2) and finite(p3)))
      ERROR("Nonfinite longitudinal scalar at " << p);
    const C root = std::sqrt(9. * p2 * p2 - 8. * p1 * p3);
    const C plus = (-p2 + root) / 2., minus = (-p2 - root) / 2.;
    const C c = std::abs(plus - p2) <= std::abs(minus - p2) ? plus : minus;
    const C d = c + 2. * p2;
    const double scale = std::max({std::abs(p1), std::abs(p2), std::abs(p3)});
    if (scale == 0.) {
      ERROR("Singular longitudinal seed: all retained scalars vanish at " << p);
    }
    const double margin = std::min(std::abs(d), std::abs(2. * c + p2)) / scale;
    if (not(finite(c) and std::isfinite(margin) and margin > 1.e-10))
      ERROR("Singular longitudinal seed branch at " << p);
    result.minimum_relative_denominator =
        std::min(result.minimum_relative_denominator, margin);
    get(result.coulomb)[p] = c;
    result.psi.get(0)[p] = 2. * p1 * p1 / d;
    result.psi.get(1)[p] = p1;
    result.psi.get(2)[p] = p2;
    result.psi.get(3)[p] = p3;
    result.psi.get(4)[p] = 2. * p3 * p3 / d;
  }
  return result;
}

WeylScalars longitudinal_seed_response(const LongitudinalSeed& seed,
                                       const WeylScalars& column) {
  auto result = column;
  for (size_t p = 0; p < get(seed.coulomb).size(); ++p) {
    const C c = get(seed.coulomb)[p], p1 = seed.psi.get(1)[p],
            p2 = seed.psi.get(2)[p], p3 = seed.psi.get(3)[p], d = c + 2. * p2;
    const C dc = ((4. * p2 - c) * column.get(2)[p] -
                  2. * p3 * column.get(1)[p] - 2. * p1 * column.get(3)[p]) /
                 (2. * c + p2);
    result.get(0)[p] = 4. * p1 * column.get(1)[p] / d -
                       2. * p1 * p1 * (dc + 2. * column.get(2)[p]) / (d * d);
    result.get(4)[p] = 4. * p3 * column.get(3)[p] / d -
                       2. * p3 * p3 * (dc + 2. * column.get(2)[p]) / (d * d);
  }
  return result;
}

LongitudinalTideEvaluation evaluate_longitudinal_tide(
    const WeylScalars& measured, const RealMatrix& rotation,
    const tnsr::ii<DataVector, 3>& metric, const TriadVector& labels,
    const DataVector& weights, const double mass, const size_t lm,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& dk, const Scalar<DataVector>& dtk) {
  const size_t n = weights.size();
  if (n < 4 or measured.get(1).size() != n or measured.get(2).size() != n or
      measured.get(3).size() != n or not(std::isfinite(mass) and mass > 0.))
    ERROR("Invalid longitudinal tide data");
  const auto seed = longitudinal_seed(measured);
  LongitudinalTideEvaluation result{};
  result.minimum_relative_denominator = seed.minimum_relative_denominator;
  result.registration = register_seed(seed, rotation, mass);
  auto& r = result.registration;
  result.rapidity =
      invariant_rapidity(r.member, rotation, metric, lapse, shift, dk, dtk);
  result.frame = geometric_frame(r, result.rapidity, rotation, metric, labels,
                                 weights, mass, lm);
  const auto& frame = result.frame;
  const Scalar<DataVector> zero(n, 0.), boost{frame.boost};
  const auto rest_seed = type_iii(r.pulled_back, boost, zero);
  auto grad_rotation = RealMatrix(n, 0.), aux_rotation = RealMatrix(n, 0.);
  const auto inverse_lower = inverse_lower_triangular(cholesky_factor(metric));
  AdaptedFourVector ell{}, kay{};
  principal_null_pair(make_not_null(&ell), make_not_null(&kay),
                      r.rotation.a_bar, r.rotation.b);
  DataVector g_e(n), sigma(n);
  for (size_t p = 0; p < n; ++p) {
    std::array<double, 4> cov{{get(dtk)[p], 0., 0., 0.}};
    for (size_t i = 0; i < 3; ++i) {
      cov[0] -= shift.get(i)[p] * dk.get(i)[p];
      for (size_t j = 0; j < 3; ++j)
        for (size_t k = 0; k < 3; ++k)
          cov[i + 1] +=
              rotation.get(i, j)[p] * inverse_lower.get(j, k)[p] * dk.get(k)[p];
    }
    cov[0] /= get(lapse)[p];
    const auto m = rest_dyad(r, p);
    std::array<double, 3> g{};
    for (size_t j = 0; j < 4; ++j) {
      g[0] += (exp(frame.boost[p]) * ell.get(j)[p] -
               exp(-frame.boost[p]) * kay.get(j)[p]) /
              sqrt(2.) * cov[j];
      g[1] += sqrt(2.) * m[j].real() * cov[j];
      g[2] += sqrt(2.) * m[j].imag() * cov[j];
    }
    const double norm = std::hypot(g[0], g[1], g[2]);
    if (not(std::isfinite(norm) and norm > 0. and
            std::abs(g[0]) / norm > 1.e-10))
      ERROR("Longitudinal gradient direction is singular at " << p);
    const double orient = std::copysign(1., g[0]);
    for (double& v : g)
      v *= orient / norm;
    result.maximum_gradient_tilt = std::max(
        result.maximum_gradient_tilt, std::acos(std::clamp(g[0], -1., 1.)));
    // Shortest SO(3) rotation carrying e_radial to the oriented gradient.
    const std::array<double, 3> u{{0., -g[2], g[1]}};
    const double cross[3][3] = {
        {0., -u[2], u[1]}, {u[2], 0., -u[0]}, {-u[1], u[0], 0.}};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j) {
        double v = (i == j ? 1. : 0.) + cross[i][j];
        for (size_t k = 0; k < 3; ++k)
          v += cross[i][k] * cross[k][j] / (1. + g[0]);
        grad_rotation.get(j, i)[p] = v;
      }
    double orientation = 0.;
    for (size_t i = 0; i < 3; ++i) {
      const size_t j = (i + 1) % 3, k = (i + 2) % 3;
      orientation +=
          (C{0., 1.} * (frame.dyad[j][p] * std::conj(frame.dyad[k][p]) -
                        frame.dyad[k][p] * std::conj(frame.dyad[j][p])))
              .real() *
          frame.map.direction.get(i)[p];
    }
    if (std::abs(std::abs(orientation) - 1.) > 1.e-8)
      ERROR("Invalid longitudinal map dyad orientation at " << p);
    sigma[p] = std::copysign(1., orientation);
    for (size_t i = 0; i < 3; ++i) {
      aux_rotation.get(0, i)[p] = sigma[p] * frame.map.direction.get(i)[p];
      aux_rotation.get(1, i)[p] = sqrt(2.) * frame.dyad[i][p].real();
      aux_rotation.get(2, i)[p] = sqrt(2.) * frame.dyad[i][p].imag();
    }
    const double f = 1. - 2. * mass / get(r.measured_radius)[p];
    g_e[p] =
        sqrt(f) * (1. + 2. * mass / get(r.measured_radius)[p]) + 1. / sqrt(f);
  }
  const auto gradient_seed = weyl_scalars_from_tidal_tensor(rotate_symmetric(
      tidal_tensor_from_weyl_scalars(rest_seed), grad_rotation));
  const ComplexDataVector sensor =
      gradient_seed.get(3) - conj(gradient_seed.get(1));
  std::array<WeylScalars, 10> completion_columns;
  std::vector<double> design(3 * n * 10), rhs(3 * n);
  for (size_t p = 0; p < n; ++p) {
    const double w = frame.map.weights[p];
    if (not(std::isfinite(w) and w > 0.))
      ERROR("Invalid longitudinal quadrature weight");
    rhs[3 * p] = sqrt(w) * sensor[p].real() / (sqrt(2.) * g_e[p]);
    rhs[3 * p + 1] = sqrt(w) * sensor[p].imag() / (sqrt(2.) * g_e[p]);
    rhs[3 * p + 2] = sqrt(w) * 2. * get(seed.coulomb)[p].imag();
  }
  for (size_t a = 0; a < 10; ++a) {
    std::array<double, 5> unit{}, empty{};
    unit[a % 5] = 1.;
    const auto basis = stf_from_components(unit),
               none = stf_from_components(empty);
    const auto q = quadrupole_tide_tensor(
        a < 5 ? basis : none, a < 5 ? none : basis, frame.map.direction,
        r.measured_radius, mass, false);
    const auto rest_column =
        weyl_scalars_from_tidal_tensor(rotate_symmetric(q, aux_rotation));
    const auto column = push_forward(
        type_iii(rest_column, Scalar<DataVector>{-frame.boost}, zero),
        r.rotation.a_bar, r.rotation.b);
    const auto ds = longitudinal_seed_response(seed, column);
    auto correction = ds;
    for (size_t k = 0; k < 5; ++k)
      correction.get(k) -= column.get(k);
    const auto z = type_iii(
        pull_back(correction, r.rotation.a_bar, r.rotation.b), boost, zero);
    completion_columns[a] = column;
    for (size_t k = 0; k < 5; ++k)
      completion_columns[a].get(k) -= ds.get(k);
    for (size_t p = 0; p < n; ++p) {
      C electric = 0.;
      if (a < 5)
        for (size_t i = 0; i < 3; ++i)
          for (size_t j = 0; j < 3; ++j)
            electric += sqrt(2.) * sigma[p] * g_e[p] * basis.get(i, j) *
                        frame.map.direction.get(i)[p] *
                        std::conj(frame.dyad[j][p]);
      electric += z.get(3)[p] - std::conj(z.get(1)[p]);
      const C c = get(seed.coulomb)[p], p1 = seed.psi.get(1)[p],
              p2 = seed.psi.get(2)[p], p3 = seed.psi.get(3)[p];
      const C dc = ((4. * p2 - c) * column.get(2)[p] -
                    2. * p3 * column.get(1)[p] - 2. * p1 * column.get(3)[p]) /
                   (2. * c + p2);
      const double w = sqrt(frame.map.weights[p]);
      design[(3 * p) * 10 + a] = w * electric.real() / (sqrt(2.) * g_e[p]);
      design[(3 * p + 1) * 10 + a] = w * electric.imag() / (sqrt(2.) * g_e[p]);
      design[(3 * p + 2) * 10 + a] = w * 2. * dc.imag();
    }
  }
  std::array<double, 10> norms{};
  auto normalized = design;
  for (size_t a = 0; a < 10; ++a) {
    for (size_t p = 0; p < 3 * n; ++p)
      norms[a] = std::hypot(norms[a], design[p * 10 + a]);
    if (not(std::isfinite(norms[a]) and norms[a] > 0.))
      ERROR("Empty longitudinal design column");
    for (size_t p = 0; p < 3 * n; ++p)
      normalized[p * 10 + a] /= norms[a];
  }
  std::array<double, 100> v{};
  std::array<double, 10> singular{}, work{}, solution{};
  auto av = gsl_matrix_view_array(normalized.data(), 3 * n, 10);
  auto vv = gsl_matrix_view_array(v.data(), 10, 10);
  auto sv = gsl_vector_view_array(singular.data(), 10),
       wv = gsl_vector_view_array(work.data(), 10),
       bv = gsl_vector_view_array(rhs.data(), 3 * n),
       xv = gsl_vector_view_array(solution.data(), 10);
  const int status =
      gsl_linalg_SV_decomp(&av.matrix, &vv.matrix, &sv.vector, &wv.vector);
  if (status != 0 or not(singular[9] > 1.e-12 * singular[0]))
    ERROR("Rank-deficient longitudinal fit");
  result.normalized_condition_number = singular[0] / singular[9];
  if (gsl_linalg_SV_solve(&av.matrix, &vv.matrix, &sv.vector, &bv.vector,
                          &xv.vector) != 0)
    ERROR("Longitudinal SVD solve failed");
  for (size_t a = 0; a < 10; ++a) {
    solution[a] /= norms[a];
    if (not std::isfinite(solution[a]))
      ERROR("Nonfinite longitudinal moment");
  }
  double residual = 0., data_norm = 0.;
  for (size_t p = 0; p < 3 * n; ++p) {
    double value = -rhs[p];
    for (size_t a = 0; a < 10; ++a)
      value += design[p * 10 + a] * solution[a];
    residual = std::hypot(residual, value);
    data_norm = std::hypot(data_norm, rhs[p]);
  }
  result.relative_residual = data_norm > 0. ? residual / data_norm : 0.;
  result.completed = seed.psi;
  for (size_t a = 0; a < 5; ++a)
    result.moments[a] = C{solution[a], solution[a + 5]};
  for (size_t a = 0; a < 10; ++a)
    for (size_t k = 0; k < 5; ++k)
      result.completed.get(k) += solution[a] * completion_columns[a].get(k);
  for (const auto& field : result.completed)
    for (const auto value : field)
      if (not finite(value))
        ERROR("Nonfinite longitudinal target");
  return result;
}
}  // namespace gr::np
