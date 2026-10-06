// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>
#include <limits>
#include <memory>
#include <vector>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gr::np {
namespace {
using V3 = std::array<double, 3>;
using V4 = std::array<double, 4>;
using C4 = std::array<std::complex<double>, 4>;
using M3 = std::array<V3, 3>;
double dot(const V3& a, const V3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
V3 cross(const V3& a, const V3& b) {
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}
double dot4(const V4& a, const V4& b) {
  return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}
double determinant(const M3& a) { return dot(a[0], cross(a[1], a[2])); }

// Cached scalar harmonics depend only on the fixed NR labels, not geometry.
struct Basis {
  size_t lmax{}, n{};
  std::vector<V3> labels;
  std::vector<double> y, dt, dp;
  Basis(const TriadVector& dirs, size_t lm) : lmax(lm), n((lm + 1) * (lm + 1)) {
    const size_t np = dirs.get(0).size();
    labels.resize(np);
    y.resize(np * n);
    dt.resize(np * n);
    dp.resize(np * n);
    for (size_t p = 0; p < np; ++p) {
      for (size_t i = 0; i < 3; ++i)
        labels[p][i] = dirs.get(i)[p];
      const double x = std::clamp(labels[p][2], -1., 1.);
      const double st = std::hypot(labels[p][0], labels[p][1]);
      if (st < 1.e-12)
        ERROR("Eigenmap requires a pole-free angular grid");
      const double phi = std::atan2(labels[p][1], labels[p][0]);
      size_t a = 0;
      for (size_t l = 0; l <= lm; ++l)
        for (size_t m = 0; m <= l; ++m) {
          const double plm = gsl_sf_legendre_sphPlm(static_cast<int>(l),
                                                    static_cast<int>(m), x);
          const double previous =
              l > m ? gsl_sf_legendre_sphPlm(static_cast<int>(l - 1),
                                             static_cast<int>(m), x)
                    : 0.;
          const double dl = static_cast<double>(l), dm = static_cast<double>(m);
          const double factor =
              l > m
                  ? sqrt((dl * dl - dm * dm) * (2. * dl + 1.) / (2. * dl - 1.))
                  : 0.;
          const double derivative = (dl * x * plm - factor * previous) / st;
          const double scale = m == 0 ? 1. : sqrt(2.);
          const double c = std::cos(dm * phi), s = std::sin(dm * phi);
          y[p * n + a] = scale * plm * c;
          dt[p * n + a] = scale * derivative * c;
          dp[p * n + a] = -scale * plm * dm * s / st;
          ++a;
          if (m != 0) {
            y[p * n + a] = scale * plm * s;
            dt[p * n + a] = scale * derivative * s;
            dp[p * n + a] = scale * plm * dm * c / st;
            ++a;
          }
        }
    }
  }
  bool matches(const TriadVector& d, size_t lm) const {
    if (lmax != lm or labels.size() != d.get(0).size())
      return false;
    for (size_t p = 0; p < labels.size(); ++p)
      for (size_t i = 0; i < 3; ++i)
        if (std::abs(labels[p][i] - d.get(i)[p]) > 1.e-13)
          return false;
    return true;
  }
};

M3 rotation(const std::vector<V3>& raw, const std::vector<V3>& labels,
            const DataVector& area) {
  M3 c{}, v{}, r{};
  for (size_t p = 0; p < raw.size(); ++p)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        c[i][j] += area[p] * labels[p][i] * raw[p][j];
  auto u = gsl_matrix_view_array(c[0].data(), 3, 3);
  auto vv = gsl_matrix_view_array(v[0].data(), 3, 3);
  std::array<double, 3> s{}, w{};
  auto ss = gsl_vector_view_array(s.data(), 3);
  auto ww = gsl_vector_view_array(w.data(), 3);
  const int info =
      gsl_linalg_SV_decomp(&u.matrix, &vv.matrix, &ss.vector, &ww.vector);
  if (info != 0)
    ERROR("Eigenmap alignment SVD failed: " << info);
  const double sign = determinant(c) * determinant(v) < 0. ? -1. : 1.;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        r[i][j] += c[i][k] * v[j][k] * (k == 2 ? sign : 1.);
  return r;
}
V3 rotate(const M3& r, const V3& v) {
  return {{dot(r[0], v), dot(r[1], v), dot(r[2], v)}};
}
}  // namespace

SpatialRotation sphere_map_rotation(const TriadVector& old_direction,
                                    const TriadVector& new_direction,
                                    const DataVector& weights) {
  std::vector<V3> old(weights.size()), current(weights.size());
  for (size_t p = 0; p < weights.size(); ++p)
    for (size_t i = 0; i < 3; ++i) {
      old[p][i] = old_direction.get(i)[p];
      current[p][i] = new_direction.get(i)[p];
    }
  return rotation(old, current, weights);
}

EigenSphereMap laplace_eigenmap(const TriadVector& labels,
                                const DataVector& weights,
                                const std::array<DataVector, 3>& h, size_t lm) {
  const size_t np = weights.size(), n = (lm + 1) * (lm + 1);
  if (lm < 2 or n > np or labels.get(0).size() != np)
    ERROR("Invalid eigenmap basis/grid size");
  static thread_local std::unique_ptr<Basis> cache;
  if (not cache or not cache->matches(labels, lm))
    cache = std::make_unique<Basis>(labels, lm);
  const auto& basis = *cache;
  DataVector density(np), area(np);
  double label_area = 0., metric_area = 0.;
  for (size_t p = 0; p < np; ++p) {
    const double det = h[0][p] * h[2][p] - h[1][p] * h[1][p];
    if (not(h[0][p] > 0. and det > 0. and std::isfinite(det) and
            weights[p] > 0.))
      ERROR("Invalid screen metric or quadrature at " << p);
    density[p] = sqrt(det);
    label_area += weights[p];
    metric_area += weights[p] * density[p];
  }
  const double mean_area = metric_area / label_area;
  const double scale = 4. * std::acos(-1.) / label_area;
  std::vector<double> stiffness(n * n, 0.), mass(n * n, 0.), values(n),
      vectors(n * n);
  for (size_t p = 0; p < np; ++p) {
    area[p] = scale * weights[p] * density[p] / mean_area;
    // sqrt(det h) h^{-1} is unchanged by the constant area normalization.
    const double factor = scale * weights[p] / density[p];
    for (size_t a = 0; a < n; ++a) {
      const double ta = basis.dt[p * n + a], pa = basis.dp[p * n + a],
                   ya = basis.y[p * n + a];
      const double ft = factor * (h[2][p] * ta - h[1][p] * pa);
      const double fp = factor * (h[0][p] * pa - h[1][p] * ta);
      for (size_t b = 0; b <= a; ++b) {
        stiffness[a * n + b] +=
            ft * basis.dt[p * n + b] + fp * basis.dp[p * n + b];
        mass[a * n + b] += area[p] * ya * basis.y[p * n + b];
      }
    }
  }
  for (size_t a = 0; a < n; ++a)
    for (size_t b = 0; b < a; ++b) {
      stiffness[b * n + a] = stiffness[a * n + b];
      mass[b * n + a] = mass[a * n + b];
    }
  auto aa = gsl_matrix_view_array(stiffness.data(), n, n),
       bb = gsl_matrix_view_array(mass.data(), n, n);
  auto vv = gsl_matrix_view_array(vectors.data(), n, n);
  auto ev = gsl_vector_view_array(values.data(), n);
  auto* work = gsl_eigen_gensymmv_alloc(n);
  const int info =
      gsl_eigen_gensymmv(&aa.matrix, &bb.matrix, &ev.vector, &vv.matrix, work);
  gsl_eigen_gensymmv_free(work);
  if (info != 0)
    ERROR("Screen Laplacian generalized eigensolve failed: " << info);
  gsl_eigen_gensymmv_sort(&ev.vector, &vv.matrix, GSL_EIGEN_SORT_VAL_ASC);
  if (n < 5 or values[1] <= 0. or values[4] <= values[3])
    ERROR("Screen eigenmap triplet is not isolated");
  const double normalization = sqrt(4. * std::acos(-1.) / 3.);
  std::vector<V3> raw(np), d0(np), d1(np);
  double handedness = 0.;
  V3 integral_norm{};
  for (size_t p = 0; p < np; ++p) {
    for (size_t a = 0; a < n; ++a)
      for (size_t i = 0; i < 3; ++i) {
        const double c = vectors[a * n + i + 1] * normalization;
        raw[p][i] += basis.y[p * n + a] * c;
        d0[p][i] += basis.dt[p * n + a] * c;
        d1[p][i] += basis.dp[p * n + a] * c;
      }
    for (size_t i = 0; i < 3; ++i) {
      integral_norm[i] += area[p] * raw[p][i] * raw[p][i];
    }
  }
  // GSL returns coefficient vectors of Euclidean norm one, whereas the
  // embedding needs integral_h F_i F_j = (4 pi / 3) delta_ij. Unlike an
  // overall scaling, unequal component norms survive pointwise normalization.
  for (size_t i = 0; i < 3; ++i) {
    const double factor = normalization / sqrt(integral_norm[i]);
    for (size_t p = 0; p < np; ++p) {
      raw[p][i] *= factor;
      d0[p][i] *= factor;
      d1[p][i] *= factor;
    }
  }
  for (size_t p = 0; p < np; ++p) {
    handedness += area[p] * dot(raw[p], cross(d0[p], d1[p]));
  }
  if (handedness < 0.)
    for (size_t p = 0; p < np; ++p) {
      raw[p][2] *= -1.;
      d0[p][2] *= -1.;
      d1[p][2] *= -1.;
    }
  const M3 r = rotation(raw, basis.labels, area);
  EigenSphereMap result{TriadVector(np, 0.),
                        {{TriadVector(np, 0.), TriadVector(np, 0.)}},
                        DataVector(np, 0.),
                        {},
                        std::numeric_limits<double>::infinity(),
                        0.};
  std::copy_n(values.begin(), 5, result.eigenvalues.begin());
  for (size_t p = 0; p < np; ++p) {
    V3 v = rotate(r, raw[p]), a = rotate(r, d0[p]), b = rotate(r, d1[p]);
    const double length = sqrt(dot(v, v));
    if (not(length > 1.e-12))
      ERROR("Eigenmap direction vanishes at " << p);
    for (auto& x : v)
      x /= length;
    const double va = dot(v, a), vb = dot(v, b);
    for (size_t i = 0; i < 3; ++i) {
      a[i] = (a[i] - v[i] * va) / length;
      b[i] = (b[i] - v[i] * vb) / length;
    }
    const double jac = dot(v, cross(a, b));
    if (not(jac > 0. and std::isfinite(jac)))
      ERROR("Eigenmap folds at " << p << ": Jacobian " << jac);
    result.minimum_jacobian = std::min(result.minimum_jacobian, jac);
    result.weights[p] = weights[p] * jac;
    result.area_over_label_area += result.weights[p] / label_area;
    for (size_t i = 0; i < 3; ++i) {
      result.direction.get(i)[p] = v[i];
      result.derivative[0].get(i)[p] = a[i];
      result.derivative[1].get(i)[p] = b[i];
    }
  }
  return result;
}

std::pair<DataVector, double> sphere_gradient_potential(
    const TriadVector& labels, const DataVector& weights,
    const std::array<DataVector, 2>& gradient, const size_t lm) {
  const Basis basis(labels, lm);
  const size_t np = weights.size(), n = basis.n - 1;
  std::vector<double> normal(n * n, 0.), rhs(n, 0.), coefficients(n);
  for (size_t p = 0; p < np; ++p) {
    for (size_t a = 0; a < n; ++a) {
      const double ta = basis.dt[p * basis.n + a + 1];
      const double pa = basis.dp[p * basis.n + a + 1];
      rhs[a] += weights[p] * (ta * gradient[0][p] + pa * gradient[1][p]);
      for (size_t b = 0; b < n; ++b) {
        normal[a * n + b] += weights[p] * (ta * basis.dt[p * basis.n + b + 1] +
                                           pa * basis.dp[p * basis.n + b + 1]);
      }
    }
  }
  auto matrix = gsl_matrix_view_array(normal.data(), n, n);
  auto right = gsl_vector_view_array(rhs.data(), n);
  auto solution = gsl_vector_view_array(coefficients.data(), n);
  if (gsl_linalg_cholesky_decomp(&matrix.matrix) != 0 or
      gsl_linalg_cholesky_solve(&matrix.matrix, &right.vector,
                                &solution.vector) != 0) {
    ERROR("Slice-time potential solve failed");
  }
  DataVector potential(np, 0.);
  double mean = 0., area = 0., error = 0., norm = 0.;
  for (size_t p = 0; p < np; ++p) {
    std::array<double, 2> fitted{};
    for (size_t a = 0; a < n; ++a) {
      potential[p] += coefficients[a] * basis.y[p * basis.n + a + 1];
      fitted[0] += coefficients[a] * basis.dt[p * basis.n + a + 1];
      fitted[1] += coefficients[a] * basis.dp[p * basis.n + a + 1];
    }
    mean += weights[p] * potential[p];
    area += weights[p];
    for (size_t a = 0; a < 2; ++a) {
      error += weights[p] * pow(fitted[a] - gradient[a][p], 2);
      norm += weights[p] * pow(gradient[a][p], 2);
    }
  }
  potential -= mean / area;
  return {std::move(potential), sqrt(error / std::max(norm, 1.e-300))};
}

GeometricFrame geometric_frame(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& rotation_matrix,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& metric,
    const TriadVector& labels, const DataVector& weights, double mass,
    size_t lm) {
  const size_t np = weights.size();
  const double is2 = 1. / sqrt(2.);
  const auto lower = cholesky_factor(metric);
  std::vector<std::array<V4, 2>> screen(np);
  std::vector<C4> dyad(np);
  GeometricFrame result{};
  result.boost = DataVector(np);
  result.observer = AdaptedFourVector(np, 0.);
  for (auto& field : result.screen)
    field = AdaptedFourVector(np, 0.);
  for (auto& field : result.dyad)
    field = ComplexDataVector(np, 0.);
  std::array<DataVector, 3> gram{
      {DataVector(np), DataVector(np), DataVector(np)}},
      h = gram;
  for (size_t p = 0; p < np; ++p) {
    C4 ell{{is2, is2, 0., 0.}}, kay{{is2, -is2, 0., 0.}},
        emm{{0., 0., is2, {0., is2}}};
    const auto a = -std::conj(get(registration.rotation.a_bar)[p]);
    const auto b = -get(registration.rotation.b)[p];
    C4 nk{}, nm{};
    for (size_t j = 0; j < 4; ++j) {
      nk[j] = kay[j] + std::conj(a) * emm[j] + a * std::conj(emm[j]) +
              std::norm(a) * ell[j];
      nm[j] = emm[j] + a * ell[j];
    }
    kay = nk;
    emm = nm;
    for (size_t j = 0; j < 4; ++j) {
      ell[j] +=
          std::conj(b) * emm[j] + b * std::conj(emm[j]) + std::norm(b) * kay[j];
      nm[j] = emm[j] + b * kay[j];
    }
    emm = nm;
    dyad[p] = emm;
    result.boost[p] =
        .5 * log(kay[0].real() / ell[0].real()) + get(rapidity)[p];
    V4 l{}, k{};
    for (size_t j = 0; j < 4; ++j) {
      l[j] = ell[j].real();
      k[j] = kay[j].real();
      result.observer.get(j)[p] =
          is2 * (exp(result.boost[p]) * l[j] + exp(-result.boost[p]) * k[j]);
    }
    const V3 n{{labels.get(0)[p], labels.get(1)[p], labels.get(2)[p]}};
    const double st = std::hypot(n[0], n[1]);
    const V3 et{{n[2] * n[0] / st, n[2] * n[1] / st, -st}}, ep = cross(n, et);
    for (size_t A = 0; A < 2; ++A) {
      const V3 tangent = A == 0 ? et : ep;
      V4 x{};
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          for (size_t z = 0; z < 3; ++z)
            x[i + 1] +=
                rotation_matrix.get(i, j)[p] * lower.get(z, j)[p] * tangent[z];
      const double lx = dot4(l, x), kx = dot4(k, x);
      for (size_t j = 0; j < 4; ++j) {
        screen[p][A][j] = x[j] + l[j] * kx + k[j] * lx;
        result.screen[A].get(j)[p] = screen[p][A][j];
      }
    }
    gram[0][p] = dot4(screen[p][0], screen[p][0]);
    gram[1][p] = dot4(screen[p][0], screen[p][1]);
    gram[2][p] = dot4(screen[p][1], screen[p][1]);
    const double rr = get(registration.measured_radius)[p];
    if (not(mass > 0. and rr > 2. * mass and std::isfinite(rr))) {
      ERROR("Geometric tide requires positive mass and radius outside 2M");
    }
    for (size_t j = 0; j < 3; ++j)
      h[j][p] = gram[j][p] / (rr * rr);
  }
  result.map = laplace_eigenmap(labels, weights, h, lm);
  for (size_t p = 0; p < np; ++p) {
    const double g0 = sqrt(gram[0][p]), ratio = gram[1][p] / gram[0][p];
    const double g1 = sqrt(gram[2][p] - gram[1][p] * ratio);
    V3 z0{}, z1{};
    std::complex<double> m0 = 0., m1 = 0.;
    for (size_t j = 0; j < 4; ++j) {
      const double sign = j == 0 ? -1. : 1.;
      m0 += sign * screen[p][0][j] / g0 * dyad[p][j];
      m1 +=
          sign * (screen[p][1][j] - ratio * screen[p][0][j]) / g1 * dyad[p][j];
    }
    for (size_t i = 0; i < 3; ++i) {
      z0[i] = result.map.derivative[0].get(i)[p] / g0;
      z1[i] = (result.map.derivative[1].get(i)[p] -
               ratio * result.map.derivative[0].get(i)[p]) /
              g1;
    }
    const double c00 = dot(z0, z0), c01 = dot(z0, z1), c11 = dot(z1, z1);
    const double det = c00 * c11 - c01 * c01, s = sqrt(det),
                 t = sqrt(c00 + c11 + 2. * s);
    if (not(det > 0. and std::isfinite(t)))
      ERROR("Singular polar dyad at " << p);
    const double denom = (c00 + s) * (c11 + s) - c01 * c01;
    const double i00 = t * (c11 + s) / denom, i01 = -t * c01 / denom,
                 i11 = t * (c00 + s) / denom;
    std::array<std::complex<double>, 3> mu{};
    std::complex<double> null = 0.;
    double norm = 0.;
    for (size_t i = 0; i < 3; ++i) {
      mu[i] =
          (i00 * z0[i] + i01 * z1[i]) * m0 + (i01 * z0[i] + i11 * z1[i]) * m1;
      result.dyad[i][p] = mu[i];
      null += mu[i] * mu[i];
      norm += std::norm(mu[i]);
    }
    result.maximum_dyad_error = std::max(
        {result.maximum_dyad_error, std::abs(null), std::abs(norm - 1.)});
  }
  if (result.maximum_dyad_error > 1.e-9)
    ERROR("Polar dyad normalization failed: " << result.maximum_dyad_error);
  return result;
}

GeometricTideEvaluation evaluate_geometric_second_order(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& rotation_matrix,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& metric,
    const TriadVector& labels, const DataVector& weights, double mass,
    size_t lm, const std::optional<TidalMoments>& imposed) {
  const auto frame = geometric_frame(registration, rapidity, rotation_matrix,
                                     metric, labels, weights, mass, lm);
  return evaluate_geometric_second_order(registration, frame, mass, imposed);
}

GeometricTideEvaluation evaluate_geometric_second_order(
    const FrameRegistration& registration, const GeometricFrame& frame,
    const double mass, const std::optional<TidalMoments>& imposed) {
  const size_t np = frame.map.weights.size();
  GeometricTideEvaluation result{};
  result.map = frame.map;
  result.maximum_dyad_error = frame.maximum_dyad_error;
  std::array<Scalar<ComplexDataVector>, 5> c0{}, c4{};
  for (size_t a = 0; a < 5; ++a) {
    get(c0[a]) = ComplexDataVector(np, 0.);
    get(c4[a]) = ComplexDataVector(np, 0.);
  }
  for (size_t p = 0; p < np; ++p) {
    const std::array<std::complex<double>, 3> mu{
        {frame.dyad[0][p], frame.dyad[1][p], frame.dyad[2][p]}};
    const std::array<std::complex<double>, 5> contraction{
        {mu[0] * mu[0] - mu[2] * mu[2], mu[1] * mu[1] - mu[2] * mu[2],
         2. * mu[0] * mu[1], 2. * mu[0] * mu[2], 2. * mu[1] * mu[2]}};
    const double f = 1. - 2. * mass / get(registration.measured_radius)[p];
    if (not(f > 0.))
      ERROR("Geometric tide requires an outside-horizon radius");
    for (size_t a = 0; a < 5; ++a) {
      get(c0[a])[p] = f * exp(-2. * frame.boost[p]) * contraction[a];
      get(c4[a])[p] = f * exp(2. * frame.boost[p]) * std::conj(contraction[a]);
    }
  }
  auto& ev = result.second_order;
  const Scalar<ComplexDataVector> data{registration.pulled_back.get(4)};
  if (imposed) {
    ev.fit.components = *imposed;
    ev.fit.relative_residual =
        relative_fit_residual(data, c4, *imposed, result.map.weights);
  } else
    ev.fit = fit_psi4(data, c4, result.map.weights);
  auto replacement = registration.pulled_back;
  replacement.get(0) = 0.;
  for (size_t a = 0; a < 5; ++a) {
    replacement.get(0) += ev.fit.components[a] * get(c0[a]);
    WeylScalars column(np, std::complex<double>{0., 0.});
    column.get(0) = get(c0[a]);
    column.get(4) = get(c4[a]);
    ev.direct_columns[a] = push_forward(column, registration.rotation.a_bar,
                                        registration.rotation.b);
  }
  get(ev.psi0_target) = push_forward(replacement, registration.rotation.a_bar,
                                     registration.rotation.b)
                            .get(0);
  return result;
}
}  // namespace gr::np
