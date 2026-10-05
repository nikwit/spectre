// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/ThirdOrderTide.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <limits>
#include <pup.h>
#include <utility>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gr::np {
namespace {
using V3 = std::array<double, 3>;
using M3 = std::array<V3, 3>;
using T3 = std::array<M3, 3>;
using C3 = std::array<std::complex<double>, 3>;
using V4 = std::array<double, 4>;
double dot4(const V4& a, const V4& b) {
  return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}
V3 cross(const V3& a, const V3& b) {
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}
M3 rank2(const std::array<double, 5>& x) {
  return {{{{x[0], x[2], x[3]}},
           {{x[2], x[1], x[4]}},
           {{x[3], x[4], -x[0] - x[1]}}}};
}
T3 rank3(const std::array<double, 7>& x) {
  T3 tensor{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k) {
        std::array<size_t, 3> a{{i, j, k}};
        std::sort(a.begin(), a.end());
        if (a == std::array<size_t, 3>{{0, 0, 0}})
          tensor[i][j][k] = x[0];
        else if (a == std::array<size_t, 3>{{0, 0, 1}})
          tensor[i][j][k] = x[1];
        else if (a == std::array<size_t, 3>{{0, 0, 2}})
          tensor[i][j][k] = x[2];
        else if (a == std::array<size_t, 3>{{0, 1, 1}})
          tensor[i][j][k] = x[3];
        else if (a == std::array<size_t, 3>{{0, 1, 2}})
          tensor[i][j][k] = x[4];
        else if (a == std::array<size_t, 3>{{0, 2, 2}})
          tensor[i][j][k] = -x[0] - x[3];
        else if (a == std::array<size_t, 3>{{1, 1, 1}})
          tensor[i][j][k] = x[5];
        else if (a == std::array<size_t, 3>{{1, 1, 2}})
          tensor[i][j][k] = x[6];
        else if (a == std::array<size_t, 3>{{1, 2, 2}})
          tensor[i][j][k] = -x[1] - x[5];
        else
          tensor[i][j][k] = -x[2] - x[6];
      }
  return tensor;
}
std::complex<double> contract(const M3& h, const C3& m) {
  std::complex<double> value = 0.;
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      value += h[i][j] * m[i] * m[j];
  return value;
}
// eps_{kli} H_jk N_l, symmetrized. Contracting with m_i m_j performs
// the symmetrization, and m.N=m.m=0 removes all nontransverse pieces.
M3 induction(const M3& h, const V3& n) {
  M3 result{};
  for (size_t j = 0; j < 3; ++j) {
    const auto v = cross(h[j], n);
    for (size_t i = 0; i < 3; ++i)
      result[i][j] = v[i];
  }
  return result;
}
}  // namespace

ThirdOrderColumns third_order_tide_columns(const GeometricFrame& frame,
                                           const Scalar<DataVector>& radius,
                                           const double mass,
                                           const DataVector& tilt) {
  const size_t np = get(radius).size();
  ThirdOrderColumns result{};
  for (auto& column : result)
    for (auto& slot : column)
      slot = ComplexDataVector(np, 0.);
  const std::complex<double> imaginary{0., 1.};
  for (size_t p = 0; p < np; ++p) {
    const double r = get(radius)[p], c = mass / r, f = 1. - 2. * c;
    if (not(mass > 0. and f > 0. and std::isfinite(r)))
      ERROR("Third-order tide requires positive mass and radius outside 2M");
    const double e_ind = 1. - c / (4. * f) + .75 * c * f,
                 b_ind = 1. - c / (4. * f);
    const double w = mass * (f * (-92. / 15. + 2. * log(f)) +
                             (52. / 15. - 148. * c / 15. + 28. * c * c / 15. +
                              16. * pow(c, 3) / 3. + 8. * pow(c, 4) / 3.) /
                                 f);
    const double wb = mass * (f * (-76. / 15. + 2. * log(f)) +
                              (29. / 10. - 38. * c / 5. - 2. * c * c / 5. +
                               16. * pow(c, 3) / 3. + 8. * pow(c, 4) / 3.) /
                                  f);
    const V3 n{{frame.map.direction.get(0)[p], frame.map.direction.get(1)[p],
                frame.map.direction.get(2)[p]}};
    for (size_t slot = 0; slot < 2; ++slot) {
      C3 m{};
      for (size_t i = 0; i < 3; ++i)
        m[i] = slot == 0 ? frame.dyad[i][p] : conj(frame.dyad[i][p]);
      const double boost = exp((slot == 0 ? -2. : 2.) * frame.boost[p]);
      for (size_t a = 0; a < 5; ++a) {
        std::array<double, 5> basis{};
        basis[a] = 1.;
        const auto h = rank2(basis);
        const auto q = contract(h, m), ind = contract(induction(h, n), m);
        result[a][slot][p] = boost * f * q;
        result[a + 5][slot][p] = imaginary * result[a][slot][p];
        result[a + 24][slot][p] =
            boost *
            ((w + tilt[p] * f) * q + imaginary * (2. * r / 3.) * e_ind * ind);
        result[a + 29][slot][p] = boost * (imaginary * (wb + tilt[p] * f) * q -
                                           (2. * r / 3.) * b_ind * ind);
      }
      for (size_t a = 0; a < 7; ++a) {
        std::array<double, 7> basis{};
        basis[a] = 1.;
        const auto h = rank3(basis);
        M3 contracted{};
        for (size_t i = 0; i < 3; ++i)
          for (size_t j = 0; j < 3; ++j)
            for (size_t k = 0; k < 3; ++k)
              contracted[i][j] += h[i][j][k] * n[k];
        result[a + 10][slot][p] =
            boost * r * f * (1. - c) * contract(contracted, m);
        result[a + 17][slot][p] =
            (4. / 3.) * imaginary * result[a + 10][slot][p];
      }
    }
  }
  return result;
}

ThirdOrderFit fit_third_order_tide(const FrameRegistration& registration,
                                   const ThirdOrderColumns& columns,
                                   const DataVector& weights,
                                   const DottedTidalMoments& dots) {
  constexpr size_t nc = 24;
  const size_t np = weights.size();
  if (2 * np < nc)
    ERROR("Third-order fit needs at least 12 complex data points");
  std::vector<double> design(2 * np * nc), v(nc * nc), s(nc), work(nc),
      rhs(2 * np), x(nc);
  std::array<double, nc> scale{};
  for (size_t p = 0; p < np; ++p) {
    if (not(weights[p] > 0.))
      ERROR("Invalid third-order quadrature");
    auto data = registration.pulled_back.get(4)[p];
    for (size_t a = 0; a < 10; ++a)
      data -= dots[a] * columns[a + 24][1][p];
    rhs[p] = sqrt(weights[p]) * data.real();
    rhs[p + np] = sqrt(weights[p]) * data.imag();
    for (size_t a = 0; a < nc; ++a) {
      const auto value = sqrt(weights[p]) * columns[a][1][p];
      design[p * nc + a] = value.real();
      design[(p + np) * nc + a] = value.imag();
      scale[a] += std::norm(value);
    }
  }
  for (size_t a = 0; a < nc; ++a) {
    scale[a] = sqrt(scale[a]);
    if (not(scale[a] > 0. and std::isfinite(scale[a])))
      ERROR("Empty or nonfinite third-order design column");
    for (size_t p = 0; p < 2 * np; ++p)
      design[p * nc + a] /= scale[a];
  }
  auto u = gsl_matrix_view_array(design.data(), 2 * np, nc);
  auto vv = gsl_matrix_view_array(v.data(), nc, nc);
  auto ss = gsl_vector_view_array(s.data(), nc);
  auto ww = gsl_vector_view_array(work.data(), nc);
  auto bb = gsl_vector_view_array(rhs.data(), 2 * np);
  auto xx = gsl_vector_view_array(x.data(), nc);
  if (gsl_linalg_SV_decomp(&u.matrix, &vv.matrix, &ss.vector, &ww.vector) != 0)
    ERROR("Third-order design SVD failed");
  if (not(s.back() > 1.e-12 * s.front()))
    ERROR("Third-order Psi4 fit is rank deficient");
  if (gsl_linalg_SV_solve(&u.matrix, &vv.matrix, &ss.vector, &bb.vector,
                          &xx.vector) != 0)
    ERROR("Third-order Psi4 solve failed");
  ThirdOrderFit result{};
  result.condition_number = s.front() / s.back();
  for (size_t a = 0; a < nc; ++a)
    result.moments[a] = x[a] / scale[a];
  auto replacement = registration.pulled_back;
  replacement.get(0) = 0.;
  double error = 0., norm = 0.;
  for (size_t p = 0; p < np; ++p) {
    std::complex<double> outgoing = 0.;
    for (size_t a = 0; a < 34; ++a) {
      const double coefficient = a < 24 ? result.moments[a] : dots[a - 24];
      replacement.get(0)[p] += coefficient * columns[a][0][p];
      outgoing += coefficient * columns[a][1][p];
    }
    error +=
        weights[p] * std::norm(outgoing - registration.pulled_back.get(4)[p]);
    norm += weights[p] * std::norm(registration.pulled_back.get(4)[p]);
  }
  result.relative_residual = sqrt(error / std::max(norm, 1.e-300));
  get(result.psi0_target) =
      push_forward(replacement, registration.rotation.a_bar,
                   registration.rotation.b)
          .get(0);
  return result;
}

GeometricTimeData geometric_time_data(
    const GeometricFrame& frame, const TriadVector& labels, const double radius,
    const Scalar<DataVector>& measured_radius, const double mass,
    const RealMatrix& rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& metric,
    const Scalar<DataVector>& lapse, const TriadVector& shift,
    const TriadVector& velocity) {
  const size_t np = frame.map.weights.size();
  const auto lower = cholesky_factor(metric);
  GeometricTimeData result{};
  result.minimum_clock_rate = std::numeric_limits<double>::infinity();
  for (auto& component : result.flow)
    component = DataVector(np, 0.);
  std::array<DataVector, 2> gradient{{DataVector(np, 0.), DataVector(np, 0.)}};
  double area = 0.;
  for (size_t p = 0; p < np; ++p) {
    const double f = 1. - 2. * mass / get(measured_radius)[p];
    if (not(f > 0. and std::isfinite(f) and radius > 0.))
      ERROR("Invalid geometric model clock geometry");
    const double root_f = sqrt(f);
    M3 transform{};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        for (size_t k = 0; k < 3; ++k)
          transform[i][j] += rotation.get(i, k)[p] * lower.get(j, k)[p];
    V4 xi{{get(lapse)[p], 0., 0., 0.}}, observer{};
    for (size_t i = 0; i < 4; ++i)
      observer[i] = frame.observer.get(i)[p];
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        xi[i + 1] += transform[i][j] * (shift.get(j)[p] + velocity.get(j)[p]);
    const V3 n{{labels.get(0)[p], labels.get(1)[p], labels.get(2)[p]}};
    const double st = hypot(n[0], n[1]);
    const V3 et{{n[2] * n[0] / st, n[2] * n[1] / st, -st}}, ep = cross(n, et);
    std::array<V4, 2> screen{};
    std::array<double, 2> covector{};
    for (size_t a = 0; a < 2; ++a) {
      const auto& tangent = a == 0 ? et : ep;
      V4 original{};
      for (size_t i = 0; i < 4; ++i)
        screen[a][i] = radius * frame.screen[a].get(i)[p];
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          original[i + 1] += radius * transform[i][j] * tangent[j];
      gradient[a][p] = -dot4(observer, original) / root_f;
      covector[a] = dot4(screen[a], xi);
    }
    const double g00 = dot4(screen[0], screen[0]),
                 g01 = dot4(screen[0], screen[1]),
                 g11 = dot4(screen[1], screen[1]), det = g00 * g11 - g01 * g01;
    if (not(det > 0.))
      ERROR("Singular no-screen-flow metric");
    result.flow[0][p] = (-g11 * covector[0] + g01 * covector[1]) / det;
    result.flow[1][p] = (g01 * covector[0] - g00 * covector[1]) / det;
    const double clock = -dot4(observer, xi) / root_f +
                         result.flow[0][p] * gradient[0][p] +
                         result.flow[1][p] * gradient[1][p];
    if (not(clock > 0. and std::isfinite(clock)))
      ERROR("Third-order matching requires a positive finite model clock");
    result.minimum_clock_rate = std::min(result.minimum_clock_rate, clock);
    result.clock_rate += frame.map.weights[p] * clock;
    area += frame.map.weights[p];
  }
  result.clock_rate /= area;
  auto potential =
      sphere_gradient_potential(labels, frame.map.weights, gradient, 2);
  result.slice_tilt = std::move(potential.first);
  result.tilt_residual = potential.second;
  return result;
}

void GeometricTideSample::pup(PUP::er& p) {
  p | time;
  p | quadrupole;
  p | direction;
}
bool operator==(const GeometricTideSample& a, const GeometricTideSample& b) {
  return a.time == b.time and a.quadrupole == b.quadrupole and
         a.direction == b.direction;
}
void GeometricTideHistory::pup(PUP::er& p) { p | samples; }
bool operator==(const GeometricTideHistory& a, const GeometricTideHistory& b) {
  return a.samples == b.samples;
}

CausalTidalDerivative causal_tidal_derivative(
    const gsl::not_null<GeometricTideHistory*> history,
    const GeometricTideSample& current, const EigenSphereMap& map,
    const GeometricTimeData& time_data, const double step_start,
    const bool full_step) {
  auto& past = history->samples;
  const size_t np = map.weights.size();
  if (not(std::isfinite(current.time) and time_data.clock_rate > 0. and
          std::isfinite(time_data.clock_rate) and current.time >= step_start))
    ERROR("Third-order history requires forward time and a positive clock");
  if (not past.empty() and (past.back().time > step_start or
                            past.back().direction.get(0).size() != np))
    past.clear();
  // Repeated full-step calls replace the endpoint. They must use the same
  // strictly past stencil, including when the RHS is recomputed after a retry.
  if (full_step) {
    while (not past.empty() and past.back().time >= current.time)
      past.pop_back();
  }
  const size_t available = static_cast<size_t>(std::count_if(
      past.begin(), past.end(),
      [&current](const auto& sample) { return sample.time < current.time; }));
  const size_t count = std::min(size_t{4}, available);
  CausalTidalDerivative result{};
  result.derivative_order = count;
  if (count > 0) {
    std::vector<const GeometricTideSample*> nodes;
    for (size_t a = available - count; a < available; ++a)
      nodes.push_back(&past[a]);
    nodes.push_back(&current);
    const double scale = current.time - nodes.front()->time;
    std::vector<double> x(nodes.size()), weights(nodes.size(), 0.);
    for (size_t a = 0; a < nodes.size(); ++a)
      x[a] = (nodes[a]->time - current.time) / scale;
    // Derivatives at x=0 of the Lagrange polynomials. Scaling keeps
    // nonuniform small timesteps away from powers of an absolute time.
    for (size_t a = 0; a < nodes.size(); ++a)
      for (size_t k = 0; k < nodes.size(); ++k) {
        if (k == a)
          continue;
        double term = 1. / (x[a] - x[k]);
        for (size_t j = 0; j < nodes.size(); ++j)
          if (j != a and j != k)
            term *= -x[j] / (x[a] - x[j]);
        weights[a] += term / scale;
      }
    DottedTidalMoments raw{};
    TriadVector direction_dt(np, 0.);
    for (size_t a = 0; a < nodes.size(); ++a) {
      for (size_t i = 0; i < 10; ++i)
        raw[i] +=
            weights[a] * (nodes[a]->quadrupole[i] - current.quadrupole[i]);
      for (size_t i = 0; i < 3; ++i)
        direction_dt.get(i) += weights[a] * (nodes[a]->direction.get(i) -
                                             current.direction.get(i));
    }
    M3 inertia{};
    V3 rhs{};
    for (size_t p = 0; p < np; ++p) {
      V3 n{}, dn{};
      for (size_t i = 0; i < 3; ++i) {
        n[i] = current.direction.get(i)[p];
        dn[i] = direction_dt.get(i)[p];
        for (size_t a = 0; a < 2; ++a)
          dn[i] += time_data.flow[a][p] * map.derivative[a].get(i)[p];
      }
      const auto angular = cross(n, dn);
      for (size_t i = 0; i < 3; ++i) {
        rhs[i] += map.weights[p] * angular[i];
        for (size_t j = 0; j < 3; ++j)
          inertia[i][j] += map.weights[p] * ((i == j ? 1. : 0.) - n[i] * n[j]);
      }
    }
    auto aa = gsl_matrix_view_array(inertia[0].data(), 3, 3);
    auto bb = gsl_vector_view_array(rhs.data(), 3);
    auto xx = gsl_vector_view_array(result.angular_velocity.data(), 3);
    if (gsl_linalg_cholesky_decomp(&aa.matrix) != 0 or
        gsl_linalg_cholesky_solve(&aa.matrix, &bb.vector, &xx.vector) != 0)
      ERROR("Geometric rotation connection solve failed");
    M3 omega{};
    const auto& w = result.angular_velocity;
    omega[0][1] = -w[2];
    omega[0][2] = w[1];
    omega[1][0] = w[2];
    omega[1][2] = -w[0];
    omega[2][0] = -w[1];
    omega[2][1] = w[0];
    const std::array<std::array<size_t, 2>, 5> indices{
        {{{0, 0}}, {{1, 1}}, {{0, 1}}, {{0, 2}}, {{1, 2}}}};
    for (size_t parity = 0; parity < 2; ++parity) {
      std::array<double, 5> values{};
      std::copy_n(current.quadrupole.begin() + 5 * parity, 5, values.begin());
      const auto h = rank2(values);
      for (size_t a = 0; a < 5; ++a) {
        const auto i = indices[a][0], j = indices[a][1];
        double connection = 0.;
        for (size_t k = 0; k < 3; ++k)
          connection += omega[i][k] * h[k][j] - h[i][k] * omega[k][j];
        result.dots[a + 5 * parity] =
            (raw[a + 5 * parity] - connection) / time_data.clock_rate;
      }
    }
  }
  if (full_step) {
    past.push_back(current);
    if (past.size() > 5)
      past.erase(past.begin());
  }
  return result;
}
}  // namespace gr::np
