// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTideTransport.hpp"

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
using V4 = std::array<double, 4>;
using M3 = SpatialRotation;
double dot4(const V4& a, const V4& b) {
  return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}
V3 cross(const V3& a, const V3& b) {
  return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0]}};
}
M3 multiply(const M3& a, const M3& b) {
  M3 result{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        result[i][j] += a[i][k] * b[k][j];
  return result;
}
M3 rotation_step(const V3& rate, const double dt) {
  const V3 v{{dt * rate[0], dt * rate[1], dt * rate[2]}};
  const double a2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  const double a = sqrt(a2);
  const double s = a < 1.e-4 ? 1. - a2 / 6. + a2 * a2 / 120. : sin(a) / a;
  const double c =
      a < 1.e-4 ? .5 - a2 / 24. + a2 * a2 / 720. : (1. - cos(a)) / a2;
  const M3 skew{
      {{{0., -v[2], v[1]}}, {{v[2], 0., -v[0]}}, {{-v[1], v[0], 0.}}}};
  auto result = multiply(skew, skew);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      result[i][j] = c * result[i][j] + s * skew[i][j] + (i == j ? 1. : 0.);
  return result;
}
TidalMoments rotate(const M3& r, const TidalMoments& x) {
  const std::array<std::array<std::complex<double>, 3>, 3> h{
      {{{x[0], x[2], x[3]}},
       {{x[2], x[1], x[4]}},
       {{x[3], x[4], -x[0] - x[1]}}}};
  TidalMoments result{};
  const std::array<std::array<size_t, 2>, 5> indices{
      {{{0, 0}}, {{1, 1}}, {{0, 1}}, {{0, 2}}, {{1, 2}}}};
  for (size_t a = 0; a < 5; ++a)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        result[a] += r[indices[a][0]][i] * h[i][j] * r[indices[a][1]][j];
  return result;
}
}  // namespace

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
      ERROR("Geometric matching requires a positive finite model clock");
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

std::array<double, 3> geometric_rotation_rate(
    const EigenSphereMap& map, const TriadVector& direction_dt,
    const std::array<DataVector, 2>& flow) {
  const size_t np = map.weights.size();
  std::array<double, 3> result{};
  M3 inertia{};
  V3 rhs{};
  for (size_t p = 0; p < np; ++p) {
    V3 n{}, dn{};
    for (size_t i = 0; i < 3; ++i) {
      n[i] = map.direction.get(i)[p];
      dn[i] = direction_dt.get(i)[p];
      for (size_t a = 0; a < 2; ++a)
        dn[i] += flow[a][p] * map.derivative[a].get(i)[p];
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
  auto xx = gsl_vector_view_array(result.data(), 3);
  if (gsl_linalg_cholesky_decomp(&aa.matrix) != 0 or
      gsl_linalg_cholesky_solve(&aa.matrix, &bb.vector, &xx.vector) != 0)
    ERROR("Geometric rotation connection solve failed");
  return result;
}

void GeometricRelaxationSample::pup(PUP::er& p) {
  p | time;
  p | clock_rate;
  p | direction;
  p | weights;
  p | flow_rotation;
  p | raw;
  p | filtered;
}
bool operator==(const GeometricRelaxationSample& a,
                const GeometricRelaxationSample& b) {
  return a.time == b.time and a.clock_rate == b.clock_rate and
         a.direction == b.direction and a.weights == b.weights and
         a.flow_rotation == b.flow_rotation and a.raw == b.raw and
         a.filtered == b.filtered;
}
void GeometricRelaxationHistory::pup(PUP::er& p) { p | samples; }
bool operator==(const GeometricRelaxationHistory& a,
                const GeometricRelaxationHistory& b) {
  return a.samples == b.samples;
}

TidalMoments relax_geometric_tidal_moments(
    const gsl::not_null<GeometricRelaxationHistory*> history,
    const TidalMoments& raw, const EigenSphereMap& map,
    const GeometricTimeData& temporal, const double tau, const bool model_time,
    const double time, const double step_start, const bool full_step) {
  if (not(std::isfinite(tau) and tau > 0. and std::isfinite(time) and
          std::isfinite(step_start) and time >= step_start and
          std::isfinite(temporal.clock_rate) and temporal.clock_rate > 0.))
    ERROR(
        "Geometric relaxation requires finite forward time, positive tau and "
        "clock");
  const size_t np = map.weights.size();
  auto& past = history->samples;
  const bool compatible =
      past.empty() or past.back().direction.get(0).size() == np;
  // Stages only read history. A full-step evaluation discards superseded
  // endpoints, so a retry/recomputed endpoint starts from the same anchor.
  if (full_step) {
    if (not compatible)
      past.clear();
    while (not past.empty() and past.back().time >= time)
      past.pop_back();
  }
  GeometricRelaxationSample current{};
  current.time = time;
  current.clock_rate = temporal.clock_rate;
  current.direction = map.direction;
  current.weights = map.weights;
  current.raw = raw;
  current.filtered = raw;
  current.flow_rotation =
      geometric_rotation_rate(map, TriadVector(np, 0.), temporal.flow);
  const GeometricRelaxationSample* anchor = nullptr;
  if (compatible)
    for (const auto& sample : past)
      if (sample.time < time and sample.time <= step_start)
        anchor = &sample;
  if (anchor != nullptr) {
    const double dt = time - anchor->time;
    const auto q = multiply(
        rotation_step(current.flow_rotation, .5 * dt),
        multiply(sphere_map_rotation(anchor->direction, map.direction,
                                     0.5 * (anchor->weights + map.weights)),
                 rotation_step(anchor->flow_rotation, .5 * dt)));
    const auto old = rotate(q, anchor->filtered);
    const auto old_raw = rotate(q, anchor->raw);
    const double x =
        dt *
        (model_time ? .5 * (anchor->clock_rate + temporal.clock_rate) : 1.) /
        tau;
    const double decay = exp(-x), one_minus_decay = -expm1(-x);
    // Exact integral for linearly interpolated forcing; series avoids
    // cancellation when dt/tau is small.
    const double new_weight =
        x < 1.e-4 ? x * (.5 + x * (-1. / 6. + x * (1. / 24. - x / 120.)))
                  : 1. - one_minus_decay / x;
    const double old_weight = one_minus_decay - new_weight;
    for (size_t a = 0; a < 5; ++a)
      current.filtered[a] =
          decay * old[a] + old_weight * old_raw[a] + new_weight * raw[a];
  }
  const auto result = current.filtered;
  if (full_step) {
    past.push_back(std::move(current));
    if (past.size() > 6)
      past.erase(past.begin());
  }
  return result;
}
}  // namespace gr::np
