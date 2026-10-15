// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/ThirdOrderTide.hpp"
#include "ThirdOrderReference.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace gr::np {
namespace {
GeometricFrame reference_frame() {
  namespace ref = third_order_reference;
  const size_t np = ref::weights.size();
  GeometricFrame frame{};
  frame.map.direction = TriadVector(np, 0.);
  for (auto& d : frame.map.derivative)
    d = TriadVector(np, 0.);
  frame.map.weights = DataVector(np);
  frame.boost = DataVector(np);
  for (auto& m : frame.dyad)
    m = ComplexDataVector(np, 0.);
  for (size_t p = 0; p < np; ++p) {
    frame.map.weights[p] = ref::weights[p];
    frame.boost[p] = ref::boost[p];
    for (size_t i = 0; i < 3; ++i) {
      frame.map.direction.get(i)[p] = ref::labels[3 * p + i];
      frame.dyad[i][p] = {ref::dyad[6 * p + 2 * i],
                          ref::dyad[6 * p + 2 * i + 1]};
    }
    const double x = ref::labels[3 * p], y = ref::labels[3 * p + 1],
                 z = ref::labels[3 * p + 2], st = hypot(x, y);
    const std::array<double, 3> et{{z * x / st, z * y / st, -st}},
        ep{{-y / st, x / st, 0.}};
    for (size_t i = 0; i < 3; ++i) {
      frame.map.derivative[0].get(i)[p] = et[i];
      frame.map.derivative[1].get(i)[p] = ep[i];
    }
  }
  return frame;
}

void test_profiles_and_fit() {
  namespace ref = third_order_reference;
  auto frame = reference_frame();
  const size_t np = frame.map.weights.size();
  Scalar<DataVector> radius{DataVector(np)};
  DataVector tilt(np);
  for (size_t p = 0; p < np; ++p) {
    get(radius)[p] = ref::radius[p];
    tilt[p] = ref::tilt[p];
  }
  const auto columns = third_order_tide_columns(frame, radius, .7, tilt);
  const auto tol = Approx::custom().epsilon(3.e-12).scale(1.);
  for (size_t a = 0; a < 34; ++a)
    for (size_t slot = 0; slot < 2; ++slot)
      for (size_t p = 0; p < 7; ++p) {
        const size_t offset = 2 * ((a * 2 + slot) * 7 + p);
        CHECK(real(columns[a][slot][p]) == tol(ref::columns[offset]));
        CHECK(imag(columns[a][slot][p]) == tol(ref::columns[offset + 1]));
      }
  FrameRegistration registration{};
  registration.rotation.a_bar = Scalar<ComplexDataVector>(np, 0.);
  registration.rotation.b = Scalar<ComplexDataVector>(np, 0.);
  registration.pulled_back = WeylScalars(np, std::complex<double>{0., 0.});
  registration.pulled_back.get(2) = -.01;
  for (size_t p = 0; p < np; ++p)
    registration.pulled_back.get(4)[p] = {ref::psi[2 * (np + p)],
                                          ref::psi[2 * (np + p) + 1]};
  DottedTidalMoments dots{};
  std::copy_n(ref::coefficients.begin() + 24, 10, dots.begin());
  const auto fit =
      fit_third_order_tide(registration, columns, frame.map.weights, dots);
  for (size_t a = 0; a < 24; ++a)
    CHECK(fit.moments[a] == tol(ref::coefficients[a]));
  double error = 0., norm = 0.;
  for (size_t p = 0; p < np; ++p) {
    const std::complex<double> truth{ref::psi[2 * p], ref::psi[2 * p + 1]};
    error += std::norm(get(fit.psi0_target)[p] - truth);
    norm += std::norm(truth);
  }
  CHECK(sqrt(error / norm) < 1.e-10);
  CHECK(fit.relative_residual < 1.e-11);
  // End-to-end causal prediction: the true quadrupoles vary linearly in
  // model time. Preliminary undotted fits contain a constant derivative
  // bias, which drops out of their causal derivative. The held-out incoming
  // data still come exclusively from the independent Python fixture.
  GeometricTideHistory history{};
  GeometricTimeData temporal{};
  temporal.clock_rate = 2.5;
  temporal.flow = {{DataVector(np, 0.), DataVector(np, 0.)}};
  for (const double time : {-.41, -.25, -.13, -.08, 0.}) {
    auto slice = registration;
    for (size_t p = 0; p < np; ++p) {
      for (size_t a = 0; a < 10; ++a) {
        slice.pulled_back.get(4)[p] +=
            time * temporal.clock_rate * dots[a] * columns[a][1][p];
      }
    }
    const auto raw =
        fit_third_order_tide(slice, columns, frame.map.weights, {});
    GeometricTideSample sample{};
    sample.time = time;
    sample.direction = frame.map.direction;
    std::copy_n(raw.moments.begin(), 10, sample.quadrupole.begin());
    const auto causal = causal_tidal_derivative(
        make_not_null(&history), sample, frame.map, temporal, time, true);
    if (time == 0.) {
      CHECK_ITERABLE_CUSTOM_APPROX(causal.dots, dots, tol);
      const auto prediction =
          fit_third_order_tide(slice, columns, frame.map.weights, causal.dots);
      CHECK_ITERABLE_CUSTOM_APPROX(get(prediction.psi0_target),
                                   get(fit.psi0_target), tol);
    }
  }
  // A change of all model axes must not alter the target. Refit the
  // outgoing data; do not rotate the answer with the implementation's basis.
  for (size_t p = 0; p < np; ++p) {
    const double x = frame.map.direction.get(0)[p],
                 z = frame.map.direction.get(2)[p];
    frame.map.direction.get(0)[p] = cos(.41) * x + sin(.41) * z;
    frame.map.direction.get(2)[p] = -sin(.41) * x + cos(.41) * z;
    const auto mx = frame.dyad[0][p], mz = frame.dyad[2][p];
    frame.dyad[0][p] = cos(.41) * mx + sin(.41) * mz;
    frame.dyad[2][p] = -sin(.41) * mx + cos(.41) * mz;
  }
  // Use undotted data for this independent octupole rotation check.
  for (size_t p = 0; p < np; ++p)
    for (size_t a = 0; a < 10; ++a)
      registration.pulled_back.get(4)[p] -= dots[a] * columns[a + 24][1][p];
  const auto original =
      fit_third_order_tide(registration, columns, frame.map.weights, {});
  const auto rotated = fit_third_order_tide(
      registration, third_order_tide_columns(frame, radius, .7, tilt),
      frame.map.weights, {});
  CHECK_ITERABLE_CUSTOM_APPROX(get(original.psi0_target),
                               get(rotated.psi0_target), tol);
  // Nonzero measured longitudinal slots survive incoming-slot replacement.
  registration.rotation.a_bar =
      Scalar<ComplexDataVector>(np, std::complex<double>{.02, -.03});
  registration.rotation.b =
      Scalar<ComplexDataVector>(np, std::complex<double>{-.01, .04});
  registration.pulled_back.get(1) = std::complex<double>{.007, .013};
  registration.pulled_back.get(3) = std::complex<double>{-.019, .002};
  auto replaced = registration.pulled_back;
  replaced.get(0) = get(original.psi0_target);
  const auto expected = push_forward(replaced, registration.rotation.a_bar,
                                     registration.rotation.b);
  const auto pushed =
      fit_third_order_tide(registration, columns, frame.map.weights, {});
  CHECK_ITERABLE_CUSTOM_APPROX(get(pushed.psi0_target), expected.get(0), tol);
}

void test_clock_and_potential() {
  const auto reference = reference_frame();
  const auto& labels = reference.map.direction;
  const auto& weights = reference.map.weights;
  const size_t np = weights.size();
  tnsr::ii<DataVector, 3, Frame::Inertial> metric(np, 0.);
  for (size_t i = 0; i < 3; ++i)
    metric.get(i, i) = 1.;
  TriadVector inward = labels;
  for (auto& c : inward)
    c *= -1.;
  const auto rotation = adapted_triad(metric, inward);
  FrameRegistration registration{};
  registration.rotation.a_bar = Scalar<ComplexDataVector>(np, 0.);
  registration.rotation.b = Scalar<ComplexDataVector>(np, 0.);
  registration.measured_radius = Scalar<DataVector>(np, 4.);
  const auto frame = geometric_frame(registration, Scalar<DataVector>(np, 0.),
                                     rotation, metric, labels, weights, 1., 4);
  TriadVector zero(np, 0.);
  const auto temporal = geometric_time_data(
      frame, labels, 4., registration.measured_radius, 1., rotation, metric,
      Scalar<DataVector>(np, sqrt(.5)), zero, zero);
  CHECK(temporal.clock_rate == approx(1.));
  CHECK_ITERABLE_APPROX(temporal.slice_tilt, DataVector(np, 0.));
  CHECK_ITERABLE_APPROX(temporal.flow[0], DataVector(np, 0.));
  CHECK_ITERABLE_APPROX(temporal.flow[1], DataVector(np, 0.));
  // Prescribed rotating grid: no-screen flow cancels its motion.
  auto velocity = zero;
  velocity.get(0) = -.2 * 4. * labels.get(1);
  velocity.get(1) = .2 * 4. * labels.get(0);
  const auto moving = geometric_time_data(
      frame, labels, 4., registration.measured_radius, 1., rotation, metric,
      Scalar<DataVector>(np, sqrt(.5)), zero, velocity);
  for (size_t i = 0; i < 3; ++i) {
    const DataVector total = velocity.get(i) / 4. +
                             moving.flow[0] * frame.map.derivative[0].get(i) +
                             moving.flow[1] * frame.map.derivative[1].get(i);
    CHECK_ITERABLE_APPROX(total, DataVector(np, 0.));
  }
  CHECK(moving.clock_rate == approx(1.));
  // An analytic linear-plus-quadratic model-time tilt and its gradient.
  std::array<DataVector, 2> gradient{{DataVector(np, 0.), DataVector(np, 0.)}};
  DataVector truth = .3 * labels.get(0) + .2 * labels.get(1) * labels.get(2);
  for (size_t a = 0; a < 2; ++a)
    gradient[a] = .3 * reference.map.derivative[a].get(0) +
                  .2 * (reference.map.derivative[a].get(1) * labels.get(2) +
                        labels.get(1) * reference.map.derivative[a].get(2));
  const auto potential =
      sphere_gradient_potential(labels, weights, gradient, 2);
  CHECK_ITERABLE_APPROX(potential.first, truth);
  CHECK(potential.second < 1.e-12);
}

void test_causal_history() {
  const auto frame = reference_frame();
  const size_t np = frame.map.weights.size();
  GeometricTimeData temporal{};
  temporal.clock_rate = 2.5;
  temporal.flow = {{DataVector(np, 0.), DataVector(np, 0.)}};
  GeometricTideHistory history{};
  const auto sample = [&frame](double t) {
    GeometricTideSample s{};
    s.time = t;
    s.direction = frame.map.direction;
    for (size_t a = 0; a < 10; ++a)
      s.quadrupole[a] = (.2 + .03 * a) * (.7 + .1 * t + .02 * t * t -
                                          .004 * pow(t, 3) + .002 * pow(t, 4));
    return s;
  };
  for (const double t : {0., .08, .13, .25, .41}) {
    const auto d = causal_tidal_derivative(make_not_null(&history), sample(t),
                                           frame.map, temporal, t, true);
    if (t == 0.) {
      CHECK(d.derivative_order == 0);
      CHECK(d.dots == DottedTidalMoments{});
    }
  }
  const auto current = sample(.41);
  const auto d = causal_tidal_derivative(make_not_null(&history), current,
                                         frame.map, temporal, .41, true);
  CHECK(d.derivative_order == 4);
  for (size_t a = 0; a < 10; ++a)
    CHECK(d.dots[a] ==
          approx((.2 + .03 * a) *
                 (.1 + .04 * .41 - .012 * pow(.41, 2) + .008 * pow(.41, 3)) /
                 2.5));
  const auto before = history;
  const auto repeated = causal_tidal_derivative(
      make_not_null(&history), current, frame.map, temporal, .41, true);
  CHECK(history == before);
  CHECK_ITERABLE_APPROX(repeated.dots, d.dots);
  history = serialize_and_deserialize(history);
  CHECK(history == before);
  // An intermediate stage is evaluated but cannot contaminate the history.
  const auto stage = causal_tidal_derivative(
      make_not_null(&history), sample(.47), frame.map, temporal, .41, false);
  CHECK(history == before);
  for (size_t a = 0; a < 10; ++a)
    CHECK(stage.dots[a] ==
          approx((.2 + .03 * a) *
                 (.1 + .04 * .47 - .012 * pow(.47, 2) + .008 * pow(.47, 3)) /
                 2.5));
  causal_tidal_derivative(make_not_null(&history), sample(.41), frame.map,
                          temporal, .41, false);
  CHECK(history == before);
  causal_tidal_derivative(make_not_null(&history), sample(.43), frame.map,
                          temporal, .41, false);
  CHECK(history == before);
  // A self-start/time rollback discards the previous attempt.
  const auto reset = causal_tidal_derivative(
      make_not_null(&history), sample(.1), frame.map, temporal, .1, true);
  CHECK(reset.derivative_order == 0);
  CHECK(history.samples.size() == 1);
  history.samples.back().direction = TriadVector(np + 1, 0.);
  const auto resized = causal_tidal_derivative(
      make_not_null(&history), sample(.2), frame.map, temporal, .2, true);
  CHECK(resized.derivative_order == 0);
}

// Rotating auxiliary axes and physical change are deliberately independent.
// The transported derivative should recover R dH R^T, not d(R H R^T)/dT.
double rotation_error(const double dt, const bool rotating_grid) {
  const auto base = reference_frame();
  const size_t np = base.map.weights.size();
  GeometricTimeData temporal{};
  temporal.clock_rate = 1.;
  temporal.flow = {{DataVector(np, 0.), DataVector(np, 0.)}};
  GeometricTideHistory history{};
  CausalTidalDerivative last{};
  for (int step = -4; step <= 0; ++step) {
    const double t = dt * step, c = cos(.7 * t), s = sin(.7 * t);
    GeometricTideSample current{};
    current.time = t;
    auto map = base.map;
    for (size_t p = 0; p < np; ++p) {
      const double x = base.map.direction.get(0)[p],
                   y = base.map.direction.get(1)[p];
      map.direction.get(0)[p] = c * x - s * y;
      map.direction.get(1)[p] = s * x + c * y;
      for (size_t a = 0; a < 2; ++a) {
        const double dx = base.map.derivative[a].get(0)[p],
                     dy = base.map.derivative[a].get(1)[p];
        map.derivative[a].get(0)[p] = c * dx - s * dy;
        map.derivative[a].get(1)[p] = s * dx + c * dy;
      }
      if (rotating_grid) {
        const std::array<double, 3> counter{
            {.7 * map.direction.get(1)[p], -.7 * map.direction.get(0)[p], 0.}};
        for (size_t a = 0; a < 2; ++a) {
          temporal.flow[a][p] = 0.;
          for (size_t i = 0; i < 3; ++i)
            temporal.flow[a][p] += counter[i] * map.derivative[a].get(i)[p];
        }
      }
    }
    current.direction = map.direction;
    // H=diag(1,-1,0), with dH=diag(.2,-.2,0).
    current.quadrupole[0] = rotating_grid ? 1. : (1. + .2 * t) * cos(1.4 * t);
    current.quadrupole[1] = -current.quadrupole[0];
    current.quadrupole[2] = rotating_grid ? 0. : (1. + .2 * t) * sin(1.4 * t);
    last = causal_tidal_derivative(make_not_null(&history), current, map,
                                   temporal, t, true);
  }
  double error = 0.;
  for (size_t a = 0; a < 10; ++a) {
    const double truth =
        rotating_grid ? 0. : (a == 0 ? .2 : (a == 1 ? -.2 : 0.));
    error += pow(last.dots[a] - truth, 2);
  }
  return sqrt(error);
}
}  // namespace
}  // namespace gr::np

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.GeneralRelativity.NP.ThirdOrderTide",
                  "[Unit][PointwiseFunctions]") {
  gr::np::test_profiles_and_fit();
  gr::np::test_clock_and_potential();
  gr::np::test_causal_history();
  for (bool moving : {false, true}) {
    const double coarse = gr::np::rotation_error(.04, moving);
    const double fine = gr::np::rotation_error(.02, moving);
    CHECK(fine < coarse / 12.);
    CHECK(fine < 1.e-6);
  }
}
