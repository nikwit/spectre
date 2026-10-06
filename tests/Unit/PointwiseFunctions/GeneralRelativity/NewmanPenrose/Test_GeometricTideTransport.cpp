// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTideTransport.hpp"
#include "ThirdOrderReference.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace {
using namespace gr::np;
using C = std::complex<double>;
using M = SpatialRotation;
M rz(const double t) {
  return {{{{cos(t), -sin(t), 0.}}, {{sin(t), cos(t), 0.}}, {{0., 0., 1.}}}};
}
M rx(const double t) {
  return {{{{1., 0., 0.}}, {{0., cos(t), -sin(t)}}, {{0., sin(t), cos(t)}}}};
}
M product(const M& a, const M& b) {
  M c{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        c[i][j] += a[i][k] * b[k][j];
  return c;
}
TidalMoments rotated(const TidalMoments& a, const M& r) {
  const std::array<std::array<C, 3>, 3> h{{{{a[0], a[2], a[3]}},
                                           {{a[2], a[1], a[4]}},
                                           {{a[3], a[4], -a[0] - a[1]}}}};
  std::array<std::array<C, 3>, 3> out{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          out[i][j] += r[i][k] * r[j][l] * h[k][l];
  return {{out[0][0], out[1][1], out[0][1], out[0][2], out[1][2]}};
}
std::pair<EigenSphereMap, GeometricTimeData> geometry(const double time,
                                                      const bool moving_axes,
                                                      const bool moving_grid,
                                                      const double clock = 1.) {
  namespace ref = third_order_reference;
  const size_t n = ref::weights.size();
  EigenSphereMap map{};
  map.direction = TriadVector(n, 0.);
  for (auto& d : map.derivative)
    d = TriadVector(n, 0.);
  map.weights = DataVector(n);
  GeometricTimeData temporal{};
  temporal.clock_rate = clock;
  temporal.flow = {{DataVector(n, 0.), DataVector(n, 0.)}};
  const auto r = product(rx(moving_axes ? .9 * time : 0.),
                         rz(moving_grid ? .7 * time : 0.));
  for (size_t p = 0; p < n; ++p) {
    map.weights[p] = ref::weights[p];
    const double x = ref::labels[3 * p], y = ref::labels[3 * p + 1],
                 z = ref::labels[3 * p + 2];
    const double st = hypot(x, y);
    const std::array<std::array<double, 3>, 3> basis{
        {{{x, y, z}},
         {{z * x / st, z * y / st, -st}},
         {{-y / st, x / st, 0.}}}};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j) {
        map.direction.get(i)[p] += r[i][j] * basis[0][j];
        for (size_t a = 0; a < 2; ++a)
          map.derivative[a].get(i)[p] += r[i][j] * basis[a + 1][j];
      }
    if (moving_grid)
      temporal.flow[1][p] = -.7 * st;
  }
  return {map, temporal};
}
TidalMoments tide(const double t) {
  return {{C(1. + .3 * t, .2 - .1 * t), C(-.3 + .2 * t, -.4),
           C(.1, -.2 + .1 * t), C(.2, .1), C(-.1, .3)}};
}
void test_constant_and_linear() {
  for (bool axes : {false, true})
    for (bool grid : {false, true})
      for (bool model : {false, true}) {
        CAPTURE(axes, grid, model);
        GeometricRelaxationHistory h{};
        const double tau = .8;
        // Nonuniform steps; a linear physical tide has a closed-form response.
        for (double t : {0., .02, .09, .15, .31, .46, .61}) {
          const double clock = 2.3 + .4 * t;
          const auto [map, temporal] = geometry(t, axes, grid, clock);
          const double physical_t = model ? 2.3 * t + .2 * t * t : t;
          const auto r = rx(axes ? .9 * t : 0.);
          const auto raw = rotated(tide(physical_t), r);
          const auto out = relax_geometric_tidal_moments(
              make_not_null(&h), raw, map, temporal, tau, model, t, t, true);
          const double filtered_t =
              physical_t - tau * (1. - exp(-physical_t / tau));
          const auto expected = rotated(tide(filtered_t), r);
          for (size_t a = 0; a < 5; ++a) {
            CHECK(real(out[a]) == approx(real(expected[a])));
            CHECK(imag(out[a]) == approx(imag(expected[a])));
          }
          CHECK(serialize_and_deserialize(h) == h);
        }
      }
  // Constant physical tide under independent, noncommuting axis/grid rotations.
  GeometricRelaxationHistory h{};
  for (double t : {0., .07, .18, .4, .7}) {
    const auto [map, temporal] = geometry(t, true, true);
    const auto raw = rotated(tide(0.), rx(.9 * t));
    const auto out = relax_geometric_tidal_moments(
        make_not_null(&h), raw, map, temporal, 1., false, t, t, true);
    for (size_t a = 0; a < 5; ++a)
      CHECK(abs(out[a] - raw[a]) < 2.e-14);
  }
}
void test_lifecycle() {
  GeometricRelaxationHistory h{};
  const auto evaluate = [&h](double t, double start, bool full) {
    const auto [map, temporal] = geometry(t, true, true);
    return relax_geometric_tidal_moments(make_not_null(&h),
                                         rotated(tide(t), rx(.9 * t)), map,
                                         temporal, .4, false, t, start, full);
  };
  evaluate(0., 0., true);
  evaluate(.1, .1, true);
  const auto committed = h;
  const auto stage = evaluate(.2, .1, false);
  CHECK(h == committed);
  evaluate(.3, .1,
           false);  // Later stage evaluated first must not contaminate .2.
  CHECK(evaluate(.2, .1, false) == stage);
  CHECK(h == committed);
  const auto endpoint = evaluate(.2, .2, true);
  CHECK(endpoint == stage);
  const auto saved = h;
  CHECK(evaluate(.2, .2, true) == endpoint);
  CHECK(h == saved);
  h = serialize_and_deserialize(h);
  const auto restarted = evaluate(.25, .2, false);
  h = saved;
  CHECK(evaluate(.25, .2, false) == restarted);
  evaluate(.3, .3, true);
  // Retry from .1: retained .1 anchor is restored; future .2/.3 are discarded.
  CHECK(evaluate(.2, .2, true) == endpoint);
  CHECK(h == saved);
  // Recomputed endpoint uses the revised raw forcing, not old endpoint state.
  auto [map, temporal] = geometry(.2, true, true);
  const auto revised = rotated(tide(.5), rx(.18));
  const auto revised_out = relax_geometric_tidal_moments(
      make_not_null(&h), revised, map, temporal, .4, false, .2, .2, true);
  auto fresh = committed;
  CHECK(revised_out == relax_geometric_tidal_moments(make_not_null(&fresh),
                                                     revised, map, temporal, .4,
                                                     false, .2, .2, true));
  // Angular resolution change resets from raw without touching history in
  // stages.
  auto changed = map;
  changed.weights = DataVector(3, 1.);
  changed.direction = TriadVector(size_t{3}, 0.);
  for (size_t i = 0; i < 3; ++i)
    changed.direction.get(i)[i] = 1.;
  for (auto& d : changed.derivative)
    d = TriadVector(size_t{3}, 0.);
  temporal.flow = {{DataVector(3, 0.), DataVector(3, 0.)}};
  const auto before = h;
  CHECK(relax_geometric_tidal_moments(make_not_null(&h), revised, changed,
                                      temporal, .4, false, .25, .2,
                                      false) == revised);
  CHECK(h == before);
  CHECK(relax_geometric_tidal_moments(make_not_null(&h), revised, changed,
                                      temporal, .4, false, .3, .3,
                                      true) == revised);
  CHECK(h.samples.size() == 1);
  // Startup/self-start before the retained history is reinitialized explicitly.
  CHECK(relax_geometric_tidal_moments(make_not_null(&h), revised, changed,
                                      temporal, .4, false, 0., 0.,
                                      true) == revised);
  CHECK(h.samples.size() == 1);
  for (size_t n = 1; n <= 8; ++n) {
    const double t = .1 * n;
    relax_geometric_tidal_moments(make_not_null(&h), revised, changed, temporal,
                                  .4, false, t, t, true);
  }
  CHECK(h.samples.size() == 6);
  CHECK(h.samples.front().time > 0.);
  CHECK(relax_geometric_tidal_moments(make_not_null(&h), revised, changed,
                                      temporal, .4, false, 0., 0.,
                                      true) == revised);
  CHECK(h.samples.size() == 1);
  CHECK_THROWS_WITH(
      (relax_geometric_tidal_moments(make_not_null(&h), revised, changed,
                                     temporal, 0., false, .1, .1, true)),
      Catch::Matchers::ContainsSubstring("positive tau"));
}
void test_held_forcing() {
  const auto [map, temporal] = geometry(0., false, false);
  for (const double x : {1.e-9, .01, 1., 10., 1.e4}) {
    GeometricRelaxationHistory h{};
    const auto raw = tide(0.);
    relax_geometric_tidal_moments(make_not_null(&h), raw, map, temporal, 1.,
                                  false, 0., 0., true);
    const auto initial = tide(1.);
    h.samples.back().filtered = initial;
    const auto out = relax_geometric_tidal_moments(
        make_not_null(&h), raw, map, temporal, 1., false, x, x, true);
    for (size_t a = 0; a < 5; ++a) {
      CHECK(abs(out[a] - (raw[a] + exp(-x) * (initial[a] - raw[a]))) < 2.e-14);
    }
  }
}
double sinusoid_error(const double dt) {
  // A rotating physical electric quadrupole (orbital frequency omega/2),
  // independent of auxiliary-axis and grid motion.
  GeometricRelaxationHistory h{};
  constexpr double omega = 1.7, tau = .6;
  double error = 0.;
  for (size_t n = 0; n <= static_cast<size_t>(1. / dt + .5); ++n) {
    const double t = n * dt;
    const auto [map, temporal] = geometry(t, true, true);
    const C signal = exp(C(0., omega * t));
    const TidalMoments raw0{
        {real(signal), -real(signal), imag(signal), 0., 0.}};
    const auto out = relax_geometric_tidal_moments(
        make_not_null(&h), rotated(raw0, rx(.9 * t)), map, temporal, tau, false,
        t, t, true);
    const C truth =
        (signal + C(0., omega * tau) * exp(-t / tau)) / (C(1., omega * tau));
    const auto expected =
        rotated(TidalMoments{{real(truth), -real(truth), imag(truth), 0., 0.}},
                rx(.9 * t));
    for (size_t a = 0; a < 5; ++a)
      error = std::max(error, abs(out[a] - expected[a]));
  }
  return error;
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.GeneralRelativity.NP.GeometricTideTransport",
    "[Unit][PointwiseFunctions]") {
  test_constant_and_linear();
  test_lifecycle();
  test_held_forcing();
  const double coarse = sinusoid_error(.02), fine = sinusoid_error(.01);
  CHECK(coarse / fine > 3.8);
  CHECK(fine < 3.e-5);
}
