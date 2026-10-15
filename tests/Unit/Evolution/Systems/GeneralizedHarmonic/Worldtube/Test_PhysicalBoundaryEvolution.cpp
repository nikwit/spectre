// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeEvolvedPhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Actions/EvolvePhysicalTrace.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryEvolution.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/Rk4Owren.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Tensor = tnsr::aa<DataVector, 3>;
using Phi = tnsr::iaa<DataVector, 3>;
using gh::worldtube::PhysicalTraceGeometry;

void test_empty_history_projection() {
  using Tag = gh::worldtube::Tags::PhysicalTrace;
  const Slab slab(0., 1.);
  const TimeStepId id(true, 0, slab.start());
  TimeSteppers::History<DataVector> bulk(4), surface(1);
  bulk.insert(id, DataVector(8, 1.), DataVector(8, 2.));
  DataVector state, dt_state;
  gh::worldtube::ProjectPhysicalTrace<Tag>::apply(
      make_not_null(&state), make_not_null(&dt_state), make_not_null(&surface),
      bulk, std::pair<Mesh<3>, Element<3>>{});
  CHECK(surface.integration_order() == 4);
  REQUIRE(surface.size() == 1);
  CHECK(surface.back().time_step_id == id);
  CHECK(surface.back().derivative.size() == 0);
  CHECK(state.size() == 0);
  CHECK_THROWS_WITH(
      gh::worldtube::ProjectPhysicalTrace<Tag>::check_empty(DataVector(2, 0.)),
      Catch::Matchers::ContainsSubstring(
          "AMR of an evolved physical boundary"));
}

// Exact vacuum Rosen wave: ds^2=-dt^2+dx^2+a(u)^2 dY^2+b(u)^2 dZ^2,
// u=t-x, a=cos(k u), b=cosh(k u). R_uu=-(a''/a+b''/b)=0.
// Rotate (Y,Z) by a constant angle to excite both stored polarizations.
// This is a nonlinear solution, so there is no linearization error floor.
// Rosen curvature convention: https://arxiv.org/abs/1004.2734, Sec. 3.
struct Wave {
  static constexpr size_t n = 3;
  Tensor g{n, 0.}, dg{n, 0.}, ddg{n, 0.}, pi{n, 0.}, target{n, 0.};
  Phi phi{n, 0.};
  tnsr::a<DataVector, 3> h{n, 0.};
  tnsr::ab<DataVector, 3> dh{n, 0.};
  PhysicalTraceGeometry geometry;

  static std::array<tnsr::I<DataVector, 3>, 2> seeds() {
    std::array<tnsr::I<DataVector, 3>, 2> result{tnsr::I<DataVector, 3>(n, 0.),
                                                 tnsr::I<DataVector, 3>(n, 0.)};
    get<1>(result[0]) = 1.;
    get<2>(result[1]) = 1.;
    return result;
  }
  static tnsr::i<DataVector, 3> normal() {
    tnsr::i<DataVector, 3> result(n, 0.);
    get<0>(result) = -1.;
    return result;
  }
  static Tensor flat() {
    Tensor result(n, 0.);
    get<0, 0>(result) = -1.;
    get<1, 1>(result) = get<2, 2>(result) = get<3, 3>(result) = 1.;
    return result;
  }
  explicit Wave(const double time)
      : geometry(flat(), Tensor(n, 0.), normal(), seeds()) {
    get<0, 0>(g) = -1.;
    get<1, 1>(g) = 1.;
    for (size_t p = 0; p < n; ++p) {
      const double k = 0.5 + 0.1 * static_cast<double>(p);
      const double a = cos(k * time), b = cosh(k * time);
      const double da = -k * sin(k * time), db = k * sinh(k * time);
      const double dda = -k * k * a, ddb = k * k * b;
      const double angle = 0.17 + 0.21 * static_cast<double>(p);
      const std::array<double, 2> e{cos(angle), sin(angle)};
      const std::array<double, 2> f{-sin(angle), cos(angle)};
      for (size_t i = 0; i < 2; ++i) {
        for (size_t j = i; j < 2; ++j) {
          g.get(i + 2, j + 2)[p] = a * a * e[i] * e[j] + b * b * f[i] * f[j];
          dg.get(i + 2, j + 2)[p] =
              2. * a * da * e[i] * e[j] + 2. * b * db * f[i] * f[j];
          ddg.get(i + 2, j + 2)[p] = 2. * (da * da + a * dda) * e[i] * e[j] +
                                     2. * (db * db + b * ddb) * f[i] * f[j];
          target.get(i + 2, j + 2)[p] =
              -4. * a * dda * e[i] * e[j] - 4. * b * ddb * f[i] * f[j];
        }
      }
      // H=-Gamma is a prescribed nonzero gauge for these coordinates.
      const double P = da / a + db / b;
      const double dP = dda / a - square(da / a) + ddb / b - square(db / b);
      get<0>(h)[p] = P;
      get<1>(h)[p] = -P;
      get<0, 0>(dh)[p] = dP;
      get<0, 1>(dh)[p] = get<1, 0>(dh)[p] = -dP;
      get<1, 1>(dh)[p] = dP;
    }
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        pi.get(a, b) = -dg.get(a, b);
        phi.get(0, a, b) = -dg.get(a, b);
      }
    }
    geometry = PhysicalTraceGeometry(g, dg, normal(), seeds());
  }
  DataVector amplitudes() const {
    return gh::worldtube::physical_trace_amplitudes(pi, phi, geometry);
  }
  DataVector rhs(const DataVector& state, const bool full_jet = false) const {
    Tensor ghost_pi(n);
    Phi ghost_phi(n), dmetric(n, 0.), dpi(n, 0.);
    tnsr::ijaa<DataVector, 3> dphi(n, 0.);
    const gh::worldtube::PhysicalScreen screen(g, geometry.inverse_metric,
                                               geometry.time_normal,
                                               geometry.normal_vector);
    gh::worldtube::physical_ghost(
        make_not_null(&ghost_pi), make_not_null(&ghost_phi), pi, phi,
        gh::worldtube::physical_trace_tensor(state, geometry),
        geometry.normal_covector, geometry.normal_vector, screen);
    if (full_jet) {
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = a; b < 4; ++b) {
          dmetric.get(0, a, b) = -dg.get(a, b);
          dpi.get(0, a, b) = ddg.get(a, b);
          dphi.get(0, 0, a, b) = ddg.get(a, b);
        }
      }
    }
    return gh::worldtube::physical_boundary_rhs(
        g, ghost_pi, ghost_phi, dmetric, dpi, dphi, Scalar<DataVector>(n, 0.7),
        Scalar<DataVector>(n, 1.3), h, dh, target, geometry);
  }
};

void test_geometry_and_rhs() {
  for (double t : {0.1, 0.4, 0.9}) {
    const Wave w(t), plus(t + 1.e-5), minus(t - 1.e-5);
    auto tolerance = Approx::custom().epsilon(3.e-9).scale(1.);
    for (size_t p = 0; p < 2; ++p) {
      for (size_t c = 0; c < w.geometry.polarization[p].size(); ++c) {
        CHECK_ITERABLE_CUSTOM_APPROX(w.geometry.dt_polarization[p][c],
                                     (plus.geometry.polarization[p][c] -
                                      minus.geometry.polarization[p][c]) /
                                         2.e-5,
                                     tolerance);
      }
      for (size_t q = 0; q < 2; ++q) {
        DataVector norm(Wave::n, 0.);
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = 0; b < 4; ++b) {
            norm += w.geometry.polarization[p].get(a, b) *
                    w.geometry.lower_polarization[q].get(a, b);
          }
        }
        CHECK_ITERABLE_APPROX(norm, DataVector(Wave::n, p == q ? 1. : 0.));
      }
    }
    const DataVector state = w.amplitudes();
    CHECK(min(abs(state)) > 1.e-3);
    CHECK_ITERABLE_CUSTOM_APPROX(
        w.rhs(state), (plus.amplitudes() - minus.amplitudes()) / 2.e-5,
        tolerance);
    CHECK_ITERABLE_APPROX(w.rhs(state), w.rhs(state, true));
    // This is also checked away from the analytic trace. Normal derivatives
    // cannot enter the surface equation even when reduction constraints differ.
    CHECK_ITERABLE_APPROX(w.rhs(state + 0.023), w.rhs(state + 0.023, true));
  }
  // An oblique normal and a changing receiver metric exercise dot(s), which
  // vanishes in the normal-incidence Rosen-wave example.
  Tensor g = Wave::flat(), dg(Wave::n, 0.);
  get<1, 2>(dg) = 0.12;
  get<1, 1>(dg) = 0.31;
  get<2, 3>(dg) = -0.08;
  tnsr::i<DataVector, 3> normal(Wave::n, 0.);
  get<0>(normal) = 1.;
  get<1>(normal) = 0.3;
  auto seeds = Wave::seeds();
  get<0>(seeds[0]) = -0.3;
  const PhysicalTraceGeometry geometry(g, dg, normal, seeds);
  auto gp = g, gm = g;
  for (size_t c = 0; c < g.size(); ++c) {
    gp[c] += 1.e-5 * dg[c];
    gm[c] -= 1.e-5 * dg[c];
  }
  const PhysicalTraceGeometry plus(gp, dg, normal, seeds),
      minus(gm, dg, normal, seeds);
  auto tolerance = Approx::custom().epsilon(1.e-9).scale(1.);
  for (size_t i = 0; i < 3; ++i) {
    CHECK_ITERABLE_CUSTOM_APPROX(
        geometry.dt_normal_vector.get(i),
        (plus.normal_vector.get(i) - minus.normal_vector.get(i)) / 2.e-5,
        tolerance);
  }
}

// Exercise both motion of the coordinate surface and receiver metric changes.
void test_moving_geometry() {
  const ylm::Spherepack sphere(6, 6);
  const auto angles = sphere.theta_phi_points();
  const size_t n = sphere.physical_size();
  const auto evaluate = [&](const double t) {
    tnsr::I<DataVector, 3> x(n), velocity(n);
    const DataVector xx = sin(angles[0]) * cos(angles[1]);
    const DataVector yy = sin(angles[0]) * sin(angles[1]);
    get<0>(x) = 2. + .13 * t + cos(t) * xx - sin(t) * yy;
    get<1>(x) = -1. - .07 * t + sin(t) * xx + cos(t) * yy;
    get<2>(x) = .2 * t + (1. + .1 * t) * cos(angles[0]);
    get<0>(velocity) = .13 - sin(t) * xx - cos(t) * yy;
    get<1>(velocity) = -.07 + cos(t) * xx - sin(t) * yy;
    get<2>(velocity) = .2 + .1 * cos(angles[0]);
    Tensor g(n, 0.), dtg(n, 0.);
    get<0, 0>(g) = -1.;
    get<1, 1>(g) = 1. + .3 * t;
    get<2, 2>(g) = 1. + .2 * t;
    get<3, 3>(g) = 1. - .1 * t;
    get<1, 2>(g) = .08 * t;
    get<1, 1>(dtg) = .3;
    get<2, 2>(dtg) = .2;
    get<3, 3>(dtg) = -.1;
    get<1, 2>(dtg) = .08;
    return gh::worldtube::mapped_sphere_trace_geometry(g, dtg, x, velocity, 6);
  };
  const double t = .23, h = 1.e-5;
  const auto actual = evaluate(t), plus = evaluate(t + h),
             minus = evaluate(t - h);
  const auto tolerance = Approx::custom().epsilon(2.e-8).scale(1.);
  for (size_t i = 0; i < 3; ++i) {
    CHECK_ITERABLE_CUSTOM_APPROX(
        actual.dt_normal_vector.get(i),
        (plus.normal_vector.get(i) - minus.normal_vector.get(i)) / (2. * h),
        tolerance);
  }
  for (size_t p = 0; p < 2; ++p) {
    for (size_t c = 0; c < actual.polarization[p].size(); ++c) {
      CHECK_ITERABLE_CUSTOM_APPROX(
          actual.dt_polarization[p][c],
          (plus.polarization[p][c] - minus.polarization[p][c]) / (2. * h),
          tolerance);
    }
  }
}

void test_ale_and_damped_harmonic() {
  const Wave w(.4);
  const gh::gauges::DampedHarmonic gauge(17., {1., 0., 1.}, {2, 2, 2});
  const size_t n = Wave::n;
  tnsr::I<DataVector, 3> x(n, .3), velocity(n);
  get<0>(velocity) = .11;
  get<1>(velocity) = .07;
  get<2>(velocity) = -.04;
  Tensor pi(n);
  Phi phi(n);
  const gh::worldtube::PhysicalScreen screen(w.g, w.geometry.inverse_metric,
                                             w.geometry.time_normal,
                                             w.geometry.normal_vector);
  gh::worldtube::physical_ghost(
      make_not_null(&pi), make_not_null(&phi), w.pi, w.phi,
      gh::worldtube::physical_trace_tensor(w.amplitudes() + .023, w.geometry),
      w.geometry.normal_covector, w.geometry.normal_vector, screen);
  // Normal extensions must cancel even off the exact trace, with nonzero
  // mesh advection and native damped-harmonic derivatives on the ghost state.
  for (const bool dh : {false, true}) {
    Phi dg(n, 0.), dpi(n, 0.);
    tnsr::ijaa<DataVector, 3> dphi(n, 0.);
    const auto rhs = [&]() {
      return gh::worldtube::physical_boundary_rhs(
          w.g, pi, phi, dg, dpi, dphi, Scalar<DataVector>(n, .7),
          Scalar<DataVector>(n, 1.3), w.h, w.dh, w.target, w.geometry, velocity,
          dh ? &gauge : nullptr, &x);
    };
    const auto tangential = rhs();
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        dg.get(0, a, b) = .017 * static_cast<double>(a + b + 1);
        dpi.get(0, a, b) = -.023 * static_cast<double>(a + 2 * b + 1);
        for (size_t j = 0; j < 3; ++j) {
          dphi.get(0, j, a, b) = .031 * static_cast<double>(a + b + j + 1);
        }
      }
    }
    CHECK_ITERABLE_APPROX(tangential, rhs());
  }
}

double evolve_wave(const size_t steps, const bool restart_every_stage) {
  const TimeSteppers::Rk4Owren stepper{};
  const Slab slab(0.2, 0.9);
  auto id = TimeStepId(true, 0, slab.start());
  const auto dt = slab.duration() / static_cast<int64_t>(steps);
  TimeSteppers::History<DataVector> history(4);
  DataVector state = Wave(0.2).amplitudes();
  for (size_t step = 0; step < steps; ++step) {
    const auto saved_state = state;
    const auto saved_history = history;
    const auto saved_id = id;
    // First perform and discard a trial step. Restore all evolved state and
    // history, including any FSAL derivative retained by the native stepper.
    for (size_t attempt = 0; attempt < 2; ++attempt) {
      if (attempt == 1) {
        state = saved_state;
        history = saved_history;
        id = saved_id;
      }
      do {
        history.insert(id, state, Wave(id.substep_time()).rhs(state));
        if (restart_every_stage) {
          state = serialize_and_deserialize(state);
          history = serialize_and_deserialize(history);
          id = serialize_and_deserialize(id);
        }
        stepper.update_u(make_not_null(&state), history, dt);
        const auto updated = state;
        // Re-evaluating an update must not integrate the boundary twice.
        stepper.update_u(make_not_null(&state), history, dt);
        CHECK_ITERABLE_APPROX(state, updated);
        stepper.clean_history(make_not_null(&history));
        id = stepper.next_time_id(id, dt);
      } while (id.substep() != 0);
    }
  }
  return max(abs(state - Wave(0.9).amplitudes()));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Gh.PhysicalBoundaryEvolution",
                  "[Unit][Evolution][PhysicalTrace]") {
  register_classes_with_charm(
      tmpl::list<gh::Solutions::WrappedGr<gr::Solutions::Minkowski<3>>>{});
  const auto boundary = TestHelpers::test_creation<
      gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost>(
      "GaugeSector: SommerfeldAbsorbing\nAmplitude: 0.001\nFrequency: 0.7");
  const auto restored = serialize_and_deserialize(boundary);
  CHECK(restored.amplitude() == boundary.amplitude());
  CHECK(restored.frequency() == boundary.frequency());
  CHECK(boundary.get_clone() != nullptr);
  test_empty_history_projection();
  test_geometry_and_rhs();
  test_moving_geometry();
  test_ale_and_damped_harmonic();
  double previous = 0.;
  for (const size_t steps : std::array<size_t, 4>{4, 8, 16, 32}) {
    const double error = evolve_wave(steps, false);
    const double restarted_error = evolve_wave(steps, true);
    CAPTURE(steps, error, restarted_error);
    CHECK(error < 1.e-4);
    CHECK(error == approx(restarted_error));
    if (previous > 0.) {
      CHECK(previous / error > 13.);
    }
    previous = error;
  }
}
