// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SelfForce.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SelfForceDerivatives.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Framework/TestHelpers.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeodesicAcceleration.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {
constexpr double charge = 0.3;
constexpr double rest_mass = 1.2;
constexpr double turn_on_time = 2.;
constexpr double turn_on_interval = 0.8;

// An exact flat-space wave, quadratic in t - n.x. All nonzero spatial
// derivatives fit in the reconstruction, so its closure is exact.
RegularFieldDerivatives plane_wave(const tnsr::I<double, 3>& x,
                                   const double time) {
  const std::array<double, 3> n{{0.6, 0.8, 0.}};
  const double phase = time - 0.6 * (get<0>(x) - 8.) - 0.8 * (get<1>(x) - 6.);
  const double slope = 0.02 + 0.006 * phase;
  RegularFieldDerivatives field{};
  field.psi = 0.1 + 0.02 * phase + 0.003 * phase * phase;
  field.dt_psi = slope;
  for (size_t i = 0; i < 3; ++i) {
    field.d_psi.get(i) = -n.at(i) * slope;
    field.d_dt_psi.get(i) = -0.006 * n.at(i);
    for (size_t j = i; j < 3; ++j) {
      field.d2_psi.get(i, j) = 0.006 * n.at(i) * n.at(j);
      field.d2_dt_psi.get(i, j) = 0.;
    }
  }
  return field;
}

#ifdef SPECTRE_AUTODIFF
using State = std::array<tnsr::I<double, 3>, 2>;
State rhs(const State& state, const double time, const bool evolve_mass) {
  const auto field = plane_wave(state[0], time);
  const double mass = rest_mass - (evolve_mass ? charge * field.psi : 0.);
  const double w = time > turn_on_time
                       ? turn_on_function(time - turn_on_time, turn_on_interval)
                       : 0.;
  double inverse_u0_squared = 1.;
  for (size_t i = 0; i < 3; ++i) {
    inverse_u0_squared -= state[1].get(i) * state[1].get(i);
  }
  State result{{state[1], tnsr::I<double, 3>{0.}}};
  for (size_t i = 0; i < 3; ++i) {
    result[1].get(i) = w * charge / mass * inverse_u0_squared *
                       (field.d_psi.get(i) + state[1].get(i) * field.dt_psi);
  }
  return result;
}
State add(const State& a, const State& b, const double scale) {
  auto result = a;
  for (size_t v = 0; v < 2; ++v) {
    for (size_t i = 0; i < 3; ++i) {
      result.at(v).get(i) += scale * b.at(v).get(i);
    }
  }
  return result;
}
State rk4(const State& state, const double time, const double h,
          const bool evolve_mass) {
  const auto k1 = rhs(state, time, evolve_mass);
  const auto k2 = rhs(add(state, k1, h / 2.), time + h / 2., evolve_mass);
  const auto k3 = rhs(add(state, k2, h / 2.), time + h / 2., evolve_mass);
  const auto k4 = rhs(add(state, k3, h), time + h, evolve_mass);
  return add(add(add(add(state, k1, h / 6.), k2, h / 3.), k3, h / 3.), k4,
             h / 6.);
}

void test_time_derivatives() {
  const State state{{tnsr::I<double, 3>{{8., 6., 0.}},
                     tnsr::I<double, 3>{{-0.04, 0.08, 0.}}}};
  const double h = 2.e-4;
  const auto local_approx = Approx::custom().epsilon(2.e-10).scale(1.);
  for (const bool evolve_mass : {false, true}) {
    for (const double time : {1.5, 2., 2.7, 5.}) {
      CAPTURE(evolve_mass, time);
      const auto evaluate_at = [&](const double dt) {
        const auto shifted = rk4(state, time, dt, evolve_mass);
        return second_order_self_force_derivatives(
            plane_wave(shifted[0], time + dt), shifted[0], shifted[1], charge,
            rest_mass, evolve_mass, 0., time + dt, turn_on_time,
            turn_on_interval);
      };
      const auto result = evaluate_at(0.);
      const auto minus2 = evaluate_at(-2. * h);
      const auto minus1 = evaluate_at(-h);
      const auto plus1 = evaluate_at(h);
      const auto plus2 = evaluate_at(2. * h);
      CHECK_ITERABLE_APPROX(result.coordinate_acceleration,
                            rhs(state, time, evolve_mass)[1]);
      for (size_t a = 0; a < 4; ++a) {
        CHECK(result.dt_f.get(a) ==
              local_approx((minus2.f.get(a) - 8. * minus1.f.get(a) +
                            8. * plus1.f.get(a) - plus2.f.get(a)) /
                           (12. * h)));
        CHECK(result.dt_Du_f.get(a) ==
              local_approx((minus2.Du_f.get(a) - 8. * minus1.Du_f.get(a) +
                            8. * plus1.Du_f.get(a) - plus2.Du_f.get(a)) /
                           (12. * h)));
        // The turn-on is C^3; the stencil across its join has reduced order
        // for the derivative of D_u^2 f. Test that join separately below.
        if (time != turn_on_time) {
          CHECK(result.dt_DuDu_f.get(a) ==
                local_approx((minus2.DuDu_f.get(a) - 8. * minus1.DuDu_f.get(a) +
                              8. * plus1.DuDu_f.get(a) - plus2.DuDu_f.get(a)) /
                             (12. * h)));
        } else {
          CHECK(result.dt_DuDu_f.get(a) == 0.);
        }
      }
    }
  }
}

void test_curved_background() {
  const tnsr::I<double, 3> position{{8., 6., 0.}};
  const tnsr::I<double, 3> velocity{{-0.04, 0.08, 0.}};
  constexpr double mass = 1.4;
  const gr::Solutions::KerrSchild ks{mass, {{0., 0., 0.}}, {{0., 0., 0.}}};
  Tags::BackgroundQuantities<3>::type background{};
  Tags::BackgroundQuantitiesCompute<3>::function(make_not_null(&background),
                                                 {position, velocity}, ks);
  const auto& metric = get<gr::Tags::SpacetimeMetric<double, 3>>(background);
  const auto& inverse =
      get<gr::Tags::InverseSpacetimeMetric<double, 3>>(background);
  const auto& connection =
      get<gr::Tags::SpacetimeChristoffelSecondKind<double, 3>>(background);
  tnsr::A<double, 3> u{
      {1., get<0>(velocity), get<1>(velocity), get<2>(velocity)}};
  double norm = 0.;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      norm -= metric.get(a, b) * u.get(a) * u.get(b);
    }
  }
  for (auto& component : u) {
    component /= sqrt(norm);
  }
  const auto dot = [&metric](const auto& a, const auto& b) {
    double result = 0.;
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        result += metric.get(i, j) * a.get(i) * b.get(j);
      }
    }
    return result;
  };
  const auto local_approx = Approx::custom().epsilon(2.e-13).scale(1.);
  for (const double q : {0., charge}) {
    for (const double time : {1.5, 2.7, 5.}) {
      CAPTURE(q, time);
      const auto field = plane_wave(position, time);
      const auto result = second_order_self_force_derivatives(
          field, position, velocity, q, rest_mass, true, mass, time,
          turn_on_time, turn_on_interval);
      const double w =
          time > turn_on_time
              ? turn_on_function(time - turn_on_time, turn_on_interval)
              : 0.;
      const tnsr::a<double, 3> gradient{{field.dt_psi, get<0>(field.d_psi),
                                         get<1>(field.d_psi),
                                         get<2>(field.d_psi)}};
      auto expected_f = self_force_per_mass(gradient, u, q,
                                            rest_mass - q * field.psi, inverse);
      auto expected_acc = gr::geodesic_acceleration(velocity, connection);
      for (size_t a = 0; a < 4; ++a) {
        expected_f.get(a) *= w;
      }
      for (size_t i = 0; i < 3; ++i) {
        expected_acc.get(i) +=
            (expected_f.get(i + 1) - velocity.get(i) * expected_f.get(0)) /
            (get<0>(u) * get<0>(u));
      }
      CHECK_ITERABLE_APPROX(result.f, expected_f);
      CHECK_ITERABLE_APPROX(result.coordinate_acceleration, expected_acc);
      CHECK_ITERABLE_APPROX(
          result.Du_f,
          Du_self_force_per_mass(result.f, result.dt_f, u, connection));
      CHECK_ITERABLE_APPROX(
          result.DuDu_f,
          Du_self_force_per_mass(result.Du_f, result.dt_Du_f, u, connection));
      // Differentiate u.u=-1 along the accelerated worldline. These check
      // velocity normalization, connection terms, and the force's projection.
      CHECK(dot(u, result.f) == local_approx(0.));
      CHECK(dot(u, result.Du_f) == local_approx(-dot(result.f, result.f)));
      CHECK(dot(u, result.DuDu_f) ==
            local_approx(-3. * dot(result.f, result.Du_f)));
      if (q == 0. or w == 0.) {
        CHECK_ITERABLE_APPROX(result.dt_DuDu_f, (tnsr::A<double, 3>{0.}));
      }
    }
  }
}

// Independently evolve the local Cauchy data to first order in time using
// the wave equation and the established GR background implementation. Spatial
// derivatives of that RHS are finite-differenced only in this test. This
// exercises the curvature contributions to the highest force derivative.
void test_curved_time_derivative() {
  const tnsr::I<double, 3> position{{8., 6., 0.}};
  const tnsr::I<double, 3> velocity{{-0.04, 0.08, 0.}};
  constexpr double mass = 1.4;
  constexpr double time = 2.7;
  auto field = plane_wave(position, time);
  field.d2_psi.get(2, 2) = 0.004;
  field.d2_dt_psi.get(0, 0) = 0.003;
  field.d2_dt_psi.get(0, 1) = -0.002;
  field.d2_dt_psi.get(2, 2) = -0.001;
  const gr::Solutions::KerrSchild ks{mass, {{0., 0., 0.}}, {{0., 0., 0.}}};
  const auto spatial_field = [&](const tnsr::I<double, 3>& x) {
    auto shifted = field;
    for (size_t i = 0; i < 3; ++i) {
      const double dx_i = x.get(i) - position.get(i);
      shifted.psi += field.d_psi.get(i) * dx_i;
      shifted.dt_psi += field.d_dt_psi.get(i) * dx_i;
      for (size_t j = 0; j < 3; ++j) {
        const double dx_j = x.get(j) - position.get(j);
        shifted.psi += 0.5 * field.d2_psi.get(i, j) * dx_i * dx_j;
        shifted.dt_psi += 0.5 * field.d2_dt_psi.get(i, j) * dx_i * dx_j;
        shifted.d_psi.get(i) += field.d2_psi.get(i, j) * dx_j;
        shifted.d_dt_psi.get(i) += field.d2_dt_psi.get(i, j) * dx_j;
      }
    }
    return shifted;
  };
  const auto wave_rhs = [&](const tnsr::I<double, 3>& x) {
    const auto shifted = spatial_field(x);
    Tags::BackgroundQuantities<3>::type background{};
    Tags::BackgroundQuantitiesCompute<3>::function(make_not_null(&background),
                                                   {x, velocity}, ks);
    const auto& inv =
        get<gr::Tags::InverseSpacetimeMetric<double, 3>>(background);
    const auto& gamma =
        get<gr::Tags::TraceSpacetimeChristoffelSecondKind<double, 3>>(
            background);
    double result = gamma.get(0) * shifted.dt_psi;
    for (size_t i = 0; i < 3; ++i) {
      result += -2. * inv.get(0, i + 1) * shifted.d_dt_psi.get(i) +
                gamma.get(i + 1) * shifted.d_psi.get(i);
      for (size_t j = 0; j < 3; ++j) {
        result -= inv.get(i + 1, j + 1) * shifted.d2_psi.get(i, j);
      }
    }
    return result / inv.get(0, 0);
  };
  const auto result = second_order_self_force_derivatives(
      field, position, velocity, charge, rest_mass, true, mass, time,
      turn_on_time, turn_on_interval);
  const auto at_time = [&](const double dt) {
    auto x = position;
    auto v = velocity;
    for (size_t i = 0; i < 3; ++i) {
      x.get(i) += dt * velocity.get(i);
      v.get(i) += dt * result.coordinate_acceleration.get(i);
    }
    auto evolved = spatial_field(x);
    const auto initial = evolved;
    evolved.psi += dt * initial.dt_psi;
    evolved.dt_psi += dt * wave_rhs(x);
    constexpr double h = 0.01;
    const std::array<double, 4> offsets{{-2., -1., 1., 2.}};
    const std::array<double, 4> weights{{1., -8., 8., -1.}};
    for (size_t i = 0; i < 3; ++i) {
      double d_rhs = 0.;
      for (size_t a = 0; a < 4; ++a) {
        auto point = x;
        point.get(i) += offsets.at(a) * h;
        d_rhs += weights.at(a) * wave_rhs(point) / (12. * h);
      }
      evolved.d_psi.get(i) += dt * initial.d_dt_psi.get(i);
      evolved.d_dt_psi.get(i) += dt * d_rhs;
      for (size_t j = i; j < 3; ++j) {
        double d2_rhs = 0.;
        // Applying two fourth-order first-derivative stencils works for
        // diagonal as well as mixed derivatives.
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = 0; b < 4; ++b) {
            auto point = x;
            point.get(i) += offsets.at(a) * h;
            point.get(j) += offsets.at(b) * h;
            d2_rhs += weights.at(a) * weights.at(b) * wave_rhs(point) /
                      (144. * h * h);
          }
        }
        evolved.d2_psi.get(i, j) += dt * initial.d2_dt_psi.get(i, j);
        evolved.d2_dt_psi.get(i, j) += dt * d2_rhs;
      }
    }
    return second_order_self_force_derivatives(evolved, x, v, charge, rest_mass,
                                               true, mass, time + dt,
                                               turn_on_time, turn_on_interval);
  };
  constexpr double dt = 1.e-5;
  const auto plus = at_time(dt);
  const auto minus = at_time(-dt);
  const auto local_approx = Approx::custom().epsilon(3.e-9).scale(1.);
  for (size_t a = 0; a < 4; ++a) {
    CHECK(result.dt_f.get(a) ==
          local_approx((plus.f.get(a) - minus.f.get(a)) / (2. * dt)));
    CHECK(result.dt_Du_f.get(a) ==
          local_approx((plus.Du_f.get(a) - minus.Du_f.get(a)) / (2. * dt)));
    CHECK(result.dt_DuDu_f.get(a) ==
          local_approx((plus.DuDu_f.get(a) - minus.DuDu_f.get(a)) / (2. * dt)));
  }
}

#endif  // SPECTRE_AUTODIFF

void test_reconstruction() {
  const tnsr::I<double, 3> position{{8., 6., 0.}};
  const tnsr::I<double, 3> velocity{{-0.04, 0.08, 0.}};
  auto expected = plane_wave(position, 2.7);
  // Include a nonzero Hessian of the time-derivative field as well.
  expected.d2_dt_psi.get(0, 0) = 0.03;
  expected.d2_dt_psi.get(1, 2) = -0.02;
  for (const double radius : {0.1, 0.7}) {
    double trace = 0., dt_trace = 0., total_dt = expected.dt_psi;
    for (size_t i = 0; i < 3; ++i) {
      trace += expected.d2_psi.get(i, i);
      dt_trace += expected.d2_dt_psi.get(i, i);
      total_dt += velocity.get(i) * expected.d_psi.get(i);
    }
    tnsr::ii<double, 3> stf{}, dt_stf{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = i; j < 3; ++j) {
        stf.get(i, j) =
            expected.d2_psi.get(i, j) / 2. - (i == j ? trace / 6. : 0.);
        dt_stf.get(i, j) =
            expected.d2_dt_psi.get(i, j) / 2. - (i == j ? dt_trace / 6. : 0.);
      }
    }
    const auto result = reconstruct_regular_field(
        expected.psi + radius * radius * trace / 6.,
        expected.dt_psi + radius * radius * dt_trace / 6., expected.d_psi,
        expected.d_dt_psi, stf, dt_stf, expected.psi, total_dt, velocity,
        radius);
    CHECK(result.psi == approx(expected.psi));
    CHECK(result.dt_psi == approx(expected.dt_psi));
    CHECK_ITERABLE_APPROX(result.d_psi, expected.d_psi);
    CHECK_ITERABLE_APPROX(result.d_dt_psi, expected.d_dt_psi);
    CHECK_ITERABLE_APPROX(result.d2_psi, expected.d2_psi);
    CHECK_ITERABLE_APPROX(result.d2_dt_psi, expected.d2_dt_psi);
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.CSW.Worldtube.SelfForceDerivatives",
                  "[Unit][Evolution]") {
  test_reconstruction();
#ifdef SPECTRE_AUTODIFF
  test_time_derivatives();
  test_curved_background();
  test_curved_time_derivative();
#else
  const tnsr::I<double, 3> position{{8., 6., 0.}};
  const tnsr::I<double, 3> velocity{{-0.04, 0.08, 0.}};
  CHECK_THROWS_WITH(
      second_order_self_force_derivatives(plane_wave(position, 2.7), position,
                                          velocity, charge, rest_mass, true, 1.,
                                          2.7, turn_on_time, turn_on_interval),
      Catch::Matchers::ContainsSubstring("require SPECTRE_AUTODIFF=ON"));
#endif
}
}  // namespace CurvedScalarWave::Worldtube
