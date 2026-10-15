// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"
#include "Framework/TestHelpers.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave {
namespace {
using deriv_psi_tag =
    ::Tags::deriv<Tags::Psi, tmpl::size_t<3>, Frame::Inertial>;
using puncture_vars =
    Variables<tmpl::list<Tags::Psi, ::Tags::dt<Tags::Psi>, deriv_psi_tag>>;

struct AccelerationData {
  std::array<double, 3> f;
  std::array<double, 3> dt_f;
  std::array<double, 3> Du_f;
  std::array<double, 3> dt_Du_f;
  std::array<double, 3> DuDu_f;
  std::array<double, 3> dt_DuDu_f;
};

AccelerationData sample_acceleration_data() {
  return {{1. / 230., -1. / 110., 1. / 170.},
          {1. / 3700., -1. / 2900., 1. / 4300.},
          {-1. / 500., 1. / 370., 1. / 630.},
          {-1. / 5300., 1. / 4700., 1. / 5900.},
          {1. / 1900., -1. / 1700., 1. / 2300.},
          {1. / 7100., -1. / 6700., 1. / 7900.}};
}

puncture_vars evaluate(const tnsr::I<DataVector, 3>& centered_coords,
                       const tnsr::I<double, 3>& position,
                       const tnsr::I<double, 3>& velocity,
                       const tnsr::I<double, 3>& coordinate_acceleration,
                       const AccelerationData& a, const double mass) {
  // Also exercise initialization when the output has the wrong size.
  puncture_vars result{2};
  Worldtube::acceleration_terms_2(
      make_not_null(&result), centered_coords, position, velocity,
      coordinate_acceleration, a.f[0], a.f[1], a.f[2], a.dt_f[0], a.dt_f[1],
      a.dt_f[2], a.Du_f[0], a.Du_f[1], a.Du_f[2], a.dt_Du_f[0], a.dt_Du_f[1],
      a.dt_Du_f[2], a.DuDu_f[0], a.DuDu_f[1], a.DuDu_f[2], a.dt_DuDu_f[0],
      a.dt_DuDu_f[1], a.dt_DuDu_f[2], mass);
  return result;
}

tnsr::I<DataVector, 3> sample_points() {
  tnsr::I<DataVector, 3> result{size_t{3}};
  get<0>(result) = DataVector{0.17, -0.23, 0.31};
  get<1>(result) = DataVector{-0.11, 0.19, 0.07};
  get<2>(result) = DataVector{0.09, -0.13, 0.21};
  return result;
}

void test_dispatch() {
  const auto coords = sample_points();
  const tnsr::I<double, 3> position{{-6., 8., 0.}};
  const tnsr::I<double, 3> velocity{{-0.04, 1. / 12., 0.}};
  const tnsr::I<double, 3> acceleration{{-1. / 700., 1. / 1100., 0.}};
  const auto a = sample_acceleration_data();
  const auto expected =
      evaluate(coords, position, velocity, acceleration, a, 1.);
  const Worldtube::PunctureField puncture{
      Worldtube::PunctureField::Schwarzschild{2, 1.}};
  puncture_vars result{};
  puncture.apply_acceleration_terms(
      make_not_null(&result), coords, position, velocity, acceleration, a.f[0],
      a.f[1], a.f[2], a.dt_f[0], a.dt_f[1], a.dt_f[2], a.Du_f[0], a.Du_f[1],
      a.Du_f[2], a.dt_Du_f[0], a.dt_Du_f[1], a.dt_Du_f[2],
      std::array<double, 6>{{a.DuDu_f[0], a.DuDu_f[1], a.DuDu_f[2],
                             a.dt_DuDu_f[0], a.dt_DuDu_f[1], a.dt_DuDu_f[2]}});
  CHECK_VARIABLES_APPROX(result, expected);
  CHECK_THROWS_WITH(
      puncture.apply_acceleration_terms(
          make_not_null(&result), coords, position, velocity, acceleration,
          a.f[0], a.f[1], a.f[2], a.dt_f[0], a.dt_f[1], a.dt_f[2], a.Du_f[0],
          a.Du_f[1], a.Du_f[2], a.dt_Du_f[0], a.dt_Du_f[1], a.dt_Du_f[2]),
      Catch::Matchers::ContainsSubstring("require D_u^2 f"));
}

void test_mathematica_reference() {
  // Evaluated at 35-digit precision from the corrected Mathematica graph,
  // before SymPy CSE, integer-power expansion, or vector-buffer reuse.
  // The input components are deliberately independent to exercise all inputs.
  const std::array<std::array<double, 5>, 4> expected{
      {{{0.00516081411271969485, 0.000355180982272231497,
         0.00239538427795082192, -0.000723022550794723062,
         -0.00615187781710198601}},
       {{0.00434874033743261853, 0.000984475607389477640,
         0.00220279951566647045, -0.00877985730166446076,
         -0.00742035032589388452}},
       {{0.00410194107951176822, 0.00102941695644192179, 0.00341245872491572217,
         -0.00884726711806267524, -0.00585810005463911629}},
       {{0.00338809119410003104, 0.000972716222534694328,
         0.00385213933337423895, -0.00828948474421676578,
         -0.00412134855801704270}}}};
  const std::array<double, 4> x_positions{{-6., 9., -15., 21.}};
  const std::array<double, 4> y_positions{{8., 12., 20., 28.}};
  const tnsr::I<double, 3> velocity{{-1. / 25., 1. / 12., 0.}};
  const tnsr::I<double, 3> coordinate_acceleration{
      {-1. / 700., 1. / 1100., 0.}};
  const auto a = sample_acceleration_data();
  const auto local_approx = Approx::custom().epsilon(2.e-13).scale(1.);
  for (size_t i = 0; i < expected.size(); ++i) {
    CAPTURE(i);
    const double n = static_cast<double>(i + 1);
    tnsr::I<DataVector, 3> coords{size_t{1}};
    get<0>(coords) = (2. + n) / 17.;
    get<1>(coords) = (-3. + n) / 19.;
    get<2>(coords) = (1. + n) / 23.;
    const tnsr::I<double, 3> position{
        {gsl::at(x_positions, i), gsl::at(y_positions, i), 0.}};
    const auto result =
        evaluate(coords, position, velocity, coordinate_acceleration, a, 1.);
    CHECK(get(get<Tags::Psi>(result))[0] ==
          local_approx(gsl::at(expected, i)[0]));
    CHECK(get(get<::Tags::dt<Tags::Psi>>(result))[0] ==
          local_approx(gsl::at(expected, i)[1]));
    for (size_t d = 0; d < 3; ++d) {
      CHECK(get<deriv_psi_tag>(result).get(d)[0] ==
            local_approx(gsl::at(gsl::at(expected, i), d + 2)));
    }
  }
}

void test_derivatives_and_zero_acceleration() {
  const auto coords = sample_points();
  const tnsr::I<double, 3> position{{-6., 8., 0.}};
  const tnsr::I<double, 3> velocity{{-1. / 25., 1. / 12., 0.}};
  const tnsr::I<double, 3> coordinate_acceleration{
      {-1. / 700., 1. / 1100., 0.}};
  const auto a = sample_acceleration_data();
  const auto result =
      evaluate(coords, position, velocity, coordinate_acceleration, a, 1.);
  const double step = 1.e-5;
  const auto local_approx = Approx::custom().epsilon(1.e-10).scale(1.);
  for (size_t d = 0; d < 3; ++d) {
    CAPTURE(d);
    auto plus_coords = coords;
    auto minus_coords = coords;
    plus_coords.get(d) += step;
    minus_coords.get(d) -= step;
    const auto plus = evaluate(plus_coords, position, velocity,
                               coordinate_acceleration, a, 1.);
    const auto minus = evaluate(minus_coords, position, velocity,
                                coordinate_acceleration, a, 1.);
    const DataVector finite_difference =
        (get(get<Tags::Psi>(plus)) - get(get<Tags::Psi>(minus))) / (2. * step);
    CHECK_ITERABLE_CUSTOM_APPROX(get<deriv_psi_tag>(result).get(d),
                                 finite_difference, local_approx);
  }

  // Move the particle and its acceleration data while holding the inertial
  // field point fixed. Thus centered_coords change by -velocity * dt. The
  // ordinary time derivatives of f, D_u f and D_u^2 f are independent inputs.
  const auto at_time = [&](const double dt) {
    auto shifted_coords = coords;
    auto shifted_position = position;
    auto shifted_velocity = velocity;
    auto shifted_a = a;
    for (size_t d = 0; d < 3; ++d) {
      shifted_coords.get(d) -= dt * velocity.get(d);
      shifted_position.get(d) += dt * velocity.get(d);
      shifted_velocity.get(d) += dt * coordinate_acceleration.get(d);
      gsl::at(shifted_a.f, d) += dt * gsl::at(a.dt_f, d);
      gsl::at(shifted_a.Du_f, d) += dt * gsl::at(a.dt_Du_f, d);
      gsl::at(shifted_a.DuDu_f, d) += dt * gsl::at(a.dt_DuDu_f, d);
    }
    return evaluate(shifted_coords, shifted_position, shifted_velocity,
                    coordinate_acceleration, shifted_a, 1.);
  };
  const auto plus = at_time(step);
  const auto minus = at_time(-step);
  const DataVector finite_difference =
      (get(get<Tags::Psi>(plus)) - get(get<Tags::Psi>(minus))) / (2. * step);
  CHECK_ITERABLE_CUSTOM_APPROX(get(get<::Tags::dt<Tags::Psi>>(result)),
                               finite_difference, local_approx);

  const auto zero = evaluate(coords, position, velocity,
                             coordinate_acceleration, AccelerationData{}, 1.);
  CHECK_VARIABLES_APPROX(zero, puncture_vars(3, 0.));
}

void test_uniform_acceleration_in_flat_spacetime() {
  // At the instantaneously stationary point of a hyperbolic worldline with
  // proper acceleration a in the x direction, the exact time-symmetric field
  // is [s sqrt(1 + a X + a^2 s^2/4)]^{-1}. Subtract 1/s and expand through
  // distance squared to obtain an independent nonlinear-acceleration check.
  const double a = 0.07;
  const double a2 = a * a;
  const double a3 = a2 * a;
  const AccelerationData data{{0., a, 0.},  {a2, 0., 0.}, {a2, 0., 0.},
                              {0., a3, 0.}, {0., a3, 0.}, {a2 * a2, 0., 0.}};
  const auto coords = sample_points();
  const auto result = evaluate(coords, tnsr::I<double, 3>{{10., 0., 0.}},
                               tnsr::I<double, 3>{{0., 0., 0.}},
                               tnsr::I<double, 3>{{a, 0., 0.}}, data, 0.);
  const auto local_approx = Approx::custom().epsilon(2.e-13).scale(1.);
  for (size_t p = 0; p < get<0>(coords).size(); ++p) {
    CAPTURE(p);
    const double x = get<0>(coords)[p];
    const double y = get<1>(coords)[p];
    const double z = get<2>(coords)[p];
    const double s = sqrt(x * x + y * y + z * z);
    const double s3 = s * s * s;
    const double expected = -a * x / (2. * s) + a2 * (3. * x * x / s - s) / 8. +
                            a3 * (3. * x * s - 5. * x * x * x / s) / 16.;
    CHECK(get(get<Tags::Psi>(result))[p] == local_approx(expected));
    CHECK(get(get<::Tags::dt<Tags::Psi>>(result))[p] == local_approx(0.));
    for (size_t d = 0; d < 3; ++d) {
      const double delta_x = d == 0 ? 1. : 0.;
      const double xi = coords.get(d)[p];
      const double expected_derivative =
          -a / 2. * (delta_x / s - x * xi / s3) +
          a2 / 8. * (6. * x * delta_x / s - 3. * x * x * xi / s3 - xi / s) +
          a3 / 16. *
              (3. * delta_x * s + 3. * x * xi / s - 15. * x * x * delta_x / s +
               5. * x * x * x * xi / s3);
      CHECK(get<deriv_psi_tag>(result).get(d)[p] ==
            local_approx(expected_derivative));
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.CurvedScalarWave.AccelerationTermsOrder2",
    "[Unit][Evolution]") {
  test_dispatch();
  test_mathematica_reference();
  test_derivatives_and_zero_acceleration();
  test_uniform_acceleration_in_flat_spacetime();
}
}  // namespace CurvedScalarWave
