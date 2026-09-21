// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <random>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/CoulombDecode.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
constexpr const char* module = "NpMatching";

double relative_l2(const ComplexDataVector& difference,
                   const ComplexDataVector& reference) {
  double diff = 0.;
  double ref = 0.;
  for (size_t p = 0; p < difference.size(); ++p) {
    diff += std::norm(difference[p]);
    ref += std::norm(reference[p]);
  }
  return sqrt(diff / std::max(ref, 1.e-300));
}

void test_background_and_radius_solve() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 12;
  std::uniform_real_distribution<> radius_dist(2.6, 12.);
  std::uniform_real_distribution<> factor_dist(0.7, 1.6);
  std::uniform_real_distribution<> noise_dist(-0.1, 0.1);
  const DataVector used_for_size(num_points);
  const double mass = 1.;
  const auto radius = make_with_random_values<DataVector>(
      make_not_null(&generator), make_not_null(&radius_dist), used_for_size);
  CHECK_ITERABLE_APPROX(
      background_radial_derivative(radius, mass),
      pypp::call<DataVector>(module, "coulomb_background_radial_derivative",
                             radius, mass));
  CHECK_ITERABLE_APPROX(
      background_radial_derivative(radius, mass),
      DataVector(sqrt(1. - 2. * mass / radius) * 3. * mass / pow<4>(radius)));

  // The solve recovers the radius from the normal derivative for any
  // (positive or negative) kinematic factor, from a perturbed start
  const auto factor = make_with_random_values<DataVector>(
      make_not_null(&generator), make_not_null(&factor_dist), used_for_size);
  const DataVector normal_derivative =
      -background_radial_derivative(radius, mass) * factor;
  const DataVector initial =
      radius * (1. + make_with_random_values<DataVector>(
                         make_not_null(&generator), make_not_null(&noise_dist),
                         used_for_size));
  const RadiusSolve solve =
      radius_from_normal_derivative(normal_derivative, factor, initial, mass);
  CHECK(solve.valid);
  Approx tight = Approx::custom().epsilon(1.e-9).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(solve.radius, radius, tight);
  CHECK_ITERABLE_CUSTOM_APPROX(
      solve.radius,
      pypp::call<DataVector>(module, "coulomb_radius_from_normal_derivative",
                             normal_derivative, factor, initial, mass),
      tight);
  CHECK(pypp::call<bool>(module, "coulomb_radius_solve_valid",
                         normal_derivative, factor, initial, mass));

  // A normal derivative above the maximum of B(r), reached at the turning
  // point 9M/4, has no solution: the solve reports it. (A point inside the
  // turning point itself has a valid outer-branch solution with the same
  // derivative, which the decode cannot distinguish.)
  DataVector inside_derivative = normal_derivative;
  inside_derivative[3] =
      1.2 * background_radial_derivative(DataVector{2.25 * mass}, mass)[0] *
      factor[3];
  const RadiusSolve failed =
      radius_from_normal_derivative(inside_derivative, factor, radius, mass);
  CHECK_FALSE(failed.valid);
  CHECK_FALSE(pypp::call<bool>(module, "coulomb_radius_solve_valid",
                               inside_derivative, factor, radius, mass));
}

// The manufactured slice (exact type-D background at areal radius 0.5 with
// a quadrupole tide of relative size 1e-3 and a boost of speed 0.22): with
// the exact normal derivative of the background Coulomb scalar the decode
// recovers the radius, the manufactured moments and the held-out Psi0, and
// agrees with the reference implementation.
void test_decode_on_manufactured_slice() {
  const size_t num_points = 240;
  const auto slice = helpers::manufactured_slice(71, num_points);
  const auto& geometry = slice.geometry;
  const Scalar<DataVector> rapidity{
      pypp::call<DataVector>(module, "manufactured_rapidity", 71, num_points)};
  const auto psi_list = helpers::pack_scalars(slice.psi);
  const auto reals = helpers::pack_reals(geometry);
  const double true_radius = 0.5;

  const FrameRegistration registration =
      register_frame(slice.psi, geometry.rotation, slice.mass);
  const DataVector factor =
      normal_derivative_factor(registration.member, rapidity);
  CHECK_ITERABLE_APPROX(
      factor,
      pypp::call<DataVector>(module, "coulomb_normal_derivative_factor",
                             psi_list, reals, slice.mass, get(rapidity)));
  const Scalar<DataVector> normal_derivative{
      background_radial_derivative(DataVector(num_points, true_radius),
                                   slice.mass) *
      factor};

  // Columns against the reference
  const auto columns = coulomb_tide_columns(
      registration.radial_direction, registration.measured_radius, slice.mass);
  for (const size_t index : {0_st, 3_st, 5_st, 9_st}) {
    CHECK_ITERABLE_APPROX(
        get(gsl::at(columns, index)),
        pypp::call<ComplexDataVector>(
            module, "coulomb_tide_column", psi_list, reals, slice.mass,
            get(registration.measured_radius), index));
  }

  const CoulombDecode decode = decode_tidal_moments_from_coulomb(
      registration, rapidity, slice.mass, normal_derivative);
  REQUIRE(decode.valid);
  Approx custom_approx = Approx::custom().epsilon(1.e-9).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get(decode.areal_radius),
      pypp::call<DataVector>(module, "coulomb_decode_radius", psi_list, reals,
                             slice.mass, get(rapidity), get(normal_derivative),
                             std::vector<double>{}),
      custom_approx);
  const auto expected_components = pypp::call<std::array<double, 10>>(
      module, "coulomb_decode_components", psi_list, reals, slice.mass,
      get(rapidity), get(normal_derivative), std::vector<double>{});
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::real(gsl::at(decode.components, a)) ==
          custom_approx(gsl::at(expected_components, a)));
    CHECK(std::imag(gsl::at(decode.components, a)) ==
          custom_approx(gsl::at(expected_components, 5 + a)));
  }
  CHECK(decode.relative_residual ==
        custom_approx(pypp::call<double>(
            module, "coulomb_decode_residual", psi_list, reals, slice.mass,
            get(rapidity), get(normal_derivative), std::vector<double>{})));

  // Physics: the radius is the manufactured one (the tide does not enter the
  // normal derivative); the moments are the manufactured ones to 0.5%: the
  // decode reads a signal of relative size 1e-4 off the Coulomb scalar, and
  // the O(v eps) aberration of the measured radial direction that the tidal
  // model identifies with the rest-frame one is what is left; the target
  // reproduces the held-out Psi0
  CHECK_ITERABLE_CUSTOM_APPROX(get(decode.areal_radius),
                               DataVector(num_points, true_radius),
                               Approx::custom().epsilon(1.e-8).scale(1.));
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::abs(std::real(gsl::at(decode.components, a)) -
                   gsl::at(slice.truth, a)) <
          1.e-2 * std::abs(gsl::at(slice.truth, a)) + 1.e-9);
    CHECK(std::abs(std::imag(gsl::at(decode.components, a)) -
                   gsl::at(slice.truth, 5 + a)) <
          1.e-2 * std::abs(gsl::at(slice.truth, 5 + a)) + 1.e-9);
  }
  const SecondOrderEvaluation evaluation =
      evaluate_second_order(registration, rapidity, geometry.rotation,
                            slice.mass, std::nullopt, decode.components);
  CHECK(relative_l2(get(evaluation.psi0_target) - slice.psi.get(0),
                    slice.psi.get(0)) < 1.e-3);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.CoulombDecode",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_background_and_radius_solve();
  test_decode_on_manufactured_slice();
}
}  // namespace gr::np
