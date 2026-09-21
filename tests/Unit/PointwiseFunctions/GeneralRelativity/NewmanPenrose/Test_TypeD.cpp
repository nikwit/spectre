// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <complex>
#include <cstddef>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
const std::string module = "NpMatching";

void test_type_d() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const double mass = 1. / 9.;
  const auto data =
      helpers::near_type_d(make_not_null(&generator), num_points, 0.1, 1.e-3);
  const Scalar<ComplexDataVector> i = invariant_i(data.psi);
  const Scalar<ComplexDataVector> j = invariant_j(data.psi);

  const Scalar<ComplexDataVector> coulomb = coulomb_scalar(i, j);
  CHECK_ITERABLE_APPROX(coulomb, pypp::call<Scalar<ComplexDataVector>>(
                                     module, "coulomb_scalar", i, j));
  CHECK_ITERABLE_APPROX(background_radius(i, j, mass),
                        pypp::call<Scalar<ComplexDataVector>>(
                            module, "background_radius", i, j, mass));
  CHECK_ITERABLE_APPROX(
      kinnersley_scalars(coulomb),
      pypp::call<WeylScalars>(module, "kinnersley_scalars", coulomb));
  CHECK_ITERABLE_APPROX(type_d_scalars(coulomb, data.a_bar, data.b),
                        pypp::call<WeylScalars>(module, "type_d_scalars",
                                                coulomb, data.a_bar, data.b));

  const TypeDRotation rotation = solve_type_d_rotation(data.psi, coulomb);
  CHECK_ITERABLE_APPROX(rotation.a_bar,
                        pypp::call<Scalar<ComplexDataVector>>(
                            module, "solve_type_d_a_bar", data.psi, coulomb));
  CHECK_ITERABLE_APPROX(rotation.b,
                        pypp::call<Scalar<ComplexDataVector>>(
                            module, "solve_type_d_b", data.psi, coulomb));
  CHECK_ITERABLE_APPROX(rotation.x,
                        pypp::call<Scalar<ComplexDataVector>>(
                            module, "solve_type_d_x", data.psi, coulomb));
  CHECK_ITERABLE_APPROX(
      rotation.predicted_psi,
      pypp::call<WeylScalars>(module, "solve_type_d_predicted_psi", data.psi,
                              coulomb));
  // Slots 1 and 2 are reproduced by construction
  CHECK_ITERABLE_APPROX(rotation.predicted_psi.get(1), data.psi.get(1));
  CHECK_ITERABLE_APPROX(rotation.predicted_psi.get(2), data.psi.get(2));

  CHECK_ITERABLE_APPROX(pull_back(data.psi, rotation.a_bar, rotation.b),
                        pypp::call<WeylScalars>(module, "pull_back", data.psi,
                                                rotation.a_bar, rotation.b));
  Approx round_trip_approx = Approx::custom().epsilon(1.e-12).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(
      push_forward(pull_back(data.psi, rotation.a_bar, rotation.b),
                   rotation.a_bar, rotation.b),
      data.psi, round_trip_approx);
  CHECK_ITERABLE_APPROX(psi0_leading(coulomb, rotation.b),
                        pypp::call<Scalar<ComplexDataVector>>(
                            module, "psi0_leading", coulomb, rotation.b));

  // On exact type-D data the solve recovers the rotation and the pull-back
  // lands on the Kinnersley scalars
  const auto exact =
      helpers::near_type_d(make_not_null(&generator), num_points, 0.1, 0.);
  const Scalar<ComplexDataVector> exact_coulomb =
      coulomb_scalar(invariant_i(exact.psi), invariant_j(exact.psi));
  Approx custom_approx = Approx::custom().epsilon(1.e-10).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(exact_coulomb, exact.coulomb, custom_approx);
  const TypeDRotation exact_rotation =
      solve_type_d_rotation(exact.psi, exact_coulomb);
  CHECK_ITERABLE_CUSTOM_APPROX(
      Scalar<ComplexDataVector>{get(exact_rotation.a_bar) *
                                get(exact_rotation.b)},
      Scalar<ComplexDataVector>{get(exact.a_bar) * get(exact.b)},
      custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(
      pull_back(exact.psi, exact_rotation.a_bar, exact_rotation.b),
      kinnersley_scalars(exact_coulomb), custom_approx);
}

// A tetrad already aligned with the principal null directions: Psi1 and
// Psi3 vanish up to roundoff, so a_bar = x / b is a ratio of roundoff-sized
// numbers. The solve must return the aligned rotation, not O(1) noise.
void test_aligned_limit() {
  const size_t num_points = 3;
  const std::complex<double> coulomb_value{-0.064, 0.};
  Scalar<ComplexDataVector> coulomb(num_points, coulomb_value);
  WeylScalars psi(num_points, std::complex<double>{0., 0.});
  // Deterministic roundoff-sized contamination: Psi2 off by 1e-14 relative,
  // Psi1 and Psi3 at 1e-15 of the Coulomb scalar
  for (size_t p = 0; p < num_points; ++p) {
    psi.get(2)[p] = coulomb_value * (1. + 1.e-14);
    psi.get(1)[p] = 1.e-15 * std::abs(coulomb_value);
    psi.get(3)[p] = -1.e-15 * std::abs(coulomb_value);
  }
  const TypeDRotation rotation = solve_type_d_rotation(psi, coulomb);
  for (size_t p = 0; p < num_points; ++p) {
    CHECK(get(rotation.a_bar)[p] == std::complex<double>{0., 0.});
    CHECK(std::abs(get(rotation.b)[p]) < 1.e-14);
  }
  CHECK_ITERABLE_CUSTOM_APPROX(rotation.predicted_psi, psi,
                               Approx::custom().epsilon(1.e-12).scale(0.1));
  // Without the alignment threshold the same input amplifies the roundoff
  const TypeDRotation unregularized = solve_type_d_rotation(psi, coulomb, 0.);
  CHECK(std::abs(get(unregularized.a_bar)[0]) > 0.1);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.TypeD",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_type_d();
  test_aligned_limit();
}
}  // namespace gr::np
