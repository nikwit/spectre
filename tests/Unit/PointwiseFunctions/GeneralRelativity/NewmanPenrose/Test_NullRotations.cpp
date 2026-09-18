// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <complex>
#include <cstddef>
#include <random>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
const std::string module = "NpMatching";

void test_rotations_and_invariants() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const WeylScalars psi =
      helpers::random_scalars(make_not_null(&generator), num_points, 1.);
  const Scalar<ComplexDataVector> a_bar{
      helpers::random_complex(make_not_null(&generator), num_points, 0.5)};
  const Scalar<ComplexDataVector> b{
      helpers::random_complex(make_not_null(&generator), num_points, 0.5)};
  std::uniform_real_distribution<> dist(-0.5, 0.5);
  const DataVector used_for_size(num_points);
  const auto eta = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), make_not_null(&dist), used_for_size);
  const auto chi = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), make_not_null(&dist), used_for_size);

  CHECK_ITERABLE_APPROX(type_i(psi, a_bar),
                        pypp::call<WeylScalars>(module, "type_i", psi, a_bar));
  CHECK_ITERABLE_APPROX(type_ii(psi, b),
                        pypp::call<WeylScalars>(module, "type_ii", psi, b));
  CHECK_ITERABLE_APPROX(
      type_iii(psi, eta, chi),
      pypp::call<WeylScalars>(module, "type_iii", psi, eta, chi));

  const Scalar<ComplexDataVector> i = invariant_i(psi);
  const Scalar<ComplexDataVector> j = invariant_j(psi);
  CHECK_ITERABLE_APPROX(
      i, pypp::call<Scalar<ComplexDataVector>>(module, "invariant_i", psi));
  CHECK_ITERABLE_APPROX(
      j, pypp::call<Scalar<ComplexDataVector>>(module, "invariant_j", psi));

  // The rotations form groups and the invariants are invariant
  const Scalar<ComplexDataVector> minus_a_bar{-get(a_bar)};
  const Scalar<ComplexDataVector> minus_b{-get(b)};
  const Scalar<DataVector> minus_eta{-get(eta)};
  const Scalar<DataVector> minus_chi{-get(chi)};
  CHECK_ITERABLE_APPROX(type_i(type_i(psi, a_bar), minus_a_bar), psi);
  CHECK_ITERABLE_APPROX(type_ii(type_ii(psi, b), minus_b), psi);
  CHECK_ITERABLE_APPROX(type_iii(type_iii(psi, eta, chi), minus_eta, minus_chi),
                        psi);
  Approx custom_approx = Approx::custom().epsilon(1.e-11).scale(1.);
  const WeylScalars transformed =
      type_iii(type_i(type_ii(psi, b), a_bar), eta, chi);
  CHECK_ITERABLE_CUSTOM_APPROX(invariant_i(transformed), i, custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(invariant_j(transformed), j, custom_approx);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.NullRotations",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_rotations_and_invariants();
}
}  // namespace gr::np
