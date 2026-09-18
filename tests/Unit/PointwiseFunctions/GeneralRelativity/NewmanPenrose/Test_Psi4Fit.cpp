// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
const std::string module = "NpMatching";

double relative_l2(const ComplexDataVector& difference,
                   const ComplexDataVector& reference) {
  double num = 0.;
  double den = 0.;
  for (size_t p = 0; p < difference.size(); ++p) {
    num += std::norm(difference[p]);
    den += std::norm(reference[p]);
  }
  return sqrt(num / den);
}

void test_fit_psi4() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 60;
  std::uniform_real_distribution<> positive(0.5, 1.5);
  const DataVector used_for_size(num_points);
  // An arbitrary complex design with known coefficients
  std::array<std::complex<double>, 5> truth{};
  std::array<Scalar<ComplexDataVector>, 5> columns{};
  std::vector<ComplexDataVector> column_list{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(truth, a) =
        helpers::random_complex(make_not_null(&generator), 1, 1.)[0];
    get(gsl::at(columns, a)) =
        helpers::random_complex(make_not_null(&generator), num_points, 1.);
    column_list.push_back(get(gsl::at(columns, a)));
  }
  Scalar<ComplexDataVector> data{
      ComplexDataVector(num_points, std::complex<double>{0., 0.})};
  for (size_t a = 0; a < 5; ++a) {
    get(data) += gsl::at(truth, a) * get(gsl::at(columns, a));
  }
  Approx custom_approx = Approx::custom().epsilon(1.e-10).scale(1.);
  {
    const DirectPsi4Fit fit = fit_psi4(data, columns);
    for (size_t a = 0; a < 5; ++a) {
      CHECK(std::real(gsl::at(fit.components, a)) ==
            custom_approx(std::real(gsl::at(truth, a))));
      CHECK(std::imag(gsl::at(fit.components, a)) ==
            custom_approx(std::imag(gsl::at(truth, a))));
    }
    CHECK(fit.relative_residual == custom_approx(0.));
    const auto expected = pypp::call<std::array<double, 10>>(
        module, "fit_psi4", get(data), column_list, std::vector<double>{});
    for (size_t a = 0; a < 5; ++a) {
      CHECK(fit.electric()[a] == custom_approx(gsl::at(expected, a)));
      CHECK(fit.magnetic()[a] == custom_approx(gsl::at(expected, 5 + a)));
    }
  }
  // Inconsistent data with weights: compare to the reference solver
  {
    get(data) +=
        helpers::random_complex(make_not_null(&generator), num_points, 0.3);
    const auto weights = make_with_random_values<DataVector>(
        make_not_null(&generator), make_not_null(&positive), used_for_size);
    const DirectPsi4Fit fit = fit_psi4(data, columns, weights);
    const auto expected = pypp::call<std::array<double, 10>>(
        module, "fit_psi4", get(data), column_list,
        std::vector<double>(weights.begin(), weights.end()));
    for (size_t a = 0; a < 5; ++a) {
      CHECK(fit.electric()[a] == custom_approx(gsl::at(expected, a)));
      CHECK(fit.magnetic()[a] == custom_approx(gsl::at(expected, 5 + a)));
    }
    // A pointwise phase on data and columns leaves the fit unchanged
    Scalar<ComplexDataVector> spun_data = data;
    auto spun_columns = columns;
    for (size_t p = 0; p < num_points; ++p) {
      const std::complex<double> phase =
          std::exp(std::complex<double>{0., 2. * positive(generator)});
      get(spun_data)[p] *= phase;
      for (size_t a = 0; a < 5; ++a) {
        get(gsl::at(spun_columns, a))[p] *= phase;
      }
    }
    const DirectPsi4Fit spun_fit = fit_psi4(spun_data, spun_columns, weights);
    for (size_t a = 0; a < 5; ++a) {
      CHECK(std::real(gsl::at(spun_fit.components, a)) ==
            custom_approx(std::real(gsl::at(fit.components, a))));
      CHECK(std::imag(gsl::at(spun_fit.components, a)) ==
            custom_approx(std::imag(gsl::at(fit.components, a))));
    }
  }
}

// The manufactured slice of the minimal study, with the exact type-III
// rapidity of its boost supplied: the full chain reproduces the reference
// implementation and recovers the manufactured moments and the held-out Psi0.
void test_second_order_on_manufactured_slice() {
  const size_t num_points = 240;
  const auto slice = helpers::manufactured_slice(71, num_points);
  const auto& geometry = slice.geometry;
  const Scalar<DataVector> rapidity{
      pypp::call<DataVector>(module, "manufactured_rapidity", 71, num_points)};

  const FrameRegistration registration =
      register_frame(slice.psi, geometry.rotation, slice.mass);
  const SecondOrderEvaluation evaluation = evaluate_second_order(
      registration, rapidity, geometry.rotation, slice.mass);

  const auto psi_list = helpers::pack_scalars(slice.psi);
  const auto reals = helpers::pack_reals(geometry);
  Approx custom_approx = Approx::custom().epsilon(1.e-9).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get(registration.measured_radius),
      pypp::call<DataVector>(module, "register_frame_radius", psi_list, reals,
                             slice.mass),
      custom_approx);
  const auto expected_components = pypp::call<std::array<double, 10>>(
      module, "evaluate_second_order_components", psi_list, reals, slice.mass,
      get(rapidity));
  for (size_t a = 0; a < 5; ++a) {
    CHECK(evaluation.fit.electric()[a] ==
          custom_approx(gsl::at(expected_components, a)));
    CHECK(evaluation.fit.magnetic()[a] ==
          custom_approx(gsl::at(expected_components, 5 + a)));
  }
  const auto expected_target =
      pypp::call<ComplexDataVector>(module, "evaluate_second_order_psi0_target",
                                    psi_list, reals, slice.mass, get(rapidity));
  CHECK(relative_l2(get(evaluation.psi0_target) - expected_target,
                    expected_target) < 1.e-10);

  // Physics: the manufactured slice is an exact type-D background with a
  // quadrupole tide of relative size 1e-3 and a boost of speed 0.22. The
  // radius is 0.5 up to the O(eps^2) tidal contamination of the leading-order
  // measurement, the fitted moments are the manufactured ones up to the
  // O(eps^4) remainder, and the held-out Psi0 is reproduced. The leading
  // type-D frame leaves O(eps^2) longitudinal residuals in the pulled-back
  // Psi1 and Psi3, which contaminate the target at 4e-4 relative here; the
  // quasi-Kinnersley Newton correction of the reference brings that to 1e-6
  // and is not part of this library.
  Approx loose = Approx::custom().epsilon(1.e-3).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(get(registration.measured_radius),
                               DataVector(num_points, 0.5), loose);
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::abs(evaluation.fit.electric()[a] - gsl::at(slice.truth, a)) <
          1.e-3 * std::abs(gsl::at(slice.truth, a)) + 1.e-9);
    CHECK(std::abs(evaluation.fit.magnetic()[a] - gsl::at(slice.truth, 5 + a)) <
          1.e-3 * std::abs(gsl::at(slice.truth, 5 + a)) + 1.e-9);
  }
  CHECK(relative_l2(get(evaluation.psi0_target) - slice.psi.get(0),
                    slice.psi.get(0)) < 1.e-3);

  // Blind fixed point: with Psi0 hidden from the frame measurement the
  // iteration Psi0 <- target converges onto the held-out value
  WeylScalars working = slice.psi;
  working.get(0) = 0.;
  double change = 1.;
  for (size_t iteration = 0; iteration < 12 and change > 1.e-10; ++iteration) {
    const SecondOrderEvaluation blind = evaluate_second_order(
        register_frame(working, geometry.rotation, slice.mass), rapidity,
        geometry.rotation, slice.mass);
    change = relative_l2(get(blind.psi0_target) - working.get(0),
                         get(blind.psi0_target));
    working.get(0) = get(blind.psi0_target);
  }
  CHECK(change < 1.e-10);
  CHECK(relative_l2(working.get(0) - slice.psi.get(0), slice.psi.get(0)) <
        1.e-3);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.Psi4Fit",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_fit_psi4();
  test_second_order_on_manufactured_slice();
}
}  // namespace gr::np
