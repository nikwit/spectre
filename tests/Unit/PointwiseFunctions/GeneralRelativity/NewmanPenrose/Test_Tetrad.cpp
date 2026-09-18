// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <complex>
#include <cstddef>
#include <random>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
using frame = Frame::Inertial;
using SymmetricTensor = tnsr::ii<DataVector, 3, frame>;
const std::string module = "NpMatching";

void test_cholesky_and_triad() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const auto geometry =
      helpers::random_geometry(make_not_null(&generator), num_points);
  const RealMatrix lower = cholesky_factor(geometry.spatial_metric);
  CHECK_ITERABLE_APPROX(lower, pypp::call<RealMatrix>(module, "cholesky_factor",
                                                      geometry.spatial_metric));
  CHECK_ITERABLE_APPROX(
      inverse_lower_triangular(lower),
      pypp::call<RealMatrix>(module, "inverse_lower_triangular", lower));
  CHECK_ITERABLE_APPROX(
      adapted_triad(geometry.spatial_metric, geometry.directions),
      pypp::call<RealMatrix>(module, "adapted_triad", geometry.spatial_metric,
                             geometry.directions));

  // The triad rows are orthonormal in the Cholesky triad and the first row is
  // the unit normal covector of the coordinate sphere: L^{-1} d, normalized
  const RealMatrix inverse = inverse_lower_triangular(lower);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      DataVector product(num_points, 0.);
      for (size_t k = 0; k < 3; ++k) {
        product += geometry.rotation.get(i, k) * geometry.rotation.get(j, k);
      }
      CHECK_ITERABLE_APPROX(product, DataVector(num_points, i == j ? 1. : 0.));
    }
  }

  // At a coordinate pole the azimuth is fixed to zero
  tnsr::I<DataVector, 3, frame> pole(num_points, 0.);
  get<2>(pole) = 1.;
  CHECK_ITERABLE_APPROX(adapted_triad(geometry.spatial_metric, pole),
                        pypp::call<RealMatrix>(module, "adapted_triad",
                                               geometry.spatial_metric, pole));
}

void test_scalars_dictionary() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const auto geometry =
      helpers::random_geometry(make_not_null(&generator), num_points);
  // A symmetric trace-free complex Q
  ComplexMatrix q(num_points, std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      q.get(i, j) =
          helpers::random_complex(make_not_null(&generator), num_points, 1.);
      q.get(j, i) = q.get(i, j);
    }
  }
  const ComplexDataVector trace = q.get(0, 0) + q.get(1, 1) + q.get(2, 2);
  for (size_t i = 0; i < 3; ++i) {
    q.get(i, i) -= trace / 3.;
  }
  CHECK_ITERABLE_APPROX(rotate_symmetric(q, geometry.rotation),
                        pypp::call<ComplexMatrix>(module, "rotate_symmetric", q,
                                                  geometry.rotation));
  const WeylScalars psi = weyl_scalars_from_tidal_tensor(q);
  CHECK_ITERABLE_APPROX(psi, pypp::call<WeylScalars>(
                                 module, "weyl_scalars_from_tidal_tensor", q));
  CHECK_ITERABLE_APPROX(
      tidal_tensor_from_weyl_scalars(psi),
      pypp::call<ComplexMatrix>(module, "tidal_tensor_from_weyl_scalars", psi));
  // Round trip
  CHECK_ITERABLE_APPROX(tidal_tensor_from_weyl_scalars(psi), q);

  // From the electric and magnetic parts in coordinate components
  std::uniform_real_distribution<> dist(-1., 1.);
  const DataVector used_for_size(num_points);
  const auto electric = make_with_random_values<tnsr::ii<DataVector, 3, frame>>(
      make_not_null(&generator), make_not_null(&dist), used_for_size);
  const auto magnetic = make_with_random_values<tnsr::ii<DataVector, 3, frame>>(
      make_not_null(&generator), make_not_null(&dist), used_for_size);
  CHECK_ITERABLE_APPROX(
      weyl_scalars_from_electric_magnetic(
          electric, magnetic, geometry.spatial_metric, geometry.directions),
      pypp::call<WeylScalars>(module, "weyl_scalars_from_electric_magnetic",
                              electric, magnetic, geometry.spatial_metric,
                              geometry.directions));
}

void test_incoming_field() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const auto geometry =
      helpers::random_geometry(make_not_null(&generator), num_points);
  const Scalar<ComplexDataVector> psi0{
      helpers::random_complex(make_not_null(&generator), num_points, 1.)};
  const auto w_minus = incoming_weyl_field(psi0, geometry.rotation);
  CHECK_ITERABLE_APPROX(
      w_minus, pypp::call<SymmetricTensor>(module, "incoming_weyl_field", psi0,
                                           geometry.rotation));
  // |w^-|_F^2 = 8 |Psi0|^2, trace free and tangent to the cut
  DataVector frobenius(num_points, 0.);
  DataVector trace(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    trace += w_minus.get(i, i);
    for (size_t j = 0; j < 3; ++j) {
      frobenius += square(w_minus.get(i, j));
    }
  }
  CHECK_ITERABLE_APPROX(frobenius, DataVector{8. * square(abs(get(psi0)))});
  CHECK_ITERABLE_APPROX(trace, DataVector(num_points, 0.));
  for (size_t i = 0; i < 3; ++i) {
    DataVector contraction(num_points, 0.);
    for (size_t j = 0; j < 3; ++j) {
      contraction += w_minus.get(i, j) * geometry.rotation.get(0, j);
    }
    CHECK_ITERABLE_APPROX(contraction, DataVector(num_points, 0.));
  }
  // Dyad-spin invariance: rotating the tangent legs by chi sends
  // m -> exp(-i chi) m and Psi0 -> exp(-2 i chi) Psi0, and leaves w^- unchanged
  const double chi = 0.7;
  RealMatrix spun = geometry.rotation;
  for (size_t i = 0; i < 3; ++i) {
    spun.get(1, i) = cos(chi) * geometry.rotation.get(1, i) +
                     sin(chi) * geometry.rotation.get(2, i);
    spun.get(2, i) = -sin(chi) * geometry.rotation.get(1, i) +
                     cos(chi) * geometry.rotation.get(2, i);
  }
  const Scalar<ComplexDataVector> spun_psi0{
      get(psi0) * std::exp(std::complex<double>{0., -2. * chi})};
  CHECK_ITERABLE_APPROX(incoming_weyl_field(spun_psi0, spun), w_minus);

  const RealMatrix lower = cholesky_factor(geometry.spatial_metric);
  CHECK_ITERABLE_APPROX(
      orthonormal_to_coordinate_covariant(w_minus, lower),
      pypp::call<SymmetricTensor>(module, "orthonormal_to_coordinate_covariant",
                                  w_minus, lower));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.Tetrad",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_cholesky_and_triad();
  test_scalars_dictionary();
  test_incoming_field();
}
}  // namespace gr::np
