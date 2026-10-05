// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPoints.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace gr::np {
namespace {
std::pair<TriadVector, DataVector> sphere(const size_t lm, const double tilt) {
  const size_t nt = lm + 1, nf = 2 * lm + 1, np = nt * nf;
  const auto& z = Spectral::collocation_points<Spectral::Basis::Legendre,
                                               Spectral::Quadrature::Gauss>(nt);
  const auto& w = Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                               Spectral::Quadrature::Gauss>(nt);
  TriadVector labels(np, 0.);
  DataVector weights(np);
  for (size_t j = 0; j < nf; ++j) {
    const double phi = 2. * acos(-1.) * static_cast<double>(j) / nf;
    for (size_t i = 0; i < nt; ++i) {
      const size_t k = i + nt * j;
      const double x = sqrt(1. - square(z[i])) * cos(phi);
      const double y = sqrt(1. - square(z[i])) * sin(phi);
      labels.get(0)[k] = cos(tilt) * x + sin(tilt) * z[i];
      labels.get(1)[k] = y;
      labels.get(2)[k] = -sin(tilt) * x + cos(tilt) * z[i];
      weights[k] = w[i] / (2. * nf);
    }
  }
  return {labels, weights};
}

void test_round_map(const double tilt, const double scale) {
  const auto [labels, weights] = sphere(8, tilt);
  const size_t np = weights.size();
  const std::array<DataVector, 3> metric{
      {DataVector(np, scale), DataVector(np, 0.), DataVector(np, scale)}};
  const auto map = laplace_eigenmap(labels, weights, metric, 8);
  const auto tolerance = Approx::custom().epsilon(2.e-10).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(map.direction, labels, tolerance);
  CHECK_ITERABLE_CUSTOM_APPROX(map.weights, weights, tolerance);
  CHECK(map.minimum_jacobian == tolerance(1.));
  CHECK(map.area_over_label_area == tolerance(1.));
  CHECK_ITERABLE_CUSTOM_APPROX(map.eigenvalues,
                               (std::array<double, 5>{{0., 2., 2., 2., 6.}}),
                               tolerance);
  for (size_t p = 0; p < np; ++p) {
    const double x = labels.get(0)[p], y = labels.get(1)[p],
                 z = labels.get(2)[p], st = hypot(x, y);
    const std::array<double, 3> et{{z * x / st, z * y / st, -st}},
        ep{{-y / st, x / st, 0.}};
    for (size_t i = 0; i < 3; ++i) {
      CHECK(map.derivative[0].get(i)[p] == tolerance(et[i]));
      CHECK(map.derivative[1].get(i)[p] == tolerance(ep[i]));
    }
  }
}

// Independent incoming target: construct Psi4 from an electric+magnetic STF
// tide on a Schwarzschild sphere and keep the measured Psi0 identically zero.
// A rotated collocation grid checks the BBH angular orientation assumption.
void test_prescribed_tide(const double tilt) {
  const auto [labels, weights] = sphere(8, tilt);
  const size_t np = weights.size();
  tnsr::ii<DataVector, 3, Frame::Inertial> metric(np, 0.);
  TriadVector inward = labels;
  for (size_t i = 0; i < 3; ++i) {
    inward.get(i) *= -1.;
    for (size_t j = i; j < 3; ++j) {
      metric.get(i, j) = .8 * labels.get(i) * labels.get(j);
      if (i == j) {
        metric.get(i, j) += 1.;
      }
    }
  }
  const auto rotation = adapted_triad(metric, inward);
  const double eta = -atanh(.8), mass = 1., radius = 2.5;
  const std::array<std::array<std::complex<double>, 3>, 3> tide{
      {{{{.4, .2}, {.12, -.17}, {.08, .1}}},
       {{{.12, -.17}, {-.1, .3}, {.11, -.03}}},
       {{{.08, .1}, {.11, -.03}, {-.3, -.5}}}}};
  WeylScalars psi(np, std::complex<double>{0., 0.});
  psi.get(2) = -mass / cube(radius);
  ComplexDataVector truth(np, std::complex<double>{0., 0.});
  for (size_t p = 0; p < np; ++p) {
    const double x = labels.get(0)[p], y = labels.get(1)[p],
                 z = labels.get(2)[p], st = hypot(x, y);
    const std::array<std::complex<double>, 3> mu{
        {{z * x / (st * sqrt(2.)), y / (st * sqrt(2.))},
         {z * y / (st * sqrt(2.)), -x / (st * sqrt(2.))},
         {-st / sqrt(2.), 0.}}};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        truth[p] += 1.e-7 * (1. - 2. * mass / radius) * exp(-2. * eta) *
                    tide[i][j] * mu[i] * mu[j];
        psi.get(4)[p] += 1.e-7 * (1. - 2. * mass / radius) * exp(2. * eta) *
                         tide[i][j] * conj(mu[i]) * conj(mu[j]);
      }
    }
  }
  const auto registration = register_frame(psi, rotation, mass);
  const auto evaluation = evaluate_geometric_second_order(
      registration, Scalar<DataVector>{DataVector(np, eta)}, rotation, metric,
      labels, weights, mass, 8);
  double error = 0., norm = 0.;
  for (size_t p = 0; p < np; ++p) {
    error += weights[p] *
             std::norm(get(evaluation.second_order.psi0_target)[p] - truth[p]);
    norm += weights[p] * std::norm(truth[p]);
  }
  CHECK(sqrt(error / norm) < 1.e-10);
  CHECK(evaluation.maximum_dyad_error < 1.e-12);
  const auto imposed = evaluate_geometric_second_order(
      registration, Scalar<DataVector>{DataVector(np, eta)}, rotation, metric,
      labels, weights, mass, 8, evaluation.second_order.fit.components);
  CHECK_ITERABLE_APPROX(get(imposed.second_order.psi0_target),
                        get(evaluation.second_order.psi0_target));
}
}  // namespace
}  // namespace gr::np

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.GeneralRelativity.NP.GeometricTide",
                  "[Unit][PointwiseFunctions]") {
  for (const double tilt : {0., .31}) {
    gr::np::test_round_map(tilt, 1.);
    gr::np::test_round_map(tilt, 7.3);
    gr::np::test_prescribed_tide(tilt);
  }
}
