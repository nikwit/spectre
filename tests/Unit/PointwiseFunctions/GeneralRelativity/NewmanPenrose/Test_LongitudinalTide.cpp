// Distributed under the MIT License.
// See LICENSE.txt for details.
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Framework/TestingFramework.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPoints.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/LongitudinalTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
namespace gr::np {
namespace {
using C = std::complex<double>;
void seed_test() {
  const size_t n = 4;
  Scalar<ComplexDataVector> c(n, C{-.03, .0001}), a(n), b(n);
  get(a) = ComplexDataVector{C{0., 0.}, C{.04, .02}, C{0., 0.}, C{.06, -.03}};
  get(b) = ComplexDataVector{C{0., 0.}, C{0., 0.}, C{-.05, .03}, C{-.03, .02}};
  auto p = type_d_scalars(c, a, b);
  const auto seed = longitudinal_seed(p);
  CHECK_ITERABLE_APPROX(seed.psi, p);
  CHECK_ITERABLE_APPROX(seed.coulomb, c);
  p.get(0) = p.get(4) = C{std::numeric_limits<double>::quiet_NaN(), 0.};
  CHECK(longitudinal_seed(p).psi == seed.psi);
  WeylScalars column(n, C{0., 0.});
  for (size_t k = 0; k < 5; ++k)
    for (size_t j = 0; j < n; ++j)
      column.get(k)[j] = C{.12 + .01 * k, -.07 + .02 * j};
  const auto ds = longitudinal_seed_response(seed, column);
  const double h = 1.e-7;
  auto plus = seed.psi, minus = seed.psi;
  for (size_t k = 0; k < 5; ++k) {
    plus.get(k) += h * column.get(k);
    minus.get(k) -= h * column.get(k);
  }
  const auto sp = longitudinal_seed(plus), sm = longitudinal_seed(minus);
  const auto tol = Approx::custom().epsilon(2.e-9).scale(1.);
  for (size_t k = 0; k < 5; ++k)
    CHECK_ITERABLE_CUSTOM_APPROX(
        ds.get(k),
        ComplexDataVector{(sp.psi.get(k) - sm.psi.get(k)) / (2. * h)}, tol);
  CHECK_THROWS_WITH(
      longitudinal_seed(WeylScalars(n, C{})),
      Catch::Matchers::ContainsSubstring("Singular longitudinal"));
}
std::pair<TriadVector, DataVector> sphere() {
  const size_t nt = 9, nf = 17, n = nt * nf;
  const auto& z = Spectral::collocation_points<Spectral::Basis::Legendre,
                                               Spectral::Quadrature::Gauss>(nt);
  const auto& w = Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                               Spectral::Quadrature::Gauss>(nt);
  TriadVector dir(n, 0.);
  DataVector weights(n);
  for (size_t j = 0; j < nf; ++j)
    for (size_t i = 0; i < nt; ++i) {
      const size_t p = i + nt * j;
      const double phi = 2. * acos(-1.) * j / nf;
      dir.get(0)[p] = sqrt(1. - z[i] * z[i]) * cos(phi);
      dir.get(1)[p] = sqrt(1. - z[i] * z[i]) * sin(phi);
      dir.get(2)[p] = z[i];
      weights[p] = w[i] / (2. * nf);
    }
  return {dir, weights};
}
double full_test(const double amplitude, const double sigma, const double eta,
                 const double radius) {
  const auto [labels, weights] = sphere();
  const size_t n = weights.size();
  const double mass = 1., c = mass / pow(radius, 3),
               root_f = sqrt(1. - 2. * mass / radius);
  auto direction = labels;
  for (auto& v : direction)
    v *= sigma;
  tnsr::ii<DataVector, 3> metric(n, 0.);
  for (size_t i = 0; i < 3; ++i)
    metric.get(i, i) = 1.;
  const auto rotation = adapted_triad(metric, direction);
  auto e =
      stf_from_components(std::array<double, 5>{{.4, -.3, .12, .08, -.11}});
  auto b =
      stf_from_components(std::array<double, 5>{{-.2, .1, -.09, .13, .05}});
  for (auto& v : e)
    v *= amplitude * c;
  for (auto& v : b)
    v *= amplitude * c;
  auto q = quadrupole_tide_tensor(e, b, labels, Scalar<DataVector>(n, radius),
                                  mass, false);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      q.get(i, j) +=
          c * ((i == j ? 1. : 0.) - 3. * labels.get(i) * labels.get(j));
  auto psi =
      type_iii(weyl_scalars_from_tidal_tensor(rotate_symmetric(q, rotation)),
               Scalar<DataVector>(n, -eta), Scalar<DataVector>(n, 0.));
  Scalar<DataVector> lapse(n, .83), dt(n, 0.);
  TriadVector shift(n, 0.);
  tnsr::i<DataVector, 3> dk(n, 0.);
  for (size_t i = 0; i < 3; ++i)
    shift.get(i) = std::array<double, 3>{{.12, -.04, .06}}[i];
  for (size_t p = 0; p < n; ++p) {
    double enn = 0.;
    std::array<double, 3> en{};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j) {
        en[i] += e.get(i, j) * labels.get(j)[p];
        enn += e.get(i, j) * labels.get(i)[p] * labels.get(j)[p];
      }
    const double radial = (-288. * c * c + 144. * c * enn) * root_f / radius;
    get(dt)[p] = get(lapse)[p] * (-sinh(eta) * sigma * radial);
    for (size_t i = 0; i < 3; ++i) {
      dk.get(i)[p] = cosh(eta) * radial * labels.get(i)[p] -
                     96. * c / radius * (en[i] - enn * labels.get(i)[p]);
      get(dt)[p] += shift.get(i)[p] * dk.get(i)[p];
    }
  }
  const auto fit = evaluate_longitudinal_tide(
      psi, rotation, metric, labels, weights, mass, 8, lapse, shift, dk, dt);
  CHECK(std::isfinite(fit.relative_residual));
  CHECK(fit.normalized_condition_number < 4.);
  for (size_t k = 1; k <= 3; ++k)
    CHECK(fit.completed.get(k) == psi.get(k));
  double error = 0.;
  for (size_t p = 0; p < n; ++p)
    error += weights[p] * std::norm(fit.completed.get(0)[p] - psi.get(0)[p]);
  error = sqrt(error) / c;
  if (amplitude == .001) {
    auto poisoned = psi;
    poisoned.get(0) = poisoned.get(4) =
        C{std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN()};
    const auto same =
        evaluate_longitudinal_tide(poisoned, rotation, metric, labels, weights,
                                   mass, 8, lapse, shift, dk, dt);
    CHECK(same.completed == fit.completed);
    CHECK(same.moments == fit.moments);
    // A point-dependent dyad spin must rotate the target with spin weight +2.
    Scalar<DataVector> chi(n);
    get(chi) = .13 * labels.get(0) - .09 * labels.get(1);
    auto spun_rotation = rotation;
    for (size_t i = 0; i < 3; ++i) {
      spun_rotation.get(1, i) = cos(get(chi)) * rotation.get(1, i) -
                                sin(get(chi)) * rotation.get(2, i);
      spun_rotation.get(2, i) = sin(get(chi)) * rotation.get(1, i) +
                                cos(get(chi)) * rotation.get(2, i);
    }
    const auto spun = evaluate_longitudinal_tide(
        type_iii(psi, Scalar<DataVector>(n, 0.), chi), spun_rotation, metric,
        labels, weights, mass, 8, lapse, shift, dk, dt);
    CHECK(max(abs(spun.completed.get(0) -
                  exp(C{0., 2.} * get(chi)) * fit.completed.get(0))) <
          1.e-11 * c);
    // Check the actual matching dispatcher and incoming tensor conversion.
    if (sigma < 0.) {
      auto inverse_rotation = rotation;
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          inverse_rotation.get(i, j) = rotation.get(j, i);
      const auto nrq = rotate_symmetric(tidal_tensor_from_weyl_scalars(psi),
                                        inverse_rotation);
      tnsr::ii<DataVector, 3> electric(n), magnetic(n);
      tnsr::i<DataVector, 3> normal(n);
      for (size_t i = 0; i < 3; ++i) {
        normal.get(i) = direction.get(i);
        for (size_t j = i; j < 3; ++j) {
          electric.get(i, j) = real(nrq.get(i, j));
          magnetic.get(i, j) = imag(nrq.get(i, j));
        }
      }
      gh::worldtube::KretschmannFaceData<3> face{};
      face.direction = Direction<3>::lower_xi();
      face.quadrature_weights = weights;
      face.d_kretschmann = dk;
      face.dt_kretschmann = dt;
      face.kretschmann = Scalar<DataVector>(n, 0.);
      const auto matched = gh::worldtube::evaluate_matching(
          gh::worldtube::PhysicalModel::QuadrupoleLongitudinal, mass, electric,
          magnetic, metric, normal, lapse, shift, &face, std::nullopt);
      CHECK(max(abs(get(matched.psi0_target) - fit.completed.get(0))) <
            1.e-12 * c);
      auto expected = incoming_weyl_field(
          Scalar<ComplexDataVector>{fit.completed.get(0)}, rotation);
      for (auto& v : expected)
        v *= .5;
      CHECK_ITERABLE_APPROX(matched.incoming_mode, expected);
    }
  }
  std::cout << "LONGITUDINAL_TEST amplitude=" << amplitude << " sigma=" << sigma
            << " eta=" << eta << " radius=" << radius << " error=" << error
            << " condition=" << fit.normalized_condition_number << '\n';
  return error;
}
}  // namespace
SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.Np.LongitudinalTide",
                  "[Unit][PointwiseFunctions]") {
  seed_test();
  for (const double sigma : {-1., 1.})
    for (const double eta : {0., .35})
      for (const double radius : {3., 5.}) {
        const double e1 = full_test(.001, sigma, eta, radius),
                     e2 = full_test(.0005, sigma, eta, radius),
                     e3 = full_test(.00025, sigma, eta, radius);
        CHECK(e1 < 2.e-6);
        CHECK(e2 / e1 < .30);
        CHECK(e2 / e1 > .20);
        CHECK(e3 / e2 < .30);
        CHECK(e3 / e2 > .20);
      }
}
}  // namespace gr::np
