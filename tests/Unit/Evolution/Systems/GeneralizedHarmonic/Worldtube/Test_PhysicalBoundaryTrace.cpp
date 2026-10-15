// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <random>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubePhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/AveragedUpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryTrace.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Tensor = tnsr::aa<DataVector, 3>;
using Phi = tnsr::iaa<DataVector, 3>;
using Dt = gh::TimeDerivative<gh::Solutions::all_solutions<3>, 3>;
using Temp = Variables<Dt::temporary_tags>;
using Solution = gh::Solutions::WrappedGr<gr::Solutions::KerrSchild>;

struct NativeRhs {
  explicit NativeRhs(size_t n) : g(n), pi(n), phi(n), temp(n) {}
  Tensor g, pi;
  Phi phi;
  Temp temp;
};

NativeRhs volume_rhs(const Tensor& g, const Tensor& pi, const Phi& phi,
                     const Phi& dg, const Phi& dpi,
                     const tnsr::ijaa<DataVector, 3>& dphi,
                     const tnsr::I<DataVector, 3>& velocity,
                     const Scalar<DataVector>& gamma2) {
  const size_t n = get<0, 0>(g).size();
  NativeRhs result(n);
  const Scalar<DataVector> gamma0(n, 0.7), gamma1(n, -1.);
  const Mesh<3> mesh(2, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto);
  const tnsr::I<DataVector, 3> coords(n, 0.);
  const InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>
      inv_jac(n, 0.);
  tmpl::as_pack<Dt::temporary_tags>(
      [&]<typename... Tags>(tmpl::type_<Tags>...) {
        Dt::apply(make_not_null(&result.g), make_not_null(&result.pi),
                  make_not_null(&result.phi),
                  make_not_null(&get<Tags>(result.temp))..., dg, dpi, dphi, g,
                  pi, phi, gamma0, gamma1, gamma2, gh::gauges::Harmonic{}, mesh,
                  0., coords, inv_jac, std::optional{velocity});
      });
  // The generic DG volume action adds this advection separately.
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t i = 0; i < 3; ++i) {
        result.g.get(a, b) += velocity.get(i) * dg.get(i, a, b);
        result.pi.get(a, b) += velocity.get(i) * dpi.get(i, a, b);
        for (size_t j = 0; j < 3; ++j) {
          result.phi.get(j, a, b) += velocity.get(i) * dphi.get(i, j, a, b);
        }
      }
    }
  }
  return result;
}

Tensor curvature(const Tensor& g, const Tensor& pi, const Phi& phi,
                 const Phi& dpi, const tnsr::ijaa<DataVector, 3>& dphi,
                 const tnsr::i<DataVector, 3>& sc,
                 const tnsr::I<DataVector, 3>& sv, const Temp& temp,
                 const gh::worldtube::PhysicalScreen& screen) {
  const size_t n = get<0, 0>(g).size();
  const auto& normal =
      get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(temp);
  const auto K = gh::extrinsic_curvature(normal, pi, phi);
  Tensor result(n, 0.);
  gh::BoundaryConditions::Bjorhus::detail::add_physical_terms_to_dt_v_minus(
      make_not_null(&result), Scalar<DataVector>(n, 0.), sc, sv, normal,
      screen.lower, screen.mixed, screen.upper,
      get<gr::Tags::InverseSpatialMetric<DataVector, 3>>(temp), K, g,
      get<gr::Tags::InverseSpacetimeMetric<DataVector, 3>>(temp), Phi(n, 0.),
      Tensor(n, 0.), phi, dphi, dpi,
      std::array<DataVector, 4>{DataVector(n, 0.), DataVector(n, 0.),
                                DataVector(n, 0.), DataVector(n, 1.)});
  return result;
}

template <typename T>
void randomize(T& tensor, std::mt19937& gen) {
  std::uniform_real_distribution<double> distribution(-0.2, 0.2);
  for (auto& component : tensor) {
    for (auto& value : component) {
      value = distribution(gen);
    }
  }
}

Tensor characteristic(const Tensor& pi, const Phi& phi,
                      const tnsr::I<DataVector, 3>& normal, double sign) {
  auto result = pi;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t i = 0; i < 3; ++i) {
        result.get(a, b) += sign * normal.get(i) * phi.get(i, a, b);
      }
    }
  }
  return result;
}

void test_native() {
  constexpr size_t n = 8;
  std::mt19937 gen(271828);
  tnsr::I<DataVector, 3> coords(n);
  for (size_t p = 0; p < n; ++p) {
    get<0>(coords)[p] = 3. + 0.1 * p;
    get<1>(coords)[p] = 0.3 + 0.2 * p;
    get<2>(coords)[p] = -0.7 + 0.1 * p;
  }
  const Solution solution(1., std::array<double, 3>{0.1, 0.2, 0.3},
                          std::array<double, 3>{0., 0., 0.});
  const auto exact = solution.variables(
      coords, 0.,
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
                 gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>{});
  const auto& g = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(exact);
  Tensor pi(n), target(n);
  Phi phi(n), dg(n), dpi(n);
  tnsr::ijaa<DataVector, 3> dphi(n);
  tnsr::I<DataVector, 3> velocity(n), sv(n), dt_sv(n, 0.);
  tnsr::i<DataVector, 3> sc(n);
  Scalar<DataVector> gamma2(n, 0.5), gamma1(n, -1.);
  const Tensor zero(n, 0.);
  const std::array<tnsr::AA<DataVector, 3>, 2> dt_E{
      tnsr::AA<DataVector, 3>(n, 0.), tnsr::AA<DataVector, 3>(n, 0.)};
  for (size_t trial = 0; trial < 8; ++trial) {
    CAPTURE(trial);
    randomize(pi, gen);
    randomize(phi, gen);
    randomize(target, gen);
    randomize(dg, gen);
    randomize(dpi, gen);
    randomize(dphi, gen);
    randomize(velocity, gen);
    randomize(sc, gen);
    const auto full = volume_rhs(g, pi, phi, dg, dpi, dphi, velocity, gamma2);
    const auto& hi =
        get<gr::Tags::InverseSpatialMetric<DataVector, 3>>(full.temp);
    DataVector norm(n, 0.);
    for (size_t i = 0; i < 3; ++i) {
      sv.get(i) = 0.;
      for (size_t j = 0; j < 3; ++j) {
        sv.get(i) += hi.get(i, j) * sc.get(j);
      }
      norm += sc.get(i) * sv.get(i);
    }
    for (size_t i = 0; i < 3; ++i) {
      sc.get(i) /= sqrt(norm);
      sv.get(i) /= sqrt(norm);
    }
    const gh::worldtube::PhysicalScreen screen(
        g, get<gr::Tags::InverseSpacetimeMetric<DataVector, 3>>(full.temp),
        get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(full.temp), sv);

    // Two independent screen probes. Orthonormalization is unnecessary for
    // the cancellation identity (it holds for any screen STF contraction).
    std::array<tnsr::AA<DataVector, 3>, 2> E{tnsr::AA<DataVector, 3>(n, 0.),
                                             tnsr::AA<DataVector, 3>(n, 0.)};
    const auto& gi =
        get<gr::Tags::InverseSpacetimeMetric<DataVector, 3>>(full.temp);
    for (auto& e : E) {
      Tensor seed(n);
      randomize(seed, gen);
      const auto cov = gh::worldtube::physical_projection(seed, screen);
      for (size_t a = 0; a < 4; ++a)
        for (size_t b = a; b < 4; ++b)
          for (size_t c = 0; c < 4; ++c)
            for (size_t d = 0; d < 4; ++d)
              e.get(a, b) += gi.get(a, c) * gi.get(b, d) * cov.get(c, d);
    }
    auto tg = dg, tp = dpi;
    auto tf = dphi;
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b)
        for (size_t i = 0; i < 3; ++i)
          for (size_t j = 0; j < 3; ++j) {
            tg.get(i, a, b) -= sc.get(i) * sv.get(j) * dg.get(j, a, b);
            tp.get(i, a, b) -= sc.get(i) * sv.get(j) * dpi.get(j, a, b);
            for (size_t k = 0; k < 3; ++k)
              tf.get(i, k, a, b) -=
                  sc.get(i) * sv.get(j) * dphi.get(j, k, a, b);
          }
    const auto tangential =
        volume_rhs(g, pi, phi, tg, tp, tf, velocity, gamma2);
    const auto W = curvature(g, pi, phi, dpi, dphi, sc, sv, full.temp, screen);
    const auto Wt = curvature(g, pi, phi, tp, tf, sc, sv, full.temp, screen);
    const auto& lapse = get<gr::Tags::Lapse<DataVector>>(full.temp);
    const auto& shift = get<gr::Tags::Shift<DataVector, 3>>(full.temp);
    DataVector lambda = -get(lapse);
    for (size_t i = 0; i < 3; ++i)
      lambda -= sc.get(i) * (shift.get(i) + velocity.get(i));
    const auto surface = [&](const NativeRhs& rhs, const Tensor& weyl) {
      return gh::worldtube::physical_trace_rhs(
          pi, phi, rhs.pi, rhs.phi,
          get<gh::Tags::ThreeIndexConstraint<DataVector, 3>>(rhs.temp), weyl,
          target, gamma2, lambda, sv, dt_sv, E, dt_E);
    };
    CHECK_ITERABLE_APPROX(surface(full, W), surface(tangential, Wt));

    // Arbitrary nonzero target, curved metric, oblique normal: all other
    // characteristic sectors must be preserved exactly up to roundoff.
    Tensor ghost_pi(n);
    Phi ghost_phi(n);
    gh::worldtube::physical_ghost(make_not_null(&ghost_pi),
                                  make_not_null(&ghost_phi), pi, phi, target,
                                  sc, sv, screen);
    const auto old_minus = characteristic(pi, phi, sv, -1.);
    const auto new_minus = characteristic(ghost_pi, ghost_phi, sv, -1.);
    CHECK_ITERABLE_APPROX(gh::worldtube::physical_projection(new_minus, screen),
                          gh::worldtube::physical_projection(target, screen));
    CHECK_ITERABLE_APPROX(characteristic(ghost_pi, ghost_phi, sv, 1.),
                          characteristic(pi, phi, sv, 1.));
    Tensor jump = new_minus;
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b)
        jump.get(a, b) -= old_minus.get(a, b);
    CHECK_ITERABLE_APPROX(jump,
                          gh::worldtube::physical_projection(jump, screen));
    auto tangent_jump = ghost_phi;
    for (size_t i = 0; i < 3; ++i)
      for (size_t a = 0; a < 4; ++a)
        for (size_t b = a; b < 4; ++b) {
          tangent_jump.get(i, a, b) -= phi.get(i, a, b);
          tangent_jump.get(i, a, b) += 0.5 * sc.get(i) * jump.get(a, b);
        }
    CHECK_ITERABLE_APPROX(tangent_jump, Phi(n, 0.));

    // Exercise the actual numerical flux with opposite ghost-side normals.
    auto opposite = sc;
    for (auto& c : opposite)
      c *= -1.;
    Tensor flux_g(n), flux_pi(n);
    Phi flux_phi(n);
    gh::BoundaryCorrections::AveragedUpwindPenalty<3>{}.dg_boundary_terms(
        make_not_null(&flux_g), make_not_null(&flux_pi),
        make_not_null(&flux_phi), g, pi, phi, gamma1, gamma2, sc, velocity, g,
        ghost_pi, ghost_phi, gamma1, gamma2, opposite, velocity,
        dg::Formulation::StrongInertial);
    CHECK_ITERABLE_APPROX(flux_g, zero);
    CHECK_ITERABLE_APPROX(characteristic(flux_pi, flux_phi, sv, 1.), zero);
    auto expected_flux = jump;
    for (auto& c : expected_flux)
      c *= lambda;
    CHECK_ITERABLE_APPROX(characteristic(flux_pi, flux_phi, sv, -1.),
                          expected_flux);

    // Check the public boundary interface, including removal of the Frozen
    // physical RHS, preservation of CP/gauge, cloning, and serialization.
    using BC = gh::BoundaryConditions::WorldtubePhysicalGhost;
    const BC bc(std::make_unique<Solution>(solution),
                gh::BoundaryConditions::detail::SectorImposition::Frozen);
    const auto restored = serialize_and_deserialize(bc);
    REQUIRE(restored.get_clone() != nullptr);
    const Domain<3> domain{};
    const Element<3> element(ElementId<3>(0), {});
    const domain::FunctionsOfTimeMap functions{};
    const gh::worldtube::KretschmannFaceData<3> face{};
    const ConstraintDamping::Constant<3, Frame::Grid> damping(0.5);
    const auto common = [&](const auto& boundary, Tensor* cg, Tensor* cp,
                            Phi* cf) {
      return boundary.dg_time_derivative(
          make_not_null(cg), make_not_null(cp), make_not_null(cf),
          std::optional{velocity}, sc, sv, g, pi, phi, coords, gamma1, gamma2,
          lapse, shift, gi,
          get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(full.temp),
          get<gh::Tags::ThreeIndexConstraint<DataVector, 3>>(full.temp),
          get<gh::Tags::GaugeH<DataVector, 3>>(full.temp),
          get<gh::Tags::SpacetimeDerivGaugeH<DataVector, 3>>(full.temp), full.g,
          full.pi, full.phi, dg, dpi, dphi, 0., domain, element, functions,
          face, damping);
    };
    REQUIRE_FALSE(common(restored, &flux_g, &flux_pi, &flux_phi).has_value());
    auto correction_minus = characteristic(flux_pi, flux_phi, sv, -1.);
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b)
        correction_minus.get(a, b) -= get(gamma2) * flux_g.get(a, b);
    CHECK_ITERABLE_APPROX(
        gh::worldtube::physical_projection(correction_minus, screen), zero);
    using Sector = gh::BoundaryConditions::detail::SectorImposition;
    const gh::BoundaryConditions::WorldtubeTypeD<3> baseline(
        Sector::Bjorhus, Sector::Frozen, Sector::Frozen,
        gh::BoundaryConditions::detail::PhysicalModel::None);
    Tensor base_g(n), base_pi(n);
    Phi base_phi(n);
    REQUIRE_FALSE(common(baseline, &base_g, &base_pi, &base_phi).has_value());
    CHECK_ITERABLE_APPROX(flux_g, base_g);
    auto base_minus = characteristic(base_pi, base_phi, sv, -1.);
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b)
        base_minus.get(a, b) -= get(gamma2) * base_g.get(a, b);
    const auto base_physical =
        gh::worldtube::physical_projection(base_minus, screen);
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b)
        base_minus.get(a, b) -= base_physical.get(a, b);
    CHECK_ITERABLE_APPROX(correction_minus, base_minus);
    CHECK_ITERABLE_APPROX(characteristic(flux_pi, flux_phi, sv, 1.),
                          characteristic(base_pi, base_phi, sv, 1.));
    Scalar<DataVector> ghost_gamma1(n), ghost_gamma2(n), ghost_lapse(n);
    tnsr::I<DataVector, 3> ghost_shift(n);
    tnsr::II<DataVector, 3> ghost_hi(n);
    REQUIRE_FALSE(
        restored
            .dg_ghost(
                make_not_null(&flux_g), make_not_null(&ghost_pi),
                make_not_null(&ghost_phi), make_not_null(&ghost_gamma1),
                make_not_null(&ghost_gamma2), make_not_null(&ghost_lapse),
                make_not_null(&ghost_shift), make_not_null(&ghost_hi),
                std::optional{velocity}, sc, sv, g, pi, phi, coords, gamma1,
                gamma2, lapse, shift, gi,
                get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(full.temp),
                get<gh::Tags::ThreeIndexConstraint<DataVector, 3>>(full.temp),
                get<gh::Tags::GaugeH<DataVector, 3>>(full.temp),
                get<gh::Tags::SpacetimeDerivGaugeH<DataVector, 3>>(full.temp),
                full.g, full.pi, full.phi, dg, dpi, dphi, 0., domain, element,
                functions, face, damping)
            .has_value());
    CHECK_ITERABLE_APPROX(flux_g, g);
    CHECK_ITERABLE_APPROX(ghost_gamma1, gamma1);
    CHECK_ITERABLE_APPROX(ghost_gamma2, gamma2);
    CHECK_ITERABLE_APPROX(ghost_lapse, lapse);
    CHECK_ITERABLE_APPROX(ghost_shift, shift);
    CHECK_ITERABLE_APPROX(ghost_hi, hi);
    const auto exact_target =
        characteristic(get<gh::Tags::Pi<DataVector, 3>>(exact),
                       get<gh::Tags::Phi<DataVector, 3>>(exact), sv, -1.);
    CHECK_ITERABLE_APPROX(
        gh::worldtube::physical_projection(
            characteristic(ghost_pi, ghost_phi, sv, -1.), screen),
        gh::worldtube::physical_projection(exact_target, screen));
    CHECK_ITERABLE_APPROX(characteristic(ghost_pi, ghost_phi, sv, 1.),
                          characteristic(pi, phi, sv, 1.));
    // Reject a spacelike tube and unsupported gamma1 before touching data.
    get(gamma1) = 0.;
    REQUIRE(common(restored, &flux_g, &flux_pi, &flux_phi).has_value());
    get(gamma1) = -1.;
    for (size_t i = 0; i < 3; ++i)
      velocity.get(i) = 10. * sv.get(i);
    REQUIRE(common(restored, &flux_g, &flux_pi, &flux_phi).has_value());
  }
}
// Two-polarization planar wave fixes the curvature sign and factor of two.
void test_plane_wave() {
  constexpr size_t n = 8;
  Tensor g(n, 0.), pi(n, 0.);
  Phi phi(n, 0.), dpi(n, 0.);
  tnsr::ijaa<DataVector, 3> dphi(n, 0.);
  get<0, 0>(g) = -1.;
  for (size_t i = 1; i < 4; ++i)
    g.get(i, i) = 1.;
  tnsr::i<DataVector, 3> sc(n, 0.);
  get<0>(sc) = 1.;
  tnsr::I<DataVector, 3> sv(n, 0.), velocity(n, 0.);
  get<0>(sv) = 1.;
  std::array<tnsr::AA<DataVector, 3>, 2> E{tnsr::AA<DataVector, 3>(n, 0.),
                                           tnsr::AA<DataVector, 3>(n, 0.)};
  get<2, 2>(E[0]) = 1. / sqrt(2.);
  get<3, 3>(E[0]) = -1. / sqrt(2.);
  get<2, 3>(E[1]) = 1. / sqrt(2.);
  const std::array<double, 2> second{0.7, -0.4};
  for (size_t a = 0; a < 4; ++a)
    for (size_t b = a; b < 4; ++b)
      for (size_t p = 0; p < 2; ++p) {
        dpi.get(0, a, b) -= second[p] * E[p].get(a, b);
        dphi.get(0, 0, a, b) += second[p] * E[p].get(a, b);
      }
  const auto rhs = volume_rhs(g, pi, phi, phi, dpi, dphi, velocity,
                              Scalar<DataVector>(n, 0.));
  const gh::worldtube::PhysicalScreen screen(
      g, get<gr::Tags::InverseSpacetimeMetric<DataVector, 3>>(rhs.temp),
      get<gr::Tags::SpacetimeNormalVector<DataVector, 3>>(rhs.temp), sv);
  const auto W = curvature(g, pi, phi, dpi, dphi, sc, sv, rhs.temp, screen);
  const auto rate = gh::worldtube::physical_trace_rhs(
      pi, phi, pi, phi, phi, pi, W, Scalar<DataVector>(n, 0.),
      DataVector(n, -1.), sv, velocity, E,
      {tnsr::AA<DataVector, 3>(n, 0.), tnsr::AA<DataVector, 3>(n, 0.)});
  for (size_t p = 0; p < 2; ++p)
    CHECK_ITERABLE_APPROX(rate[p], DataVector(n, -2. * second[p]));
}

void test_frame_transport() {
  constexpr size_t n = 8;
  std::mt19937 gen(314159);
  Tensor pi(n), rhs_pi(n);
  Phi phi(n), rhs_phi(n);
  randomize(pi, gen);
  randomize(phi, gen);
  randomize(rhs_pi, gen);
  randomize(rhs_phi, gen);
  constexpr double omega = 0.7;
  const auto frame = [=](double t) {
    std::array<tnsr::AA<DataVector, 3>, 2> E{tnsr::AA<DataVector, 3>(n, 0.),
                                             tnsr::AA<DataVector, 3>(n, 0.)};
    const std::array<double, 4> e{0., -sin(omega * t), cos(omega * t), 0.};
    const std::array<double, 4> f{0., 0., 0., 1.};
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b) {
        E[0].get(a, b) = (e[a] * e[b] - f[a] * f[b]) / sqrt(2.);
        E[1].get(a, b) = (e[a] * f[b] + f[a] * e[b]) / sqrt(2.);
      }
    return E;
  };
  const auto E = frame(0.);
  auto dt_E = E;
  for (auto& e : dt_E)
    for (auto& c : e)
      c = 0.;
  get<1, 2>(dt_E[0]) = -omega / sqrt(2.);
  get<1, 3>(dt_E[1]) = -omega / sqrt(2.);
  tnsr::I<DataVector, 3> sv(n, 0.), dt_sv(n, 0.);
  get<0>(sv) = 1.;
  get<1>(dt_sv) = omega;
  const auto rate = gh::worldtube::physical_trace_rhs(
      pi, phi, rhs_pi, rhs_phi, Phi(n, 0.), Tensor(n, 0.), Tensor(n, 0.),
      Scalar<DataVector>(n, 0.), DataVector(n, -1.), sv, dt_sv, E, dt_E);
  const auto amplitude = [&](double t) {
    tnsr::I<DataVector, 3> normal(n, 0.);
    get<0>(normal) = cos(omega * t);
    get<1>(normal) = sin(omega * t);
    auto pt = pi;
    auto ft = phi;
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = a; b < 4; ++b) {
        pt.get(a, b) += t * rhs_pi.get(a, b);
        for (size_t i = 0; i < 3; ++i)
          ft.get(i, a, b) += t * rhs_phi.get(i, a, b);
      }
    const auto v = characteristic(pt, ft, normal, -1.);
    const auto et = frame(t);
    std::array<DataVector, 2> result{DataVector(n, 0.), DataVector(n, 0.)};
    for (size_t p = 0; p < 2; ++p)
      for (size_t a = 0; a < 4; ++a)
        for (size_t b = 0; b < 4; ++b)
          result[p] += et[p].get(a, b) * v.get(a, b);
    return result;
  };
  double previous_error = 0.;
  for (const double h : {1.e-3, 5.e-4, 2.5e-4}) {
    const auto plus = amplitude(h), minus = amplitude(-h);
    double error = 0.;
    for (size_t p = 0; p < 2; ++p)
      error =
          std::max(error, max(abs((plus[p] - minus[p]) / (2. * h) - rate[p])));
    CHECK(error < 1.e-6);
    if (previous_error > 0.)
      CHECK(previous_error / error > 3.9);
    previous_error = error;
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Gh.PhysicalBoundaryTrace",
                  "[Unit][Evolution][PhysicalTrace]") {
  register_classes_with_charm(tmpl::list<Solution>{});
  test_native();
  test_plane_wave();
  test_frame_transport();
}
