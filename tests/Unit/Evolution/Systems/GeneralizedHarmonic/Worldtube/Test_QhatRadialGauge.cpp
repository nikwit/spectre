// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <cmath>

#include "DataStructures/Tensor/Slice.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/QhatRadialGauge.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/TestHelpers.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace {
void test_profile() {
  const auto setup =
      TestHelpers::gh_worldtube::shell_element(10, 5, {{0., 0., 0.}}, 2.5, 3.);
  const size_t n = setup.mesh.number_of_grid_points();
  tnsr::aa<DataVector, 3> metric(n, 0.), pi(n, 0.);
  tnsr::iaa<DataVector, 3> phi(n, 0.);
  get<0, 0>(metric) = -1.;
  for (size_t i = 1; i < 4; ++i) {
    metric.get(i, i) = 1.;
  }
  DataVector radius_squared(n, 0.);
  for (const auto& x : setup.inertial_coords) {
    radius_squared += square(x);
  }
  get<0, 0>(pi) = .03 * radius_squared;
  get<0, 1>(pi) = .01 * radius_squared;
  std::optional<gh::worldtube::QhatRadialGaugeData> initial{};
  const auto initialize = [&](const double time) {
    gh::worldtube::initialize_qhat_radial_gauge(
        make_not_null(&initial), metric, pi, phi, setup.inertial_coords,
        tnsr::I<double, 3>(0.), setup.mesh, setup.inverse_jacobian,
        Direction<3>::lower_xi(), time, {});
  };
  initialize(0.);
  REQUIRE(initial.has_value());
  gh::worldtube::KretschmannFaceData<3> holder{};
  holder.qhat_radial_gauge = initial;
  CHECK(serialize_and_deserialize(holder) == holder);
  const auto captured = *initial;
  get<0, 0>(pi) *= 2.;
  initialize(1.);
  initialize(0.);  // Self-start resets must not recapture.
  CHECK(*initial == captured);

  // Rotate the sphere and vector data by 90 degrees about z, and translate it.
  // The profile must rotate with the mesh without using a future donor trace.
  domain::FunctionsOfTimeMap functions;
  functions["Rotation"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.,
          std::array<DataVector, 3>{{DataVector{sqrt(.5), 0., 0., sqrt(.5)},
                                     DataVector(4, 0.), DataVector(4, 0.)}},
          10.);
  const tnsr::I<double, 3> center{{.4, -.2, .1}};
  auto xyz = data_on_slice(setup.inertial_coords, setup.mesh.extents(), 0, 0);
  const auto old_xyz = xyz;
  get<0>(xyz) = -get<1>(old_xyz) + get<0>(center);
  get<1>(xyz) = get<0>(old_xyz) + get<1>(center);
  get<2>(xyz) += get<2>(center);
  const size_t nf = get<0>(xyz).size();
  auto g = data_on_slice(metric, setup.mesh.extents(), 0, 0);
  tnsr::aa<DataVector, 3> p(nf, 0.), metric_correction(nf, 0.);
  tnsr::iaa<DataVector, 3> f(nf, 0.), dg(nf, 0.), dp(nf, 0.);
  tnsr::ijaa<DataVector, 3> df(nf, 0.);
  tnsr::I<DataVector, 3> radial(nf), shift(nf, 0.);
  std::optional<tnsr::I<DataVector, 3>> velocity{
      tnsr::I<DataVector, 3>(nf, 0.)};
  get<0>(*velocity) = .1;
  for (size_t i = 0; i < 3; ++i) {
    radial.get(i) = (xyz.get(i) - center.get(i)) / 2.5;
    dp.get(i, 0, 0) = .06 * 2.5 * radial.get(i);
    dp.get(i, 0, 2) = .02 * 2.5 * radial.get(i);
  }
  get<0, 0>(p) = .03 * 2.5 * 2.5;
  get<0, 2>(p) = .01 * 2.5 * 2.5;
  const Scalar<DataVector> lapse(nf, 1.), gamma2(nf, 37.);
  const std::array<double, 3> coefficients{{-2.65, -2.33, -2.1}};
  const auto correction = [&]() {
    return gh::worldtube::qhat_radial_gauge_correction(
        *initial, coefficients, g, p, f, dg, dp, df, xyz, center, lapse, shift,
        velocity, gamma2, metric_correction, 1., functions);
  };
  const auto tolerance = Approx::custom().epsilon(1.e-10).scale(1.);
  const tnsr::a<DataVector, 3> zero(nf, 0.);
  CHECK_ITERABLE_CUSTOM_APPROX(correction(), zero, tolerance);
  const auto q0 =
      gh::worldtube::qhat_radial_data(g, p, f, dg, dp, df, radial).first;
  for (auto& x : p) {
    x *= 1.2;
  }
  for (auto& x : dp) {
    x *= 1.2;
  }
  const auto actual = correction();
  DataVector delta_r(nf, 0.);
  for (size_t i = 0; i < 3; ++i) {
    delta_r += .2 * radial.get(i) * q0.get(i + 1);
  }
  const DataVector speed = 1. - .1 * get<0>(radial);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get<0>(actual),
      DataVector{speed * (2. - coefficients[0]) / 2.5 * .2 * get<0>(q0)},
      tolerance);
  for (size_t i = 0; i < 3; ++i) {
    CHECK_ITERABLE_CUSTOM_APPROX(
        actual.get(i + 1),
        DataVector{
            speed / 2.5 *
            ((2. - coefficients[2]) * .2 * q0.get(i + 1) +
             (coefficients[2] - coefficients[1]) * radial.get(i) * delta_r)},
        tolerance);
  }

  // Test analytic basis derivatives against directional finite differences
  // on a curved, off-diagonal metric with nonzero Phi.
  get<0, 0>(g) -= .2;
  get<0, 1>(g) = .03;
  get<1, 2>(g) = .02;
  for (size_t c = 0; c < f.size(); ++c) {
    f[c] = .0001 * static_cast<double>(c + 1);
  }
  tnsr::aa<DataVector, 3> drg(nf), drp(nf);
  tnsr::iaa<DataVector, 3> drf(nf);
  for (size_t c = 0; c < drg.size(); ++c) {
    drg[c] = .001 * static_cast<double>(c + 1);
    drp[c] = -.002 * static_cast<double>(c + 1);
  }
  for (size_t c = 0; c < drf.size(); ++c) {
    drf[c] = .0003 * static_cast<double>(c + 1);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        dg.get(i, a, b) = radial.get(i) * drg.get(a, b);
        dp.get(i, a, b) = radial.get(i) * drp.get(a, b);
        for (size_t j = 0; j < 3; ++j) {
          df.get(i, j, a, b) = radial.get(i) * drf.get(j, a, b);
        }
      }
    }
  }
  const auto evaluate = [&](const double eps, const bool metric_only) {
    auto gg = g;
    auto pp = p;
    auto ff = f;
    for (size_t c = 0; c < gg.size(); ++c) {
      gg[c] += eps * drg[c];
      if (not metric_only) {
        pp[c] += eps * drp[c];
      }
    }
    if (not metric_only) {
      for (size_t c = 0; c < ff.size(); ++c) {
        ff[c] += eps * drf[c];
      }
    }
    return gh::worldtube::qhat_radial_data(gg, pp, ff, dg, dp, df, radial)
        .first;
  };
  constexpr double eps = 1.e-5;
  const auto plus = evaluate(eps, false), minus = evaluate(-eps, false);
  const auto dr =
      gh::worldtube::qhat_radial_data(g, p, f, dg, dp, df, radial).second;
  for (size_t a = 0; a < 4; ++a) {
    CHECK_ITERABLE_CUSTOM_APPROX(
        dr.get(a), DataVector{(plus.get(a) - minus.get(a)) / (2. * eps)},
        Approx::custom().epsilon(1.e-8).scale(1.));
  }
  // Nonzero metric corrections from other characteristic sectors must not
  // alter the intended full-qhat evolution through their changing basis.
  const auto before = correction();
  metric_correction = drg;
  const auto after = correction();
  const auto mp = evaluate(eps, true), mm = evaluate(-eps, true);
  // qhat(g,Pi+gamma2*dg,Phi)-qhat(g,Pi,Phi) is gamma2 ell.dg.
  auto shifted_p = p;
  for (size_t c = 0; c < p.size(); ++c) {
    shifted_p[c] += 37. * drg[c];
  }
  const auto qp =
      gh::worldtube::qhat_radial_data(g, shifted_p, f, dg, dp, df, radial)
          .first;
  const auto q =
      gh::worldtube::qhat_radial_data(g, p, f, dg, dp, df, radial).first;
  for (size_t a = 0; a < 4; ++a) {
    CHECK_ITERABLE_CUSTOM_APPROX(
        DataVector{after.get(a) - before.get(a)},
        DataVector{-(mp.get(a) - mm.get(a)) / (2. * eps) - qp.get(a) +
                   q.get(a)},
        Approx::custom().epsilon(1.e-8).scale(1.));
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Gh.QhatRadialGauge",
                  "[Unit][Evolution]") {
  domain::FunctionsOfTime::register_derived_with_charm();
  test_profile();
}
