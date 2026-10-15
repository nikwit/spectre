// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <limits>

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeEvolvedPhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryEvolution.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Tensor = tnsr::aa<DataVector, 3>;
using Phi = tnsr::iaa<DataVector, 3>;
using Target =
    gh::BoundaryConditions::detail::PrescribedSchwarzschildQuadrupole;
using Pulse = gh::BoundaryConditions::detail::RecordedMomentPulse;

tnsr::I<DataVector, 3> coordinates() {
  tnsr::I<DataVector, 3> x(size_t{7});
  for (size_t p = 0; p < 7; ++p) {
    const double theta = 0.4 + 0.35 * static_cast<double>(p);
    const double phi = 0.21 + 0.77 * static_cast<double>(p);
    get<0>(x)[p] = 3.25 * sin(theta) * cos(phi);
    get<1>(x)[p] = 3.25 * sin(theta) * sin(phi);
    get<2>(x)[p] = 3.25 * cos(theta);
  }
  return x;
}

void test_prescribed_np() {
  const auto x = coordinates();
  const size_t n = get<0>(x).size();
  const std::array<double, 5> electric{0.001, -0.002, 0.0003, -0.0004, 0.0005};
  const std::array<double, 5> magnetic{-0.0004, 0.0002, 0.0007, 0.0006,
                                       -0.0003};
  const Target target(1., Pulse(2., 8., electric, magnetic));
  const auto restored = serialize_and_deserialize(target);
  const auto actual = target.incoming_mode(6., x);
  CHECK_ITERABLE_APPROX(actual, restored.incoming_mode(6., x));
  for (const double time : {-3., 2., 10., 12.}) {
    const auto zero = target.incoming_mode(time, x);
    for (const auto& component : zero) {
      CHECK(max(abs(component)) < 1.e-25);
    }
  }
  const auto quarter = target.incoming_mode(4., x);
  for (size_t c = 0; c < actual.size(); ++c) {
    CHECK_ITERABLE_APPROX(4. * quarter[c], actual[c]);
  }

  // Independent known frame: the inward radial leg in Cholesky components is
  // row zero of the adapted triad, and a static KS observer has
  // tanh(eta)=-2M/r. This bypasses evaluate_matching, its invariant boost and
  // frame recovery.
  const gr::Solutions::KerrSchild background(
      1., std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0., 0., 0.});
  const auto data = background.variables(
      x, 0., tmpl::list<gr::Tags::SpatialMetric<DataVector, 3>>{});
  const auto& metric = get<gr::Tags::SpatialMetric<DataVector, 3>>(data);
  tnsr::I<DataVector, 3> directions(n);
  for (size_t i = 0; i < 3; ++i) {
    directions.get(i) = -x.get(i) / 3.25;
  }
  const auto rotation = gr::np::adapted_triad(metric, directions);
  gr::np::TriadVector radial(n), velocity(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    radial.get(i) = rotation.get(0, i);
  }
  const auto columns = gr::np::direct_tide_scalar_columns(
      radial, velocity, Scalar<DataVector>(n, -atanh(2. / 3.25)),
      Scalar<DataVector>(n, 3.25), rotation, 1.);
  Scalar<ComplexDataVector> psi0(n, 0.);
  for (size_t a = 0; a < 5; ++a) {
    get(psi0) +=
        std::complex<double>{electric[a], magnetic[a]} * columns[a].get(0);
  }
  auto expected = gr::np::orthonormal_to_coordinate_covariant(
      gr::np::incoming_weyl_field(psi0, rotation),
      gr::np::cholesky_factor(metric));
  for (auto& component : expected) {
    component *= 0.5;
  }
  CHECK_ITERABLE_APPROX(actual, expected);
  CHECK(max(abs(get(psi0))) > 1.e-5);
  CHECK_THROWS_WITH(Target(0., Pulse{}),
                    Catch::Matchers::ContainsSubstring("positive finite mass"));
  CHECK_THROWS_WITH(target.incoming_mode(6., tnsr::I<DataVector, 3>(n, 0.5)),
                    Catch::Matchers::ContainsSubstring("outside r=2M"));
}

void test_receiver_conversion() {
  const auto x = coordinates();
  const size_t n = get<0>(x).size();
  // Spin gives the receiver a tangential shift and a different screen metric
  // from the prescribed Schwarzschild target. This exercises temporal lift,
  // screen detracing, both polarizations, and nonzero normal components.
  const gh::Solutions::WrappedGr<gr::Solutions::KerrSchild> receiver(
      0.9, std::array<double, 3>{0.2, -0.3, 0.1},
      std::array<double, 3>{0., 0., 0.});
  const auto fields = receiver.variables(
      x, 0.,
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
                 gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>{});
  const auto& g = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(fields);
  const auto& pi = get<gh::Tags::Pi<DataVector, 3>>(fields);
  const auto& phi = get<gh::Tags::Phi<DataVector, 3>>(fields);
  const auto geometry =
      gh::worldtube::spherical_trace_geometry(g, Tensor(n, 0.), x);
  const gh::worldtube::PhysicalScreen screen(
      g, geometry.inverse_metric, geometry.time_normal, geometry.normal_vector);
  tnsr::ii<DataVector, 3> incoming(n);
  for (size_t c = 0; c < incoming.size(); ++c) {
    for (size_t p = 0; p < n; ++p) {
      incoming[c][p] = 0.01 * sin(static_cast<double>(3 * c + p + 1));
    }
  }
  const auto target =
      gh::worldtube::incoming_weyl_to_physical_target(incoming, screen);
  const auto K = gh::extrinsic_curvature(geometry.time_normal, pi, phi);
  const DataVector lambda(n, -0.7);
  const auto bjorhus = [&](const tnsr::ii<DataVector, 3>* prescribed) {
    Tensor correction(n, 0.);
    gh::BoundaryConditions::Bjorhus::detail::add_physical_terms_to_dt_v_minus(
        make_not_null(&correction), Scalar<DataVector>(n, 0.4),
        geometry.normal_covector, geometry.normal_vector, geometry.time_normal,
        screen.lower, screen.mixed, screen.upper,
        geometry.inverse_spatial_metric, K, g, geometry.inverse_metric,
        Phi(n, 0.), Tensor(n, 0.), phi, tnsr::ijaa<DataVector, 3>(n, 0.),
        Phi(n, 0.),
        std::array<DataVector, 4>{DataVector(n, 0.), DataVector(n, 0.),
                                  DataVector(n, 0.), lambda},
        std::numeric_limits<double>::signaling_NaN(), nullptr,
        gh::BoundaryConditions::Bjorhus::default_incoming_wave_components,
        prescribed);
    return correction;
  };
  const auto off = bjorhus(nullptr);
  const auto on = bjorhus(&incoming);
  for (size_t c = 0; c < target.size(); ++c) {
    CHECK_ITERABLE_APPROX(DataVector(on[c] - off[c]),
                          DataVector(-lambda * target[c]));
  }
  CHECK(max(abs(get<0, 1>(target))) > 1.e-6);
  CHECK_ITERABLE_APPROX(target,
                        gh::worldtube::physical_projection(target, screen));

  const std::string options =
      "GaugeSector: SommerfeldAbsorbing\nAmplitude: 0\nFrequency: 0\n"
      "PrescribedQuadrupole:\n  Mass: 1\n  MomentPulse:\n"
      "    StartTime: 2\n    Duration: 8\n"
      "    Electric: [0.001, -0.002, 0.0003, -0.0004, 0.0005]\n"
      "    Magnetic: [-0.0004, 0.0002, 0.0007, 0.0006, -0.0003]\n";
  const auto boundary = TestHelpers::test_creation<
      gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost>(options);
  const auto restored = serialize_and_deserialize(boundary);
  const auto clone = boundary.get_clone();
  const auto& cloned = dynamic_cast<
      const gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost&>(*clone);
  const auto w = boundary.curvature_target(6., x, screen);
  const Target prescribed(
      1., Pulse(2., 8., {0.001, -0.002, 0.0003, -0.0004, 0.0005},
                {-0.0004, 0.0002, 0.0007, 0.0006, -0.0003}));
  CHECK_ITERABLE_APPROX(w, gh::worldtube::incoming_weyl_to_physical_target(
                               prescribed.incoming_mode(6., x), screen));
  CHECK(max(abs(get<1, 1>(w))) > 1.e-5);
  CHECK_ITERABLE_APPROX(w, restored.curvature_target(6., x, screen));
  CHECK_ITERABLE_APPROX(w, cloned.curvature_target(6., x, screen));
  CHECK_THROWS_WITH(
      (gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost(
          gh::BoundaryConditions::detail::SectorImposition::SommerfeldAbsorbing,
          1., 1., Target{})),
      Catch::Matchers::ContainsSubstring("requires Amplitude=0"));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.Gh.PhysicalBoundaryTarget",
                  "[Unit][Evolution][PhysicalTrace]") {
  register_classes_with_charm(
      tmpl::list<gh::Solutions::WrappedGr<gr::Solutions::Minkowski<3>>>{});
  test_prescribed_np();
  test_receiver_conversion();
}
