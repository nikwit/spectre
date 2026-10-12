// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/ModalVector.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ChalmersWithSpecCoarsening.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Metric = gr::Tags::SpacetimeMetric<DataVector, 3>;
using Pi = gh::Tags::Pi<DataVector, 3>;
using Phi = gh::Tags::Phi<DataVector, 3>;
using Inverse =
    domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>;
using Coordinates = domain::Tags::Coordinates<3, Frame::Inertial>;
using Criterion = gh::ChalmersWithSpecCoarsening;
struct Metavariables {
  static constexpr size_t volume_dim = 3;
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<amr::Criterion, tmpl::list<Criterion>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.Spectrum",
                  "[Unit][Evolution]") {
  using gh::spec_amr::diagnose_power;
  DataVector power(16, 0.0);
  for (size_t k = 0; k < 15; ++k) {
    power[k] = pow(10.0, -0.7 * k);
  }
  const auto exponential = diagnose_power(power);
  CHECK(exponential.convergence == approx(0.7));
  CHECK(exponential.pile_up < 1.e-10);
  CHECK(exponential.next_error / exponential.error == approx(pow(10., 0.7)));
  power[15] = 1.e100;  // the killed mode never drives a decision
  CHECK(diagnose_power(power).error == exponential.error);
  CHECK(diagnose_power(DataVector(16, 0.0)).degenerate);
  CHECK(diagnose_power(DataVector(16, 1.0)).pile_up == 0.0);
  // A true flattened tail after a well-resolved exponential prefix.
  power = DataVector(20, 0.0);
  for (size_t k = 0; k < 19; ++k) {
    power[k] = pow(10., -0.7 * std::min(k, size_t{10}));
  }
  CHECK(diagnose_power(power).pile_up > 3.0);
  // Parity zeros must not manufacture a plateau or an artificially tiny tail.
  for (size_t k = 0; k < 19; ++k) {
    power[k] = k % 2 == 0 ? pow(10., -0.5 * k) : 0.;
  }
  CHECK(diagnose_power(power).convergence == approx(0.5));
  CHECK(diagnose_power(power).pile_up < 1.e-10);
  CHECK(diagnose_power(power).next_error >= diagnose_power(power).error);
  power = DataVector(12, 0.0);
  power[7] = 1.0;
  CHECK(diagnose_power(power).next_error == 0.0);
  power[10] = 1.0;
  CHECK(diagnose_power(power).error == 1.0);
  CHECK_FALSE(diagnose_power(power).degenerate);
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.Policy",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  Settings settings{1.e-4, {{1.e-3, 1.e-3}}, false, 5, 2};
  Diagnostics data{};
  for (auto& group : data.spectra) {
    group.fill(Spectrum{2.e-5, 4.e-5, 0.5, 0.0, false});
  }
  for (auto& group : data.constraints) {
    group.fill(1.e-4);
  }
  const std::array<size_t, 3> extents{{12, 12, 8}};
  const auto none = make_array<3>(amr::Flag::DoNothing);
  const auto down = make_array<3>(amr::Flag::DecreaseResolution);
  History history{};
  CHECK(decide(make_not_null(&history), data, extents, 0., none, settings) ==
        none);
  CHECK(decide(make_not_null(&history), data, extents, 1., none, settings) ==
        none);
  // Current constraints increased 20%, yet the predicted coarsened norm
  // remains below tolerance. Monotonic constraints are not required.
  for (auto& group : data.constraints) {
    group.fill(1.2e-4);
  }
  CHECK(decide(make_not_null(&history), data, extents, 2., none, settings) ==
        down);
  const auto copies = history.calls;
  CHECK(decide(make_not_null(&history), data, extents, 2., none, settings) ==
        down);
  CHECK(history.calls == copies);
  // P addition is governed by the same tensor budgets. A componentwise
  // Chalmers request cannot undo an otherwise justified p reduction.
  auto add_history = history;
  CHECK(decide(make_not_null(&add_history), data, extents, 3.,
               make_array<3>(amr::Flag::IncreaseResolution), settings) == down);
  // An accepted sibling join must not acquire p decreases in other axes:
  // the prospective-parent test was performed at the existing child extents.
  auto join_history = history;
  auto join = none;
  join[1] = amr::Flag::Join;
  CHECK(decide(make_not_null(&join_history), data, extents, 3., join,
               settings) == join);
  auto restored = serialize_and_deserialize(history);
  auto reduced = extents;
  --reduced[0];
  CHECK(decide(make_not_null(&restored), data, reduced, 3., none,
               settings)[0] == amr::Flag::DoNothing);
  CHECK(decide(make_not_null(&restored), data, reduced, 4., none,
               settings)[0] == amr::Flag::DecreaseResolution);

  // A rising constraint forecast prevents ordinary coarsening.
  data.constraints[0][0] = 5.e-4;
  CHECK(decide(make_not_null(&history), data, extents, 3., none, settings)[0] ==
        amr::Flag::DoNothing);
  // Forced pile-up can override both error budgets and p-addition, but only
  // if BOTH monitored tensor groups show pile-up.
  for (auto& group : data.spectra) {
    group[0] = Spectrum{0.01, 0.02, 0.4, 5.0, false};
  }
  settings.force_pile_up = true;
  auto up = none;
  up[0] = amr::Flag::IncreaseResolution;
  CHECK(decide(make_not_null(&history), data, extents, 4., up, settings)[0] ==
        amr::Flag::DecreaseResolution);
  data.spectra[1][0].pile_up = 0.0;
  CHECK(decide(make_not_null(&history), data, extents, 5., up, settings)[0] ==
        amr::Flag::IncreaseResolution);
  data.spectra[1][0].pile_up = 3.5;
  CHECK(decide(make_not_null(&history), data, extents, 6., up, settings)[0] ==
        amr::Flag::DoNothing);  // avoid immediately adding the point back
  data.spectra[1][0].pile_up = 5.;
  auto minimum = extents;
  minimum[0] = 5;
  decide(make_not_null(&history), data, minimum, 7., up, settings);
  CHECK(decide(make_not_null(&history), data, minimum, 8., up, settings)[0] ==
        amr::Flag::DoNothing);
  up[0] = amr::Flag::Split;
  CHECK(decide(make_not_null(&history), data, minimum, 9., up, settings)[0] ==
        amr::Flag::Split);
  up[0] = amr::Flag::Join;
  data.constraints[1][2] = 1.;
  CHECK(decide(make_not_null(&history), data, minimum, 10., up, settings)[0] ==
        amr::Flag::DoNothing);
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.TargetShift",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  Settings settings{1.e-4, {{1.e-3, 1.e-3}}, true, 7, 2, -1.};
  Diagnostics data{};
  for (auto& group : data.spectra) {
    group.fill(Spectrum{4.e-4, 5.e-4, 0.5, 0., false});
  }
  for (auto& group : data.constraints) {
    group.fill(1.2e-3);
  }
  History history{};
  history.calls = 3;
  history.extents = {{12, 12, 8}};
  history.time = 0.;
  for (auto& group : history.constraints) {
    group.fill(1.e-3);
  }
  const auto check_shift = [&](const double expected) {
    CHECK_ITERABLE_APPROX(updated_target_shifts(history, data, settings),
                          make_array<3>(expected));
  };
  // A sufficiently unresolved spectrum already demands refinement and does
  // not need additional tightening from the constraints.
  check_shift(0.);
  for (auto& group : data.spectra) {
    for (auto& spectrum : group) {
      spectrum.error = 2.e-4;
    }
  }
  check_shift(log10(2.) - 0.4);  // 0.8-mode placement
  history.log10_target_shift.fill(-0.2);
  check_shift(-0.45);  // rising constraints: half a mode is more restrictive
  for (auto& group : data.spectra) {
    for (auto& spectrum : group) {
      spectrum.error = 8.e-5;
    }
  }
  check_shift(log10(0.8) - 0.4);  // now placement is more restrictive
  for (auto& group : history.constraints) {
    group.fill(1.3e-3);
  }
  check_shift(-0.45);  // decreasing, but forecast still exceeds the target
  history.log10_target_shift.fill(-0.9);
  check_shift(-1.);  // bounded tightening
  history.log10_target_shift.fill(-1.);
  check_shift(-1.);

  // Both constraints must be small AND decreasing to relax. Recovery takes
  // 0.1 mode per call, not the 0.5-mode tightening step.
  history.log10_target_shift.fill(-0.5);
  for (auto& group : data.constraints) {
    group.fill(4.e-4);
  }
  for (auto& group : history.constraints) {
    group.fill(5.e-4);
  }
  check_shift(-0.45);
  history.constraints[1].fill(4.e-4);
  check_shift(-0.5);  // one stationary constraint prevents relaxation
  history.constraints[1].fill(6.e-4);
  data.constraints[1].fill(5.e-4);
  check_shift(-0.5);  // strict half-tolerance recovery threshold
  data.constraints[1].fill(4.e-4);
  history.log10_target_shift.fill(-0.01);
  check_shift(0.);  // no relaxation beyond the original target

  // The forecast can tighten even when the current value is below tolerance.
  history.log10_target_shift.fill(-0.2);
  for (auto& group : history.constraints) {
    group.fill(0.);
  }
  for (auto& group : data.constraints) {
    group.fill(7.5e-4);
  }
  check_shift(log10(0.8) - 0.4);  // tiny-old-value fallback predicts 1.5e-3
  history.calls = 0;
  check_shift(-0.2);  // without history, current value alone is below target
  data.constraints[1].fill(1.2e-3);
  check_shift(-0.45);  // first-call tightening by half a mode
  history.calls = 3;

  // Directions and tensor groups are independent. The smaller usable slope
  // sets the step, while a nondecaying or absent fit cannot define one.
  history.log10_target_shift.fill(-0.2);
  for (auto& group : history.constraints) {
    group.fill(1.3e-3);
  }
  for (auto& group : data.constraints) {
    group.fill(1.2e-3);
  }
  data.spectra[1][0].convergence = 0.2;
  data.spectra[1][1].convergence = -0.1;
  data.spectra[0][2].degenerate = true;
  data.spectra[1][2].degenerate = true;
  CHECK_ITERABLE_APPROX(updated_target_shifts(history, data, settings),
                        (std::array{-0.3, -0.2, -0.2}));
  CHECK(history.log10_target_shift == make_array<3>(-0.2));
  settings.adjustment_limit = 0.;
  check_shift(0.);  // disabled controller retains legacy targets
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.SharedTarget",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  Settings settings{1.e-4, {{1.e-3, 1.e-3}}, true, 7, 2, -1.};
  Diagnostics data{};
  for (auto& group : data.spectra) {
    group.fill(Spectrum{5.e-5, 7.e-5, 0.5, 3.3, false});
  }
  for (auto& group : data.constraints) {
    group.fill(1.2e-3);
  }
  History history{};
  history.calls = 3;
  history.extents = {{12, 12, 8}};
  history.time = 0.;
  for (auto& group : history.constraints) {
    group.fill(1.e-3);
  }
  const auto extents = history.extents;
  const auto none = make_array<3>(amr::Flag::DoNothing);
  const auto up = make_array<3>(amr::Flag::IncreaseResolution);
  const auto down = make_array<3>(amr::Flag::DecreaseResolution);
  auto legacy = settings;
  legacy.adjustment_limit = 0.;
  auto legacy_history = history;
  // This is the original remove-then-add failure: the next error fits the
  // unshifted target, but not the constraint-adjusted target.
  CHECK(decide(make_not_null(&legacy_history), data, extents, 1., up, legacy) ==
        down);
  CHECK(decide(make_not_null(&history), data, extents, 1., up, settings) ==
        none);
  CHECK_ITERABLE_APPROX(history.log10_target_shift,
                        make_array<3>(log10(0.5) - 0.4));
  const auto once = history;
  CHECK(decide(make_not_null(&history), data, extents, 1., up, settings) ==
        none);
  CHECK(history.log10_target_shift == once.log10_target_shift);
  CHECK(history.calls == once.calls);
  for (auto& group : data.spectra) {
    for (auto& spectrum : group) {
      spectrum.pile_up = 0.;
    }
  }
  // Refinement uses the shifted threshold too, even below the base target.
  CHECK(decide(make_not_null(&history), data, extents, 2., up, settings) == up);
  // SpEC additionally requires the base spectral or constraint permission.
  // A lingering tight target alone must not add points after recovery.
  auto quiet_history = history;
  auto quiet_data = data;
  for (auto& group : quiet_history.constraints) {
    group.fill(1.e-4);
  }
  for (auto& group : quiet_data.constraints) {
    group.fill(1.e-4);
  }
  CHECK(decide(make_not_null(&quiet_history), quiet_data, extents, 3., up,
               settings) == none);
  // High constraints no longer bypass that shared threshold by themselves.
  for (auto& group : data.spectra) {
    for (auto& spectrum : group) {
      spectrum.error = 1.e-6;
      spectrum.next_error = 2.e-6;
    }
  }
  CHECK(decide(make_not_null(&history), data, extents, 3., up, settings) ==
        none);
  CHECK(history.log10_target_shift == make_array<3>(-1.));
  // Preserve Chalmers' h decision and the strong-pile-up escape route.
  CHECK(decide(make_not_null(&history), data, extents, 4.,
               make_array<3>(amr::Flag::Split),
               settings) == make_array<3>(amr::Flag::Split));
  for (auto& group : data.spectra) {
    group.fill(Spectrum{1.e-2, 2.e-2, 0.5, 4.2, false});
  }
  CHECK(decide(make_not_null(&history), data, extents, 5., up, settings) ==
        down);

  // Checkpointing preserves the shifted target and constraint trend. A p
  // change triggers the cooldown but does not erase or freeze target history.
  history.log10_target_shift.fill(-0.2);
  for (auto& group : data.spectra) {
    group.fill(Spectrum{5.e-5, 7.e-5, 0.5, 0., false});
  }
  auto restored = serialize_and_deserialize(history);
  CHECK(restored.log10_target_shift == history.log10_target_shift);
  CHECK(restored.constraints == history.constraints);
  auto reduced = extents;
  --reduced[0];
  const auto expected_shift = updated_target_shifts(restored, data, settings);
  CHECK(decide(make_not_null(&restored), data, reduced, 6., up, settings)[0] ==
        amr::Flag::DoNothing);
  CHECK(restored.log10_target_shift == expected_shift);
  CHECK(restored.log10_target_shift[0] < -0.2);
  ProjectHistory::apply(make_not_null(&restored),
                        std::pair{Mesh<3>{8, Spectral::Basis::Legendre,
                                          Spectral::Quadrature::GaussLobatto},
                                  Element<3>{ElementId<3>{0}, {}}});
  CHECK(restored.log10_target_shift == expected_shift);
  ProjectHistory::apply(make_not_null(&restored), 0);
  CHECK(restored.calls == 0);
  CHECK(restored.log10_target_shift == make_array<3>(0.));

  // Recovery eventually permits accuracy reduction again. Keep the mesh
  // fixed to isolate the controller's slow release of its stored target.
  history.log10_target_shift.fill(-1.);
  for (auto& group : history.constraints) {
    group.fill(4.e-4);
  }
  for (auto& group : data.spectra) {
    group.fill(Spectrum{5.e-5, 7.e-5, 1., 0., false});
  }
  auto decision = none;
  for (size_t step = 0; step < 10; ++step) {
    for (size_t c = 0; c < 2; ++c) {
      for (size_t d = 0; d < 3; ++d) {
        data.constraints[c][d] = 0.8 * history.constraints[c][d];
      }
    }
    decision = decide(make_not_null(&history), data, extents, 6. + step, none,
                      settings);
    if (step == 0) {
      CHECK(decision == none);
      CHECK_ITERABLE_APPROX(history.log10_target_shift, make_array<3>(-0.9));
    }
  }
  CHECK(decision == down);
  CHECK_ITERABLE_APPROX(history.log10_target_shift, make_array<3>(0.));
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.DiagnosticsAndHistory",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  const Mesh<3> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const size_t n = mesh.number_of_grid_points();
  Metric::type metric(n, 0.0);
  metric.get(0, 0) = -1.;
  for (size_t i = 1; i < 4; ++i) {
    metric.get(i, i) = 1.;
  }
  Pi::type pi(n, 0.0);
  Phi::type phi(n, 0.0);
  Inverse::type inverse(n, 0.0);
  Coordinates::type coords(n, 0.0);
  const auto xi = logical_coordinates(mesh);
  for (size_t d = 0; d < 3; ++d) {
    inverse.get(d, d) = 1. / (d + 1.);
    coords.get(d) = (d + 1.) * xi.get(d);
  }
  const gh::gauges::Harmonic harmonic{};
  const Scalar<DataVector> gamma2(n, 1.0);
  const auto minkowski = gh::spec_amr::diagnose(metric, pi, phi, mesh, inverse,
                                                coords, 0., harmonic, gamma2);
  for (const auto& constraint : minkowski.constraints) {
    CHECK(max(abs(DataVector(constraint))) < 1.e-12);
  }
  pi.get(0, 0) = 2.e-3 * coords.get(0);
  pi.get(0, 1) = 3.e-3 * coords.get(1);
  const auto manufactured = gh::spec_amr::diagnose(
      metric, pi, phi, mesh, inverse, coords, 0., harmonic, gamma2);
  CHECK_ITERABLE_APPROX(manufactured.constraints[0],
                        (std::array{1.e-3, 3.e-3, 0.}));
  CHECK_ITERABLE_APPROX(manufactured.constraints[1], (std::array{0., 0., 0.}));
  phi.get(2, 0, 1) = 4.e-3;
  const auto reduction = gh::spec_amr::diagnose(metric, pi, phi, mesh, inverse,
                                                coords, 0., harmonic, gamma2);
  CHECK_ITERABLE_APPROX(reduction.constraints[1], (std::array{0., 0., 4.e-3}));
  // Restore Minkowski for a wrapper and polymorphic-dispatch test.
  pi = Pi::type(n, 0.0);
  phi = Phi::type(n, 0.0);
  const gh::ChalmersWithConstraintFeedback<3>::FieldCriterion field{
      {"SpacetimeMetric", "Pi", "Phi"}, 1.e-4, 0., 0.5, false, true, 0., false};
  const gh::ChalmersWithConstraintFeedback<3> refinement{field, 1.e-3, 1.,
                                                         false, false};
  const Criterion criterion{
      refinement, 1.e-4, {{1.e-3, 1.e-3}}, true, true, 7, 2, 0., false};
  auto box = db::create<
      tmpl::list<domain::Tags::Mesh<3>, Metric, Pi, Phi, Inverse, Coordinates,
                 ::Tags::Time, gh::gauges::Tags::GaugeCondition,
                 gh::Tags::ConstraintGamma2, gh::spec_amr::HistoryTag>>(
      mesh, metric, pi, phi, inverse, coords, 0.,
      std::unique_ptr<gh::gauges::GaugeCondition>(
          std::make_unique<gh::gauges::Harmonic>()),
      gamma2, gh::spec_amr::History{});
  auto observation = make_observation_box<tmpl::list<>>(make_not_null(&box));
  Parallel::GlobalCache<Metavariables> cache{};
  const ElementId<3> id{0};
  const amr::Criterion& base = criterion;
  base.evaluate(observation, cache, id);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).calls == 0);
  base.evaluate_and_update(make_not_null(&observation), cache, id);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).calls == 1);
  base.evaluate(observation, cache, id);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).calls == 1);
  CHECK(serialize_and_deserialize(criterion)(box, cache, id) ==
        criterion(box, cache, id));
  CHECK_THROWS_WITH(
      (Criterion{refinement, 0., {{1., 1.}}, true, true, 7, 2, 0., false}),
      Catch::Matchers::ContainsSubstring("Invalid spectral"));
  CHECK_THROWS_WITH(
      (Criterion{refinement, 1.e-4, {{1., 1.}}, true, true, 6, 2, 0., false}),
      Catch::Matchers::ContainsSubstring("minimum modes"));
  auto history = db::get<gh::spec_amr::HistoryTag>(box);
  gh::spec_amr::ProjectHistory::apply(make_not_null(&history),
                                      std::pair{mesh, Element<3>{id, {}}});
  CHECK(history.calls == 1);
  gh::spec_amr::ProjectHistory::apply(make_not_null(&history), 0);
  CHECK(history.calls == 0);
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.WrapperTargetShift",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  const Mesh<3> mesh{12, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const size_t n = mesh.number_of_grid_points();
  Metric::type metric(n, 0.0);
  for (size_t i = 1; i < 4; ++i) {
    metric.get(i, i) = 1.;
  }
  ModalVector modes(n, 0.0);
  for (size_t k = 0; k < 11; ++k) {
    modes[k] = 0.01 * exp(-0.8 * static_cast<double>(k));
  }
  metric.get(0, 0) = -1. + to_nodal_coefficients(modes, mesh);
  Pi::type pi(n, 0.0);
  Inverse::type inverse(n, 0.0);
  Coordinates::type coords(n, 0.0);
  const auto xi = logical_coordinates(mesh);
  for (size_t d = 0; d < 3; ++d) {
    inverse.get(d, d) = 1.;
    coords.get(d) = xi.get(d);
  }
  auto phi = partial_derivative(metric, mesh, inverse);
  phi.get(0, 0, 0) += 0.01;
  const Scalar<DataVector> gamma2(n, 1.);
  const auto diagnostics =
      gh::spec_amr::diagnose(metric, pi, phi, mesh, inverse, coords, 0.,
                             gh::gauges::Harmonic{}, gamma2);
  const double tensor_error = std::max(diagnostics.spectra[0][0].error,
                                       diagnostics.spectra[1][0].error);
  REQUIRE(tensor_error > 0.);
  const double component_error = ::Chalmers::estimate_directional_energy(
                                     metric.get(0, 0), mesh, 0., 1, true)[0]
                                     .error;
  const gh::ChalmersWithConstraintFeedback<3>::FieldCriterion field{
      {"SpacetimeMetric"},
      1.2 * component_error,
      0.,
      0.5,
      false,
      true,
      0.,
      false};
  const gh::ChalmersWithConstraintFeedback<3> refinement{field, 1.e9, 1., false,
                                                         false};
  const auto make_criterion = [&](const double start, const bool enabled) {
    return Criterion{refinement,
                     1.2 * tensor_error,
                     {{1.e9, 1.e-4}},
                     enabled,
                     true,
                     7,
                     2,
                     start,
                     false,
                     -1.};
  };
  const auto criterion = make_criterion(0., true);
  gh::spec_amr::History initial_history{};
  initial_history.calls = 3;
  initial_history.time = -1.;
  initial_history.extents = mesh.extents().indices();
  initial_history.constraints = diagnostics.constraints;
  auto box = db::create<
      tmpl::list<domain::Tags::Mesh<3>, Metric, Pi, Phi, Inverse, Coordinates,
                 ::Tags::Time, gh::gauges::Tags::GaugeCondition,
                 gh::Tags::ConstraintGamma2, gh::spec_amr::HistoryTag>>(
      mesh, metric, pi, phi, inverse, coords, 0.,
      std::unique_ptr<gh::gauges::GaugeCondition>(
          std::make_unique<gh::gauges::Harmonic>()),
      gamma2, initial_history);
  auto observation = make_observation_box<tmpl::list<>>(make_not_null(&box));
  Parallel::GlobalCache<Metavariables> cache{};
  const ElementId<3> id{0};
  CHECK(refinement(box, cache, id)[0] == amr::Flag::DoNothing);
  const amr::Criterion& base = criterion;
  // This call must actually request refinement of the smooth field, not only
  // veto a coarsening request. Observation must not accumulate target shifts.
  CHECK(base.evaluate(observation, cache, id)[0] ==
        amr::Flag::IncreaseResolution);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).calls == 3);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).log10_target_shift ==
        make_array<3>(0.));
  CHECK(serialize_and_deserialize(criterion)(box, cache, id)[0] ==
        amr::Flag::IncreaseResolution);
  CHECK(base.evaluate_and_update(make_not_null(&observation), cache, id)[0] ==
        amr::Flag::IncreaseResolution);
  const auto after = db::get<gh::spec_amr::HistoryTag>(box);
  CHECK(after.calls == 4);
  CHECK(after.log10_target_shift[0] < 0.);
  base.evaluate_and_update(make_not_null(&observation), cache, id);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).calls == after.calls);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).log10_target_shift ==
        after.log10_target_shift);

  const auto reset_history = [&]() {
    db::mutate<gh::spec_amr::HistoryTag>(
        [&](const auto h) { *h = initial_history; }, make_not_null(&box));
    observation.reset();
  };
  reset_history();
  const auto delayed = make_criterion(2., true);
  CHECK(delayed.evaluate_with_history(make_not_null(&observation), cache,
                                      id)[0] == amr::Flag::DoNothing);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).log10_target_shift ==
        make_array<3>(0.));
  db::mutate<::Tags::Time>([](const auto t) { *t = 2.; }, make_not_null(&box));
  observation.reset();
  CHECK(delayed.evaluate_with_history(make_not_null(&observation), cache,
                                      id)[0] == amr::Flag::IncreaseResolution);
  reset_history();
  const auto disabled = make_criterion(0., false);
  CHECK(disabled.evaluate_with_history(make_not_null(&observation), cache,
                                       id)[0] == amr::Flag::DoNothing);
  CHECK(db::get<gh::spec_amr::HistoryTag>(box).log10_target_shift ==
        make_array<3>(0.));

  const auto parsed =
      TestHelpers::test_factory_creation<amr::Criterion, Criterion>(
          "ChalmersWithSpecCoarsening:\n"
          "  RefinementCriterion:\n"
          "    Chalmers:\n"
          "      VariablesToMonitor: [SpacetimeMetric, Pi, Phi]\n"
          "      AbsoluteTarget: 1.e-4\n"
          "      RelativeTarget: 0\n"
          "      SmoothnessThreshold: 0.5\n"
          "      PRefineOnly: false\n"
          "      UseDirectionalEnergy: true\n"
          "      CoarseningFactor: 0.05\n"
          "      HCoarsening: false\n"
          "    ConstraintTarget: 1.e-4\n"
          "    MinimumTargetFactor: 0.1\n"
          "    EnableCoarseningVeto: true\n"
          "    Verbose: false\n"
          "  RelativeTruncationTarget: 1.e-4\n"
          "  ConstraintTargets: [1.e-3, 1.e-3]\n"
          "  EnableCoarsening: true\n"
          "  ForceDropPileUpModes: true\n"
          "  MinimumModes: 7\n"
          "  CallsBetweenChanges: 2\n"
          "  StartTime: 0\n"
          "  Verbose: false\n"
          "  TruncationErrorAdjustmentLimit: -1\n");
  CHECK(parsed->type() == amr::Criteria::Type::hp);
  for (const double limit :
       {0.1, -1000., std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_WITH(
        (Criterion{refinement,
                   1.e-4,
                   {{1.e-3, 1.e-3}},
                   true,
                   true,
                   7,
                   2,
                   0.,
                   false,
                   limit}),
        Catch::Matchers::ContainsSubstring("TruncationErrorAdjustmentLimit"));
  }
}
