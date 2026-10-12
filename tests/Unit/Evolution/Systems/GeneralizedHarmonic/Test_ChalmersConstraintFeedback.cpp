// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/ModalVector.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ChalmersConstraintFeedback.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "ParallelAlgorithms/Amr/Policies/Policies.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using Feedback = gh::ChalmersWithConstraintFeedback<3>;
struct Metavariables {
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<amr::Criterion, tmpl::list<Feedback>>>;
  };
};
using Metric = gr::Tags::SpacetimeMetric<DataVector, 3>;
using Pi = gh::Tags::Pi<DataVector, 3>;
using Phi = gh::Tags::Phi<DataVector, 3>;
using Inverse =
    domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>;
}  // namespace

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.ChalmersConstraintFeedback",
                  "[Unit][Evolution]") {
  register_factory_classes_with_charm<Metavariables>();
  const Mesh<3> mesh{12, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const auto n = mesh.number_of_grid_points();
  Metric::type metric(n, 0.0);
  Pi::type pi(n, 0.0);
  Inverse::type inverse(n, 0.0);
  // Orthogonal directions rotated in x/y, with anisotropic and varying scales.
  const auto xi = logical_coordinates(mesh);
  inverse.get(0, 1) = 2.0 + 0.2 * get<0>(xi);
  inverse.get(1, 0) = -3.0;
  inverse.get(2, 2) = 0.5;
  metric.get(0, 0) = square(get<0>(xi)) + get<2>(xi);
  auto phi = partial_derivative(metric, mesh, inverse);
  CHECK_ITERABLE_APPROX(
      gh::chalmers_directional_constraints(metric, phi, mesh, inverse),
      (std::array{0.0, 0.0, 0.0}));
  // C_y00=2 and C_x01=C_x10=3: normalization removes map stretching,
  // rotation assigns directions correctly, and off-diagonals count twice.
  phi.get(1, 0, 0) -= 2.0;
  phi.get(0, 0, 1) -= 3.0;
  CHECK_ITERABLE_APPROX(
      gh::chalmers_directional_constraints(metric, phi, mesh, inverse),
      (std::array{2.0, std::sqrt(18.0), 0.0}));

  // Constant fields are spectrally resolved, but can violate reduction
  // constraints. Veto only the bad p direction and retain other coarsening.
  metric = Metric::type(n, 0.0);
  phi = Phi::type(n, 0.0);
  phi.get(1, 0, 0) = 2.e-3;
  const Feedback::FieldCriterion fields{{"SpacetimeMetric", "Pi", "Phi"},
                                        1.e-4,
                                        0.0,
                                        0.5,
                                        false,
                                        true,
                                        0.05,
                                        false};
  Parallel::GlobalCache<Metavariables> cache{};
  auto box = db::create<
      tmpl::list<domain::Tags::Mesh<3>, Metric, Pi, Phi, Inverse, ::Tags::Time,
                 amr::Tags::CriterionNeighborData<3>, amr::Tags::Policies>>(
      mesh, metric, pi, phi, inverse, 0.0, amr::CriterionNeighborData<3>{},
      amr::Policies{amr::Isotropy::Anisotropic, amr::Limits{}, true, true});
  const ElementId<3> id{
      0, std::array{SegmentId{1, 0}, SegmentId{0, 0}, SegmentId{0, 0}}};
  const auto field_flags = fields(box, cache, id);
  CHECK(field_flags == make_array<3>(amr::Flag::DecreaseResolution));
  const Feedback baseline{fields, 1.e-3, 1.0, false, false};
  CHECK(baseline(box, cache, id) == field_flags);
  const Feedback veto{fields, 1.e-3, 1.0, true, false};
  CHECK(veto(box, cache, id) ==
        (std::array{amr::Flag::DoNothing, amr::Flag::DecreaseResolution,
                    amr::Flag::DecreaseResolution}));
  CHECK(serialize_and_deserialize(veto)(box, cache, id) ==
        veto(box, cache, id));

  // The ordinary sibling machinery still handles joins. A bad transverse
  // constraint vetoes a join, even when the joining direction is good.
  const Feedback::FieldCriterion joins{{"SpacetimeMetric", "Pi", "Phi"},
                                       1.e-4,
                                       0.0,
                                       0.5,
                                       false,
                                       true,
                                       0.05,
                                       true};
  const Feedback join_guard{joins, 1.e-3, 1.0, true, false};
  CHECK(join_guard.needs_sibling_data());
  db::mutate<Phi>(
      [](const auto p) {
        *p = Phi::type(p->begin()->size(), 0.0);
        p->get(0, 0, 0) = 2.e-3;
      },
      make_not_null(&box));
  const auto packed = join_guard.pack_sibling_data(box);
  db::mutate<amr::Tags::CriterionNeighborData<3>>(
      [&packed, &mesh](const auto state) {
        state->round = 1;
        state->sibling.emplace(
            ElementId<3>{0, std::array{SegmentId{1, 1}, SegmentId{0, 0},
                                       SegmentId{0, 0}}},
            amr::CriterionElementData<3>{1, mesh, packed});
      },
      make_not_null(&box));
  CHECK(joins(box, cache, id)[0] == amr::Flag::Join);
  CHECK(join_guard(box, cache, id)[0] == amr::Flag::DoNothing);
  CHECK((Feedback{joins, 1.e-2, 1.0, true, false}(box, cache, id)[0]) ==
        amr::Flag::Join);

  // A smooth spectrum just below target becomes p-refined when only its
  // directional target is tightened. The bounded factor is checked by
  // comparison with the same field criterion at explicit factors.
  ModalVector modes(n, 0.0);
  for (size_t k = 0; k < 11; ++k) {
    modes[k] = std::exp(-0.8 * static_cast<double>(k));
  }
  metric.get(0, 0) = to_nodal_coefficients(modes, mesh);
  phi = partial_derivative(metric, mesh, inverse);
  phi.get(1, 0, 0) += 1.0;
  db::mutate<Metric, Phi>(
      [&](const auto g, const auto p) {
        *g = metric;
        *p = phi;
      },
      make_not_null(&box));
  const double error = ::Chalmers::estimate_directional_energy(
                           metric.get(0, 0), mesh, 0.0, 1, true)[0]
                           .error;
  const Feedback::FieldCriterion smooth{
      {"SpacetimeMetric"}, 1.2 * error, 0.0, 0.5, false, true, 0.0, false};
  CHECK(smooth(box, cache, id)[0] == amr::Flag::DoNothing);
  const Feedback tightened{smooth, 1.e-3, 0.1, true, false};
  CHECK(tightened(box, cache, id)[0] == amr::Flag::IncreaseResolution);
  CHECK(tightened(box, cache, id) ==
        smooth.evaluate_with_target_factors(box, cache, id,
                                            std::array{0.1, 1.0, 1.0}));
  const Feedback::FieldCriterion slow{
      {"SpacetimeMetric"}, 0.5 * error, 0.0, 1.0, false, true, 0.0, false};
  CHECK((Feedback{slow, 1.e-3, 0.1, true, false}(box, cache, id)[0]) ==
        amr::Flag::Split);
  // Once the constraint disappears, target tightening also disappears.
  db::mutate<Phi>(
      [&](const auto p) { *p = partial_derivative(metric, mesh, inverse); },
      make_not_null(&box));
  CHECK(tightened(box, cache, id) == smooth(box, cache, id));

  const auto parsed =
      TestHelpers::test_factory_creation<amr::Criterion, Feedback>(
          "ChalmersWithConstraintFeedback:\n"
          "  Chalmers:\n"
          "    VariablesToMonitor: [SpacetimeMetric, Pi, Phi]\n"
          "    AbsoluteTarget: 1.e-4\n"
          "    RelativeTarget: 0\n"
          "    SmoothnessThreshold: 0.5\n"
          "    PRefineOnly: false\n"
          "    UseDirectionalEnergy: true\n"
          "    CoarseningFactor: 0.05\n"
          "    HCoarsening: true\n"
          "  ConstraintTarget: 1.e-4\n"
          "  MinimumTargetFactor: 0.1\n"
          "  EnableCoarseningVeto: true\n"
          "  Verbose: false\n");
  CHECK(parsed->needs_sibling_data());
  CHECK(parsed->type() == amr::Criteria::Type::hp);
  CHECK_THROWS_WITH((Feedback{fields, 0.0, 0.1, true, false}),
                    Catch::Matchers::ContainsSubstring("finite and positive"));
  CHECK_THROWS_WITH((Feedback{fields, 1.e-3, 0.0, true, false}),
                    Catch::Matchers::ContainsSubstring("(0, 1]"));
}
