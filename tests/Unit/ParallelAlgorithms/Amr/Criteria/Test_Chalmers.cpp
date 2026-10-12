// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Chalmers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
struct Field : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct OtherField : db::SimpleTag {
  using type = Scalar<DataVector>;
};
using Criterion = amr::Criteria::Chalmers<1, tmpl::list<Field>>;
using FilteredCriterion = amr::Criteria::Chalmers<1, tmpl::list<Field>, true>;
using ParityCriterion =
    amr::Criteria::Chalmers<1, tmpl::list<Field>, true, true>;
struct Metavariables {
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<amr::Criterion,
                   tmpl::list<Criterion, FilteredCriterion, ParityCriterion>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Amr.Criteria.ChalmersParity",
                  "[Unit][ParallelAlgorithms]") {
  register_factory_classes_with_charm<Metavariables>();
  const auto criterion =
      TestHelpers::test_factory_creation<amr::Criterion, ParityCriterion>(
          "ChalmersParityIgnoringHighestMode:\n"
          "  VariablesToMonitor: [Field]\n"
          "  AbsoluteTarget: 3.e-5\n"
          "  RelativeTarget: 0\n"
          "  SmoothnessThreshold: 0.5\n"
          "  PRefineOnly: false\n"
          "  CoarseningFactor: 0.0\n"
          "  HCoarsening: false\n");
  CHECK(criterion->observation_name() == "ChalmersParityIgnoringHighestMode");
  const auto serialized = serialize_and_deserialize(
      dynamic_cast<const ParityCriterion&>(*criterion));
  const Mesh<1> mesh{12, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  Parallel::GlobalCache<Metavariables> cache{};
  for (const double noise : {0.0, 2.e-15, 1.e-13}) {
    ModalVector modes(12, noise);
    for (size_t k = 1; k < 11; k += 2) {
      modes[k] = 0.1 * std::exp(-0.8 * k);
    }
    modes[11] = 0.0;
    const auto box = db::create<tmpl::list<domain::Tags::Mesh<1>, Field>>(
        mesh, Scalar<DataVector>{to_nodal_coefficients(modes, mesh)});
    CHECK(serialized(box, cache, ElementId<1>{0})[0] == amr::Flag::DoNothing);
  }
}

SPECTRE_TEST_CASE("Unit.Amr.Criteria.Chalmers", "[Unit][ParallelAlgorithms]") {
  register_factory_classes_with_charm<Metavariables>();
  const auto criterion =
      TestHelpers::test_factory_creation<amr::Criterion, Criterion>(
          "Chalmers:\n"
          "  VariablesToMonitor: [Field]\n"
          "  AbsoluteTarget: 1.e-8\n"
          "  RelativeTarget: 0.0\n"
          "  SmoothnessThreshold: 1.0\n"
          "  PRefineOnly: false\n"
          "  CoarseningFactor: 0.0\n"
          "  HCoarsening: false\n");
  const auto filtered =
      TestHelpers::test_factory_creation<amr::Criterion, FilteredCriterion>(
          "ChalmersIgnoringHighestMode:\n"
          "  VariablesToMonitor: [Field]\n"
          "  AbsoluteTarget: 1.e-8\n"
          "  RelativeTarget: 0\n"
          "  SmoothnessThreshold: 0.5\n"
          "  PRefineOnly: false\n"
          "  CoarseningFactor: 0.1\n"
          "  HCoarsening: false\n");
  CHECK(filtered->observation_name() == "ChalmersIgnoringHighestMode");
  CHECK(criterion->type() == amr::Criteria::Type::hp);
  const Criterion low_threshold{{"Field"}, 1.e-8, 0.0, 0.5, false};
  const Criterion p_only{{"Field"}, 1.e-8, 0.0, 1.0, true};
  const auto serialized = serialize_and_deserialize(low_threshold);
  const Mesh<1> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  Parallel::GlobalCache<Metavariables> cache{};
  const auto evaluate = [&](const auto& c, const double sigma,
                            const double amplitude) {
    ModalVector modes(8);
    for (size_t k = 0; k < modes.size(); ++k) {
      modes[k] = amplitude * std::exp(-sigma * static_cast<double>(k));
    }
    const auto box = db::create<tmpl::list<domain::Tags::Mesh<1>, Field>>(
        mesh, Scalar<DataVector>{to_nodal_coefficients(modes, mesh)});
    return c(box, cache, ElementId<1>{0})[0];
  };
  CHECK(evaluate(serialize_and_deserialize(
                     dynamic_cast<const FilteredCriterion&>(*filtered)),
                 0.7, 1.0) == amr::Flag::IncreaseResolution);
  const auto& parsed = dynamic_cast<const Criterion&>(*criterion);
  CHECK(evaluate(parsed, 0.7, 1.0) == amr::Flag::Split);
  CHECK(evaluate(low_threshold, 0.7, 1.0) == amr::Flag::IncreaseResolution);
  CHECK(evaluate(serialized, 0.7, 1.0) == amr::Flag::IncreaseResolution);
  CHECK(evaluate(p_only, 0.2, 1.0) == amr::Flag::IncreaseResolution);
  CHECK(evaluate(parsed, 1.3, 1.0) == amr::Flag::IncreaseResolution);
  CHECK(evaluate(parsed, 0.7, 1.e-12) == amr::Flag::DoNothing);
  CHECK(evaluate(parsed, 0.7, 0.0) == amr::Flag::DoNothing);
  // Nondecaying spectra must not chase negligible far tails. Above the
  // amplitude guard, they still request refinement.
  CHECK(evaluate(parsed, -0.1, 1.e-12) == amr::Flag::DoNothing);
  CHECK(evaluate(parsed, -0.1, 1.e-6) == amr::Flag::Split);
  // The evolution criterion must pass its absolute budget to the
  // multidimensional estimator, while retaining resolved polynomial data.
  const Mesh<3> cube_mesh{8, Spectral::Basis::Legendre,
                          Spectral::Quadrature::GaussLobatto};
  ModalVector cube_modes(cube_mesh.number_of_grid_points(), 0.0);
  cube_modes[0] = 1.0;
  cube_modes[1] = 0.5;
  cube_modes[2] = 0.1;
  cube_modes[7 + 8 * 7] = 1.e-10;
  const auto cube_box = db::create<tmpl::list<domain::Tags::Mesh<3>, Field>>(
      cube_mesh,
      Scalar<DataVector>{to_nodal_coefficients(cube_modes, cube_mesh)});
  const amr::Criteria::Chalmers<3, tmpl::list<Field>> cube_criterion{
      {"Field"}, 1.e-7, 0.0, 0.5, false};
  CHECK(cube_criterion(cube_box, cache, ElementId<3>{0}) ==
        std::array{amr::Flag::DoNothing, amr::Flag::DoNothing,
                   amr::Flag::DoNothing});
  const Mesh<3> transverse_mesh{{{8, 3, 3}},
                                Spectral::Basis::Legendre,
                                Spectral::Quadrature::GaussLobatto};
  ModalVector transverse_modes(transverse_mesh.number_of_grid_points(), 0.0);
  transverse_modes[0] = 1.0;
  transverse_modes[8 * 2] = 1.e-10;
  const auto transverse_box =
      db::create<tmpl::list<domain::Tags::Mesh<3>, Field>>(
          transverse_mesh, Scalar<DataVector>{to_nodal_coefficients(
                               transverse_modes, transverse_mesh)});
  const amr::Criteria::Chalmers<3, tmpl::list<Field>> relative_criterion{
      {"Field"}, 0.0, 1.e-7, 0.5, false};
  CHECK(relative_criterion(transverse_box, cache, ElementId<3>{0}) ==
        std::array{amr::Flag::DoNothing, amr::Flag::DoNothing,
                   amr::Flag::DoNothing});
  using CubeCriterion = amr::Criteria::Chalmers<3, tmpl::list<Field>>;
  const auto energy_criterion = TestHelpers::test_creation<CubeCriterion>(
      "VariablesToMonitor: [Field]\n"
      "AbsoluteTarget: 1.e-4\n"
      "RelativeTarget: 0.0\n"
      "SmoothnessThreshold: 1.0\n"
      "PRefineOnly: false\n"
      "UseDirectionalEnergy: true\n"
      "CoarseningFactor: 0.0\n"
      "HCoarsening: false\n");
  const CubeCriterion stripe_criterion{{"Field"}, 1.e-4, 0.0,
                                       1.0,       false, false};
  const Mesh<3> mixed_mesh{{{8, 7, 3}},
                           Spectral::Basis::Legendre,
                           Spectral::Quadrature::GaussLobatto};
  ModalVector mixed_modes(mixed_mesh.number_of_grid_points(), 0.0);
  for (size_t k = 0; k < 8; ++k) {
    mixed_modes[k] = 2.0 * std::exp(-1.2 * k);
    mixed_modes[k + 8] = 1.e-4 * std::exp(-0.1 * k);
  }
  const auto mixed_box = db::create<tmpl::list<domain::Tags::Mesh<3>, Field>>(
      mixed_mesh,
      Scalar<DataVector>{to_nodal_coefficients(mixed_modes, mixed_mesh)});
  CHECK(energy_criterion(mixed_box, cache, ElementId<3>{0})[0] ==
        amr::Flag::IncreaseResolution);
  CHECK(serialize_and_deserialize(energy_criterion)(mixed_box, cache,
                                                    ElementId<3>{0})[0] ==
        amr::Flag::IncreaseResolution);
  CHECK(stripe_criterion(mixed_box, cache, ElementId<3>{0})[0] ==
        amr::Flag::Split);
  CHECK(serialize_and_deserialize(stripe_criterion)(
            mixed_box, cache, ElementId<3>{0})[0] == amr::Flag::Split);
  const auto coarsening = TestHelpers::test_creation<Criterion>(
      "VariablesToMonitor: [Field]\n"
      "AbsoluteTarget: 1.e-3\n"
      "RelativeTarget: 0.0\n"
      "SmoothnessThreshold: 0.5\n"
      "PRefineOnly: false\n"
      "CoarseningFactor: 0.1\n"
      "HCoarsening: false\n");
  CHECK(evaluate(coarsening, 1.3, 0.1) == amr::Flag::DecreaseResolution);
  // Between the two thresholds, keep the resolution: no immediate reversal.
  CHECK(evaluate(coarsening, 1.3, 1.0) == amr::Flag::DoNothing);
  CHECK(evaluate(coarsening, 0.7, 1.0) == amr::Flag::IncreaseResolution);
  CHECK(evaluate(coarsening, 0.2, 1.0) == amr::Flag::Split);
  CHECK(evaluate(coarsening, 0.2, 0.0) == amr::Flag::DecreaseResolution);
  CHECK(evaluate(serialize_and_deserialize(coarsening), 1.3, 0.1) ==
        amr::Flag::DecreaseResolution);
  // A resolved polynomial must retain enough modes to certify its tail.
  // A second monitored field can veto the first field's coarsening request.
  using TwoFields = amr::Criteria::Chalmers<1, tmpl::list<Field, OtherField>>;
  const TwoFields two_fields{
      {"Field", "OtherField"}, 1.e-3, 0.0, 0.5, false, 0.1};
  ModalVector polynomial(8, 0.0);
  polynomial[0] = 1.0;
  polynomial[5] = 0.2;
  auto two_box =
      db::create<tmpl::list<domain::Tags::Mesh<1>, Field, OtherField>>(
          mesh, Scalar<DataVector>{DataVector(8, 0.0)},
          Scalar<DataVector>{to_nodal_coefficients(polynomial, mesh)});
  CHECK(two_fields(two_box, cache, ElementId<1>{0})[0] == amr::Flag::DoNothing);
  db::mutate<OtherField>([](const auto field) { get(*field) = 0.0; },
                         make_not_null(&two_box));
  CHECK(two_fields(two_box, cache, ElementId<1>{0})[0] ==
        amr::Flag::DecreaseResolution);
  for (const double factor : {-0.1, 1.0, 1.1}) {
    CHECK_THROWS_WITH((Criterion{{"Field"}, 1.e-8, 0.0, 0.5, false, factor}),
                      Catch::Matchers::ContainsSubstring("CoarseningFactor"));
  }
  CHECK_THROWS_WITH((Criterion{{"Field"}, 1.e-8, 0.0, 0.0, false}),
                    Catch::Matchers::ContainsSubstring("finite and positive"));
}
