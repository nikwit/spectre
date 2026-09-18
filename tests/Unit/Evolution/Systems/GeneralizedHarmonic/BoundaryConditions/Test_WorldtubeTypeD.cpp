// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace {
using frame = Frame::Inertial;
constexpr size_t Dim = 3;
using Imposition = gh::BoundaryConditions::detail::SectorImposition;
using Worldtube = gh::BoundaryConditions::WorldtubeTypeD<Dim>;

struct Metavariables {
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<gh::BoundaryConditions::BoundaryCondition<Dim>,
                             tmpl::list<Worldtube>>>;
  };
};

// The face data of a boundary condition call. Only the relations that the
// corrections rely on are enforced: the spacetime metric, its inverse, the
// normal vector and the lapse and shift belong together, and the interface
// normal is a unit covector. Everything else is random.
struct FaceData {
  tnsr::i<DataVector, Dim, frame> normal_covector;
  tnsr::I<DataVector, Dim, frame> normal_vector;
  tnsr::aa<DataVector, Dim, frame> spacetime_metric;
  tnsr::aa<DataVector, Dim, frame> pi;
  tnsr::iaa<DataVector, Dim, frame> phi;
  tnsr::I<DataVector, Dim, frame> coords;
  Scalar<DataVector> gamma1;
  Scalar<DataVector> gamma2;
  Scalar<DataVector> lapse;
  tnsr::I<DataVector, Dim, frame> shift;
  tnsr::AA<DataVector, Dim, frame> inverse_spacetime_metric;
  tnsr::A<DataVector, Dim, frame> spacetime_unit_normal_vector;
  tnsr::iaa<DataVector, Dim, frame> three_index_constraint;
  tnsr::a<DataVector, Dim, frame> gauge_source;
  tnsr::ab<DataVector, Dim, frame> spacetime_deriv_gauge_source;
  tnsr::aa<DataVector, Dim, frame> dt_spacetime_metric;
  tnsr::aa<DataVector, Dim, frame> dt_pi;
  tnsr::iaa<DataVector, Dim, frame> dt_phi;
  tnsr::iaa<DataVector, Dim, frame> d_spacetime_metric;
  tnsr::iaa<DataVector, Dim, frame> d_pi;
  tnsr::ijaa<DataVector, Dim, frame> d_phi;
};

struct Corrections {
  tnsr::aa<DataVector, Dim, frame> dt_spacetime_metric;
  tnsr::aa<DataVector, Dim, frame> dt_pi;
  tnsr::iaa<DataVector, Dim, frame> dt_phi;
};

// `shift_along_normal` sets the normal component of the shift. Above the
// lapse (about one) every characteristic field is incoming at every point, so
// no correction is zeroed by the characteristic-speed test and the identities
// below hold exactly.
template <typename Generator>
FaceData make_face_data(const gsl::not_null<Generator*> generator,
                        const size_t num_points,
                        const double shift_along_normal) {
  std::uniform_real_distribution<> small(-0.1, 0.1);
  std::uniform_real_distribution<> positive(0.2, 1.0);
  std::uniform_real_distribution<> radius(1.0, 2.0);
  const DataVector used_for_size(num_points);
  const auto random = [&generator, &small, &used_for_size](auto tensor_type) {
    using T = decltype(tensor_type);
    return make_with_random_values<T>(generator, make_not_null(&small),
                                      used_for_size);
  };

  FaceData data{};
  // Spatial metric: flat plus a small perturbation
  auto spatial_metric = random(tnsr::ii<DataVector, Dim, frame>{});
  for (size_t i = 0; i < Dim; ++i) {
    spatial_metric.get(i, i) += 1.;
  }
  const auto inverse_spatial_metric =
      determinant_and_inverse(spatial_metric).second;
  data.lapse = make_with_value<Scalar<DataVector>>(used_for_size, 1.);
  get(data.lapse) += make_with_random_values<DataVector>(
      generator, make_not_null(&small), used_for_size);
  data.shift = random(tnsr::I<DataVector, Dim, frame>{});
  get<0>(data.shift) += shift_along_normal;
  data.spacetime_metric =
      gr::spacetime_metric(data.lapse, data.shift, spatial_metric);
  data.inverse_spacetime_metric = gr::inverse_spacetime_metric(
      data.lapse, data.shift, inverse_spatial_metric);
  data.spacetime_unit_normal_vector =
      gr::spacetime_normal_vector(data.lapse, data.shift);

  // Unit normal, mostly along x
  data.normal_covector = random(tnsr::i<DataVector, Dim, frame>{});
  get<0>(data.normal_covector) += 1.;
  DataVector norm(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      norm += inverse_spatial_metric.get(i, j) * data.normal_covector.get(i) *
              data.normal_covector.get(j);
    }
  }
  norm = sqrt(norm);
  for (size_t i = 0; i < Dim; ++i) {
    data.normal_covector.get(i) /= norm;
  }
  data.normal_vector =
      make_with_value<tnsr::I<DataVector, Dim, frame>>(used_for_size, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      data.normal_vector.get(i) +=
          inverse_spatial_metric.get(i, j) * data.normal_covector.get(j);
    }
  }

  data.pi = random(tnsr::aa<DataVector, Dim, frame>{});
  data.phi = random(tnsr::iaa<DataVector, Dim, frame>{});
  data.coords = make_with_random_values<tnsr::I<DataVector, Dim, frame>>(
      generator, make_not_null(&radius), used_for_size);
  data.gamma1 = make_with_random_values<Scalar<DataVector>>(
      generator, make_not_null(&positive), used_for_size);
  data.gamma2 = make_with_random_values<Scalar<DataVector>>(
      generator, make_not_null(&positive), used_for_size);
  data.three_index_constraint = random(tnsr::iaa<DataVector, Dim, frame>{});
  data.gauge_source = random(tnsr::a<DataVector, Dim, frame>{});
  data.spacetime_deriv_gauge_source =
      random(tnsr::ab<DataVector, Dim, frame>{});
  data.dt_spacetime_metric = random(tnsr::aa<DataVector, Dim, frame>{});
  data.dt_pi = random(tnsr::aa<DataVector, Dim, frame>{});
  data.dt_phi = random(tnsr::iaa<DataVector, Dim, frame>{});
  data.d_spacetime_metric = random(tnsr::iaa<DataVector, Dim, frame>{});
  data.d_pi = random(tnsr::iaa<DataVector, Dim, frame>{});
  data.d_phi = random(tnsr::ijaa<DataVector, Dim, frame>{});
  return data;
}

Corrections make_corrections(const size_t num_points) {
  const DataVector used_for_size(num_points);
  return {
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(used_for_size, 0.),
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(used_for_size, 0.),
      make_with_value<tnsr::iaa<DataVector, Dim, frame>>(used_for_size, 0.)};
}

Corrections apply_worldtube(
    const Worldtube& boundary_condition, const FaceData& data,
    const double time, const Domain<Dim>& domain, const Element<Dim>& element,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  auto corrections = make_corrections(get_size(get(data.lapse)));
  const auto error = boundary_condition.dg_time_derivative(
      make_not_null(&corrections.dt_spacetime_metric),
      make_not_null(&corrections.dt_pi), make_not_null(&corrections.dt_phi),
      std::nullopt, data.normal_covector, data.normal_vector,
      data.spacetime_metric, data.pi, data.phi, data.coords, data.gamma1,
      data.gamma2, data.lapse, data.shift, data.inverse_spacetime_metric,
      data.spacetime_unit_normal_vector, data.three_index_constraint,
      data.gauge_source, data.spacetime_deriv_gauge_source,
      data.dt_spacetime_metric, data.dt_pi, data.dt_phi,
      data.d_spacetime_metric, data.d_pi, data.d_phi, time, domain, element,
      functions_of_time);
  CHECK_FALSE(error.has_value());
  return corrections;
}

Corrections apply_constraint_preserving_bjorhus(const FaceData& data,
                                                const double time) {
  const gh::BoundaryConditions::ConstraintPreservingBjorhus<Dim>
      boundary_condition{
          gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType::
              ConstraintPreservingPhysical};
  auto corrections = make_corrections(get_size(get(data.lapse)));
  const auto error = boundary_condition.dg_time_derivative(
      make_not_null(&corrections.dt_spacetime_metric),
      make_not_null(&corrections.dt_pi), make_not_null(&corrections.dt_phi),
      std::nullopt, data.normal_covector, data.normal_vector,
      data.spacetime_metric, data.pi, data.phi, data.coords, data.gamma1,
      data.gamma2, data.lapse, data.shift, data.inverse_spacetime_metric,
      data.spacetime_unit_normal_vector, data.three_index_constraint,
      data.gauge_source, data.spacetime_deriv_gauge_source,
      data.dt_spacetime_metric, data.dt_pi, data.dt_phi,
      data.d_spacetime_metric, data.d_pi, data.d_phi, time);
  CHECK_FALSE(error.has_value());
  return corrections;
}

template <typename T>
void add_to(const gsl::not_null<T*> result, const T& summand,
            const double factor = 1.) {
  for (size_t i = 0; i < result->size(); ++i) {
    (*result)[i] += factor * summand[i];
  }
}

void check_corrections_equal(const Corrections& lhs, const Corrections& rhs) {
  // The two boundary conditions assemble the same total through differently
  // ordered sums of the three sector projections, so allow roundoff.
  Approx custom_approx = Approx::custom().epsilon(1.e-11).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(lhs.dt_spacetime_metric, rhs.dt_spacetime_metric,
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(lhs.dt_pi, rhs.dt_pi, custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(lhs.dt_phi, rhs.dt_phi, custom_approx);
}

// A domain whose only excision sphere is centered at the origin, and an
// element on its boundary
struct ExcisedSphere {
  ExcisedSphere()
      : domain(domain::creators::Sphere{
            1.0, 3.0, domain::creators::Sphere::Excision{}, 0_st, 3_st, true}
                   .create_domain()),
        element(ElementId<Dim>{0}, {}) {
    REQUIRE(domain.excision_spheres().size() == 1);
    REQUIRE(domain.excision_spheres()
                .begin()
                ->second.abutting_direction(element.id())
                .has_value());
  }
  Domain<Dim> domain;
  Element<Dim> element;
  domain::FunctionsOfTimeMap functions_of_time{};
};

void test_option_parsing_and_serialization() {
  {
    const auto created = TestHelpers::test_creation<
        std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
        Metavariables>(
        "WorldtubeTypeD:\n"
        "  ConstraintPreservingSector: Bjorhus\n"
        "  PhysicalSector: Frozen\n"
        "  GaugeSector: SommerfeldAbsorbing");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->constraint_v_psi() == Imposition::Bjorhus);
    CHECK(worldtube->constraint_v_zero() == Imposition::Bjorhus);
    CHECK(worldtube->constraint_v_minus() == Imposition::Bjorhus);
    CHECK(worldtube->physical_sector() == Imposition::Frozen);
    CHECK(worldtube->gauge_sector() == Imposition::SommerfeldAbsorbing);
    CHECK(*worldtube == Worldtube{Imposition::Bjorhus, Imposition::Frozen,
                                  Imposition::SommerfeldAbsorbing});
    CHECK(*worldtube != Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::SommerfeldAbsorbing});
    const auto cloned = worldtube->get_clone();
    CHECK(*dynamic_cast<const Worldtube*>(cloned.get()) == *worldtube);
  }
  {
    const auto created = TestHelpers::test_creation<
        std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
        Metavariables>(
        "WorldtubeTypeD:\n"
        "  ConstraintPreservingSector:\n"
        "    VPsi: Frozen\n"
        "    VZero: Bjorhus\n"
        "    VMinus: Frozen\n"
        "  PhysicalSector: Bjorhus\n"
        "  GaugeSector: Frozen");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->constraint_v_psi() == Imposition::Frozen);
    CHECK(worldtube->constraint_v_zero() == Imposition::Bjorhus);
    CHECK(worldtube->constraint_v_minus() == Imposition::Frozen);
    CHECK(worldtube->physical_sector() == Imposition::Bjorhus);
    CHECK(worldtube->gauge_sector() == Imposition::Frozen);
  }
  test_serialization_via_base<gh::BoundaryConditions::BoundaryCondition<Dim>,
                              Worldtube>(
      gh::BoundaryConditions::detail::PerFieldConstraintSectors{
          Imposition::Frozen, Imposition::Bjorhus, Imposition::Bjorhus},
      Imposition::Frozen, Imposition::SommerfeldOutgoing);

  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Bjorhus")),
      Catch::Matchers::ContainsSubstring("GaugeSector: Bjorhus is not"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: SommerfeldAbsorbing\n"
                         "  GaugeSector: Frozen")),
      Catch::Matchers::ContainsSubstring(
          "PhysicalSector: the Sommerfeld conditions apply only"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Ghost\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen")),
      Catch::Matchers::ContainsSubstring("Must be one of Bjorhus, Frozen"));
}

// With every sector imposed by its Bjorhus correction and the outgoing
// Sommerfeld sign, the worldtube condition on an excision centered at the
// origin is the outer-boundary ConstraintPreservingPhysical condition.
void test_reproduces_constraint_preserving_bjorhus() {
  MAKE_GENERATOR(generator);
  const ExcisedSphere sphere{};
  const Worldtube worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                            Imposition::SommerfeldOutgoing};
  for (const double shift_along_normal : {0.0, 2.0}) {
    CAPTURE(shift_along_normal);
    const auto data =
        make_face_data(make_not_null(&generator), 5, shift_along_normal);
    check_corrections_equal(
        apply_worldtube(worldtube, data, 1.3, sphere.domain, sphere.element,
                        sphere.functions_of_time),
        apply_constraint_preserving_bjorhus(data, 1.3));
  }
}

// With every sector frozen, and every characteristic field incoming, the
// corrections cancel the time derivatives exactly.
void test_all_frozen_cancels_time_derivatives() {
  MAKE_GENERATOR(generator);
  const ExcisedSphere sphere{};
  const Worldtube worldtube{Imposition::Frozen, Imposition::Frozen,
                            Imposition::Frozen};
  const auto data = make_face_data(make_not_null(&generator), 5, 2.0);
  const auto corrections =
      apply_worldtube(worldtube, data, 0.0, sphere.domain, sphere.element,
                      sphere.functions_of_time);
  Approx custom_approx = Approx::custom().epsilon(1.e-11).scale(1.);
  auto expected = make_corrections(5);
  add_to(make_not_null(&expected.dt_spacetime_metric), data.dt_spacetime_metric,
         -1.);
  add_to(make_not_null(&expected.dt_pi), data.dt_pi, -1.);
  add_to(make_not_null(&expected.dt_phi), data.dt_phi, -1.);
  CHECK_ITERABLE_CUSTOM_APPROX(corrections.dt_spacetime_metric,
                               expected.dt_spacetime_metric, custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(corrections.dt_pi, expected.dt_pi,
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(corrections.dt_phi, expected.dt_phi,
                               custom_approx);
}

// Switching a sector on adds that sector's terms and nothing else: the
// corrections are additive over the sectors relative to the all-frozen
// baseline, for both Sommerfeld signs.
void test_sectors_are_independent() {
  MAKE_GENERATOR(generator);
  const ExcisedSphere sphere{};
  const auto data = make_face_data(make_not_null(&generator), 5, 2.0);
  const auto apply = [&](const Imposition constraint, const Imposition physical,
                         const Imposition gauge) {
    return apply_worldtube(Worldtube{constraint, physical, gauge}, data, 0.7,
                           sphere.domain, sphere.element,
                           sphere.functions_of_time);
  };
  const auto frozen =
      apply(Imposition::Frozen, Imposition::Frozen, Imposition::Frozen);
  for (const auto gauge :
       {Imposition::SommerfeldAbsorbing, Imposition::SommerfeldOutgoing}) {
    CAPTURE(gauge);
    const auto all = apply(Imposition::Bjorhus, Imposition::Bjorhus, gauge);
    const auto only_constraint =
        apply(Imposition::Bjorhus, Imposition::Frozen, Imposition::Frozen);
    const auto only_physical =
        apply(Imposition::Frozen, Imposition::Bjorhus, Imposition::Frozen);
    const auto only_gauge =
        apply(Imposition::Frozen, Imposition::Frozen, gauge);

    // all - frozen == sum over sectors of (only_sector - frozen)
    Corrections expected = frozen;
    for (const auto* const only :
         {&only_constraint, &only_physical, &only_gauge}) {
      add_to(make_not_null(&expected.dt_spacetime_metric),
             only->dt_spacetime_metric);
      add_to(make_not_null(&expected.dt_pi), only->dt_pi);
      add_to(make_not_null(&expected.dt_phi), only->dt_phi);
      add_to(make_not_null(&expected.dt_spacetime_metric),
             frozen.dt_spacetime_metric, -1.);
      add_to(make_not_null(&expected.dt_pi), frozen.dt_pi, -1.);
      add_to(make_not_null(&expected.dt_phi), frozen.dt_phi, -1.);
    }
    check_corrections_equal(all, expected);
  }
  // The two Sommerfeld signs differ
  const auto absorbing = apply(Imposition::Frozen, Imposition::Frozen,
                               Imposition::SommerfeldAbsorbing);
  const auto outgoing = apply(Imposition::Frozen, Imposition::Frozen,
                              Imposition::SommerfeldOutgoing);
  CHECK(max(abs(get<0, 0>(absorbing.dt_pi) - get<0, 0>(outgoing.dt_pi))) >
        1.e-8);
}

void test_excision_sphere_center() {
  const ElementId<Dim> abutting_element{2};
  const ElementId<Dim> other_element{5};
  domain::FunctionsOfTimeMap functions_of_time{};
  functions_of_time["Translation"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0,
          std::array<DataVector, 3>{
              {{1.0, -2.0, 0.5}, {0.0, 3.0, 0.0}, {0.0, 0.0, 0.0}}},
          10.0);

  std::unordered_map<std::string, ExcisionSphere<Dim>> excision_spheres{};
  excision_spheres.emplace(
      "Static",
      ExcisionSphere<Dim>{1.5,
                          tnsr::I<double, Dim, Frame::Grid>{{{1., 2., 3.}}},
                          {{2, Direction<Dim>::lower_zeta()}}});
  ExcisionSphere<Dim> moving{0.5,
                             tnsr::I<double, Dim, Frame::Grid>{{{-1., 0., 4.}}},
                             {{7, Direction<Dim>::lower_zeta()}}};
  moving.inject_time_dependent_maps(
      std::make_unique<domain::CoordinateMap<
          Frame::Grid, Frame::Inertial,
          domain::CoordinateMaps::TimeDependent::Translation<Dim>>>(
          domain::CoordinateMaps::TimeDependent::Translation<Dim>{
              "Translation"}));
  excision_spheres.emplace("Moving", std::move(moving));

  const auto static_center =
      gh::BoundaryConditions::detail::excision_sphere_center(
          excision_spheres, abutting_element, 2.0, functions_of_time);
  CHECK(static_center == tnsr::I<double, Dim, frame>{{{1., 2., 3.}}});

  const auto moving_center =
      gh::BoundaryConditions::detail::excision_sphere_center(
          excision_spheres, ElementId<Dim>{7}, 2.0, functions_of_time);
  // Center translated by (1 + 0 * 2, -2 + 3 * 2, 0.5) at t = 2
  CHECK_ITERABLE_APPROX(moving_center,
                        (tnsr::I<double, Dim, frame>{{{0., 4., 4.5}}}));

  CHECK_THROWS_WITH(
      gh::BoundaryConditions::detail::excision_sphere_center(
          excision_spheres, other_element, 2.0, functions_of_time),
      Catch::Matchers::ContainsSubstring(
          "abuts none of the domain's excision spheres"));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.BoundaryConditions.Worldtube",
    "[Unit][Evolution]") {
  test_option_parsing_and_serialization();
  test_reproduces_constraint_preserving_bjorhus();
  test_all_frozen_cancels_time_derivatives();
  test_sectors_are_independent();
  test_excision_sphere_center();
}
