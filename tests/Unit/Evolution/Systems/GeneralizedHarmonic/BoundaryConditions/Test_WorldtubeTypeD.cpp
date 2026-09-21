// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
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
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/CovariantDerivOfExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Phi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Pi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/TimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylPropagating.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace {
using frame = Frame::Inertial;
constexpr size_t Dim = 3;
using Imposition = gh::BoundaryConditions::detail::SectorImposition;
using Model = gh::BoundaryConditions::detail::PhysicalModel;
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
    const domain::FunctionsOfTimeMap& functions_of_time,
    const gh::worldtube::KretschmannFaceData<Dim>& face_data = {}) {
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
      functions_of_time, face_data);
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
        "  GaugeSector: SommerfeldAbsorbing\n"
        "  PhysicalModel: None\n"
        "  Mass: None\n"
        "  MomentRelaxationTime: None");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->constraint_v_psi() == Imposition::Bjorhus);
    CHECK(worldtube->constraint_v_zero() == Imposition::Bjorhus);
    CHECK(worldtube->constraint_v_minus() == Imposition::Bjorhus);
    CHECK(worldtube->physical_sector() == Imposition::Frozen);
    CHECK(worldtube->gauge_sector() == Imposition::SommerfeldAbsorbing);
    CHECK(*worldtube == Worldtube{Imposition::Bjorhus, Imposition::Frozen,
                                  Imposition::SommerfeldAbsorbing,
                                  Model::None});
    CHECK(*worldtube != Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::SommerfeldAbsorbing,
                                  Model::None});
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
        "  GaugeSector: Frozen\n"
        "  PhysicalModel: None\n"
        "  Mass: None\n"
        "  MomentRelaxationTime: None");
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
      Imposition::Frozen, Imposition::SommerfeldOutgoing, Model::None,
      std::nullopt);
  test_serialization_via_base<gh::BoundaryConditions::BoundaryCondition<Dim>,
                              Worldtube>(
      Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
      Model::Quadrupole, std::optional<double>{1.0});
  test_serialization_via_base<gh::BoundaryConditions::BoundaryCondition<Dim>,
                              Worldtube>(
      Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
      Model::Quadrupole, std::optional<double>{1.0},
      std::optional<double>{10.0});
  {
    const auto created = TestHelpers::test_creation<
        std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
        Metavariables>(
        "WorldtubeTypeD:\n"
        "  ConstraintPreservingSector: Bjorhus\n"
        "  PhysicalSector: Bjorhus\n"
        "  GaugeSector: Frozen\n"
        "  PhysicalModel: Quadrupole\n"
        "  Mass: 1.0\n"
        "  MomentRelaxationTime: 10.0");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->moment_relaxation_time() == std::optional<double>{10.0});
    CHECK(*worldtube == Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::Frozen, Model::Quadrupole, 1.0,
                                  10.0});
    CHECK(*worldtube != Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::Frozen, Model::Quadrupole, 1.0});
  }
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: TypeD\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: 10.0")),
      Catch::Matchers::ContainsSubstring(
          "MomentRelaxationTime is only used by"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: Quadrupole\n"
                         "  Mass: 1.0\n"
                         "  MomentRelaxationTime: 0.0")),
      Catch::Matchers::ContainsSubstring(
          "MomentRelaxationTime must be positive"));

  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Bjorhus\n"
                         "  PhysicalModel: None\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("GaugeSector: Bjorhus is not"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: SommerfeldAbsorbing\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: None\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring(
          "PhysicalSector: the Sommerfeld conditions apply only"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Ghost\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: None\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("Must be one of Bjorhus, Frozen"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Frozen\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: TypeD\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("Use PhysicalSector: Bjorhus"));
  {
    const auto created = TestHelpers::test_creation<
        std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
        Metavariables>(
        "WorldtubeTypeD:\n"
        "  ConstraintPreservingSector: Bjorhus\n"
        "  PhysicalSector: Bjorhus\n"
        "  GaugeSector: SommerfeldAbsorbing\n"
        "  PhysicalModel: TypeD\n"
        "  Mass: None\n"
        "  MomentRelaxationTime: None");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->physical_model() == Model::TypeD);
    CHECK_FALSE(worldtube->mass().has_value());
  }
  {
    const auto created = TestHelpers::test_creation<
        std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
        Metavariables>(
        "WorldtubeTypeD:\n"
        "  ConstraintPreservingSector: Bjorhus\n"
        "  PhysicalSector: Bjorhus\n"
        "  GaugeSector: SommerfeldAbsorbing\n"
        "  PhysicalModel: Quadrupole\n"
        "  Mass: 1.0\n"
        "  MomentRelaxationTime: None");
    const auto* const worldtube = dynamic_cast<const Worldtube*>(created.get());
    REQUIRE(worldtube != nullptr);
    CHECK(worldtube->physical_model() == Model::Quadrupole);
    CHECK(worldtube->mass() == std::optional<double>{1.0});
    CHECK(*worldtube == Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::SommerfeldAbsorbing,
                                  Model::Quadrupole, 1.0});
    CHECK(*worldtube != Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                  Imposition::SommerfeldAbsorbing,
                                  Model::Quadrupole, 2.0});
  }
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: Quadrupole\n"
                         "  Mass: None\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("Quadrupole needs the mass"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: TypeD\n"
                         "  Mass: 1.0\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("Mass is only used by"));
  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<
          std::unique_ptr<gh::BoundaryConditions::BoundaryCondition<Dim>>,
          Metavariables>("WorldtubeTypeD:\n"
                         "  ConstraintPreservingSector: Bjorhus\n"
                         "  PhysicalSector: Bjorhus\n"
                         "  GaugeSector: Frozen\n"
                         "  PhysicalModel: Quadrupole\n"
                         "  Mass: -1.0\n"
                         "  MomentRelaxationTime: None")),
      Catch::Matchers::ContainsSubstring("Mass must be positive"));
}

// With every sector imposed by its Bjorhus correction and the outgoing
// Sommerfeld sign, the worldtube condition on an excision centered at the
// origin is the outer-boundary ConstraintPreservingPhysical condition.
void test_reproduces_constraint_preserving_bjorhus() {
  MAKE_GENERATOR(generator);
  const ExcisedSphere sphere{};
  const Worldtube worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                            Imposition::SommerfeldOutgoing, Model::None};
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
                            Imposition::Frozen, Model::None};
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
    return apply_worldtube(Worldtube{constraint, physical, gauge, Model::None},
                           data, 0.7, sphere.domain, sphere.element,
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

// The face normal of the boundary condition points out of the domain, into
// the hole. With the worldtube-adapted tetrad built along that normal, the
// incoming Weyl mode of Kidder-Scheel-Teukolsky is half the note's w^- built
// from Psi0 -- and would be Psi4's field with the opposite orientation, which
// is the orientation the offline study used. Checked on a flat slice against
// gr::weyl_propagating with a constraint-consistent covariant derivative of
// the extrinsic curvature, for which the two forms of the mode agree.
void test_incoming_mode_normalization() {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> dist(-1., 1.);
  const size_t num_points = 6;
  const DataVector used_for_size(num_points);
  const auto random_symmetric = [&]() {
    return make_with_random_values<tnsr::ii<DataVector, Dim, frame>>(
        make_not_null(&generator), make_not_null(&dist), used_for_size);
  };
  const auto ricci = random_symmetric();
  const auto extrinsic_curvature = random_symmetric();
  // nabla_k K_ij with the momentum constraint D^j K_ij = D_i K enforced:
  // add delta_ij u_k with u = v/2, where v_i is the constraint defect
  auto cov_deriv_k = make_with_random_values<tnsr::ijj<DataVector, Dim, frame>>(
      make_not_null(&generator), make_not_null(&dist), used_for_size);
  for (size_t i = 0; i < Dim; ++i) {
    DataVector defect(num_points, 0.);
    for (size_t j = 0; j < Dim; ++j) {
      defect += cov_deriv_k.get(j, i, j) - cov_deriv_k.get(i, j, j);
    }
    for (size_t j = 0; j < Dim; ++j) {
      cov_deriv_k.get(i, j, j) += 0.5 * defect;
    }
  }
  tnsr::ii<DataVector, Dim, frame> spatial_metric(num_points, 0.);
  tnsr::II<DataVector, Dim, frame> inverse_spatial_metric(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    spatial_metric.get(i, i) = 1.;
    inverse_spatial_metric.get(i, i) = 1.;
  }
  auto normal_covector =
      make_with_random_values<tnsr::i<DataVector, Dim, frame>>(
          make_not_null(&generator), make_not_null(&dist), used_for_size);
  const DataVector norm =
      sqrt(square(get<0>(normal_covector)) + square(get<1>(normal_covector)) +
           square(get<2>(normal_covector)));
  tnsr::I<DataVector, Dim, frame> normal_vector(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    normal_covector.get(i) /= norm;
    normal_vector.get(i) = normal_covector.get(i);
  }
  const auto projection_IJ =
      gr::transverse_projection_operator(inverse_spatial_metric, normal_vector);
  const auto projection_ij =
      gr::transverse_projection_operator(spatial_metric, normal_covector);
  const auto projection_Ij =
      gr::transverse_projection_operator(normal_vector, normal_covector);
  const auto u8_minus = gr::weyl_propagating(
      ricci, extrinsic_curvature, inverse_spatial_metric, cov_deriv_k,
      normal_vector, projection_IJ, projection_ij, projection_Ij, -1.);

  const auto electric =
      gr::weyl_electric(ricci, extrinsic_curvature, inverse_spatial_metric);
  const auto magnetic = gr::weyl_magnetic(cov_deriv_k, spatial_metric,
                                          Scalar<DataVector>(num_points, 1.));

  // Tetrad along the face normal: U8- = w^-(Psi0) / 2
  const gr::np::WeylScalars psi = gr::np::weyl_scalars_from_electric_magnetic(
      electric, magnetic, spatial_metric, normal_vector);
  auto half_w_minus = gr::np::incoming_weyl_field(
      Scalar<ComplexDataVector>{psi.get(0)},
      gr::np::adapted_triad(spatial_metric, normal_vector));
  for (auto& component : half_w_minus) {
    component *= 0.5;
  }
  Approx custom_approx = Approx::custom().epsilon(1.e-11).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(u8_minus, half_w_minus, custom_approx);

  // Tetrad along -normal (the study's orientation): the same mode is the
  // field of Psi4, conjugated
  tnsr::I<DataVector, Dim, frame> flipped(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    flipped.get(i) = -normal_vector.get(i);
  }
  const gr::np::WeylScalars psi_flipped =
      gr::np::weyl_scalars_from_electric_magnetic(electric, magnetic,
                                                  spatial_metric, flipped);
  auto half_w_flipped = gr::np::incoming_weyl_field(
      Scalar<ComplexDataVector>{conj(psi_flipped.get(4))},
      gr::np::adapted_triad(spatial_metric, flipped));
  for (auto& component : half_w_flipped) {
    component *= 0.5;
  }
  CHECK_ITERABLE_CUSTOM_APPROX(u8_minus, half_w_flipped, custom_approx);
}

// Face data of a boosted Kerr-Schild hole on a coordinate sphere around its
// center, with the derivatives of Phi and Pi by fourth-order finite
// differences of the analytic solution. Fields the tests do not compare are
// random.
template <typename Generator>
FaceData kerr_schild_face_data(const gsl::not_null<Generator*> generator,
                               const size_t num_points) {
  const gr::Solutions::KerrSchild solution{
      1., {{0., 0., 0.}}, {{0., 0., 0.}}, {{0.2, -0.1, 0.15}}};
  const double radius = 4.;
  const double step = 1.e-3;
  std::uniform_real_distribution<> unit(-1., 1.);
  std::uniform_real_distribution<> small(-0.1, 0.1);
  std::uniform_real_distribution<> positive(0.2, 1.0);
  const DataVector used_for_size(num_points);
  auto coords = make_with_random_values<tnsr::I<DataVector, Dim, frame>>(
      generator, make_not_null(&unit), used_for_size);
  const DataVector coord_norm = sqrt(
      square(get<0>(coords)) + square(get<1>(coords)) + square(get<2>(coords)));
  for (size_t i = 0; i < Dim; ++i) {
    coords.get(i) *= radius / coord_norm;
  }

  using tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, ::Tags::dt<gr::Tags::Lapse<DataVector>>,
      gr::Solutions::KerrSchild::DerivLapse<DataVector, frame>,
      gr::Tags::Shift<DataVector, Dim, frame>,
      ::Tags::dt<gr::Tags::Shift<DataVector, Dim, frame>>,
      gr::Solutions::KerrSchild::DerivShift<DataVector, frame>,
      gr::Tags::SpatialMetric<DataVector, Dim, frame>,
      ::Tags::dt<gr::Tags::SpatialMetric<DataVector, Dim, frame>>,
      gr::Solutions::KerrSchild::DerivSpatialMetric<DataVector, frame>,
      gr::Tags::ExtrinsicCurvature<DataVector, Dim, frame>,
      gr::Tags::InverseSpatialMetric<DataVector, Dim, frame>>;
  const auto phi_and_pi_at = [&solution](
                                 const tnsr::I<DataVector, Dim, frame>& x) {
    const auto vars = solution.variables(x, 0., tags{});
    const auto phi = gh::phi(
        get<gr::Tags::Lapse<DataVector>>(vars),
        get<gr::Solutions::KerrSchild::DerivLapse<DataVector, frame>>(vars),
        get<gr::Tags::Shift<DataVector, Dim, frame>>(vars),
        get<gr::Solutions::KerrSchild::DerivShift<DataVector, frame>>(vars),
        get<gr::Tags::SpatialMetric<DataVector, Dim, frame>>(vars),
        get<gr::Solutions::KerrSchild::DerivSpatialMetric<DataVector, frame>>(
            vars));
    const auto pi = gh::pi(
        get<gr::Tags::Lapse<DataVector>>(vars),
        get<::Tags::dt<gr::Tags::Lapse<DataVector>>>(vars),
        get<gr::Tags::Shift<DataVector, Dim, frame>>(vars),
        get<::Tags::dt<gr::Tags::Shift<DataVector, Dim, frame>>>(vars),
        get<gr::Tags::SpatialMetric<DataVector, Dim, frame>>(vars),
        get<::Tags::dt<gr::Tags::SpatialMetric<DataVector, Dim, frame>>>(vars),
        phi);
    return std::pair{phi, pi};
  };

  FaceData data{};
  const auto vars = solution.variables(coords, 0., tags{});
  data.coords = coords;
  data.lapse = get<gr::Tags::Lapse<DataVector>>(vars);
  data.shift = get<gr::Tags::Shift<DataVector, Dim, frame>>(vars);
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataVector, Dim, frame>>(vars);
  const auto& inverse_spatial_metric =
      get<gr::Tags::InverseSpatialMetric<DataVector, Dim, frame>>(vars);
  data.spacetime_metric =
      gr::spacetime_metric(data.lapse, data.shift, spatial_metric);
  data.inverse_spacetime_metric = gr::inverse_spacetime_metric(
      data.lapse, data.shift, inverse_spatial_metric);
  data.spacetime_unit_normal_vector =
      gr::spacetime_normal_vector(data.lapse, data.shift);
  std::tie(data.phi, data.pi) = phi_and_pi_at(coords);
  // Phi_iab = d_i g_ab exactly for the analytic solution
  data.d_spacetime_metric = data.phi;
  data.three_index_constraint =
      make_with_value<tnsr::iaa<DataVector, Dim, frame>>(used_for_size, 0.);
  // Fourth-order central differences for d_k Phi and d_k Pi
  data.d_phi =
      make_with_value<tnsr::ijaa<DataVector, Dim, frame>>(used_for_size, 0.);
  data.d_pi =
      make_with_value<tnsr::iaa<DataVector, Dim, frame>>(used_for_size, 0.);
  for (size_t k = 0; k < Dim; ++k) {
    for (const auto& [multiple, weight] :
         {std::pair{-2., 1. / 12.}, std::pair{-1., -8. / 12.},
          std::pair{1., 8. / 12.}, std::pair{2., -1. / 12.}}) {
      auto shifted = coords;
      shifted.get(k) += multiple * step;
      const auto [phi_shifted, pi_shifted] = phi_and_pi_at(shifted);
      for (size_t a = 0; a <= Dim; ++a) {
        for (size_t b = a; b <= Dim; ++b) {
          data.d_pi.get(k, a, b) += weight / step * pi_shifted.get(a, b);
          for (size_t i = 0; i < Dim; ++i) {
            data.d_phi.get(k, i, a, b) +=
                weight / step * phi_shifted.get(i, a, b);
          }
        }
      }
    }
  }

  // Unit normal into the hole: covector along -x, normalized with gamma
  data.normal_covector =
      make_with_value<tnsr::i<DataVector, Dim, frame>>(used_for_size, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    data.normal_covector.get(i) = -coords.get(i);
  }
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

  // The rest does not enter what the tests compare
  const auto random = [&](auto tensor_type) {
    using T = decltype(tensor_type);
    return make_with_random_values<T>(generator, make_not_null(&small),
                                      used_for_size);
  };
  data.gamma1 = make_with_random_values<Scalar<DataVector>>(
      generator, make_not_null(&positive), used_for_size);
  data.gamma2 = make_with_random_values<Scalar<DataVector>>(
      generator, make_not_null(&positive), used_for_size);
  data.gauge_source = random(tnsr::a<DataVector, Dim, frame>{});
  data.spacetime_deriv_gauge_source =
      random(tnsr::ab<DataVector, Dim, frame>{});
  data.dt_spacetime_metric = random(tnsr::aa<DataVector, Dim, frame>{});
  data.dt_pi = random(tnsr::aa<DataVector, Dim, frame>{});
  data.dt_phi = random(tnsr::iaa<DataVector, Dim, frame>{});
  return data;
}

// A boosted Kerr-Schild hole is exactly type D, so the type-D model must
// reproduce the incoming mode the face data carry, and the boundary condition
// with the model must then differ from the model-free one by exactly the
// projected model mode.
void test_type_d_model_on_boosted_kerr_schild() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 8;
  const auto data =
      kerr_schild_face_data(make_not_null(&generator), num_points);
  tnsr::ii<DataVector, Dim, frame> spatial_metric(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      spatial_metric.get(i, j) = data.spacetime_metric.get(i + 1, j + 1);
    }
  }
  const auto det_and_inverse = determinant_and_inverse(spatial_metric);
  const auto& inverse_spatial_metric = det_and_inverse.second;
  const auto extrinsic_curvature = gh::extrinsic_curvature(
      data.spacetime_unit_normal_vector, data.pi, data.phi);

  tnsr::ii<DataVector, Dim, frame> electric{};
  tnsr::ii<DataVector, Dim, frame> magnetic{};
  gh::worldtube::weyl_electric_magnetic(
      make_not_null(&electric), make_not_null(&magnetic), data.phi, data.d_phi,
      data.d_pi, data.spacetime_unit_normal_vector, spatial_metric,
      inverse_spatial_metric, extrinsic_curvature,
      data.inverse_spacetime_metric);

  // Exact type D: the type-D solve reproduces all five scalars, and the
  // Coulomb scalar is real and -M/r^3 with r the rest-frame areal radius. The
  // lab-frame coordinate sphere of radius 4 is Lorentz-contracted, so that
  // radius lies between 4 and 4 gamma with gamma the boost factor.
  const gr::np::WeylScalars psi = gr::np::weyl_scalars_from_electric_magnetic(
      electric, magnetic, spatial_metric, data.normal_vector);
  const auto coulomb = gr::np::coulomb_scalar(gr::np::invariant_i(psi),
                                              gr::np::invariant_j(psi));
  Approx fd_approx = Approx::custom().epsilon(1.e-7).scale(1.);
  const double lorentz_factor =
      1. / sqrt(1. - square(0.2) - square(0.1) - square(0.15));
  for (size_t p = 0; p < num_points; ++p) {
    CHECK(std::abs(std::imag(get(coulomb)[p])) < 1.e-7);
    CHECK(std::real(get(coulomb)[p]) <
          -1. / 64. / cube(lorentz_factor) + 1.e-7);
    CHECK(std::real(get(coulomb)[p]) > -1. / 64. - 1.e-7);
  }
  const auto rotation = gr::np::solve_type_d_rotation(psi, coulomb);
  CHECK_ITERABLE_CUSTOM_APPROX(rotation.predicted_psi, psi, fd_approx);
  // The boost misaligns the tetrad, so the incoming mode is not zero
  CHECK(max(abs(psi.get(0))) > 1.e-5);

  // The model's incoming mode equals the one the face data carry
  const auto projection_IJ = gr::transverse_projection_operator(
      inverse_spatial_metric, data.normal_vector);
  const auto projection_ij =
      gr::transverse_projection_operator(spatial_metric, data.normal_covector);
  const auto projection_Ij = gr::transverse_projection_operator(
      data.normal_vector, data.normal_covector);
  tnsr::ijj<DataVector, Dim, frame> d_spatial_metric(num_points);
  for (size_t k = 0; k < Dim; ++k) {
    for (size_t i = 0; i < Dim; ++i) {
      for (size_t j = i; j < Dim; ++j) {
        d_spatial_metric.get(k, i, j) = data.phi.get(k, i + 1, j + 1);
      }
    }
  }
  const auto cov_deriv_k = gh::covariant_deriv_of_extrinsic_curvature(
      extrinsic_curvature, data.spacetime_unit_normal_vector,
      raise_or_lower_first_index(gr::christoffel_first_kind(d_spatial_metric),
                                 inverse_spatial_metric),
      data.inverse_spacetime_metric, data.phi, data.d_pi, data.d_phi);
  const auto ricci =
      gh::spatial_ricci_tensor(data.phi, data.d_phi, inverse_spatial_metric);
  const auto u8_minus = gr::weyl_propagating(
      ricci, extrinsic_curvature, inverse_spatial_metric, cov_deriv_k,
      data.normal_vector, projection_IJ, projection_ij, projection_Ij, -1.);
  const auto model_mode = gh::BoundaryConditions::detail::type_d_incoming_mode(
      electric, magnetic, spatial_metric, data.normal_covector);
  CHECK_ITERABLE_CUSTOM_APPROX(model_mode, u8_minus, fd_approx);

  // Through the boundary condition: the model changes the correction by
  // -lambda_- P_TT(2 P^i_a P^j_b U8_ij) relative to no model
  const ExcisedSphere sphere{};
  const auto with_model = apply_worldtube(
      Worldtube{Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
                Model::TypeD},
      data, 0., sphere.domain, sphere.element, sphere.functions_of_time);
  const auto without_model = apply_worldtube(
      Worldtube{Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
                Model::None},
      data, 0., sphere.domain, sphere.element, sphere.functions_of_time);
  gh::BoundaryConditions::Bjorhus::IntermediateVariables<Dim> vars{num_points};
  gh::BoundaryConditions::Bjorhus::compute_intermediate_variables<Dim>(
      make_not_null(&vars), std::nullopt, data.normal_covector,
      data.spacetime_metric, data.pi, data.phi, data.gamma1, data.gamma2,
      data.lapse, data.shift, data.inverse_spacetime_metric,
      data.spacetime_unit_normal_vector, data.three_index_constraint,
      data.gauge_source, data.spacetime_deriv_gauge_source,
      data.dt_spacetime_metric, data.dt_pi, data.dt_phi,
      data.d_spacetime_metric, data.d_pi, data.d_phi);
  auto u3_of_mode =
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(get(data.lapse), 0.);
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = a; b <= Dim; ++b) {
      for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = 0; j < Dim; ++j) {
          u3_of_mode.get(a, b) += 2. * vars.projection_Ab.get(i + 1, a) *
                                  vars.projection_Ab.get(j + 1, b) *
                                  model_mode.get(i, j);
        }
      }
    }
  }
  auto bc_dt_v_minus =
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(get(data.lapse), 0.);
  gh::BoundaryConditions::Bjorhus::detail::add_physical_sector_projection(
      make_not_null(&bc_dt_v_minus), DataVector{-vars.char_speeds[3]},
      vars.projection_ab, vars.projection_Ab, vars.projection_AB, u3_of_mode);
  auto zero_psi =
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(get(data.lapse), 0.);
  auto zero_zero =
      make_with_value<tnsr::iaa<DataVector, Dim, frame>>(get(data.lapse), 0.);
  auto zero_plus =
      make_with_value<tnsr::aa<DataVector, Dim, frame>>(get(data.lapse), 0.);
  Corrections expected_difference = make_corrections(num_points);
  gh::BoundaryConditions::Bjorhus::project_corrections_onto_evolved_variables(
      make_not_null(&expected_difference.dt_spacetime_metric),
      make_not_null(&expected_difference.dt_pi),
      make_not_null(&expected_difference.dt_phi), make_not_null(&zero_psi),
      make_not_null(&zero_zero), make_not_null(&zero_plus),
      make_not_null(&bc_dt_v_minus), vars.char_speeds, data.gamma2,
      data.normal_covector);
  Corrections difference = with_model;
  add_to(make_not_null(&difference.dt_spacetime_metric),
         without_model.dt_spacetime_metric, -1.);
  add_to(make_not_null(&difference.dt_pi), without_model.dt_pi, -1.);
  add_to(make_not_null(&difference.dt_phi), without_model.dt_phi, -1.);
  check_corrections_equal(difference, expected_difference);
}

// The Kerr-Schild hole carries no tide, so the order-two model must find no
// tidal moments and reproduce the order-zero target; the invariant rapidity
// must boost the tangent member to the rest frame of the hole.
void test_quadrupole_model_on_boosted_kerr_schild() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 8;
  const std::array<double, 3> velocity{{0.2, -0.1, 0.15}};
  const auto data =
      kerr_schild_face_data(make_not_null(&generator), num_points);
  tnsr::ii<DataVector, Dim, frame> spatial_metric(num_points, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      spatial_metric.get(i, j) = data.spacetime_metric.get(i + 1, j + 1);
    }
  }
  const auto det_and_inverse = determinant_and_inverse(spatial_metric);
  const auto& inverse_spatial_metric = det_and_inverse.second;
  const auto extrinsic_curvature = gh::extrinsic_curvature(
      data.spacetime_unit_normal_vector, data.pi, data.phi);
  tnsr::ii<DataVector, Dim, frame> electric{};
  tnsr::ii<DataVector, Dim, frame> magnetic{};
  gh::worldtube::weyl_electric_magnetic(
      make_not_null(&electric), make_not_null(&magnetic), data.phi, data.d_phi,
      data.d_pi, data.spacetime_unit_normal_vector, spatial_metric,
      inverse_spatial_metric, extrinsic_curvature,
      data.inverse_spacetime_metric);

  // The Kretschmann scalar of the face curvature is 48 M^2 / r'^6 with r' the
  // rest-frame radius
  const auto analytic =
      TestHelpers::gh_worldtube::boosted_kretschmann(data.coords, velocity, 0.);
  CHECK_ITERABLE_CUSTOM_APPROX(get(gh::worldtube::kretschmann_scalar(
                                   electric, magnetic, inverse_spatial_metric)),
                               get(analytic.kretschmann),
                               Approx::custom().epsilon(1.e-6).scale(
                                   max(abs(get(analytic.kretschmann)))));

  gh::worldtube::KretschmannFaceData<Dim> face_data{};
  face_data.direction = Direction<Dim>::lower_zeta();
  face_data.time = 0.;
  face_data.kretschmann = analytic.kretschmann;
  face_data.d_kretschmann = analytic.d_kretschmann;
  face_data.dt_kretschmann = analytic.dt_kretschmann;
  face_data.quadrature_weights = DataVector(num_points, 1. / num_points);

  const auto quadrupole = gh::worldtube::evaluate_matching(
      Model::Quadrupole, 1.0, electric, magnetic, spatial_metric,
      data.normal_covector, data.lapse, data.shift, &face_data);
  const auto type_d = gh::worldtube::evaluate_matching(
      Model::TypeD, std::nullopt, electric, magnetic, spatial_metric,
      data.normal_covector, data.lapse, data.shift, nullptr);
  REQUIRE(quadrupole.rapidity.has_value());
  REQUIRE(quadrupole.registration.has_value());
  REQUIRE(quadrupole.second_order.has_value());

  // Rest frame: u = cosh(eta) Gamma (n + w) + sinh(eta) r_hat, with the
  // tangent member (Gamma, w, r_hat) in the adapted triad, moves with the
  // boost velocity of the hole
  const auto& member = quadrupole.registration->member;
  const auto& rotation = quadrupole.adapted_rotation;
  const auto cholesky_inverse =
      gr::np::inverse_lower_triangular(gr::np::cholesky_factor(spatial_metric));
  Approx velocity_approx = Approx::custom().epsilon(1.e-5).scale(1.);
  for (size_t p = 0; p < num_points; ++p) {
    const double eta = get(*quadrupole.rapidity)[p];
    const double gamma = get(member.lorentz_factor)[p];
    std::array<double, 3> u_triad{};
    for (size_t a = 0; a < 3; ++a) {
      gsl::at(u_triad, a) =
          std::cosh(eta) * gamma * member.transverse_velocity.get(a)[p] +
          std::sinh(eta) * member.radial_direction.get(a)[p];
    }
    const double u_0 = std::cosh(eta) * gamma / get(data.lapse)[p];
    for (size_t i = 0; i < 3; ++i) {
      double u_i = -data.shift.get(i)[p] * u_0;
      for (size_t a = 0; a < 3; ++a) {
        for (size_t k = 0; k < 3; ++k) {
          u_i += gsl::at(u_triad, a) * rotation.get(a, k)[p] *
                 cholesky_inverse.get(k, i)[p];
        }
      }
      CHECK(u_i / u_0 == velocity_approx(gsl::at(velocity, i)));
    }
  }

  // No tide: the fitted moments vanish and the target is the type-D one
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::abs(gsl::at(quadrupole.second_order->fit.components, a)) <
          1.e-5);
  }
  Approx fd_approx = Approx::custom().epsilon(1.e-7).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(get(quadrupole.psi0_target),
                               get(type_d.psi0_target), fd_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(quadrupole.incoming_mode, type_d.incoming_mode,
                               fd_approx);
  const double lorentz_factor =
      1. / sqrt(1. - square(0.2) - square(0.1) - square(0.15));
  for (size_t p = 0; p < num_points; ++p) {
    const double radius = get(quadrupole.registration->measured_radius)[p];
    CHECK(radius > 4. - 1.e-6);
    CHECK(radius < 4. * lorentz_factor + 1.e-6);
  }

  // Through the boundary condition the two models agree to the same accuracy
  const ExcisedSphere sphere{};
  const auto with_quadrupole =
      apply_worldtube(Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                Imposition::Frozen, Model::Quadrupole, 1.0},
                      data, 0., sphere.domain, sphere.element,
                      sphere.functions_of_time, face_data);
  const auto with_type_d = apply_worldtube(
      Worldtube{Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
                Model::TypeD},
      data, 0., sphere.domain, sphere.element, sphere.functions_of_time);
  CHECK_ITERABLE_CUSTOM_APPROX(with_quadrupole.dt_spacetime_metric,
                               with_type_d.dt_spacetime_metric, fd_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(with_quadrupole.dt_pi, with_type_d.dt_pi,
                               fd_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(with_quadrupole.dt_phi, with_type_d.dt_phi,
                               fd_approx);

  // With the fitted moments stored as the relaxed moments of the face data,
  // the condition skips the fit and imposes them: the same corrections
  {
    gh::worldtube::KretschmannFaceData<Dim> with_moments = face_data;
    with_moments.filtered_moments = quadrupole.second_order->fit.components;
    with_moments.filtered_moments_time = 0.;
    const auto imposed = gh::worldtube::evaluate_matching(
        Model::Quadrupole, 1.0, electric, magnetic, spatial_metric,
        data.normal_covector, data.lapse, data.shift, &with_moments,
        with_moments.filtered_moments);
    CHECK(imposed.second_order->fit.components ==
          quadrupole.second_order->fit.components);
    CHECK(imposed.second_order->fit.relative_residual ==
          approx(quadrupole.second_order->fit.relative_residual));
    CHECK_ITERABLE_APPROX(get(imposed.psi0_target),
                          get(quadrupole.psi0_target));
    const auto with_stored_moments = apply_worldtube(
        Worldtube{Imposition::Bjorhus, Imposition::Bjorhus, Imposition::Frozen,
                  Model::Quadrupole, 1.0, 10.0},
        data, 0., sphere.domain, sphere.element, sphere.functions_of_time,
        with_moments);
    check_corrections_equal(with_stored_moments, with_quadrupole);
    // Zero moments impose the type-D target
    with_moments.filtered_moments = gr::np::TidalMoments{};
    const auto zero_moments = gh::worldtube::evaluate_matching(
        Model::Quadrupole, 1.0, electric, magnetic, spatial_metric,
        data.normal_covector, data.lapse, data.shift, &with_moments,
        with_moments.filtered_moments);
    CHECK_ITERABLE_CUSTOM_APPROX(get(zero_moments.psi0_target),
                                 get(type_d.psi0_target), fd_approx);
  }

  // The order-two model needs the face data of the element
  CHECK_THROWS_WITH(
      apply_worldtube(Worldtube{Imposition::Bjorhus, Imposition::Bjorhus,
                                Imposition::Frozen, Model::Quadrupole, 1.0},
                      data, 0., sphere.domain, sphere.element,
                      sphere.functions_of_time),
      Catch::Matchers::ContainsSubstring("needs the Kretschmann face data"));
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
  test_incoming_mode_normalization();
  test_type_d_model_on_boosted_kerr_schild();
  test_quadrupole_model_on_boosted_kerr_schild();
}
