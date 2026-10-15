// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeEvolvedPhysicalGhost.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

#include <cmath>
#include <future>
#include <memory>
#include <numbers>
#include <string>
#include "DataStructures/Variables.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Shape.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/Creators/TimeDependentOptions/ExpansionMap.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Domain/Creators/TimeDependentOptions/ShapeMap.hpp"
#include "Domain/FunctionsOfTime/FixedInertialRadius.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplay.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/AveragedUpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/SphereReplay.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Serialization/Serialize.hpp"

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Domain/InterfaceLogicalCoordinates.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplayTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Characteristics.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeRecordingData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/LiftFlux.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"

namespace {
using domain::FunctionsOfTime::PiecewisePolynomial;
using gh::worldtube::ReplayTags;

double value(const size_t c, const double t, const double nx, const double ny,
             const double nz) {
  const double perturbation =
      static_cast<double>(c + 1) *
      (1. + .1 * t * nx + .02 * t * t * ny * nz + .001 * pow(t, 5) * nz);
  if (c < 10) {
    tnsr::aa<double, 3> eta(0.);
    get<0, 0>(eta) = -1.;
    get<1, 1>(eta) = 1.;
    get<2, 2>(eta) = 1.;
    get<3, 3>(eta) = 1.;
    return 3. * eta[c] + .001 * perturbation;
  }
  return perturbation;
}
tnsr::I<double, 3> linear_center(const double time) {
  tnsr::I<double, 3> center;
  for (size_t d = 0; d < 3; ++d) {
    center.get(d) = .2 * (d + 1) * time;
  }
  return center;
}
void write_tape(const std::string& filename, const bool theta_fast,
                const gh::worldtube::SphereReplay::CenterTrajectory&
                    center_at_time = linear_center,
                const bool physical = false) {
  const ylm::Spherepack sp(4, 4);
  const auto& angles = sp.theta_phi_points();
  const size_t nt = 5, np = 9, n = nt * np;
  h5::H5File<h5::AccessType::ReadWrite> file(filename);
  auto& volume = file.insert<h5::VolumeData>("GhWorldtube", 0);
  // Deliberately write observations out of chronological order.
  for (size_t ii = 0; ii < 8; ++ii) {
    const size_t i = 7 - ii;
    const double time = .1 * i;
    std::vector<ElementVolumeData> grids;
    for (const double radius : {.3, .4}) {
      std::vector<TensorComponent> components;
      Variables<ReplayTags> vars(n);
      for (size_t q = 0; q < n; ++q) {
        const double nx = sin(angles[0][q]) * cos(angles[1][q]);
        const double ny = sin(angles[0][q]) * sin(angles[1][q]);
        const double nz = cos(angles[0][q]);
        for (size_t c = 0; c < 50; ++c)
          vars.data()[c * n + q] = radius * value(c, time, nx, ny, nz);
      }
      const auto reorder = [&](const DataVector& data) {
        DataVector result(n);
        for (size_t j = 0; j < np; ++j)
          for (size_t k = 0; k < nt; ++k)
            result[theta_fast ? k + nt * j : j + np * k] = data[k + nt * j];
        return result;
      };
      const auto append = [&](const auto& tensor, const std::string& name) {
        for (size_t c = 0; c < tensor.size(); ++c)
          components.emplace_back(name + tensor.component_suffix(c),
                                  reorder(tensor[c]));
      };
      append(get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
             "SpacetimeMetric");
      append(get<gh::Tags::Pi<DataVector, 3>>(vars), "Pi");
      append(get<gh::Tags::Phi<DataVector, 3>>(vars), "Phi");
      if (physical) {
        tnsr::iaa<DataVector, 3> dpi(n);
        tnsr::ijaa<DataVector, 3> dphi(n);
        for (size_t c = 0; c < dpi.size(); ++c) {
          dpi[c] = DataVector(n, (c + 1.) * (1. + pow(time, 5)));
        }
        for (size_t c = 0; c < dphi.size(); ++c) {
          dphi[c] = DataVector(n, (c + 1.) * (2. - pow(time, 4)));
        }
        append(dpi, "deriv(Pi)");
        append(dphi, "deriv(Phi)");
        Scalar<DataVector> k(n, 2. + pow(time, 5));
        tnsr::i<DataVector, 3> dk(n);
        for (size_t d = 0; d < 3; ++d) {
          dk.get(d) = d + 1.;
          for (size_t q = 0; q < n; ++q) {
            const double direction =
                d == 2 ? cos(angles[0][q])
                       : sin(angles[0][q]) *
                             (d == 0 ? cos(angles[1][q]) : sin(angles[1][q]));
            get(k)[q] +=
                (d + 1.) * (center_at_time(time).get(d) + radius * direction);
          }
        }
        append(k, "Kretschmann");
        append(dk, "deriv(Kretschmann)");
      }
      components.emplace_back("InertialRadius", DataVector(n, radius));
      for (size_t d = 0; d < 3; ++d) {
        const double center = center_at_time(time).get(d);
        DataVector xyz(n);
        for (size_t q = 0; q < n; ++q)
          xyz[q] =
              center + radius * (d == 2 ? cos(angles[0][q])
                                        : sin(angles[0][q]) *
                                              (d == 0 ? cos(angles[1][q])
                                                      : sin(angles[1][q])));
        components.emplace_back(
            "InertialCoordinates_" + std::string(1, "xyz"[d]), reorder(xyz));
        components.emplace_back("InertialCenter_" + std::string(1, "xyz"[d]),
                                DataVector(n, center));
      }
      grids.emplace_back(
          radius == .3 ? "Sphere0" : "Sphere1", std::move(components),
          theta_fast ? std::vector<size_t>{nt, np}
                     : std::vector<size_t>{np, nt},
          std::vector{Spectral::Basis::SphericalHarmonic,
                      Spectral::Basis::SphericalHarmonic},
          theta_fast ? std::vector{Spectral::Quadrature::Gauss,
                                   Spectral::Quadrature::Equiangular}
                     : std::vector{Spectral::Quadrature::Equiangular,
                                   Spectral::Quadrature::Gauss});
    }
    volume.write_volume_data(100 - i, time, grids);
  }
}

void test_physical_reader(const bool theta_fast) {
  const std::string filename =
      theta_fast ? "PhysicalReplayTheta.h5" : "PhysicalReplayPhi.h5";
  if (file_system::check_if_file_exists(filename))
    file_system::rm(filename, true);
  write_tape(filename, theta_fast, linear_center, true);
  const gh::worldtube::SphereReplay tape(filename, "GhWorldtube", .4, true);
  const gh::worldtube::SphereReplay fields_only(filename, "GhWorldtube", .4);
  tnsr::I<double, 3> velocity;
  for (size_t d = 0; d < 3; ++d)
    velocity.get(d) = .2 * (d + 1.);
  for (const double t : {0., .3, .315, .7}) {
    tnsr::I<DataVector, 3> coords(3, 0.);
    const auto center = linear_center(t);
    for (size_t d = 0; d < 3; ++d) {
      coords.get(d) = center.get(d);
      coords.get(d)[d] += .4;
    }
    const auto data =
        tape.evaluate_physical(t, coords, linear_center, velocity);
    const auto expected_fields = fields_only.evaluate(t, coords, linear_center);
    tmpl::for_each<ReplayTags>([&](auto tag_v) {
      using Tag = tmpl::type_from<decltype(tag_v)>;
      CHECK_ITERABLE_APPROX(get<Tag>(data.fields), get<Tag>(expected_fields));
    });
    for (size_t c = 0; c < data.d_pi.size(); ++c) {
      CHECK_ITERABLE_APPROX(data.d_pi[c],
                            DataVector(3, (c + 1.) * (1. + pow(t, 5))));
    }
    for (size_t c = 0; c < data.d_phi.size(); ++c) {
      CHECK_ITERABLE_APPROX(data.d_phi[c],
                            DataVector(3, (c + 1.) * (2. - pow(t, 4))));
    }
    DataVector k(3, 2. + pow(t, 5));
    for (size_t d = 0; d < 3; ++d) {
      k += (d + 1.) * coords.get(d);
      CHECK_ITERABLE_APPROX(data.d_kretschmann.get(d), DataVector(3, d + 1.));
    }
    CHECK_ITERABLE_APPROX(get(data.kretschmann), k);
    // Differentiation divides SH transform roundoff by the tape cadence.
    CHECK_ITERABLE_CUSTOM_APPROX(get(data.dt_kretschmann),
                                 DataVector(3, 5. * pow(t, 4)),
                                 Approx::custom().epsilon(1.e-12).scale(1.));
    CHECK_THROWS_WITH(
        fields_only.evaluate_physical(t, coords, linear_center, velocity),
        Catch::Matchers::ContainsSubstring(
            "Physical sphere replay data was not loaded"));
  }
  file_system::rm(filename, true);
}

void write_curvature_tape(const std::string& filename,
                          const bool moving = false) {
  // Compute derivatives in a three-dimensional donor volume before slicing.
  const auto setup =
      TestHelpers::gh_worldtube::shell_element(12, 12, {{0., 0., 0.}}, .4, 1.0);
  const auto logical = interface_logical_coordinates(setup.mesh.slice_away(0),
                                                     Direction<3>::lower_xi());
  const auto coords = setup.element_map(logical);
  const gr::Solutions::KerrSchild solution(.11, {{.1, .15, -.2}},
                                           {{0., 0., 0.}});
  const auto evolved = TestHelpers::gh_worldtube::evolved_variables(
      solution, setup.inertial_coords, 0.);
  Variables<gh::worldtube::recording_tags> volume;
  gh::worldtube::compute_recording_data(
      make_not_null(&volume),
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
      get<gh::Tags::Pi<DataVector, 3>>(evolved),
      get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
      setup.inverse_jacobian);
  const auto face =
      intrp::Irregular<3>(setup.mesh, logical).interpolate(volume);
  const size_t n = get<0>(coords).size();
  std::vector<TensorComponent> components;
  tmpl::for_each<gh::worldtube::recording_tags>([&](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    const auto& tensor = get<Tag>(face);
    for (size_t c = 0; c < tensor.size(); ++c) {
      components.emplace_back(db::tag_name<Tag>() + tensor.component_suffix(c),
                              tensor[c]);
    }
  });
  components.emplace_back("InertialRadius", DataVector(n, .4));
  for (size_t d = 0; d < 3; ++d) {
    components.emplace_back("InertialCoordinates_" + std::string(1, "xyz"[d]),
                            coords.get(d));
    components.emplace_back("InertialCenter_" + std::string(1, "xyz"[d]),
                            DataVector(n, 0.));
  }
  h5::H5File<h5::AccessType::ReadWrite> file(filename);
  auto& tape = file.insert<h5::VolumeData>("GhWorldtube", 0);
  for (size_t i = 0; i < 8; ++i) {
    auto sample_components = components;
    if (moving) {
      for (size_t d = 0; d < 3; ++d) {
        const double center = linear_center(.01 * i).get(d);
        for (auto& c : sample_components) {
          if (c.name == "InertialCoordinates_" + std::string(1, "xyz"[d]) or
              c.name == "InertialCenter_" + std::string(1, "xyz"[d])) {
            std::get<DataVector>(c.data) += center;
          }
        }
      }
    }
    tape.write_volume_data(
        i, .1 * i,
        {ElementVolumeData{
            "Sphere",
            sample_components,
            {13, 25},
            {Spectral::Basis::SphericalHarmonic,
             Spectral::Basis::SphericalHarmonic},
            {Spectral::Quadrature::Gauss, Spectral::Quadrature::Equiangular}}});
  }
}

void test_recorded_target(const std::string& filename,
                          const bool moving = false) {
  using gh::BoundaryConditions::WorldtubeReplayTypeD;
  auto setup =
      TestHelpers::gh_worldtube::shell_element(8, 8, {{0., 0., 0.}}, .4, 1.4);
  auto coords = setup.element_map(interface_logical_coordinates(
      setup.mesh.slice_away(0), Direction<3>::lower_xi()));
  const size_t n = get<0>(coords).size();
  domain::FunctionsOfTimeMap fots;
  tnsr::I<double, 3> velocity(0.);
  if (moving) {
    fots["Translation"] = std::make_unique<PiecewisePolynomial<2>>(
        0.,
        std::array{DataVector{0., 0., 0.}, DataVector{.02, .04, .06},
                   DataVector{0., 0., 0.}},
        1.);
    setup.domain.inject_time_dependent_map_for_excision_sphere(
        setup.domain.excision_spheres().begin()->first,
        domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
            domain::CoordinateMaps::TimeDependent::Translation<3>{
                "Translation"}));
    for (size_t d = 0; d < 3; ++d) {
      coords.get(d) += linear_center(.0315).get(d);
      velocity.get(d) = .02 * (d + 1.);
    }
  }
  const auto center = [moving](const double t) {
    return moving ? linear_center(.1 * t) : tnsr::I<double, 3>(0.);
  };
  const auto data =
      gh::worldtube::SphereReplay(filename, "GhWorldtube", .4, true)
          .evaluate_physical(.315, coords, center, velocity);
  const auto curvature = gh::worldtube::weyl_curvature(
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(data.fields),
      get<gh::Tags::Pi<DataVector, 3>>(data.fields),
      get<gh::Tags::Phi<DataVector, 3>>(data.fields), data.d_pi, data.d_phi);
  tnsr::i<DataVector, 3> normal(n);
  for (size_t d = 0; d < 3; ++d)
    normal.get(d) = center(.315).get(d) - coords.get(d);
  const auto norm = magnitude(normal, curvature.inverse_spatial_metric);
  for (auto& x : normal)
    x /= get(norm);
  gh::worldtube::KretschmannFaceData<3> face;
  face.direction = Direction<3>::lower_xi();
  face.time = .315;
  face.kretschmann = data.kretschmann;
  face.d_kretschmann = data.d_kretschmann;
  // Exact advected-stationary derivative, independent of the interpolant's
  // cancellation and the BC's extraction of center velocity from the map.
  face.dt_kretschmann = Scalar<DataVector>(n, 0.);
  for (size_t d = 0; d < 3; ++d) {
    get(face.dt_kretschmann) -= velocity.get(d) * data.d_kretschmann.get(d);
  }
  face.quadrature_weights =
      gh::worldtube::face_quadrature_weights<3>(setup.mesh.slice_away(0));
  const auto expected = gh::worldtube::evaluate_matching(
      gh::worldtube::PhysicalModel::Quadrupole, .1, curvature.electric,
      curvature.magnetic, curvature.spatial_metric, normal, curvature.lapse,
      curvature.shift, &face, std::nullopt);
  const auto ghost = serialize_and_deserialize(
      TestHelpers::test_creation<
          gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost>(
          "RecordedWeyl:\n  H5Filename: " + filename +
          "\n  SubfileName: GhWorldtube\n  Radius: 0.4"));
  REQUIRE(ghost.uses_recorded_weyl());
  CHECK(ghost.get_clone() != nullptr);
  const WorldtubeReplayTypeD weyl_provider(
      filename, "GhWorldtube", .4, gh::worldtube::PhysicalModel::None,
      std::nullopt, std::nullopt, true, false, std::nullopt, true,
      "RecordedWeyl");
  tnsr::I<DataVector, 3> normal_vector(n, 0.);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      normal_vector.get(i) +=
          curvature.inverse_spatial_metric.get(i, j) * normal.get(j);
  const gh::worldtube::PhysicalScreen screen(
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(data.fields),
      gr::inverse_spacetime_metric(curvature.lapse, curvature.shift,
                                   curvature.inverse_spatial_metric),
      gr::spacetime_normal_vector(curvature.lapse, curvature.shift),
      normal_vector);
  CHECK_ITERABLE_APPROX(
      ghost.curvature_target(.315, coords, screen, &setup.domain,
                             &setup.element, &fots, &setup.mesh),
      gh::worldtube::incoming_weyl_to_physical_target(
          weyl_provider.recorded_incoming_mode(.315, coords, setup.domain,
                                               setup.element, fots, setup.mesh),
          screen));
  const WorldtubeReplayTypeD model(
      filename, "GhWorldtube", .4, gh::worldtube::PhysicalModel::Quadrupole, .1,
      std::nullopt, true, false, std::nullopt, true, "RecordedModel");
  CHECK_ITERABLE_APPROX(
      model.recorded_incoming_mode(.315, coords, setup.domain, setup.element,
                                   fots, setup.mesh),
      expected.incoming_mode);
  const WorldtubeReplayTypeD moments(
      filename, "GhWorldtube", .4, gh::worldtube::PhysicalModel::Quadrupole, .1,
      std::nullopt, true, false, std::nullopt, true, "RecordedMoments");
  CHECK_ITERABLE_APPROX(
      moments.recorded_tidal_moments(.315, coords, setup.domain, setup.element,
                                     fots, setup.mesh),
      expected.second_order->fit.components);
  Scalar<ComplexDataVector> psi0;
  get(psi0) = expected.psi.get(0);
  auto weyl = gr::np::orthonormal_to_coordinate_covariant(
      gr::np::incoming_weyl_field(psi0, expected.adapted_rotation),
      gr::np::cholesky_factor(curvature.spatial_metric));
  for (auto& c : weyl)
    c *= .5;
  const WorldtubeReplayTypeD control(
      filename, "GhWorldtube", .4, gh::worldtube::PhysicalModel::Quadrupole, .1,
      std::nullopt, true, false, std::nullopt, true, "RecordedWeyl");
  CHECK_ITERABLE_APPROX(
      control.recorded_incoming_mode(.315, coords, setup.domain, setup.element,
                                     fots, setup.mesh),
      weyl);
  const std::string options =
      "H5Filename: " + filename +
      "\nSubfileName: GhWorldtube\nRadius: 0.4\nPhysicalModel: Quadrupole"
      "\nMass: 0.1\nConstraintPreserving: true\nPhysicalTarget: ";
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          options + "RecordedModel\nMomentRelaxationTime: 1."),
      Catch::Matchers::ContainsSubstring("MomentRelaxationTime: None"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          options + "RecordedMoments\nMomentRelaxationTime: 1."),
      Catch::Matchers::ContainsSubstring("MomentRelaxationTime: None"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          "H5Filename: " + filename +
          "\nSubfileName: GhWorldtube\nRadius: 0.4\nPhysicalModel: TypeD"
          "\nMass: None\nMomentRelaxationTime: None"
          "\nPhysicalTarget: RecordedMoments"),
      Catch::Matchers::ContainsSubstring(
          "requires an order-two physical model"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          options + "RecordedMoments\nMomentRelaxationTime: None"
                    "\nReplayBoundaryData: false\nSommerfeldGauge: true"),
      Catch::Matchers::ContainsSubstring("require ReplayBoundaryData"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          options + "Unknown\nMomentRelaxationTime: None"),
      Catch::Matchers::ContainsSubstring("PhysicalTarget must be"));
}

void test_physical_replay(
    const std::string& filename, const bool quadrupole,
    const bool constraint_preserving = false,
    const bool sommerfeld_gauge = false, const bool qhat_radial_gauge = false,
    const bool replay_boundary_data = true,
    const std::string& physical_target = "Live", const bool geometric = false,
    const bool third_order = false, const bool pulse = false,
    const bool fixed_point = false,
    const std::optional<double> derivative_strength = std::nullopt) {
  using namespace gh::BoundaryConditions;
  const auto parsed = TestHelpers::test_creation<WorldtubeReplayTypeD>(
      "H5Filename: " + filename + "\nSubfileName: GhWorldtube\nRadius: 0.4" +
      (quadrupole ? std::string{"\nPhysicalModel: "} +
                        (fixed_point   ? "\n  QuadrupoleGeometricFixedPoint:"
                                         "\n    RelativeTolerance: 1.e-10"
                                         "\n    AbsoluteTolerance: 1.e-14"
                                         "\n    MaxIterations: 20"
                                         "\n    Damping: 1.0"
                         : third_order ? "ThirdOrderGeometric"
                                       : (geometric ? "QuadrupoleGeometric"
                                                    : "Quadrupole")) +
                        "\nMass: 0.1\nMomentRelaxationTime: None"
                  : "") +
      (constraint_preserving ? "\nConstraintPreserving: true" : "") +
      (sommerfeld_gauge ? "\nSommerfeldGauge: true" : "") +
      (qhat_radial_gauge ? "\nQhatRadialResponse: [-2.65, -2.33, -2.1]" : "") +
      (replay_boundary_data ? "" : "\nReplayBoundaryData: false") +
      "\nPhysicalTarget: " + physical_target +
      (pulse ? "\nMomentPulse:\n  StartTime: -1.\n  Duration: 2."
               "\n  Electric: [0.000001, 0., 0., 0., 0.]"
               "\n  Magnetic: [0., 0., 0.000002, 0., 0.]"
               "\nMomentDiagnosticsInterval: 0.01"
             : "") +
      (derivative_strength ? "\nDerivativeFitReplay:\n  LiveStrength: " +
                                 std::to_string(*derivative_strength) +
                                 "\n  RampStart: -2.\n  RampDuration: 1."
                           : ""));
  const auto hybrid = serialize_and_deserialize(parsed);
  const auto physical_model =
      third_order ? gh::worldtube::PhysicalModel::ThirdOrderGeometric
      : geometric ? gh::worldtube::PhysicalModel::QuadrupoleGeometric
                  : (quadrupole ? gh::worldtube::PhysicalModel::Quadrupole
                                : gh::worldtube::PhysicalModel::TypeD);
  const auto mass = quadrupole ? std::optional<double>{.1} : std::nullopt;
  CHECK(hybrid.physical_target() == physical_target);
  CHECK(hybrid.derivative_fit_replay().has_value() ==
        derivative_strength.has_value());
  if (derivative_strength) {
    CHECK((*hybrid.derivative_fit_replay())(0.) == *derivative_strength);
  }
  CHECK(hybrid.physical_model() == physical_model);
  CHECK(hybrid.psi0_fixed_point().has_value() == fixed_point);
  CHECK(hybrid.constraint_preserving() == constraint_preserving);
  CHECK(hybrid.sommerfeld_gauge() == sommerfeld_gauge);
  CHECK(hybrid.qhat_radial_response().has_value() == qhat_radial_gauge);
  CHECK(hybrid.replay_boundary_data() == replay_boundary_data);
  CHECK(hybrid.mass() == mass);
  CHECK_FALSE(hybrid.moment_relaxation_time().has_value());
  CHECK(dynamic_cast<WorldtubeReplayTypeD*>(hybrid.get_clone().get()) !=
        nullptr);
  const auto replay =
      replay_boundary_data
          ? std::make_optional<WorldtubeReplay>(filename, "GhWorldtube", .4)
          : std::nullopt;
  // A radially and angularly refined wedge: the lift must use the element's
  // Jacobian, not the block's or the inertial coordinate sphere radius.
  const auto setup =
      physical_target == "Live" and not(geometric or third_order)
          ? TestHelpers::gh_worldtube::wedge_element(
                6, {{0., 0., 0.}}, .4, 1.4, 0,
                {{SegmentId{1, 0}, SegmentId{1, 1}, SegmentId{2, 0}}})
          : TestHelpers::gh_worldtube::shell_element(8, 8, {{0., 0., 0.}}, .4,
                                                     1.4);
  const auto direction =
      setup.domain.excision_spheres().begin()->second.abutting_directions().at(
          setup.element.id().block_id());
  if (geometric or third_order) {
    gh::worldtube::validate_geometric_matching_face(
        setup.inertial_coords, setup.mesh, setup.inverse_jacobian, direction);
  }
  const auto logical = interface_logical_coordinates(
      setup.mesh.slice_away(direction.dimension()), direction);
  const auto coords = setup.element_map(logical);
  const auto invjac = setup.element_map.inv_jacobian(logical);
  const size_t n = get<0>(coords).size();
  const gr::Solutions::KerrSchild solution(.1, {{.2, -.1, .15}},
                                           {{0., 0., 0.}});
  const auto evolved = TestHelpers::gh_worldtube::evolved_variables(
      solution, setup.inertial_coords, 0.);
  Variables<gh::worldtube::recording_tags> volume;
  gh::worldtube::compute_recording_data(
      make_not_null(&volume),
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
      get<gh::Tags::Pi<DataVector, 3>>(evolved),
      get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
      setup.inverse_jacobian);
  const auto face =
      intrp::Irregular<3>(setup.mesh, logical).interpolate(volume);
  const auto& g = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(face);
  const auto& pi = get<gh::Tags::Pi<DataVector, 3>>(face);
  const auto& phi = get<gh::Tags::Phi<DataVector, 3>>(face);
  const auto& dg = get<
      gh::worldtube::recording_deriv<gr::Tags::SpacetimeMetric<DataVector, 3>>>(
      face);
  const auto& dpi =
      get<gh::worldtube::recording_deriv<gh::Tags::Pi<DataVector, 3>>>(face);
  const auto& dphi =
      get<gh::worldtube::recording_deriv<gh::Tags::Phi<DataVector, 3>>>(face);
  const auto inverse = determinant_and_inverse(gr::spatial_metric(g)).second;
  const auto shift = gr::shift(g, inverse);
  const auto lapse = gr::lapse(shift, g);
  const auto inverse4 = gr::inverse_spacetime_metric(lapse, shift, inverse);
  const auto time_normal = gr::spacetime_normal_vector(lapse, shift);
  tnsr::i<DataVector, 3> raw_normal(n);
  for (size_t i = 0; i < 3; ++i)
    raw_normal.get(i) = direction.sign() * invjac.get(direction.dimension(), i);
  const auto mag = magnitude(raw_normal, inverse);
  auto normal = raw_normal;
  tnsr::I<DataVector, 3> nv(n, 0.);
  for (size_t i = 0; i < 3; ++i)
    normal.get(i) /= get(mag);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      nv.get(i) += inverse.get(i, j) * normal.get(j);
  auto c3 = dg;
  for (size_t c = 0; c < c3.size(); ++c)
    c3[c] -= phi[c];
  Scalar<DataVector> gamma1(n, -.5);
  get(gamma1)[2] =
      -1.;  // A zero-speed metric field, including on moving meshes.
  const Scalar<DataVector> gamma2(n, 2.);
  const tnsr::a<DataVector, 3> gauge(n, .01);
  const tnsr::ab<DataVector, 3> dgauge(n, .02);
  const tnsr::aa<DataVector, 3> dtg(n, .03), dtp(n, -.02);
  const tnsr::iaa<DataVector, 3> dtf(n, .04);
  const std::unique_ptr<evolution::BoundaryCorrection> flux =
      std::make_unique<gh::BoundaryCorrections::AveragedUpwindPenalty<3>>();
  gh::worldtube::KretschmannFaceData<3> face_data;
  if (quadrupole) {
    gh::worldtube::update_kretschmann_face_data(
        make_not_null(&face_data),
        get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
        get<gh::Tags::Pi<DataVector, 3>>(evolved),
        get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
        setup.inverse_jacobian, setup.element, setup.domain.excision_spheres(),
        std::optional<tnsr::I<DataVector, 3>>{}, 0.);
    if (third_order) {
      gh::worldtube::update_third_order_matching(
          make_not_null(&face_data),
          get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
          get<gh::Tags::Pi<DataVector, 3>>(evolved),
          get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
          setup.inverse_jacobian, setup.inertial_coords, std::nullopt, *mass,
          0., TimeStepId{true, 0, Slab{0., 1.}.start()});
      REQUIRE(face_data.third_order.has_value());
      CHECK(face_data.third_order->derivative_order == 0);
      face_data = serialize_and_deserialize(face_data);
    } else if (fixed_point) {
      gh::worldtube::update_geometric_fixed_point(
          make_not_null(&face_data),
          get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
          get<gh::Tags::Pi<DataVector, 3>>(evolved),
          get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
          setup.inverse_jacobian, setup.inertial_coords, *mass, 0.,
          *hybrid.psi0_fixed_point());
      REQUIRE(face_data.psi0_fixed_point.has_value());
      CHECK_FALSE(face_data.filtered_moments.has_value());
      face_data = serialize_and_deserialize(face_data);
    } else {
      gh::worldtube::update_filtered_tidal_moments(
          make_not_null(&face_data),
          get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
          get<gh::Tags::Pi<DataVector, 3>>(evolved),
          get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.mesh,
          setup.inverse_jacobian, physical_model, *mass, std::nullopt, 0.);
      REQUIRE(face_data.filtered_moments.has_value());
    }
  }
  if (qhat_radial_gauge) {
    tnsr::I<DataVector, 3> radial = coords;
    for (auto& x : radial) {
      x /= .4;
    }
    const auto [q, drq] =
        gh::worldtube::qhat_radial_data(g, pi, phi, dg, dpi, dphi, radial);
    face_data.qhat_radial_gauge =
        gh::worldtube::QhatRadialGaugeData{q, drq, radial, {}, .4};
    // Nonzero response as well as native-lift cancellation is exercised.
    for (auto& x : face_data.qhat_radial_gauge->initial_q) {
      x *= .9;
    }
  }
  for (const bool moving : {false, true}) {
    CAPTURE(moving);
    std::optional<tnsr::I<DataVector, 3>> velocity;
    if (moving) {
      velocity = tnsr::I<DataVector, 3>(n);
      for (size_t i = 0; i < 3; ++i) {
        velocity->get(i) = .02 * static_cast<double>(i + 1);
        // Exercise outgoing-point masking on the same moving face.
        velocity->get(i)[0] = -5. * nv.get(i)[0];
        // Also exercise all incoming fields on this face.
        // For the model-only case keep u-plus outgoing while making the
        // advected fields incoming.
        velocity->get(i)[1] = replay_boundary_data
                                  ? 5. * nv.get(i)[1]
                                  : -shift.get(i)[1] + .2 * nv.get(i)[1];
        velocity->get(i)[3] = -shift.get(i)[3];  // Zero-speed u-zero.
      }
    }
    if (physical_target == "RecordedMoments") {
      // The donor moments must not freeze the live boost. Change its input
      // independently of the tape and leave stale live moments in face_data.
      get(face_data.dt_kretschmann) =
          (moving ? .02 : .01) * face_data.d_kretschmann.get(0);
      if (pulse) {
        const auto before = face_data;
        hybrid.observe_moment_response(
            0., get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
            get<gh::Tags::Pi<DataVector, 3>>(evolved),
            get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.inertial_coords,
            setup.mesh, setup.inverse_jacobian, setup.domain, setup.element, {},
            face_data);
        CHECK(face_data == before);
      }
    }
    tnsr::aa<DataVector, 3> extra_g(n), extra_p(n), model_g(n), model_p(n);
    tnsr::iaa<DataVector, 3> extra_f(n), model_f(n);
    CHECK_FALSE(hybrid
                    .dg_time_derivative(
                        make_not_null(&extra_g), make_not_null(&extra_p),
                        make_not_null(&extra_f), velocity, normal, nv, g, pi,
                        phi, coords, gamma1, gamma2, lapse, shift, inverse4,
                        time_normal, c3, gauge, dgauge, dtg, dtp, dtf, dg, dpi,
                        dphi, 0., setup.domain, setup.element, {}, setup.mesh,
                        ::dg::Formulation::StrongInertial, *flux,
                        setup.inverse_jacobian, face_data)
                    .has_value());
    const WorldtubeTypeD<3> model{
        constraint_preserving
            ? gh::BoundaryConditions::detail::SectorImposition::Bjorhus
            : gh::BoundaryConditions::detail::SectorImposition::Frozen,
        gh::BoundaryConditions::detail::SectorImposition::Bjorhus,
        sommerfeld_gauge
            ? gh::BoundaryConditions::detail::SectorImposition::
                  SommerfeldAbsorbing
            : gh::BoundaryConditions::detail::SectorImposition::Frozen,
        fixed_point
            ? WorldtubeReplayTypeD::PhysicalModel::
                  type{gh::BoundaryConditions::detail::GeometricFixedPoint{}}
            : WorldtubeReplayTypeD::PhysicalModel::type{physical_model},
        mass,
        hybrid.moment_relaxation_time()};
    std::optional<tnsr::ii<DataVector, 3>> prescribed{};
    if (physical_target == "RecordedMoments") {
      auto donor_moments = hybrid.recorded_tidal_moments(
          0., coords, setup.domain, setup.element, {}, setup.mesh);
      if (pulse) {
        // Independent expected input at the center of the pulse. Verify
        // injection through the moving-mesh Bjorhus/interface calculation.
        donor_moments[0] += 1.e-6;
        donor_moments[2] += std::complex<double>{0., 2.e-6};
      }
      const auto curvature =
          gh::worldtube::weyl_curvature(g, pi, phi, dpi, dphi);
      const auto expected = gh::worldtube::evaluate_matching(
          physical_model, mass, curvature.electric, curvature.magnetic,
          curvature.spatial_metric, normal, lapse, shift, &face_data,
          donor_moments);
      prescribed = expected.incoming_mode;
      // A full recorded target would be a different boundary condition.
      const auto full_recorded = hybrid.recorded_incoming_mode(
          0., coords, setup.domain, setup.element, {}, setup.mesh);
      double largest_difference = 0.;
      for (size_t c = 0; c < prescribed->size(); ++c) {
        largest_difference = std::max(
            largest_difference, max(abs((*prescribed)[c] - full_recorded[c])));
      }
      CHECK(largest_difference > 1.e-8);
    } else if (derivative_strength) {
      const auto live = gh::worldtube::evaluate_matching(
          physical_model, mass,
          gh::worldtube::weyl_curvature(g, pi, phi, dpi, dphi).electric,
          gh::worldtube::weyl_curvature(g, pi, phi, dpi, dphi).magnetic,
          gr::spatial_metric(g), normal, lapse, shift, &face_data);
      const auto reference_data =
          gh::worldtube::SphereReplay(filename, "GhWorldtube", .4, true)
              .evaluate_physical(
                  0., coords, [](double) { return tnsr::I<double, 3>(0.); },
                  tnsr::I<double, 3>(0.));
      auto mixed_dpi = dpi;
      auto mixed_dphi = dphi;
      for (size_t c = 0; c < dpi.size(); ++c) {
        mixed_dpi[c] = *derivative_strength * dpi[c] +
                       (1. - *derivative_strength) * reference_data.d_pi[c];
      }
      for (size_t c = 0; c < dphi.size(); ++c) {
        mixed_dphi[c] = *derivative_strength * dphi[c] +
                        (1. - *derivative_strength) * reference_data.d_phi[c];
      }
      const auto mixed =
          gh::worldtube::weyl_curvature(g, pi, phi, mixed_dpi, mixed_dphi);
      auto labels = coords;
      for (auto& x : labels) {
        x /= .4;
      }
      auto inward = labels;
      for (auto& x : inward) {
        x *= -1.;
      }
      const auto fit_psi = gr::np::weyl_scalars_from_electric_magnetic(
          mixed.electric, mixed.magnetic, mixed.spatial_metric, inward);
      const auto frame = gr::np::geometric_frame(
          *live.registration, *live.rapidity, live.adapted_rotation,
          mixed.spatial_metric, labels, face_data.quadrature_weights, *mass, 8);
      auto fit_registration = *live.registration;
      fit_registration.pulled_back =
          gr::np::pull_back(fit_psi, live.registration->rotation.a_bar,
                            live.registration->rotation.b);
      auto expected_moments = gr::np::evaluate_geometric_second_order(
                                  fit_registration, frame, *mass)
                                  .second_order.fit.components;
      const auto actual = hybrid.derivative_fit_evaluation(
          0., g, pi, phi, dpi, dphi, normal, coords, setup.domain,
          setup.element, {}, face_data);
      CHECK_ITERABLE_APPROX(actual.second_order->fit.components,
                            expected_moments);
      CHECK_ITERABLE_APPROX(actual.psi, live.psi);
      CHECK_ITERABLE_APPROX(actual.registration->pulled_back,
                            live.registration->pulled_back);
      CHECK_ITERABLE_APPROX(actual.registration->measured_radius,
                            live.registration->measured_radius);
      CHECK_ITERABLE_APPROX(*actual.rapidity, *live.rapidity);
      CHECK_ITERABLE_APPROX(actual.adapted_rotation, live.adapted_rotation);
      for (size_t c = 0; c < 5; ++c) {
        CHECK_ITERABLE_APPROX(actual.second_order->direct_columns[c],
                              live.second_order->direct_columns[c]);
      }
      if (*derivative_strength == 1.) {
        CHECK_ITERABLE_APPROX(actual.psi0_target, live.psi0_target);
        CHECK_ITERABLE_APPROX(actual.second_order->fit.components,
                              live.second_order->fit.components);
      }
      if (pulse) {
        expected_moments[0] += 1.e-6;
        expected_moments[2] += std::complex<double>{0., 2.e-6};
      }
      const auto expected = gh::worldtube::evaluate_matching(
          physical_model, mass,
          gh::worldtube::weyl_curvature(g, pi, phi, dpi, dphi).electric,
          gh::worldtube::weyl_curvature(g, pi, phi, dpi, dphi).magnetic,
          mixed.spatial_metric, normal, lapse, shift, &face_data,
          expected_moments);
      prescribed = expected.incoming_mode;
      if (not pulse) {
        CHECK_ITERABLE_APPROX(actual.psi0_target, expected.psi0_target);
      }
      const auto before = face_data;
      hybrid.observe_moment_response(
          0., get<gr::Tags::SpacetimeMetric<DataVector, 3>>(evolved),
          get<gh::Tags::Pi<DataVector, 3>>(evolved),
          get<gh::Tags::Phi<DataVector, 3>>(evolved), setup.inertial_coords,
          setup.mesh, setup.inverse_jacobian, setup.domain, setup.element, {},
          face_data);
      CHECK(face_data == before);
    } else if (physical_target != "Live") {
      prescribed = hybrid.recorded_incoming_mode(0., coords, setup.domain,
                                                 setup.element, {}, setup.mesh);
    }
    Bjorhus::IntermediateVariables<3> model_vars(n);
    CHECK_FALSE(model
                    .dg_time_derivative_with_intermediates(
                        make_not_null(&model_vars), make_not_null(&model_g),
                        make_not_null(&model_p), make_not_null(&model_f),
                        velocity, normal, nv, g, pi, phi, coords, gamma1,
                        gamma2, lapse, shift, inverse4, time_normal, c3, gauge,
                        dgauge, dtg, dtp, dtf, dg, dpi, dphi, 0., setup.domain,
                        setup.element, {}, face_data,
                        ConstraintDamping::Constant<3, Frame::Grid>{0.},
                        prescribed.has_value() ? &*prescribed : nullptr)
                    .has_value());
    Variables<ReplayTags> ghost(n), reference(n), hybrid_ghost(n);
    auto& gm = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(ghost);
    auto& gp = get<gh::Tags::Pi<DataVector, 3>>(ghost);
    auto& gf = get<gh::Tags::Phi<DataVector, 3>>(ghost);
    Scalar<DataVector> out1, out2, out_lapse;
    tnsr::I<DataVector, 3> out_shift;
    tnsr::II<DataVector, 3> out_inverse;
    if (replay_boundary_data) {
      CHECK_FALSE(
          replay
              ->dg_ghost(make_not_null(&gm), make_not_null(&gp),
                         make_not_null(&gf), make_not_null(&out1),
                         make_not_null(&out2), make_not_null(&out_lapse),
                         make_not_null(&out_shift), make_not_null(&out_inverse),
                         velocity, normal, nv, coords, gamma1, gamma2, 0.,
                         setup.domain, setup.element, {})
              .has_value());
    } else {
      gm = g;
      gp = pi;
      gf = phi;
    }
    CHECK_FALSE(
        hybrid
            .dg_ghost(
                make_not_null(&get<gr::Tags::SpacetimeMetric<DataVector, 3>>(
                    hybrid_ghost)),
                make_not_null(&get<gh::Tags::Pi<DataVector, 3>>(hybrid_ghost)),
                make_not_null(&get<gh::Tags::Phi<DataVector, 3>>(hybrid_ghost)),
                make_not_null(&out1), make_not_null(&out2),
                make_not_null(&out_lapse), make_not_null(&out_shift),
                make_not_null(&out_inverse), velocity, normal, nv, g, pi, phi,
                coords, gamma1, gamma2, lapse, shift, inverse4, time_normal, c3,
                gauge, dgauge, dtg, dtp, dtf, dg, dpi, dphi, 0., setup.domain,
                setup.element, {}, setup.mesh,
                ::dg::Formulation::StrongInertial, *flux,
                setup.inverse_jacobian, face_data)
            .has_value());
    if (not replay_boundary_data) {
      CHECK_ITERABLE_APPROX(out1, gamma1);
      CHECK_ITERABLE_APPROX(out2, gamma2);
      CHECK_ITERABLE_APPROX(out_lapse, lapse);
      CHECK_ITERABLE_APPROX(out_shift, shift);
      CHECK_ITERABLE_APPROX(out_inverse, inverse);
    }
    tmpl::for_each<ReplayTags>([&](auto tag_v) {
      using Tag = tmpl::type_from<decltype(tag_v)>;
      CHECK_ITERABLE_APPROX(get<Tag>(hybrid_ghost), get<Tag>(ghost));
    });
    auto donor_normal = normal;
    const auto donor_mag = magnitude(normal, out_inverse);
    for (auto& c : donor_normal)
      c /= -get(donor_mag);
    auto& rg = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(reference);
    auto& rp = get<gh::Tags::Pi<DataVector, 3>>(reference);
    auto& rf = get<gh::Tags::Phi<DataVector, 3>>(reference);
    const auto v = velocity.value_or(tnsr::I<DataVector, 3>(n, 0.));
    gh::BoundaryCorrections::AveragedUpwindPenalty<3>{}.dg_boundary_terms(
        make_not_null(&rg), make_not_null(&rp), make_not_null(&rf), g, pi, phi,
        gamma1, gamma2, normal, v, gm, gp, gf, gamma1, gamma2, donor_normal, v,
        ::dg::Formulation::StrongInertial);
    ::dg::lift_flux(make_not_null(&reference),
                    setup.mesh.extents(direction.dimension()), mag);
    if (not replay_boundary_data) {
      // Identical receiver/ghost states must give exactly zero DG penalty.
      for (size_t c = 0; c < reference.size(); ++c) {
        CHECK(reference.data()[c] == 0.);
      }
    }
    Bjorhus::IntermediateVariables<3> vars(n);
    Bjorhus::compute_intermediate_variables(
        make_not_null(&vars), velocity, normal, g, pi, phi, gamma1, gamma2,
        lapse, shift, inverse4, time_normal, c3, gauge, dgauge, dtg, dtp, dtf,
        dg, dpi, dphi);
    const auto extra = gh::characteristic_fields(gamma2, inverse, extra_g,
                                                 extra_p, extra_f, normal);
    const auto ref =
        gh::characteristic_fields(gamma2, inverse, rg, rp, rf, normal);
    const auto expected = gh::characteristic_fields(gamma2, inverse, model_g,
                                                    model_p, model_f, normal);
    // The large donor/live jump is amplified by the refined face lift;
    // subtracting its selected projections loses a few digits.
    const auto tolerance = Approx::custom().epsilon(1.e-11).scale(1.);
    const auto check_field = [&](const auto& e, const auto& r, const auto& m,
                                 const size_t speed, const bool replace) {
      for (size_t c = 0; c < e.size(); ++c) {
        for (size_t q = 0; q < n; ++q) {
          if (replace and vars.char_speeds[speed][q] < 0.) {
            CHECK(e[c][q] + r[c][q] == tolerance(m[c][q]));
          } else {
            CHECK(e[c][q] == tolerance(0.));
          }
        }
      }
    };
    check_field(get<gh::Tags::VSpacetimeMetric<DataVector, 3>>(extra),
                get<gh::Tags::VSpacetimeMetric<DataVector, 3>>(ref),
                get<gh::Tags::VSpacetimeMetric<DataVector, 3>>(expected), 0,
                constraint_preserving);
    check_field(get<gh::Tags::VZero<DataVector, 3>>(extra),
                get<gh::Tags::VZero<DataVector, 3>>(ref),
                get<gh::Tags::VZero<DataVector, 3>>(expected), 1,
                constraint_preserving);
    check_field(get<gh::Tags::VPlus<DataVector, 3>>(extra),
                get<gh::Tags::VPlus<DataVector, 3>>(ref),
                get<gh::Tags::VPlus<DataVector, 3>>(expected), 2, false);
    const auto physical = [&](const auto& x) {
      tnsr::aa<DataVector, 3> r(n, 0.);
      Bjorhus::detail::add_physical_sector_projection(
          make_not_null(&r), DataVector(n, 1.), vars.projection_ab,
          vars.projection_Ab, vars.projection_AB, x);
      return r;
    };
    const auto constraint = [&](const auto& x) {
      tnsr::aa<DataVector, 3> r(n, 0.);
      Bjorhus::detail::add_constraint_sector_projection(
          make_not_null(&r), DataVector(n, 1.), vars.outgoing_null_one_form,
          vars.incoming_null_vector, vars.projection_ab, vars.projection_Ab,
          vars.projection_AB, x);
      return r;
    };
    const auto gauge_sector = [&](const auto& x) {
      tnsr::aa<DataVector, 3> r(n, 0.);
      Bjorhus::detail::add_gauge_sector_projection(
          make_not_null(&r), DataVector(n, 1.), vars.incoming_null_one_form,
          vars.outgoing_null_one_form, vars.incoming_null_vector,
          vars.outgoing_null_vector, vars.projection_Ab, x);
      return r;
    };
    const auto& em = get<gh::Tags::VMinus<DataVector, 3>>(extra);
    const auto& rm = get<gh::Tags::VMinus<DataVector, 3>>(ref);
    const auto& mm = get<gh::Tags::VMinus<DataVector, 3>>(expected);
    check_field(physical(em), physical(rm), physical(mm), 3, true);
    check_field(constraint(em), constraint(rm), constraint(mm), 3,
                constraint_preserving);
    if (qhat_radial_gauge) {
      auto metric_correction = rg;
      if (constraint_preserving) {
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = a; b < 4; ++b) {
            for (size_t p = 0; p < n; ++p) {
              if (vars.char_speeds[0][p] < 0.) {
                metric_correction.get(a, b)[p] = model_g.get(a, b)[p];
              }
            }
          }
        }
      }
      const auto target = gh::worldtube::qhat_radial_gauge_correction(
          *face_data.qhat_radial_gauge, *hybrid.qhat_radial_response(), g, pi,
          phi, dg, dpi, dphi, coords, tnsr::I<double, 3>(0.), lapse, shift,
          velocity, gamma2, metric_correction, 0., {});
      const auto extra_gauge = gauge_sector(em);
      for (size_t a = 0; a < 4; ++a) {
        DataVector actual(n, 0.);
        for (size_t b = 0; b < 4; ++b) {
          actual +=
              vars.outgoing_null_vector.get(b) * (em.get(a, b) + rm.get(a, b));
        }
        for (size_t p = 0; p < n; ++p) {
          if (vars.char_speeds[3][p] < 0.) {
            CHECK(actual[p] == tolerance(target.get(a)[p]));
          } else {
            for (const auto& component : extra_gauge) {
              CHECK(component[p] == tolerance(0.));
            }
          }
        }
      }
    } else {
      check_field(gauge_sector(em), gauge_sector(rm), gauge_sector(mm), 3,
                  sommerfeld_gauge);
    }
    // The evolved-trace boundary retains the full original gauge/CP/advected
    // sectors and supplies only its physical incoming sector by a ghost flux.
    tnsr::aa<DataVector, 3> nonphys_g(n), nonphys_p(n);
    tnsr::iaa<DataVector, 3> nonphys_f(n);
    CHECK_FALSE(hybrid
                    .dg_time_derivative(
                        make_not_null(&nonphys_g), make_not_null(&nonphys_p),
                        make_not_null(&nonphys_f), velocity, normal, nv, g, pi,
                        phi, coords, gamma1, gamma2, lapse, shift, inverse4,
                        time_normal, c3, gauge, dgauge, dtg, dtp, dtf, dg, dpi,
                        dphi, 0., setup.domain, setup.element, {}, setup.mesh,
                        ::dg::Formulation::StrongInertial, *flux,
                        setup.inverse_jacobian, face_data, true)
                    .has_value());
    const auto nonphys = gh::characteristic_fields(
        gamma2, inverse, nonphys_g, nonphys_p, nonphys_f, normal);
    tmpl::for_each<tmpl::list<gh::Tags::VSpacetimeMetric<DataVector, 3>,
                              gh::Tags::VPlus<DataVector, 3>,
                              gh::Tags::VZero<DataVector, 3>>>([&](auto tag) {
      using Tag = tmpl::type_from<decltype(tag)>;
      const auto& actual = get<Tag>(nonphys);
      for (size_t c = 0; c < actual.size(); ++c) {
        CHECK_ITERABLE_CUSTOM_APPROX(
            actual[c], get<Tag>(extra)[c] + get<Tag>(ref)[c], tolerance);
      }
    });
    const auto& nonphys_minus = get<gh::Tags::VMinus<DataVector, 3>>(nonphys);
    CHECK_ITERABLE_CUSTOM_APPROX(physical(nonphys_minus),
                                 (tnsr::aa<DataVector, 3>(n, 0.)), tolerance);
    const auto check_nonphysical = [&](const auto& projector) {
      const auto actual = projector(nonphys_minus);
      const auto expected_extra = projector(em), expected_ref = projector(rm);
      for (size_t c = 0; c < actual.size(); ++c) {
        CHECK_ITERABLE_CUSTOM_APPROX(
            actual[c], expected_extra[c] + expected_ref[c], tolerance);
      }
    };
    check_nonphysical(constraint);
    check_nonphysical(gauge_sector);
    if (moving) {
      for (size_t c = 0; c < vars.char_speeds.size(); ++c) {
        CHECK(vars.char_speeds[c][0] > 0.);
        if (replay_boundary_data or c != 2) {
          CHECK(vars.char_speeds[c][1] < 0.);
        } else {
          CHECK(vars.char_speeds[c][1] > 0.);
        }
      }
      CHECK(vars.char_speeds[1][3] == 0.);
      if (not replay_boundary_data) {
        for (size_t i = 0; i < 3; ++i) {
          velocity->get(i)[1] = 5. * nv.get(i)[1];
        }
        const auto error = hybrid.dg_time_derivative(
            make_not_null(&extra_g), make_not_null(&extra_p),
            make_not_null(&extra_f), velocity, normal, nv, g, pi, phi, coords,
            gamma1, gamma2, lapse, shift, inverse4, time_normal, c3, gauge,
            dgauge, dtg, dtp, dtf, dg, dpi, dphi, 0., setup.domain,
            setup.element, {}, setup.mesh, ::dg::Formulation::StrongInertial,
            *flux, setup.inverse_jacobian, face_data);
        REQUIRE(error.has_value());
        CHECK_THAT(*error, Catch::Matchers::ContainsSubstring(
                               "no model for incoming u-plus"));
      }
    }
    CHECK(vars.char_speeds[0][2] == 0.);
    // Check that this exercised a real replacement, not a zero correction.
    CHECK(max(abs(get<1, 1>(extra_p))) > 1.e-8);
  }
}

void test_curved_center() {
  const std::string filename = "SphereReplayCurvedCenter.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  const auto trajectory = [](const double time) {
    auto center = linear_center(time);
    get<0>(center) += .01 * sin(20. * time);
    return center;
  };
  write_tape(filename, true, trajectory);
  const gh::worldtube::SphereReplay reader(filename, "GhWorldtube", .4);
  const double time = .015;
  const auto center = trajectory(time);
  // Demonstrate that interpolating these six sampled centers would trigger
  // the old 1e-7 center check, although the exact trajectory is correct.
  double interpolated_x = 0.;
  for (size_t i = 0; i < 6; ++i) {
    double weight = 1.;
    for (size_t j = 0; j < 6; ++j) {
      if (i != j) {
        weight *= (time - .1 * j) / (.1 * i - .1 * j);
      }
    }
    interpolated_x += weight * get<0>(trajectory(.1 * i));
  }
  CHECK(std::abs(interpolated_x - get<0>(center)) > 1.e-7);
  tnsr::I<DataVector, 3> xyz(1_st);
  for (size_t d = 0; d < 3; ++d) {
    xyz.get(d)[0] = center.get(d) + (d == 0 ? .4 : 0.);
  }
  const auto result = reader.evaluate(time, xyz, trajectory);
  for (size_t c = 0; c < 50; ++c) {
    CHECK(result.data()[c] == approx(.4 * value(c, time, 1., 0., 0.)));
  }
  // Repeated queries return independent copies of the cached fields.
  auto cached = reader.evaluate(time, xyz, trajectory);
  CHECK(cached == result);
  cached.data()[0] += 1.;
  CHECK(reader.evaluate(time, xyz, trajectory) == result);

  // Change time, angles and point count, including revisiting an older query.
  // These all need distinct values even though they use the same tape.
  const auto check_query = [&](const double query_time, const size_t count,
                               const bool y_axis) {
    const auto query_center = trajectory(query_time);
    tnsr::I<DataVector, 3> query_coords(count);
    for (size_t d = 0; d < 3; ++d) {
      query_coords.get(d) =
          query_center.get(d) + (d == (y_axis ? 1_st : 0_st) ? .4 : 0.);
    }
    return reader.evaluate(query_time, query_coords, trajectory);
  };
  for (const double query_time : {.02, time, .02}) {
    for (const size_t count : {1_st, 2_st}) {
      for (const bool y_axis : {false, true}) {
        const auto values = check_query(query_time, count, y_axis);
        for (size_t c = 0; c < 50; ++c) {
          for (size_t q = 0; q < count; ++q) {
            CHECK(values.data()[c * count + q] ==
                  approx(.4 * value(c, query_time, y_axis ? 0. : 1.,
                                    y_axis ? 1. : 0., 0.)));
          }
        }
      }
    }
  }
  // Boundary objects share the tape between worker threads. Exercise
  // simultaneous hits and misses with distinct geometry and time.
  std::vector<std::future<Variables<gh::worldtube::ReplayTags>>> queries;
  for (size_t i = 0; i < 8; ++i) {
    queries.push_back(std::async(std::launch::async, [&, i]() {
      return check_query(i % 2 == 0 ? time : .02, 1_st, i % 2 != 0);
    }));
  }
  for (size_t i = 0; i < queries.size(); ++i) {
    const auto values = queries[i].get();
    const bool y_axis = i % 2 != 0;
    for (size_t c = 0; c < 50; ++c) {
      CHECK(values.data()[c] ==
            approx(.4 * value(c, y_axis ? .02 : time, y_axis ? 0. : 1.,
                              y_axis ? 1. : 0., 0.)));
    }
  }
  // Refill the original key before testing trajectory/radius validation on
  // a would-be cache hit. A matching query must not bypass these checks.
  CHECK(reader.evaluate(time, xyz, trajectory) == result);
  // This trajectory matches the true center at the evaluation time, but not
  // at the tape samples. Keep detecting that genuine trajectory mismatch.
  const auto wrong_trajectory = [&](const double sample_time) {
    auto wrong = trajectory(sample_time);
    get<1>(wrong) += .001 * (sample_time - time);
    return wrong;
  };
  CHECK_THROWS_WITH(reader.evaluate(time, xyz, wrong_trajectory),
                    Catch::Matchers::ContainsSubstring(
                        "Receiver center does not follow the recorded sphere"));
  const auto wrong_current_center = [&](const double query_time) {
    auto current = trajectory(query_time);
    if (query_time == time) {
      get<0>(current) += .01;
    }
    return current;
  };
  CHECK_THROWS_WITH(reader.evaluate(time, xyz, wrong_current_center),
                    Catch::Matchers::ContainsSubstring("not on the recorded"));
  get<0>(xyz)[0] += .01;
  CHECK_THROWS_WITH(reader.evaluate(time, xyz, trajectory),
                    Catch::Matchers::ContainsSubstring("not on the recorded"));
  file_system::rm(filename, true);
}

void test_reader(const bool theta_fast) {
  const std::string filename =
      theta_fast ? "SphereReplayStrahl.h5" : "SphereReplayCce.h5";
  if (file_system::check_if_file_exists(filename))
    file_system::rm(filename, true);
  write_tape(filename, theta_fast);
  // Boundary-condition fixtures have a stationary domain, so their tape must
  // have a stationary center at all stencil times, not only at t=0.
  const std::string static_filename = "Static" + filename;
  if (file_system::check_if_file_exists(static_filename)) {
    file_system::rm(static_filename, true);
  }
  write_tape(static_filename, theta_fast,
             [](double) { return tnsr::I<double, 3>(0.); });
  test_physical_replay(static_filename, false);
  test_physical_replay(static_filename, true);
  test_physical_replay(static_filename, false, true);
  test_physical_replay(static_filename, true, true);
  test_physical_replay(static_filename, false, false, true);
  test_physical_replay(static_filename, false, true, true);
  test_physical_replay(static_filename, false, false, false, true);
  test_physical_replay(static_filename, true, true, false, true);
  // No recording exists at this path: both ghost packaging and the complete
  // model correction must work without constructing or reading a tape.
  const std::string missing_tape = "NonexistentModelBoundaryRecording.h5";
  REQUIRE_FALSE(file_system::check_if_file_exists(missing_tape));
  test_physical_replay(missing_tape, true, true, false, true, false);
  test_physical_replay(missing_tape, false, true, true, false, false);
  const auto no_tape_options =
      TestHelpers::test_creation<gh::BoundaryConditions::WorldtubeReplayTypeD>(
          "Radius: 0.4\nReplayBoundaryData: false\nConstraintPreserving: true"
          "\nQhatRadialResponse: [-2.65, -2.33, -2.1]");
  CHECK_FALSE(no_tape_options.replay_boundary_data());
  for (const std::string& incomplete :
       {std::string{""}, std::string{"\nConstraintPreserving: true"},
        std::string{"\nQhatRadialResponse: [-2.65, -2.33, -2.1]"}}) {
    CHECK_THROWS_WITH(
        TestHelpers::test_creation<
            gh::BoundaryConditions::WorldtubeReplayTypeD>(
            "Radius: 0.4\nReplayBoundaryData: false" + incomplete),
        Catch::Matchers::ContainsSubstring(
            "requires ConstraintPreserving: true and either"));
  }
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<gh::BoundaryConditions::WorldtubeReplayTypeD>(
          "H5Filename: " + filename +
          "\nSubfileName: GhWorldtube\nRadius: 0.4\nSommerfeldGauge: true"
          "\nQhatRadialResponse: [-2.65, -2.33, -2.1]"),
      Catch::Matchers::ContainsSubstring("mutually exclusive"));
  const gh::worldtube::SphereReplay reader(filename, "GhWorldtube", .4);
  CHECK(reader.time_bounds()[0] == 0.);
  CHECK(reader.time_bounds()[1] == approx(.7));
  for (const double time : {0., .015, .3, .347, .699, .7}) {
    tnsr::I<double, 3> center;
    tnsr::I<DataVector, 3> xyz(3);
    const std::array<std::array<double, 3>, 3> directions{
        {{{1., 0., 0.}}, {{0., 0., 1.}}, {{.36, .48, .8}}}};
    for (size_t d = 0; d < 3; ++d) {
      center.get(d) = .2 * (d + 1) * time;
      for (size_t q = 0; q < 3; ++q)
        xyz.get(d)[q] = center.get(d) + .4 * directions[q][d];
    }
    const auto result = reader.evaluate(time, xyz, linear_center);
    for (size_t c = 0; c < 50; ++c)
      for (size_t q = 0; q < 3; ++q)
        CHECK(result.data()[c * 3 + q] ==
              approx(.4 * value(c, time, directions[q][0], directions[q][1],
                                directions[q][2])));
    CHECK_THROWS_WITH(reader.evaluate(-.01, xyz, linear_center),
                      Catch::Matchers::ContainsSubstring("extrapolate"));
    CHECK_THROWS_WITH(reader.evaluate(.71, xyz, linear_center),
                      Catch::Matchers::ContainsSubstring("extrapolate"));
    xyz.get(0)[0] += .01;
    CHECK_THROWS_WITH(
        reader.evaluate(time, xyz, linear_center),
        Catch::Matchers::ContainsSubstring("not on the recorded"));
  }
  const auto first =
      gh::worldtube::cached_sphere_replay(filename, "GhWorldtube", .4);
  CHECK(first ==
        gh::worldtube::cached_sphere_replay(filename, "GhWorldtube", .4));
  const auto bc =
      TestHelpers::test_creation<gh::BoundaryConditions::WorldtubeReplay>(
          "H5Filename: " + static_filename +
          "\nSubfileName: GhWorldtube\nRadius: 0.4");
  const auto migrated = serialize_and_deserialize(bc);
  CHECK(migrated.get_clone() != nullptr);
  // Exercise the native ghost contract and moving-mesh numerical flux.
  const auto setup =
      TestHelpers::gh_worldtube::wedge_element(4, {{0., 0., 0.}}, .4, 2.);
  tnsr::I<DataVector, 3> xyz(1_st, 0.), velocity(1_st, .2), shift(1_st, 0.);
  get<0>(xyz) = .4;
  tnsr::i<DataVector, 3> normal(1_st, 0.);
  get<0>(normal) = -1.;
  tnsr::I<DataVector, 3> normal_vector(1_st, 0.);
  tnsr::aa<DataVector, 3> metric(1_st, 0.), pi(1_st, 0.);
  tnsr::iaa<DataVector, 3> phi(1_st, 0.);
  tnsr::II<DataVector, 3> inverse(1_st, 0.);
  Scalar<DataVector> gamma1(1_st, -.5), gamma2(1_st, 2.), lapse(1_st, 0.), out1,
      out2;
  CHECK_FALSE(migrated
                  .dg_ghost(make_not_null(&metric), make_not_null(&pi),
                            make_not_null(&phi), make_not_null(&out1),
                            make_not_null(&out2), make_not_null(&lapse),
                            make_not_null(&shift), make_not_null(&inverse),
                            std::optional{velocity}, normal, normal_vector, xyz,
                            gamma1, gamma2, 0., setup.domain, setup.element, {})
                  .has_value());
  CHECK_ITERABLE_APPROX(out1, gamma1);
  CHECK_ITERABLE_APPROX(out2, gamma2);
  CHECK(get(lapse)[0] > 0.);
  const auto reference = reader.evaluate(0., xyz, linear_center);
  CHECK_ITERABLE_APPROX(
      metric, (get<gr::Tags::SpacetimeMetric<DataVector, 3>>(reference)));
  CHECK_ITERABLE_APPROX(pi, (get<gh::Tags::Pi<DataVector, 3>>(reference)));
  CHECK_ITERABLE_APPROX(phi, (get<gh::Tags::Phi<DataVector, 3>>(reference)));
  get<0>(normal) = -1. / sqrt(get<0, 0>(inverse));
  auto opposite = normal;
  for (auto& c : opposite)
    c = -c;
  tnsr::aa<DataVector, 3> dg(1_st, 0.), dp(1_st, 0.);
  tnsr::iaa<DataVector, 3> df(1_st, 0.);
  gh::BoundaryCorrections::AveragedUpwindPenalty<3>{}.dg_boundary_terms(
      make_not_null(&dg), make_not_null(&dp), make_not_null(&df), metric, pi,
      phi, gamma1, gamma2, normal, velocity, metric, pi, phi, gamma1, gamma2,
      opposite, velocity, dg::Formulation::StrongInertial);
  for (const auto& c : dg)
    CHECK(max(abs(c)) == 0.);
  for (const auto& c : dp)
    CHECK(max(abs(c)) == 0.);
  for (const auto& c : df)
    CHECK(max(abs(c)) == 0.);

  CHECK_THROWS_WITH((gh::worldtube::SphereReplay(filename, "GhWorldtube", .35)),
                    Catch::Matchers::ContainsSubstring("absent"));
  file_system::rm(filename, true);
  file_system::rm(static_filename, true);
}
void test_map_replay_options() {
  domain::FunctionsOfTime::register_derived_with_charm();
  const std::string filename = "SphereReplayMaps.h5";
  if (file_system::check_if_file_exists(filename))
    file_system::rm(filename, true);
  domain::FunctionsOfTimeMap global, local;
  global["Expansion"] = std::make_unique<PiecewisePolynomial<2>>(
      0., std::array{DataVector{.8}, DataVector{.03}, DataVector{0.}}, 2.);
  local["Expansion"] = global.at("Expansion")->create_at_time(0., .1);
  global["ExpansionOuterBoundary"] =
      std::make_unique<domain::FunctionsOfTime::FixedSpeedCubic>(1., 0., .01,
                                                                 10.);
  const DataVector zero(ylm::Spherepack::spectral_size(4, 4), 0.);
  global["ShapeB"] = std::make_unique<PiecewisePolynomial<2>>(
      0., std::array{zero, zero, zero}, 2.);
  global["SizeB"] = std::make_unique<PiecewisePolynomial<3>>(
      0.,
      std::array{DataVector{0.}, DataVector{0.}, DataVector{0.},
                 DataVector{0.}},
      2.);
  {
    h5::H5File<h5::AccessType::ReadWrite> file(filename);
    auto& vol = file.insert<h5::VolumeData>("ReplayVolume", 0);
    vol.write_volume_data(
        17, 0.,
        {ElementVolumeData{"dummy",
                           {TensorComponent{"x", DataVector{0., 1.}}},
                           {2},
                           {Spectral::Basis::Legendre},
                           {Spectral::Quadrature::GaussLobatto}}},
        std::nullopt, serialize(local), serialize(global));
  }
  using namespace domain::creators::time_dependent_options;
  const auto replay = TestHelpers::test_creation<FromVolumeFile>(
      "H5Filename: " + filename + "\nSubfileName: ReplayVolume\nReplay: true");
  const auto ordinary = TestHelpers::test_creation<FromVolumeFile>(
      "H5Filename: " + filename + "\nSubfileName: ReplayVolume");
  CHECK_FALSE(ordinary.replay());
  CHECK(ordinary.retrieve_function_of_time({"Expansion"}, 0.)
            .at("Expansion")
            ->time_bounds()[1] == .1);
  // No observation exists at 1.5, and the local history expires at .1.
  auto read = replay.retrieve_function_of_time({"Expansion"}, 1.5);
  CHECK(read.at("Expansion")->time_bounds()[1] == 2.);
  CHECK_ITERABLE_APPROX(read.at("Expansion")->func_and_2_derivs(1.5),
                        global.at("Expansion")->func_and_2_derivs(1.5));
  const auto expansions = get_expansion<false>(
      std::variant<ExpansionMapOptions<false>, FromVolumeFile>{replay}, 1.5,
      1.6);
  REQUIRE(expansions.size() == 2);
  CHECK_ITERABLE_APPROX(expansions.at("Expansion")->func(1.5),
                        global.at("Expansion")->func(1.5));
  CHECK_ITERABLE_APPROX(expansions.at("ExpansionOuterBoundary")->func(1.5),
                        global.at("ExpansionOuterBoundary")->func(1.5));
  const auto shape = TestHelpers::test_creation<
      FromVolumeFileShapeSize<domain::ObjectLabel::B>>(
      "LMax: Auto\nCoefficientTruncationLimit: 0.\nTransitionEndsAtCube: "
      "true\nH5Filename: " +
      filename +
      "\nSubfileName: ReplayVolume\nReplay: true\nFixedInertialRadius: 0.3");
  const auto result = get_shape_and_size<true, domain::ObjectLabel::B>(
      std::variant<ShapeMapOptions<true, domain::ObjectLabel::B>,
                   FromVolumeFileShapeSize<domain::ObjectLabel::B>>{shape},
      0., .1, .1, .375);
  CHECK(dynamic_cast<const domain::FunctionsOfTime::FixedInertialRadius*>(
            result.at("SizeB").get()) != nullptr);
  CHECK(max(abs(result.at("ShapeB")->func(1.5)[0])) == 0.);
  CHECK(result.at("SizeB")->time_bounds()[1] == 2.);
  file_system::rm(filename, true);
}
void test_fixed_radius() {
  domain::FunctionsOfTime::register_derived_with_charm();
  domain::FunctionsOfTimeMap fots;
  auto expansion = std::make_unique<PiecewisePolynomial<2>>(
      0., std::array{DataVector{.8}, DataVector{.03}, DataVector{.004}}, 1.);
  expansion->update(1., DataVector{-.002}, 2.);
  fots["Size"] = std::make_unique<domain::FunctionsOfTime::FixedInertialRadius>(
      .375, .3, expansion->get_clone());
  const DataVector zero(ylm::Spherepack::spectral_size(4, 4), 0.);
  fots["Shape"] = std::make_unique<PiecewisePolynomial<2>>(
      0., std::array{zero, zero, zero}, 2.);
  const domain::CoordinateMaps::TimeDependent::Shape shape(
      {{0., 0., 0.}}, 0.,
      std::make_unique<domain::CoordinateMaps::ShapeMapTransitionFunctions::
                           SphereTransition>(.375, 2.),
      "Shape", "Size");
  const auto copy = serialize_and_deserialize(fots["Size"]);
  for (const double time : {.1, .7, 1., 1.1, 1.9}) {
    const auto a = expansion->func_and_2_derivs(time);
    const auto size = fots.at("Size")->func_and_2_derivs(time);
    CHECK_ITERABLE_APPROX(copy->func_and_2_derivs(time), size);
    const double norm = sqrt(4. * std::numbers::pi);
    const double r = .375 - size[0][0] / norm;
    CHECK(a[0][0] * r == approx(.3));
    CHECK(a[1][0] * r - a[0][0] * size[1][0] / norm == approx(0.));
    CHECK(a[2][0] * r - 2. * a[1][0] * size[1][0] / norm -
              a[0][0] * size[2][0] / norm ==
          approx(0.));
    const std::array<double, 3> x{{.375, 0., 0.}};
    const auto mapped = shape(x, time, fots);
    const auto velocity = shape.frame_velocity(x, time, fots);
    CHECK(a[0][0] * mapped[0] == approx(.3));
    CHECK(a[1][0] * mapped[0] + a[0][0] * velocity[0] == approx(0.));
    const auto jac = shape.jacobian(x, time, fots);
    const auto inv = shape.inv_jacobian(x, time, fots);
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j) {
        double sum = 0.;
        for (size_t k = 0; k < 3; ++k)
          sum += jac.get(i, k) * inv.get(k, j);
        CHECK(sum == approx(i == j ? 1. : 0.));
      }
  }
}
}  // namespace
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.SphereReplay", "[Unit][Evolution]") {
  test_reader(true);
  test_reader(false);
  test_curved_center();
  test_fixed_radius();
  test_map_replay_options();
}

// [[TimeOut, 15]]
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.RecordedPhysicalTarget",
                  "[Unit][Evolution]") {
  test_physical_reader(true);
  test_physical_reader(false);
  const std::string filename = "RecordedCurvatureReplay.h5";
  if (file_system::check_if_file_exists(filename))
    file_system::rm(filename, true);
  write_curvature_tape(filename);
  test_recorded_target(filename);
  test_physical_replay(filename, true, true, false, false, true,
                       "RecordedModel");
  test_physical_replay(filename, true, true, false, false, true,
                       "RecordedWeyl");
  file_system::rm(filename, true);
  const std::string moving_filename = "RecordedMovingCurvatureReplay.h5";
  if (file_system::check_if_file_exists(moving_filename))
    file_system::rm(moving_filename, true);
  write_curvature_tape(moving_filename, true);
  test_recorded_target(moving_filename, true);
  file_system::rm(moving_filename, true);
}

// [[TimeOut, 15]]
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.RecordedMoments", "[Unit][Evolution]") {
  const std::string filename = "RecordedMomentsReplay.h5";
  if (file_system::check_if_file_exists(filename))
    file_system::rm(filename, true);
  write_curvature_tape(filename);
  test_recorded_target(filename);
  test_physical_replay(filename, true, true, false, false, true,
                       "RecordedMoments");
  test_physical_replay(filename, true, true, false, false, true,
                       "RecordedMoments", false, false, true);
  file_system::rm(filename, true);
}

SPECTRE_TEST_CASE("Unit.Gh.Worldtube.MomentPulse", "[Unit][Evolution]") {
  using gh::BoundaryConditions::detail::RecordedMomentPulse;
  const auto pulse =
      serialize_and_deserialize(TestHelpers::test_creation<RecordedMomentPulse>(
          "StartTime: 20.\nDuration: 2.\nElectric: [1., 2., 3., 4., 5.]"
          "\nMagnetic: [-1., 0., 1., 2., 3.]"));
  const gr::np::TidalMoments zero{};
  CHECK(pulse(19.) == zero);
  CHECK(pulse(20.) == zero);
  CHECK(pulse(22.) == zero);
  CHECK(pulse(23.) == zero);
  for (size_t a = 0; a < 5; ++a) {
    const std::complex<double> peak{1. + a, -1. + a};
    CHECK(pulse(21.)[a].real() == approx(peak.real()));
    CHECK(pulse(21.)[a].imag() == approx(peak.imag()));
    CHECK(pulse(20.5)[a].real() == approx(.25 * peak.real()));
    CHECK(pulse(20.5)[a].imag() == approx(.25 * peak.imag()));
    CHECK(pulse(21.5)[a].real() == approx(.25 * peak.real()));
    CHECK(pulse(21.5)[a].imag() == approx(.25 * peak.imag()));
  }
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<RecordedMomentPulse>(
          "StartTime: 20.\nDuration: 0.\nElectric: [0.,0.,0.,0.,0.]"
          "\nMagnetic: [0.,0.,0.,0.,0.]"),
      Catch::Matchers::ContainsSubstring("positive finite Duration"));
  const std::string filename = "MomentPulseReplay.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  write_curvature_tape(filename);
  const std::string options =
      "H5Filename: " + filename +
      "\nSubfileName: GhWorldtube\nRadius: 0.4"
      "\nPhysicalModel: Quadrupole\nMass: 0.1\nMomentRelaxationTime: None";
  using gh::BoundaryConditions::WorldtubeReplayTypeD;
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<WorldtubeReplayTypeD>(
          options + "\nMomentDiagnosticsInterval: 0.01"),
      Catch::Matchers::ContainsSubstring("PhysicalTarget: RecordedMoments"));
  CHECK_THROWS_WITH(TestHelpers::test_creation<WorldtubeReplayTypeD>(
                        options + "\nPhysicalTarget: RecordedMoments"
                                  "\nMomentDiagnosticsInterval: -0.1"),
                    Catch::Matchers::ContainsSubstring("positive and finite"));
  file_system::rm(filename, true);
}

// [[TimeOut, 30]]
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.GeometricReplay", "[Unit][Evolution]") {
  const std::string filename = "GeometricCurvatureReplay.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  write_curvature_tape(filename);
  // Exercise the live geometric physical sector with CP, recorded gauge,
  // and both moving and stationary faces. The replay penalty must cancel
  // in the physical/CP sectors while retaining the recorded gauge sector.
  test_physical_replay(filename, true, true, false, false, true, "Live", true);
  test_physical_replay(filename, true, true, false, false, true, "Live", true,
                       false, false, true);
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<gh::BoundaryConditions::WorldtubeReplayTypeD>(
          "H5Filename: " + filename +
          "\nSubfileName: GhWorldtube\nRadius: 0.4\nMass: 0.1"
          "\nMomentRelaxationTime: None\nPhysicalTarget: RecordedMoments"
          "\nPhysicalModel:\n  QuadrupoleGeometricFixedPoint:"
          "\n    RelativeTolerance: 1.e-10\n    AbsoluteTolerance: 1.e-14"
          "\n    MaxIterations: 20\n    Damping: 1.0"),
      Catch::Matchers::ContainsSubstring("PhysicalTarget: Live"));
  for (const bool model_time : {false, true}) {
    const std::string options =
        "H5Filename: " + filename +
        "\nSubfileName: GhWorldtube\nRadius: 0.4"
        "\nPhysicalModel: QuadrupoleGeometric\nMass: 0.1"
        "\nConstraintPreserving: true\nMomentRelaxationTime: " +
        (model_time ? "{ModelTime: {Timescale: 1.0}}" : "1.0");
    const auto relaxed = serialize_and_deserialize(
        TestHelpers::test_creation<
            gh::BoundaryConditions::WorldtubeReplayTypeD>(options));
    REQUIRE(relaxed.moment_relaxation_time().has_value());
    CHECK(*relaxed.moment_relaxation_time() == 1.0);
    CHECK(relaxed.moment_relaxation_uses_model_time() == model_time);
    CHECK_THROWS_WITH(
        TestHelpers::test_creation<
            gh::BoundaryConditions::WorldtubeReplayTypeD>(
            options + "\nPhysicalTarget: RecordedMoments"),
        Catch::Matchers::ContainsSubstring("MomentRelaxationTime: None"));
  }
  file_system::rm(filename, true);
}

// [[TimeOut, 30]]
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.ThirdOrderReplay", "[Unit][Evolution]") {
  const std::string filename = "ThirdOrderCurvatureReplay.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  write_curvature_tape(filename);
  // The serialized causal target must reach the physical Bjorhus sector,
  // with live CP and recorded gauge, including moving-face speed masking.
  test_physical_replay(filename, true, true, false, false, true, "Live", false,
                       true);
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<gh::BoundaryConditions::WorldtubeReplayTypeD>(
          "H5Filename: " + filename +
          "\nSubfileName: GhWorldtube\nRadius: 0.4"
          "\nPhysicalModel: ThirdOrderGeometric\nMass: 0.1"
          "\nMomentRelaxationTime: None\nPhysicalTarget: RecordedModel"),
      Catch::Matchers::ContainsSubstring("needs live causal history"));
  file_system::rm(filename, true);
}

// [[TimeOut, 180]]
SPECTRE_TEST_CASE("Unit.Gh.Worldtube.DerivativeFitReplay",
                  "[Unit][Evolution]") {
  using gh::BoundaryConditions::detail::DerivativeFitReplay;
  for (double strength : {0., .5, 1.}) {
    const auto ramp =
        serialize_and_deserialize(DerivativeFitReplay{strength, 2., 4.});
    CHECK(ramp(1.) == 1.);
    CHECK(ramp(2.) == 1.);
    CHECK(ramp(4.) == approx(.5 * (1. + strength)));
    CHECK(ramp(6.) == strength);
    CHECK(ramp(7.) == strength);
  }
  CHECK_THROWS_WITH((DerivativeFitReplay{-1., 0., 1.}),
                    Catch::Matchers::ContainsSubstring("LiveStrength"));
  CHECK_THROWS_WITH((DerivativeFitReplay{0., 0., 0.}),
                    Catch::Matchers::ContainsSubstring("RampDuration"));
  const std::string filename = "DerivativeFitReplay.h5";
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
  write_curvature_tape(filename);
  for (double strength : {0., .5, 1.}) {
    test_physical_replay(filename, true, true, false, false, true, "Live", true,
                         false, false, false, strength);
  }
  test_physical_replay(filename, true, true, false, false, true, "Live", true,
                       false, true, false, 0.);
  file_system::rm(filename, true);
}
