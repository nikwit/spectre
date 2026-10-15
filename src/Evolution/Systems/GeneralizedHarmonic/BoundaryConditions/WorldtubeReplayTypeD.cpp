// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplayTypeD.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <numbers>
#include <sstream>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Slice.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/AveragedUpwindPenalty.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/LiftFlux.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/ParseError.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/LongitudinalTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"

namespace gh::BoundaryConditions {
detail::DerivativeFitReplay::DerivativeFitReplay(
    const double strength, const double start, const double duration,
    const Options::Context& context)
    : live_strength(strength), ramp_start(start), ramp_duration(duration) {
  if (not std::isfinite(strength) or strength < 0. or strength > 1. or
      not std::isfinite(start) or not std::isfinite(duration) or
      duration <= 0. or not std::isfinite(start + duration)) {
    PARSE_ERROR(context,
                "DerivativeFitReplay requires finite LiveStrength in "
                "[0,1], finite RampStart and positive RampDuration.");
  }
}
double detail::DerivativeFitReplay::operator()(const double time) const {
  if (time <= ramp_start or live_strength == 1.) {
    return 1.;
  }
  if (time >= ramp_start + ramp_duration) {
    return live_strength;
  }
  const double x = (time - ramp_start) / ramp_duration;
  const double step = x * x * x * (10. + x * (-15. + 6. * x));
  return 1. + (live_strength - 1.) * step;
}
void detail::DerivativeFitReplay::pup(PUP::er& p) {
  p | live_strength;
  p | ramp_start;
  p | ramp_duration;
}

detail::RecordedMomentPulse::RecordedMomentPulse(
    const double start_in, const double duration_in,
    std::array<double, 5> electric_in, std::array<double, 5> magnetic_in,
    const Options::Context& context)
    : start(start_in),
      duration(duration_in),
      electric(electric_in),
      magnetic(magnetic_in) {
  if (not std::isfinite(start) or not std::isfinite(duration) or
      duration <= 0. or not std::isfinite(start + duration)) {
    PARSE_ERROR(context,
                "MomentPulse requires finite StartTime and positive "
                "finite Duration.");
  }
  for (size_t a = 0; a < 5; ++a) {
    if (not std::isfinite(electric[a]) or not std::isfinite(magnetic[a])) {
      PARSE_ERROR(context, "MomentPulse components must be finite.");
    }
  }
}
gr::np::TidalMoments detail::RecordedMomentPulse::operator()(
    const double time) const {
  gr::np::TidalMoments result{};
  if (time <= start or time >= start + duration) {
    return result;
  }
  const double envelope =
      pow(sin(std::numbers::pi * (time - start) / duration), 4);
  for (size_t a = 0; a < 5; ++a) {
    result[a] = envelope * std::complex<double>{electric[a], magnetic[a]};
  }
  return result;
}
void detail::RecordedMomentPulse::pup(PUP::er& p) {
  p | start;
  p | duration;
  p | electric;
  p | magnetic;
}
WorldtubeReplayTypeD WorldtubeReplayTypeD::from_yaml(
    const Options::Option& input) {
  // Retain the original three-option TypeD input. A selected model uses
  // explicit mass and relaxation options, including None where appropriate.
  auto node = YAML::Clone(input.node());
  for (const char* key : {"MomentPulse", "MomentDiagnosticsInterval",
                          "MomentFeedback", "DerivativeFitReplay"}) {
    if (node.IsMap() and not node[key]) {
      node[key] = "None";
    }
  }
  if (node.IsMap() and not node["PhysicalTarget"]) {
    node["PhysicalTarget"] = "Live";
  }
  if (node.IsMap() and not node["ReplayBoundaryData"]) {
    node["ReplayBoundaryData"] = true;
  }
  if (node.IsMap() and not node["ReplayBoundaryData"].as<bool>()) {
    if (not node["H5Filename"]) {
      node["H5Filename"] = "";
    }
    if (not node["SubfileName"]) {
      node["SubfileName"] = "";
    }
  }
  if (node.IsMap() and not node["ConstraintPreserving"]) {
    node["ConstraintPreserving"] = false;
  }
  if (node.IsMap() and not node["SommerfeldGauge"]) {
    node["SommerfeldGauge"] = false;
  }
  if (node.IsMap() and not node["QhatRadialResponse"]) {
    node["QhatRadialResponse"] = "None";
  }
  if (node.IsMap() and not node["PhysicalModel"] and not node["Mass"] and
      not node["MomentRelaxationTime"]) {
    node["PhysicalModel"] = "TypeD";
    node["Mass"] = "None";
    node["MomentRelaxationTime"] = "None";
  }
  Options::Parser<options> parser(help);
  parser.parse(Options::Option(std::move(node), input.context()));
  return parser.apply_all([&input](auto /*option_list*/, auto&&... args) {
    return WorldtubeReplayTypeD(std::forward<decltype(args)>(args)...,
                                input.context());
  });
}
WorldtubeReplayTypeD::WorldtubeReplayTypeD(std::string filename,
                                           std::string subfile,
                                           const double radius)
    : WorldtubeReplayTypeD(std::move(filename), std::move(subfile), radius,
                           worldtube::PhysicalModel::TypeD, std::nullopt,
                           std::nullopt) {}
WorldtubeReplayTypeD::WorldtubeReplayTypeD(
    std::string filename, std::string subfile, const double radius,
    const PhysicalModel::type physical_model, const std::optional<double> mass,
    const detail::MomentRelaxation moment_relaxation_time,
    const bool constraint_preserving, const bool sommerfeld_gauge,
    const std::optional<std::array<double, 3>> qhat_radial_response,
    const bool replay_boundary_data, std::string physical_target,
    std::optional<detail::RecordedMomentPulse> moment_pulse,
    const std::optional<double> moment_diagnostics_interval,
    std::optional<detail::MomentFeedback> moment_feedback,
    std::optional<detail::DerivativeFitReplay> derivative_fit_replay,
    const Options::Context& context)
    : replay_(replay_boundary_data
                  ? std::make_optional<WorldtubeReplay>(
                        std::move(filename), std::move(subfile), radius,
                        physical_target != "Live" or
                            derivative_fit_replay.has_value())
                  : std::nullopt),
      model_(constraint_preserving ? detail::SectorImposition::Bjorhus
                                   : detail::SectorImposition::Frozen,
             detail::SectorImposition::Bjorhus,
             sommerfeld_gauge ? detail::SectorImposition::SommerfeldAbsorbing
                              : detail::SectorImposition::Frozen,
             physical_model, mass, moment_relaxation_time, context),
      qhat_radial_response_(qhat_radial_response),
      replay_boundary_data_(replay_boundary_data),
      physical_target_(std::move(physical_target)),
      moment_pulse_(std::move(moment_pulse)),
      moment_diagnostics_interval_(moment_diagnostics_interval),
      moment_feedback_(std::move(moment_feedback)),
      derivative_fit_replay_(std::move(derivative_fit_replay)) {
  if (model_.physical_model() ==
          worldtube::PhysicalModel::QuadrupoleLongitudinal and
      (physical_target_ != "Live" or moment_pulse_ or moment_feedback_ or
       derivative_fit_replay_)) {
    PARSE_ERROR(context,
                "QuadrupoleLongitudinal requires Live with no MomentPulse, "
                "MomentFeedback or DerivativeFitReplay");
  }
  if (derivative_fit_replay_ and
      (physical_target_ != "Live" or not replay_boundary_data_ or
       model_.physical_model() !=
           worldtube::PhysicalModel::QuadrupoleGeometric or
       model_.psi0_fixed_point() or moment_relaxation_time.has_value() or
       moment_feedback_)) {
    PARSE_ERROR(context,
                "DerivativeFitReplay requires Live QuadrupoleGeometric, "
                "ReplayBoundaryData and no relaxation, fixed point or "
                "MomentFeedback.");
  }
  if (moment_feedback_.has_value() and
      (physical_target_ != "RecordedMoments" or
       model_.physical_model() != worldtube::PhysicalModel::Quadrupole)) {
    PARSE_ERROR(
        context,
        "MomentFeedback requires RecordedMoments with legacy Quadrupole.");
  }
  if (model_.psi0_fixed_point() and physical_target_ != "Live") {
    PARSE_ERROR(context,
                "QuadrupoleGeometricFixedPoint requires PhysicalTarget: Live.");
  }
  if ((moment_pulse_.has_value() or
       moment_diagnostics_interval_.has_value()) and
      physical_target_ != "RecordedMoments" and not derivative_fit_replay_ and
      model_.physical_model() !=
          worldtube::PhysicalModel::QuadrupoleLongitudinal) {
    PARSE_ERROR(context,
                "MomentPulse and MomentDiagnosticsInterval require "
                "PhysicalTarget: RecordedMoments or DerivativeFitReplay.");
  }
  if (moment_diagnostics_interval_.has_value() and
      (not std::isfinite(*moment_diagnostics_interval_) or
       *moment_diagnostics_interval_ <= 0.)) {
    PARSE_ERROR(context,
                "MomentDiagnosticsInterval must be positive and finite.");
  }
  if (physical_target_ != "Live" and physical_target_ != "RecordedModel" and
      physical_target_ != "RecordedWeyl" and
      physical_target_ != "RecordedMoments") {
    PARSE_ERROR(context,
                "PhysicalTarget must be Live, RecordedModel, RecordedMoments "
                "or RecordedWeyl.");
  }
  if (physical_target_ != "Live" and
      (not replay_boundary_data_ or moment_relaxation_time.has_value())) {
    PARSE_ERROR(context,
                "Recorded physical targets require ReplayBoundaryData: "
                "true and MomentRelaxationTime: None.");
  }
  if (physical_target_ == "RecordedModel" and
      model_.physical_model() ==
          worldtube::PhysicalModel::ThirdOrderGeometric) {
    PARSE_ERROR(context,
                "ThirdOrderGeometric needs live causal history; use "
                "PhysicalTarget: Live.");
  }
  if (physical_target_ == "RecordedModel" and
      model_.physical_model() == worldtube::PhysicalModel::None) {
    PARSE_ERROR(context, "RecordedModel requires a physical model.");
  }
  if (physical_target_ == "RecordedMoments" and
      not worldtube::is_order_two(model_.physical_model())) {
    PARSE_ERROR(context,
                "RecordedMoments requires an order-two physical model.");
  }
  if (not replay_boundary_data_ and
      (not constraint_preserving or
       (not sommerfeld_gauge and not qhat_radial_response.has_value()))) {
    PARSE_ERROR(context,
                "ReplayBoundaryData: false requires ConstraintPreserving: "
                "true and either SommerfeldGauge or QhatRadialResponse.");
  }
  if (sommerfeld_gauge and qhat_radial_response_.has_value()) {
    PARSE_ERROR(
        context,
        "SommerfeldGauge and QhatRadialResponse are mutually exclusive.");
  }
  if (qhat_radial_response_.has_value()) {
    for (const double coefficient : *qhat_radial_response_) {
      if (not std::isfinite(coefficient)) {
        PARSE_ERROR(context, "QhatRadialResponse coefficients must be finite.");
      }
    }
  }
}
WorldtubeReplayTypeD::WorldtubeReplayTypeD(CkMigrateMessage* msg)
    : BoundaryCondition<3>(msg) {}
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubeReplayTypeD::get_clone() const {
  return std::make_unique<WorldtubeReplayTypeD>(*this);
}
void WorldtubeReplayTypeD::pup(PUP::er& p) {
  BoundaryCondition<3>::pup(p);
  p | replay_;
  p | model_;
  p | qhat_radial_response_;
  p | replay_boundary_data_;
  p | physical_target_;
  p | moment_pulse_;
  p | moment_diagnostics_interval_;
  p | moment_feedback_;
  p | derivative_fit_replay_;
}
gr::np::TidalMoments WorldtubeReplayTypeD::moment_perturbation(
    const double time) const {
  return moment_pulse_.has_value() ? (*moment_pulse_)(time)
                                   : gr::np::TidalMoments{};
}

worldtube::MatchingEvaluation WorldtubeReplayTypeD::derivative_fit_evaluation(
    const double time, const tnsr::aa<DataVector, 3>& g,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi,
    const tnsr::i<DataVector, 3>& normal, const tnsr::I<DataVector, 3>& coords,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const worldtube::KretschmannFaceData<3>& face_data) const {
  ASSERT(derivative_fit_replay_ and replay_, "Derivative replay is disabled.");
  if (face_data.time != time) {
    ERROR("DerivativeFitReplay requires current Kretschmann face data.");
  }
  const auto live = worldtube::weyl_curvature(g, pi, phi, d_pi, d_phi);
  const double strength = (*derivative_fit_replay_)(time);
  std::optional<gr::np::WeylScalars> fit_psi{};
  if (strength != 1.) {
    const auto center_at_time = [&](const double t) {
      return detail::excision_sphere_center(domain.excision_spheres(),
                                            element.id(), t, functions_of_time);
    };
    tnsr::I<double, 3> velocity(0.);
    for (const auto& [name, sphere] : domain.excision_spheres()) {
      (void)name;
      if (sphere.abutting_directions().count(element.id().block_id()) != 0 and
          sphere.is_time_dependent()) {
        velocity =
            std::get<3>(sphere.moving_mesh_grid_to_inertial_map()
                            .coords_frame_velocity_jacobians(
                                sphere.center(), time, functions_of_time));
      }
    }
    const auto reference = replay_->tape().evaluate_physical(
        time, coords, center_at_time, velocity);
    auto mixed_dpi = d_pi;
    auto mixed_dphi = d_phi;
    for (size_t a = 0; a < mixed_dpi.size(); ++a) {
      mixed_dpi[a] =
          reference.d_pi[a] + strength * (d_pi[a] - reference.d_pi[a]);
    }
    for (size_t a = 0; a < mixed_dphi.size(); ++a) {
      mixed_dphi[a] =
          reference.d_phi[a] + strength * (d_phi[a] - reference.d_phi[a]);
    }
    const auto mixed =
        worldtube::weyl_curvature(g, pi, phi, mixed_dpi, mixed_dphi);
    const size_t n = get<0>(normal).size();
    tnsr::I<DataVector, 3> directions(n);
    const DataVector length =
        sqrt(square(get<0>(normal)) + square(get<1>(normal)) +
             square(get<2>(normal)));
    for (size_t i = 0; i < 3; ++i) {
      directions.get(i) = normal.get(i) / length;
    }
    fit_psi = gr::np::weyl_scalars_from_electric_magnetic(
        mixed.electric, mixed.magnetic, live.spatial_metric, directions);
  }
  return worldtube::evaluate_matching_with_fit_input(
      model_.physical_model(), model_.mass(), live.electric, live.magnetic,
      live.spatial_metric, normal, live.lapse, live.shift, &face_data,
      std::nullopt, false, fit_psi ? &*fit_psi : nullptr);
}

void WorldtubeReplayTypeD::observe_moment_response(
    const double time, const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& coordinates, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const worldtube::KretschmannFaceData<3>& face_data) const {
  ASSERT((physical_target_ == "RecordedMoments" or derivative_fit_replay_ or
          model_.physical_model() ==
              worldtube::PhysicalModel::QuadrupoleLongitudinal) and
             face_data.direction,
         "Moment response diagnostics require a recorded-moments face.");
  const size_t dim = face_data.direction->dimension();
  const size_t index = index_to_slice_at(mesh.extents(), *face_data.direction);
  const auto coords = data_on_slice(coordinates, mesh.extents(), dim, index);
  const auto g = data_on_slice(spacetime_metric, mesh.extents(), dim, index);
  const auto p = data_on_slice(pi, mesh.extents(), dim, index);
  const auto f = data_on_slice(phi, mesh.extents(), dim, index);
  const auto dpi = data_on_slice(partial_derivative(pi, mesh, inverse_jacobian),
                                 mesh.extents(), dim, index);
  const auto dphi =
      data_on_slice(partial_derivative(phi, mesh, inverse_jacobian),
                    mesh.extents(), dim, index);
  const auto curvature = worldtube::weyl_curvature(g, p, f, dpi, dphi);
  // Same outward domain normal used by the DG face (inward at excision).
  const auto inv_jac =
      data_on_slice(inverse_jacobian, mesh.extents(), dim, index);
  const size_t n = get<0>(coords).size();
  tnsr::i<DataVector, 3> normal(n);
  DataVector norm(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    normal.get(i) = face_data.direction->sign() * inv_jac.get(dim, i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      norm += normal.get(i) * normal.get(j) *
              curvature.inverse_spatial_metric.get(i, j);
    }
  }
  norm = sqrt(norm);
  for (auto& component : normal) {
    component /= norm;
  }
  if (model_.physical_model() ==
      worldtube::PhysicalModel::QuadrupoleLongitudinal) {
    gr::np::TriadVector direction(n), labels(n);
    const DataVector length = sqrt(
        square(normal.get(0)) + square(normal.get(1)) + square(normal.get(2)));
    for (size_t i = 0; i < 3; ++i) {
      direction.get(i) = normal.get(i) / length;
      labels.get(i) = -direction.get(i);
    }
    const auto rotation =
        gr::np::adapted_triad(curvature.spatial_metric, direction);
    const auto psi = gr::np::weyl_scalars_from_electric_magnetic(
        curvature.electric, curvature.magnetic, curvature.spatial_metric,
        direction);
    const size_t grid_lmax =
        static_cast<size_t>((sqrt(8. * n + 1.) - 3.) / 4. + .5);
    const char* setting = std::getenv("NP_GEOMETRIC_LMAX");
    const size_t lm = setting == nullptr ? std::min(size_t{8}, grid_lmax)
                                         : std::stoul(setting);
    const auto fit = gr::np::evaluate_longitudinal_tide(
        psi, rotation, curvature.spatial_metric, labels,
        face_data.quadrature_weights, *model_.mass(), lm, curvature.lapse,
        curvature.shift, face_data.d_kretschmann, face_data.dt_kretschmann);
    const auto axes = gr::np::sphere_map_rotation(
        fit.frame.map.direction, labels, face_data.quadrature_weights);
    double target = 0., area = 0.;
    for (size_t i = 0; i < n; ++i) {
      target +=
          face_data.quadrature_weights[i] * std::norm(fit.completed.get(0)[i]);
      area += face_data.quadrature_weights[i];
    }
    std::ostringstream out;
    out << std::setprecision(17) << "LONGITUDINAL_FIT time=" << time
        << " element=" << element.id() << " moments=";
    for (size_t a = 0; a < 5; ++a)
      out << (a == 0 ? "" : ",") << fit.moments[a].real();
    for (const auto& value : fit.moments)
      out << ',' << value.imag();
    out << " residual=" << fit.relative_residual
        << " condition=" << fit.normalized_condition_number
        << " seed_margin=" << fit.minimum_relative_denominator
        << " gradient_tilt=" << fit.maximum_gradient_tilt
        << " target_rms=" << sqrt(target / area) << " axes_to_inertial=";
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        out << (i == 0 and j == 0 ? "" : ",") << axes[i][j];
    std::fprintf(stdout, "%s\n", out.str().c_str());
    std::fflush(stdout);
    return;
  }
  const auto raw = worldtube::evaluate_matching(
      model_.physical_model(), model_.mass(), curvature.electric,
      curvature.magnetic, curvature.spatial_metric, normal, curvature.lapse,
      curvature.shift, &face_data, std::nullopt);
  if (derivative_fit_replay_) {
    const auto fitted = derivative_fit_evaluation(
        time, g, p, f, dpi, dphi, normal, coords, domain, element,
        functions_of_time, face_data);
    const auto pulse = moment_perturbation(time);
    auto imposed = fitted.second_order->fit.components;
    for (size_t a = 0; a < 5; ++a) {
      imposed[a] += pulse[a];
    }
    const auto applied = worldtube::evaluate_matching(
        model_.physical_model(), model_.mass(), curvature.electric,
        curvature.magnetic, curvature.spatial_metric, normal, curvature.lapse,
        curvature.shift, &face_data, imposed);
    auto retained = raw.registration->pulled_back;
    retained.get(0) = 0.;
    const auto keep =
        gr::np::push_forward(retained, raw.registration->rotation.a_bar,
                             raw.registration->rotation.b)
            .get(0);
    // Reconstruct the whole derivative contribution jointly, with live
    // coefficients. Zero derivative inputs remove exactly P from curvature.
    const auto rest =
        worldtube::weyl_curvature(g, p, f, tnsr::iaa<DataVector, 3>(n, 0.),
                                  tnsr::ijaa<DataVector, 3>(n, 0.));
    tnsr::I<DataVector, 3> directions(n);
    const DataVector length =
        sqrt(square(get<0>(normal)) + square(get<1>(normal)) +
             square(get<2>(normal)));
    for (size_t i = 0; i < 3; ++i) {
      directions.get(i) = normal.get(i) / length;
    }
    auto principal = raw.psi;
    const auto rest_psi = gr::np::weyl_scalars_from_electric_magnetic(
        rest.electric, rest.magnetic, curvature.spatial_metric, directions);
    for (size_t a = 0; a < 5; ++a) {
      principal.get(a) -= rest_psi.get(a);
    }
    const auto principal_fit = worldtube::evaluate_matching_with_fit_input(
        model_.physical_model(), model_.mass(), curvature.electric,
        curvature.magnetic, curvature.spatial_metric, normal, curvature.lapse,
        curvature.shift, &face_data, std::nullopt, false, &principal);
    for (auto& component : directions) {
      component *= -1.;
    }
    const size_t grid_lmax = static_cast<size_t>(
        (sqrt(8. * static_cast<double>(n) + 1.) - 3.) / 4. + 0.5);
    const char* setting = std::getenv("NP_GEOMETRIC_LMAX");
    const size_t lm = setting == nullptr
                          ? std::min(size_t{8}, grid_lmax)
                          : static_cast<size_t>(std::stoul(setting));
    const auto frame = gr::np::geometric_frame(
        *raw.registration, *raw.rapidity, raw.adapted_rotation,
        curvature.spatial_metric, directions, face_data.quadrature_weights,
        *model_.mass(), lm);
    const auto axes_to_inertial = gr::np::sphere_map_rotation(
        frame.map.direction, directions, face_data.quadrature_weights);
    const auto rms = [&face_data](const ComplexDataVector& v) {
      double value = 0.;
      double area = 0.;
      for (size_t i = 0; i < v.size(); ++i) {
        value += face_data.quadrature_weights[i] * std::norm(v[i]);
        area += face_data.quadrature_weights[i];
      }
      return sqrt(value / area);
    };
    std::ostringstream out;
    out << std::setprecision(17) << "DERIVATIVE_FIT time=" << time
        << " element=" << element.id()
        << " strength=" << (*derivative_fit_replay_)(time);
    const auto append = [&out](const char* name,
                               const gr::np::TidalMoments& q) {
      out << ' ' << name << '=';
      for (size_t a = 0; a < 5; ++a) {
        out << (a == 0 ? "" : ",") << q[a].real();
      }
      for (const auto& x : q) {
        out << ',' << x.imag();
      }
    };
    append("shadow", raw.second_order->fit.components);
    append("fitted", fitted.second_order->fit.components);
    append("imposed", imposed);
    append("pulse", pulse);
    append("principal", principal_fit.second_order->fit.components);
    out << " axes_to_inertial=";
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        out << (i == 0 and j == 0 ? "" : ",") << axes_to_inertial[i][j];
      }
    }
    out << " shadow_residual=" << raw.second_order->fit.relative_residual
        << " fitted_residual=" << fitted.second_order->fit.relative_residual
        << " keep_rms=" << rms(keep) << " fit_rms="
        << rms(ComplexDataVector{get(applied.psi0_target) - keep})
        << " target_rms=" << rms(get(applied.psi0_target))
        << " shadow_target_rms=" << rms(get(raw.psi0_target))
        << " psi4_rms=" << rms(raw.psi.get(4))
        << " principal_psi4_rms=" << rms(principal.get(4));
    std::fprintf(stdout, "%s\n", out.str().c_str());
    std::fflush(stdout);
    return;
  }
  const auto recorded = recorded_tidal_moments(time, coords, domain, element,
                                               functions_of_time, mesh);
  const auto pulse = moment_perturbation(time);
  auto imposed = recorded;
  if (moment_feedback_) {
    imposed = moment_feedback_->imposed(
        time, raw.second_order->fit.components, recorded, pulse,
        curvature.spatial_metric, coords, face_data.quadrature_weights);
  } else {
    for (size_t a = 0; a < 5; ++a) {
      imposed[a] += pulse[a];
    }
  }
  std::ostringstream out;
  out << std::setprecision(17) << "MOMENT_RESPONSE time=" << time
      << " element=" << element.id();
  const auto append = [&out](const char* name, const gr::np::TidalMoments& q) {
    out << ' ' << name << '=';
    for (size_t a = 0; a < 5; ++a) {
      out << (a == 0 ? "" : ",") << q[a].real();
    }
    for (const auto& component : q) {
      out << ',' << component.imag();
    }
  };
  append("recorded", recorded);
  append("pulse", pulse);
  append("imposed", imposed);
  append("inferred", raw.second_order->fit.components);
  if (moment_feedback_) {
    out << " feedback_strength=" << moment_feedback_->strength(time);
    const auto gains = moment_feedback_->strengths(time);
    out << " feedback_electric_strength=" << gains[0]
        << " feedback_magnetic_strength=" << gains[1];
  }
  out << " residual=" << raw.second_order->fit.relative_residual;
  // A single stdio call prevents interleaving among threaded elements.
  std::fprintf(stdout, "%s\n", out.str().c_str());
  std::fflush(stdout);
}
tnsr::ii<DataVector, 3> WorldtubeReplayTypeD::recorded_incoming_mode(
    const double time, const tnsr::I<DataVector, 3>& coords,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const Mesh<3>& mesh) const {
  return recorded_target_data(time, coords, domain, element, functions_of_time,
                              mesh)
      .first;
}

gr::np::TidalMoments WorldtubeReplayTypeD::recorded_tidal_moments(
    const double time, const tnsr::I<DataVector, 3>& coords,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const Mesh<3>& mesh) const {
  auto result = recorded_target_data(time, coords, domain, element,
                                     functions_of_time, mesh);
  if (not result.second.has_value()) {
    ERROR("Recorded tidal moments require a recorded order-two model fit.");
  }
  return *result.second;
}

std::pair<tnsr::ii<DataVector, 3>, std::optional<gr::np::TidalMoments>>
WorldtubeReplayTypeD::recorded_target_data(
    const double time, const tnsr::I<DataVector, 3>& coords,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const Mesh<3>& mesh) const {
  if (physical_target_ == "Live" or not replay_.has_value()) {
    ERROR("recorded_incoming_mode requires a recorded physical target.");
  }
  const auto center_at_time = [&](const double t) {
    return detail::excision_sphere_center(domain.excision_spheres(),
                                          element.id(), t, functions_of_time);
  };
  tnsr::I<double, 3> center_velocity(0.);
  const auto direction =
      worldtube::excision_face_direction(domain.excision_spheres(), element);
  if (not direction.has_value()) {
    ERROR("Recorded physical targets require an excision face.");
  }
  for (const auto& [name, sphere] : domain.excision_spheres()) {
    (void)name;
    if (sphere.abutting_directions().count(element.id().block_id()) != 0 and
        sphere.is_time_dependent()) {
      center_velocity =
          std::get<3>(sphere.moving_mesh_grid_to_inertial_map()
                          .coords_frame_velocity_jacobians(
                              sphere.center(), time, functions_of_time));
    }
  }
  const auto data = replay_->tape().evaluate_physical(
      time, coords, center_at_time, center_velocity);
  const auto curvature = worldtube::weyl_curvature(
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(data.fields),
      get<gh::Tags::Pi<DataVector, 3>>(data.fields),
      get<gh::Tags::Phi<DataVector, 3>>(data.fields), data.d_pi, data.d_phi);
  const auto center = center_at_time(time);
  const size_t n = get<0>(coords).size();
  tnsr::i<DataVector, 3> normal(n);
  DataVector norm(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    normal.get(i) = center.get(i) - coords.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      norm += normal.get(i) * normal.get(j) *
              curvature.inverse_spatial_metric.get(i, j);
    }
  }
  norm = sqrt(norm);
  for (auto& component : normal) {
    component /= norm;
  }
  if (physical_target_ != "RecordedWeyl") {
    const auto face_mesh = mesh.slice_away(direction->dimension());
    if (worldtube::is_order_two(model_.physical_model()) and
        (face_mesh.basis(0) != Spectral::Basis::SphericalHarmonic or
         face_mesh.basis(1) != Spectral::Basis::SphericalHarmonic)) {
      ERROR(
          "Recorded quadrupole targets require a full spherical-harmonic "
          "face.");
    }
    worldtube::KretschmannFaceData<3> recorded_face{};
    recorded_face.direction = direction;
    recorded_face.time = time;
    recorded_face.kretschmann = data.kretschmann;
    recorded_face.d_kretschmann = data.d_kretschmann;
    recorded_face.dt_kretschmann = data.dt_kretschmann;
    recorded_face.quadrature_weights =
        worldtube::face_quadrature_weights<3>(face_mesh);
    // In particular, do not pass the receiver's filtered_moments: even with
    // relaxation disabled they contain the instantaneous *live* fit.
    auto evaluation = worldtube::evaluate_matching(
        model_.physical_model(), model_.mass(), curvature.electric,
        curvature.magnetic, curvature.spatial_metric, normal, curvature.lapse,
        curvature.shift, &recorded_face, std::nullopt);
    return {std::move(evaluation.incoming_mode),
            evaluation.second_order.has_value()
                ? std::optional{evaluation.second_order->fit.components}
                : std::nullopt};
  }

  tnsr::I<DataVector, 3> directions(n);
  norm = sqrt(square(get<0>(normal)) + square(get<1>(normal)) +
              square(get<2>(normal)));
  for (size_t i = 0; i < 3; ++i) {
    directions.get(i) = normal.get(i) / norm;
  }
  const auto rotation =
      gr::np::adapted_triad(curvature.spatial_metric, directions);
  const auto psi = gr::np::weyl_scalars_from_electric_magnetic(
      curvature.electric, curvature.magnetic, curvature.spatial_metric,
      directions);
  Scalar<ComplexDataVector> psi0;
  get(psi0) = psi.get(0);
  auto incoming = gr::np::orthonormal_to_coordinate_covariant(
      gr::np::incoming_weyl_field(psi0, rotation),
      gr::np::cholesky_factor(curvature.spatial_metric));
  // Match evaluate_matching and the Bjorhus physical-sector convention.
  for (auto& component : incoming) {
    component *= 0.5;
  }
  return {std::move(incoming), std::nullopt};
}
std::optional<std::string> WorldtubeReplayTypeD::dg_ghost(
    gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> ghost_phi,
    gsl::not_null<Scalar<DataVector>*> ghost_gamma1,
    gsl::not_null<Scalar<DataVector>*> ghost_gamma2,
    gsl::not_null<Scalar<DataVector>*> ghost_lapse,
    gsl::not_null<tnsr::I<DataVector, 3>*> ghost_shift,
    gsl::not_null<tnsr::II<DataVector, 3>*> ghost_inverse_spatial_metric,
    const std::optional<tnsr::I<DataVector, 3>>& face_mesh_velocity,
    const tnsr::i<DataVector, 3>& normal_covector,
    const tnsr::I<DataVector, 3>& normal_vector,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& coords, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::AA<DataVector, 3>& /*inverse_spacetime_metric*/,
    const tnsr::A<DataVector, 3>& /*spacetime_unit_normal_vector*/,
    const tnsr::iaa<DataVector, 3>& /*three_index_constraint*/,
    const tnsr::a<DataVector, 3>& /*gauge_source*/,
    const tnsr::ab<DataVector, 3>& /*spacetime_deriv_gauge_source*/,
    const tnsr::aa<DataVector, 3>& /*logical_dt_spacetime_metric*/,
    const tnsr::aa<DataVector, 3>& /*logical_dt_pi*/,
    const tnsr::iaa<DataVector, 3>& /*logical_dt_phi*/,
    const tnsr::iaa<DataVector, 3>& /*d_spacetime_metric*/,
    const tnsr::iaa<DataVector, 3>& /*d_pi*/,
    const tnsr::ijaa<DataVector, 3>& /*d_phi*/, double time,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const Mesh<3>& /*mesh*/, ::dg::Formulation /*formulation*/,
    const evolution::BoundaryCorrection& /*boundary_correction*/,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& /*volume_inverse_jacobian*/,
    const worldtube::KretschmannFaceData<3>& /*face_data*/) const {
  if (not replay_boundary_data_) {
    // Identical states have exactly zero AveragedUpwindPenalty, so the
    // time-derivative method supplies the entire boundary correction.
    *ghost_metric = spacetime_metric;
    *ghost_pi = pi;
    *ghost_phi = phi;
    *ghost_gamma1 = gamma1;
    *ghost_gamma2 = gamma2;
    *ghost_lapse = lapse;
    *ghost_shift = shift;
    *ghost_inverse_spatial_metric =
        determinant_and_inverse(gr::spatial_metric(spacetime_metric)).second;
    return {};
  }
  return replay_->dg_ghost(ghost_metric, ghost_pi, ghost_phi, ghost_gamma1,
                           ghost_gamma2, ghost_lapse, ghost_shift,
                           ghost_inverse_spatial_metric, face_mesh_velocity,
                           normal_covector, normal_vector, coords, gamma1,
                           gamma2, time, domain, element, functions_of_time);
}
std::optional<std::string> WorldtubeReplayTypeD::dg_time_derivative(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> dt_phi,
    const std::optional<tnsr::I<DataVector, 3>>& face_mesh_velocity,
    const tnsr::i<DataVector, 3>& normal_covector,
    const tnsr::I<DataVector, 3>& normal_vector,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& coords, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::AA<DataVector, 3>& inverse_spacetime_metric,
    const tnsr::A<DataVector, 3>& spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, 3>& three_index_constraint,
    const tnsr::a<DataVector, 3>& gauge_source,
    const tnsr::ab<DataVector, 3>& spacetime_deriv_gauge_source,
    const tnsr::aa<DataVector, 3>& logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& logical_dt_pi,
    const tnsr::iaa<DataVector, 3>& logical_dt_phi,
    const tnsr::iaa<DataVector, 3>& d_spacetime_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi, double time,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time, const Mesh<3>& mesh,
    ::dg::Formulation formulation,
    const evolution::BoundaryCorrection& boundary_correction,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& volume_inverse_jacobian,
    const worldtube::KretschmannFaceData<3>& face_data,
    const bool nonphysical_only) const {
  const auto* flux =
      dynamic_cast<const gh::BoundaryCorrections::AveragedUpwindPenalty<3>*>(
          &boundary_correction);
  if (flux == nullptr) {
    return "WorldtubeReplayTypeD requires AveragedUpwindPenalty.";
  }
  std::optional<Direction<3>> direction;
  for (const auto& [name, sphere] : domain.excision_spheres()) {
    (void)name;
    const auto found =
        sphere.abutting_directions().find(element.id().block_id());
    if (found != sphere.abutting_directions().end()) {
      direction = found->second;
    }
  }
  if (not direction.has_value()) {
    return "WorldtubeReplayTypeD requires an excision face.";
  }
  const size_t dim = direction->dimension();
  if (mesh.quadrature(dim) != Spectral::Quadrature::GaussLobatto) {
    return "WorldtubeReplayTypeD requires GaussLobatto in the face-normal "
           "direction.";
  }

  const size_t n = get<0>(coords).size();
  // Use the selected curvature prescription and moving-mesh Bjorhus
  // correction. Only the requested characteristic sectors are retained below.
  std::optional<tnsr::ii<DataVector, 3>> recorded_target{};
  std::optional<gr::np::TidalMoments> recorded_moments{};
  const auto diagnostic_pulse = moment_perturbation(time);
  bool has_diagnostic_pulse = false;
  for (const auto& component : diagnostic_pulse) {
    has_diagnostic_pulse = has_diagnostic_pulse or component != 0.;
  }
  // Preserve the original stage path exactly at unit strength without a pulse,
  // including its live prefit. The full-step shadow diagnostics are read-only.
  if (derivative_fit_replay_ and
      ((*derivative_fit_replay_)(time) != 1. or has_diagnostic_pulse)) {
    const auto fitted = derivative_fit_evaluation(
        time, spacetime_metric, pi, phi, d_pi, d_phi, normal_covector, coords,
        domain, element, functions_of_time, face_data);
    recorded_moments = fitted.second_order->fit.components;
    const auto pulse = moment_perturbation(time);
    for (size_t a = 0; a < 5; ++a) {
      (*recorded_moments)[a] += pulse[a];
    }
  } else if (physical_target_ == "RecordedMoments") {
    recorded_moments = recorded_tidal_moments(time, coords, domain, element,
                                              functions_of_time, mesh);
    const auto pulse = moment_perturbation(time);
    if (moment_feedback_) {
      if (face_data.time != time) {
        ERROR("MomentFeedback requires current Kretschmann face data.");
      }
      // RecordedMoments deliberately bypasses the action's live prefit.
      // Fit the current face here from the same volume derivatives passed
      // to Bjorhus, and retain the current action's invariant-gradient data.
      const auto curvature =
          worldtube::weyl_curvature(spacetime_metric, pi, phi, d_pi, d_phi);
      const auto raw = worldtube::evaluate_matching(
          model_.physical_model(), model_.mass(), curvature.electric,
          curvature.magnetic, curvature.spatial_metric, normal_covector,
          curvature.lapse, curvature.shift, &face_data, std::nullopt);
      *recorded_moments = moment_feedback_->imposed(
          time, raw.second_order->fit.components, *recorded_moments, pulse,
          curvature.spatial_metric, coords, face_data.quadrature_weights);
    } else {
      for (size_t a = 0; a < 5; ++a) {
        (*recorded_moments)[a] += pulse[a];
      }
    }
  } else if (physical_target_ != "Live") {
    recorded_target = recorded_incoming_mode(time, coords, domain, element,
                                             functions_of_time, mesh);
  }
  Bjorhus::IntermediateVariables<3> vars{n};
  if (auto error = model_.dg_time_derivative_with_intermediates(
          make_not_null(&vars), dt_metric, dt_pi, dt_phi, face_mesh_velocity,
          normal_covector, normal_vector, spacetime_metric, pi, phi, coords,
          gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
          spacetime_unit_normal_vector, three_index_constraint, gauge_source,
          spacetime_deriv_gauge_source, logical_dt_spacetime_metric,
          logical_dt_pi, logical_dt_phi, d_spacetime_metric, d_pi, d_phi, time,
          domain, element, functions_of_time, face_data,
          ConstraintDamping::Constant<3, Frame::Grid>{0.},
          recorded_target.has_value() ? &*recorded_target : nullptr,
          recorded_moments);
      error.has_value()) {
    return error;
  }

  if (not replay_boundary_data_ and min(vars.char_speeds[2]) < 0.) {
    return "ReplayBoundaryData: false requires nonnegative u-plus speeds; "
           "there is no model for incoming u-plus.";
  }

  Variables<worldtube::ReplayTags> reference(n, 0.);
  auto& rg = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(reference);
  auto& rp = get<gh::Tags::Pi<DataVector, 3>>(reference);
  auto& rf = get<gh::Tags::Phi<DataVector, 3>>(reference);
  if (replay_boundary_data_) {
    // Reproduce the native replay lift, including distinct interior/exterior
    // normal normalization, receiver mesh velocity, and the *element* inverse
    // Jacobian (so radial h-refinement contributes its scale).
    Variables<worldtube::ReplayTags> ghost(n);
    auto& gm = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(ghost);
    auto& gp = get<gh::Tags::Pi<DataVector, 3>>(ghost);
    auto& gf = get<gh::Tags::Phi<DataVector, 3>>(ghost);
    Scalar<DataVector> out_gamma1, out_gamma2, out_lapse;
    tnsr::I<DataVector, 3> out_shift;
    tnsr::II<DataVector, 3> out_inverse;
    if (auto error = replay_->dg_ghost(
            make_not_null(&gm), make_not_null(&gp), make_not_null(&gf),
            make_not_null(&out_gamma1), make_not_null(&out_gamma2),
            make_not_null(&out_lapse), make_not_null(&out_shift),
            make_not_null(&out_inverse), face_mesh_velocity, normal_covector,
            normal_vector, coords, gamma1, gamma2, time, domain, element,
            functions_of_time);
        error.has_value()) {
      return error;
    }
    const auto inverse_jacobian = data_on_slice(
        volume_inverse_jacobian, mesh.extents(), dim,
        direction->side() == Side::Lower ? 0 : mesh.extents(dim) - 1);
    Scalar<DataVector> magnitude(n, 0.);
    DataVector donor_magnitude(n, 0.);
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        get(magnitude) += inverse_jacobian.get(dim, i) *
                          vars.inverse_spatial_metric.get(i, j) *
                          inverse_jacobian.get(dim, j);
        donor_magnitude += normal_covector.get(i) * out_inverse.get(i, j) *
                           normal_covector.get(j);
      }
    }
    get(magnitude) = sqrt(get(magnitude));
    donor_magnitude = sqrt(donor_magnitude);
    tnsr::i<DataVector, 3> donor_normal(n);
    for (size_t i = 0; i < 3; ++i) {
      donor_normal.get(i) = -normal_covector.get(i) / donor_magnitude;
    }
    const auto velocity =
        face_mesh_velocity.value_or(tnsr::I<DataVector, 3>(n, 0.));
    flux->dg_boundary_terms(make_not_null(&rg), make_not_null(&rp),
                            make_not_null(&rf), spacetime_metric, pi, phi,
                            gamma1, gamma2, normal_covector, velocity, gm, gp,
                            gf, out_gamma1, out_gamma2, donor_normal, velocity,
                            formulation);
    ::dg::lift_flux(make_not_null(&reference), mesh.extents(dim), magnitude,
                    mesh.basis(dim));
  }

  // Native DG adds the reference lift (zero without boundary replay).
  // Supply the selected incoming projections of (C_model - C_reference).
  tnsr::aa<DataVector, 3> difference(n), selected(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      difference.get(a, b) =
          dt_pi->get(a, b) - rp.get(a, b) -
          get(gamma2) * (dt_metric->get(a, b) - rg.get(a, b));
      for (size_t i = 0; i < 3; ++i) {
        difference.get(a, b) -= vars.unit_interface_normal_vector.get(i) *
                                (dt_phi->get(i, a, b) - rf.get(i, a, b));
      }
    }
  }
  Bjorhus::detail::add_physical_sector_projection(
      make_not_null(&selected), DataVector(n, 1.), vars.projection_ab,
      vars.projection_Ab, vars.projection_AB, difference);
  if (sommerfeld_gauge()) {
    Bjorhus::detail::add_gauge_sector_projection(
        make_not_null(&selected), DataVector(n, 1.),
        vars.incoming_null_one_form, vars.outgoing_null_one_form,
        vars.incoming_null_vector, vars.outgoing_null_vector,
        vars.projection_Ab, difference);
  }
  if (qhat_radial_response_.has_value()) {
    if (not face_data.qhat_radial_gauge.has_value()) {
      return "QhatRadialResponse requires initialized receiver face data.";
    }
    // This is the metric correction that native replay plus our selected CP
    // replacement will actually impose. Account for its null-basis variation.
    auto total_metric_correction = rg;
    if (constraint_preserving()) {
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = a; b < 4; ++b) {
          for (size_t p = 0; p < n; ++p) {
            if (vars.char_speeds[0][p] < 0.) {
              total_metric_correction.get(a, b)[p] = dt_metric->get(a, b)[p];
            }
          }
        }
      }
    }
    const auto center = detail::excision_sphere_center(
        domain.excision_spheres(), element.id(), time, functions_of_time);
    auto correction = worldtube::qhat_radial_gauge_correction(
        *face_data.qhat_radial_gauge, *qhat_radial_response_, spacetime_metric,
        pi, phi, d_spacetime_metric, d_pi, d_phi, coords, center, lapse, shift,
        face_mesh_velocity, gamma2, total_metric_correction, time,
        functions_of_time);
    // Remove only the gauge contraction of the native reference lift.
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        DataVector reference_minus = rp.get(a, b) - get(gamma2) * rg.get(a, b);
        for (size_t i = 0; i < 3; ++i) {
          reference_minus -=
              vars.unit_interface_normal_vector.get(i) * rf.get(i, a, b);
        }
        correction.get(a) -= vars.outgoing_null_vector.get(b) * reference_minus;
      }
    }
    DataVector ell_correction(n, 0.);
    for (size_t a = 0; a < 4; ++a) {
      ell_correction += vars.outgoing_null_vector.get(a) * correction.get(a);
    }
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        selected.get(a, b) -=
            vars.incoming_null_one_form.get(a) * correction.get(b) +
            vars.incoming_null_one_form.get(b) * correction.get(a) +
            vars.incoming_null_one_form.get(a) *
                vars.incoming_null_one_form.get(b) * ell_correction;
      }
    }
  }
  tnsr::aa<DataVector, 3> selected_psi(n, 0.), selected_plus(n, 0.);
  tnsr::iaa<DataVector, 3> selected_zero(n, 0.);
  if (constraint_preserving()) {
    Bjorhus::detail::add_constraint_sector_projection(
        make_not_null(&selected), DataVector(n, 1.),
        vars.outgoing_null_one_form, vars.incoming_null_vector,
        vars.projection_ab, vars.projection_Ab, vars.projection_AB, difference);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        selected_psi.get(a, b) = dt_metric->get(a, b) - rg.get(a, b);
        for (size_t i = 0; i < 3; ++i) {
          selected_zero.get(i, a, b) = dt_phi->get(i, a, b) - rf.get(i, a, b);
          for (size_t j = 0; j < 3; ++j) {
            selected_zero.get(i, a, b) -=
                normal_covector.get(i) *
                vars.unit_interface_normal_vector.get(j) *
                (dt_phi->get(j, a, b) - rf.get(j, a, b));
          }
        }
      }
    }
  }
  // Zero-speed fields receive no model correction. The reference vanishes
  // without boundary replay. The Bjorhus helper below masks
  // strictly positive speeds and performs the characteristic inversion.
  const auto mask_zero_speed = [](const auto field, const DataVector& speed) {
    for (auto& component : *field) {
      for (size_t p = 0; p < speed.size(); ++p) {
        if (speed[p] == 0.) {
          component[p] = 0.;
        }
      }
    }
  };
  mask_zero_speed(make_not_null(&selected_psi), vars.char_speeds[0]);
  mask_zero_speed(make_not_null(&selected_zero), vars.char_speeds[1]);
  mask_zero_speed(make_not_null(&selected), vars.char_speeds[3]);
  Bjorhus::project_corrections_onto_evolved_variables(
      dt_metric, dt_pi, dt_phi, make_not_null(&selected_psi),
      make_not_null(&selected_zero), make_not_null(&selected_plus),
      make_not_null(&selected), vars.char_speeds, gamma2, normal_covector);
  if (nonphysical_only) {
    // This caller supplies a purely physical ghost instead of the replay
    // ghost. Restore the full native replay lift, then remove its total
    // incoming physical projection. All other sectors retain R + C_model.
    for (size_t c = 0; c < dt_metric->size(); ++c) {
      (*dt_metric)[c] += rg[c];
      (*dt_pi)[c] += rp[c];
    }
    for (size_t c = 0; c < dt_phi->size(); ++c) {
      (*dt_phi)[c] += rf[c];
    }
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        difference.get(a, b) =
            dt_pi->get(a, b) - get(gamma2) * dt_metric->get(a, b);
        for (size_t i = 0; i < 3; ++i) {
          difference.get(a, b) -= normal_vector.get(i) * dt_phi->get(i, a, b);
        }
      }
    }
    selected = tnsr::aa<DataVector, 3>(n, 0.);
    Bjorhus::detail::add_physical_sector_projection(
        make_not_null(&selected), DataVector(n, 1.), vars.projection_ab,
        vars.projection_Ab, vars.projection_AB, difference);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        dt_pi->get(a, b) -= 0.5 * selected.get(a, b);
        for (size_t i = 0; i < 3; ++i) {
          dt_phi->get(i, a, b) +=
              0.5 * normal_covector.get(i) * selected.get(a, b);
        }
      }
    }
  }
  return {};
}
PUP::able::PUP_ID WorldtubeReplayTypeD::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::BoundaryConditions
