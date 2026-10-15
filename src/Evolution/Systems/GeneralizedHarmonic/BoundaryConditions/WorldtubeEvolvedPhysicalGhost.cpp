// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeEvolvedPhysicalGhost.hpp"
#include <cmath>
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Options/ParseError.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"

namespace gh::BoundaryConditions {
detail::RecordedWeylTarget::RecordedWeylTarget(std::string filename,
                                               std::string subfile,
                                               const double radius)
    : boundary(std::move(filename), std::move(subfile), radius,
               worldtube::PhysicalModel::None, std::nullopt, std::nullopt, true,
               false, std::nullopt, true, "RecordedWeyl") {}
void detail::RecordedWeylTarget::pup(PUP::er& p) { p | boundary; }

detail::PrescribedSchwarzschildQuadrupole::PrescribedSchwarzschildQuadrupole(
    const double mass_in, RecordedMomentPulse pulse_in,
    const Options::Context& context)
    : mass(mass_in), pulse(std::move(pulse_in)) {
  if (not std::isfinite(mass) or mass <= 0.) {
    PARSE_ERROR(context,
                "Prescribed quadrupole requires a positive finite mass.");
  }
}

void detail::PrescribedSchwarzschildQuadrupole::pup(PUP::er& p) {
  p | mass;
  p | pulse;
}

tnsr::ii<DataVector, 3>
detail::PrescribedSchwarzschildQuadrupole::incoming_mode(
    const double time, const tnsr::I<DataVector, 3>& coordinates) const {
  const size_t n = get<0>(coordinates).size();
  const DataVector radius =
      sqrt(square(get<0>(coordinates)) + square(get<1>(coordinates)) +
           square(get<2>(coordinates)));
  if (min(radius) <= 2. * mass) {
    ERROR("Prescribed quadrupole requires a sphere outside r=2M.");
  }
  const gr::Solutions::KerrSchild background(mass,
                                             std::array<double, 3>{0., 0., 0.},
                                             std::array<double, 3>{0., 0., 0.});
  const auto fields =
      background.variables(coordinates, time,
                           tmpl::list<gr::Tags::SpatialMetric<DataVector, 3>,
                                      gr::Tags::Lapse<DataVector>,
                                      gr::Tags::Shift<DataVector, 3>>{});
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataVector, 3>>(fields);
  // Exact Schwarzschild curvature in the spherical Kerr-Schild slice:
  // E_ij=M/r^3 (gamma_ij-3 s_i s_j), B_ij=0, K=48 M^2/r^6.
  // The domain normal points INTO the hole, matching evaluate_matching.
  tnsr::i<DataVector, 3> normal(n);
  for (size_t i = 0; i < 3; ++i) {
    normal.get(i) =
        -sqrt(1. + 2. * mass / radius) * coordinates.get(i) / radius;
  }
  tnsr::ii<DataVector, 3> electric(n), magnetic(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      electric.get(i, j) =
          mass / pow(radius, 3) *
          (spatial_metric.get(i, j) - 3. * normal.get(i) * normal.get(j));
    }
  }
  worldtube::KretschmannFaceData<3> face{};
  face.direction = Direction<3>::lower_xi();
  face.time = time;
  get(face.kretschmann) = 48. * square(mass) / pow(radius, 6);
  get(face.dt_kretschmann) = DataVector(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    face.d_kretschmann.get(i) =
        -6. * get(face.kretschmann) * coordinates.get(i) / square(radius);
  }
  // Reuse the production NP model without changing its radial response,
  // boost, or w^-/2 normalization. The imposed moments bypass its fit.
  return worldtube::evaluate_matching(
             worldtube::PhysicalModel::Quadrupole, mass, electric, magnetic,
             spatial_metric, normal, get<gr::Tags::Lapse<DataVector>>(fields),
             get<gr::Tags::Shift<DataVector, 3>>(fields), &face, pulse(time))
      .incoming_mode;
}

WorldtubeEvolvedPhysicalGhost WorldtubeEvolvedPhysicalGhost::from_yaml(
    const Options::Option& input) {
  auto node = YAML::Clone(input.node());
  if (node.IsMap() and not node["PrescribedQuadrupole"]) {
    node["PrescribedQuadrupole"] = "None";
  }
  if (node.IsMap()) {
    if (not node["RecordedWeyl"])
      node["RecordedWeyl"] = "None";
    if (not node["GaugeSector"])
      node["GaugeSector"] = "Frozen";
    if (not node["Amplitude"])
      node["Amplitude"] = 0.;
    if (not node["Frequency"])
      node["Frequency"] = 0.;
  }
  Options::Parser<options> parser(help);
  parser.parse(Options::Option(std::move(node), input.context()));
  return parser.apply_all([&input](auto /*option_list*/, auto&&... args) {
    return WorldtubeEvolvedPhysicalGhost(std::forward<decltype(args)>(args)...,
                                         input.context());
  });
}

WorldtubeEvolvedPhysicalGhost::WorldtubeEvolvedPhysicalGhost(
    const detail::SectorImposition gauge_sector, const double amplitude,
    const double frequency,
    std::optional<detail::PrescribedSchwarzschildQuadrupole>
        prescribed_quadrupole,
    std::optional<detail::RecordedWeylTarget> recorded_weyl,
    const Options::Context& context)
    // RecordedWeyl supplies its own nonphysical sectors. Do not retain an
    // unused analytic boundary: binary runs need not register that analytic
    // solution for Charm serialization.
    : nonphysical_(recorded_weyl.has_value()
                       ? std::nullopt
                       : std::make_optional<WorldtubePhysicalGhost>(
                             std::make_unique<gh::Solutions::WrappedGr<
                                 gr::Solutions::Minkowski<3>>>(),
                             gauge_sector, context)),
      amplitude_(amplitude),
      frequency_(frequency),
      prescribed_quadrupole_(std::move(prescribed_quadrupole)),
      recorded_weyl_(std::move(recorded_weyl)) {
  if (not std::isfinite(amplitude) or not std::isfinite(frequency)) {
    PARSE_ERROR(context, "Physical trace forcing parameters must be finite.");
  }
  if (recorded_weyl_.has_value() and
      (amplitude != 0. or prescribed_quadrupole_.has_value() or
       gauge_sector != detail::SectorImposition::Frozen)) {
    PARSE_ERROR(context,
                "RecordedWeyl requires Amplitude=0, PrescribedQuadrupole=None "
                "and GaugeSector=Frozen (the gauge is supplied by replay).");
  }
  if (prescribed_quadrupole_.has_value() and amplitude != 0.) {
    PARSE_ERROR(
        context,
        "PrescribedQuadrupole requires Amplitude=0 to select one target.");
  }
}

tnsr::aa<DataVector, 3> WorldtubeEvolvedPhysicalGhost::curvature_target(
    const double time, const tnsr::I<DataVector, 3>& coordinates,
    const worldtube::PhysicalScreen& receiver_screen, const Domain<3>* domain,
    const Element<3>* element,
    const domain::FunctionsOfTimeMap* functions_of_time,
    const Mesh<3>* mesh) const {
  if (recorded_weyl_.has_value()) {
    if (domain == nullptr or element == nullptr or
        functions_of_time == nullptr or mesh == nullptr) {
      ERROR(
          "RecordedWeyl requires the receiver domain, element, maps and mesh.");
    }
    return worldtube::incoming_weyl_to_physical_target(
        recorded_weyl_->boundary.recorded_incoming_mode(
            time, coordinates, *domain, *element, *functions_of_time, *mesh),
        receiver_screen);
  }
  if (prescribed_quadrupole_.has_value()) {
    return worldtube::incoming_weyl_to_physical_target(
        prescribed_quadrupole_->incoming_mode(time, coordinates),
        receiver_screen);
  }
  tnsr::aa<DataVector, 3> target(get<0>(coordinates).size(), 0.);
  get<1, 1>(target) = amplitude_ * sin(frequency_ * time);
  get<2, 2>(target) = -get<1, 1>(target);
  return worldtube::physical_projection(target, receiver_screen);
}
WorldtubeEvolvedPhysicalGhost::WorldtubeEvolvedPhysicalGhost(
    CkMigrateMessage* msg)
    : BoundaryCondition<3>(msg) {}
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubeEvolvedPhysicalGhost::get_clone() const {
  return std::make_unique<WorldtubeEvolvedPhysicalGhost>(*this);
}
void WorldtubeEvolvedPhysicalGhost::pup(PUP::er& p) {
  BoundaryCondition<3>::pup(p);
  p | nonphysical_;
  p | amplitude_;
  p | frequency_;
  p | prescribed_quadrupole_;
  p | recorded_weyl_;
}
PUP::able::PUP_ID WorldtubeEvolvedPhysicalGhost::my_PUP_ID = 0;
std::optional<std::string> WorldtubeEvolvedPhysicalGhost::dg_ghost(
    gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> ghost_phi,
    gsl::not_null<Scalar<DataVector>*> ghost_gamma1,
    gsl::not_null<Scalar<DataVector>*> ghost_gamma2,
    gsl::not_null<Scalar<DataVector>*> ghost_lapse,
    gsl::not_null<tnsr::I<DataVector, 3>*> ghost_shift,
    gsl::not_null<tnsr::II<DataVector, 3>*> ghost_inverse_spatial_metric,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
    // c.f. dg_interior_evolved_variables_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    // c.f. dg_interior_temporary_tags
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, 3, Frame::Inertial>& spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, 3, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, 3, Frame::Inertial>&
        spacetime_deriv_gauge_source,
    // c.f. dg_interior_dt_vars_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& logical_dt_phi,
    // c.f. dg_interior_deriv_vars_tags
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi,
    // c.f. dg_gridless_tags
    double time, const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const worldtube::KretschmannFaceData<3>& face_data,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2,
    const DataVector& physical_trace, const Mesh<3>& mesh,
    const ::dg::Formulation formulation,
    const evolution::BoundaryCorrection& boundary_correction,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& volume_inverse_jacobian) const {
  (void)three_index_constraint;
  (void)gauge_source;
  (void)spacetime_deriv_gauge_source;
  (void)logical_dt_spacetime_metric;
  (void)logical_dt_pi;
  (void)logical_dt_phi;
  (void)d_spacetime_metric;
  (void)d_pi;
  (void)d_phi;
  (void)time;
  (void)domain;
  (void)element;
  (void)functions_of_time;
  (void)face_data;
  (void)damping_gamma2;

  if (min(get(gamma1)) != -1. or max(get(gamma1)) != -1.) {
    return "Evolved physical trace requires gamma1=-1.";
  }
  const auto geometry =
      uses_recorded_weyl()
          ? worldtube::mapped_sphere_trace_geometry(
                spacetime_metric,
                tnsr::aa<DataVector, 3>(get(lapse).size(), 0.), coords,
                face_mesh_velocity.value_or(
                    tnsr::I<DataVector, 3>(get(lapse).size(), 0.)),
                mesh.extents(1) - 1)
          : worldtube::spherical_trace_geometry(
                spacetime_metric,
                tnsr::aa<DataVector, 3>(get(lapse).size(), 0.), coords);
  (void)formulation;
  (void)boundary_correction;
  (void)volume_inverse_jacobian;
  DataVector speed = -get(lapse);
  for (size_t i = 0; i < 3; ++i) {
    if (max(abs(normal_vector.get(i) - geometry.normal_vector.get(i))) >
        1.e-10) {
      return "Evolved physical trace requires the lower radial face normal "
             "to agree with the mapped surface geometry.";
    }
    speed -= shift.get(i) * normal_covector.get(i);
    if (face_mesh_velocity.has_value()) {
      speed -= face_mesh_velocity->get(i) * normal_covector.get(i);
    }
  }
  if (max(speed) >= 0. or min(speed + 2. * get(lapse)) <= 0.) {
    return "Evolved physical trace requires lambda_minus<0<lambda_plus.";
  }
  const auto target =
      worldtube::physical_trace_tensor(physical_trace, geometry);
  const worldtube::PhysicalScreen screen(
      spacetime_metric, inverse_spacetime_metric, spacetime_unit_normal_vector,
      normal_vector);
  worldtube::physical_ghost(ghost_pi, ghost_phi, pi, phi, target,
                            normal_covector, normal_vector, screen);
  *ghost_metric = spacetime_metric;
  *ghost_gamma1 = gamma1;
  *ghost_gamma2 = gamma2;
  *ghost_lapse = lapse;
  *ghost_shift = shift;
  *ghost_inverse_spatial_metric = geometry.inverse_spatial_metric;
  return {};
}
std::optional<std::string> WorldtubeEvolvedPhysicalGhost::dg_time_derivative(
    gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
        dt_spacetime_metric_correction,
    gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> dt_pi_correction,
    gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> dt_phi_correction,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
    // c.f. dg_interior_evolved_variables_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    // c.f. dg_interior_temporary_tags
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, 3, Frame::Inertial>& spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, 3, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, 3, Frame::Inertial>&
        spacetime_deriv_gauge_source,
    // c.f. dg_interior_dt_vars_tags
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& logical_dt_phi,
    // c.f. dg_interior_deriv_vars_tags
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi,
    // c.f. dg_gridless_tags
    double time, const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const worldtube::KretschmannFaceData<3>& face_data,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2,
    const DataVector& physical_trace, const Mesh<3>& mesh,
    const ::dg::Formulation formulation,
    const evolution::BoundaryCorrection& boundary_correction,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& volume_inverse_jacobian) const {
  (void)physical_trace;
  if (recorded_weyl_.has_value()) {
    return recorded_weyl_->boundary.dg_time_derivative(
        dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
        face_mesh_velocity, normal_covector, normal_vector, spacetime_metric,
        pi, phi, coords, gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
        spacetime_unit_normal_vector, three_index_constraint, gauge_source,
        spacetime_deriv_gauge_source, logical_dt_spacetime_metric,
        logical_dt_pi, logical_dt_phi, d_spacetime_metric, d_pi, d_phi, time,
        domain, element, functions_of_time, mesh, formulation,
        boundary_correction, volume_inverse_jacobian, face_data, true);
  }
  return nonphysical_->dg_time_derivative(
      dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
      face_mesh_velocity, normal_covector, normal_vector, spacetime_metric, pi,
      phi, coords, gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
      spacetime_unit_normal_vector, three_index_constraint, gauge_source,
      spacetime_deriv_gauge_source, logical_dt_spacetime_metric, logical_dt_pi,
      logical_dt_phi, d_spacetime_metric, d_pi, d_phi, time, domain, element,
      functions_of_time, face_data, damping_gamma2);
}
}  // namespace gh::BoundaryConditions
