// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubePhysicalGhost.hpp"

#include <cmath>
#include <utility>

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryTrace.hpp"
#include "Options/ParseError.hpp"

namespace gh::BoundaryConditions {
namespace {
std::optional<std::string> check_characteristics(
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift, const tnsr::i<DataVector, 3>& normal,
    const std::optional<tnsr::I<DataVector, 3>>& velocity) {
  for (size_t p = 0; p < get(lapse).size(); ++p) {
    if (get(gamma1)[p] != -1.) {
      return "WorldtubePhysicalGhost requires gamma1=-1.";
    }
    double advection = 0.;
    for (size_t i = 0; i < 3; ++i) {
      advection += normal.get(i)[p] *
                   (shift.get(i)[p] + (velocity ? velocity->get(i)[p] : 0.));
    }
    if (not(-get(lapse)[p] - advection < 0. and
            get(lapse)[p] - advection > 0.)) {
      return "WorldtubePhysicalGhost requires lambda_minus<0<lambda_plus.";
    }
  }
  return {};
}
}  // namespace
WorldtubePhysicalGhost::WorldtubePhysicalGhost(
    std::unique_ptr<evolution::initial_data::InitialData> analytic_prescription,
    const detail::SectorImposition gauge_sector,
    const Options::Context& context)
    : analytic_(std::move(analytic_prescription)),
      nonphysical_(detail::SectorImposition::Bjorhus,
                   detail::SectorImposition::Frozen, gauge_sector,
                   detail::PhysicalModel::None, std::nullopt, std::nullopt,
                   context) {
  if (gauge_sector != detail::SectorImposition::SommerfeldAbsorbing and
      gauge_sector != detail::SectorImposition::SommerfeldOutgoing and
      gauge_sector != detail::SectorImposition::Frozen) {
    PARSE_ERROR(context,
                "WorldtubePhysicalGhost supports SommerfeldAbsorbing, "
                "SommerfeldOutgoing, or Frozen gauge sectors.");
  }
}
WorldtubePhysicalGhost::WorldtubePhysicalGhost(CkMigrateMessage* msg)
    : BoundaryCondition<3>(msg) {}
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubePhysicalGhost::get_clone() const {
  return std::make_unique<WorldtubePhysicalGhost>(*this);
}
void WorldtubePhysicalGhost::pup(PUP::er& p) {
  BoundaryCondition<3>::pup(p);
  p | analytic_;
  p | nonphysical_;
}
PUP::able::PUP_ID WorldtubePhysicalGhost::my_PUP_ID = 0;
std::optional<std::string> WorldtubePhysicalGhost::dg_ghost(
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
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2)
    const {
  (void)three_index_constraint;
  (void)gauge_source;
  (void)spacetime_deriv_gauge_source;
  (void)logical_dt_spacetime_metric;
  (void)logical_dt_pi;
  (void)logical_dt_phi;
  (void)d_spacetime_metric;
  (void)d_pi;
  (void)d_phi;
  (void)domain;
  (void)element;
  (void)functions_of_time;
  (void)face_data;
  (void)damping_gamma2;
  if (const auto error = check_characteristics(
          gamma1, lapse, shift, normal_covector, face_mesh_velocity)) {
    return error;
  }
  // Obtain the prescribed V=Pi-s.Phi using the receiver's normal. A metric
  // mismatch must not inject an artificial -gamma2*g_analytic trace.
  if (const auto error = analytic_.dg_ghost(
          ghost_metric, ghost_pi, ghost_phi, ghost_gamma1, ghost_gamma2,
          ghost_lapse, ghost_shift, ghost_inverse_spatial_metric,
          face_mesh_velocity, normal_covector, normal_vector, coords, gamma1,
          gamma2, time)) {
    return error;
  }
  auto target = *ghost_pi;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t i = 0; i < 3; ++i) {
        target.get(a, b) -= normal_vector.get(i) * ghost_phi->get(i, a, b);
      }
    }
  }
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
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      ghost_inverse_spatial_metric->get(i, j) =
          inverse_spacetime_metric.get(i + 1, j + 1) +
          spacetime_unit_normal_vector.get(i + 1) *
              spacetime_unit_normal_vector.get(j + 1);
    }
  }
  return {};
}
std::optional<std::string> WorldtubePhysicalGhost::dg_time_derivative(
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
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2)
    const {
  if (const auto error = check_characteristics(
          gamma1, lapse, shift, normal_covector, face_mesh_velocity)) {
    return error;
  }
  if (const auto error = nonphysical_.dg_time_derivative(
          dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
          face_mesh_velocity, normal_covector, normal_vector, spacetime_metric,
          pi, phi, coords, gamma1, gamma2, lapse, shift,
          inverse_spacetime_metric, spacetime_unit_normal_vector,
          three_index_constraint, gauge_source, spacetime_deriv_gauge_source,
          logical_dt_spacetime_metric, logical_dt_pi, logical_dt_phi,
          d_spacetime_metric, d_pi, d_phi, time, domain, element,
          functions_of_time, face_data, damping_gamma2)) {
    return error;
  }
  // Frozen in WorldtubeTypeD means cancelling the physical volume RHS.
  // Remove that correction here: the ghost flux alone supplies this sector.
  auto minus_correction = *dt_pi_correction;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      minus_correction.get(a, b) -=
          get(gamma2) * dt_spacetime_metric_correction->get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        minus_correction.get(a, b) -=
            normal_vector.get(i) * dt_phi_correction->get(i, a, b);
      }
    }
  }
  const worldtube::PhysicalScreen screen(
      spacetime_metric, inverse_spacetime_metric, spacetime_unit_normal_vector,
      normal_vector);
  const auto physical =
      worldtube::physical_projection(minus_correction, screen);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_pi_correction->get(a, b) -= 0.5 * physical.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        dt_phi_correction->get(i, a, b) +=
            0.5 * normal_covector.get(i) * physical.get(a, b);
      }
    }
  }
  return {};
}
}  // namespace gh::BoundaryConditions
