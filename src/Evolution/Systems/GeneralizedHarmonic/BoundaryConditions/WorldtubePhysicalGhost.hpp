// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletAnalytic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"

namespace gh::BoundaryConditions {
/// Prescribed physical trace coupled through the upwind flux, with live
/// constraint-preserving and gauge time-derivative corrections.
///
/// This first-stage boundary uses an analytic prescription for the screen
/// projection of Pi-s.Phi. The projection and s are evaluated with receiver
/// geometry. It supplies only the two physical components of incoming u-minus;
/// every other ghost characteristic is copied from the receiver. The physical
/// time-derivative correction is zero. No independent surface state or NP
/// curvature-to-trace evolution is enabled by this boundary.
///
/// Requires gamma1=-1 and a timelike tube (lambda_minus<0<lambda_plus).
class WorldtubePhysicalGhost final : public BoundaryCondition<3> {
 public:
  using AnalyticPrescription = DirichletAnalytic<3>::AnalyticPrescription;
  struct GaugeSector {
    using type = detail::SectorImposition;
    static constexpr Options::String help =
        "SommerfeldAbsorbing, SommerfeldOutgoing, or Frozen gauge sector.";
  };
  using options = tmpl::list<AnalyticPrescription, GaugeSector>;
  static constexpr Options::String help =
      "Analytic physical incoming trace through the ghost flux with live "
      "constraint preservation. Prototype for gamma1=-1 timelike tubes.";
  WorldtubePhysicalGhost() = default;
  WorldtubePhysicalGhost(std::unique_ptr<evolution::initial_data::InitialData>
                             analytic_prescription,
                         detail::SectorImposition gauge_sector,
                         const Options::Context& context = {});
  explicit WorldtubePhysicalGhost(CkMigrateMessage* msg);
  WRAPPED_PUPable_decl_template(WorldtubePhysicalGhost);
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override;
  void pup(PUP::er& p) override;
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::GhostAndTimeDerivative;
  using dg_interior_evolved_variables_tags =
      WorldtubeTypeD<3>::dg_interior_evolved_variables_tags;
  using dg_interior_temporary_tags =
      WorldtubeTypeD<3>::dg_interior_temporary_tags;
  using dg_interior_dt_vars_tags = WorldtubeTypeD<3>::dg_interior_dt_vars_tags;
  using dg_interior_deriv_vars_tags =
      WorldtubeTypeD<3>::dg_interior_deriv_vars_tags;
  using dg_gridless_tags = WorldtubeTypeD<3>::dg_gridless_tags;

  std::optional<std::string> dg_ghost(
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
      const tnsr::A<DataVector, 3, Frame::Inertial>&
          spacetime_unit_normal_vector,
      const tnsr::iaa<DataVector, 3, Frame::Inertial>& three_index_constraint,
      const tnsr::a<DataVector, 3, Frame::Inertial>& gauge_source,
      const tnsr::ab<DataVector, 3, Frame::Inertial>&
          spacetime_deriv_gauge_source,
      // c.f. dg_interior_dt_vars_tags
      const tnsr::aa<DataVector, 3, Frame::Inertial>&
          logical_dt_spacetime_metric,
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
      const;

  std::optional<std::string> dg_time_derivative(
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*>
          dt_spacetime_metric_correction,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> dt_pi_correction,
      gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*>
          dt_phi_correction,
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
      const tnsr::A<DataVector, 3, Frame::Inertial>&
          spacetime_unit_normal_vector,
      const tnsr::iaa<DataVector, 3, Frame::Inertial>& three_index_constraint,
      const tnsr::a<DataVector, 3, Frame::Inertial>& gauge_source,
      const tnsr::ab<DataVector, 3, Frame::Inertial>&
          spacetime_deriv_gauge_source,
      // c.f. dg_interior_dt_vars_tags
      const tnsr::aa<DataVector, 3, Frame::Inertial>&
          logical_dt_spacetime_metric,
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
      const;

 private:
  DirichletAnalytic<3> analytic_;
  WorldtubeTypeD<3> nonphysical_;
};
}  // namespace gh::BoundaryConditions
