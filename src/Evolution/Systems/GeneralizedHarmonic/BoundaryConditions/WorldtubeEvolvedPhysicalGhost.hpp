// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubePhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplayTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryEvolution.hpp"

namespace gh::BoundaryConditions {
namespace detail {
/// The existing Quadrupole model with prescribed moments and a known, static
/// Schwarzschild Kerr-Schild frame. Coefficients use the model's Cholesky
/// triad, just as RecordedMoments. No receiver field or curvature fit enters
/// U^{8-}.
struct PrescribedSchwarzschildQuadrupole {
  struct Mass {
    using type = double;
    static constexpr Options::String help = "Positive Schwarzschild mass.";
  };
  struct MomentPulse {
    using type = RecordedMomentPulse;
    static constexpr Options::String help =
        "Prescribed electric and magnetic moments in the legacy model axes.";
  };
  using options = tmpl::list<Mass, MomentPulse>;
  static constexpr Options::String help =
      "Existing NP Quadrupole target on a fixed Schwarzschild frame. "
      "Requires a centered sphere outside r=2M. The target is prescribed; "
      "only its screen projection uses the receiver metric.";
  PrescribedSchwarzschildQuadrupole() = default;
  PrescribedSchwarzschildQuadrupole(double mass, RecordedMomentPulse pulse,
                                    const Options::Context& context = {});
  tnsr::ii<DataVector, 3> incoming_mode(
      double time, const tnsr::I<DataVector, 3>& coordinates) const;
  void pup(PUP::er& p);
  double mass{1.};
  RecordedMomentPulse pulse{};
};
/// Recorded true incoming Weyl curvature, with recorded gauge and live CP.
struct RecordedWeylTarget {
  using H5Filename = WorldtubeReplay::H5Filename;
  using SubfileName = WorldtubeReplay::SubfileName;
  using Radius = WorldtubeReplay::Radius;
  using options = tmpl::list<H5Filename, SubfileName, Radius>;
  static constexpr Options::String help =
      "Recorded Weyl curvature on a moving sphere.";
  RecordedWeylTarget() = default;
  RecordedWeylTarget(std::string filename, std::string subfile, double radius);
  void pup(PUP::er& p);
  WorldtubeReplayTypeD boundary;
};
}  // namespace detail

/// Fixed-sphere prototype. The two physical amplitudes are evolved in the
/// element DataBox with the bulk time stepper. Prescribed W forcing or the
/// existing NP quadrupole on a prescribed Schwarzschild frame is supported.
/// No live NP fit is used.
class WorldtubeEvolvedPhysicalGhost final : public BoundaryCondition<3> {
 public:
  using GaugeSector = WorldtubePhysicalGhost::GaugeSector;
  struct Amplitude {
    using type = double;
    static constexpr Options::String help =
        "Amplitude of W_target=A sin(omega*t) Q[diag(1,-1,0)].";
  };
  struct Frequency {
    using type = double;
    static constexpr Options::String help =
        "Angular frequency of the W target.";
  };
  struct PrescribedQuadrupole {
    using type = Options::Auto<detail::PrescribedSchwarzschildQuadrupole,
                               Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Prescribed NP quadrupole pulse, or None (default). Requires "
        "Amplitude=0; otherwise the oscillating W target is used.";
  };
  struct RecordedWeyl {
    using type =
        Options::Auto<detail::RecordedWeylTarget, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "True recorded incoming Weyl target. Retains recorded gauge and live "
        "CP; "
        "requires Amplitude=0, GaugeSector=Frozen and "
        "PrescribedQuadrupole=None.";
  };
  using options = tmpl::list<GaugeSector, Amplitude, Frequency,
                             PrescribedQuadrupole, RecordedWeyl>;
  static WorldtubeEvolvedPhysicalGhost from_yaml(const Options::Option& input);
  static constexpr Options::String help =
      "Evolved physical trace with prescribed curvature on a full-sphere face. "
      "Requires global stepping, no surface AMR and gamma1=-1. Supports "
      "Harmonic, AnalyticChristoffel and DampedHarmonic. RecordedWeyl also "
      "supports moving maps and retains recorded gauge plus live CP.";
  WorldtubeEvolvedPhysicalGhost() = default;
  WorldtubeEvolvedPhysicalGhost(
      detail::SectorImposition gauge_sector, double amplitude, double frequency,
      std::optional<detail::PrescribedSchwarzschildQuadrupole>
          prescribed_quadrupole = std::nullopt,
      std::optional<detail::RecordedWeylTarget> recorded_weyl = std::nullopt,
      const Options::Context& context = {});
  explicit WorldtubeEvolvedPhysicalGhost(CkMigrateMessage* msg);
  WRAPPED_PUPable_decl_template(WorldtubeEvolvedPhysicalGhost);
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override;
  void pup(PUP::er& p) override;
  double amplitude() const { return amplitude_; }
  double frequency() const { return frequency_; }
  bool uses_recorded_weyl() const { return recorded_weyl_.has_value(); }
  tnsr::aa<DataVector, 3> curvature_target(
      double time, const tnsr::I<DataVector, 3>& coordinates,
      const worldtube::PhysicalScreen& receiver_screen,
      const Domain<3>* domain = nullptr, const Element<3>* element = nullptr,
      const domain::FunctionsOfTimeMap* functions_of_time = nullptr,
      const Mesh<3>* mesh = nullptr) const;
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::GhostAndTimeDerivative;
  using dg_interior_evolved_variables_tags =
      WorldtubeTypeD<3>::dg_interior_evolved_variables_tags;
  using dg_interior_temporary_tags =
      WorldtubeTypeD<3>::dg_interior_temporary_tags;
  using dg_interior_dt_vars_tags = WorldtubeTypeD<3>::dg_interior_dt_vars_tags;
  using dg_interior_deriv_vars_tags =
      WorldtubeTypeD<3>::dg_interior_deriv_vars_tags;
  using dg_gridless_tags = tmpl::push_back<
      WorldtubeTypeD<3>::dg_gridless_tags, worldtube::Tags::PhysicalTrace,
      domain::Tags::Mesh<3>, ::dg::Tags::Formulation,
      evolution::Tags::BoundaryCorrection,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>>;

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
      const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2,
      const DataVector& physical_trace, const Mesh<3>& mesh,
      ::dg::Formulation formulation,
      const evolution::BoundaryCorrection& boundary_correction,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& volume_inverse_jacobian) const;

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
      const ConstraintDamping::DampingFunction<3, Frame::Grid>& damping_gamma2,
      const DataVector& physical_trace, const Mesh<3>& mesh,
      ::dg::Formulation formulation,
      const evolution::BoundaryCorrection& boundary_correction,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& volume_inverse_jacobian) const;

 private:
  std::optional<WorldtubePhysicalGhost> nonphysical_;
  double amplitude_ = 0.;
  double frequency_ = 0.;
  std::optional<detail::PrescribedSchwarzschildQuadrupole>
      prescribed_quadrupole_{};
  std::optional<detail::RecordedWeylTarget> recorded_weyl_{};
};
}  // namespace gh::BoundaryConditions

template <>
struct Options::create_from_yaml<
    gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost> {
  template <typename Metavariables>
  static auto create(const Options::Option& input) {
    return gh::BoundaryConditions::WorldtubeEvolvedPhysicalGhost::from_yaml(
        input);
  }
};
