// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <optional>
#include <utility>

#include "Evolution/BoundaryCorrectionTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/MomentFeedback.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplay.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags/Formulation.hpp"

namespace gh::BoundaryConditions {
namespace detail {
/// Additive compact pulse in the model's STF component basis. The last
/// diagonal entry is minus the sum of the first two diagonal entries.
struct RecordedMomentPulse {
  struct StartTime {
    using type = double;
    static constexpr Options::String help = "Pulse start in simulation time.";
  };
  struct Duration {
    using type = double;
    static constexpr Options::String help = "Positive pulse duration.";
  };
  struct Electric {
    using type = std::array<double, 5>;
    static constexpr Options::String help =
        "Peak additive [E11,E22,E12,E13,E23], in inverse length squared.";
  };
  struct Magnetic {
    using type = std::array<double, 5>;
    static constexpr Options::String help =
        "Peak additive [B11,B22,B12,B13,B23], in inverse length squared.";
  };
  using options = tmpl::list<StartTime, Duration, Electric, Magnetic>;
  static constexpr Options::String help =
      "Prescribed sin^4 pulse, exactly zero outside its time interval.";
  RecordedMomentPulse() = default;
  RecordedMomentPulse(double start, double duration,
                      std::array<double, 5> electric,
                      std::array<double, 5> magnetic,
                      const Options::Context& context = {});
  gr::np::TidalMoments operator()(double time) const;
  void pup(PUP::er& p);
  double start{0.};
  double duration{1.};
  std::array<double, 5> electric{};
  std::array<double, 5> magnetic{};
};
/// Diagnostic strength of the joint dPi/dPhi input to the moment fit.
struct DerivativeFitReplay {
  struct LiveStrength {
    using type = double;
    static constexpr Options::String help = "Final live strength in [0,1].";
  };
  struct RampStart {
    using type = double;
    static constexpr Options::String help = "Start time of the smooth switch.";
  };
  struct RampDuration {
    using type = double;
    static constexpr Options::String help = "Positive duration of the switch.";
  };
  using options = tmpl::list<LiveStrength, RampStart, RampDuration>;
  static constexpr Options::String help =
      "Joint recorded derivative substitution in the geometric moment fit. "
      "Strength starts at one and ramps to LiveStrength with a quintic step.";
  DerivativeFitReplay() = default;
  DerivativeFitReplay(double strength, double start, double duration,
                      const Options::Context& context = {});
  double operator()(double time) const;
  void pup(PUP::er& p);
  double live_strength{1.};
  double ramp_start{0.};
  double ramp_duration{1.};
};
}  // namespace detail
/// Replay the native GH numerical interface, replacing only the incoming
/// physical u-minus RHS projection by a live curvature-model Bjorhus
/// correction. The historical name is retained for existing TypeD inputs.
/// Optionally replace the incoming constraint sectors by live Bjorhus terms.
/// Optionally use inner-boundary Sommerfeld gauge terms. By default gauge
/// retains replay; u-plus always retains the reference interface contribution.
/// With ReplayBoundaryData: false, use identical interior/ghost fields and
/// only live incoming corrections. This requires constraint preservation and
/// a model gauge, and rejects incoming u-plus characteristics.
/// Supports moving/refined excision faces with
/// Gauss-Lobatto points and AveragedUpwindPenalty.
class WorldtubeReplayTypeD final : public BoundaryCondition<3> {
 public:
  using PhysicalModel = WorldtubeTypeD<3>::PhysicalModel;
  using Mass = WorldtubeTypeD<3>::Mass;
  using MomentRelaxationTime = WorldtubeTypeD<3>::MomentRelaxationTime;
  struct ConstraintPreserving {
    using type = bool;
    static constexpr Options::String help =
        "Use live Bjorhus constraint corrections for incoming u-psi, u-zero "
        "and the constraint projection of u-minus. Otherwise retain replay.";
  };
  struct SommerfeldGauge {
    using type = bool;
    static constexpr Options::String help =
        "Replace the incoming u-minus gauge projection by the live "
        "SommerfeldAbsorbing Bjorhus terms about the excision center. "
        "Otherwise retain replay.";
  };
  struct QhatRadialResponse {
    using type = Options::Auto<std::array<double, 3>, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Dimensionless [time, radial, tangential] responses R*K for "
        "qhat=ell.(Pi-n.Phi), or None. Captures the receiver's initial "
        "profile and radial derivative, transporting spatial components "
        "with Rotation. Fixed-radius spherical-harmonic faces and fixed "
        "resolution only. Mutually exclusive with SommerfeldGauge.";
  };
  struct ReplayBoundaryData {
    using type = bool;
    static constexpr Options::String help =
        "Use recorded ghost fields (default true). If false, no boundary "
        "recording is read and H5Filename/SubfileName may be omitted. "
        "Requires ConstraintPreserving and a model gauge; u-plus must "
        "remain outgoing or zero-speed.";
  };
  struct PhysicalTarget {
    using type = std::string;
    static constexpr Options::String help =
        "Live (default), RecordedModel (fit the model entirely from the "
        "recording), RecordedMoments (recorded tidal fit with live frame and "
        "boost), or RecordedWeyl (recorded true incoming Weyl tensor). "
        "Recorded targets require ReplayBoundaryData and no moment relaxation. "
        "CP and the physical Bjorhus application still use live fields.";
  };
  struct MomentPulse {
    using type =
        Options::Auto<detail::RecordedMomentPulse, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Add a prescribed pulse to RecordedMoments or DerivativeFitReplay.";
  };
  struct MomentDiagnosticsInterval {
    using type = Options::Auto<double, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Minimum simulation-time spacing of full-step MOMENT_RESPONSE log "
        "records, or None. Records donor, imposed, and independently fitted "
        "live moments. Requires RecordedMoments or DerivativeFitReplay.";
  };
  struct MomentFeedback {
    using type =
        Options::Auto<detail::MomentFeedback, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Optional reference-centered feedback diagnostic for RecordedMoments "
        "with legacy Quadrupole. None preserves the existing prescription.";
  };
  struct DerivativeFitReplay {
    using type =
        Options::Auto<detail::DerivativeFitReplay, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "Fit-only joint derivative replay diagnostic, or None (default). "
        "Requires Live QuadrupoleGeometric, replay and no relaxation.";
  };
  using options = tmpl::append<
      WorldtubeReplay::options,
      tmpl::list<PhysicalModel, Mass, MomentRelaxationTime,
                 ConstraintPreserving, SommerfeldGauge, QhatRadialResponse,
                 ReplayBoundaryData, PhysicalTarget, MomentPulse,
                 MomentDiagnosticsInterval, MomentFeedback,
                 DerivativeFitReplay>>;
  static WorldtubeReplayTypeD from_yaml(const Options::Option& options);
  static constexpr Options::String help =
      "Sphere replay with the physical incoming sector supplied by the "
      "selected Bjorhus model and optional live constraint preservation. "
      "Optional Sommerfeld or damping-subtracted radial gauge; "
      "otherwise gauge retains replay. "
      "ReplayBoundaryData: false disables the recorded interface and "
      "requires models for all incoming sectors. Otherwise u-plus retains "
      "replay. Omitting the "
      "model options selects zeroth-order TypeD.";
  WorldtubeReplayTypeD() = default;
  WorldtubeReplayTypeD(std::string filename, std::string subfile,
                       double radius);
  WorldtubeReplayTypeD(
      std::string filename, std::string subfile, double radius,
      PhysicalModel::type physical_model, std::optional<double> mass,
      detail::MomentRelaxation moment_relaxation_time,
      bool constraint_preserving = false, bool sommerfeld_gauge = false,
      std::optional<std::array<double, 3>> qhat_radial_response = std::nullopt,
      bool replay_boundary_data = true, std::string physical_target = "Live",
      std::optional<detail::RecordedMomentPulse> moment_pulse = std::nullopt,
      std::optional<double> moment_diagnostics_interval = std::nullopt,
      std::optional<detail::MomentFeedback> moment_feedback = std::nullopt,
      std::optional<detail::DerivativeFitReplay> derivative_fit_replay =
          std::nullopt,
      const Options::Context& context = {});
  explicit WorldtubeReplayTypeD(CkMigrateMessage* msg);
  worldtube::PhysicalModel physical_model() const {
    return model_.physical_model();
  }
  const std::optional<double>& mass() const { return model_.mass(); }
  bool constraint_preserving() const {
    return model_.constraint_v_psi() == detail::SectorImposition::Bjorhus;
  }
  bool sommerfeld_gauge() const {
    return model_.gauge_sector() ==
           detail::SectorImposition::SommerfeldAbsorbing;
  }
  const std::optional<std::array<double, 3>>& qhat_radial_response() const {
    return qhat_radial_response_;
  }
  bool replay_boundary_data() const { return replay_boundary_data_; }
  const std::optional<double>& moment_relaxation_time() const {
    return model_.moment_relaxation_time();
  }
  const std::optional<gr::np::GeometricFixedPointOptions>& psi0_fixed_point()
      const {
    return model_.psi0_fixed_point();
  }
  bool moment_relaxation_uses_model_time() const {
    return model_.moment_relaxation_uses_model_time();
  }
  const std::string& physical_target() const { return physical_target_; }
  const std::optional<double>& moment_diagnostics_interval() const {
    return moment_diagnostics_interval_;
  }
  gr::np::TidalMoments moment_perturbation(double time) const;
  const std::optional<detail::DerivativeFitReplay>& derivative_fit_replay()
      const {
    return derivative_fit_replay_;
  }
  /// Same path used by stage boundary application and full-step diagnostics.
  worldtube::MatchingEvaluation derivative_fit_evaluation(
      double time, const tnsr::aa<DataVector, 3>& g,
      const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
      const tnsr::iaa<DataVector, 3>& d_pi,
      const tnsr::ijaa<DataVector, 3>& d_phi,
      const tnsr::i<DataVector, 3>& normal,
      const tnsr::I<DataVector, 3>& coords, const Domain<3>& domain,
      const Element<3>& element,
      const domain::FunctionsOfTimeMap& functions_of_time,
      const worldtube::KretschmannFaceData<3>& face_data) const;
  /// Called on full steps by the face-data action. Uses volume derivatives
  /// before restriction to the face, and never updates the prescribed target.
  void observe_moment_response(
      double time, const tnsr::aa<DataVector, 3>& spacetime_metric,
      const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
      const tnsr::I<DataVector, 3>& coordinates, const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const Domain<3>& domain, const Element<3>& element,
      const domain::FunctionsOfTimeMap& functions_of_time,
      const worldtube::KretschmannFaceData<3>& face_data) const;
  /// Incoming Weyl target constructed without any receiver field or history.
  /// Coordinates and maps identify the recorded sphere; the mesh supplies
  /// only its quadrature weights.
  tnsr::ii<DataVector, 3> recorded_incoming_mode(
      double time, const tnsr::I<DataVector, 3>& coords,
      const Domain<3>& domain, const Element<3>& element,
      const domain::FunctionsOfTimeMap& functions_of_time,
      const Mesh<3>& mesh) const;
  /// Fit the five complex tidal moments entirely from the recording.
  /// Requires an order-two model and a full spherical-harmonic face.
  gr::np::TidalMoments recorded_tidal_moments(
      double time, const tnsr::I<DataVector, 3>& coords,
      const Domain<3>& domain, const Element<3>& element,
      const domain::FunctionsOfTimeMap& functions_of_time,
      const Mesh<3>& mesh) const;
  WRAPPED_PUPable_decl_template(WorldtubeReplayTypeD);
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
  using dg_gridless_tags = tmpl::list<
      ::Tags::Time, domain::Tags::Domain<3>, domain::Tags::Element<3>,
      domain::Tags::FunctionsOfTime, domain::Tags::Mesh<3>,
      ::dg::Tags::Formulation, evolution::Tags::BoundaryCorrection,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>,
      worldtube::Tags::KretschmannFaceData<3>>;
  std::optional<std::string> dg_ghost(
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
      const worldtube::KretschmannFaceData<3>& face_data = {}) const;
  std::optional<std::string> dg_time_derivative(
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
      const worldtube::KretschmannFaceData<3>& face_data = {},
      bool nonphysical_only = false) const;

 private:
  std::pair<tnsr::ii<DataVector, 3>, std::optional<gr::np::TidalMoments>>
  recorded_target_data(double time, const tnsr::I<DataVector, 3>& coords,
                       const Domain<3>& domain, const Element<3>& element,
                       const domain::FunctionsOfTimeMap& functions_of_time,
                       const Mesh<3>& mesh) const;
  std::optional<WorldtubeReplay> replay_{};
  WorldtubeTypeD<3> model_{};
  std::optional<std::array<double, 3>> qhat_radial_response_{};
  bool replay_boundary_data_{true};
  std::string physical_target_{"Live"};
  std::optional<detail::RecordedMomentPulse> moment_pulse_{};
  std::optional<double> moment_diagnostics_interval_{};
  std::optional<detail::MomentFeedback> moment_feedback_{};
  std::optional<detail::DerivativeFitReplay> derivative_fit_replay_{};
};
}  // namespace gh::BoundaryConditions

template <>
struct Options::create_from_yaml<gh::BoundaryConditions::WorldtubeReplayTypeD> {
  template <typename Metavariables>
  static gh::BoundaryConditions::WorldtubeReplayTypeD create(
      const Options::Option& options) {
    return gh::BoundaryConditions::WorldtubeReplayTypeD::from_yaml(options);
  }
};
