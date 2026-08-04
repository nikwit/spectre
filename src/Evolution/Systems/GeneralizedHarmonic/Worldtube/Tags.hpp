// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Options/Auto.hpp"
#include "Options/String.hpp"
#include "Time/History.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

/// Online first-order worldtube matching for the generalized harmonic
/// system: fit the affine-map parameters from the evolved fields on the
/// excision sphere as the run progresses, to drive the ghost gauge boundary
/// condition (findings §15b/§15c prescription, made self-consistent).
namespace gh::Worldtube {

/// Number of first-order map parameters
/// (qdot0, beta_i, qdot^i, sigma_ij).
static constexpr size_t num_map_parameters = 13;

/// Option-created configuration of the online matcher.
struct MatcherConfig {
  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of the black hole"};
    static type lower_bound() { return 0.; }
  };
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "The worldtube center, fixed in time"};
  };
  struct CenterVelocity {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Kinematic pin for the map velocity: qdot^i = (1 + qdot^0) v_center. "
        "Zero for a hole at rest."};
  };
  struct TraceStrainPin {
    using type = double;
    static constexpr Options::String help = {
        "Pin for tr(sigma)/3. Zero for an unperturbed single hole; for a "
        "tidally perturbed hole use the invariant-radius value (findings "
        "§12)."};
  };
  struct FitLMax {
    using type = size_t;
    static constexpr Options::String help = {
        "The fit residual keeps spherical-harmonic modes with l <= FitLMax "
        "of the gauge components. The first-order model space has content "
        "up to l = 4."};
  };
  struct FitInterval {
    using type = double;
    static constexpr Options::String help = {
        "Minimum simulation-time interval between fits. The fit runs at the "
        "first step whose time exceeds the previous fit time by this "
        "amount. In the default strict first-order value mode the fitted "
        "p_(1) state is held constant while the zeroth-order model center is "
        "advanced with its fitted O(epsilon) velocity; derivative ODE modes "
        "retain their higher-order extrapolation."};
    static type lower_bound() { return 0.; }
  };
  struct FitRadialIndex {
    using type = size_t;
    static constexpr Options::String help = {
        "Radial collocation index of the worldtube element at which the fit "
        "reads the fields: 0 is the excision face itself. In the closed "
        "loop the boundary condition clamps the gauge components of u^- at "
        "index 0 (the fit then sees no signal and the loop settles on the "
        "trivial frozen fixed point), so fit a few points off the boundary "
        "(index 4 is r ~ 2.55 for the standard Shell0) and let the BC "
        "impose the model at the face."};
  };
  struct RateOde {
    using type = bool;
    static constexpr Options::String help = {
        "HIGHER-ORDER EXPERIMENT (not part of the strict Dhesi slow-time "
        "first-order system): obtain p(t) by fitting the parameter rates "
        "pdot from the "
        "time-derivative content of the boundary data (dt g = beta.Phi - "
        "alpha.Pi, projected on the covariant response columns, a linear "
        "solve) and integrating the first-order ODE dp/dt = pdot by the "
        "trapezoid rule from p = 0 (Schwarzschild) at the first fit. The "
        "algebraic value fit is bypassed. Complementary to the default "
        "mode: uses the Pi channel instead of the u^- gauge components, "
        "smooth by construction, but integration accumulates rate bias "
        "with no restoring force."};
  };
  struct SecondOrderOde {
    using type = bool;
    static constexpr Options::String help = {
        "HIGHER-ORDER EXPERIMENT (not part of the strict Dhesi slow-time "
        "first-order system): obtain p(t) by fitting the parameter "
        "accelerations pddot from the "
        "evolution equations (d2t g assembled from the GH right-hand sides "
        "stored in the DataBox, projected on the covariant response columns, "
        "a linear solve) and integrating the second-order ODE by velocity "
        "Verlet from (p, pdot) = (0, 0) (Schwarzschild) at the first fit. "
        "Produces a smooth, dynamically consistent (p, pdot) pair - the "
        "scalar-worldtube architecture. Mutually exclusive with RateOde and "
        "FitCenterOffset."};
  };
  struct StepperOde {
    using type = bool;
    static constexpr Options::String help = {
        "HIGHER-ORDER EXPERIMENT (not part of the strict Dhesi slow-time "
        "first-order system): integrate (p, pdot) through the element's own "
        "TimeStepper: the "
        "acceleration is measured from the just-computed right-hand sides "
        "at every (sub)step and recorded in a TimeSteppers::History, with "
        "the integration order slaved to the system history, dormancy "
        "during self-start, and TimeStepId-keyed rewind on step rejection. "
        "FitInterval is ignored. Mutually exclusive with the other modes."};
  };
  struct GaugeDamping {
    using type = double;
    static constexpr Options::String help = {
        "Restoring rate gamma for the ODE modes: the integrated "
        "acceleration is pddot_meas - 2*gamma*pdot - gamma^2*(p - pin), "
        "with the pin the kinematic pin values (zero for the unpinned "
        "parameters). The affine map is pure gauge, so (p, pdot) is a "
        "double zero root of the raw equations-of-motion loop; any "
        "measurement bias epsilon splits it into +-sqrt(epsilon), one "
        "growing. gamma^2 > epsilon turns the pair into a damped "
        "oscillator, i.e. gauge-fixes the flat direction. Zero disables."};
    static double lower_bound() { return 0.; }
  };
  struct UPlusAnchor {
    using type = double;
    static constexpr Options::String help = {
        "Anchor rate kappa for the StepperOde mode: every FitInterval the "
        "map parameters are value-fitted to the gauge projection of the "
        "OUTGOING characteristic u^+ = Pi - n Phi - gamma2 g (same null "
        "projectors as u^-), the one channel the ghost BC does not set, "
        "and the ODE acceleration gains -2 kappa (pdot - pdot_anchor) "
        "- kappa^2 (p - p_anchor) on the free directions, unfolded "
        "through the pins. u^+ is the data the ambient evolution feeds "
        "the excision, so this anchors the map to the ambient chart and "
        "lifts the double zero root of the self-referential loop "
        "(kappa^2 must exceed the loop bias epsilon). Zero disables."};
    static double lower_bound() { return 0.; }
  };
  struct FitUPlus {
    using type = bool;
    static constexpr Options::String help = {
        "Target the value fit at the gauge projection of the OUTGOING "
        "characteristic u^+ = Pi - n Phi - gamma2 g instead of u^-. u^+ "
        "is the data the ambient evolution feeds the excision and the one "
        "channel the ghost BC does not set, so the fit reads a quantity "
        "the closed loop cannot manufacture (the u^- value fit at the "
        "face measures the penalty-driven combination and is "
        "tautological). Combined with StepperOde, the acceleration "
        "sensor targets dt of the u^+ gauge projection instead of the "
        "all-components d2t g projection."};
  };
  struct KretschmannTracePin {
    using type = bool;
    static constexpr Options::String help = {
        "Set the trace-strain pin from the curvature each fit instead of "
        "the TraceStrainPin constant: the vacuum Gauss-Bonnet scalar "
        "(= Kretschmann) gives the invariant harmonic radius rho_GB = "
        "(48 M^2/K)^{1/6} - M, and tr sigma / 3 = 1 - <rho_GB>/<R> on "
        "the face (findings 12: the fit-independent measurement of the "
        "one strain direction the metric fit cannot own). Applies to the "
        "value-fit modes; the ODE modes keep the constant."};
  };
  struct TracePinInterval {
    using type = double;
    static constexpr Options::String help = {
        "Recompute cadence of the Kretschmann trace pin (held constant "
        "in between). The pin is secular physics and must NOT be an "
        "instantaneous closure: recomputing it every substep couples the "
        "fit to the fields through the second-derivative operator at "
        "unit gain and blows up in ~2 M regardless of the value-fit "
        "cadence (measured). Ignored unless KretschmannTracePin."};
    static double lower_bound() { return 0.; }
  };
  struct FitTraceStrain {
    using type = bool;
    static constexpr Options::String help = {
        "Free the trace of the strain as a tenth fit parameter instead "
        "of pinning it: the trace then rides the same ambient anchor as "
        "the other directions (meaningful for the u^+ target, which the "
        "ghost BC cannot drag). TraceStrainPin/KretschmannTracePin then "
        "affect nothing in the solve; the Kretschmann measurement, if "
        "enabled, is logged as an open-loop diagnostic only (closing it "
        "as feedback measured a sampled-loop gain of -2, findings 15t)."};
  };
  struct SpatialMonopoleWeight {
    using type = double;
    static constexpr Options::String help = {
        "Relative weight of the l = 0 modes of the spatial components in "
        "the rate/acceleration least squares. The q8 block analysis found "
        "the spatial monopole absorbs unmodeled content into the trace "
        "strain and clock rate without degrading the condition number; "
        "set to 0 to exclude it from the solve. The demoted rows are "
        "still evaluated in the held-out closure diagnostics."};
    static double lower_bound() { return 0.; }
  };
  struct UPlusBlockWeights {
    using type = std::array<double, 15>;
    static constexpr Options::String help = {
        "Multiplicative residual weights for the FitUPlus value and "
        "StepperOde acceleration solves, ordered [A_l0..A_l4, "
        "C_l0..C_l4, V_l0..V_l4]. A zero excludes a block from the solve, "
        "but its unweighted residual is still evaluated in the diagnostics. "
        "These weights do not apply to the rate or metric-acceleration "
        "fits."};
  };
  struct FitCenterOffset {
    using type = bool;
    static constexpr Options::String help = {
        "Also fit the zeroth-order spatial offset q^i of the map: the model "
        "is evaluated about Center + q with q fitted (exactly, not "
        "linearized), correcting an imperfect worldtube center online. The "
        "time offset q^0 is an exact zero mode of the static-in-time model "
        "and is never fitted."};
  };

  struct FitVelocity {
    using type = bool;
    static constexpr Options::String help = {
        "Free the three map velocities qdot^i instead of pinning them "
        "kinematically to (1 + qdot^0) CenterVelocity. All thirteen map "
        "parameters are then fitted (with FitTraceStrain), which is what "
        "the binary needs: the hole's centre moves relative to the "
        "excision centre and its velocity is wanted as an input to the "
        "control system rather than an output of it. Only the value fit "
        "supports this; combining it with a rate/ODE mode is an error. "
        "findings 9 measured qdot^i = (1 + qdot^0) dz/dT offline to 2.3% "
        "-- freeing it tests that relation online."};
  };
  struct CentreAdvection {
    using type = bool;
    static constexpr Options::String help = {
        "Include the strict first-order centre-motion term -qdot^k "
        "Phi^(0)_kab in the model's d_t g and advance the zeroth-order "
        "analytic model center between value fits. This does not move the "
        "mesh. It is the only channel that puts qdot^i into Pi; with it off, "
        "qdot^i is visible only in g and Phi. Physically it should always be "
        "on -- off is for the A/B that demonstrates the difference. "
        "Identically zero when qdot^i = 0, so it cannot alter a run whose "
        "velocity is pinned to a vanishing CenterVelocity."};
  };

  using options =
      tmpl::list<Mass, Center, CenterVelocity, TraceStrainPin, FitLMax,
                 FitInterval, FitCenterOffset, RateOde, SecondOrderOde,
                 StepperOde, GaugeDamping, UPlusAnchor, FitUPlus,
                 KretschmannTracePin, TracePinInterval, FitTraceStrain,
                 FitVelocity, CentreAdvection, SpatialMonopoleWeight,
                 UPlusBlockWeights, FitRadialIndex>;
  static constexpr Options::String help = {
      "Online worldtube matching. By default, algebraically fit the "
      "instantaneous center and 13 affine-map coefficients using a strict "
      "Dhesi slow-time first-order model: no coefficient rates, no "
      "between-fit extrapolation, and no nonlinear resummation. RateOde, "
      "SecondOrderOde, and StepperOde select explicitly higher-order "
      "experimental systems instead."};

  MatcherConfig() = default;
  MatcherConfig(double mass, const std::array<double, 3>& center,
                const std::array<double, 3>& center_velocity,
                double trace_strain_pin, size_t fit_l_max, double fit_interval,
                bool fit_center_offset, bool rate_ode, bool second_order_ode,
                bool stepper_ode, double gauge_damping, double uplus_anchor,
                bool fit_uplus, bool kretschmann_trace_pin,
                double trace_pin_interval, bool fit_trace_strain,
                bool fit_velocity, bool centre_advection,
                double spatial_monopole_weight,
                const std::array<double, 15>& uplus_block_weights,
                size_t fit_radial_index)
      : mass(mass),
        center(center),
        center_velocity(center_velocity),
        trace_strain_pin(trace_strain_pin),
        fit_l_max(fit_l_max),
        fit_interval(fit_interval),
        fit_center_offset(fit_center_offset),
        rate_ode(rate_ode),
        second_order_ode(second_order_ode),
        stepper_ode(stepper_ode),
        gauge_damping(gauge_damping),
        uplus_anchor(uplus_anchor),
        fit_uplus(fit_uplus),
        kretschmann_trace_pin(kretschmann_trace_pin),
        trace_pin_interval(trace_pin_interval),
        fit_trace_strain(fit_trace_strain),
        fit_velocity(fit_velocity),
        centre_advection(centre_advection),
        spatial_monopole_weight(spatial_monopole_weight),
        uplus_block_weights(uplus_block_weights),
        fit_radial_index(fit_radial_index) {
    if (fit_velocity and (rate_ode or second_order_ode or stepper_ode)) {
      ERROR(
          "FitVelocity is implemented for the value fit only: the rate and "
          "acceleration solves carry their own velocity pin (qddot^i = "
          "qddot^0 v_centre) and a fixed nine-column layout, so enabling "
          "both would silently keep the velocity pinned. Set RateOde, "
          "SecondOrderOde and StepperOde false.");
    }
    for (const double weight : uplus_block_weights) {
      if (weight < 0.) {
        ERROR("UPlusBlockWeights entries must be non-negative, got " << weight
                                                                     << ".");
      }
    }
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  double mass = 1.0;
  std::array<double, 3> center{{0., 0., 0.}};
  std::array<double, 3> center_velocity{{0., 0., 0.}};
  double trace_strain_pin = 0.;
  size_t fit_l_max = 4;
  double fit_interval = 0.;
  bool fit_center_offset = false;
  bool rate_ode = false;
  bool second_order_ode = false;
  bool stepper_ode = false;
  double gauge_damping = 0.;
  double uplus_anchor = 0.;
  bool fit_uplus = false;
  bool kretschmann_trace_pin = false;
  double trace_pin_interval = 0.5;
  bool fit_trace_strain = false;
  bool fit_velocity = false;
  bool centre_advection = true;
  double spatial_monopole_weight = 1.;
  std::array<double, 15> uplus_block_weights{
      {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}};
  size_t fit_radial_index = 0;
};

bool operator==(const MatcherConfig& lhs, const MatcherConfig& rhs);
bool operator!=(const MatcherConfig& lhs, const MatcherConfig& rhs);

/// The per-element state of the online matcher: the latest fitted
/// parameters, their estimated rates, and the previous fit for the
/// backward-difference drives. Meaningful only on the element owning the
/// excision face; default (invalid) everywhere else.
struct MapParameterData {
  double last_fit_time = std::numeric_limits<double>::lowest();
  double previous_fit_time = std::numeric_limits<double>::lowest();
  std::array<double, num_map_parameters> p{};
  std::array<double, num_map_parameters> p_previous{};
  std::array<double, num_map_parameters> pdot{};
  std::array<double, num_map_parameters> pddot{};
  std::array<double, 3> center_offset{};
  /// Stepper-integrated mode: the 26-component state (p, pdot) and its
  /// time-stepper history
  DataVector ode_state{};
  TimeSteppers::History<DataVector> ode_history{};
  /// Stepper-integrated mode: the state at the start of the current step
  /// and that step's substep-0 id, so a rejected step restarts from the
  /// correct value instead of the rejected end-of-step value
  DataVector ode_step_start{};
  TimeStepId ode_step_id{};
  /// Kretschmann trace pin: the held value and its measurement time
  double trace_pin_value = 0.;
  double trace_pin_time = std::numeric_limits<double>::lowest();
  /// Open-loop centre measurement from the l = 1 content of the
  /// Kretschmann radius on the excision sphere, with the previous sample so
  /// its backward difference can be logged as an independent velocity. A
  /// diagnostic only: never an input to the solve (findings 15t).
  std::array<double, 3> gb_dipole{};
  std::array<double, 3> gb_dipole_previous{};
  std::array<double, 3> gb_dipole_velocity{};
  double gb_dipole_time = std::numeric_limits<double>::lowest();
  double gb_dipole_time_previous = std::numeric_limits<double>::lowest();
  /// u^+ anchor state: the latest and previous value fits of the map to
  /// the gauge projection of the outgoing characteristic, their times,
  /// and the warm-start vector of the Gauss-Newton solve
  std::array<double, num_map_parameters> anchor_p{};
  std::array<double, num_map_parameters> anchor_p_previous{};
  double anchor_time = std::numeric_limits<double>::lowest();
  double anchor_time_previous = std::numeric_limits<double>::lowest();
  bool anchor_valid = false;
  bool valid = false;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

namespace OptionTags {
struct WorldtubeMatcher {
  using type = Options::Auto<MatcherConfig, Options::AutoLabel::None>;
  static constexpr Options::String help = {
      "Online worldtube matching of the first-order affine map from the "
      "evolved fields on the excision sphere. Set to None to disable."};
};
}  // namespace OptionTags

namespace Tags {
/// Global-cache tag holding the matcher configuration (nullopt = disabled).
struct Matcher : db::SimpleTag {
  using type = std::optional<MatcherConfig>;
  using option_tags = tmpl::list<OptionTags::WorldtubeMatcher>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// DataBox tag with the latest fitted map parameters (see
/// `MapParameterData`).
struct MapParameters : db::SimpleTag {
  using type = MapParameterData;
};
}  // namespace Tags

namespace Initialization {
/// Initializes the matcher state to its default (invalid) value.
struct InitializeMapParameters {
  using return_tags = tmpl::list<Tags::MapParameters>;
  using argument_tags = tmpl::list<>;
  using simple_tags = tmpl::list<Tags::MapParameters>;
  using compute_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  static void apply(const gsl::not_null<MapParameterData*> data) {
    *data = MapParameterData{};
  }
};
}  // namespace Initialization
}  // namespace gh::Worldtube
