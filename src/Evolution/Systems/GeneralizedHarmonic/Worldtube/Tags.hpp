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
        "Fallback inertial center used by direct matcher calls. During an "
        "online evolution the matcher measures the instantaneous excision-"
        "sphere center from its inertial collocation coordinates instead, so "
        "a time-dependent domain map remains independent of the fit."};
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
        "advanced with the fitted finite bulk velocity when FitBulkBoost is "
        "active, or otherwise its fitted O(epsilon) velocity; derivative ODE "
        "modes retain their higher-order extrapolation."};
    static type lower_bound() { return 0.; }
  };
  struct FitRadialIndex {
    using type = size_t;
    static constexpr Options::String help = {
        "Radial collocation index of the worldtube element at which the fit "
        "reads the fields: 0 is the excision face itself. With FitUPlus, "
        "index 0 is valid because the boundary condition does not prescribe "
        "u+; this is the direct face-matching configuration. If fitting u- "
        "instead, index 0 is tautological because the boundary condition "
        "clamps that gauge characteristic, so an off-face index is then "
        "required for an independent signal."};
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
  struct FitBulkBoost {
    using type = bool;
    static constexpr Options::String help = {
        "Before fitting the first-order affine residual, determine a separate "
        "finite Lorentz-boost velocity from the covariant metric on the "
        "sampling sphere. The boost is exact in velocity and the residual "
        "remains strict first order. This is implemented only for the value "
        "fit and is mutually exclusive with FitVelocity: otherwise the "
        "tangent of the finite boost would be fitted a second time as beta_i "
        "and qdot^i."};
  };
  struct FitExactFrame {
    using type = bool;
    static constexpr Options::String help = {
        "Replace the linear 13-parameter value fit with the zeroth-order "
        "exact-frame model: one nonlinear Gauss-Newton fit of the finite "
        "frame map L = B(rapidity) S(s0, sigma, s_ij) -- the exact "
        "pushforward of harmonic Schwarzschild, no linear residual stage -- "
        "to the l <= FitLMax gauge projection of u+. All 13 frame parameters "
        "are free. CenterVelocity and the trace pin act only as cold-start "
        "priors (V ~ CenterVelocity, s0 ~ -pin, s_ij ~ pin delta_ij), never "
        "as constraints. Requires FitUPlus; mutually exclusive with "
        "FitBulkBoost, FitVelocity, FitTraceStrain, FitCenterOffset, and the "
        "ODE modes. Between fits the model centre advances with the exact "
        "coordinate centre velocity V_c = L^i_0/L^0_0 regardless of "
        "CentreAdvection (the exact model contains its advection by "
        "construction, so that A/B switch does not apply)."};
  };
  struct FitVelocitySeparately {
    using type = bool;
    static constexpr Options::String help = {
        "Exact-frame mode only: split the 13-parameter Gauss-Newton solve "
        "into two stages per fit instant -- first the boost rapidity alone "
        "with the symmetric factor S frozen at its warm start, then the ten "
        "S parameters with the rapidity frozen. Each stage sees only its own "
        "residual response, so the near-degenerate combined direction "
        "delta V = delta sigma (which leaves the coordinate centre velocity "
        "V_c unchanged) cannot be traversed within a single fit cycle. "
        "Requires FitExactFrame."};
  };
  struct PinSymmetricFactor {
    using type = bool;
    static constexpr Options::String help = {
        "Exact-frame mode only: pin the eta-symmetric frame factor to the "
        "identity, S = 1 (s0 = sigma_i = s_ij = 0), and fit only the three "
        "boost rapidity components. Physically exact for an isolated "
        "(companion-free) hole, where the symmetric sector carries no "
        "content; removes the near-degenerate delta V = delta sigma "
        "direction by construction. Do not use with a companion present: "
        "S = 1 then discards the O(m1/a) uniform-potential sector. Requires "
        "FitExactFrame; mutually exclusive with FitVelocitySeparately."};
  };
  struct FitRadialDerivative {
    using type = bool;
    static constexpr Options::String help = {
        "Exact-frame mode only: append the radial derivative of u^+ at the "
        "fit shell to the residual, formed on both the data and model side "
        "with the same Lagrange differentiation row over all radial "
        "collocation shells of the boundary element. The derivative rows "
        "carry the first-order radial-profile information that "
        "distinguishes the boost from the simultaneity mixing sigma, which "
        "the single-shell values only separate at second order. Requires "
        "FitExactFrame."};
  };
  struct RadialDerivativeWeight {
    using type = double;
    static constexpr Options::String help = {
        "Overall weight multiplying the radial-derivative rows relative to "
        "the value rows (dimensionally a length scale; 1.0 weights d_r u^+ "
        "in units of the mass). Only used with FitRadialDerivative."};
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

  struct ExcisionSphereName {
    using type = std::string;
    static constexpr Options::String help = {
        "Name of the excision sphere in the domain whose abutting block "
        "hosts the worldtube boundary element the matcher reads. "
        "'ExcisionSphere' for the single-hole Sphere/SphericalShells "
        "domains; 'ExcisionSphereB' for the small hole in a "
        "BinaryCompactObject domain."};
    static type suggested_value() { return "ExcisionSphere"; }
  };

  using options =
      tmpl::list<Mass, Center, CenterVelocity, TraceStrainPin, FitLMax,
                 FitInterval, FitCenterOffset, RateOde, SecondOrderOde,
                 StepperOde, GaugeDamping, UPlusAnchor, FitUPlus,
                 KretschmannTracePin, TracePinInterval, FitTraceStrain,
                 FitVelocity, FitBulkBoost, FitExactFrame,
                 FitVelocitySeparately, PinSymmetricFactor, FitRadialDerivative,
                 RadialDerivativeWeight, CentreAdvection, SpatialMonopoleWeight,
                 UPlusBlockWeights, FitRadialIndex, ExcisionSphereName>;
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
                bool fit_velocity, bool fit_bulk_boost, bool fit_exact_frame,
                bool fit_velocity_separately, bool pin_symmetric_factor,
                bool fit_radial_derivative, double radial_derivative_weight,
                bool centre_advection, double spatial_monopole_weight,
                const std::array<double, 15>& uplus_block_weights,
                size_t fit_radial_index,
                std::string excision_sphere_name = "ExcisionSphere")
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
        fit_bulk_boost(fit_bulk_boost),
        fit_exact_frame(fit_exact_frame),
        fit_velocity_separately(fit_velocity_separately),
        pin_symmetric_factor(pin_symmetric_factor),
        fit_radial_derivative(fit_radial_derivative),
        radial_derivative_weight(radial_derivative_weight),
        centre_advection(centre_advection),
        spatial_monopole_weight(spatial_monopole_weight),
        uplus_block_weights(uplus_block_weights),
        fit_radial_index(fit_radial_index),
        excision_sphere_name(std::move(excision_sphere_name)) {
    if (fit_velocity and (rate_ode or second_order_ode or stepper_ode)) {
      ERROR(
          "FitVelocity is implemented for the value fit only: the rate and "
          "acceleration solves carry their own velocity pin (qddot^i = "
          "qddot^0 v_centre) and a fixed nine-column layout, so enabling "
          "both would silently keep the velocity pinned. Set RateOde, "
          "SecondOrderOde and StepperOde false.");
    }
    if (fit_bulk_boost and (rate_ode or second_order_ode or stepper_ode)) {
      ERROR(
          "FitBulkBoost is implemented for the strict value fit only. Set "
          "RateOde, SecondOrderOde and StepperOde false.");
    }
    if (fit_bulk_boost and fit_velocity) {
      ERROR(
          "FitBulkBoost and FitVelocity are mutually exclusive: the latter "
          "would duplicate the tangent of the finite Lorentz boost in the "
          "first-order residual.");
    }
    if (fit_bulk_boost and
        center_velocity != std::array<double, 3>{{0., 0., 0.}}) {
      ERROR(
          "FitBulkBoost requires CenterVelocity = [0, 0, 0]. The finite "
          "boost owns the bulk center motion; a separate velocity pin would "
          "double count it.");
    }
    if (fit_exact_frame and (rate_ode or second_order_ode or stepper_ode)) {
      ERROR(
          "FitExactFrame is a value fit. Set RateOde, SecondOrderOde and "
          "StepperOde false.");
    }
    if (fit_exact_frame and (fit_bulk_boost or fit_velocity)) {
      ERROR(
          "FitExactFrame owns the full frame including the boost, so "
          "FitBulkBoost and FitVelocity would fit the same velocity a second "
          "time. Set both false.");
    }
    if (fit_exact_frame and not fit_uplus) {
      ERROR(
          "FitExactFrame fits the gauge projection of the outgoing "
          "characteristic; set FitUPlus true.");
    }
    if (fit_exact_frame and fit_center_offset) {
      ERROR(
          "FitExactFrame takes the centre from the tracked worldtube (the "
          "Gauss-Bonnet track); a fitted centre offset is not part of the "
          "zeroth-order model. Set FitCenterOffset false.");
    }
    if (fit_exact_frame and fit_trace_strain) {
      ERROR(
          "FitExactFrame always fits the full symmetric strain including its "
          "trace; FitTraceStrain configures the linear path only. Set it "
          "false.");
    }
    if (fit_velocity_separately and not fit_exact_frame) {
      ERROR(
          "FitVelocitySeparately splits the exact-frame Gauss-Newton solve "
          "and requires FitExactFrame.");
    }
    if (pin_symmetric_factor and not fit_exact_frame) {
      ERROR(
          "PinSymmetricFactor pins the exact-frame symmetric factor and "
          "requires FitExactFrame.");
    }
    if (pin_symmetric_factor and fit_velocity_separately) {
      ERROR(
          "PinSymmetricFactor leaves only the boost stage, so "
          "FitVelocitySeparately has nothing to split. Set one of them "
          "false.");
    }
    if (fit_radial_derivative and not fit_exact_frame) {
      ERROR(
          "FitRadialDerivative extends the exact-frame fit and requires "
          "FitExactFrame.");
    }
    if (fit_radial_derivative and radial_derivative_weight <= 0.) {
      ERROR("RadialDerivativeWeight must be positive, got "
            << radial_derivative_weight << ".");
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
  bool fit_bulk_boost = false;
  bool fit_exact_frame = false;
  bool fit_velocity_separately = false;
  bool pin_symmetric_factor = false;
  bool fit_radial_derivative = false;
  double radial_derivative_weight = 1.0;
  bool centre_advection = true;
  double spatial_monopole_weight = 1.;
  std::array<double, 15> uplus_block_weights{
      {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}};
  size_t fit_radial_index = 0;
  std::string excision_sphere_name{"ExcisionSphere"};
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
  /// Finite, nonperturbative Lorentz velocity of the Schwarzschild
  /// background. This is separate from the epsilon-order affine residual.
  std::array<double, 3> bulk_velocity{};
  /// Instantaneous inertial center of the worldtube sampling sphere, measured
  /// from the l=0 modes of its inertial collocation coordinates. This is
  /// numerical-domain geometry, not a black-hole-center measurement.
  std::array<double, 3> worldtube_center{};
  /// Worldtube center at `last_fit_time`. Together with `worldtube_center`,
  /// this removes the domain-map displacement from the advected fit offset.
  std::array<double, 3> worldtube_center_at_last_fit{};
  bool worldtube_center_valid = false;
  /// Black-hole center minus `worldtube_center_at_last_fit`. Thus this stays
  /// small when the control system tracks the hole; it is not the absolute
  /// inertial black-hole position.
  std::array<double, 3> center_offset{};
  /// Zeroth-order exact-frame state (FitExactFrame mode): the fitted frame
  /// parameters ordered (rapidity[3], s0, sigma[3], s_(ij)[6]) and the
  /// derived coordinate centre velocity V_c = L^i_0/L^0_0, which advects the
  /// model centre between fits (spec Eq. Z4).
  std::array<double, num_map_parameters> exact_frame_theta{};
  std::array<double, 3> exact_frame_center_velocity{};
  bool exact_frame_valid = false;
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
  /// its backward difference can be logged as the hole--worldtube relative
  /// velocity. Adding the independently logged worldtube-center velocity
  /// gives an inertial velocity estimate. A diagnostic only: never an input
  /// to the solve (findings 15t).
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
