// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"

/// \cond
namespace ylm {
class Spherepack;
}  // namespace ylm
/// \endcond

namespace gh::Worldtube {

namespace detail {
/// Result of a block-weighted linear fit and an unweighted held-out test.
struct WeightedModeFit {
  std::vector<double> coefficients{};
  std::vector<double> fitted_residual{};
  std::vector<double> held_out_residual{};
  /// 2-norm condition number of the weighted design matrix.
  double condition_number = 0.;
};

/*!
 * \brief Solve a linear modal fit with per-block residual weights.
 *
 * `columns` is stored column-major. `row_blocks` maps every row to one entry
 * of `block_weights`. The same fitted coefficients are applied to
 * `held_out_target`, but the held-out rows never enter the solve.
 */
WeightedModeFit fit_weighted_modes(
    const std::vector<double>& target,
    const std::vector<double>& held_out_target,
    const std::vector<std::vector<double>>& columns,
    const std::vector<size_t>& row_blocks,
    const std::array<double, 15>& block_weights);

/*!
 * \brief Gauge-mode content of one GH characteristic tensor on a sphere.
 *
 * `normal_sign` is -1 for \f$u^+\f$ and +1 for \f$u^-\f$.  This helper is
 * shared by the matcher and its derivative-convention regression test.
 */
std::vector<double> characteristic_gauge_modes(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const std::array<double, 3>& center, const ylm::Spherepack& ylm_transform,
    size_t fit_l_max, double normal_sign);

/*!
 * \brief Fixed-sphere time derivative of the gauge-mode content of a GH
 * characteristic tensor.
 *
 * The inertial coordinates and sphere center are held fixed. `dt_gamma2` is
 * required because \f$\gamma_2\f$ can depend explicitly on time through its
 * functions of time even on a fixed grid.
 */
std::vector<double> characteristic_gauge_time_derivative_modes(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& dt_gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const std::array<double, 3>& center, const ylm::Spherepack& ylm_transform,
    size_t fit_l_max, double normal_sign);
}  // namespace detail

/// Diagnostics and result of one online fit.
struct FitResult {
  std::array<double, num_map_parameters> p{};
  /// Fitted zeroth-order spatial offset of the map: the model is centred at
  /// `config.center + center_offset`. Zero unless `FitCenterOffset` is on.
  std::array<double, 3> center_offset{};
  double residual_initial = 0.;
  double residual_final = 0.;
  /// Correction-relative denominator evaluated at zero free coefficients.
  double baseline_residual = 0.;
  /// Held-out opposite-characteristic baseline and final residual norms.
  double minus_baseline_residual = 0.;
  double minus_residual_final = 0.;
  /// Unweighted closure and absolute RMS per {A,C,V} x ell block.
  std::array<double, 15> block_closure{};
  std::array<double, 15> block_target_rms{};
  std::array<double, 15> block_residual_rms{};
  std::array<double, 15> block_minus_closure{};
  std::array<double, 15> block_minus_target_rms{};
  std::array<double, 15> block_minus_residual_rms{};
  /// Collocation-space covariant metric and Phi residual norms.
  double metric_baseline_residual = 0.;
  double metric_residual_final = 0.;
  double phi_baseline_residual = 0.;
  double phi_residual_final = 0.;
  /// 2-norm condition number of the weighted final design matrix.
  double condition_number = 0.;
  size_t iterations = 0;
};

/*!
 * \brief Fit the 13 first-order affine-map parameters from the evolved
 * fields on the excision sphere.
 *
 * Implements the findings-§15b prescription: the fit target is the
 * spherical-harmonic content (l <= `config.fit_l_max`) of the null-basis
 * gauge components {A, C, V} of the selected characteristic field. The
 * production value path uses the boundary-uncontrolled outgoing \f$u^+\f$;
 * \f$u^-\f$ remains available as a diagnostic convention. The model side is
 * `gh::Solutions::affine_map_model::first_order_evolved_variables` combined
 * with the
 * *data* normal, frame, and \f$\gamma_2\f$, mirroring what the ghost
 * boundary condition applies. The velocity can be pinned at strict first
 * order to `config.center_velocity` or fitted, and the trace strain can
 * likewise be pinned or fitted. Coefficient rates are absent by slow-time power
 * counting; center advection is retained through \f$\dot q^i\f$. Depending on
 * those choices and center fitting, the system has 9--16 free values, solved by
 * warm-started Gauss--Newton. `config.uplus_block_weights` weights the
 * {A,C,V} x ell residual blocks in the solve; all reported closures are
 * unweighted, and the opposite characteristic is never used in the solve.
 *
 * All face tensors must be in the Spherepack collocation order of `ylm`
 * (the natural order of an excision-face slice on the Ylm-basis domain).
 * `center_offset_start` is the Gauss--Newton initial guess when center fitting
 * is enabled and the fixed instantaneous center offset when it is disabled.
 */
FitResult fit_map_parameters(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const std::array<double, num_map_parameters>& p_start,
    const std::array<double, 3>& center_offset_start, double normal_sign,
    double trace_pin);

/// Result of one linear rate fit (the `RateOde` mode).
struct RateFitResult {
  std::array<double, num_map_parameters> pdot{};
  double residual_initial = 0.;
  double residual_final = 0.;
  /// Held-out closure: relative residual per (component class, l) block,
  /// indexed 5*class + l with class 0 = TT, 1 = Ti, 2 = ij and l = 0..4.
  /// Evaluated on the isotropy-weighted rows regardless of solve weights.
  std::array<double, 15> block_closure{};
  /// Absolute RMS target and residual for the same blocks as `block_closure`.
  std::array<double, 15> block_target_rms{};
  std::array<double, 15> block_residual_rms{};
  /// Held-out D_T u^- closure from a D_T u^+ acceleration fit. These remain
  /// zero for the metric-component rate and acceleration fits.
  std::array<double, 15> block_minus_closure{};
  std::array<double, 15> block_minus_target_rms{};
  std::array<double, 15> block_minus_residual_rms{};
  double minus_residual_initial = 0.;
  double minus_residual_final = 0.;
  /// 2-norm condition number of the weighted design matrix and norm of the
  /// unfolded fitted rate/acceleration vector.
  double condition_number = 0.;
  double parameter_derivative_norm = 0.;
  /// Largest parameter shift when the boosts are re-solved on the
  /// audit-preferred V1 blocks (TT l=1, Ti l=0) with everything else
  /// frozen at the global solution.
  double spread_vector = 0.;
  /// Same for clock rate and strain on the C3 blocks (TT l=0, ij l=2).
  double spread_clock_strain = 0.;
};

/*!
 * \brief Fit the rates \f$\dot p_A\f$ of the map parameters from the
 * time-derivative content of the boundary data.
 *
 * The data-side time derivative is \f$\partial_t g_{ab} = \beta^k
 * \Phi_{kab} - \alpha \Pi_{ab}\f$ (all data quantities); the model relation
 * is \f$\partial_t g = -\left(g\, \sum_A \dot p_A R_A\, g\right)\f$ with the
 * *data* metric in the sandwich (first-order equivalent to the model
 * metric). Both sides are projected onto spherical-harmonic modes with
 * l <= `config.fit_l_max` and the nine unpinned rates solved by linear
 * least squares; the rate pins are the time derivatives of the value pins
 * (\f$\ddot q^i = \ddot q^0 v_\text{center}\f$, trace rate zero for a
 * constant trace pin). This is the online analogue of the offline
 * "drives direct from Pi, Phi" measurement (findings §9).
 */
RateFitResult fit_map_parameter_rates(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config);

/*!
 * \brief Fit the accelerations \f$\ddot p_A\f$ of the map parameters from
 * the evolution equations (the `SecondOrderOde` mode).
 *
 * The data-side second time derivative is assembled from the GH right-hand
 * sides: \f$\partial_t^2 g = \partial_t\beta\,\Phi + \beta\,\partial_t\Phi
 * - \partial_t\alpha\,\Pi - \alpha\,\partial_t\Pi\f$, with
 * \f$\partial_t\Pi\f$, \f$\partial_t\Phi\f$, \f$\partial_t g\f$ the sliced
 * `dt` variables (the actual equations of motion) and \f$\partial_t\alpha,
 * \partial_t\beta\f$ derived from \f$\partial_t g\f$ by the 3+1 algebra.
 * The model relation, exact within the amplitude-linear map, is
 * \f$\partial_t^2 g = -\left(g\, \sum_A \ddot p_A R_A\, g\right)
 * + 2\, g S g S g\f$ with \f$S = \sum_A \dot p_A R_A\f$ built from the
 * current state `pdot_state`; the (pdot-quadratic) Hessian term is
 * subtracted from the data-side target, so the solve stays linear in
 * \f$\ddot p\f$ and the projection and acceleration-level pins are
 * identical in structure to the rate fit. The returned `pdot` member
 * holds \f$\ddot p\f$.
 */
RateFitResult fit_map_parameter_accelerations(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi,
    const std::array<double, num_map_parameters>& pdot_state,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config);

/*!
 * \brief Fit \f$\ddot p_A\f$ against the time derivative of the gauge
 * projection of the OUTGOING characteristic (the `StepperOde` +
 * `FitUPlus` mode).
 *
 * The target is \f$\partial_t\{A^+, C^+, V^+\}\f$ of
 * \f$u^+ = \Pi - n^k\Phi_k - \gamma_2 g\f$ under the same null
 * projectors as \f$u^-\f$, including the frame-motion terms (the data
 * frame and its time derivative are used on both sides). The penalty
 * drags only the \f$u^-\f$ projection of the fields, so this sensor
 * does not echo the map's own drift along gauge directions — the
 * self-consistency that produced the §15l double zero root is absent by
 * construction. The model side is evaluated at the current
 * (`p_state`, `pdot_state`) with all pieces analytic;
 * \f$\ddot p\f$ enters only through \f$\partial_t\Pi_{\rm model}\f$,
 * linearly. The solve applies `config.uplus_block_weights`; diagnostics are
 * always evaluated without those weights. `block_closure` reports fitted
 * \f$D_Tu^+\f$ {A, C, V} x l in the 15 slots, and the `block_minus_*`
 * members report the held-out \f$D_Tu^-\f$ prediction made with the same
 * fitted acceleration. The estimator-spread members are not defined for this
 * target and stay zero.
 */
RateFitResult fit_map_parameter_accelerations_uplus(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& dt_gamma2,
    const std::array<double, num_map_parameters>& p_state,
    const std::array<double, num_map_parameters>& pdot_state,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    bool evaluate_held_out_uminus);
}  // namespace gh::Worldtube
