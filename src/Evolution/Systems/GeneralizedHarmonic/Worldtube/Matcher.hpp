// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"

/// \cond
namespace ylm {
class Spherepack;
}  // namespace ylm
/// \endcond

namespace gh::Worldtube {

/// Diagnostics and result of one online fit.
struct FitResult {
  std::array<double, num_map_parameters> p{};
  /// Fitted zeroth-order spatial offset of the map: the model is centred at
  /// `config.center + center_offset`. Zero unless `FitCenterOffset` is on.
  std::array<double, 3> center_offset{};
  double residual_initial = 0.;
  double residual_final = 0.;
  size_t iterations = 0;
};

/*!
 * \brief Fit the 13 first-order affine-map parameters from the evolved
 * fields on the excision sphere.
 *
 * Implements the findings-§15b prescription: the fit target is the
 * spherical-harmonic content (l <= `config.fit_l_max`) of the null-basis
 * gauge components {A, C, V} of the incoming characteristic field
 * \f$u^-_{ab} = \Pi_{ab} + n^k \Phi_{kab} - \gamma_2 g_{ab}\f$, with uniform
 * (absolute) weights. The model side is
 * `gh::Solutions::affine_map_model::evolved_variables` combined with the
 * *data* normal, frame, and \f$\gamma_2\f$, mirroring what the ghost
 * boundary condition applies. Pins: velocity kinematic to
 * `config.center_velocity`, trace strain to `config.trace_strain_pin`,
 * drives held at `pdot_estimate` (backward differences of the fit history).
 * Nine free parameters, solved by warm-started Gauss-Newton.
 *
 * All face tensors must be in the Spherepack collocation order of `ylm`
 * (the natural order of an excision-face slice on the Ylm-basis domain).
 */
FitResult fit_map_parameters(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const std::array<double, num_map_parameters>& p_start,
    const std::array<double, 3>& center_offset_start,
    const std::array<double, num_map_parameters>& pdot_estimate);

/// Result of one linear rate fit (the `RateOde` mode).
struct RateFitResult {
  std::array<double, num_map_parameters> pdot{};
  double residual_initial = 0.;
  double residual_final = 0.;
  /// Held-out closure: relative residual per (component class, l) block,
  /// indexed 5*class + l with class 0 = TT, 1 = Ti, 2 = ij and l = 0..4.
  /// Evaluated on the isotropy-weighted rows regardless of solve weights.
  std::array<double, 15> block_closure{};
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
}  // namespace gh::Worldtube
