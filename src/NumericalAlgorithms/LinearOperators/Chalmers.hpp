// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>

class DataVector;
template <size_t Dim>
class Mesh;

/// Spectral-decay estimator used in the Chalmers hp-refinement experiment.
namespace Chalmers {
enum class FitStatus {
  BelowFloor,
  TerminatedTail,
  InsufficientModes,
  Nondecaying,
  Fit
};

struct SequenceEstimate {
  double error{0.0};
  double sigma{std::numeric_limits<double>::quiet_NaN()};
  double log_fit_rms{std::numeric_limits<double>::quiet_NaN()};
  size_t significant_modes{0};
  FitStatus status{FitStatus::BelowFloor};
};

struct DirectionEstimate {
  double error{0.0};
  double sigma{std::numeric_limits<double>::quiet_NaN()};
  double max_log_fit_rms{0.0};
  size_t insufficient_sequences{0};
  size_t ignored_sequences{0};
};

/*!
 * \brief Fit a real Legendre coefficient sequence to A exp(-sigma k).
 *
 * Implements the conventions of tools/Amr/Chalmers.py: equal-weight natural-log
 * regression on all significant coefficients, retaining original mode indices;
 * a default floor of 100 epsilon times the largest coefficient; at least three
 * significant modes; and two absent final modes interpreted as a terminated
 * tail. This last convention cannot detect aliasing of an unresolved function.
 *
 * The squared error indicator is the squared norm of the last retained mode
 * plus the extrapolated infinite tail, with Legendre norms 2/(2k+1). It is not
 * a rigorous error bound. Insufficient or nondecaying fits have infinite error.
 */
SequenceEstimate estimate_sequence(
    const DataVector& coefficients,
    std::optional<double> coefficient_floor = std::nullopt);

/*!
 * \brief Independently fit even and odd degrees with their original indices.
 *
 * Each active parity guards its last retained mode and extrapolates only modes
 * of that parity; the two indicators combine in quadrature. Two absent final
 * modes of the same parity indicate a terminated polynomial. Otherwise at
 * least three significant coefficients are required. Insufficient fits return
 * infinite error and undetermined sigma. The slowest active fit sets sigma.
 *
 * A parity whose resolved nonconstant L2 norm is at most
 * 0.01 * variation_budget / sqrt(2) is omitted from regression, but its norm
 * is retained in the indicator. Thus omitted resolved variation uses at most
 * 1% of the budget. This does not bound unsampled structure. The constant
 * coefficient is excluded only from this negligibility test.
 */
DirectionEstimate estimate_sequence_parity(const DataVector& coefficients,
                                           double variation_budget = 0.0);

/// The optional ignored_highest_modes removes filtered modes from the fit and
/// treats them as unresolved tail, rather than evidence of spectral decay.
/// Fit individual modal stripes, returning max error and min sigma per
/// direction. Only tensor-product Legendre meshes are supported. Any active
/// insufficient stripe makes the directional sigma undetermined (quiet NaN).
/// If absolute_stripe_budget is positive, omit stripes whose nonconstant part
/// has reference L2 norm below budget / sqrt(2^(Dim-1) * number_of_stripes).
/// The total norm of the omitted varying part of the interpolant is then
/// bounded by the budget. This guards against extrapolating negligible
/// transverse noise to an infinite tail; it is not an error bound for the
/// unsampled function. Zero preserves the raw fit.
template <size_t Dim>
std::array<DirectionEstimate, Dim> estimate(const DataVector& values,
                                            const Mesh<Dim>& mesh,
                                            double absolute_stripe_budget = 0.0,
                                            size_t ignored_highest_modes = 0);

/*!
 * \brief Fit the transverse L2 norm of each directional modal coefficient.
 *
 * In direction d, the sequence is
 * S_k^2 = sum_{j_perp} c_{k,j_perp}^2 prod_{i!=d} 2/(2 j_i+1).
 * The error and sigma come from the same fit to S_k. Thus a weak, poorly
 * decaying stripe cannot independently determine the directional smoothness.
 * This is an experimental extension of the Chalmers rule, distinct from the
 * individual-stripe fits in estimate(). No coefficients cancel in this norm.
 *
 * If the reference L2 norm of the nonconstant part in this direction is at
 * most absolute_variation_budget, no fit is needed. The constant mode cannot
 * require directional refinement. This bounds only the sampled interpolant's
 * variation, not the error of an unsampled function. The returned error uses
 * the full reference-element L2 normalization. There is one sequence per
 * direction, so ignored_sequences and insufficient_sequences are zero or one.
 * If separate_parity is true, uses estimate_sequence_parity instead; its
 * counts refer to the two parity subsequences.
 */
template <size_t Dim>
std::array<DirectionEstimate, Dim> estimate_directional_energy(
    const DataVector& values, const Mesh<Dim>& mesh,
    double absolute_variation_budget = 0.0, size_t ignored_highest_modes = 0,
    bool separate_parity = false);

/*!
 * \brief Estimate the error after removing one mode in each direction.
 *
 * With ignored_highest_modes, also budgets the modes that the filter will
 * remove on the smaller candidate mesh, and keeps six usable modes.
 * Each candidate is considered separately. Adds the exact reference L2 norm
 * of the removed modal slab to the directional-energy estimate on the retained
 * spectrum. This includes resolved content that would otherwise disappear
 * from the fit. A variation below absolute_variation_budget returns that
 * variation norm, rather than fitting negligible tails. Indeterminate fits
 * and meshes with six or fewer modes veto coarsening by returning infinity.
 * Keeping six modes leaves at least three modes of each parity for a later
 * refinement fit when a previously quiet direction becomes active.
 * Like the refinement indicator, this is not an error bound for an unsampled
 * function. The AMR projector uses the same modal truncation for p reduction.
 * The separate_parity switch applies to the retained-spectrum fit.
 */
template <size_t Dim>
std::array<double, Dim> p_coarsening_error(
    const DataVector& values, const Mesh<Dim>& mesh,
    double absolute_variation_budget = 0.0, size_t ignored_highest_modes = 0,
    bool separate_parity = false);

struct JoinEstimate {
  double projection_error;
  double error;
  double parent_amplitude;
};

/*!
 * \brief Check the representation obtained by joining two face siblings.
 *
 * Uses the native child-to-parent projector at the maximum child extents.
 * Measures the exact reference-parent L2 difference from the piecewise child
 * polynomials, including their interface mismatch, then adds the sum of the
 * parent's directional-energy indicators. The difference is evaluated at the
 * larger extents even when the children have unequal polynomial degrees, so
 * parent modes absent from a smaller child's mesh cannot be hidden.
 * The variation budget is divided between the parent's directions. This is a
 * sampled-polynomial loss plus a spectral indicator, not a bound on unsampled
 * data or on physical derivatives. No p reduction is included in this join.
 * The separate_parity switch applies to the parent's directional estimates.
 */
template <size_t Dim>
JoinEstimate h_coarsening_error(const DataVector& lower_values,
                                const Mesh<Dim>& lower_mesh,
                                const DataVector& upper_values,
                                const Mesh<Dim>& upper_mesh, size_t dimension,
                                double absolute_variation_budget = 0.0,
                                size_t ignored_highest_modes = 0,
                                bool separate_parity = false);
}  // namespace Chalmers
