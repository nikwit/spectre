// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::Solutions::order_by_order_worldtube {

/// Number of first-order time-dependent affine-map profiles.
static constexpr size_t num_affine_rates = 16;

/// Coefficients of the first-order profiles, ordered as
///
/// \f$(\ddot q^0, \dot\beta_i, \ddot q^i,
///       \dot\sigma_{xx}, \dot\sigma_{xy}, \dot\sigma_{xz},
///       \dot\sigma_{yy}, \dot\sigma_{yz}, \dot\sigma_{zz},
///       \dot\omega_{xy}, \dot\omega_{xz}, \dot\omega_{yz})\f$.
using AffineRates = std::array<double, num_affine_rates>;

/*!
 * \brief First-order inverse-metric response to time-dependent affine-map
 * parameters in the local harmonic-Schwarzschild chart.
 *
 * This is the 16-column order-one model used by the q=8 matching analysis:
 * one clock acceleration, three time-gradient rates, three center
 * accelerations, six symmetric spatial-strain rates, and three
 * rotation/precession rates. The response is exactly linear in `rates`.
 * `local_coords` are spatial harmonic coordinates relative to the black-hole
 * center and must lie outside the harmonic horizon, \f$\rho>M\f$.
 */
void local_affine_rate_inverse_metric_response(
    gsl::not_null<tnsr::AA<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& local_coords, double mass,
    const AffineRates& rates);

/*!
 * \brief The local first-order rate response composed through a finite exact
 * frame map.
 *
 * For the event \f$X-z=(T,\mathbf X-\mathbf z)\f$, this computes
 *
 * \f[
 * \delta G^{AB}(X) = L^A{}_{\mu} L^B{}_{\nu}
 *   \delta g^{\mu\nu}((L^{-1}(X-z))^i).
 * \f]
 *
 * The finite frame `frame_map_matrix` is not expanded. This function returns
 * only the first-order response, not the zeroth-order background.
 */
void affine_rate_inverse_metric_response(
    gsl::not_null<tnsr::AA<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& x, double time, double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix, const AffineRates& rates);

/*!
 * \brief Exact finite-frame harmonic-Schwarzschild inverse metric plus the
 * linear first-order affine-rate response.
 *
 * This is the value-level order-zero-plus-order-one model. It intentionally
 * contains neither second-order profiles nor a nonlinear resummation of the
 * first-order coefficients.
 */
void inverse_metric(gsl::not_null<tnsr::AA<DataVector, 3>*> result,
                    const tnsr::I<DataVector, 3>& x, double time, double mass,
                    const std::array<double, 3>& center,
                    const exact_frame::FrameMatrix& frame_map_matrix,
                    const AffineRates& rates);

/*!
 * \brief Strict first-order GH-variable response to the affine rates on the
 * finite exact-frame background.
 *
 * Returns only \f$(\delta g_{ab},\delta\Pi_{ab},\delta\Phi_{iab})\f$.
 * Spatial and time derivatives of the inverse-metric response are analytic;
 * the response is linearized once about the exact zeroth-order metric, so no
 * products of affine rates are retained.
 */
void affine_rate_evolved_variables_response(
    gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric_response,
    gsl::not_null<tnsr::aa<DataVector, 3>*> pi_response,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> phi_response,
    const tnsr::I<DataVector, 3>& x, double time, double mass,
    const std::array<double, 3>& center,
    const exact_frame::FrameMatrix& frame_map_matrix, const AffineRates& rates);

/// Exact zeroth-order finite-frame GH variables plus the strict linear
/// first-order affine-rate response.
void evolved_variables(gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
                       gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
                       gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
                       const tnsr::I<DataVector, 3>& x, double time,
                       double mass, const std::array<double, 3>& center,
                       const exact_frame::FrameMatrix& frame_map_matrix,
                       const AffineRates& rates);

}  // namespace gh::Solutions::order_by_order_worldtube
