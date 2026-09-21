// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

/// \cond
namespace gsl {
template <class T>
class not_null;
}  // namespace gsl
/// \endcond

namespace gr::np {
/*!
 * \brief The principal null pair of the type-D solve in NR adapted components
 * \f$(n, s, \hat\theta, \hat\phi)\f$, up to the type-III boost,
 * eq. `kinnersley-uN`. Mirrors `restframe.principal_null_pair`.
 *
 * \details The Kinnersley pair in NR components is the type-I rotation
 * (\f$-a\f$) followed by the type-II rotation (\f$-b\f$) applied to the NR
 * tetrad, with the tetrad parameter \f$a = \overline{\bar a}\f$. The legs are
 * exactly real by the conjugate pairing of the rotation terms.
 */
void principal_null_pair(gsl::not_null<AdaptedFourVector*> outgoing,
                         gsl::not_null<AdaptedFourVector*> incoming,
                         const Scalar<ComplexDataVector>& a_bar,
                         const Scalar<ComplexDataVector>& b);

/*!
 * \brief The boost-family member whose radial leg is tangent to the slice,
 * eq. `transverse-member`. Mirrors `restframe.TangentBoostMember`.
 *
 * \details Rescaling the null legs so their time components agree makes
 * \f$r^0 = 0\f$; then \f$\Gamma = t^0\f$, the Eulerian velocity of the member
 * is \f$w = \vec t/t^0\f$, and the measured radial direction is
 * \f$\hat r = \vec r\f$, all in adapted triad components.
 */
struct TangentBoostMember {
  TriadVector radial_direction;
  TriadVector transverse_velocity;
  Scalar<DataVector> lorentz_factor;
};

/// \brief Select the tangent member from the type-D rotation parameters.
/// Errors if the null legs are not future directed or the decoded velocity
/// is not subluminal. Mirrors `restframe.tangent_boost_member`.
TangentBoostMember tangent_boost_member(const Scalar<ComplexDataVector>& a_bar,
                                        const Scalar<ComplexDataVector>& b);

/*!
 * \brief The pointwise rapidity from the curvature gradient,
 * eq. `invariant-boost-condition`. Mirrors `restframe.invariant_rapidity`.
 *
 * \details The rest-frame member of the boost family annihilates the
 * spacetime gradient of the Kretschmann scalar \f$K\f$: with
 * \f$t(\eta) = \cosh\eta\, t_{0p} + \sinh\eta\, r_{0p}\f$ anchored at the
 * tangent member,
 * \f$\tanh\eta_p = -\Gamma(n\cdot\nabla K + w\cdot DK)/(\hat r\cdot DK)\f$
 * with \f$n\cdot\nabla K = (\partial_t K - \beta^i D_i K)/\alpha\f$.
 * `spatial_gradient` is \f$D_i K\f$ in coordinate covector components and
 * `time_derivative` is \f$\partial_t K\f$ at fixed inertial coordinates.
 * Errors if the rapidity is not physical, i.e. if \f$|\tanh\eta| \geq 1\f$
 * at any point; `invariant_tanh_rapidity()` returns \f$\tanh\eta\f$
 * without that check, for diagnostics.
 */
/// @{
Scalar<DataVector> invariant_tanh_rapidity(
    const TangentBoostMember& member, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& spatial_gradient,
    const Scalar<DataVector>& time_derivative);

Scalar<DataVector> invariant_rapidity(
    const TangentBoostMember& member, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& spatial_gradient,
    const Scalar<DataVector>& time_derivative);
/// @}
}  // namespace gr::np
