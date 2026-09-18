// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <complex>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

namespace gr::np {
/*!
 * \brief The six quadrupole radial profiles of Table `tab:responses`, exact
 * in \f$M/r\f$ at first order in the moments. Mirrors
 * `profiles.quadrupole_profiles`.
 *
 * \details With \f$f = 1 - 2M/r\f$: \f$e_L = b_L = 1\f$,
 * \f$e_V = \sqrt f (1 + 2M/r)\f$, \f$b_V = \sqrt f\f$, \f$e_T = b_T = f\f$.
 * Errors if the radius lies inside the horizon.
 */
struct QuadrupoleProfiles {
  DataVector e_l;
  DataVector e_v;
  DataVector e_t;
  DataVector b_l;
  DataVector b_v;
  DataVector b_t;
};
QuadrupoleProfiles quadrupole_profiles(const Scalar<DataVector>& radius,
                                       double mass);

/// @{
/// \brief A symmetric trace-free tensor from its five direct Cartesian
/// components \f$(X_{11}, X_{22}, X_{12}, X_{13}, X_{23})\f$, with
/// \f$X_{33} = -X_{11} - X_{22}\f$. Mirrors `tides.stf_from_components`.
tnsr::ii<double, 3, Frame::Inertial> stf_from_components(
    const std::array<double, 5>& components);
tnsr::ii<std::complex<double>, 3, Frame::Inertial> stf_from_components(
    const std::array<std::complex<double>, 5>& components);
/// @}

/*!
 * \brief The model tidal tensor \f$Q^{\rm model} = E + iB\f$ at first order
 * in the moments, eq. `tidal-tensor-eps2`. Mirrors
 * `profiles.quadrupole_tide_tensor`.
 *
 * \details `electric` and `magnetic` are constant real STF tensors in the
 * same orthonormal triad as `direction`, the measured radial direction, and
 * `radius` is the measured background radius. Each \f$\hat r\f$-irreducible
 * channel (longitudinal, vector, transverse) carries its own profile. With
 * `transverse_only` the longitudinal and vector channels are omitted: that is
 * the piece entering the \f$\Psi_0\f$ and \f$\Psi_4\f$ slots of the aligned
 * frame (eq. `second-order-psi0`), the other channels being carried by the
 * measured Coulomb shift and frame rotations.
 */
ComplexMatrix quadrupole_tide_tensor(
    const tnsr::ii<double, 3, Frame::Inertial>& electric,
    const tnsr::ii<double, 3, Frame::Inertial>& magnetic,
    const TriadVector& direction, const Scalar<DataVector>& radius, double mass,
    bool transverse_only);

/*!
 * \brief SO(3,C) representation of a pure boost acting on the self-dual
 * \f$Q = E + iB\f$:
 * \f$a a^T + \cosh\eta\,(1 - a a^T) + i\sinh\eta\,[a]_\times\f$ for the unit
 * axis \f$a\f$. Mirrors `psi0.self_dual_boost`.
 */
ComplexMatrix self_dual_boost(const TriadVector& axis,
                              const Scalar<DataVector>& rapidity);

/// \brief The SO(3,C) image of the rest-frame registration: a radial boost by
/// `rapidity` along `radial_direction` followed by the transverse boost with
/// velocity `transverse_velocity`. Mirrors `psi0.rest_to_slice_map`.
ComplexMatrix rest_to_slice_map(const TriadVector& radial_direction,
                                const TriadVector& transverse_velocity,
                                const Scalar<DataVector>& rapidity);

/*!
 * \brief The five NR-frame scalar columns for the direct complex components
 * \f$(H_{11}, H_{22}, H_{12}, H_{13}, H_{23})\f$ of \f$H = E + iB\f$.
 * Mirrors `psi0.direct_tide_scalar_columns` with `transverse_only=True`.
 *
 * \details Each column is the transverse quadrupole response to a real unit
 * value of one component, built in the hole rest frame, carried to the slice
 * with `rest_to_slice_map()`, rotated into the adapted triad and read off as
 * Weyl scalars. Complex least-squares coefficients then supply the electric
 * and magnetic parts at once, because \f$e_T = b_T\f$ and the self-dual
 * action is complex linear. All directions are Cholesky-triad components.
 */
std::array<WeylScalars, 5> direct_tide_scalar_columns(
    const TriadVector& radial_direction, const TriadVector& transverse_velocity,
    const Scalar<DataVector>& rapidity, const Scalar<DataVector>& radius,
    const RealMatrix& adapted_rotation, double mass);
}  // namespace gr::np
