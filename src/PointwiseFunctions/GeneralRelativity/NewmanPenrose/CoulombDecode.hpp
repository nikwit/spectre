// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <optional>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

namespace gr::np {
/*!
 * \brief The tide of one sphere decoded from the Coulomb channel, see
 * `decode_tidal_moments_from_coulomb()`. Mirrors `coulomb.CoulombDecode`.
 *
 * `components` are the five complex direct components of \f$E + iB\f$ in the
 * convention of `direct_tide_scalar_columns()`, usable as the imposed
 * components of `evaluate_second_order()`. `valid` is false when a point
 * reached the turning point \f$r = 9M/4\f$ of the radius solve or the solve
 * did not converge; the other members are then meaningless.
 */
struct CoulombDecode {
  bool valid;
  Scalar<DataVector> areal_radius;
  TidalMoments components;
  /// The fitted \f$l = 0\f$ and \f$l = 1\f$ Coulomb offsets (a mass
  /// mismatch and a residual dipole)
  std::array<std::complex<double>, 4> nuisance;
  double relative_residual;
  double newton_residual;
};

/// The radius of the turning point of \f$B(r)\f$ in units of the mass
constexpr double coulomb_decode_turning_point_over_mass = 2.25;

/// \brief \f$B(r) = \sqrt{1 - 2M/r}\, 3M/r^4\f$, the proper radial
/// derivative of the Coulomb scalar of a hole at rest. Mirrors
/// `coulomb.background_radial_derivative`.
DataVector background_radial_derivative(const DataVector& radius, double mass);

/// \brief \f$\cosh\eta\, (\hat r\cdot s) + \sinh\eta\, \Gamma\, (w\cdot s)\f$:
/// the component along the sphere normal of the rest-frame radial unit
/// vector, for the tangent member in adapted-triad components. Mirrors
/// `coulomb.normal_derivative_factor`.
DataVector normal_derivative_factor(const TangentBoostMember& member,
                                    const Scalar<DataVector>& rapidity);

/*!
 * \brief Solve \f$B(r)\,|\text{factor}| = |\partial_s \mathrm{Re}\,\Psi_2^K|\f$
 * for the areal radius \f$r\f$ of every point on the outer branch of
 * \f$B\f$, by Newton's method from `initial_radius`. Mirrors
 * `coulomb.radius_from_normal_derivative`.
 *
 * Returns the radius, whether every point stayed above the turning point
 * and the solve converged, and the largest relative residual of the solve.
 */
struct RadiusSolve {
  DataVector radius;
  bool valid;
  double newton_residual;
};
RadiusSolve radius_from_normal_derivative(
    const DataVector& normal_derivative_of_coulomb_real,
    const DataVector& factor, const DataVector& initial_radius, double mass,
    size_t iterations = 40);

/*!
 * \brief The response of the Coulomb scalar \f$-3J/I\f$ to the ten real
 * unit moments: the five electric direct components, then the five magnetic
 * ones. Mirrors `coulomb.coulomb_tide_columns`.
 *
 * \details \f$-3J/I\f$ is a tetrad invariant, so the tide's contribution to
 * it is the rest-frame value \f$\tfrac12 Q(\hat R, \hat R) = \tfrac12 e_L
 * H(\hat R, \hat R)\f$ along the rest-frame radial direction, with no boost
 * mixing: a boost rearranges the scalars between the slots of a tetrad, not
 * the invariant. `radial_direction` is the measured radial direction in the
 * Cholesky frame, which the tidal model identifies with the rest-frame
 * radial direction; `radius` enters only through the profile \f$e_L\f$
 * (which is 1 for the quadrupole).
 */
std::array<Scalar<ComplexDataVector>, 10> coulomb_tide_columns(
    const TriadVector& radial_direction, const Scalar<DataVector>& radius,
    double mass);

/*!
 * \brief Decode the ten quadrupole moments of a sphere from the Coulomb
 * scalar and the real part of its derivative along the sphere normal.
 * Mirrors `coulomb.decode_tidal_moments_from_coulomb`.
 *
 * \details On a single sphere the electric tide is degenerate with a
 * displacement of the sphere in the Coulomb scalar alone,
 * \f$\Psi_2^K = -M/r^3 + \tfrac12 E_{nn}\f$ with an unknown areal-radius
 * field. The normal derivative breaks the degeneracy: a static tide does not
 * change it at first order (profile \f$e_L = 1\f$), while for the type-D
 * background it is \f$B(r)\f$ times the factor of
 * `normal_derivative_factor()`, exactly and pointwise. The radius solve
 * therefore yields the areal radius of every face point free of the tide,
 * whatever the offset, shape and motion of the coordinate sphere, and the
 * tide is the pointwise excess \f$\Psi_2^K + M/r^3\f$, fitted to
 * `coulomb_tide_columns()` with \f$l = 0, 1\f$ nuisance terms by weighted
 * least squares. The radiative modes enter only through the invariants,
 * quadratically, which is what makes this decode usable in a boundary
 * condition where the \f$\Psi_4\f$ fit feeds back on itself.
 *
 * The map \f$r \to B(r)\f$ turns at \f$r = 9M/4\f$, so the sphere must lie
 * outside about \f$2.4M\f$; the sensitivity is \f$1/|(M/r)/f - 4|\f$.
 */
CoulombDecode decode_tidal_moments_from_coulomb(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    double mass, const Scalar<DataVector>& normal_derivative_of_coulomb_real,
    const std::optional<DataVector>& point_weights = std::nullopt);
}  // namespace gr::np
