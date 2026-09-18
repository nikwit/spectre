// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

namespace gr::np {
/*!
 * \brief The Kinnersley Coulomb scalar \f$\Psi_2^K = -3J/I\f$,
 * eq. `typeD-IJ`. Mirrors `kinnersley.coulomb_scalar`.
 *
 * \details Exact on a type-D background and correct to \f$O(\varepsilon^2)\f$
 * with the tide. For Schwarzschild \f$\Psi_2^K = -M/r^3 < 0\f$.
 */
Scalar<ComplexDataVector> coulomb_scalar(
    const Scalar<ComplexDataVector>& invariant_i,
    const Scalar<ComplexDataVector>& invariant_j);

/// \brief The background areal radius \f$r = (M I / 3J)^{1/3}\f$,
/// eq. `radius`, on the principal branch: the real part is the radius and
/// the imaginary part an \f$O(\varepsilon^2)\f$ diagnostic. Mirrors
/// `kinnersley.background_radius`.
Scalar<ComplexDataVector> background_radius(
    const Scalar<ComplexDataVector>& invariant_i,
    const Scalar<ComplexDataVector>& invariant_j, double mass);

/// The Kinnersley-frame scalars \f$(0, 0, \Psi_2^K, 0, 0)\f$,
/// eq. `kinnersley-scalars`. Mirrors `kinnersley.kinnersley_scalars`.
WeylScalars kinnersley_scalars(const Scalar<ComplexDataVector>& coulomb);

/// \brief The type-D scalars in a generic tetrad, eq. `typeD-generic`: a
/// type-II rotation (\f$b\f$) followed by a type-I rotation (\f$\bar a\f$)
/// applied to the Kinnersley scalars. Mirrors `kinnersley.type_d_scalars`.
WeylScalars type_d_scalars(const Scalar<ComplexDataVector>& coulomb,
                           const Scalar<ComplexDataVector>& a_bar,
                           const Scalar<ComplexDataVector>& b);

/*!
 * \brief The null rotations from the Kinnersley frame to the NR tetrad,
 * fixed by the \f$\Psi_1\f$ and \f$\Psi_2\f$ equations of eq. `typeD-generic`.
 *
 * \details `a_bar` and `b` carry the undetermined type-III factor; only
 * `x` \f$= \bar a b\f$ and `predicted_psi` are convention free. The latter
 * are the five type-D scalars in the NR tetrad: slots 1 and 2 match the
 * measured scalars by construction, and slots 0, 3, 4 differ from them at
 * \f$O(\varepsilon^2)\f$ -- the tide.
 */
struct TypeDRotation {
  Scalar<ComplexDataVector> a_bar;
  Scalar<ComplexDataVector> b;
  Scalar<ComplexDataVector> x;
  WeylScalars predicted_psi;
};

/*!
 * \brief Fix \f$(\bar a, b)\f$ from the \f$\Psi_1\f$ and \f$\Psi_2\f$
 * equations. Mirrors `kinnersley.solve_type_d_rotation`.
 *
 * \details The \f$\Psi_2\f$ equation is the quadratic
 * \f$\bar a b = (-6 \pm \sqrt{36 - 24(1 - \Psi_2/\Psi_2^K)})/12\f$, whose
 * near-zero root aligns the ingoing and outgoing null legs with those of the
 * NR tetrad. The \f$\Psi_1\f$ equation,
 * \f$\Psi_1 = 3b(1 + 2\bar a b)\Psi_2^K\f$, then gives \f$b\f$, and
 * \f$\bar a = x/b\f$. Errors if the aligning root exceeds 0.5 in magnitude
 * (the tetrad is longitudinally exchanged) or if \f$\Psi_1 = 0\f$ while
 * \f$x \neq 0\f$.
 */
TypeDRotation solve_type_d_rotation(const WeylScalars& psi,
                                    const Scalar<ComplexDataVector>& coulomb);

/// \brief Apply the inverse of the Kinnersley rotations: type I
/// (\f$-\bar a\f$) followed by type II (\f$-b\f$). On type-D scalars this
/// lands on `kinnersley_scalars()`; on data it exposes the
/// \f$O(\varepsilon^2)\f$ longitudinal residuals in slots 1 and 3. Mirrors
/// `kinnersley.pull_back`.
WeylScalars pull_back(const WeylScalars& psi,
                      const Scalar<ComplexDataVector>& a_bar,
                      const Scalar<ComplexDataVector>& b);

/// \brief The inverse of `pull_back()`: type II (\f$b\f$) followed by type I
/// (\f$\bar a\f$), carrying scalars from the Kinnersley frame to the NR
/// tetrad. `push_forward(pull_back(psi)) == psi`.
WeylScalars push_forward(const WeylScalars& psi,
                         const Scalar<ComplexDataVector>& a_bar,
                         const Scalar<ComplexDataVector>& b);

/// The leading-order boundary value \f$\Psi_0^{\rm BC} = 6 b^2 \Psi_2^K\f$,
/// eq. `psi0-leading`. Mirrors `kinnersley.psi0_leading`.
Scalar<ComplexDataVector> psi0_leading(const Scalar<ComplexDataVector>& coulomb,
                                       const Scalar<ComplexDataVector>& b);
}  // namespace gr::np
