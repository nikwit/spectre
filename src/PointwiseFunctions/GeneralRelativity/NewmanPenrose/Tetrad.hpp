// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

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
/// @{
/// \brief The lower Cholesky factor \f$L\f$ of the spatial metric,
/// \f$\gamma = L L^T\f$. Its columns' inverse transpose are the vectors of
/// the Cholesky-orthonormal triad.
void cholesky_factor(gsl::not_null<RealMatrix*> lower,
                     const tnsr::ii<DataVector, 3, Frame::Inertial>& metric);
RealMatrix cholesky_factor(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& metric);
/// @}

/// The inverse of a lower-triangular 3x3 matrix
RealMatrix inverse_lower_triangular(const RealMatrix& lower);

/*!
 * \brief Rows \f$(s, \hat\theta, \hat\phi)\f$ of the worldtube-adapted triad
 * in Cholesky-orthonormal components. Mirrors
 * `scalars.adapted_tetrad_rotation`.
 *
 * \details `directions` is the Euclidean unit normal of the coordinate sphere.
 * The sphere normal is the covector proportional to it, with orthonormal
 * components \f$L^{-1} d\f$, normalized. The tangent legs start from the
 * gauge-coordinate spherical vectors \f$\partial_\theta x\f$ and
 * \f$\partial_\phi x\f$, are carried to the orthonormal triad with
 * \f$L^T\f$ and orthonormalized there with the orientation
 * \f$s \times \hat\theta = \hat\phi\f$. At a coordinate pole the azimuth is
 * fixed to zero.
 */
RealMatrix adapted_triad(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame::Inertial>& directions);

/// Components \f$R T R^T\f$ of a rank-2 tensor in the triad whose basis
/// vectors are the rows of `rotation`. Mirrors `scalars.rotate_symmetric`.
ComplexMatrix rotate_symmetric(const ComplexMatrix& tensor,
                               const RealMatrix& rotation);

/*!
 * \brief Weyl scalars from the tidal tensor \f$Q = E + iB\f$ in the adapted
 * triad (rows \f$s, t_1, t_2\f$), eq. `scalars-from-Q`. Mirrors
 * `scalars.psi_from_q`.
 *
 * \details The dictionary is
 * \f{align*}{
 * \Psi_0 &= Q(m, m), & \Psi_1 &= -Q(s, m)/\sqrt 2, & \Psi_2 &= Q(s, s)/2,\\
 * \Psi_3 &= Q(s, \bar m)/\sqrt 2, & \Psi_4 &= Q(\bar m, \bar m),
 * \f}
 * with \f$m = (t_1 + i t_2)/\sqrt 2\f$.
 */
WeylScalars weyl_scalars_from_tidal_tensor(const ComplexMatrix& q);

/// The inverse of `weyl_scalars_from_tidal_tensor()`: the symmetric
/// trace-free \f$Q\f$ in the adapted triad. Mirrors `scalars.q_from_psi`.
ComplexMatrix tidal_tensor_from_weyl_scalars(const WeylScalars& psi);

/*!
 * \brief The Weyl scalars in the worldtube-adapted NR tetrad from the
 * electric and magnetic parts of the Weyl tensor in coordinate components.
 * Mirrors the derivation in `data.load_slice`.
 *
 * The tensor \f$E_{ij} + i B_{ij}\f$ is carried to the Cholesky triad as
 * \f$L^{-1}(E + iB)L^{-T}\f$, symmetrized and de-traced (the hygiene
 * projection of the reference), rotated into the adapted triad of
 * `adapted_triad()`, and read off with `weyl_scalars_from_tidal_tensor()`.
 */
WeylScalars weyl_scalars_from_electric_magnetic(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame::Inertial>& directions);

/*!
 * \brief The incoming characteristic field of the Weyl tensor,
 * eq. `incoming-weyl`,
 * \f$w^-_{ij} = 2(\bar\Psi_0 m_i m_j + \Psi_0 \bar m_i \bar m_j)\f$,
 * in Cholesky-orthonormal components. Mirrors `scalars.w_minus`.
 *
 * \details Real, symmetric, trace free and tangent to the cut, and invariant
 * under the dyad spin \f$m \to e^{i\chi} m\f$. Pointwise
 * \f$|w^-|_F^2 = 8|\Psi_0|^2\f$.
 */
tnsr::ii<DataVector, 3, Frame::Inertial> incoming_weyl_field(
    const Scalar<ComplexDataVector>& psi0, const RealMatrix& rotation);

/// Covariant coordinate components \f$T_{ij} = L_{ia} L_{jb} T_{ab}\f$ of a
/// symmetric tensor given in Cholesky-orthonormal components.
tnsr::ii<DataVector, 3, Frame::Inertial> orthonormal_to_coordinate_covariant(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& orthonormal,
    const RealMatrix& cholesky);
}  // namespace gr::np
