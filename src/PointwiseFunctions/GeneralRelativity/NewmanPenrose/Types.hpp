// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"

/*!
 * \ingroup GeneralRelativityGroup
 * \brief Pointwise Newman-Penrose algebra for the curvature matching at a
 * worldtube boundary.
 *
 * \details Implements, function by function, the pure-NumPy reference
 * `npmatch` of the worldtube NP-matching study: the worldtube-adapted tetrad
 * and Weyl scalars, the null rotations and curvature invariants, the type-D
 * (Kinnersley) frame solve, the hole rest frame from the curvature gradient,
 * the tidal response profiles, and the least-squares fit of the
 * quadrupole moments through \f$\Psi_4\f$ that supplies the boundary value of
 * \f$\Psi_0\f$. Every function names the reference function it mirrors and is
 * tested against it through pypp.
 *
 * Conventions follow the reference: spatial components are taken in
 * orthonormal triads. The Cholesky triad of the spatial metric
 * \f$\gamma = L L^T\f$ is the working triad, and the worldtube-adapted triad
 * has rows \f$(s, \hat\theta, \hat\phi)\f$ with \f$s\f$ the unit outward
 * normal of the excision sphere. The NR tetrad is
 * \f$l = (n + s)/\sqrt 2\f$, \f$k = (n - s)/\sqrt 2\f$,
 * \f$m = (\hat\theta + i\hat\phi)/\sqrt 2\f$, and the five complex Weyl
 * scalars are stored with index \f$A = 0..4\f$ for \f$\Psi_A\f$.
 */
namespace gr::np {
/// The five Weyl scalars \f$\Psi_0, \ldots, \Psi_4\f$ at each point
using WeylScalars = tnsr::a<ComplexDataVector, 4, Frame::Inertial>;
/// A real 3x3 matrix at each point, e.g. the rows of a triad
using RealMatrix = tnsr::ij<DataVector, 3, Frame::Inertial>;
/// A complex 3x3 matrix at each point, e.g. the tidal tensor \f$E + iB\f$
using ComplexMatrix = tnsr::ij<ComplexDataVector, 3, Frame::Inertial>;
/// A spatial vector in triad components at each point
using TriadVector = tnsr::I<DataVector, 3, Frame::Inertial>;
/// A four-vector in the adapted components \f$(n, s, \hat\theta, \hat\phi)\f$
using AdaptedFourVector = tnsr::A<DataVector, 3, Frame::Inertial>;
}  // namespace gr::np
