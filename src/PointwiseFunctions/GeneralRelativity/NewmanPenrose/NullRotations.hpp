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
/*!
 * \brief Type-I null rotation of the Weyl scalars, eq. `type-I-scalars`
 * (\f$l\f$ fixed). Mirrors `scalars.type_i`.
 *
 * \details The parameter is \f$\bar a\f$ because the scalars transform with
 * the conjugate of the tetrad parameter \f$a\f$:
 * \f$\Psi_1 \to \Psi_1 + \bar a \Psi_0\f$,
 * \f$\Psi_2 \to \Psi_2 + 2\bar a\Psi_1 + \bar a^2\Psi_0\f$, and so on.
 */
WeylScalars type_i(const WeylScalars& psi,
                   const Scalar<ComplexDataVector>& a_bar);

/// \brief Type-II null rotation of the Weyl scalars, eq. `type-II-scalars`
/// (\f$k\f$ fixed): \f$\Psi_3 \to \Psi_3 + b\Psi_4\f$, ... Mirrors
/// `scalars.type_ii`.
WeylScalars type_ii(const WeylScalars& psi, const Scalar<ComplexDataVector>& b);

/// \brief Type-III boost and spin of the Weyl scalars, eq. `type-III-scalars`:
/// \f$\Psi_A \to e^{(2 - A)(\eta + i\chi)}\Psi_A\f$. Mirrors
/// `scalars.type_iii`.
WeylScalars type_iii(const WeylScalars& psi, const Scalar<DataVector>& eta,
                     const Scalar<DataVector>& chi);

/// @{
/*!
 * \brief The curvature invariants
 * \f$I = \Psi_0\Psi_4 - 4\Psi_1\Psi_3 + 3\Psi_2^2\f$ and
 * \f$J = \det\begin{pmatrix}\Psi_0&\Psi_1&\Psi_2\\ \Psi_1&\Psi_2&\Psi_3\\
 * \Psi_2&\Psi_3&\Psi_4\end{pmatrix}\f$, eq. `invariants`. Both are unchanged
 * under every tetrad transformation. Mirrors `scalars.invariants`.
 */
void invariants(gsl::not_null<Scalar<ComplexDataVector>*> invariant_i,
                gsl::not_null<Scalar<ComplexDataVector>*> invariant_j,
                const WeylScalars& psi);
Scalar<ComplexDataVector> invariant_i(const WeylScalars& psi);
Scalar<ComplexDataVector> invariant_j(const WeylScalars& psi);
/// @}
}  // namespace gr::np
