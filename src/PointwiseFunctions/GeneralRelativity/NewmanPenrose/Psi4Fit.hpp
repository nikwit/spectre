// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <complex>
#include <optional>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

namespace gr::np {
/*!
 * \brief The five complex components \f$(H_{11}, H_{22}, H_{12}, H_{13},
 * H_{23})\f$ of \f$H = E + iB\f$ fitted from the pulled-back
 * \f$\Psi_4\f$, eq. `psi4-fit`. Mirrors `wplus.DirectPsi4Fit`.
 */
struct DirectPsi4Fit {
  std::array<std::complex<double>, 5> components;
  /// \f$|Z x - d| / |d|\f$ of the (weighted) design
  double relative_residual;
  std::array<double, 5> electric() const;
  std::array<double, 5> magnetic() const;
};

/*!
 * \brief Solve the complex least-squares problem
 * \f$\hat x = \arg\min \sum_p w_p |d_p - Z_{pa} x_a|^2\f$ for the five
 * complex components. Mirrors `wplus.fit_psi4` /
 * `wplus._solve_complex_design`.
 *
 * \details `measured_psi4` is \f$d_p\f$, the \f$\Psi_4\f$ pulled back to the
 * Kinnersley frame, and `projected_columns` are the \f$\Psi_4\f$ slots of the
 * five tide columns pulled back the same way. Optional `point_weights` are
 * quadrature weights; rows are scaled by their square root so the collocation
 * sum becomes the \f$L^2(S^2)\f$ inner product. The columns are normalized
 * before the normal equations are solved with LAPACK, purely as a
 * preconditioner. Under a pointwise dyad spin the data and every column
 * acquire the same phase, so the solution is independent of that spin.
 */
DirectPsi4Fit fit_psi4(
    const Scalar<ComplexDataVector>& measured_psi4,
    const std::array<Scalar<ComplexDataVector>, 5>& projected_columns,
    const std::optional<DataVector>& point_weights = std::nullopt);

/*!
 * \brief The leading frame: the type-D rotation, the scalars pulled back to
 * the Kinnersley frame, the tangent boost member, the measured radius, and
 * the measured directions in Cholesky-triad components.
 *
 * \details The type-III rapidity that completes the rest frame is not part
 * of the registration; it comes from `invariant_rapidity()` and is passed to
 * `evaluate_second_order()`. `radial_direction` and `transverse_velocity`
 * are the members of `member` rotated with \f$R^T\f$ into the Cholesky
 * triad, as the tidal response functions expect.
 */
struct FrameRegistration {
  Scalar<ComplexDataVector> coulomb;
  TypeDRotation rotation;
  WeylScalars pulled_back;
  TangentBoostMember member;
  Scalar<DataVector> measured_radius;
  TriadVector radial_direction;
  TriadVector transverse_velocity;
};

/// \brief Measure the frame from the scalars `psi` in the adapted NR
/// tetrad. Mirrors `wplus.register_frame` without the Newton corrections and
/// without a boost.
FrameRegistration register_frame(const WeylScalars& psi,
                                 const RealMatrix& adapted_rotation,
                                 double mass);

/*!
 * \brief One evaluation of the minimal second-order construction. Mirrors
 * `wplus.SecondOrderEvaluation`.
 *
 * \details `psi0_target` is the boundary value of \f$\Psi_0\f$ in the NR
 * tetrad, eq. `psi0-nr` in its exact-composition form: the Kinnersley
 * scalars pushed forward to the NR tetrad, plus the fitted transverse tide
 * read off in the NR tetrad.
 */
struct SecondOrderEvaluation {
  std::array<WeylScalars, 5> direct_columns;
  DirectPsi4Fit fit;
  Scalar<ComplexDataVector> psi0_target;
};

/// \brief Fit the pulled-back \f$\Psi_4\f$ and construct the NR
/// \f$\Psi_0\f$ target for a registered frame and the type-III `rapidity`
/// of its rest frame. Mirrors `wplus.evaluate_second_order`.
SecondOrderEvaluation evaluate_second_order(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& adapted_rotation, double mass,
    const std::optional<DataVector>& fit_point_weights = std::nullopt);
}  // namespace gr::np
