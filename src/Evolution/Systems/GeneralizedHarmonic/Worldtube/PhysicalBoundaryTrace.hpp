// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::worldtube {
/// Screen geometry shared by the physical ghost and surface evolution.
/// The interface normal points OUT of the numerical domain (into the hole).
struct PhysicalScreen {
  PhysicalScreen(const tnsr::aa<DataVector, 3>& metric,
                 const tnsr::AA<DataVector, 3>& inverse_metric,
                 const tnsr::A<DataVector, 3>& time_normal,
                 const tnsr::I<DataVector, 3>& interface_normal);
  tnsr::aa<DataVector, 3> lower;
  tnsr::Ab<DataVector, 3> mixed;
  tnsr::AA<DataVector, 3> upper;
};

/// Apply the rank-two, screen transverse-traceless projector.
tnsr::aa<DataVector, 3> physical_projection(
    const tnsr::aa<DataVector, 3>& tensor, const PhysicalScreen& screen);

/// Convert the existing NP matcher's U^{8-}=w^-/2 in covariant coordinate
/// components into the surface equation's W target. Uses the RECEIVER screen,
/// including its temporal components, exactly as the physical Bjorhus term:
/// W_ab = Q_ab^{cd} 2 P^i_c P^j_d U^{8-}_{ij}.
tnsr::aa<DataVector, 3> incoming_weyl_to_physical_target(
    const tnsr::ii<DataVector, 3>& incoming_mode,
    const PhysicalScreen& receiver_screen);

/// Reconstruct only the physical incoming characteristic. The metric is copied
/// by the caller. `target` prescribes the screen projection of Pi - s.Phi;
/// its complementary components are ignored. All geometry is the receiver's.
void physical_ghost(gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_pi,
                    gsl::not_null<tnsr::iaa<DataVector, 3>*> ghost_phi,
                    const tnsr::aa<DataVector, 3>& pi,
                    const tnsr::iaa<DataVector, 3>& phi,
                    const tnsr::aa<DataVector, 3>& target,
                    const tnsr::i<DataVector, 3>& normal_covector,
                    const tnsr::I<DataVector, 3>& normal_vector,
                    const PhysicalScreen& screen);

/// Assemble the two physical trace derivatives for gamma1=-1.
///
/// `rhs_pi` and `rhs_phi` are the volume GH RHS along the moving tube, BEFORE
/// any DG or boundary correction. `weyl` is twice the lifted incoming Weyl
/// field, INCLUDING the Bjorhus C4 additions; `target_weyl` has the same
/// normalization. The three-index constraint and the volume/curvature RHS
/// must use the SAME derivative jet. With a tangential jet the normal
/// contraction of this constraint is -s.Phi, even on an exact solution.
///
/// The polarizations are contravariant, screen tracefree and orthonormal.
/// Their derivatives and the normal derivative use the actual receiver metric
/// RHS. A fixed mesh does not in general make these derivatives vanish.
/// This operation does not evolve state or evaluate the live NP target.
std::array<DataVector, 2> physical_trace_rhs(
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& rhs_pi,
    const tnsr::iaa<DataVector, 3>& rhs_phi,
    const tnsr::iaa<DataVector, 3>& three_index_constraint,
    const tnsr::aa<DataVector, 3>& weyl,
    const tnsr::aa<DataVector, 3>& target_weyl,
    const Scalar<DataVector>& gamma2, const DataVector& lambda_minus,
    const tnsr::I<DataVector, 3>& normal_vector,
    const tnsr::I<DataVector, 3>& dt_normal_vector,
    const std::array<tnsr::AA<DataVector, 3>, 2>& polarization,
    const std::array<tnsr::AA<DataVector, 3>, 2>& dt_polarization);
}  // namespace gh::worldtube
