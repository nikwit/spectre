// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

/// Items shared by the worldtube boundary condition, the action supplying its
/// curvature-gradient data and the diagnostics event
namespace gh::worldtube {
/*!
 * \brief The electric and magnetic parts of the Weyl tensor from the
 * generalized harmonic variables and their spatial derivatives.
 *
 * \details The spatial Ricci tensor comes from \f$\Phi\f$ and its derivative
 * and the covariant derivative of the extrinsic curvature from
 * `gh::covariant_deriv_of_extrinsic_curvature()`; then
 * `gr::weyl_electric()` and `gr::weyl_magnetic()`. Unlike the physical
 * Bjorhus term, no four-index-constraint terms are added to the Ricci tensor:
 * the curvature is consumed as it is.
 */
void weyl_electric_magnetic(
    gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*> electric,
    gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*> magnetic,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::A<DataVector, 3, Frame::Inertial>& spacetime_unit_normal_vector,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& extrinsic_curvature,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric);

/// The 3+1 quantities and the Weyl curvature derived from the evolved
/// generalized harmonic variables at a set of points
struct WeylCurvature {
  tnsr::ii<DataVector, 3, Frame::Inertial> spatial_metric;
  tnsr::II<DataVector, 3, Frame::Inertial> inverse_spatial_metric;
  Scalar<DataVector> lapse;
  tnsr::I<DataVector, 3, Frame::Inertial> shift;
  tnsr::ii<DataVector, 3, Frame::Inertial> electric;
  tnsr::ii<DataVector, 3, Frame::Inertial> magnetic;
};

/// \brief `weyl_electric_magnetic()` starting from the evolved variables
/// \f$g_{ab}, \Pi_{ab}, \Phi_{iab}\f$ and the spatial derivatives of
/// \f$\Pi\f$ and \f$\Phi\f$, deriving the 3+1 quantities on the way.
WeylCurvature weyl_curvature(
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi);

/*!
 * \brief The Kretschmann scalar \f$K = C_{abcd} C^{abcd}\f$ of a vacuum
 * spacetime from the electric and magnetic parts of the Weyl tensor,
 * \f$K = 8 (E_{ij} E^{ij} - B_{ij} B^{ij})\f$.
 *
 * \details This is \f$16\,\mathrm{Re}\, I\f$ with \f$I\f$ the quadratic Weyl
 * invariant of the NP-matching note, the field whose spacetime gradient
 * fixes the invariant rest frame (`gr::np::invariant_rapidity()`).
 */
Scalar<DataVector> kretschmann_scalar(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inverse_spatial_metric);
}  // namespace gh::worldtube
