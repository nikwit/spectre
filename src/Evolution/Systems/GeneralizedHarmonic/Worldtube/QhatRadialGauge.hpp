// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <optional>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

namespace PUP {
class er;
}

namespace gh::worldtube {
/// Startup qhat and its radial derivative, with spatial components expressed
/// in the frame of the Rotation function of time. Angular labels follow that
/// rotation. Supports fixed-radius spherical faces at fixed resolution,
/// including initial radial refinement. A field-only restart resets the datum.
struct QhatRadialGaugeData {
  tnsr::a<DataVector, 3> initial_q{}, initial_dr_q{};
  tnsr::I<DataVector, 3> corotating_radial_direction{};
  std::array<size_t, 3> extents{};
  double radius{0.};
  void pup(PUP::er& p);
};
bool operator==(const QhatRadialGaugeData& a, const QhatRadialGaugeData& b);
bool operator!=(const QhatRadialGaugeData& a, const QhatRadialGaugeData& b);

/// Compute qhat=ell.(Pi-n.Phi) and its full radial derivative, including
/// metric derivatives of both n and ell. The normal points into the hole.
std::pair<tnsr::a<DataVector, 3>, tnsr::a<DataVector, 3>> qhat_radial_data(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi,
    const tnsr::I<DataVector, 3>& radial);

void initialize_qhat_radial_gauge(
    gsl::not_null<std::optional<QhatRadialGaugeData>*> data,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Direction<3>& direction, double time,
    const domain::FunctionsOfTimeMap& functions_of_time);

/// Desired ell contraction of the frozen-characteristic gauge correction.
/// Implements Delta(D_c qhat)=s[dr qhat-dr qhat0-(K/R)(qhat-qhat0)] in all
/// four channels. K contains dimensionless time/radial/tangential responses.
/// Accounts for the metric correction imposed by the other boundary sectors,
/// so its changes to the null basis do not contaminate the modeled qhat RHS.
tnsr::a<DataVector, 3> qhat_radial_gauge_correction(
    const QhatRadialGaugeData& initial,
    const std::array<double, 3>& coefficients,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi,
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity,
    const Scalar<DataVector>& gamma2,
    const tnsr::aa<DataVector, 3>& metric_correction, double time,
    const domain::FunctionsOfTimeMap& functions_of_time);
}  // namespace gh::worldtube
