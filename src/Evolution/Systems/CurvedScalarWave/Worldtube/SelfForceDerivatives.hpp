// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/Tensor/Tensor.hpp"

namespace CurvedScalarWave::Worldtube {

/// Values and spatial derivatives of the two quadratic regular-field
/// expansions at the particle. All derivatives are inertial partial
/// derivatives; the Hessians are twice the quadratic Taylor coefficients.
struct RegularFieldDerivatives {
  double psi;
  double dt_psi;
  tnsr::i<double, 3> d_psi;
  tnsr::i<double, 3> d_dt_psi;
  tnsr::ii<double, 3> d2_psi;
  tnsr::ii<double, 3> d2_dt_psi;
};

/// Reconstruct the full quadratic fields from the boundary STF coefficients
/// and evolved central value and its total time derivative. In particular,
/// boundary monopoles are not central values at expansion order two.
RegularFieldDerivatives reconstruct_regular_field(
    double boundary_psi_monopole, double boundary_dt_psi_monopole,
    const tnsr::i<double, 3>& psi_dipole,
    const tnsr::i<double, 3>& dt_psi_dipole,
    const tnsr::ii<double, 3>& psi_quadrupole,
    const tnsr::ii<double, 3>& dt_psi_quadrupole, double psi0,
    double total_dt_psi0, const tnsr::I<double, 3>& velocity,
    double worldtube_radius);

/// Acceleration data for the puncture, with the turn-on function included.
/// dt denotes the total coordinate-time derivative along the worldline.
struct SelfForceDerivatives {
  tnsr::I<double, 3> coordinate_acceleration;
  tnsr::A<double, 3> f;
  tnsr::A<double, 3> dt_f;
  tnsr::A<double, 3> Du_f;
  tnsr::A<double, 3> dt_Du_f;
  tnsr::A<double, 3> DuDu_f;
  tnsr::A<double, 3> dt_DuDu_f;
};

/*!
 * \brief Compute the acceleration inputs through second puncture order on a
 * Schwarzschild background in Cartesian Kerr-Schild coordinates.
 *
 * Extends the local closure of Sec. V.B of \cite Wittek:2024gxn to quadratic
 * regular fields. Third and higher spatial derivatives of the supplied
 * `psi` and `dt_psi` expansions are truncated, but their Hessians are retained.
 * The stationary-background wave equation and its derivatives determine
 * mixed/time derivatives through fourth differential order. These are local
 * closure estimates, not independently measured higher field derivatives.
 *
 * Automatic differentiation of the normalized equations of motion then
 * supplies the first three total time derivatives of the self-acceleration.
 * Covariant differentiation gives `Du_f`, `DuDu_f` and their time derivatives.
 * The differentiated expressions include the turn-on function and, when
 * `evolve_mass` is true, \f$\mu=\mu_0-q\Psi^R\f$. No numerical time
 * differencing or independent evolution of higher force derivatives is used.
 * Requires the existing SpECTRE option `SPECTRE_AUTODIFF=ON`.
 */
SelfForceDerivatives second_order_self_force_derivatives(
    const RegularFieldDerivatives& field, const tnsr::I<double, 3>& position,
    const tnsr::I<double, 3>& velocity, double charge, double rest_mass,
    bool evolve_mass, double black_hole_mass, double time, double turn_on_time,
    double turn_on_interval);

}  // namespace CurvedScalarWave::Worldtube
