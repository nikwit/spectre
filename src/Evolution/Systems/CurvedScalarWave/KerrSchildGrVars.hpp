// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
template <typename X, typename Symm, typename IndexList>
class Tensor;
/// \endcond

namespace CurvedScalarWave {

/*!
 * \brief Computes the background spacetime quantities of the
 * `CurvedScalarWave::System` for a Kerr-Schild spacetime with zero spin and
 * zero velocity in a single pass.
 *
 * This is a performance optimization of
 * `CurvedScalarWave::Actions::CalculateGrVars`, which evaluates the analytic
 * background on the moving grid at every time substep. The generic
 * `gr::Solutions::KerrSchild` solution computes the requested quantities
 * through a graph of cached intermediates, each evaluated in a separate pass
 * over the data, and copies the results out. For zero spin all requested
 * quantities have short closed forms, which this function evaluates together
 * in a single loop over the grid points, writing directly into the
 * preallocated output buffers.
 *
 * With \f$\phi = 2M/r\f$ and \f$n_i = (x_i - c_i)/r\f$, where \f$M\f$ is the
 * mass, \f$c_i\f$ the center, and \f$r\f$ the coordinate distance to the
 * center, the quantities are:
 *
 * \f{align}
 *   \alpha &= (1+\phi)^{-1/2}, \\
 *   \partial_i \alpha &= \frac{\phi}{2r} \alpha^3 n_i, \\
 *   \beta^i &= \frac{\phi}{1+\phi} n^i, \\
 *   \partial_i \beta^j &= \frac{\phi}{r(1+\phi)} \delta_i^j
 *     - \frac{\phi (2+\phi)}{r(1+\phi)^2} n_i n^j, \\
 *   \gamma_{ij} &= \delta_{ij} + \phi\, n_i n_j, \\
 *   \gamma^{ij} &= \delta^{ij} - \frac{\phi}{1+\phi} n^i n^j, \\
 *   \gamma^{jk}\Gamma^i_{jk} &= \frac{\phi(3+4\phi)}{2r(1+\phi)^2} n^i, \\
 *   K &= \frac{\phi}{r}\left(1 + \frac{3\phi}{2}\right) \alpha^3.
 * \f}
 */
void zero_spin_kerr_schild_gr_vars(
    gsl::not_null<Scalar<DataVector>*> lapse,
    gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> deriv_lapse,
    gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
    gsl::not_null<tnsr::iJ<DataVector, 3, Frame::Inertial>*> deriv_shift,
    gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*> spatial_metric,
    gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
        inverse_spatial_metric,
    gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        trace_spatial_christoffel_second_kind,
    gsl::not_null<Scalar<DataVector>*> trace_extrinsic_curvature,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords, double mass,
    const std::array<double, 3>& center);

}  // namespace CurvedScalarWave
