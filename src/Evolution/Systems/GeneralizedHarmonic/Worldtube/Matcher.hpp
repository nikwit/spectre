// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"

/// \cond
namespace ylm {
class Spherepack;
}  // namespace ylm
/// \endcond

namespace gh::Worldtube {

/// Diagnostics and result of one online fit.
struct FitResult {
  std::array<double, num_map_parameters> p{};
  double residual_initial = 0.;
  double residual_final = 0.;
  size_t iterations = 0;
};

/*!
 * \brief Fit the 13 first-order affine-map parameters from the evolved
 * fields on the excision sphere.
 *
 * Implements the findings-§15b prescription: the fit target is the
 * spherical-harmonic content (l <= `config.fit_l_max`) of the null-basis
 * gauge components {A, C, V} of the incoming characteristic field
 * \f$u^-_{ab} = \Pi_{ab} + n^k \Phi_{kab} - \gamma_2 g_{ab}\f$, with uniform
 * (absolute) weights. The model side is
 * `gh::Solutions::affine_map_model::evolved_variables` combined with the
 * *data* normal, frame, and \f$\gamma_2\f$, mirroring what the ghost
 * boundary condition applies. Pins: velocity kinematic to
 * `config.center_velocity`, trace strain to `config.trace_strain_pin`,
 * drives held at `pdot_estimate` (backward differences of the fit history).
 * Nine free parameters, solved by warm-started Gauss-Newton.
 *
 * All face tensors must be in the Spherepack collocation order of `ylm`
 * (the natural order of an excision-face slice on the Ylm-basis domain).
 */
FitResult fit_map_parameters(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const std::array<double, num_map_parameters>& p_start,
    const std::array<double, num_map_parameters>& pdot_estimate);
}  // namespace gh::Worldtube
