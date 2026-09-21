// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"

/// \cond
template <size_t Dim>
class Element;
template <size_t Dim>
class ExcisionSphere;
template <size_t Dim>
class Mesh;
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gh::worldtube {
/*!
 * \brief The Kretschmann scalar, its spacetime gradient and the quadrature
 * weights on the excision face of an element, for the second-order worldtube
 * matching.
 *
 * \details The invariant rest frame of the curvature matching
 * (`gr::np::invariant_rapidity()`) is the frame in which the Kretschmann
 * scalar \f$K\f$ is momentarily stationary, so the boundary condition needs
 * \f$D_i K\f$ and \f$\partial_t K\f$ on the face. Neither is available from
 * the face data a boundary condition receives: \f$D_i K\f$ contains the
 * normal derivative, i.e. third derivatives of the metric, so \f$K\f$ is
 * differentiated in the element volume and sliced
 * (`update_kretschmann_face_data()`);
 * \f$\partial_t K\f$ at fixed inertial coordinates is a backward difference
 * of the face values between successive evaluations, corrected for the mesh
 * velocity, \f$\partial_t K = dK/dt|_{\rm grid} - v^i D_i K\f$. The first
 * evaluation, and any evaluation at the time of the previous one, has no
 * difference to take; it uses the grid-comoving estimate
 * \f$\partial_t K = -v^i D_i K\f$, or keeps the previous estimate.
 *
 * `direction` is the direction of the excision face and is unset when the
 * element abuts no excision sphere, in which case the remaining members are
 * empty. `quadrature_weights` are the weights of the face collocation grid
 * normalized to unit sum; on a spherical-harmonic face they are the exact
 * weights of the unit-sphere measure, so the least-squares fit of the
 * matching becomes a modal fit.
 */
template <size_t Dim>
struct KretschmannFaceData {
  std::optional<Direction<Dim>> direction{};
  double time{std::numeric_limits<double>::signaling_NaN()};
  Scalar<DataVector> kretschmann{};
  tnsr::i<DataVector, Dim, Frame::Inertial> d_kretschmann{};
  Scalar<DataVector> dt_kretschmann{};
  DataVector quadrature_weights{};
  /// Face values and time of the previous evaluation, for the backward
  /// difference
  double previous_time{std::numeric_limits<double>::signaling_NaN()};
  Scalar<DataVector> previous_kretschmann{};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

template <size_t Dim>
bool operator==(const KretschmannFaceData<Dim>& lhs,
                const KretschmannFaceData<Dim>& rhs);
template <size_t Dim>
bool operator!=(const KretschmannFaceData<Dim>& lhs,
                const KretschmannFaceData<Dim>& rhs);

/// \brief The direction of the excision face of the element, if it abuts an
/// excision sphere.
template <size_t Dim>
std::optional<Direction<Dim>> excision_face_direction(
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const Element<Dim>& element);

/*!
 * \brief Quadrature weights of a face collocation grid, normalized to unit
 * sum.
 *
 * \details The weights are the tensor product of the one-dimensional weights.
 * A `Spectral::Basis::SphericalHarmonic` dimension with Gauss quadrature
 * carries the Gauss-Legendre weights in \f$\cos\theta\f$ and an equiangular
 * one the uniform weights \f$2\pi/n_\phi\f$, so that the product is the
 * unit-sphere area element on a spherical-harmonic face.
 */
template <size_t Dim>
DataVector face_quadrature_weights(const Mesh<Dim - 1>& face_mesh);

/*!
 * \brief Update the Kretschmann face data of an element from its evolved
 * variables at time `time`.
 *
 * \details See `KretschmannFaceData`. Differentiates \f$\Pi\f$ and \f$\Phi\f$
 * in the volume, builds the Weyl curvature and \f$K\f$ there
 * (`kretschmann_scalar()`), differentiates \f$K\f$ and slices \f$K\f$ and
 * \f$D_i K\f$ to the excision face. The `mesh_velocity` is the volume mesh
 * velocity, if the domain has one. Only implemented for `Dim == 3`; in other
 * dimensions only the direction and the weights are set.
 */
template <size_t Dim>
void update_kretschmann_face_data(
    gsl::not_null<KretschmannFaceData<Dim>*> data,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Element<Dim>& element,
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        mesh_velocity,
    double time);

namespace Tags {
/// The `gh::worldtube::KretschmannFaceData` of the element
template <size_t Dim>
struct KretschmannFaceData : db::SimpleTag {
  using type = worldtube::KretschmannFaceData<Dim>;
};
}  // namespace Tags
}  // namespace gh::worldtube
