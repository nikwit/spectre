// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <ostream>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/CoulombDecode.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"

/// \cond
namespace gh::worldtube {
template <size_t Dim>
struct KretschmannFaceData;
}  // namespace gh::worldtube
/// \endcond

namespace gh::worldtube {
/*!
 * \brief Which model supplies the radiation entering the domain through the
 * physical sector of the worldtube boundary condition.
 */
enum class PhysicalModel {
  /// No model: the physical Bjorhus term forbids incoming radiation, as at an
  /// outer boundary.
  None,
  /// Order zero of the curvature matching: the excised region holds a type-D
  /// (Kerr) hole. The Kinnersley Coulomb scalar \f$\Psi_2^K = -3J/I\f$ and
  /// the null rotations aligning the NR tetrad with the principal null
  /// directions are read off the face curvature pointwise, and the incoming
  /// mode is \f$\Psi_0 = 6 b^2 \Psi_2^K\f$ (eq. `psi0-leading` of the
  /// worldtube NP-matching note). Needs no mass, boost or radius.
  TypeD,
  /// Order two of the curvature matching: the type-D hole carries a
  /// quadrupolar tide. The frame is registered as for `TypeD`, boosted to
  /// the invariant rest frame fixed by the spacetime gradient of the
  /// Kretschmann scalar (`gr::np::invariant_rapidity()`), the five complex
  /// tidal moments are fitted to the pulled-back \f$\Psi_4\f$ over the face,
  /// and the incoming mode is the tide's \f$\Psi_0\f$ pushed forward to the
  /// NR tetrad (`gr::np::evaluate_second_order()`). Needs the mass of the
  /// hole and the `KretschmannFaceData` of the element.
  Quadrupole,
  /// Order two with the tide read from the Coulomb channel instead of the
  /// leaving mode: the areal radius of every face point follows from the
  /// normal derivative of the Coulomb scalar with the kinematics of the
  /// type-D background treated exactly, and the moments are the fit of the
  /// Coulomb excess over \f$-M/r^3\f$
  /// (`gr::np::decode_tidal_moments_from_coulomb()`). The radiative modes
  /// enter only through the invariants, quadratically, so the condition does
  /// not feed back on itself through \f$\Psi_4\f$. Needs the same inputs as
  /// `Quadrupole` and the excision outside about \f$2.4M\f$.
  QuadrupoleCoulomb
};

/// Whether the model is one of the order-two models, which need the mass
/// and the Kretschmann face data
bool is_order_two(PhysicalModel model);

PhysicalModel convert_physical_model_from_yaml(const Options::Option& options);

std::ostream& operator<<(std::ostream& os, PhysicalModel model);

/// The result of `evaluate_matching()`
struct MatchingEvaluation {
  /// Rows: the adapted triad \f$(s, \hat\theta, \hat\phi)\f$ in the
  /// Cholesky-orthonormal frame of the spatial metric
  gr::np::RealMatrix adapted_rotation;
  /// The Weyl scalars of the face curvature in the adapted tetrad, with
  /// \f$s\f$ along the outward normal of the domain
  gr::np::WeylScalars psi;
  Scalar<ComplexDataVector> coulomb;
  gr::np::TypeDRotation type_d_rotation;
  /// The model's \f$\Psi_0\f$, the mode entering the domain
  Scalar<ComplexDataVector> psi0_target;
  /// \f$U^{8-}_{ij} = \tfrac12 w^-_{ij}(\Psi_0)\f$ in covariant coordinate
  /// components, the normalization of `gr::weyl_propagating()` with sign
  /// \f$-1\f$
  tnsr::ii<DataVector, 3, Frame::Inertial> incoming_mode;
  /// Set for the order-two models only
  std::optional<Scalar<DataVector>> rapidity;
  std::optional<gr::np::FrameRegistration> registration;
  std::optional<gr::np::SecondOrderEvaluation> second_order;
  /// Set for `PhysicalModel::QuadrupoleCoulomb` without imposed moments
  std::optional<gr::np::CoulombDecode> coulomb_decode;
};

/*!
 * \brief The derivative of the real part of the Coulomb scalar along the
 * sphere normal from the Kretschmann data on the face.
 *
 * \details We have \f$K = 16\,\mathrm{Re}\, I = 48\,\mathrm{Re}\,\Psi_2^2\f$
 * up to terms quadratic in the tide, so
 * \f$\partial_s K = 96\,\mathrm{Re}(\Psi_2 \partial_s \Psi_2)\f$ and, for a
 * Coulomb scalar with a small imaginary part (a slowly spinning hole),
 * \f$\partial_s \mathrm{Re}\,\Psi_2 = \partial_s K /
 * (96\,\mathrm{Re}\,\Psi_2)\f$. `d_kretschmann` is the coordinate gradient of
 * \f$K\f$ and `unit_normal_vector` the unit normal vector \f$s^i\f$ of the
 * face.
 */
Scalar<DataVector> normal_derivative_of_coulomb(
    const Scalar<ComplexDataVector>& coulomb,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_kretschmann,
    const tnsr::I<DataVector, 3, Frame::Inertial>& unit_normal_vector);

/*!
 * \brief The curvature matching on the excision face of the worldtube.
 *
 * \details The worldtube-adapted tetrad is built with \f$s\f$ along
 * `unit_normal_covector`, the outward normal of the domain, which at the
 * excision points into the hole. In that tetrad the mode entering the domain
 * is \f$\Psi_0\f$ and
 * \f$U^{8-}_{ij} = \tfrac12 w^-_{ij}
 *  = \bar\Psi_0 m_i m_j + \Psi_0 \bar m_i \bar m_j\f$
 * (checked against `gr::weyl_propagating()` in the unit tests). The
 * `model` selects the value of \f$\Psi_0\f$, see `PhysicalModel`. For
 * the order-two models the `mass`, the `lapse`, the `shift` and the
 * `face_data` (with the point count of the face) are required; they are
 * ignored otherwise. With `imposed_moments` the tidal fit is skipped and the
 * target is built from these moments, e.g. the relaxed moments of
 * `KretschmannFaceData::filtered_moments` the boundary condition imposes.
 * `PhysicalModel::None` is not a model and is rejected.
 *
 * The construction consumes the interior \f$\Psi_0\f$ through the invariants
 * and the type-D solve, so an evolution imposing the result is the fixed
 * point iteration of the note's circularity discussion.
 */
MatchingEvaluation evaluate_matching(
    PhysicalModel model, std::optional<double> mass,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial>& unit_normal_covector,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const KretschmannFaceData<3>* face_data,
    const std::optional<gr::np::TidalMoments>& imposed_moments = std::nullopt);
}  // namespace gh::worldtube

template <>
struct Options::create_from_yaml<gh::worldtube::PhysicalModel> {
  template <typename Metavariables>
  static gh::worldtube::PhysicalModel create(const Options::Option& options) {
    return gh::worldtube::convert_physical_model_from_yaml(options);
  }
};
