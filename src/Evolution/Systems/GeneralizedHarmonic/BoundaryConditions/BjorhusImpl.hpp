// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/ContainerHelpers.hpp"

/// \cond
template <size_t VolumeDim, typename Fr>
class MathFunction;
namespace gsl {
template <class T>
class not_null;
}  // namespace gsl
/// \endcond

namespace gh::BoundaryConditions {
/// \brief Detailed implementation of Bjorhus-type boundary corrections
namespace Bjorhus {
/*!
 * \brief Computes the expression needed to set boundary conditions on the time
 * derivative of the characteristic field \f$v^{g}_{ab}\f$
 *
 * \details In the Bjorhus scheme, the time derivatives of evolved variables are
 * characteristic projected. A constraint-preserving correction term is added
 * here to the resulting characteristic (time-derivative) field:
 *
 * \f{align}
 * \Delta \partial_t v^{g}_{ab} = \lambda_{g} n^i C_{iab}
 * \f}
 *
 * where \f$n^i\f$ is the local unit normal to the external boundary,
 * \f$C_{iab} = \partial_i g_{ab} - \Phi_{iab}\f$ is the three-index
 * constraint, and \f$\lambda_{g}\f$ is the characteristic speed of the field
 * \f$v^{g}_{ab}\f$.
 */
template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_psi(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> bc_dt_v_psi,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const std::array<DataType, 4>& char_speeds);

/*!
 * \brief Computes the expression needed to set boundary conditions on the time
 * derivative of the characteristic field \f$v^{0}_{iab}\f$
 *
 * \details In the Bjorhus scheme, the time derivatives of evolved variables are
 * characteristic projected. A constraint-preserving correction term is added
 * here to the resulting characteristic (time-derivative) field \f$v^0_{iab}\f$:
 *
 * \f{align}
 * \Delta \partial_t v^{0}_{iab} = \lambda_{0} n^j C_{jiab}
 * \f}
 *
 * where \f$n^i\f$ is the local unit normal to the external boundary,
 * \f$C_{ijab} = \partial_i\Phi_{jab} - \partial_j\Phi_{iab}\f$ is the
 * four-index constraint, and \f$\lambda_{0}\f$ is the characteristic speed of
 * the field \f$v^0_{iab}\f$.
 *
 * \note In 3D, only the non-zero and unique components of the four-index
 * constraint are stored, as \f$\hat{C}_{iab} = \epsilon_i^{jk} C_{jkab}\f$. In
 * 2D the input expected here is \f$\hat{C}_{0ab} = C_{01ab}, \hat{C}_{1ab} =
 * C_{10ab}\f$, and in 1D \f$\hat{C}_{0ab} = 0\f$.
 */
template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_zero(
    gsl::not_null<tnsr::iaa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_zero,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        four_index_constraint,
    const std::array<DataType, 4>& char_speeds);

/// @{
/*!
 * \brief Computes the expression needed to set boundary conditions on the time
 * derivative of the characteristic field \f$v^{-}_{ab}\f$
 *
 * \details In the Bjorhus scheme, the time derivatives of evolved variables are
 * characteristic projected. Constraint-preserving correction terms
 * \f$T^{\mathrm C}_{ab}\f$, Sommerfeld condition terms for gauge degrees of
 * freedom \f$T^{\mathrm G}_{ab}\f$, and terms constraining physical degrees of
 * freedom \f$T^{\mathrm P}_{ab}\f$ are added here to the resulting
 * characteristic (time-derivative) field:
 *
 * \f{align}
 * \Delta \partial_t v^{-}_{ab} = T^{\mathrm C}_{ab} + T^{\mathrm G}_{ab}
 *                              + T^{\mathrm P}_{ab}.
 * \f}
 *
 * These terms are given by Eq. (64) of \cite Lindblom2005qh :
 *
 * \f{eqnarray}
 * T^{\mathrm C}_{ab} &=& \frac{1}{2}
 *      \left(2 k^c k^d l_a l_b - k^c l_b P^d_a - k^c l_a P^d_b - k^d l_b P^c_a
 *            - k^d l_a P^c_b + P^{cd} P_{ab}\right) \partial_t v^{-}_{cd}\\
 *     &&+ \frac{1}{\sqrt{2}} \lambda_{-}c^{\hat{0}-}_c \left(l_a l_b k^c +
 *          P_{ab} l^c - P^c_b l_a - P^c_a l_b \right) \nonumber
 * \f}
 *
 * where \f$l^a (l_a)\f$ is the outgoing null vector (one-form), \f$k^a (k_a)\f$
 * is the incoming null vector (one-form), \f$P_{ab}, P^{ab}, P^a_b\f$ are
 * spacetime projection operators defined in `transverse_projection_operator()`,
 * \f$\partial_t v^{-}_{ab}\f$ is the characteristic projected time derivative
 * of evolved variables (corresponding to the \f$v^{-}\f$ field), and
 * \f$c^{\hat{0}\pm}_a\f$ are characteristic modes of the constraint evolution
 * system:
 *
 * \f{align}\nonumber
 * c^{\hat{0}\pm}_a = F_a \mp n^k C_{ka},
 * \f}
 *
 * where \f$F_a\f$ is the generalized-harmonic (GH) F constraint [Eq. (43) of
 * \cite Lindblom2005qh], \f$C_{ka}\f$ is the GH 2-index constraint [Eq. (44)
 * of \cite Lindblom2005qh], and \f$n^k\f$ is the unit spatial normal to the
 * outer boundary. Boundary correction terms that prevent strong reflections of
 * gauge perturbations are given by Eq. (25) of \cite Rinne2007ui :
 *
 * \f{align}
 * T^{\mathrm G}_{ab} =
 *     \left(k_a P^c_b l^d + k_b P^c_a l^d -
 *          \left(k_a l_b k^c l^d + k_b l_a k_c l^d + k_a k_b l^c l^d
 *          \right) \right) \left(\gamma_2 - \frac{1}{r}
 *                          \right) \partial_t v^{g}_{cd}
 * \f}
 *
 * where \f$r\f$ is the radial coordinate at the outer boundary, which is
 * assumed to be spherical, \f$\gamma_2\f$ is a GH constraint damping parameter,
 * and \f$\partial_t v^{g}_{ab}\f$ is the characteristic projected time
 * derivative of evolved variables (corresponding to the \f$v^{g}\f$ field).
 * Finally, we constrain physical degrees of freedom using corrections from
 * Eq. (68) of \cite Lindblom2005qh :
 *
 * \f{align}
 * T^{\mathrm P}_{ab} = \left( P^c_a P^d_b - \frac{1}{2} P_{ab} P^{cd} \right)
 *     \left(\partial_t v^{-}_{cd} + \lambda_{-}
 *           \left(U^{3-}_{cd} - \gamma_2 n^i C_{icd}\right)\right),
 * \f}
 *
 * where \f$C_{icd}\f$ is the GH 3-index constraint [c.f. Eq. (26) of
 * \cite Lindblom2005qh, see also `three_index_constraint()`], and
 *
 * \f{align}\nonumber
 * U^{3-}_{ab} = 2 P^i_a P^j_b U^{8-}_{ij}
 * \f}
 *
 * is the inward propagating characteristic mode of the Weyl tensor evolution
 * \f$U^{8-}_{ab}\f$ [c.f. Eq. 75 of \cite Kidder2004rw ] projected onto the
 * outer boundary using the spatial-spacetime projection operators \f$P^i_a\f$.
 * Note that (A) the covariant derivative of extrinsic curvature needed to
 * get the Weyl propagating modes is calculated substituting the evolved
 * variable \f$\Phi_{iab}\f$ for spatial derivatives of the spacetime metric;
 * and (B) the spatial Ricci tensor used in the same calculation is also
 * calculated using the same substituion, and also includes corrections
 * proportional to the GH 4-index constraint \f$C_{ijab}\f$ [c.f. Eq. (45) of
 * \cite Lindblom2005qh , see also `four_index_constraint()`]:
 *
 * \f{align}\nonumber
 * R_{ij} \rightarrow R_{ij} + C_{(iklj)} + \frac{1}{2} n^k q^a C_{(ikj)a},
 * \f}
 *
 * where \f$q^a\f$ is the future-directed spacetime normal vector.
 */
template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_minus(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const std::array<DataType, 4>& char_speeds);

template <size_t VolumeDim, typename DataType>
void constraint_preserving_gauge_corrections_dt_v_minus(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const std::array<DataType, 4>& char_speeds);

template <size_t VolumeDim, typename DataType>
void constraint_preserving_gauge_physical_corrections_dt_v_minus(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    double time,
    const tnsr::i<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::II<DataType, VolumeDim, Frame::Inertial>&
        inverse_spatial_metric,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>& extrinsic_curvature,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& spacetime_metric,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>&
        inverse_spacetime_metric,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& phi,
    const tnsr::ijaa<DataType, VolumeDim, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& d_pi,
    const std::array<DataType, 4>& char_speeds,
    const MathFunction<1, Frame::Inertial>* incoming_wave_profile = nullptr);
/// @}

/*!
 * \brief The quantities on an external face that every Bjorhus-type boundary
 * condition of the generalized harmonic system assembles from the interior
 * data before adding its own corrections.
 *
 * \details Fill with `compute_intermediate_variables()`. The time derivatives
 * are in the inertial frame: when the face moves, the advective mesh-velocity
 * term is subtracted from the logical time derivatives supplied by the DG
 * scheme. The characteristic speeds already account for the mesh velocity, see
 * `gh::characteristic_speeds()`.
 *
 * A boundary condition then computes its corrections to the four
 * characteristic time derivatives and calls
 * `project_corrections_onto_evolved_variables()` to obtain the corrections to
 * the evolved variables. See `ConstraintPreservingBjorhus` and
 * `WorldtubeTypeD` for the two users.
 */
template <size_t Dim>
struct IntermediateVariables {
  explicit IntermediateVariables(size_t num_points);

  tnsr::a<DataVector, Dim, Frame::Inertial> spacetime_unit_normal_one_form;
  tnsr::I<DataVector, Dim, Frame::Inertial> unit_interface_normal_vector;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> four_index_constraint;
  tnsr::II<DataVector, Dim, Frame::Inertial> inverse_spatial_metric;
  tnsr::ii<DataVector, Dim, Frame::Inertial> extrinsic_curvature;
  tnsr::a<DataVector, Dim, Frame::Inertial> incoming_null_one_form;
  tnsr::a<DataVector, Dim, Frame::Inertial> outgoing_null_one_form;
  tnsr::A<DataVector, Dim, Frame::Inertial> incoming_null_vector;
  tnsr::A<DataVector, Dim, Frame::Inertial> outgoing_null_vector;
  tnsr::aa<DataVector, Dim, Frame::Inertial> projection_ab;
  tnsr::Ab<DataVector, Dim, Frame::Inertial> projection_Ab;
  tnsr::AA<DataVector, Dim, Frame::Inertial> projection_AB;
  /// Inertial-frame time derivatives of the evolved variables
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_spacetime_metric;
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_pi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> dt_phi;
  /// Characteristic projections of the inertial-frame time derivatives
  tnsr::aa<DataVector, Dim, Frame::Inertial> char_projected_rhs_dt_v_psi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> char_projected_rhs_dt_v_zero;
  tnsr::aa<DataVector, Dim, Frame::Inertial> char_projected_rhs_dt_v_plus;
  tnsr::aa<DataVector, Dim, Frame::Inertial> char_projected_rhs_dt_v_minus;
  /// \f$c^{\hat 0 \pm}_a = F_a \mp n^k C_{ka}\f$
  tnsr::a<DataVector, Dim, Frame::Inertial> constraint_char_zero_plus;
  tnsr::a<DataVector, Dim, Frame::Inertial> constraint_char_zero_minus;
  std::array<DataVector, 4> char_speeds;
};

/// \brief Fill an `IntermediateVariables` from the data on the face.
///
/// The arguments are the ones a boundary condition receives in its
/// `dg_time_derivative` for the tags of `ConstraintPreservingBjorhus`.
template <size_t Dim>
void compute_intermediate_variables(
    gsl::not_null<IntermediateVariables<Dim>*> vars,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
    const tnsr::AA<DataVector, Dim, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, Dim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, Dim, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, Dim, Frame::Inertial>&
        spacetime_deriv_gauge_source,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>&
        logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& logical_dt_phi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi);

/// \brief The smallest characteristic speed on the face. When it is
/// non-negative no characteristic field is incoming anywhere on the face and a
/// Bjorhus-type boundary condition has nothing to correct.
double min_characteristic_speed(const std::array<DataVector, 4>& char_speeds);

/*!
 * \brief Turn corrections to the characteristic time derivatives into
 * corrections to the time derivatives of the evolved variables.
 *
 * \details The corrections are first set to zero at every grid point where
 * the corresponding characteristic speed is positive, since only incoming
 * fields receive boundary conditions. They are then transformed with
 * `gh::evolved_fields_from_characteristic_fields()`, i.e. inverted from the
 * characteristic projection, c.f. Eqs. (63) - (65) of \cite Lindblom2005qh.
 * The `bc_dt_v_*` arguments are modified in place by the zeroing.
 */
template <size_t Dim>
void project_corrections_onto_evolved_variables(
    gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        dt_spacetime_metric_correction,
    gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> dt_pi_correction,
    gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        dt_phi_correction,
    gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> bc_dt_v_psi,
    gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> bc_dt_v_zero,
    gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> bc_dt_v_plus,
    gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> bc_dt_v_minus,
    const std::array<DataVector, 4>& char_speeds,
    const Scalar<DataVector>& gamma2,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector);

namespace detail {
/// @{
/*!
 * \brief Add \f$c\,X^{\rm sector}_{ab}\f$ to `result`, where
 * \f$X^{\rm sector}\f$ is the projection of the symmetric source
 * \f$X_{cd}\f$ onto one of the three characteristic sectors of \f$u^-\f$ and
 * \f$c\f$ is a scalar coefficient.
 *
 * \details The three sectors are the constraint-preserving, physical and
 * gauge parts of Eq. (64) of \cite Lindblom2005qh, i.e. the tensor structures
 * multiplying \f$\partial_t v^-_{cd}\f$ in \f$T^{\mathrm C}_{ab}\f$,
 * \f$T^{\mathrm P}_{ab}\f$ and \f$T^{\mathrm G}_{ab}\f$ (see
 * `constraint_preserving_gauge_physical_corrections_dt_v_minus()`). The gauge
 * sector is defined with the sign that makes the three projections partition
 * the identity on symmetric tensors: they are idempotent, annihilate one
 * another and sum to the identity. This is what allows a boundary condition
 * to impose the sectors independently, see `WorldtubeTypeD`.
 *
 * \note `add_gauge_sommerfeld_terms_to_dt_v_minus()` applies the negative of
 * the gauge projection, so that a positive coefficient there reads as the
 * damping rate of \f$\partial_t u^-|_{\rm gauge}\f$.
 */
template <size_t VolumeDim, typename DataType>
void add_constraint_sector_projection(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source);

template <size_t VolumeDim, typename DataType>
void add_physical_sector_projection(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source);

template <size_t VolumeDim, typename DataType>
void add_gauge_sector_projection(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source);
/// @}

/*!
 * \brief Add the gauge-sector correction \f$-c\,X^{\rm gauge}_{ab}\f$ to
 * `bc_dt_v_minus` for a pointwise coefficient \f$c\f$ and a symmetric source
 * \f$X_{cd}\f$.
 *
 * \details With \f$c = \gamma_2 - 1/r\f$ and \f$X = \partial_t v^g\f$ this is
 * the Sommerfeld gauge condition \f$T^{\mathrm G}_{ab}\f$ of
 * `constraint_preserving_gauge_physical_corrections_dt_v_minus()`;
 * `add_gauge_sommerfeld_terms_to_dt_v_minus()` forms that coefficient with the
 * radius measured from the origin of the inertial coordinates. A boundary
 * condition whose boundary is not centered at the origin forms the coefficient
 * itself and calls this function.
 */
template <size_t VolumeDim, typename DataType>
void add_gauge_sector_terms_to_dt_v_minus(
    gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source);

template <size_t VolumeDim, typename DataType>
void add_gauge_sommerfeld_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi);

template <size_t VolumeDim, typename DataType>
void add_constraint_dependent_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const std::array<DataType, 4>& char_speeds);

/*!
 * \brief Add the physical correction \f$T^{\mathrm P}_{ab}\f$ of
 * `constraint_preserving_gauge_physical_corrections_dt_v_minus()` to
 * `bc_dt_v_minus`.
 *
 * \details The incoming Weyl mode \f$U^{8-}_{ij}\f$ is computed from the
 * interior data on the face. By default the correction drives it to zero (no
 * incoming gravitational waves), optionally plus the injected wave given by
 * `incoming_wave_profile`. If `incoming_weyl_propagating_minus` is not null,
 * it is subtracted from the interior \f$U^{8-}_{ij}\f$ before the correction
 * is formed, so that the incoming field is driven towards the supplied value
 * instead. It must be given in the normalization and projection of
 * `gr::weyl_propagating()` with sign \f$-1\f$, i.e. spatial, symmetric, and
 * transverse to the face normal. This is how a worldtube boundary condition
 * supplies the incoming radiation of the excised region.
 */
template <size_t VolumeDim, typename DataType>
void add_physical_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::i<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::II<DataType, VolumeDim, Frame::Inertial>&
        inverse_spatial_metric,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>& extrinsic_curvature,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& spacetime_metric,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>&
        inverse_spacetime_metric,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& phi,
    const tnsr::ijaa<DataType, VolumeDim, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& d_pi,
    const std::array<DataType, 4>& char_speeds,
    double time = std::numeric_limits<double>::signaling_NaN(),
    const MathFunction<1, Frame::Inertial>* incoming_wave_profile = nullptr,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>*
        incoming_weyl_propagating_minus = nullptr);
}  // namespace detail
}  // namespace Bjorhus
}  // namespace gh::BoundaryConditions
