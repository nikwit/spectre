// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace CurvedScalarWave::BoundaryCorrections {
/*!
 * \brief Computes the upwind multipenalty boundary correction for scalar wave
 * in curved spacetime.
 *
 * \details This implements the upwind multipenalty boundary correction term,
 * given by \cite Teukolsky2015ega :
 *
 * \f[
 * G = S\left(\Lambda^+ S^{-1} U^{\rm int}
 *           + \Lambda^- S^{-1} U^{\rm ext}\right)
 * = S \left(\Lambda^+ \hat{U}^{\rm int}
 *           + \Lambda^- \hat{U}^{\rm ext}\right),
 * \f]
 *
 * where
 *
 *  - \f$G\f$ is the numerical upwind flux dotted with the interface normal;
 *  - \f$U\f$ is a vector of all evolved variables;
 *  - \f$S\f$ is a matrix whose columns are the eigenvectors of the
 *       characteristic matrix for the evolution system. It maps the
 *       evolved variables to characteristic variables \f$\hat{U}\f$, s.t.
 *       \f$\hat{U} := S^{-1}\cdot U\f$; and
 *  - \f$\Lambda^\pm\f$ are diagonal matrices containing
 *       positive / negative eigenvalues of the same matrix as its elements.
 *
 * The superscripts \f${\rm int}\f$ and \f${\rm ext}\f$ on \f$U\f$ indicate
 * that the corresponding set of variables at the element interface have been
 * taken from the _interior_ or _exterior_ of the element. Exterior of the
 * element is naturally the interior of its neighboring element. Therefore,
 * \f$\hat{U}^{\rm int}:= S^{-1}U^{\rm int}\f$ are the characteristic variables
 * at the element interface computed using evolved variables from the interior
 * of the element, and \f$\hat{U}^{\rm ext}= S^{-1}U^{\rm ext}\f$ are the
 * same computed from evolved variables taken from the element exterior, i.e.
 * the neighboring element.
 * The sign of characteristic speeds indicates the direction of propagation of
 * the corresponding characteristic field with respect to the interface normal
 * that the field has been computed along (with negative speeds indicating
 * incoming characteristics, and with positive speeds indicating outgoing
 * characteristics). Therefore, \f$\Lambda^+\f$ contains characterstic speeds
 * for variables that are outgoing from the element at the interface, and
 * \f$\Lambda^-\f$ contains characteristic speeds for variables that are
 * incoming to the element at the interface. An ambiguity naturally arises
 * as to which set of evolved variables (''interior'' or ''exterior'') to use
 * when computing these speeds. We compute both and use their average, i.e.
 * \f$\lambda^{\rm avg} = (\lambda^{\rm int} + \lambda^{\rm ext})/2\f$, to
 * populate \f$\Lambda^\pm\f$.
 *
 *
 * This function computes the upwind flux as follows:
 *  -# Computes internal and external characteristic variables and speeds using
 * evolved variables from
 *     both the element interior and exterior.
 *  -# Computes the average characteristic speeds and constructs
 *      \f$\Lambda^{\pm}\f$.
 *  -# Computes the upwind flux as a weigted sum of external and internal
 *     characteristic variables
 *
 * \f[
 * G = S\left(w_{\rm ext} \hat{U}^{\rm ext}
 *       + w_{\rm int} \hat{U}^{\rm int}\right),
 * \f]
 *
 * with weights \f$w_{\rm ext} = \Theta(-\Lambda)\cdot\Lambda\f$, and
 * \f$w_{\rm int} = \Theta(\Lambda)\cdot\Lambda\f$, where \f$\Theta\f$ is the
 * step function centered at zero, \f$\Lambda = \Lambda^+ + \Lambda^-\f$, and
 * the dot operator \f$(\cdot)\f$ indicates an element-wise product.
 *
 * \warning With the averaging of characteristic speeds, this flux does not
 * satisfy the generalized Rankine-Hugoniot conditions
 * \f${G}^{\rm int}(\hat{U}^{\rm int}, \hat{U}^{\rm ext})
 * = - {G}^{\rm ext}(\hat{U}^{\rm ext}, \hat{U}^{\rm int})\f$,
 * which enforces that the flux leaving the element is equal to the flux
 * entering the neighboring element, and vice versa. This condition is
 * important for a well-balanced scheme, and so please use this flux with
 * caution.
 */
template <size_t Dim>
class UpwindPenalty final : public BoundaryCorrection<Dim> {
 private:
  struct CharSpeedsTensor : db::SimpleTag {
    using type = tnsr::a<DataVector, 3, Frame::Inertial>;
  };

 public:
  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "Computes the UpwindPenalty boundary correction term for the scalar wave "
      "system."};

  UpwindPenalty() = default;
  UpwindPenalty(const UpwindPenalty&) = default;
  UpwindPenalty& operator=(const UpwindPenalty&) = default;
  UpwindPenalty(UpwindPenalty&&) = default;
  UpwindPenalty& operator=(UpwindPenalty&&) = default;
  ~UpwindPenalty() override = default;

  /// \cond
  explicit UpwindPenalty(CkMigrateMessage* msg) noexcept;
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(UpwindPenalty);  // NOLINT
  /// \endcond
  void pup(PUP::er& p) override;  // NOLINT

  std::unique_ptr<BoundaryCorrection<Dim>> get_clone() const noexcept override;

  using dg_package_field_tags =
      tmpl::list<Tags::VPsi, Tags::VZero<Dim>, Tags::VPlus, Tags::VMinus,
                 Tags::ConstraintGamma2,
                 ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<
                     Dim, Frame::Inertial>>,
                 CharSpeedsTensor>;
  using dg_package_data_temporary_tags = tmpl::list<
      gr::Tags::Lapse<DataVector>,
      gr::Tags::Shift<Dim, Frame::Inertial, DataVector>,
      gr::Tags::InverseSpatialMetric<Dim, Frame::Inertial, DataVector>,
      Tags::ConstraintGamma1, Tags::ConstraintGamma2>;
  using dg_package_data_volume_tags = tmpl::list<>;

  double dg_package_data(
      gsl::not_null<Scalar<DataVector>*> packaged_v_psi,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*> packaged_v_zero,
      gsl::not_null<Scalar<DataVector>*> packaged_v_plus,
      gsl::not_null<Scalar<DataVector>*> packaged_v_minus,
      gsl::not_null<Scalar<DataVector>*> packaged_gamma2,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          packaged_interface_unit_normal,
      gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
          packaged_char_speeds,

      const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,
      const Scalar<DataVector>& psi,

      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
      const tnsr::II<DataVector, Dim, Frame::Inertial>& inverse_spatial_metric,
      const Scalar<DataVector>& constraint_gamma1,
      const Scalar<DataVector>& constraint_gamma2,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& interface_unit_normal,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*mesh_velocity*/,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity)
      const noexcept;

  void dg_boundary_terms(
      gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          phi_boundary_correction,
      gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,

      const Scalar<DataVector>& v_psi_int,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& v_zero_int,
      const Scalar<DataVector>& v_plus_int,
      const Scalar<DataVector>& v_minus_int,
      const Scalar<DataVector>& gamma2_int,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          interface_unit_normal_int,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,

      const Scalar<DataVector>& v_psi_ext,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& v_zero_ext,
      const Scalar<DataVector>& v_plus_ext,
      const Scalar<DataVector>& v_minus_ext,
      const Scalar<DataVector>& gamma2_ext,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          interface_unit_normal_ext,
      const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
      dg::Formulation /*dg_formulation*/) const noexcept;
};
}  // namespace CurvedScalarWave::BoundaryCorrections
