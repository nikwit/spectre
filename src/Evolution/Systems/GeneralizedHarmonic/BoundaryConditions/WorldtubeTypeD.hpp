// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <unordered_map>
#include <variant>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/Context.hpp"
#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t Dim>
class Domain;
template <size_t Dim>
class Element;
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace gh::BoundaryConditions::detail {
/// How one characteristic sector of the worldtube boundary condition is
/// imposed.
enum class SectorImposition {
  /// Through a Bjorhus time-derivative correction on this sector.
  Bjorhus,
  /// Not at all: the sector is frozen, \f$\partial_t u^-|_{\rm sector} = 0\f$.
  /// A diagnostic control that deliberately imposes no condition.
  Frozen,
  /// Gauge sector only. Bayliss-Turkel \f$L=0\f$ condition in the form that
  /// annihilates an *ingoing* wave \f$g(t+r)/r\f$, i.e. one the hole absorbs:
  /// \f$(\partial_t - \partial_r - 1/r)f = 0\f$. The radius is measured from
  /// the center of the excision sphere the boundary belongs to.
  SommerfeldAbsorbing,
  /// Gauge sector only. The outer-boundary form of the same condition, which
  /// annihilates an *outgoing* wave \f$g(t-r)/r\f$:
  /// \f$(\partial_t + \partial_r + 1/r)f = 0\f$. Provided so the two signs can
  /// be compared; at an inner boundary `SommerfeldAbsorbing` is the physically
  /// motivated one.
  SommerfeldOutgoing
};

/// Whether an imposition is one of the two Sommerfeld gauge conditions.
bool is_sommerfeld(SectorImposition imposition);

/// Sign of the \f$1/r\f$ term as it enters the coefficient
/// \f$\gamma_2 - c/r\f$: \f$-1\f$ absorbs into the hole, \f$+1\f$ is the
/// outer-boundary form.
double sommerfeld_one_over_r_sign(SectorImposition imposition);

SectorImposition convert_sector_imposition_from_yaml(
    const Options::Option& options);

std::ostream& operator<<(std::ostream& os, SectorImposition imposition);

/*!
 * \brief Which model supplies the radiation entering the domain through the
 * physical sector.
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
  TypeD
};

PhysicalModel convert_physical_model_from_yaml(const Options::Option& options);

std::ostream& operator<<(std::ostream& os, PhysicalModel model);

/// \brief Per-field imposition of the constraint-preserving sector.
///
/// The sector bundles three structurally different corrections: the
/// three-index-constraint term on \f$\partial_t v_\psi\f$, the
/// four-index-constraint term on \f$\partial_t v_0\f$, and the constraint
/// projection of \f$\partial_t u^-\f$. Imposing them separately lets a
/// constraint influx be attributed to one term rather than to the sector as a
/// whole.
struct PerFieldConstraintSectors {
  struct VPsi {
    using type = SectorImposition;
    static constexpr Options::String help{
        "Imposition of the three-index-constraint term on dt v_psi. NOTE that "
        "Frozen here raises the three-index constraint by construction, "
        "because v_psi = psi_ab and C_iab = d_i psi_ab - Phi_iab: holding "
        "psi_ab fixed while Phi_iab evolves makes C_iab grow at the face, and "
        "the Bjorhus term being replaced is the term that damps it there."};
  };
  struct VZero {
    using type = SectorImposition;
    static constexpr Options::String help{
        "Imposition of the four-index-constraint term on dt v_zero."};
  };
  struct VMinus {
    using type = SectorImposition;
    static constexpr Options::String help{
        "Imposition of the constraint projection of v_minus."};
  };

  using options = tmpl::list<VPsi, VZero, VMinus>;
  static constexpr Options::String help{
      "Impose the three terms of the constraint-preserving sector separately. "
      "Diagnostic form of ConstraintPreservingSector."};

  PerFieldConstraintSectors() = default;
  PerFieldConstraintSectors(SectorImposition v_psi_in,
                            SectorImposition v_zero_in,
                            SectorImposition v_minus_in);

  SectorImposition v_psi{SectorImposition::Bjorhus};
  SectorImposition v_zero{SectorImposition::Bjorhus};
  SectorImposition v_minus{SectorImposition::Bjorhus};
};

/*!
 * \brief The inertial-frame center of the excision sphere that the element
 * abuts.
 *
 * \details Looks the element up in the abutting directions of every excision
 * sphere and maps that sphere's grid-frame center to the inertial frame with
 * its time-dependent map, if it has one. Errors if the element abuts no
 * excision sphere, since the worldtube boundary condition is only meaningful
 * on an excision boundary.
 */
template <size_t Dim>
tnsr::I<double, Dim, Frame::Inertial> excision_sphere_center(
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const ElementId<Dim>& element_id, double time,
    const domain::FunctionsOfTimeMap& functions_of_time);

/*!
 * \brief The electric and magnetic parts of the Weyl tensor on the face from
 * the generalized harmonic variables.
 *
 * \details The spatial Ricci tensor comes from \f$\Phi\f$ and its derivative
 * and the covariant derivative of the extrinsic curvature from
 * `gh::covariant_deriv_of_extrinsic_curvature()`; then
 * `gr::weyl_electric()` and `gr::weyl_magnetic()`. Unlike the physical
 * Bjorhus term, no four-index-constraint terms are added to the Ricci tensor:
 * the model consumes the curvature as it is.
 */
void face_weyl_electric_magnetic(
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

/*!
 * \brief The incoming Weyl mode \f$U^{8-}_{ij}\f$ of the type-D model, in
 * the normalization of `gr::weyl_propagating()` with sign \f$-1\f$.
 *
 * \details The worldtube-adapted tetrad is built with \f$s\f$ along
 * `unit_normal_covector`, the outward normal of the domain, which at the
 * excision points into the hole. In that tetrad the mode entering the domain
 * is \f$\Psi_0\f$ and
 * \f$U^{8-}_{ij} = \tfrac12 w^-_{ij}
 *  = \bar\Psi_0 m_i m_j + \Psi_0 \bar m_i \bar m_j\f$
 * (checked against `gr::weyl_propagating()` in the unit tests). \f$\Psi_0\f$
 * is the type-D value \f$6 b^2 \Psi_2^K\f$ measured from the face curvature
 * by `gr::np::solve_type_d_rotation()`. The result is in covariant coordinate
 * components.
 */
tnsr::ii<DataVector, 3, Frame::Inertial> type_d_incoming_mode(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial>& unit_normal_covector);
}  // namespace gh::BoundaryConditions::detail

namespace gh::BoundaryConditions {
/*!
 * \brief Bjorhus-type boundary condition for the inner boundary of a worldtube
 * around a black hole, with each characteristic sector imposed independently.
 *
 * \details The worldtube excision scheme excises a region around the smaller
 * black hole of a binary that is larger than its horizon, so the excision
 * boundary needs boundary conditions. This class imposes them on the
 * generalized harmonic system with Bjorhus' method \cite Bjorhus1995 as
 * corrections to the characteristic projections of the time derivatives of
 * the evolved variables, Eqs. (63) - (65) of \cite Lindblom2005qh, in the same
 * way as the outer-boundary condition `ConstraintPreservingBjorhus`, with
 * which it shares its implementation (see `Bjorhus::IntermediateVariables`).
 *
 * The incoming field \f$u^-\f$ splits exactly into three sectors, see
 * `Bjorhus::detail::add_constraint_sector_projection()` and its siblings:
 *
 * - the constraint-preserving sector, which together with the
 *   constraint-preserving corrections to \f$\partial_t v_\psi\f$ and
 *   \f$\partial_t v_0\f$ prevents the influx of constraint violations;
 * - the physical sector, the transverse-traceless part that carries the
 *   incoming gravitational radiation, i.e. the Newman-Penrose scalar
 *   \f$\Psi_0\f$;
 * - the gauge sector.
 *
 * Each sector is imposed independently, as selected at runtime by the options
 * `ConstraintPreservingSector`, `PhysicalSector` and `GaugeSector`, either
 * through its Bjorhus correction (`Bjorhus`), or not at all (`Frozen`, the
 * sector's time derivative is set to zero; a diagnostic that lets whatever
 * is in the volume flow in unchecked). For the gauge sector the model-free
 * Sommerfeld conditions `SommerfeldAbsorbing` and `SommerfeldOutgoing`
 * are available instead of `Bjorhus`; their \f$1/r\f$ term uses the distance
 * from the center of the excision sphere, which the boundary condition looks
 * up in the domain, so they are correct for a hole that is off-center or
 * moving. The constraint-preserving sector can also be imposed per field, see
 * `detail::PerFieldConstraintSectors`. The field \f$v^+\f$ is always frozen,
 * as in `ConstraintPreservingBjorhus`.
 *
 * The physical sector's Bjorhus term drives the incoming Weyl mode
 * \f$U^{8-}\f$ towards a value supplied by the `PhysicalModel`. With `None`
 * that value is zero, i.e. no incoming radiation, exactly like the
 * outer-boundary condition. With `TypeD` it is the order-zero curvature
 * matching of the worldtube NP-matching note, see
 * `detail::type_d_incoming_mode()`: the face curvature is read as a type-D
 * hole seen from a misaligned tetrad and the radiation such a hole sends into
 * the domain is supplied. The construction consumes the interior \f$\Psi_0\f$
 * through the invariants, so the evolution itself is the fixed-point
 * iteration of the note's circularity discussion.
 *
 * \note Unlike `ConstraintPreservingBjorhus`, this condition does not reject a
 * mesh velocity along the outward normal. At an inner boundary the outward
 * normal points into the excision, so an excision tracking a moving hole has
 * such a mesh velocity over half of the sphere by construction. The mesh
 * velocity is accounted for in the inertial time derivatives and in the
 * characteristic speeds.
 */
template <size_t Dim>
class WorldtubeTypeD final : public BoundaryCondition<Dim> {
 public:
  /// \brief How the constraint-preserving sector is imposed.
  ///
  /// Either one imposition for the whole sector, or a
  /// `detail::PerFieldConstraintSectors` map imposing its three terms
  /// separately.
  struct ConstraintPreservingSector {
    using type = std::variant<detail::SectorImposition,
                              detail::PerFieldConstraintSectors>;
    static constexpr Options::String help{
        "Imposition of the constraint-preserving sector, which comprises "
        "v_psi, v_zero and the constraint projection of v_minus. Either "
        "Bjorhus/Frozen for all three at once, or a map with the keys VPsi, "
        "VZero and VMinus to impose them separately. Frozen sets time "
        "derivatives to zero and is a DIAGNOSTIC ONLY: it deliberately lets "
        "constraint violations enter, which is what makes it useful for "
        "isolating where an influx comes from."};
  };
  /// \brief How the physical sector -- the transverse-traceless projection of
  /// \f$u^-\f$ carrying the incoming gravitational radiation -- is imposed.
  struct PhysicalSector {
    using type = detail::SectorImposition;
    static constexpr Options::String help{
        "Bjorhus or Frozen imposition of the physical sector of v_minus. "
        "Frozen sets dt of that projection to zero; like the constraint "
        "sector's Frozen it is a diagnostic, not a production setting."};
  };
  /// \brief How the gauge sector of \f$u^-\f$ is imposed.
  struct GaugeSector {
    using type = detail::SectorImposition;
    static constexpr Options::String help{
        "Imposition of the gauge sector of v_minus: Frozen (no condition at "
        "all), or SommerfeldAbsorbing / SommerfeldOutgoing for a model-free "
        "Bayliss-Turkel L=0 condition with the radius measured from the "
        "excision center. SommerfeldAbsorbing annihilates a wave falling into "
        "the hole and is the physically motivated sign at an inner boundary; "
        "SommerfeldOutgoing is the outer-boundary sign, provided for "
        "comparison. Bjorhus is not implemented."};
  };
  /// \brief Which model supplies the incoming radiation to the physical
  /// sector's Bjorhus term.
  struct PhysicalModel {
    using type = detail::PhysicalModel;
    static constexpr Options::String help{
        "The incoming radiation the physical Bjorhus term drives the boundary "
        "towards: None (no incoming radiation, the outer-boundary condition) "
        "or TypeD (order zero of the curvature matching: the excised region "
        "is a type-D hole, Psi0 = 6 b^2 Psi2^K from the face curvature). "
        "Requires PhysicalSector: Bjorhus."};
  };
  using options = tmpl::list<ConstraintPreservingSector, PhysicalSector,
                             GaugeSector, PhysicalModel>;
  static constexpr Options::String help{
      "Bjorhus-type boundary condition for the inner boundary of a worldtube "
      "around a black hole. Each of the three characteristic sectors of "
      "v_minus is imposed independently, through its Bjorhus time-derivative "
      "correction or not at all (Frozen); the gauge sector can use a "
      "Sommerfeld condition instead. The sectors partition v_minus exactly, "
      "so the choices are independent. The physical sector's Bjorhus term "
      "drives the incoming radiation towards the value of PhysicalModel."};
  static std::string name() { return "WorldtubeTypeD"; }

  WorldtubeTypeD(
      std::variant<detail::SectorImposition, detail::PerFieldConstraintSectors>
          constraint_preserving_sector,
      detail::SectorImposition physical_sector,
      detail::SectorImposition gauge_sector,
      detail::PhysicalModel physical_model,
      const Options::Context& context = {});

  WorldtubeTypeD() = default;
  /// \cond
  WorldtubeTypeD(WorldtubeTypeD&&) = default;
  WorldtubeTypeD& operator=(WorldtubeTypeD&&) = default;
  WorldtubeTypeD(const WorldtubeTypeD&) = default;
  WorldtubeTypeD& operator=(const WorldtubeTypeD&) = default;
  /// \endcond
  ~WorldtubeTypeD() override = default;

  explicit WorldtubeTypeD(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, WorldtubeTypeD);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::TimeDerivative;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, Dim>,
                 Tags::Pi<DataVector, Dim>, Tags::Phi<DataVector, Dim>>;
  using dg_interior_temporary_tags =
      tmpl::list<domain::Tags::Coordinates<Dim, Frame::Inertial>,
                 Tags::ConstraintGamma1, Tags::ConstraintGamma2,
                 gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
                 gr::Tags::InverseSpacetimeMetric<DataVector, Dim>,
                 gr::Tags::SpacetimeNormalVector<DataVector, Dim>,
                 Tags::ThreeIndexConstraint<DataVector, Dim>,
                 Tags::GaugeH<DataVector, Dim>,
                 Tags::SpacetimeDerivGaugeH<DataVector, Dim>>;
  using dg_interior_dt_vars_tags =
      tmpl::list<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>,
                 ::Tags::dt<Tags::Pi<DataVector, Dim>>,
                 ::Tags::dt<Tags::Phi<DataVector, Dim>>>;
  using dg_interior_deriv_vars_tags =
      tmpl::list<::Tags::deriv<gr::Tags::SpacetimeMetric<DataVector, Dim>,
                               tmpl::size_t<Dim>, Frame::Inertial>,
                 ::Tags::deriv<Tags::Pi<DataVector, Dim>, tmpl::size_t<Dim>,
                               Frame::Inertial>,
                 ::Tags::deriv<Tags::Phi<DataVector, Dim>, tmpl::size_t<Dim>,
                               Frame::Inertial>>;
  using dg_gridless_tags =
      tmpl::list<::Tags::Time, domain::Tags::Domain<Dim>,
                 domain::Tags::Element<Dim>, domain::Tags::FunctionsOfTime>;

  std::optional<std::string> dg_time_derivative(
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          dt_spacetime_metric_correction,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          dt_pi_correction,
      gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
          dt_phi_correction,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& /*normal_vector*/,
      // c.f. dg_interior_evolved_variables_tags
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
      // c.f. dg_interior_temporary_tags
      const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
      const tnsr::AA<DataVector, Dim, Frame::Inertial>&
          inverse_spacetime_metric,
      const tnsr::A<DataVector, Dim, Frame::Inertial>&
          spacetime_unit_normal_vector,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& three_index_constraint,
      const tnsr::a<DataVector, Dim, Frame::Inertial>& gauge_source,
      const tnsr::ab<DataVector, Dim, Frame::Inertial>&
          spacetime_deriv_gauge_source,
      // c.f. dg_interior_dt_vars_tags
      const tnsr::aa<DataVector, Dim, Frame::Inertial>&
          logical_dt_spacetime_metric,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& logical_dt_pi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& logical_dt_phi,
      // c.f. dg_interior_deriv_vars_tags
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_spacetime_metric,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
      const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi,
      // c.f. dg_gridless_tags
      double time, const Domain<Dim>& domain, const Element<Dim>& element,
      const domain::FunctionsOfTimeMap& functions_of_time) const;

  detail::SectorImposition constraint_v_psi() const {
    return constraint_v_psi_;
  }
  detail::SectorImposition constraint_v_zero() const {
    return constraint_v_zero_;
  }
  detail::SectorImposition constraint_v_minus() const {
    return constraint_v_minus_;
  }
  detail::SectorImposition physical_sector() const { return physical_sector_; }
  detail::SectorImposition gauge_sector() const { return gauge_sector_; }
  detail::PhysicalModel physical_model() const { return physical_model_; }

 private:
  detail::SectorImposition constraint_v_psi_{detail::SectorImposition::Bjorhus};
  detail::SectorImposition constraint_v_zero_{
      detail::SectorImposition::Bjorhus};
  detail::SectorImposition constraint_v_minus_{
      detail::SectorImposition::Bjorhus};
  detail::SectorImposition physical_sector_{detail::SectorImposition::Bjorhus};
  detail::SectorImposition gauge_sector_{
      detail::SectorImposition::SommerfeldAbsorbing};
  detail::PhysicalModel physical_model_{detail::PhysicalModel::None};
};

template <size_t Dim>
bool operator==(const WorldtubeTypeD<Dim>& lhs, const WorldtubeTypeD<Dim>& rhs);
template <size_t Dim>
bool operator!=(const WorldtubeTypeD<Dim>& lhs, const WorldtubeTypeD<Dim>& rhs);
}  // namespace gh::BoundaryConditions

template <>
struct Options::create_from_yaml<
    gh::BoundaryConditions::detail::SectorImposition> {
  template <typename Metavariables>
  static typename gh::BoundaryConditions::detail::SectorImposition create(
      const Options::Option& options) {
    return gh::BoundaryConditions::detail::convert_sector_imposition_from_yaml(
        options);
  }
};

template <>
struct Options::create_from_yaml<
    gh::BoundaryConditions::detail::PhysicalModel> {
  template <typename Metavariables>
  static typename gh::BoundaryConditions::detail::PhysicalModel create(
      const Options::Option& options) {
    return gh::BoundaryConditions::detail::convert_physical_model_from_yaml(
        options);
  }
};
