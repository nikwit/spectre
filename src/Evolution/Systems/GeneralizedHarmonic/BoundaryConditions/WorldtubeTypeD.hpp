// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace gh::BoundaryConditions::detail {
enum class WorldtubeTypeDType {
  ConstraintPreserving,
  /// Constraint-preserving, physical (inferred \f$\Psi_0\f$), and the
  /// Bayliss-Turkel Sommerfeld gauge condition.
  ConstraintPreservingPhysical,
  /// As `ConstraintPreservingPhysical`, but with the gauge sector frozen,
  /// \f$\partial_t u^-_{ab}|_{\rm gauge} = 0\f$, instead of the Sommerfeld
  /// condition. The Sommerfeld condition is an outer-boundary prescription
  /// whose \f$1/r\f$ is measured from the origin of the inertial coordinates,
  /// so it is correct at a worldtube only when the worldtube is centred on that
  /// origin; freezing needs no model input and isolates the effect of the gauge
  /// condition.
  ConstraintPreservingPhysicalFrozenGauge,
  /// As `ConstraintPreservingPhysical`, but the gauge sector relaxes towards
  /// the value of \f$u^-_{ab}\f$ computed from an analytic prescription:
  /// \f$\partial_t u^-_{ab}|_{\rm gauge}
  ///    = -\kappa\,(u^-_{ab} - u^{-,\rm model}_{ab})|_{\rm gauge}\f$.
  ///
  /// This is the "ghost" form of the gauge condition: it needs only
  /// \f$u^{-,\rm model}\f$ and not its time derivative, which would require a
  /// second time derivative of the model. Requires
  /// `AnalyticGaugePrescription`.
  ConstraintPreservingPhysicalAnalyticGhostGauge,
  /// Both gauge terms together:
  /// \f$\partial_t u^-_{ab}|_{\rm gauge}
  ///    = -(\gamma_2 - 1/r)(\partial_t u^g_{ab})|_{\rm gauge}
  ///      - \kappa\,(u^-_{ab} - u^{-,\rm model}_{ab})|_{\rm gauge}\f$.
  ///
  /// The two terms do different jobs and are not alternatives. The Sommerfeld
  /// term is a *radiation* condition: it lets whatever gauge dynamics the
  /// interior generates leave without strong reflection, and it is the term
  /// that matters when \f$\partial_t g\f$ at the boundary is large. The
  /// relaxation is *Dirichlet-like*: it supplies the gauge information that
  /// flows in from the excised region, and carries content when the
  /// configuration is nearly stationary and \f$\partial_t g \to 0\f$ makes the
  /// Sommerfeld term vanish.
  ///
  /// Requires `AnalyticGaugePrescription`.
  ConstraintPreservingPhysicalSommerfeldGhostGauge,
  /// As `ConstraintPreservingPhysicalAnalyticGhostGauge`, but the model
  /// \f$u^{-,\rm model}\f$ is built from the *online* worldtube matcher's
  /// latest fitted map parameters (`gh::Worldtube::Tags::MapParameters`,
  /// extrapolated linearly within the fit interval with the fitted rates)
  /// instead of an input-file prescription — the closed loop. Requires the
  /// `WorldtubeMatcher` option to be active; while no fit exists yet the
  /// gauge sector stays frozen. 3D only.
  ConstraintPreservingPhysicalOnlineGhostGauge,
  /// Constraint-preserving and physical sectors by Bjorhus time-derivative
  /// corrections as in `ConstraintPreservingPhysical*`; the **gauge sector
  /// is imposed weakly through ghost data**: the exterior state equals the
  /// interior except that the gauge projection of \f$u^-\f$ is replaced by
  /// the online matcher's model, and the upwind penalty drives the jump. No
  /// relaxation rate: the driving strength is set by the characteristic
  /// speed and the DG lifting. The gauge sector receives no Bjorhus
  /// correction (its volume dynamics stays free). Requires the
  /// `WorldtubeMatcher` option; while no fit exists the ghost state equals
  /// the interior (no driving). 3D only.
  ConstraintPreservingPhysicalGhostGauge
};

WorldtubeTypeDType convert_worldtube_type_d_type_from_yaml(
    const Options::Option& options);
}  // namespace gh::BoundaryConditions::detail

namespace gh::BoundaryConditions {
/*!
 * \brief Sets constraint preserving boundary conditions using the Bjorhus
 * method.
 *
 * \details Boundary conditions for the generalized harmonic evolution system
 * can be divided in to three parts, constraint-preserving, physical and gauge
 * boundary conditions.
 *
 * The generalized harmonic (GH) evolution system is a first-order reduction of
 * Einstein equations brought about by the imposition of GH gauge. This
 * introduces constraints on the free (evolved) variables in addition to the
 * standard Hamiltonian and momentum constraints. The constraint-preserving
 * portion of the boundary conditions is designed to prevent the influx of
 * constraint violations from external faces of the evolution domain, by damping
 * them away on a controlled and short time-scale. These conditions are imposed
 * as corrections to the characteristic projections of the right-hand-sides of
 * the GH evolution equations (i.e. using Bjorhus' method \cite Bjorhus1995),
 * as written down in Eq. (63) - (65) of \cite Lindblom2005qh . In addition to
 * these equations, the fourth projection is simply frozen in the unlikely case
 * its coordinate speed becomes negative, i.e. \f$d_t u^{\hat{1}+}{}_{ab}=0\f$
 * (in the notation of \cite Lindblom2005qh). The gauge degrees
 * of freedom are controlled by imposing a Sommerfeld-type condition (\f$L=0\f$
 * member of the hierarchy derived in \cite BaylissTurkel) that allow gauge
 * perturbations to pass through the boundary without strong reflections. These
 * assume a spherical outer boundary, and can be written down as in Eq. (25) of
 * \cite Rinne2007ui . Finally, the physical boundary conditions control the
 * influx of inward propagating gravitational-wave solutions from the external
 * boundaries. These are derived by considering the evolution system of the Weyl
 * curvature tensor, and controlling the inward propagating characteristics of
 * the system that are proportional to the Newman-Penrose curvature spinor
 * components \f$\Psi_4\f$ and \f$\Psi_0\f$. Here we use Eq. (68) of
 * \cite Lindblom2005qh to disallow any incoming waves. It is to be noted that
 * all the above conditions are also imposed on characteristic modes with speeds
 * exactly zero.
 *
 * See `detail::WorldtubeTypeDType` for the available combinations. They differ
 * only in the gauge sector: Sommerfeld, frozen, relaxation towards a model
 * \f$u^-\f$, or Sommerfeld and the relaxation together.
 *
 * We refer to `Bjorhus::constraint_preserving_corrections_dt_v_psi()`,
 * `Bjorhus::constraint_preserving_corrections_dt_v_zero()`, and
 * and the `_worldtube` variant of
 * `Bjorhus::constraint_preserving_gauge_physical_corrections_dt_v_minus`,
 * for the further details on implementation.
 *
 * \note Unlike `ConstraintPreservingBjorhus`, from which much of the above was
 * inherited, this is an **inner** boundary condition applied at a worldtube
 * around a hole, and several statements above have to be read with that in
 * mind. In particular the physical sector does not "disallow any incoming
 * waves": the incoming radiative field is *supplied*, from the \f$\Psi_0\f$
 * inferred by the type-D inversion, because at an excision that field
 * carries physical information out of the region that was removed. Likewise the
 * Sommerfeld condition's \f$1/r\f$ is the distance from the origin of the
 * inertial coordinates, which is the local worldtube radius only when the
 * worldtube is centred there: true for a single hole at the origin,
 * false for a drifting or binary hole.
 */
template <size_t Dim>
class WorldtubeTypeD final : public BoundaryCondition<Dim> {
 public:
  struct TypeOptionTag {
    using type = detail::WorldtubeTypeDType;
    static std::string name() { return "Type"; }
    static constexpr Options::String help{
        "Whether to impose ConstraintPreserving, with or without physical "
        "terms for VMinus."};
  };

  /// \brief Analytic prescription supplying the model value of
  /// \f$u^-_{ab}|_{\rm gauge}\f$.
  ///
  /// Required for `Type: ConstraintPreservingPhysicalAnalyticGhostGauge`, and
  /// `None` otherwise. In the worldtube matching scheme this is where the
  /// fitted tidally-perturbed metric would enter instead; an analytic solution
  /// is used here because it makes the model exactly known, so any deviation is
  /// attributable to the boundary condition.
  struct AnalyticGaugePrescription {
    using type = Options::Auto<
        std::unique_ptr<evolution::initial_data::InitialData>,
        Options::AutoLabel::None>;
    static std::string name() { return "AnalyticGaugePrescription"; }
    static constexpr Options::String help{
        "Analytic solution supplying the model value of u^-|gauge. Required "
        "for Type: ConstraintPreservingPhysicalAnalyticGhostGauge, else None."};
  };

  /// \brief Relaxation rate \f$\kappa\f$, in inverse mass units, towards the
  /// model value of the gauge sector.
  struct GaugeRelaxationRate {
    using type = double;
    static std::string name() { return "GaugeRelaxationRate"; }
    static constexpr Options::String help{
        "Relaxation rate kappa towards the model gauge sector. Only used for "
        "Type: ConstraintPreservingPhysicalAnalyticGhostGauge. kappa = 0 "
        "reduces to freezing the gauge sector."};
  };

  using options = tmpl::list<TypeOptionTag, AnalyticGaugePrescription,
                             GaugeRelaxationRate>;
  static constexpr Options::String help{
      "WorldtubeTypeD boundary conditions setting the value of the time "
      "derivatives of the spacetime metric, Phi and Pi to expressions that "
      "prevent the influx of constraint violations and reflections."};
  static std::string name() { return "WorldtubeTypeD"; }

  WorldtubeTypeD(
      detail::WorldtubeTypeDType type,
      std::optional<std::unique_ptr<evolution::initial_data::InitialData>>
          analytic_gauge_prescription,
      double gauge_relaxation_rate, const Options::Context& context = {});

  WorldtubeTypeD() = default;
  /// \cond
  WorldtubeTypeD(WorldtubeTypeD&&) = default;
  WorldtubeTypeD& operator=(WorldtubeTypeD&&) = default;
  WorldtubeTypeD(const WorldtubeTypeD&);
  WorldtubeTypeD& operator=(const WorldtubeTypeD&);
  /// \endcond
  ~WorldtubeTypeD() override = default;

  explicit WorldtubeTypeD(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition,
      WorldtubeTypeD);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  // GhostAndTimeDerivative: the constraint-preserving and physical sectors
  // are imposed by Bjorhus time-derivative corrections; for
  // `ConstraintPreservingPhysicalGhostGauge` the gauge sector is imposed
  // weakly through ghost data (the exterior state's u^- has its gauge
  // projection replaced by the online-matcher model) driven by the upwind
  // penalty. For all other Types the ghost state is a copy of the interior,
  // so the penalty contributes exactly zero and the behavior is unchanged.
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::GhostAndTimeDerivative;

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
      tmpl::list<::Tags::Time, gh::Worldtube::Tags::Matcher,
                 gh::Worldtube::Tags::MapParameters>;

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
      double time,
      const std::optional<gh::Worldtube::MatcherConfig>& matcher_config,
      const gh::Worldtube::MapParameterData& map_parameters) const;

  std::optional<std::string> dg_ghost(
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          spacetime_metric_ghost,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> pi_ghost,
      gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> phi_ghost,
      gsl::not_null<Scalar<DataVector>*> gamma1_ghost,
      gsl::not_null<Scalar<DataVector>*> gamma2_ghost,
      gsl::not_null<Scalar<DataVector>*> lapse_ghost,
      gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> shift_ghost,
      gsl::not_null<tnsr::II<DataVector, Dim, Frame::Inertial>*>
          inv_spatial_metric_ghost,

      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& normal_vector,
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
      double time,
      const std::optional<gh::Worldtube::MatcherConfig>& matcher_config,
      const gh::Worldtube::MapParameterData& map_parameters) const;

 private:
  void compute_intermediate_vars(
      gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          unit_interface_normal_vector,
      gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
          four_index_constraint,
      gsl::not_null<tnsr::II<DataVector, Dim, Frame::Inertial>*>
          inverse_spatial_metric,
      gsl::not_null<tnsr::ii<DataVector, Dim, Frame::Inertial>*>
          extrinsic_curvature,
      gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
          incoming_null_one_form,
      gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
          outgoing_null_one_form,
      gsl::not_null<tnsr::A<DataVector, Dim, Frame::Inertial>*>
          incoming_null_vector,
      gsl::not_null<tnsr::A<DataVector, Dim, Frame::Inertial>*>
          outgoing_null_vector,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> projection_ab,
      gsl::not_null<tnsr::Ab<DataVector, Dim, Frame::Inertial>*> projection_Ab,
      gsl::not_null<tnsr::AA<DataVector, Dim, Frame::Inertial>*> projection_AB,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          char_projected_rhs_dt_v_psi,
      gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
          char_projected_rhs_dt_v_zero,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          char_projected_rhs_dt_v_plus,
      gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
          char_projected_rhs_dt_v_minus,
      gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
          constraint_char_zero_plus,
      gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
          constraint_char_zero_minus,
      gsl::not_null<std::array<DataVector, 4>*> char_speeds,

      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
      const tnsr::AA<DataVector, Dim, Frame::Inertial>&
          inverse_spacetime_metric,
      const tnsr::A<DataVector, Dim, Frame::Inertial>&
          spacetime_unit_normal_vector,
      const tnsr::a<DataVector, Dim, Frame::Inertial>&
          spacetime_unit_normal_one_form,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& three_index_constraint,
      const tnsr::a<DataVector, Dim, Frame::Inertial>& gauge_source,
      const tnsr::ab<DataVector, Dim, Frame::Inertial>&
          spacetime_deriv_gauge_source,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& dt_pi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& dt_phi,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& dt_spacetime_metric,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
      const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_spacetime_metric)
      const;

  detail::WorldtubeTypeDType type_{
      detail::WorldtubeTypeDType::ConstraintPreservingPhysical};
  std::unique_ptr<evolution::initial_data::InitialData>
      analytic_gauge_prescription_{nullptr};
  double gauge_relaxation_rate_{0.};
};
}  // namespace gh::BoundaryConditions

template <>
struct Options::create_from_yaml<
    gh::BoundaryConditions::detail::WorldtubeTypeDType> {
  template <typename Metavariables>
  static typename gh::BoundaryConditions::detail::WorldtubeTypeDType create(
      const Options::Option& options) {
    return gh::BoundaryConditions::detail::
        convert_worldtube_type_d_type_from_yaml(options);
  }
};
