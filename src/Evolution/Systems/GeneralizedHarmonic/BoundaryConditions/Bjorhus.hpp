// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
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
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
namespace Tags {
struct Time;
}  // namespace Tags
/// \endcond

namespace gh::BoundaryConditions::detail {
enum class ConstraintPreservingBjorhusType {
  ConstraintPreserving,
  ConstraintPreservingPhysical
};

ConstraintPreservingBjorhusType
convert_constraint_preserving_bjorhus_type_from_yaml(
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
 * An optional injected incoming-wave contribution can be specified through
 * `IncomingWaveProfile`. The injected strain-rate tensor is
 * \f[
 *   \dot{h}_{ij} = f'(t)\, h_{ij},
 * \f]
 * where the configured profile is the envelope \f$f(t)\f$ and the injection
 * uses its *derivative*. With a Gaussian envelope
 * \f$f(t) = A \exp[-(t - t_p)^2 / w^2]\f$ this makes the pulse bipolar, so the
 * strain returns to zero once the pulse has passed; injecting \f$f\f$ itself
 * would leave a permanent offset, a step rather than a pulse. Choose
 * \f$t_p \gtrsim 4 w\f$ so the profile is negligible at the initial time.
 * \f$h_{ij}\f$ is the constant symmetric spatial tensor given componentwise by
 * `IncomingWaveComponents` in the order \f$(xx, xy, xz, yy, yz, zz)\f$. Only
 * the spatial block is set; the correction enters \f$\partial_t u^-_{ab}\f$
 * through the transverse-traceless projection
 * \f$P^c{}_a P^d{}_b - \frac{1}{2} P_{ab} P^{cd}\f$, so \f$h_{ij}\f$ need be
 * neither transverse nor trace free -- whatever is not transverse and trace
 * free with respect to the boundary normal is projected away.
 *
 * A *constant* Cartesian tensor can only excite \f$\ell = 2\f$; reaching
 * higher \f$\ell\f$ would need an angle-dependent amplitude on the sphere. The
 * default \f$\mathrm{diag}(1, 1, -2)\f$ is proportional to the \f$(2, 0)\f$
 * tensor of \cite Lindblom2005qh and so reaches essentially only \f$m = 0\f$,
 * while populating all five independent components spans the whole
 * \f$\ell = 2\f$ multiplet -- the better choice when the point is to probe a
 * boundary condition rather than to reproduce that reference.
 *
 * When comparing pulse *shapes*, note that the peak of \f$\dot{f}\f$ scales as
 * \f$A/w\f$ for a Gaussian of amplitude \f$A\f$ and width \f$w\f$: hold
 * \f$A/w\f$ and the Frobenius norm of \f$h_{ij}\f$ fixed, or the comparison
 * confounds strength with shape.
 *
 * It should be considered an approximate perturbation rather than an exact
 * gravitational wave, which would have to be constructed at null infinity.
 * Note that the profile of the injected wave is specified at the outer
 * boundary and its amplitude should be adjusted according to the usual 1/r
 * scaling.
 *
 * This class provides two choices of combinations of the above corrections:
 *  - `ConstraintPreserving` : this imposes the constraint-preserving and
 * gauge-controlling corrections;
 *  - `ConstraintPreservingPhysical` : this additionally restricts the influx of
 * any physical gravitational waves from the outer boundary, in addition to
 * preventing the influx of constraint violations and gauge perturbations.
 *
 * We refer to `Bjorhus::constraint_preserving_corrections_dt_v_psi()`,
 * `Bjorhus::constraint_preserving_corrections_dt_v_zero()`,
 * `Bjorhus::constraint_preserving_gauge_corrections_dt_v_minus()`, and
 * `Bjorhus::constraint_preserving_gauge_physical_corrections_dt_v_minus()`
 * for the further details on implementation.
 *
 * \note These boundary conditions assume a spherical outer boundary.
 */
template <size_t Dim>
class ConstraintPreservingBjorhus final : public BoundaryCondition<Dim> {
 public:
  struct TypeOptionTag {
    using type = detail::ConstraintPreservingBjorhusType;
    static std::string name() { return "Type"; }
    static constexpr Options::String help{
        "Whether to impose ConstraintPreserving, with or without physical "
        "terms for VMinus."};
  };

  struct IncomingWaveProfileOptionTag {
    using type =
        Options::Auto<std::unique_ptr<::MathFunction<1, Frame::Inertial>>,
                      Options::AutoLabel::None>;
    static std::string name() { return "IncomingWaveProfile"; }
    static constexpr Options::String help{
        "Optional envelope f(t) for the injected physical wave. The injected "
        "strain rate is its DERIVATIVE, f'(t), so a Gaussian envelope gives a "
        "bipolar pulse whose strain returns to zero; injecting f itself would "
        "leave a permanent offset. Choose the peak time >= 4 widths so the "
        "profile is negligible at the initial time. See the "
        "ConstraintPreservingBjorhus class documentation for the injected-wave "
        "formula. Specify `None` to disable injection. This option is only "
        "supported in 3D."};
  };

  struct IncomingWaveComponentsOptionTag {
    using type = std::array<double, 6>;
    static std::string name() { return "IncomingWaveComponents"; }
    static constexpr Options::String help{
        "Components (xx, xy, xz, yy, yz, zz) of the constant symmetric spatial "
        "tensor the injected wave is proportional to. They need be neither "
        "transverse nor trace free: the transverse-traceless projection keeps "
        "only the part that is, relative to the boundary normal. A constant "
        "tensor cannot reach beyond l = 2. The default diag(1, 1, -2) is "
        "proportional to the (2, 0) tensor and so reaches essentially only "
        "m = 0; populating all five independent components spans the whole "
        "l = 2 multiplet. Ignored when IncomingWaveProfile is None."};
    static type default_value() {
      return Bjorhus::default_incoming_wave_components;
    }
  };

  using options = tmpl::flatten<tmpl::list<
      TypeOptionTag,
      tmpl::conditional_t<Dim == 3,
                          tmpl::list<IncomingWaveProfileOptionTag,
                                     IncomingWaveComponentsOptionTag>,
                          tmpl::list<>>>>;
  static constexpr Options::String help{
      "ConstraintPreservingBjorhus boundary conditions setting the value of the"
      "time derivatives of the spacetime metric, Phi and Pi to expressions that"
      "prevent the influx of constraint violations and reflections."};
  static std::string name() { return "ConstraintPreservingBjorhus"; }

  explicit ConstraintPreservingBjorhus(
      detail::ConstraintPreservingBjorhusType type,
      std::optional<std::unique_ptr<::MathFunction<1, Frame::Inertial>>>
          incoming_wave_profile = std::nullopt,
      const std::array<double, 6>& incoming_wave_components =
          Bjorhus::default_incoming_wave_components);

  ConstraintPreservingBjorhus() = default;
  /// \cond
  ConstraintPreservingBjorhus(ConstraintPreservingBjorhus&&) = default;
  ConstraintPreservingBjorhus& operator=(ConstraintPreservingBjorhus&&) =
      default;
  ConstraintPreservingBjorhus(const ConstraintPreservingBjorhus&);
  ConstraintPreservingBjorhus& operator=(const ConstraintPreservingBjorhus&);
  /// \endcond
  ~ConstraintPreservingBjorhus() override = default;

  explicit ConstraintPreservingBjorhus(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition,
      ConstraintPreservingBjorhus);

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
  using dg_gridless_tags = tmpl::list<::Tags::Time>;

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
      double time = std::numeric_limits<double>::signaling_NaN()) const;

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

  detail::ConstraintPreservingBjorhusType type_{
      detail::ConstraintPreservingBjorhusType::ConstraintPreservingPhysical};
  std::unique_ptr<::MathFunction<1, Frame::Inertial>> incoming_wave_profile_{};
  std::array<double, 6> incoming_wave_components_{
      Bjorhus::default_incoming_wave_components};
};
}  // namespace gh::BoundaryConditions

template <>
struct Options::create_from_yaml<
    gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType> {
  template <typename Metavariables>
  static
      typename gh::BoundaryConditions::detail::ConstraintPreservingBjorhusType
      create(const Options::Option& options) {
    return gh::BoundaryConditions::detail::
        convert_constraint_preserving_bjorhus_type_from_yaml(options);
  }
};
