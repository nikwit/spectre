// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Interpolation/CubicSpline.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gh::Solutions {
/// Free-function core of the first-order affine-map model, shared by the
/// tabulated solution class below and the online worldtube matcher.
namespace affine_map_model {
/// out^{ab} = background_weight * G_0^{ab}(y) + sum_A c_A R_A^{ab}(y), with
/// G_0 the harmonic-Schwarzschild inverse metric, R_A the thirteen
/// first-order response columns of the matching scheme, and y coordinates
/// relative to the centre. Transcribes
/// worldtube_matching/responses.py.
void inverse_metric_combination(gsl::not_null<tnsr::AA<DataVector, 3>*> out,
                                const std::array<DataVector, 3>& y, double mass,
                                double background_weight,
                                const std::array<double, 13>& c);

/// The analytic spatial derivative \f$\partial_k\f$ of
/// `inverse_metric_combination` with respect to the local coordinates y.
void spatial_derivative_of_inverse_metric_combination(
    gsl::not_null<tnsr::iAA<DataVector, 3>*> out,
    const std::array<DataVector, 3>& y, double mass, double background_weight,
    const std::array<double, 13>& c);

/*!
 * \brief Strict first-order slow-time GH variables of the affine-map model.
 *
 * With \f$\widetilde t=\epsilon t\f$ and
 * \f$p_A=p_A(\widetilde t)\f$, this function retains exactly the
 * \f$O(\epsilon)\f$ coefficient of the covariant metric and its spatial
 * derivative.  It does not resum the inverse metric and it does not include
 * \f$D_t p_A=O(\epsilon)\f$, whose contribution to the metric derivative is
 * \f$O(\epsilon^2)\f$.  The one retained time derivative is the kinematic
 * motion of the zeroth-order center,
 * \f$\partial_t g^{(0)}_{ab}=-\dot q^k\Phi^{(0)}_{kab}\f$; multiplying the
 * perturbed \f$\Phi\f$ would again introduce \f$O(\epsilon^2)\f$ terms.
 *
 * The returned tensors have the bookkeeping factor \f$\epsilon\f$ set to
 * one after truncation:
 * \f$g=g^{(0)}+g^{(1)}\f$,
 * \f$\Phi=\Phi^{(0)}+\Phi^{(1)}\f$, and
 * \f$\Pi=\Pi^{(0)}+\Pi^{(1)}\f$.
 */
void first_order_evolved_variables(
    gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    bool centre_advection = true);

/*!
 * \brief Strict first-order affine-map variables on an exactly Lorentz-
 * boosted Schwarzschild background.
 *
 * `boost_velocity` is not order expanded.  The supplied `center` is the
 * instantaneous lab-frame center of the hole on the slice being evaluated.
 * The Schwarzschild background and its derivatives are transformed exactly,
 * while the affine-map perturbation and the lab-frame 3+1 reconstruction are
 * expanded once in `p`.  Consequently the result contains all powers of the
 * bulk velocity but no products of affine-map coefficients.
 *
 * A zero boost delegates to `first_order_evolved_variables`.
 */
void first_order_boosted_evolved_variables(
    gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const std::array<double, 3>& boost_velocity, bool centre_advection = true);

/// A rate-resummed extension of the affine-map variables at
/// points `x`: the metric is the inverse of the combination above, Phi its
/// analytic spatial derivative \f$\Phi_{kab} = -(g\, \partial_k G^{-1}\,
/// g)_{ab}\f$, and \f$\partial_t g\f$ the sum of the coefficient drift
/// \f$-(g \sum_A \dot p_A R_A g)\f$ and, unless `centre_advection` is
/// false, the motion of the centre \f$-\dot q^k \Phi_{kab}\f$. Pi is formed
/// with the model's own lapse and shift. Setting `centre_advection` false
/// recovers the static-centre form used before findings 15w; it changes
/// nothing when \f$\dot q^i = 0\f$. This routine is useful for higher-order
/// experiments, but is not the strict Dhesi slow-time first-order model; use
/// `first_order_evolved_variables` for that system.
void evolved_variables(gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
                       gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
                       gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
                       const tnsr::I<DataVector, 3>& x, double mass,
                       const std::array<double, 3>& center,
                       const std::array<double, 13>& p,
                       const std::array<double, 13>& pdot,
                       bool centre_advection = true);

/*!
 * \brief The same evolved variables for a hole moving with constant
 * `boost_velocity`, as an exact Lorentz boost of the static solution.
 *
 * A Lorentz boost is a constant linear coordinate transformation, so
 * \f$\Box x'^\mu = \Lambda^\mu{}_\nu \Box x^\nu = 0\f$: the boosted
 * solution is still exactly harmonic, and still an exact solution of the
 * Einstein equations at every order in \f$v\f$ — unlike prescribing a
 * velocity through the first-order map, which is exact only to
 * \f$O(v)\f$. With \f$M^{\bar a}{}_a = \partial x^{\bar a}/\partial x^a\f$
 * the matrix into the hole's rest frame (`lorentz_boost_matrix(-v)`),
 *
 * \f{align}{
 * g_{ab}(x) &= M^{\bar c}{}_a M^{\bar d}{}_b\, \hat g_{\bar c\bar d}(\bar x),
 * \\
 * \partial_\mu g_{ab}(x) &= M^{\bar e}{}_\mu M^{\bar c}{}_a M^{\bar d}{}_b\,
 *   \partial_{\bar e} \hat g_{\bar c \bar d}(\bar x),
 * \f}
 *
 * evaluated at the rest-frame image \f$\bar x^{\bar a} = M^{\bar a}{}_b
 * (t, x - c)^b\f$. The rest-frame \f$\partial_{\bar e}\hat g\f$ needs only
 * the analytic spatial derivative above plus \f$\partial_t \hat g =
 * \beta^k \Phi_k - \alpha \Pi\f$, so no new derivatives and no finite
 * differences are involved; the lab \f$\Pi\f$ and \f$\Phi\f$ are then
 * rebuilt from the transformed derivatives with the lab lapse and shift.
 *
 * To first order in \f$v\f$ a boost is the map direction
 * \f$\beta_i = \dot q^i = +v_i\f$ with \f$\dot q^0 = \sigma_{ij} = 0\f$;
 * at second order it carries the known \f$\gamma - 1 = v^2/2\f$ clock
 * offset and longitudinal contraction, so the residual error of a
 * first-order fit against this data is a *predicted* \f$O(v^2)\f$ rather
 * than an uncontrolled one. The map parameters are evaluated at the lab
 * time, which is exact for `p = pdot = 0` and only first-order consistent
 * otherwise. `boost_velocity` zero delegates to `evolved_variables`.
 */
void boosted_evolved_variables(
    gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, 3>*> pi,
    gsl::not_null<tnsr::iaa<DataVector, 3>*> phi,
    const tnsr::I<DataVector, 3>& x, double time, double mass,
    const std::array<double, 3>& center, const std::array<double, 13>& p,
    const std::array<double, 13>& pdot,
    const std::array<double, 3>& boost_velocity, bool centre_advection = true);
}  // namespace affine_map_model

/*!
 * \brief Rate-resummed affine-map extension of harmonic Schwarzschild, with
 * the thirteen map parameters supplied as tabulated functions of time.
 *
 * This analytic-solution class retains the legacy exact inversion and the
 * supplied coefficient rates. It is useful for manufactured higher-order
 * experiments and for the exact-boost ground truth at `p = pdot = 0`, but it
 * is not the strict Dhesi slow-time first-order model used by the online value
 * matcher.
 *
 * The model inverse metric in simulation coordinates is
 *
 * \f{align}{
 *   G^{ab}_\text{model}(x, t) = G_0^{ab}(y) + \sum_A p_A(t)\, R_A^{ab}(y),
 *   \qquad y = x - x_\text{center},
 * \f}
 *
 * where \f$G_0\f$ is the harmonic-Schwarzschild inverse metric and the
 * \f$R_A\f$ are the thirteen first-order response columns of the matching
 * scheme, ordered
 *
 * \f{align}{
 *   p = (\dot q^0,\; \beta_x, \beta_y, \beta_z,\; \dot q^x, \dot q^y,
 *        \dot q^z,\; \sigma_{xx}, \sigma_{xy}, \sigma_{xz}, \sigma_{yy},
 *        \sigma_{yz}, \sigma_{zz}).
 * \f}
 *
 * The responses are derivatives at fixed simulation-frame points (the strain
 * columns include the transport term), and the map direction is local
 * harmonic \f$\to\f$ simulation: \f$x_\text{sim} = (1 + \sigma)
 * x_\text{local}\f$. Each off-diagonal strain coefficient multiplies the
 * full symmetric unit perturbation (both index placements), so the six
 * independent components each appear once. This is the convention of the
 * reference implementation
 * `worldtube_matching/responses.py::first_order_response_columns`, which the
 * formulas here transcribe.
 *
 * Provides the generalized-harmonic evolved variables: the spacetime metric
 * is the inverse of \f$G^{ab}\f$; \f$\Phi_{iab}\f$ is its spatial derivative
 * by second-order central differences with step \f$10^{-4}\f$ (matching the
 * reference implementation); \f$\Pi_{ab} = (\beta^k \Phi_{kab} -
 * \partial_t g_{ab})/\alpha\f$ with \f$\partial_t g_{ab} = -\left(g\,
 * \sum_A \dot p_A R_A\, g\right)_{ab}\f$ the coefficient drift at a fixed
 * point. The parameter rates \f$\dot p_A\f$ are supplied as their own table
 * rather than differentiated internally, so a run reproduces exactly the
 * rates used by the offline fit.
 *
 * The parameter tables are interpolated in time with natural cubic splines
 * (linear for two rows, constant for one). Times outside the table are
 * clamped to its ends. With all parameters and rates zero this solution is
 * identical to `HarmonicSchwarzschild`.
 *
 * \note The center is fixed in time: there is no transport term from center
 * motion in \f$\partial_t g\f$. Suitable for a hole at rest (the single-hole
 * test ladder); a moving worldtube needs the center trajectory added.
 *
 * \warning Registered as initial data and as a `WorldtubeTypeD`
 * `AnalyticGaugePrescription` for `EvolveGhSingleBlackHole` only; other code
 * paths that dispatch over `gh::Solutions::all_solutions` do not know this
 * class.
 */
class AffineMappedHarmonicSchwarzschild
    : public virtual evolution::initial_data::InitialData,
      public MarkAsAnalyticSolution {
 public:
  static constexpr size_t volume_dim = 3;
  static constexpr size_t number_of_parameters = 13;

  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of the black hole"};
    static type lower_bound() { return 0.; }
  };
  struct Center {
    using type = std::array<double, volume_dim>;
    static constexpr Options::String help = {
        "The worldtube center, fixed in time"};
  };
  struct ParameterTimes {
    using type = std::vector<double>;
    static constexpr Options::String help = {
        "Strictly increasing times at which the map parameters are "
        "tabulated. A single entry makes the parameters constant. Evaluation "
        "times outside the table are clamped to its ends."};
  };
  struct ParameterValues {
    using type = std::vector<std::array<double, number_of_parameters>>;
    static constexpr Options::String help = {
        "Map parameters at each tabulated time, ordered (qdot0, beta_x, "
        "beta_y, beta_z, qdot_x, qdot_y, qdot_z, sigma_xx, sigma_xy, "
        "sigma_xz, sigma_yy, sigma_yz, sigma_zz)"};
  };
  struct ParameterRates {
    using type = std::vector<std::array<double, number_of_parameters>>;
    static constexpr Options::String help = {
        "HIGHER-ORDER/RESUMMED extension: time derivatives of the map "
        "parameters at each tabulated time, in the same order. Used for the "
        "time derivative of the model metric entering Pi. Supply the rates "
        "the offline fit used; they are not differentiated internally."};
  };
  struct Velocity {
    using type = std::array<double, volume_dim>;
    static constexpr Options::String help = {
        "Constant boost velocity of the hole, applied as an EXACT Lorentz "
        "boost of the (mapped) static solution rather than through the "
        "first-order map. A boost is a constant linear transformation, so "
        "the result is still exactly harmonic and still an exact solution "
        "at every order in v -- which makes it the controlled ground truth "
        "for testing whether the matcher recovers a velocity. Its "
        "first-order content is beta_i = qdot^i = +v_i; the O(v^2) "
        "remainder (gamma - 1 = v^2/2 and the longitudinal contraction) is "
        "analytically known, so a first-order fit's error against this "
        "data is predicted rather than uncontrolled. Zero for a hole at "
        "rest. Note that without a tracking excision the hole walks toward "
        "the boundary at speed v, which bounds the usable run length."};
  };
  using options = tmpl::list<Mass, Center, ParameterTimes, ParameterValues,
                             ParameterRates, Velocity>;
  static constexpr Options::String help = {
      "Rate-resummed affine-map extension of harmonic Schwarzschild, with "
      "the 13 map parameters tabulated in time, optionally boosted exactly. "
      "This class is not the strict slow-time first-order online model. With "
      "all parameters and the velocity zero this is HarmonicSchwarzschild."};

  AffineMappedHarmonicSchwarzschild() = default;
  AffineMappedHarmonicSchwarzschild(
      const AffineMappedHarmonicSchwarzschild& /*rhs*/) = default;
  AffineMappedHarmonicSchwarzschild& operator=(
      const AffineMappedHarmonicSchwarzschild& /*rhs*/) = default;
  AffineMappedHarmonicSchwarzschild(
      AffineMappedHarmonicSchwarzschild&& /*rhs*/) = default;
  AffineMappedHarmonicSchwarzschild& operator=(
      AffineMappedHarmonicSchwarzschild&& /*rhs*/) = default;
  ~AffineMappedHarmonicSchwarzschild() override = default;

  AffineMappedHarmonicSchwarzschild(
      double mass, const std::array<double, volume_dim>& center,
      std::vector<double> parameter_times,
      std::vector<std::array<double, number_of_parameters>> parameter_values,
      std::vector<std::array<double, number_of_parameters>> parameter_rates,
      const std::array<double, volume_dim>& velocity = {{0., 0., 0.}},
      const Options::Context& context = {});

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit AffineMappedHarmonicSchwarzschild(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(AffineMappedHarmonicSchwarzschild);
  /// \endcond

  // The evolved GH variables (what the ghost boundary condition and the
  // analytic-solution observers request) plus the ADM variables that
  // gh::Actions::SetInitialData requests from analytic initial data.
  template <typename DataType>
  using tags =
      tmpl::list<gr::Tags::SpacetimeMetric<DataType, volume_dim>,
                 gh::Tags::Pi<DataType, volume_dim>,
                 gh::Tags::Phi<DataType, volume_dim>, gr::Tags::Lapse<DataType>,
                 gr::Tags::Shift<DataType, volume_dim>,
                 gr::Tags::SpatialMetric<DataType, volume_dim>,
                 gr::Tags::ExtrinsicCurvature<DataType, volume_dim>>;

  using AllVars = tuples::tagged_tuple_from_typelist<tags<DataVector>>;

  template <typename... Tags>
  tuples::TaggedTuple<Tags...> variables(
      const tnsr::I<DataVector, volume_dim>& x, const double time,
      tmpl::list<Tags...> /*meta*/) const {
    AllVars vars = all_variables(x, time);
    return {std::move(get<Tags>(vars))...};
  }

  /// The interpolated map parameters and rates at `time` (clamped to the
  /// table). Exposed for tests and diagnostics.
  void map_parameters(
      gsl::not_null<std::array<double, number_of_parameters>*> values,
      gsl::not_null<std::array<double, number_of_parameters>*> rates,
      double time) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

  friend bool operator==(const AffineMappedHarmonicSchwarzschild& lhs,
                         const AffineMappedHarmonicSchwarzschild& rhs);

 private:
  AllVars all_variables(const tnsr::I<DataVector, volume_dim>& x,
                        double time) const;

  void build_splines();

  double mass_ = 1.0;
  std::array<double, volume_dim> center_{{0., 0., 0.}};
  std::array<double, volume_dim> velocity_{{0., 0., 0.}};
  std::vector<double> parameter_times_{};
  std::vector<std::array<double, number_of_parameters>> parameter_values_{};
  std::vector<std::array<double, number_of_parameters>> parameter_rates_{};
  // built from the tables when there are at least three rows; not pupped,
  // rebuilt after unpacking
  std::vector<intrp::CubicSpline> value_splines_{};
  std::vector<intrp::CubicSpline> rate_splines_{};
};

bool operator!=(const AffineMappedHarmonicSchwarzschild& lhs,
                const AffineMappedHarmonicSchwarzschild& rhs);
}  // namespace gh::Solutions
