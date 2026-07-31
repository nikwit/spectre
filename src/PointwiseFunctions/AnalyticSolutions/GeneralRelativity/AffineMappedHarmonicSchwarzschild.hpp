// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "DataStructures/DataVector.hpp"
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
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gh::Solutions {
/*!
 * \brief Harmonic Schwarzschild pushed through the first-order affine
 * worldtube-matching map, with the thirteen map parameters supplied as
 * tabulated functions of time.
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
    using type =
        std::vector<std::array<double, number_of_parameters>>;
    static constexpr Options::String help = {
        "Map parameters at each tabulated time, ordered (qdot0, beta_x, "
        "beta_y, beta_z, qdot_x, qdot_y, qdot_z, sigma_xx, sigma_xy, "
        "sigma_xz, sigma_yy, sigma_yz, sigma_zz)"};
  };
  struct ParameterRates {
    using type =
        std::vector<std::array<double, number_of_parameters>>;
    static constexpr Options::String help = {
        "Time derivatives of the map parameters at each tabulated time, in "
        "the same order. Used for the time derivative of the model metric "
        "entering Pi. Supply the rates the offline fit used; they are not "
        "differentiated internally."};
  };
  using options =
      tmpl::list<Mass, Center, ParameterTimes, ParameterValues,
                 ParameterRates>;
  static constexpr Options::String help = {
      "Harmonic Schwarzschild pushed through the first-order affine "
      "worldtube-matching map, with the 13 map parameters tabulated in "
      "time. With all parameters zero this is HarmonicSchwarzschild."};

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
                 gh::Tags::Phi<DataType, volume_dim>,
                 gr::Tags::Lapse<DataType>,
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
