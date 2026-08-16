// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "Time/History.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

/// Online first-order worldtube matching for the generalized harmonic
/// system: fit the affine-map parameters from the evolved fields on the
/// excision sphere as the run progresses, to drive the ghost gauge boundary
/// condition (findings §15b/§15c prescription, made self-consistent).
namespace gh::Worldtube {

/// Number of first-order map parameters
/// (qdot0, beta_i, qdot^i, sigma_ij).
static constexpr size_t num_map_parameters = 13;

/// Number of first-order affine-rate profiles in the order-by-order model.
static constexpr size_t num_order_one_rates = 16;

/// How the adopted q8 order-one affine-rate model is used online.
enum class OrderOneMode {
  /// Do not fit or apply the order-one rates.
  Off,
  /// Fit and observe the rates, but keep the boundary frame-only.
  Shadow,
  /// Fit the rates and add their strict linear response to the boundary model.
  Apply,
};

/// Option-created configuration of the online matcher.
struct MatcherConfig {
  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of the black hole"};
    static type lower_bound() { return 0.; }
  };
  struct FitInterval {
    using type = double;
    static constexpr Options::String help = {
        "Minimum simulation-time interval between fits. The fit runs at the "
        "first step whose time exceeds the previous fit time by this "
        "amount. Between fits the exact frame is held fixed while the model "
        "center advances with its fitted coordinate velocity."};
    static type lower_bound() { return 0.; }
  };
  struct OrderOne {
    using type = OrderOneMode;
    static constexpr Options::String help = {
        "Use of the q8-selected order-one affine-rate model. 'Off' runs only "
        "the clean-sector order-zero exact-frame fit. 'Shadow' also solves "
        "the 13 identifiable rates from ell=0 of D_T u+ stacked with "
        "D_R(D_T u+), pins rotation to zero, performs two held-out-gated "
        "feedback iterations, and observes the rates without applying them. "
        "'Apply' runs the identical fit and adds its strict linear rate "
        "response to the exact-frame model supplied to the boundary whenever "
        "the rate solve is valid. The held-out u- result remains diagnostic."};
  };
  struct ExcisionSphereName {
    using type = std::string;
    static constexpr Options::String help = {
        "Name of the excision sphere in the domain whose abutting block "
        "hosts the worldtube boundary element the matcher reads. "
        "'ExcisionSphere' for the single-hole Sphere/SphericalShells "
        "domains; 'ExcisionSphereB' for the small hole in a "
        "BinaryCompactObject domain."};
  };

  using options = tmpl::list<Mass, FitInterval, OrderOne, ExcisionSphereName>;
  static constexpr Options::String help = {
      "Online exact-frame worldtube matching. Fits the 13-parameter finite "
      "frame map to outgoing u+ at the excision face using l<=4 and the "
      "element-local radial stencil. It can run the q8-selected projected "
      "order-zero/order-one alternation either as a shadow diagnostic or as "
      "the strict first-order boundary model. "
      "Historical and A/B controls are deliberately fixed to the q8-validated "
      "production choices."};

  MatcherConfig() = default;
  MatcherConfig(double mass, double fit_interval, OrderOneMode order_one,
                std::string excision_sphere_name)
      : mass(mass),
        fit_interval(fit_interval),
        fit_uplus(true),
        fit_exact_frame(true),
        fit_radial_derivative(true),
        order_one(order_one),
        excision_sphere_name(std::move(excision_sphere_name)) {}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  double mass = 1.0;
  // Internal implementation and regression-test controls below are fixed by
  // the option constructor above; they are not input-file options.
  std::array<double, 3> center{{0., 0., 0.}};
  std::array<double, 3> center_velocity{{0., 0., 0.}};
  double trace_strain_pin = 0.;
  size_t fit_l_max = 4;
  double fit_interval = 0.;
  bool fit_center_offset = false;
  bool rate_ode = false;
  bool second_order_ode = false;
  bool stepper_ode = false;
  double gauge_damping = 0.;
  double uplus_anchor = 0.;
  bool fit_uplus = false;
  bool kretschmann_trace_pin = false;
  double trace_pin_interval = 0.5;
  bool fit_trace_strain = false;
  bool fit_velocity = false;
  bool fit_bulk_boost = false;
  bool fit_exact_frame = false;
  bool fit_velocity_separately = false;
  bool pin_symmetric_factor = false;
  bool fit_radial_derivative = false;
  double radial_derivative_weight = 1.0;
  OrderOneMode order_one = OrderOneMode::Off;
  bool centre_advection = true;
  double spatial_monopole_weight = 1.;
  std::array<double, 15> uplus_block_weights{
      {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}};
  size_t fit_radial_index = 0;
  std::string excision_sphere_name{"ExcisionSphere"};
};

bool operator==(const MatcherConfig& lhs, const MatcherConfig& rhs);
bool operator!=(const MatcherConfig& lhs, const MatcherConfig& rhs);

/// The per-element state of the online matcher: the latest fitted
/// parameters, their estimated rates, and the previous fit for the
/// backward-difference drives. Meaningful only on the element owning the
/// excision face; default (invalid) everywhere else.
struct MapParameterData {
  double last_fit_time = std::numeric_limits<double>::lowest();
  double previous_fit_time = std::numeric_limits<double>::lowest();
  std::array<double, num_map_parameters> p{};
  std::array<double, num_map_parameters> p_previous{};
  std::array<double, num_map_parameters> pdot{};
  std::array<double, num_map_parameters> pddot{};
  /// Finite, nonperturbative Lorentz velocity of the Schwarzschild
  /// background. This is separate from the epsilon-order affine residual.
  std::array<double, 3> bulk_velocity{};
  /// Instantaneous inertial center of the worldtube sampling sphere, measured
  /// from the l=0 modes of its inertial collocation coordinates. This is
  /// numerical-domain geometry, not a black-hole-center measurement.
  std::array<double, 3> worldtube_center{};
  /// Worldtube center at `last_fit_time`. Together with `worldtube_center`,
  /// this removes the domain-map displacement from the advected fit offset.
  std::array<double, 3> worldtube_center_at_last_fit{};
  bool worldtube_center_valid = false;
  /// Black-hole center minus `worldtube_center_at_last_fit`. Thus this stays
  /// small when the control system tracks the hole; it is not the absolute
  /// inertial black-hole position.
  std::array<double, 3> center_offset{};
  /// Zeroth-order exact-frame state (FitExactFrame mode): the fitted frame
  /// parameters ordered (rapidity[3], s0, sigma[3], s_(ij)[6]) and the
  /// derived coordinate centre velocity V_c = L^i_0/L^0_0, which advects the
  /// model centre between fits (spec Eq. Z4).
  std::array<double, num_map_parameters> exact_frame_theta{};
  std::array<double, 3> exact_frame_center_velocity{};
  bool exact_frame_valid = false;
  /// Adopted order-one affine rates. Persisted for restart continuity and
  /// diagnostics; consumed by the boundary in `OrderOneMode::Apply` whenever
  /// the rate solve is valid.
  std::array<double, num_order_one_rates> order_one_rates{};
  bool order_one_valid = false;
  /// Stepper-integrated mode: the 26-component state (p, pdot) and its
  /// time-stepper history
  DataVector ode_state{};
  TimeSteppers::History<DataVector> ode_history{};
  /// Stepper-integrated mode: the state at the start of the current step
  /// and that step's substep-0 id, so a rejected step restarts from the
  /// correct value instead of the rejected end-of-step value
  DataVector ode_step_start{};
  TimeStepId ode_step_id{};
  /// Kretschmann trace pin: the held value and its measurement time
  double trace_pin_value = 0.;
  double trace_pin_time = std::numeric_limits<double>::lowest();
  /// Open-loop centre measurement from the l = 1 content of the
  /// Kretschmann radius on the excision sphere, with the previous sample so
  /// its backward difference can be logged as the hole--worldtube relative
  /// velocity. Adding the independently logged worldtube-center velocity
  /// gives an inertial velocity estimate. A diagnostic only: never an input
  /// to the solve (findings 15t).
  std::array<double, 3> gb_dipole{};
  std::array<double, 3> gb_dipole_previous{};
  std::array<double, 3> gb_dipole_velocity{};
  double gb_dipole_time = std::numeric_limits<double>::lowest();
  double gb_dipole_time_previous = std::numeric_limits<double>::lowest();
  /// u^+ anchor state: the latest and previous value fits of the map to
  /// the gauge projection of the outgoing characteristic, their times,
  /// and the warm-start vector of the Gauss-Newton solve
  std::array<double, num_map_parameters> anchor_p{};
  std::array<double, num_map_parameters> anchor_p_previous{};
  double anchor_time = std::numeric_limits<double>::lowest();
  double anchor_time_previous = std::numeric_limits<double>::lowest();
  bool anchor_valid = false;
  bool valid = false;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

namespace OptionTags {
struct WorldtubeMatcher {
  using type = Options::Auto<MatcherConfig, Options::AutoLabel::None>;
  static constexpr Options::String help = {
      "Online exact-frame matching from the evolved fields on the excision "
      "sphere, with optional shadow or applied first-order affine rates. Set "
      "to None to disable all matching."};
};
}  // namespace OptionTags

namespace Tags {
/// Global-cache tag holding the matcher configuration (nullopt = disabled).
struct Matcher : db::SimpleTag {
  using type = std::optional<MatcherConfig>;
  using option_tags = tmpl::list<OptionTags::WorldtubeMatcher>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// DataBox tag with the latest fitted map parameters (see
/// `MapParameterData`).
struct MapParameters : db::SimpleTag {
  using type = MapParameterData;
};
}  // namespace Tags

namespace Initialization {
/// Initializes the matcher state to its default (invalid) value.
struct InitializeMapParameters {
  using return_tags = tmpl::list<Tags::MapParameters>;
  using argument_tags = tmpl::list<>;
  using simple_tags = tmpl::list<Tags::MapParameters>;
  using compute_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  static void apply(const gsl::not_null<MapParameterData*> data) {
    *data = MapParameterData{};
  }
};
}  // namespace Initialization
}  // namespace gh::Worldtube

template <>
struct Options::create_from_yaml<gh::Worldtube::OrderOneMode> {
  template <typename Metavariables>
  static gh::Worldtube::OrderOneMode create(const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
gh::Worldtube::OrderOneMode
Options::create_from_yaml<gh::Worldtube::OrderOneMode>::create<void>(
    const Options::Option& options);
