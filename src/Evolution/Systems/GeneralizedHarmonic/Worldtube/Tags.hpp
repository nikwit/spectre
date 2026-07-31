// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/Auto.hpp"
#include "Options/String.hpp"
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

/// Option-created configuration of the online matcher.
struct MatcherConfig {
  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of the black hole"};
    static type lower_bound() { return 0.; }
  };
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "The worldtube center, fixed in time"};
  };
  struct CenterVelocity {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Kinematic pin for the map velocity: qdot^i = (1 + qdot^0) v_center. "
        "Zero for a hole at rest."};
  };
  struct TraceStrainPin {
    using type = double;
    static constexpr Options::String help = {
        "Pin for tr(sigma)/3. Zero for an unperturbed single hole; for a "
        "tidally perturbed hole use the invariant-radius value (findings "
        "§12)."};
  };
  struct FitLMax {
    using type = size_t;
    static constexpr Options::String help = {
        "The fit residual keeps spherical-harmonic modes with l <= FitLMax "
        "of the gauge components. The first-order model space has content "
        "up to l = 4."};
  };
  struct FitInterval {
    using type = double;
    static constexpr Options::String help = {
        "Minimum simulation-time interval between fits. The fit runs at the "
        "first step whose time exceeds the previous fit time by this "
        "amount; between fits the boundary condition extrapolates linearly "
        "with the fitted rates."};
    static type lower_bound() { return 0.; }
  };

  using options = tmpl::list<Mass, Center, CenterVelocity, TraceStrainPin,
                             FitLMax, FitInterval>;
  static constexpr Options::String help = {
      "Online worldtube matching: fit the 13 first-order affine-map "
      "parameters from the evolved fields on the excision sphere."};

  MatcherConfig() = default;
  MatcherConfig(double mass, const std::array<double, 3>& center,
                const std::array<double, 3>& center_velocity,
                double trace_strain_pin, size_t fit_l_max,
                double fit_interval)
      : mass(mass),
        center(center),
        center_velocity(center_velocity),
        trace_strain_pin(trace_strain_pin),
        fit_l_max(fit_l_max),
        fit_interval(fit_interval) {}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  double mass = 1.0;
  std::array<double, 3> center{{0., 0., 0.}};
  std::array<double, 3> center_velocity{{0., 0., 0.}};
  double trace_strain_pin = 0.;
  size_t fit_l_max = 4;
  double fit_interval = 0.;
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
  bool valid = false;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

namespace OptionTags {
struct WorldtubeMatcher {
  using type = Options::Auto<MatcherConfig, Options::AutoLabel::None>;
  static constexpr Options::String help = {
      "Online worldtube matching of the first-order affine map from the "
      "evolved fields on the excision sphere. Set to None to disable."};
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
