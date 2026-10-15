// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <memory>
#include <string>
#include <variant>
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
class er;
}
namespace gh::BoundaryConditions::detail {
struct MomentFeedbackReference;

/// Reference-assisted diagnostic, restricted to legacy Quadrupole.
/// The reference is an independently evolved prescribed-moment control.
/// Align its Cholesky axes to the live axes by the weighted mean of the
/// proper polar rotations of L_live^{-1} L_control. No boost is transported.
class MomentFeedback {
 public:
  struct Strength {
    using type = std::variant<double, std::array<double, 2>>;
    static constexpr Options::String help =
        "Feedback strength in [0,1]: a scalar for both sectors, or "
        "[Electric, Magnetic] for independent fitted-moment gains.";
  };
  struct StartTime {
    using type = double;
    static constexpr Options::String help =
        "Simulation time to begin the feedback ramp.";
  };
  struct RampDuration {
    using type = double;
    static constexpr Options::String help =
        "Positive duration of the sin-squared ramp.";
  };
  struct MomentsFile {
    using type = std::string;
    static constexpr Options::String help =
        "Control CSV: time,E11,E22,E12,E13,E23,B11,B22,B12,B13,B23.";
  };
  struct GeometryFile {
    using type = std::string;
    static constexpr Options::String help =
        "NPFBG001/2 control geometry; version 2 includes the exact coordinate "
        "map.";
  };
  using options =
      tmpl::list<Strength, StartTime, RampDuration, MomentsFile, GeometryFile>;
  static constexpr Options::String help =
      "Diagnostic Q_recorded + alpha_E(t) Re(delta Q) + i alpha_B(t) "
      "Im(delta Q) + R pulse R^T, where delta Q=Q_live-R Q_control R^T. "
      "Recorded coefficients retain the existing prescribed-control "
      "convention. "
      "Requires identical moving face geometry; rejects time extrapolation.";
  MomentFeedback() = default;
  MomentFeedback(Strength::type strength, double start, double duration,
                 std::string moments_file, std::string geometry_file,
                 const Options::Context& context = {});
  /// Maximum of the two gains; retains the scalar diagnostic convention.
  double strength(double time) const;
  /// Electric and magnetic gains, after the common ramp.
  std::array<double, 2> strengths(double time) const;
  gr::np::TidalMoments imposed(double time, const gr::np::TidalMoments& live,
                               const gr::np::TidalMoments& recorded,
                               const gr::np::TidalMoments& pulse,
                               const tnsr::ii<DataVector, 3>& metric,
                               const tnsr::I<DataVector, 3>& coordinates,
                               const DataVector& weights) const;
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  void load();
  std::array<double, 2> strengths_{{0., 0.}};
  double start_{0.}, duration_{1.};
  std::string moments_file_, geometry_file_;
  std::shared_ptr<const MomentFeedbackReference> reference_;
};
}  // namespace gh::BoundaryConditions::detail
