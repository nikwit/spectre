// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <vector>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "Utilities/Gsl.hpp"

namespace PUP {
class er;
}

namespace gr::np {
struct GeometricTimeData {
  /// Angular shift along the no-screen flow, in the label tangent basis.
  std::array<DataVector, 2> flow;
  DataVector slice_tilt;
  double clock_rate{};
  double minimum_clock_rate{};
  double tilt_residual{};
};

/// Reconstruct tau = -u_flat/sqrt(f). Integrate its tangential part in
/// l<=2 scalar harmonics (the offline linear/quadratic potential), fixing
/// its mean with geometric weights. Clock rate follows the no-screen flow.
GeometricTimeData geometric_time_data(
    const GeometricFrame& frame, const TriadVector& labels,
    double coordinate_radius, const Scalar<DataVector>& measured_radius,
    double mass, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Scalar<DataVector>& lapse, const TriadVector& shift,
    const TriadVector& mesh_velocity);

/// Rigid part of dN/dT + V.dN; removes grid advection from axis motion.
std::array<double, 3> geometric_rotation_rate(
    const EigenSphereMap& map, const TriadVector& direction_dt,
    const std::array<DataVector, 2>& flow);

/// One full-step filter anchor. The STF components belong to direction's axes.
struct GeometricRelaxationSample {
  double time{};
  double clock_rate{};
  TriadVector direction;
  DataVector weights;
  std::array<double, 3> flow_rotation{};
  TidalMoments raw{};
  TidalMoments filtered{};
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
bool operator==(const GeometricRelaxationSample& a,
                const GeometricRelaxationSample& b);

/// Six full-step anchors, serialized for migration/restart. RHS stages do not
/// mutate them. Rollback replaces newer anchors; rollback before retained
/// history or a changed angular resolution reinitializes from the raw fit.
struct GeometricRelaxationHistory {
  std::vector<GeometricRelaxationSample> samples;
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
bool operator==(const GeometricRelaxationHistory& a,
                const GeometricRelaxationHistory& b);

/// Exponential relaxation with linear forcing in a common, transported frame.
/// Finite map alignment is corrected by the no-screen flow at both endpoints.
/// tau is in NR time unless model_time is true; the latter uses the trapezoidal
/// geometric clock integral. Smooth forcing/transport is second-order in time.
/// full_step means TimeStepId::substep()==0; step_start is its step_time.
/// The first evaluation initializes from raw. Same-time full-step calls are
/// recomputed from strictly earlier anchors, never repeatedly relaxed.
TidalMoments relax_geometric_tidal_moments(
    gsl::not_null<GeometricRelaxationHistory*> history, const TidalMoments& raw,
    const EigenSphereMap& map, const GeometricTimeData& temporal, double tau,
    bool model_time, double time, double step_start, bool full_step);
}  // namespace gr::np
