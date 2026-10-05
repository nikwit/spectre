// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "Utilities/Gsl.hpp"

namespace PUP {
class er;
}

namespace gr::np {
/// Real STF components: E2(5), B2(5), E3(7), B3(7).
/// Rank two uses xx, yy, xy, xz, yz, with zz = -xx-yy.
/// Rank three uses xxx, xxy, xxz, xyy, xyz, yyy, yyz; symmetry and
/// vanishing traces determine the other components. Magnetic octupoles
/// are physical B_ijk (the 4/3 factor belongs to the radial profile).
using ThirdOrderMoments = std::array<double, 24>;
using DottedTidalMoments = std::array<double, 10>;
using ThirdOrderColumns = std::array<std::array<ComplexDataVector, 2>, 34>;

/// Columns in the registered (not rest-boosted) Psi0/Psi4 slots:
/// E2, B2, E3, B3, Edot, Bdot. Dotted columns contain induction,
/// horizon-calibrated near-zone response and delta_t times quadrupole.
ThirdOrderColumns third_order_tide_columns(const GeometricFrame& frame,
                                           const Scalar<DataVector>& radius,
                                           double mass,
                                           const DataVector& slice_tilt);

struct ThirdOrderFit {
  ThirdOrderMoments moments{};
  Scalar<ComplexDataVector> psi0_target;
  double relative_residual{};
  double condition_number{};
};

/// Fit only undotted moments to measured Psi4 with fixed dotted moments.
/// Psi1..4 are retained when replacing Psi0 and pushing to the NR tetrad.
ThirdOrderFit fit_third_order_tide(const FrameRegistration& registration,
                                   const ThirdOrderColumns& columns,
                                   const DataVector& weights,
                                   const DottedTidalMoments& dots);

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

struct GeometricTideSample {
  double time{};
  DottedTidalMoments quadrupole{};
  TriadVector direction;
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
bool operator==(const GeometricTideSample& a, const GeometricTideSample& b);

/// Only full-step RHS states are retained, at most five. Substeps may use
/// history but do not enter it. A backward step or a resolution change clears
/// the history. Same-time calls replace a sample and never duplicate it.
struct GeometricTideHistory {
  std::vector<GeometricTideSample> samples;
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
bool operator==(const GeometricTideHistory& a, const GeometricTideHistory& b);

struct CausalTidalDerivative {
  DottedTidalMoments dots{};
  std::array<double, 3> angular_velocity{};
  size_t derivative_order{};
};

/// Backward polynomial derivative through current and up to four strictly
/// past full-step samples, for nonuniform timesteps. Remove [Omega,H] with
/// Omega measured from dN/dT + V^A d_A N; then divide by the model clock.
/// At startup return zero dots until a past sample exists. step_start is
/// the start of the current integration step, not an intermediate RK time.
CausalTidalDerivative causal_tidal_derivative(
    gsl::not_null<GeometricTideHistory*> history,
    const GeometricTideSample& current, const EigenSphereMap& map,
    const GeometricTimeData& time_data, double step_start, bool full_step);
}  // namespace gr::np
