// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <cstddef>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"

namespace PUP {
class er;
}

namespace gr::np {
struct GeometricFixedPointOptions {
  double relative_tolerance{1.e-10};
  /// Absolute tolerance in curvature units, on the label-weighted RMS.
  double absolute_tolerance{1.e-14};
  size_t max_iterations{20};
  double damping{1.};
  void pup(PUP::er& p);
};
bool operator==(const GeometricFixedPointOptions& a,
                const GeometricFixedPointOptions& b);

struct GeometricFixedPointDiagnostics {
  size_t iterations{};
  double absolute_residual{};
  double relative_residual{};
  /// Last / previous undamped residual norm; zero on the first evaluation.
  double residual_ratio{};
  void pup(PUP::er& p);
};
bool operator==(const GeometricFixedPointDiagnostics& a,
                const GeometricFixedPointDiagnostics& b);

struct GeometricFixedPointEvaluation {
  FrameRegistration registration;
  Scalar<DataVector> rapidity;
  GeometricTideEvaluation geometric;
  GeometricFixedPointDiagnostics diagnostics;
};

/// Solve z = F(z; Psi1..4, metric, dK, dtK) with the geometric quadrupole
/// evaluator. Starts from zero on every call, without reading measured Psi0.
/// Every iteration refits the registration, invariant boost, map and moments.
/// The measured curvature gradients stay fixed: this removes the direct
/// algebraic Psi0 dependence, not all dependence on the evolved interior.
/// Stops on ||F(z)-z|| <= atol + rtol max(||F(z)||, ||z||), using fixed
/// label-weighted RMS norms (not iteration-dependent map weights). The
/// residual is undamped, including when damping < 1. Returns F(z) and its
/// matching frame/fit; errors on nonconvergence or invalid/nonfinite data.
/// No evolution or relaxation history is changed by this function.
GeometricFixedPointEvaluation geometric_tide_fixed_point(
    const WeylScalars& measured_psi, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const TriadVector& labels, const DataVector& weights, double mass,
    size_t l_max, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_kretschmann,
    const Scalar<DataVector>& dt_kretschmann,
    const GeometricFixedPointOptions& options);
}  // namespace gr::np
