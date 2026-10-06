// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTideFixedPoint.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <pup.h>
#include <utility>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gr::np {
void GeometricFixedPointOptions::pup(PUP::er& p) {
  p | relative_tolerance;
  p | absolute_tolerance;
  p | max_iterations;
  p | damping;
}
bool operator==(const GeometricFixedPointOptions& a,
                const GeometricFixedPointOptions& b) {
  return a.relative_tolerance == b.relative_tolerance and
         a.absolute_tolerance == b.absolute_tolerance and
         a.max_iterations == b.max_iterations and a.damping == b.damping;
}
void GeometricFixedPointDiagnostics::pup(PUP::er& p) {
  p | iterations;
  p | absolute_residual;
  p | relative_residual;
  p | residual_ratio;
}
bool operator==(const GeometricFixedPointDiagnostics& a,
                const GeometricFixedPointDiagnostics& b) {
  return a.iterations == b.iterations and
         a.absolute_residual == b.absolute_residual and
         a.relative_residual == b.relative_residual and
         a.residual_ratio == b.residual_ratio;
}

GeometricFixedPointEvaluation geometric_tide_fixed_point(
    const WeylScalars& measured_psi, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const TriadVector& labels, const DataVector& weights, const double mass,
    const size_t l_max, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_kretschmann,
    const Scalar<DataVector>& dt_kretschmann,
    const GeometricFixedPointOptions& options) {
  if (not(std::isfinite(options.relative_tolerance) and
          options.relative_tolerance >= 0. and
          std::isfinite(options.absolute_tolerance) and
          options.absolute_tolerance >= 0. and
          (options.relative_tolerance > 0. or
           options.absolute_tolerance > 0.) and
          options.max_iterations > 0 and std::isfinite(options.damping) and
          options.damping > 0. and options.damping <= 1.)) {
    ERROR("Invalid geometric Psi0 fixed-point settings");
  }
  const size_t np = weights.size();
  if (np == 0 or measured_psi.get(1).size() != np or get(lapse).size() != np or
      get(dt_kretschmann).size() != np or
      not(std::isfinite(mass) and mass > 0.)) {
    ERROR("Invalid geometric Psi0 fixed-point face data");
  }
  double weight_sum = 0.;
  for (const double weight : weights) {
    if (not(std::isfinite(weight) and weight > 0.)) {
      ERROR("Geometric Psi0 fixed point needs positive finite weights");
    }
    weight_sum += weight;
  }
  if (not std::isfinite(weight_sum)) {
    ERROR("Geometric Psi0 fixed point requires finite total quadrature weight");
  }
  const auto rms = [&weights, weight_sum](const ComplexDataVector& values) {
    double norm = 0.;
    for (size_t p = 0; p < weights.size(); ++p) {
      norm =
          std::hypot(norm, sqrt(weights[p] / weight_sum) * std::abs(values[p]));
    }
    if (not std::isfinite(norm)) {
      ERROR("Nonfinite geometric Psi0 fixed-point residual or target");
    }
    return norm;
  };
  // Copy only the known slots. Even a NaN in measured Psi0 is irrelevant.
  WeylScalars working(np, std::complex<double>{0., 0.});
  for (size_t a = 1; a < 5; ++a) {
    working.get(a) = measured_psi.get(a);
    (void)rms(working.get(a));
  }
  GeometricFixedPointEvaluation result{};
  double previous_residual = 0.;
  for (size_t iteration = 1; iteration <= options.max_iterations; ++iteration) {
    result.registration = register_frame(working, adapted_rotation, mass);
    result.rapidity = invariant_rapidity(
        result.registration.member, adapted_rotation, spatial_metric, lapse,
        shift, d_kretschmann, dt_kretschmann);
    result.geometric = evaluate_geometric_second_order(
        result.registration, result.rapidity, adapted_rotation, spatial_metric,
        labels, weights, mass, l_max);
    const auto& target = get(result.geometric.second_order.psi0_target);
    const double residual = rms(target - working.get(0));
    const double scale = std::max(rms(target), rms(working.get(0)));
    result.diagnostics = {
        iteration, residual, scale > 0. ? residual / scale : 0.,
        previous_residual > 0. ? residual / previous_residual : 0.};
    if (residual <=
        options.absolute_tolerance + options.relative_tolerance * scale) {
      return result;
    }
    previous_residual = residual;
    working.get(0) += options.damping * (target - working.get(0));
  }
  ERROR("Geometric Psi0 fixed point failed to converge after "
        << result.diagnostics.iterations << " evaluations; absolute residual "
        << result.diagnostics.absolute_residual << ", relative residual "
        << result.diagnostics.relative_residual << ", residual ratio "
        << result.diagnostics.residual_ratio
        << ". No unconverged target is imposed.");
}
}  // namespace gr::np
