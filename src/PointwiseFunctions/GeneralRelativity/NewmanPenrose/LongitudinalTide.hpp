// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"

namespace gr::np {
/// Exact type-D completion of slots 1,2,3. Slots 0 and 4 are never read.
struct LongitudinalSeed {
  WeylScalars psi;
  Scalar<ComplexDataVector> coulomb;
  double minimum_relative_denominator{1.};
};
LongitudinalSeed longitudinal_seed(const WeylScalars& measured);
/// Analytic differential of the seed, including its missing transverse slots.
WeylScalars longitudinal_seed_response(const LongitudinalSeed& seed,
                                       const WeylScalars& column);
struct LongitudinalTideEvaluation {
  FrameRegistration registration;
  Scalar<DataVector> rapidity;
  GeometricFrame frame;
  TidalMoments moments;
  WeylScalars completed;
  double relative_residual{};
  double normalized_condition_number{};
  double minimum_relative_denominator{};
  double maximum_gradient_tilt{};
};
/// Single real ten-column fit of the corrected gradient-vector and imaginary
/// Coulomb sensors. Implements longitudinal_response.tex, including dS[T].
/// Only the completed slot zero is an incoming boundary target. The outgoing
/// field is not prescribed. All four gradient components are live inputs;
/// independence from transverse scalar values holds with that gradient fixed.
LongitudinalTideEvaluation evaluate_longitudinal_tide(
    const WeylScalars& measured, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const TriadVector& outward_labels, const DataVector& label_weights,
    double mass, size_t l_max, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& spatial_gradient,
    const Scalar<DataVector>& time_derivative);
}  // namespace gr::np
