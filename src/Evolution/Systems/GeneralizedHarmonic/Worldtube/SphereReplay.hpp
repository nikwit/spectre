// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube {
using ReplayTags =
    tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
               gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>;
/// Recorded fields and volume-computed derivatives for a prescribed physical
/// target. The time derivative of K is at fixed inertial coordinates.
struct RecordedPhysicalData {
  Variables<ReplayTags> fields;
  tnsr::iaa<DataVector, 3> d_pi;
  tnsr::ijaa<DataVector, 3> d_phi;
  Scalar<DataVector> kretschmann;
  tnsr::i<DataVector, 3> d_kretschmann;
  Scalar<DataVector> dt_kretschmann;
};
/// Immutable, single-radius GhWorldtube HDF5 tape. Inertial Cartesian tensor
/// components are interpolated in angle and time, never in radius.
/// A synchronized, single-entry cache reuses the last evaluation when time,
/// coordinates and center agree exactly. Geometry is validated on every call.
class SphereReplay {
 public:
  SphereReplay(const std::string& filename, const std::string& subfile,
               double radius, bool read_physical_data = false);
  /// The exact mapped center must be available at the evaluation time and
  /// all six tape times in its interpolation stencil. Validate the trajectory
  /// at those sample times, without interpolating the center itself.
  using CenterTrajectory = std::function<tnsr::I<double, 3>(double)>;
  Variables<ReplayTags> evaluate(double time,
                                 const tnsr::I<DataVector, 3>& coordinates,
                                 const CenterTrajectory& center_at_time) const;
  /// The recording points have fixed inertial angular axes and radius, and
  /// translate with the mapped center. Differentiate the six-point temporal
  /// interpolant at fixed recording angle and subtract v_center^i d_i K.
  RecordedPhysicalData evaluate_physical(
      double time, const tnsr::I<DataVector, 3>& coordinates,
      const CenterTrajectory& center_at_time,
      const tnsr::I<double, 3>& center_velocity) const;
  std::array<double, 2> time_bounds() const;
  double radius() const { return radius_; }

 private:
  struct TimeStencil {
    size_t start;
    std::array<double, 6> weights;
    std::array<double, 6> derivative_weights;
    tnsr::I<double, 3> center;
  };
  TimeStencil time_stencil(double time,
                           const CenterTrajectory& center_at_time) const;
  ylm::Spherepack::InterpolationInfo<DataVector> angular_info(
      const tnsr::I<DataVector, 3>& coordinates,
      const tnsr::I<double, 3>& center) const;
  DataVector interpolate_component(
      size_t component, const TimeStencil& stencil,
      const ylm::Spherepack::InterpolationInfo<DataVector>& info, size_t points,
      bool time_derivative = false) const;
  struct Evaluation {
    double time;
    tnsr::I<DataVector, 3> coordinates;
    tnsr::I<double, 3> center;
    Variables<ReplayTags> fields;
  };
  // Tapes are shared across boundary objects and worker threads. Only cache
  // lookup/publication holds the lock; independent evaluations can proceed.
  // Logical tape contents remain immutable; these members only memoize reads.
  // NOLINTNEXTLINE(spectre-mutable)
  mutable std::mutex evaluation_mutex_;
  // NOLINTNEXTLINE(spectre-mutable)
  mutable std::optional<Evaluation> last_evaluation_;
  double radius_;
  bool read_physical_data_;
  std::unique_ptr<ylm::Spherepack> spherepack_;
  std::vector<double> times_;
  std::vector<std::array<double, 3>> centers_;
  std::vector<std::vector<DataVector>> coefficients_;
};
/// Load during option construction / deserialization, before evolution I/O.
/// All boundary objects in a process share immutable data. The cache is not
/// serialized; migration refers back to the same file, subfile and radius.
std::shared_ptr<const SphereReplay> cached_sphere_replay(
    const std::string& filename, const std::string& subfile, double radius,
    bool read_physical_data = false);
}  // namespace gh::worldtube
