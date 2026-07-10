// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/IterateAccelerationTerms.hpp"

#include <cstddef>
#include <optional>
#include <tuple>

#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {
void IterateAccelerationTerms::apply(
    const gsl::not_null<Scalar<DataVector>*> /*acceleration_terms*/,
    const std::array<tnsr::I<double, Dim>, 2>& /*pos_vel*/,
    const tuples::TaggedTuple<
        gr::Tags::SpacetimeMetric<double, Dim>,
        gr::Tags::InverseSpacetimeMetric<double, Dim>,
        gr::Tags::SpacetimeChristoffelSecondKind<double, Dim>,
        gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim>,
        Tags::TimeDilationFactor>& /*background*/,
    const tnsr::I<double, Dim, Frame::Inertial>& /*geodesic_acc*/,
    const Scalar<double>& /*psi_monopole*/,
    const Scalar<double>& /*dt_psi_monopole*/,
    const tnsr::i<double, Dim, Frame::Grid>& /*psi_dipole*/,
    const tnsr::i<double, Dim, Frame::Grid>& /*dt_psi_dipole*/,
    const double /*charge*/, const std::optional<double> /*mass*/,
    const double /*time*/, const std::optional<double> /*turn_on_time*/,
    const std::optional<double> /*turn_on_interval*/,
    const size_t /*iteration*/) {
  // The iterative scheme computes the acceleration terms of the puncture
  // field due to the scalar self-force in the inertial frame. The grid-frame
  // worldtube scheme holds the regular field coefficients in the co-moving
  // grid frame, so the two cannot be combined.
  ERROR(
      "The scalar self-force is not supported by the grid-frame worldtube "
      "scheme. Select `SelfForce: None`.");
}
}  // namespace CurvedScalarWave::Worldtube
