// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"

namespace CurvedScalarWave::Worldtube {

/*!
 * \brief Evaluate the geodesic puncture once at the current element stage.
 *
 * Mutates an element-local workspace as well as the field. This action runs
 * before SendToWorldtube, outside the self-force iteration loop. The same
 * workspace is reused by IteratePunctureField, while the geodesic field is
 * retained until the stage's boundary data have been received.
 */
struct CalculateGeodesicPunctureField {
  using return_tags =
      tmpl::list<Tags::GeodesicPunctureField<3>, Tags::PunctureFieldScratch>;
  using simple_tags = return_tags;
  using argument_tags =
      tmpl::list<Tags::FaceCoordinates<3, Frame::Inertial, true>,
                 Tags::ParticlePositionVelocity<3>,
                 Tags::GeodesicAcceleration<3>, Tags::Charge,
                 Tags::PunctureFieldConfig>;

  static void apply(
      gsl::not_null<Tags::GeodesicPunctureField<3>::type*> result,
      gsl::not_null<PunctureFieldWorkspace*> workspace,
      const std::optional<tnsr::I<DataVector, 3>>& centered_face_coordinates,
      const std::array<tnsr::I<double, 3>, 2>& position_velocity,
      const tnsr::I<double, 3>& acceleration, double charge,
      const PunctureField& puncture);
};

}  // namespace CurvedScalarWave::Worldtube
