// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/CalculateGeodesicPunctureField.hpp"

namespace CurvedScalarWave::Worldtube {

void CalculateGeodesicPunctureField::apply(
    const gsl::not_null<Tags::GeodesicPunctureField<3>::type*> result,
    const gsl::not_null<PunctureFieldWorkspace*> workspace,
    const std::optional<tnsr::I<DataVector, 3>>& centered_face_coordinates,
    const std::array<tnsr::I<double, 3>, 2>& position_velocity,
    const tnsr::I<double, 3>& acceleration, const double charge,
    const PunctureField& puncture) {
  if (not centered_face_coordinates.has_value()) {
    result->reset();
    return;
  }
  if (not result->has_value()) {
    result->emplace();
  }
  puncture.apply_puncture(make_not_null(&result->value()),
                          centered_face_coordinates.value(),
                          position_velocity[0], position_velocity[1],
                          acceleration, workspace.get());
  result->value() *= charge;
}

}  // namespace CurvedScalarWave::Worldtube
