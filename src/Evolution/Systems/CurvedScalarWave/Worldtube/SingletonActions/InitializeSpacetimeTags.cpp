// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/InitializeSpacetimeTags.hpp"

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube::Initialization {

void InitializeSpacetimeTags::apply(
    const gsl::not_null<tnsr::AA<double, Dim, Frame::Grid>*>
        inverse_spacetime_metric,
    const gsl::not_null<tnsr::A<double, Dim, Frame::Grid>*>
        trace_spacetime_christoffel,
    const gsl::not_null<double*> expiration_time,
    const gsl::not_null<std::array<double, 2>*> worldtube_radius_and_velocity,
    const ExcisionSphere<Dim>& excision_sphere,
    const std::array<double, 3>& envelope_and_object_radii) {
  const double M = 1.;
  const double orbit_radius = get(magnitude(excision_sphere.center()));
  *inverse_spacetime_metric = tnsr::AA<double, Dim, Frame::Grid>(0.);
  get<0, 0>(*inverse_spacetime_metric) = -1. - 2. * M / orbit_radius;
  get<1, 1>(*inverse_spacetime_metric) = 1. - 2. * M / orbit_radius;
  get<2, 2>(*inverse_spacetime_metric) =
      (-2. * M + square(orbit_radius) - orbit_radius) / square(orbit_radius);
  get<3, 3>(*inverse_spacetime_metric) = 1.;
  get<0, 1>(*inverse_spacetime_metric) = 2. * M / orbit_radius;
  get<2, 0>(*inverse_spacetime_metric) =
      (2. * M + orbit_radius) / (orbit_radius * sqrt(orbit_radius));
  get<2, 1>(*inverse_spacetime_metric) =
      -2. * M / (orbit_radius * sqrt(orbit_radius));

  get<0>(*trace_spacetime_christoffel) = -2. * M / square(orbit_radius);
  get<1>(*trace_spacetime_christoffel) =
      -(2. * M + orbit_radius - 2. * M * orbit_radius) / cube(orbit_radius);
  get<2>(*trace_spacetime_christoffel) =
      6. * M / (square(orbit_radius) * sqrt(orbit_radius));
  get<3>(*trace_spacetime_christoffel) = 0.;
  worldtube_radius_and_velocity->at(0) = envelope_and_object_radii.at(1);
  worldtube_radius_and_velocity->at(1) = 0.;

  *expiration_time = 1e-8;
}
}  // namespace CurvedScalarWave::Worldtube::Initialization
