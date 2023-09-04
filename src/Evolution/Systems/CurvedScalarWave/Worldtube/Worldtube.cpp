// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/Worldtube.hpp"

namespace CurvedScalarWave::Worldtube {

double worldtube_shrink_factor(const double orbit_radius,
                               const double original_orbit_radius,
                               const double start_shrink_orbit,
                               const double end_shrink_orbit,
                               const double shrink_factor_at_end) {
  const double orbit_radius_fraction = orbit_radius / original_orbit_radius;
  return orbit_radius_fraction * sqrt(orbit_radius_fraction);
}

double worldtube_shrink_factor_derivative(const double orbit_radius,
                                          const double orbit_velocity,
                                          const double original_orbit_radius,
                                          const double start_shrink_orbit,
                                          const double end_shrink_orbit,
                                          const double shrink_factor_at_end) {
  const double orbit_radius_fraction = orbit_radius / original_orbit_radius;
  return 1.5 * sqrt(orbit_radius_fraction) * orbit_velocity /
         original_orbit_radius;
}
}  // namespace CurvedScalarWave::Worldtube
