// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "NumericalAlgorithms/Strahlkorper/Strahlkorper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace control_system::measurements {

inline DataVector gauss_bonnet_center(
    const ylm::Strahlkorper<Frame::Grid>& strahlkorper,
    const Scalar<DataVector>& gauss_bonnet_scalar) {
  ASSERT(strahlkorper.l_max() > 0 and strahlkorper.m_max() > 0,
         "Need l_max >= 1 and m_max >= 1 to compute a dipole.");
  const auto& ylm = strahlkorper.ylm_spherepack();
  const size_t points_per_sphere = ylm.physical_size();
  const DataVector& all_values = get(gauss_bonnet_scalar);
  ASSERT(all_values.size() >= points_per_sphere and
             all_values.size() % points_per_sphere == 0,
         "Unexpected number of Gauss-Bonnet points: "
             << all_values.size() << ", expected a multiple of "
             << points_per_sphere << ".");

  // If multiple radii are observed, use the outermost sphere.
  const size_t offset = all_values.size() - points_per_sphere;
  DataVector values_on_control_sphere{points_per_sphere};
  for (size_t i = 0; i < points_per_sphere; ++i) {
    values_on_control_sphere[i] = sqrt(all_values[offset + i] / 16.0 / 3.0);
  }

  const DataVector gb_coefs = ylm.phys_to_spec(values_on_control_sphere);
  ylm::SpherepackIterator iterator(strahlkorper.l_max(), strahlkorper.m_max());

  const double l0_coef = gb_coefs[iterator.set(0, 0)()] * sqrt(M_PI / 2.0);
  const std::array<double, 3> l1_coefs{
      gb_coefs[iterator.set(1, 1)()] * sqrt(M_PI),
      -gb_coefs[iterator.set(1, -1)()] * sqrt(M_PI),
      gb_coefs[iterator.set(1, 0)()] * sqrt(M_PI / 2.0)};

  const double dipole_magnitude =
      sqrt(l1_coefs[0] * l1_coefs[0] + l1_coefs[1] * l1_coefs[1] +
           l1_coefs[2] * l1_coefs[2]);
  const auto& sphere_center = strahlkorper.physical_center();
  DataVector center{sphere_center[0], sphere_center[1], sphere_center[2]};

  if (dipole_magnitude == 0.0) {
    return center;
  }

  ASSERT(l0_coef != 0.0,
         "Gauss-Bonnet monopole coefficient is zero, so the control "
         "error normalization is singular.");

  const double radius = strahlkorper.average_radius();
  const double delta = dipole_magnitude * radius / (sqrt(3.0) * l0_coef);
  center[0] += delta * l1_coefs[0] / dipole_magnitude;
  center[1] += delta * l1_coefs[1] / dipole_magnitude;
  center[2] += delta * l1_coefs[2] / dipole_magnitude;
  return center;
}

}  // namespace control_system::measurements
