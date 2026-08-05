// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>

#include "ControlSystem/Measurements/GaussBonnetCenter.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "NumericalAlgorithms/Strahlkorper/Strahlkorper.hpp"

namespace {

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.GaussBonnetCenterTracking",
    "[Unit][Evolution]") {
  const double radius = 2.0;
  const std::array<double, 3> sphere_center{{-1.0, 0.5, 0.25}};
  const ylm::Strahlkorper<Frame::Grid> strahlkorper_shifted{2_st, radius,
                                                            sphere_center};
  const ylm::Strahlkorper<Frame::Grid> strahlkorper_origin{
      2_st, radius, std::array<double, 3>{{0.0, 0.0, 0.0}}};

  // Build a positive sqrt(G_B / 48) profile with a nonzero dipole.
  const double l0_coef = 2.5;
  const double l1z_coef = 0.2;

  const auto& ylm = strahlkorper_origin.ylm_spherepack();
  DataVector gb_spectral{ylm.spectral_size(), 0.0};
  ylm::SpherepackIterator iterator{2, 2};
  gb_spectral[iterator.set(0, 0)()] = l0_coef / sqrt(M_PI / 2.0);
  gb_spectral[iterator.set(1, 0)()] = l1z_coef / sqrt(M_PI / 2.0);

  const DataVector sqrt_gb_over_48 = ylm.spec_to_phys(gb_spectral);
  Scalar<DataVector> gauss_bonnet_scalar{ylm.physical_size()};
  get(gauss_bonnet_scalar) = 48.0 * sqrt_gb_over_48 * sqrt_gb_over_48;

  const DataVector computed_center_shifted =
      control_system::measurements::gauss_bonnet_center(strahlkorper_shifted,
                                                        gauss_bonnet_scalar);
  const DataVector computed_center_origin =
      control_system::measurements::gauss_bonnet_center(strahlkorper_origin,
                                                        gauss_bonnet_scalar);

  CHECK(std::abs(computed_center_origin[2]) > 1.0e-6);
  CHECK(computed_center_shifted[0] ==
        approx(computed_center_origin[0] + sphere_center[0]).epsilon(1.0e-10));
  CHECK(computed_center_shifted[1] ==
        approx(computed_center_origin[1] + sphere_center[1]).epsilon(1.0e-10));
  CHECK(computed_center_shifted[2] ==
        approx(computed_center_origin[2] + sphere_center[2]).epsilon(1.0e-10));
}

}  // namespace
