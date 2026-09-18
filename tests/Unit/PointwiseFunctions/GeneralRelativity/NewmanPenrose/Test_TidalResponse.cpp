// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <random>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
using frame = Frame::Inertial;
const std::string module = "NpMatching";

void test_tidal_response() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const double mass = 1. / 9.;
  std::uniform_real_distribution<> radius_dist(0.4, 1.6);
  std::uniform_real_distribution<> unit(-1., 1.);
  std::uniform_real_distribution<> small(-0.3, 0.3);
  const DataVector used_for_size(num_points);
  const auto radius = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), make_not_null(&radius_dist), used_for_size);

  const QuadrupoleProfiles profiles = quadrupole_profiles(radius, mass);
  for (const auto& [name, values] :
       {std::pair{"e_L", profiles.e_l}, std::pair{"e_V", profiles.e_v},
        std::pair{"e_T", profiles.e_t}, std::pair{"b_L", profiles.b_l},
        std::pair{"b_V", profiles.b_v}, std::pair{"b_T", profiles.b_t}}) {
    CHECK_ITERABLE_APPROX(
        Scalar<DataVector>{values},
        pypp::call<Scalar<DataVector>>(module, "quadrupole_profile", radius,
                                       mass, std::string{name}));
  }

  // Random unit direction and constant STF moments
  auto direction = make_with_random_values<TriadVector>(
      make_not_null(&generator), make_not_null(&unit), used_for_size);
  const DataVector norm =
      sqrt(square(get<0>(direction)) + square(get<1>(direction)) +
           square(get<2>(direction)));
  for (size_t i = 0; i < 3; ++i) {
    direction.get(i) /= norm;
  }
  std::array<double, 5> e_components{};
  std::array<double, 5> b_components{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(e_components, a) = unit(generator);
    gsl::at(b_components, a) = unit(generator);
  }
  const auto electric = stf_from_components(e_components);
  const auto magnetic = stf_from_components(b_components);
  CHECK(get<2, 2>(electric) == approx(-e_components[0] - e_components[1]));
  for (const bool transverse_only : {false, true}) {
    CAPTURE(transverse_only);
    CHECK_ITERABLE_APPROX(
        quadrupole_tide_tensor(electric, magnetic, direction, radius, mass,
                               transverse_only),
        pypp::call<ComplexMatrix>(module, "quadrupole_tide_tensor", electric,
                                  magnetic, direction, radius, mass,
                                  transverse_only));
  }

  const auto rapidity = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), make_not_null(&small), used_for_size);
  CHECK_ITERABLE_APPROX(self_dual_boost(direction, rapidity),
                        pypp::call<ComplexMatrix>(module, "self_dual_boost",
                                                  direction, rapidity));

  // Transverse velocity orthogonal to the direction
  auto velocity = make_with_random_values<TriadVector>(
      make_not_null(&generator), make_not_null(&small), used_for_size);
  DataVector overlap(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    overlap += velocity.get(i) * direction.get(i);
  }
  for (size_t i = 0; i < 3; ++i) {
    velocity.get(i) -= overlap * direction.get(i);
  }
  CHECK_ITERABLE_APPROX(
      rest_to_slice_map(direction, velocity, rapidity),
      pypp::call<ComplexMatrix>(module, "rest_to_slice_map", direction,
                                velocity, rapidity));
  // With zero transverse velocity the map is the radial boost alone
  const TriadVector zero_velocity(num_points, 0.);
  CHECK_ITERABLE_APPROX(rest_to_slice_map(direction, zero_velocity, rapidity),
                        self_dual_boost(direction, rapidity));

  const auto geometry =
      helpers::random_geometry(make_not_null(&generator), num_points);
  const auto columns = direct_tide_scalar_columns(
      direction, velocity, rapidity, radius, geometry.rotation, mass);
  for (size_t component = 0; component < 5; ++component) {
    CAPTURE(component);
    CHECK_ITERABLE_APPROX(
        gsl::at(columns, component),
        pypp::call<WeylScalars>(module, "direct_tide_scalar_column", direction,
                                velocity, rapidity, radius, geometry.rotation,
                                mass, component));
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.TidalResponse",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_tidal_response();
}
}  // namespace gr::np
