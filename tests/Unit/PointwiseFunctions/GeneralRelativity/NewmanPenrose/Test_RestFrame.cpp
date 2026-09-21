// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NpTestHelpers.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
namespace helpers = TestHelpers::gr::np;
using frame = Frame::Inertial;
const std::string module = "NpMatching";

void test_null_pair_and_member() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 7;
  const Scalar<ComplexDataVector> a_bar{
      helpers::random_complex(make_not_null(&generator), num_points, 0.2)};
  const Scalar<ComplexDataVector> b{
      helpers::random_complex(make_not_null(&generator), num_points, 0.2)};
  AdaptedFourVector ell{};
  AdaptedFourVector kay{};
  principal_null_pair(make_not_null(&ell), make_not_null(&kay), a_bar, b);
  CHECK_ITERABLE_APPROX(ell, pypp::call<AdaptedFourVector>(
                                 module, "principal_null_outgoing", a_bar, b));
  CHECK_ITERABLE_APPROX(kay, pypp::call<AdaptedFourVector>(
                                 module, "principal_null_incoming", a_bar, b));
  // Null and normalized: g(l, k) = -1 with Minkowski diag(-1, 1, 1, 1)
  DataVector ell_norm = -square(get<0>(ell));
  DataVector kay_norm = -square(get<0>(kay));
  DataVector product = -get<0>(ell) * get<0>(kay);
  for (size_t i = 1; i < 4; ++i) {
    ell_norm += square(ell.get(i));
    kay_norm += square(kay.get(i));
    product += ell.get(i) * kay.get(i);
  }
  Approx custom_approx = Approx::custom().epsilon(1.e-12).scale(1.);
  CHECK_ITERABLE_CUSTOM_APPROX(ell_norm, DataVector(num_points, 0.),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(kay_norm, DataVector(num_points, 0.),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(product, DataVector(num_points, -1.),
                               custom_approx);

  const TangentBoostMember member = tangent_boost_member(a_bar, b);
  CHECK_ITERABLE_APPROX(
      member.radial_direction,
      pypp::call<TriadVector>(module, "tangent_radial_direction", a_bar, b));
  CHECK_ITERABLE_APPROX(
      member.transverse_velocity,
      pypp::call<TriadVector>(module, "tangent_transverse_velocity", a_bar, b));
  CHECK_ITERABLE_APPROX(member.lorentz_factor,
                        pypp::call<Scalar<DataVector>>(
                            module, "tangent_lorentz_factor", a_bar, b));
  // Unit radial leg orthogonal to the transverse velocity, Gamma =
  // 1/sqrt(1-w^2)
  DataVector radial_norm(num_points, 0.);
  DataVector overlap(num_points, 0.);
  DataVector speed_squared(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    radial_norm += square(member.radial_direction.get(i));
    overlap +=
        member.radial_direction.get(i) * member.transverse_velocity.get(i);
    speed_squared += square(member.transverse_velocity.get(i));
  }
  CHECK_ITERABLE_CUSTOM_APPROX(radial_norm, DataVector(num_points, 1.),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(overlap, DataVector(num_points, 0.),
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get(member.lorentz_factor),
                               DataVector{1. / sqrt(1. - speed_squared)},
                               custom_approx);
}

void test_gradient_boost() {
  MAKE_GENERATOR(generator);
  const size_t num_points = 40;
  const auto geometry =
      helpers::random_geometry(make_not_null(&generator), num_points);
  // A tangent member from a mildly misaligned type-D frame
  const auto data =
      helpers::near_type_d(make_not_null(&generator), num_points, 0.1, 1.e-3);
  const TypeDRotation rotation = solve_type_d_rotation(
      data.psi, coulomb_scalar(invariant_i(data.psi), invariant_j(data.psi)));
  const TangentBoostMember member =
      tangent_boost_member(rotation.a_bar, rotation.b);

  // Pointwise, with a random gradient that is dominantly radial so the
  // rapidity stays physical
  std::uniform_real_distribution<> small(-0.1, 0.1);
  const DataVector used_for_size(num_points);
  auto spatial_gradient =
      make_with_random_values<tnsr::i<DataVector, 3, frame>>(
          make_not_null(&generator), make_not_null(&small), used_for_size);
  for (size_t i = 0; i < 3; ++i) {
    spatial_gradient.get(i) -= 3. * geometry.directions.get(i);
  }
  const auto time_derivative = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), make_not_null(&small), used_for_size);
  CHECK_ITERABLE_APPROX(
      invariant_rapidity(member, geometry.rotation, geometry.spatial_metric,
                         geometry.lapse, geometry.shift, spatial_gradient,
                         time_derivative),
      pypp::call<Scalar<DataVector>>(
          module, "invariant_rapidity", member.radial_direction,
          member.transverse_velocity, member.lorentz_factor, geometry.rotation,
          geometry.spatial_metric, geometry.lapse, geometry.shift,
          spatial_gradient, time_derivative));
  CHECK_ITERABLE_APPROX(
      DataVector(tanh(get(invariant_rapidity(
          member, geometry.rotation, geometry.spatial_metric, geometry.lapse,
          geometry.shift, spatial_gradient, time_derivative)))),
      get(invariant_tanh_rapidity(
          member, geometry.rotation, geometry.spatial_metric, geometry.lapse,
          geometry.shift, spatial_gradient, time_derivative)));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.Gr.NewmanPenrose.RestFrame",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/GeneralRelativity/NewmanPenrose/"};
  test_null_pair_and_member();
  test_gradient_boost();
}
}  // namespace gr::np
