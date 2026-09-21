// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/SliceTensorToVariables.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace {
namespace helpers = TestHelpers::gh_worldtube;
using gh::worldtube::KretschmannFaceData;

void test_quadrature_weights() {
  // Spherical-harmonic face: Gauss-Legendre in cos(theta) times uniform phi,
  // theta varying fastest
  const size_t n_theta = 7;
  const size_t n_phi = 13;
  const Mesh<2> ylm_face{
      {{n_theta, n_phi}},
      {{Spectral::Basis::SphericalHarmonic,
        Spectral::Basis::SphericalHarmonic}},
      {{Spectral::Quadrature::Gauss, Spectral::Quadrature::Equiangular}}};
  const DataVector weights =
      gh::worldtube::face_quadrature_weights<3>(ylm_face);
  REQUIRE(weights.size() == n_theta * n_phi);
  CHECK(sum(weights) == approx(1.));
  const DataVector& gauss =
      Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                   Spectral::Quadrature::Gauss>(n_theta);
  for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
    for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
      // sum of the Gauss weights is 2, of the phi weights 2 pi
      CHECK(weights[i_theta + n_theta * i_phi] ==
            approx(gauss[i_theta] / (2. * static_cast<double>(n_phi))));
    }
  }

  // Legendre-Gauss-Lobatto face: product of the one-dimensional weights
  const Mesh<2> lgl_face{
      {{4, 5}}, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto};
  const DataVector lgl_weights =
      gh::worldtube::face_quadrature_weights<3>(lgl_face);
  const DataVector& lgl_4 =
      Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                   Spectral::Quadrature::GaussLobatto>(4);
  const DataVector& lgl_5 =
      Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                   Spectral::Quadrature::GaussLobatto>(5);
  for (size_t j = 0; j < 5; ++j) {
    for (size_t i = 0; i < 4; ++i) {
      CHECK(lgl_weights[i + 4 * j] == approx(lgl_4[i] * lgl_5[j] / 4.));
    }
  }

  // Lower dimensions
  const DataVector weights_1d =
      gh::worldtube::face_quadrature_weights<2>(Mesh<1>{
          5, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto});
  CHECK_ITERABLE_APPROX(weights_1d, DataVector(lgl_5 / 2.));
  CHECK(gh::worldtube::face_quadrature_weights<1>(Mesh<0>{}) == DataVector{1.});
}

void test_update_on_boosted_kerr_schild() {
  const std::array<double, 3> velocity{{0.2, -0.1, 0.15}};
  const auto setup = helpers::wedge_element(12, velocity);
  const auto& excision_spheres = setup.domain.excision_spheres();
  REQUIRE(excision_spheres.size() == 1);
  const auto evolved_at = [&setup](const double time) {
    return setup.evolved_variables(time);
  };
  const auto update =
      [&setup, &excision_spheres](
          const gsl::not_null<KretschmannFaceData<3>*> data,
          const helpers::EvolvedVariables& vars,
          const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
              mesh_velocity,
          const double time) {
        gh::worldtube::update_kretschmann_face_data(
            data, get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
            get<gh::Tags::Pi<DataVector, 3>>(vars),
            get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
            setup.inverse_jacobian, setup.element, excision_spheres,
            mesh_velocity, time);
      };

  KretschmannFaceData<3> data{};
  const auto vars_0 = evolved_at(0.);
  update(make_not_null(&data), vars_0, std::nullopt, 0.);
  REQUIRE(data.direction.has_value());
  CHECK(data.direction == excision_spheres.begin()->second.abutting_direction(
                              setup.element.id()));
  const Direction<3>& direction = *data.direction;
  const size_t sliced_dim = direction.dimension();
  const size_t fixed_index = index_to_slice_at(setup.mesh.extents(), direction);
  const size_t num_face_points =
      setup.mesh.slice_away(sliced_dim).number_of_grid_points();
  CHECK(get(data.kretschmann).size() == num_face_points);
  CHECK(get<0>(data.d_kretschmann).size() == num_face_points);
  CHECK(get(data.dt_kretschmann).size() == num_face_points);
  CHECK(data.quadrature_weights.size() == num_face_points);
  CHECK(sum(data.quadrature_weights) == approx(1.));
  CHECK(data.time == 0.);
  CHECK(data.previous_time == 0.);

  // The face is the inner boundary of the wedge, at coordinate radius 2.5
  const auto face_coords = get<domain::Tags::Coordinates<3, Frame::Inertial>>(
      data_on_slice<domain::Tags::Coordinates<3, Frame::Inertial>>(
          setup.mesh.extents(), sliced_dim, fixed_index,
          setup.inertial_coords));
  const DataVector face_radius =
      sqrt(square(get<0>(face_coords)) + square(get<1>(face_coords)) +
           square(get<2>(face_coords)));
  CHECK_ITERABLE_APPROX(face_radius, DataVector(num_face_points, 2.5));

  // K and its gradient against the boosted analytic values. K needs the
  // first numerical derivative (of Phi and Pi) and its gradient a second one,
  // so the gradient is the less accurate of the two at 12 points.
  const auto expected = helpers::boosted_kretschmann(face_coords, velocity, 0.);
  const double kretschmann_scale = max(abs(get(expected.kretschmann)));
  const double gradient_scale = max(abs(get<0>(expected.d_kretschmann)));
  CHECK_ITERABLE_CUSTOM_APPROX(
      get(data.kretschmann), get(expected.kretschmann),
      Approx::custom().epsilon(1.e-5).scale(kretschmann_scale));
  CHECK_ITERABLE_CUSTOM_APPROX(
      data.d_kretschmann, expected.d_kretschmann,
      Approx::custom().epsilon(1.e-3).scale(gradient_scale));
  // Without history and without a mesh velocity there is no time derivative
  CHECK(max(abs(get(data.dt_kretschmann))) == 0.);

  // Second evaluation: the backward difference is dt K to first order
  const double time_step = 1.e-3;
  const auto saved_kretschmann = data.kretschmann;
  update(make_not_null(&data), evolved_at(time_step), std::nullopt, time_step);
  CHECK(data.time == time_step);
  CHECK(data.previous_time == time_step);
  const auto expected_later =
      helpers::boosted_kretschmann(face_coords, velocity, time_step);
  CHECK_ITERABLE_CUSTOM_APPROX(
      get(data.dt_kretschmann), get(expected_later.dt_kretschmann),
      Approx::custom().epsilon(1.e-2).scale(
          max(abs(get(expected_later.dt_kretschmann)))));
  CHECK(get(data.previous_kretschmann) != get(saved_kretschmann));

  // Re-evaluation at the same time keeps the estimate and the history
  const KretschmannFaceData<3> saved = data;
  update(make_not_null(&data), evolved_at(time_step), std::nullopt, time_step);
  CHECK(data == saved);

  // A grid comoving with the hole: without history the estimate is the
  // comoving one, -v^i D_i K, which is exact here
  KretschmannFaceData<3> comoving{};
  tnsr::I<DataVector, 3, Frame::Inertial> mesh_velocity(
      setup.mesh.number_of_grid_points(), 0.);
  for (size_t i = 0; i < 3; ++i) {
    mesh_velocity.get(i) = gsl::at(velocity, i);
  }
  update(make_not_null(&comoving), vars_0, mesh_velocity, 0.);
  CHECK_ITERABLE_CUSTOM_APPROX(get(comoving.dt_kretschmann),
                               get(expected.dt_kretschmann),
                               Approx::custom().epsilon(1.e-3).scale(
                                   max(abs(get(expected.dt_kretschmann)))));

  // Time running backwards (a self-start reset) drops the history
  update(make_not_null(&data), vars_0, std::nullopt, -1.);
  CHECK(data.previous_time == -1.);
  CHECK(max(abs(get(data.dt_kretschmann))) == 0.);

  CHECK(serialize_and_deserialize(data) == data);
  CHECK(data != KretschmannFaceData<3>{});
  CHECK(KretschmannFaceData<3>{} == KretschmannFaceData<3>{});

  // The outer half of the wedge abuts no excision sphere
  const auto outer = helpers::wedge_element(
      6, velocity, 2.5, 4.0, 0,
      {{SegmentId{0, 0}, SegmentId{0, 0}, SegmentId{1, 1}}});
  KretschmannFaceData<3> outer_data = data;
  const auto outer_vars = outer.evolved_variables(0.);
  gh::worldtube::update_kretschmann_face_data(
      make_not_null(&outer_data),
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(outer_vars),
      get<gh::Tags::Pi<DataVector, 3>>(outer_vars),
      get<gh::Tags::Phi<DataVector, 3>>(outer_vars), outer.mesh,
      outer.inverse_jacobian, outer.element, outer.domain.excision_spheres(),
      std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>{}, 0.);
  CHECK_FALSE(outer_data.direction.has_value());
  CHECK(outer_data == KretschmannFaceData<3>{});
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.KretschmannFaceData",
    "[Unit][Evolution]") {
  test_quadrature_weights();
  test_update_on_boosted_kerr_schild();
}
