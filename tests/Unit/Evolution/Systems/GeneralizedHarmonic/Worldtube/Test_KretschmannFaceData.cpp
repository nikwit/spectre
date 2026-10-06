// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
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
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Actions/UpdateKretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "PointwiseFunctions/ConstraintDamping/Constant.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"

namespace {
namespace helpers = TestHelpers::gh_worldtube;
using gh::worldtube::KretschmannFaceData;

void test_geometric_face_geometry() {
  const auto setup = helpers::shell_element(10, 8, {{0., 0., 0.}}, 2.5, 3.);
  const double angle = .31, scale = 1.7;
  auto coordinates = setup.inertial_coords;
  auto inverse_jacobian = setup.inverse_jacobian;
  // x' = scale R x + center; J'^{-1} = J^{-1} R^T / scale.
  const std::array<double, 3> center{{-12., 3., .7}};
  const std::array<std::array<double, 3>, 3> rotation{
      {{{cos(angle), 0., sin(angle)}},
       {{0., 1., 0.}},
       {{-sin(angle), 0., cos(angle)}}}};
  for (size_t i = 0; i < 3; ++i) {
    coordinates.get(i) = center[i];
    for (size_t j = 0; j < 3; ++j) {
      coordinates.get(i) +=
          scale * rotation[i][j] * setup.inertial_coords.get(j);
      inverse_jacobian.get(j, i) = 0.;
      for (size_t k = 0; k < 3; ++k) {
        inverse_jacobian.get(j, i) +=
            setup.inverse_jacobian.get(j, k) * rotation[i][k] / scale;
      }
    }
  }
  const auto validate = [&]() {
    gh::worldtube::validate_geometric_matching_face(
        coordinates, setup.mesh, inverse_jacobian, Direction<3>::lower_xi());
  };
  CHECK_NOTHROW(validate());
  const auto good_coordinates = coordinates;
  coordinates.get(0) = center[0] + 1.01 * (coordinates.get(0) - center[0]);
  CHECK_THROWS_WITH(validate(), Catch::Matchers::ContainsSubstring(
                                    "round inertial-coordinate sphere"));
  coordinates = good_coordinates;
  inverse_jacobian.get(0, 0) += .01;
  CHECK_THROWS_WITH(validate(), Catch::Matchers::ContainsSubstring(
                                    "radial coordinate-normal"));
  const auto wedge = helpers::wedge_element(5, {{0., 0., 0.}});
  CHECK_THROWS_WITH(gh::worldtube::validate_geometric_matching_face(
                        wedge.inertial_coords, wedge.mesh,
                        wedge.inverse_jacobian, Direction<3>::lower_zeta()),
                    Catch::Matchers::ContainsSubstring(
                        "complete spherical-harmonic inner face"));
}

void test_geometric_relaxation_face() {
  const auto setup = helpers::shell_element(16, 8, {{0., 0., 0.}}, 2.5, 3.);
  const auto vars = setup.evolved_variables(0.);
  const auto& metric = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars);
  const auto& pi = get<gh::Tags::Pi<DataVector, 3>>(vars);
  const auto& phi = get<gh::Tags::Phi<DataVector, 3>>(vars);
  KretschmannFaceData<3> data{};
  const auto update = [&](const TimeStepId& id) {
    const double time = id.substep_time();
    gh::worldtube::update_kretschmann_face_data<3>(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.element, setup.domain.excision_spheres(),
        std::nullopt, time);
    gh::worldtube::update_geometric_quadrupole_relaxation(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.inertial_coords, std::nullopt, 1., .5,
        true, time, id);
  };
  const Slab slab{0., .1};
  const TimeStepId start{true, 0, slab.start()};
  update(start);
  REQUIRE(data.quadrupole_relaxation.has_value());
  REQUIRE(data.filtered_moments.has_value());
  CHECK_FALSE(data.third_order.has_value());
  CHECK(data.quadrupole_relaxation->samples.size() == 1);
  CHECK(data.quadrupole_relaxation->samples.back().clock_rate ==
        Approx::custom().epsilon(1.e-8).scale(1.)(1.));
  // Inject a previous physical tide. It should decay into the zero-tide
  // Schwarzschild fit, and the imposed target must consume this history.
  const auto initial_raw = *data.filtered_moments;
  auto& initial = data.quadrupole_relaxation->samples.back();
  initial.filtered[0] += .001;
  initial.filtered[1] -= .001;
  const auto anchor = *data.quadrupole_relaxation;
  const auto stage = start.next_substep(slab.duration(), .5);
  update(stage);
  CHECK(*data.quadrupole_relaxation == anchor);
  const auto staged = *data.filtered_moments;
  update(stage);
  CHECK_ITERABLE_APPROX(*data.filtered_moments, staged);
  CHECK(*data.quadrupole_relaxation == anchor);
  update(start.next_step(slab.duration()));
  CHECK(data.quadrupole_relaxation->samples.size() == 2);
  CHECK(data.filtered_moments_time == .1);
  const auto expected = initial_raw[0] + .001 * exp(-.1 / .5);
  CHECK(abs((*data.filtered_moments)[0] - expected) < 1.e-10);
  CHECK(serialize_and_deserialize(data) == data);
  const auto face = gh::worldtube::face_curvature(
      metric, pi, phi, setup.mesh, setup.inverse_jacobian, *data.direction);
  const auto target = gh::worldtube::evaluate_matching(
      gh::worldtube::PhysicalModel::QuadrupoleGeometric, 1., face.electric,
      face.magnetic, face.spatial_metric, face.unit_normal_covector, face.lapse,
      face.shift, &data, data.filtered_moments);
  CHECK(target.second_order->fit.components == *data.filtered_moments);
  CHECK(max(abs(get(target.psi0_target))) > 1.e-5);
}

void test_fixed_point_face() {
  const auto setup = helpers::shell_element(16, 8, {{0., 0., 0.}}, 2.5, 3.);
  const auto vars = setup.evolved_variables(0.);
  const auto& metric = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars);
  const auto& pi = get<gh::Tags::Pi<DataVector, 3>>(vars);
  const auto& phi = get<gh::Tags::Phi<DataVector, 3>>(vars);
  KretschmannFaceData<3> data{};
  gh::worldtube::update_kretschmann_face_data<3>(
      make_not_null(&data), metric, pi, phi, setup.mesh, setup.inverse_jacobian,
      setup.element, setup.domain.excision_spheres(), std::nullopt, 0.);
  const auto original_k = data.kretschmann;
  const auto original_dt_k = data.dt_kretschmann;
  data.filtered_moments = gr::np::TidalMoments{};
  data.quadrupole_relaxation.emplace();
  const auto update = [&]() {
    gh::worldtube::update_geometric_fixed_point(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.inertial_coords, 1., data.time, {});
  };
  update();
  REQUIRE(data.psi0_fixed_point.has_value());
  CHECK_FALSE(data.filtered_moments.has_value());
  CHECK_FALSE(data.quadrupole_relaxation.has_value());
  CHECK(data.kretschmann == original_k);
  CHECK(data.dt_kretschmann == original_dt_k);
  CHECK(max(abs(get(data.psi0_fixed_point->psi0_target))) < 1.e-10);
  CHECK(serialize_and_deserialize(data) == data);
  const auto first = *data.psi0_fixed_point;
  update();
  CHECK(*data.psi0_fixed_point == first);
  // Exercise the actual RHS action: it must bypass the old measured-Psi0
  // moment fit, and clear this cache again when another model is selected.
  using BC = gh::BoundaryConditions::WorldtubeTypeD<3>;
  using Sector = gh::BoundaryConditions::detail::SectorImposition;
  std::vector<DirectionMap<
      3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      conditions(1);
  conditions[0][*data.direction] = std::make_unique<BC>(
      Sector::Bjorhus, Sector::Bjorhus, Sector::Frozen,
      gh::BoundaryConditions::detail::GeometricFixedPoint{}, 1.);
  const Slab slab{0., .1};
  const TimeStepId id{true, 0, slab.start()};
  const auto action = [&]() {
    gh::worldtube::UpdateKretschmannFaceData<3>::apply(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.element, setup.domain, std::nullopt, 0.,
        conditions, setup.inertial_coords, id,
        ConstraintDamping::Constant<3, Frame::Grid>{.5}, {});
  };
  data.filtered_moments = gr::np::TidalMoments{};
  action();
  CHECK_FALSE(data.filtered_moments.has_value());
  CHECK(*data.psi0_fixed_point == first);
  conditions[0][*data.direction] = std::make_unique<BC>(
      Sector::Bjorhus, Sector::Bjorhus, Sector::Frozen,
      gh::worldtube::PhysicalModel::QuadrupoleGeometric, 1.);
  action();
  CHECK_FALSE(data.psi0_fixed_point.has_value());
  CHECK(data.filtered_moments.has_value());
  update();
  // Recompute at a stage and on rollback: no iteration history is reused.
  data.time = .05;
  update();
  CHECK(data.psi0_fixed_point->psi0_target == first.psi0_target);
  data.time = 0.;
  update();
  CHECK(*data.psi0_fixed_point == first);
  const auto face = gh::worldtube::face_curvature(
      metric, pi, phi, setup.mesh, setup.inverse_jacobian, *data.direction);
  const auto matching = [&]() {
    return gh::worldtube::evaluate_matching(
        gh::worldtube::PhysicalModel::QuadrupoleGeometric, 1., face.electric,
        face.magnetic, face.spatial_metric, face.unit_normal_covector,
        face.lapse, face.shift, &data, data.filtered_moments, true);
  };
  CHECK(matching().psi0_target == first.psi0_target);
  data.filtered_moments = gr::np::TidalMoments{};
  CHECK_THROWS_WITH(matching(), Catch::Matchers::ContainsSubstring(
                                    "without imposed moments"));
  data.filtered_moments.reset();
  data.time = .1;
  CHECK_THROWS_WITH(matching(), Catch::Matchers::ContainsSubstring(
                                    "current converged target"));
  data.psi0_fixed_point.reset();
  CHECK_THROWS_WITH(matching(), Catch::Matchers::ContainsSubstring(
                                    "current converged target"));
}

void test_third_order_face() {
  const auto setup = helpers::shell_element(16, 8, {{0., 0., 0.}}, 2.5, 3.);
  const auto vars = setup.evolved_variables(0.);
  const auto& metric = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars);
  const auto& pi = get<gh::Tags::Pi<DataVector, 3>>(vars);
  const auto& phi = get<gh::Tags::Phi<DataVector, 3>>(vars);
  KretschmannFaceData<3> data{};
  const auto update = [&](const int step) {
    const double time = .1 * step;
    gh::worldtube::update_kretschmann_face_data<3>(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.element, setup.domain.excision_spheres(),
        std::nullopt, time);
    const TimeStepId id{true, step, Slab{time, time + .1}.start()};
    gh::worldtube::update_third_order_matching(
        make_not_null(&data), metric, pi, phi, setup.mesh,
        setup.inverse_jacobian, setup.inertial_coords, std::nullopt, 1., time,
        id);
  };
  for (int step = 0; step < 5; ++step)
    update(step);
  REQUIRE(data.third_order.has_value());
  const auto& state = *data.third_order;
  CHECK(state.derivative_order == 4);
  CHECK(state.clock_rate == Approx::custom().epsilon(1.e-8).scale(1.)(1.));
  CHECK(max(abs(get(state.psi0_target))) < 1.e-9);
  for (const double dot : state.dots)
    CHECK(abs(dot) < 1.e-10);
  CHECK(state.fit_condition < 20.);
  CHECK(serialize_and_deserialize(data) == data);
  const auto saved = data;
  update(4);
  CHECK(data == saved);
  // Matching consumes the current cached target and retains wminus/2.
  const auto face = gh::worldtube::face_curvature(
      metric, pi, phi, setup.mesh, setup.inverse_jacobian, *data.direction);
  const auto target = gh::worldtube::evaluate_matching(
      gh::worldtube::PhysicalModel::ThirdOrderGeometric, 1., face.electric,
      face.magnetic, face.spatial_metric, face.unit_normal_covector, face.lapse,
      face.shift, &data);
  CHECK_ITERABLE_APPROX(target.psi0_target, data.third_order->psi0_target);
  data.third_order->time -= .1;
  CHECK_THROWS_WITH(
      gh::worldtube::evaluate_matching(
          gh::worldtube::PhysicalModel::ThirdOrderGeometric, 1., face.electric,
          face.magnetic, face.spatial_metric, face.unit_normal_covector,
          face.lapse, face.shift, &data),
      Catch::Matchers::ContainsSubstring("needs its current target"));
  update(0);
  CHECK(data.third_order->derivative_order == 0);
  CHECK(data.third_order->history.samples.size() == 1);
}

void test_radial_gauge_capture() {
  const auto setup = helpers::shell_element(18, 5, {{0., 0., 0.}}, 2.5, 3.);
  const size_t n = setup.mesh.number_of_grid_points();
  tnsr::aa<DataVector, 3> metric(n, 0.), pi(n, 0.);
  tnsr::iaa<DataVector, 3> phi(n, 0.);
  get<0, 0>(metric) = -1.;
  for (size_t i = 1; i < 4; ++i) {
    metric.get(i, i) = 1.;
  }
  DataVector radius(n, 0.);
  for (const auto& x : setup.inertial_coords) {
    radius += square(x);
  }
  radius = sqrt(radius);
  get<0, 0>(pi) = .03 * square(radius);
  const Scalar<DataVector> gamma{DataVector{.2 + .01 * radius}};
  std::optional<gh::worldtube::RadialGaugeFaceData> rd{};
  const auto update = [&](double t) {
    gh::worldtube::update_radial_gauge_face_data(
        make_not_null(&rd), metric, pi, phi, gamma, setup.inertial_coords,
        tnsr::I<double, 3>(0.), setup.mesh, setup.inverse_jacobian,
        Direction<3>::lower_xi(), t);
  };
  update(0.);
  REQUIRE(rd.has_value());
  const auto initial_q = rd->q, initial_dr = rd->dr_q;
  const size_t nf = setup.mesh.extents(1) * setup.mesh.extents(2);
  CHECK_ITERABLE_APPROX(get<0>(rd->q),
                        DataVector(nf, (.03 * 2.5 * 2.5 + .225) / sqrt(2.)));
  CHECK_ITERABLE_CUSTOM_APPROX(get<0>(rd->dr_q),
                               DataVector(nf, (.06 * 2.5 + .01) / sqrt(2.)),
                               Approx::custom().epsilon(1.e-10).scale(1.));
  for (size_t i = 0; i < 3; ++i) {
    CHECK_ITERABLE_APPROX(
        rd->q.get(i + 1),
        DataVector{.225 / sqrt(2.) * rd->radial_direction.get(i)});
    CHECK_ITERABLE_CUSTOM_APPROX(
        rd->dr_q.get(i + 1),
        DataVector{.01 / sqrt(2.) * rd->radial_direction.get(i)},
        Approx::custom().epsilon(1.e-10).scale(1.));
  }
  get<0, 0>(pi) *= 1.7;
  update(1.);
  update(0.);  // An AB self-start reset must not re-capture.
  CHECK_ITERABLE_APPROX(rd->initial_q, initial_q);
  CHECK_ITERABLE_APPROX(rd->initial_dr_q, initial_dr);
  CHECK(max(abs(get<0>(rd->q) - get<0>(initial_q))) > .01);
  KretschmannFaceData<3> holder{};
  holder.radial_gauge = rd;
  CHECK(serialize_and_deserialize(holder) == holder);
}

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
  REQUIRE(data.initial_gauge_difference.has_value());
  const auto initial_gauge = *data.initial_gauge_difference;
  CHECK(get<0>(initial_gauge).size() == num_face_points);
  CHECK(max(abs(get<0>(initial_gauge))) > 0.);

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

  CHECK_ITERABLE_APPROX(*data.initial_gauge_difference, initial_gauge);
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

void test_relax_tidal_moments() {
  using gh::worldtube::relax_tidal_moments;
  const gr::np::TidalMoments raw{
      {{1., 0.}, {0., 2.}, {3., 0.}, {0., 4.}, {5., 5.}}};
  const gr::np::TidalMoments other{
      {{2., 0.}, {0., 0.}, {1., 1.}, {0., 0.}, {0., 0.}}};
  std::optional<gr::np::TidalMoments> moments{};
  double moments_time = std::numeric_limits<double>::signaling_NaN();
  // First call initializes with the raw fit, whatever the relaxation time
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      raw, 1., 10.);
  REQUIRE(moments.has_value());
  CHECK(*moments == raw);
  CHECK(moments_time == 1.);
  // Forward Euler over the elapsed time
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      other, 3., 10.);
  for (size_t a = 0; a < 5; ++a) {
    CHECK(gsl::at(*moments, a) ==
          gsl::at(raw, a) + 0.2 * (gsl::at(other, a) - gsl::at(raw, a)));
  }
  CHECK(moments_time == 3.);
  // Same time: unchanged
  const auto saved = *moments;
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      raw, 3., 10.);
  CHECK(*moments == saved);
  // A step longer than the relaxation time is capped at one relaxation time
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      raw, 100., 10.);
  CHECK(*moments == raw);
  // Time running backwards resets
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      other, 50., 10.);
  CHECK(*moments == other);
  CHECK(moments_time == 50.);
  // No relaxation time: always the raw fit
  relax_tidal_moments(make_not_null(&moments), make_not_null(&moments_time),
                      raw, 51., std::nullopt);
  CHECK(*moments == raw);
}

void test_update_filtered_tidal_moments() {
  const std::array<double, 3> velocity{{0.2, -0.1, 0.15}};
  const auto setup = helpers::wedge_element(12, velocity);
  const auto vars = setup.evolved_variables(0.);
  KretschmannFaceData<3> data{};
  const auto update = [&setup, &vars, &data](const double time) {
    gh::worldtube::update_kretschmann_face_data(
        make_not_null(&data),
        get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
        get<gh::Tags::Pi<DataVector, 3>>(vars),
        get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
        setup.inverse_jacobian, setup.element, setup.domain.excision_spheres(),
        std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>{}, time);
  };
  // Static history so that dt K is the comoving estimate, as in a real run
  update(-1.e-3);
  update(0.);
  CHECK_FALSE(data.filtered_moments.has_value());
  gh::worldtube::update_filtered_tidal_moments(
      make_not_null(&data), get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
      get<gh::Tags::Pi<DataVector, 3>>(vars),
      get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
      setup.inverse_jacobian, gh::worldtube::PhysicalModel::Quadrupole, 1.0,
      10.0, 0.);
  REQUIRE(data.filtered_moments.has_value());
  CHECK(data.filtered_moments_time == 0.);
  // The moments are the instantaneous fit on the face curvature
  const gh::worldtube::FaceCurvature face = gh::worldtube::face_curvature(
      get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
      get<gh::Tags::Pi<DataVector, 3>>(vars),
      get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
      setup.inverse_jacobian, *data.direction);
  const auto raw = gh::worldtube::evaluate_matching(
      gh::worldtube::PhysicalModel::Quadrupole, 1.0, face.electric,
      face.magnetic, face.spatial_metric, face.unit_normal_covector, face.lapse,
      face.shift, &data);
  CHECK(*data.filtered_moments == raw.second_order->fit.components);
  // The exact type-D hole carries no tide: the moments are truncation error
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::abs(gsl::at(*data.filtered_moments, a)) < 1.e-3);
  }
  // A repeated update at a later time with the same state relaxes toward the
  // same fit and leaves the moments unchanged
  const auto saved = *data.filtered_moments;
  gh::worldtube::update_filtered_tidal_moments(
      make_not_null(&data), get<gr::Tags::SpacetimeMetric<DataVector, 3>>(vars),
      get<gh::Tags::Pi<DataVector, 3>>(vars),
      get<gh::Tags::Phi<DataVector, 3>>(vars), setup.mesh,
      setup.inverse_jacobian, gh::worldtube::PhysicalModel::Quadrupole, 1.0,
      10.0, 1.);
  for (size_t a = 0; a < 5; ++a) {
    CHECK(std::abs(gsl::at(*data.filtered_moments, a) - gsl::at(saved, a)) <
          1.e-12);
  }
  CHECK(data.filtered_moments_time == 1.);
  CHECK(serialize_and_deserialize(data) == data);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.Worldtube.KretschmannFaceData",
    "[Unit][Evolution]") {
  test_geometric_face_geometry();
  test_geometric_relaxation_face();
  test_fixed_point_face();
  test_third_order_face();
  test_radial_gauge_capture();
  test_quadrature_weights();
  test_update_on_boosted_kerr_schild();
  test_relax_tidal_moments();
  test_update_filtered_tidal_moments();
}
