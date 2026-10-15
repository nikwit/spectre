// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeExcisionRecordingData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeRecordingData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WorldtubeTestHelpers.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Strahlkorper/Strahlkorper.hpp"
#include "NumericalAlgorithms/Strahlkorper/StrahlkorperFunctions.hpp"
#include "Utilities/ProtocolHelpers.hpp"

namespace {
namespace wt = gh::worldtube;
const std::array<double, 3> center{{4., 0.2, -0.1}};
constexpr double radius = 0.12;
constexpr double distortion = 0.08;
constexpr double shear = 0.1;

// Nonuniform, non-diagonal logical-to-inertial map with time-dependent scale.
// Its full inverse Jacobian is needed even though the recorded tensor
// components are already expressed in inertial coordinates.
struct Geometry {
  explicit Geometry(const size_t extent, const double time)
      : mesh(extent, Spectral::Basis::Legendre,
             Spectral::Quadrature::GaussLobatto),
        coords(mesh.number_of_grid_points()),
        inverse_jacobian(mesh.number_of_grid_points(), 0.),
        scale{{0.5 + 0.1 * time, 0.7 + 0.05 * time, 0.9 + 0.2 * time}} {
    const auto logical = logical_coordinates(mesh);
    get<0>(coords) =
        center[0] +
        scale[0] * (get<0>(logical) + distortion * square(get<0>(logical))) +
        shear * get<1>(logical);
    get<1>(coords) = center[1] + scale[1] * get<1>(logical);
    get<2>(coords) = center[2] + scale[2] * get<2>(logical);
    get<0, 0>(inverse_jacobian) =
        1. / (scale[0] * (1. + 2. * distortion * get<0>(logical)));
    get<0, 1>(inverse_jacobian) =
        -shear / scale[1] * get<0, 0>(inverse_jacobian);
    get<1, 1>(inverse_jacobian) = 1. / scale[1];
    get<2, 2>(inverse_jacobian) = 1. / scale[2];
  }

  tnsr::I<DataVector, 3, Frame::ElementLogical> inverse(
      const tnsr::I<DataVector, 3>& x) const {
    tnsr::I<DataVector, 3, Frame::ElementLogical> result(get<0>(x).size());
    get<1>(result) = (get<1>(x) - center[1]) / scale[1];
    get<2>(result) = (get<2>(x) - center[2]) / scale[2];
    get<0>(result) =
        (sqrt(1. + 4. * distortion *
                       (get<0>(x) - center[0] - shear * get<1>(result)) /
                       scale[0]) -
         1.) /
        (2. * distortion);
    return result;
  }

  Mesh<3> mesh;
  tnsr::I<DataVector, 3> coords;
  InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian;
  std::array<double, 3> scale;
};

void test_volume_derivatives(const double time) {
  CAPTURE(time);
  const Geometry geometry(8, time);
  const size_t n = geometry.mesh.number_of_grid_points();
  Variables<wt::recording_source_tags> source(n, 0.);
  get<wt::recording_inverse_jacobian>(source) = geometry.inverse_jacobian;
  // This profile vanishes on the recording sphere, but has a nonzero normal
  // derivative there. Differentiating the interpolated surface values loses
  // that information and cannot pass this test.
  DataVector profile(n, -square(radius));
  for (size_t i = 0; i < 3; ++i) {
    profile += square(geometry.coords.get(i) - center[i]);
  }
  tmpl::for_each<wt::recording_evolved_tags>([&](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    auto& field = get<Tag>(source);
    for (size_t c = 0; c < field.size(); ++c) {
      field[c] = 0.001 * static_cast<double>(c + 1) * profile;
    }
  });
  auto& metric = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(source);
  get<0, 0>(metric) -= 1.;
  for (size_t i = 1; i < 4; ++i) {
    metric.get(i, i) += 1.;
  }
  Variables<wt::recording_tags> volume;
  wt::ComputeRecordingData::apply(make_not_null(&volume), source,
                                  geometry.mesh);
  Variables<wt::excision_recording_tags<Frame::Inertial>> stationary;
  wt::ComputeExcisionRecordingData::apply(make_not_null(&stationary), source,
                                          geometry.mesh);
  // A moving map x_inertial = 2 x_grid + translation(t). The GH tensor
  // components and derivatives must stay inertial, while the shift includes
  // the mesh velocity and transforms to grid components.
  Jacobian<DataVector, 3, Frame::Grid, Frame::Inertial> jac(n, 0.);
  InverseJacobian<DataVector, 3, Frame::Grid, Frame::Inertial> invjac(n, 0.);
  Jacobian<DataVector, 3, Frame::ElementLogical, Frame::Grid> logical_jac(n);
  InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Grid>
      logical_invjac(n);
  InverseJacobian<DataVector, 3, Frame::Grid, Frame::Grid> identity(n, 0.);
  const auto inertial_jac =
      determinant_and_inverse(geometry.inverse_jacobian).second;
  const tnsr::I<DataVector, 3> velocity(n, 0.03);
  const tnsr::I<DataVector, 3, Frame::Grid> grid_velocity(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    jac.get(i, i) = 2.;
    invjac.get(i, i) = 0.5;
    identity.get(i, i) = 1.;
    for (size_t j = 0; j < 3; ++j) {
      logical_jac.get(i, j) = 0.5 * inertial_jac.get(i, j);
      logical_invjac.get(i, j) = 2. * geometry.inverse_jacobian.get(i, j);
    }
  }
  Variables<wt::excision_recording_tags<Frame::Grid>> moving;
  wt::ComputeExcisionRecordingData::apply(
      make_not_null(&moving), source, geometry.mesh, jac, invjac, logical_jac,
      logical_invjac, identity, velocity, grid_velocity);
  tmpl::for_each<wt::recording_tags>([&](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    CHECK_ITERABLE_APPROX(get<Tag>(moving), get<Tag>(volume));
    CHECK_ITERABLE_APPROX(get<Tag>(stationary), get<Tag>(volume));
  });
  CHECK_ITERABLE_APPROX(get<gr::Tags::Lapse<DataVector>>(moving),
                        get<gr::Tags::Lapse<DataVector>>(stationary));
  for (size_t i = 0; i < 3; ++i) {
    CHECK_ITERABLE_APPROX(
        (get<gr::Tags::Shift<DataVector, 3, Frame::Grid>>(moving).get(i)),
        DataVector(0.5 *
                   (get<gr::Tags::Shift<DataVector, 3>>(stationary).get(i) +
                    velocity.get(i))));
  }
  const auto points = ylm::cartesian_coords(
      ylm::Strahlkorper<Frame::Inertial>(5, radius, center));
  const intrp::Irregular<3> interpolator(geometry.mesh,
                                         geometry.inverse(points));
  const auto surface = interpolator.interpolate(moving);
  const auto tolerance = Approx::custom().epsilon(2.e-11).scale(1.);
  tmpl::for_each<wt::recording_evolved_tags>([&](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    const auto& field = get<Tag>(surface);
    const auto& deriv = get<wt::recording_deriv<Tag>>(surface);
    for (size_t c = 0; c < field.size(); ++c) {
      const auto indices = field.get_tensor_index(c);
      double expected_value = 0.;
      if constexpr (std::is_same_v<Tag,
                                   gr::Tags::SpacetimeMetric<DataVector, 3>>) {
        if (indices[0] == indices[1]) {
          expected_value = indices[0] == 0 ? -1. : 1.;
        }
      }
      CHECK_ITERABLE_CUSTOM_APPROX(
          field[c], DataVector(get<0>(points).size(), expected_value),
          tolerance);
      std::array<size_t, Tag::type::rank() + 1> d_indices{};
      for (size_t j = 0; j < indices.size(); ++j) {
        d_indices[j + 1] = indices[j];
      }
      for (size_t i = 0; i < 3; ++i) {
        d_indices[0] = i;
        CHECK_ITERABLE_CUSTOM_APPROX(
            deriv.get(d_indices),
            DataVector(0.002 * static_cast<double>(c + 1) *
                       (points.get(i) - center[i])),
            tolerance);
      }
    }
  });
}

void test_curvature() {
  const Geometry geometry(14, 0.3);
  const gr::Solutions::KerrSchild solution(1., {{0., 0., 0.}}, {{0., 0., 0.}});
  const auto evolved = TestHelpers::gh_worldtube::evolved_variables(
      solution, geometry.coords, 0.);
  Variables<wt::recording_source_tags> source(
      geometry.mesh.number_of_grid_points());
  tmpl::for_each<wt::recording_evolved_tags>([&](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    get<Tag>(source) = get<Tag>(evolved);
  });
  get<wt::recording_inverse_jacobian>(source) = geometry.inverse_jacobian;
  Variables<wt::recording_tags> volume;
  wt::ComputeRecordingData::apply(make_not_null(&volume), source,
                                  geometry.mesh);
  const auto points = ylm::cartesian_coords(
      ylm::Strahlkorper<Frame::Inertial>(5, radius, center));
  const intrp::Irregular<3> interpolator(geometry.mesh,
                                         geometry.inverse(points));
  const auto surface = interpolator.interpolate(volume);
  // Schwarzschild gives K=48/r^6 and d_i K=-288 x_i/r^8. Check the
  // differentiated volume data and its interpolation against this independent
  // analytic answer, including the spatially varying inverse Jacobian.
  const auto check = [](const auto& data, const auto& x) {
    const auto exact =
        TestHelpers::gh_worldtube::boosted_kretschmann(x, {{0., 0., 0.}}, 0.);
    const auto tolerance = Approx::custom().epsilon(2.e-9).scale(1.);
    CHECK_ITERABLE_CUSTOM_APPROX(get<wt::Tags::Kretschmann>(data),
                                 exact.kretschmann, tolerance);
    CHECK_ITERABLE_CUSTOM_APPROX(
        get<wt::recording_deriv<wt::Tags::Kretschmann>>(data),
        exact.d_kretschmann, tolerance);
  };
  check(volume, geometry.coords);
  check(surface, points);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Gh.Worldtube.ComputeRecordingData",
                  "[Unit][Evolution]") {
  static_assert(
      tt::assert_conforms_to_v<wt::ComputeRecordingData,
                               intrp::protocols::ComputeVarsToInterpolate>);
  test_volume_derivatives(0.);
  test_volume_derivatives(0.7);
  test_curvature();
}
