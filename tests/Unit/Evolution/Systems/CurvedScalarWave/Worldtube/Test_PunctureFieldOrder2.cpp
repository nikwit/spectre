// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldWorkspace.hpp"
#include "Framework/TestHelpers.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {
using Deriv = ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                            Frame::Inertial>;
using Fields =
    Variables<tmpl::list<CurvedScalarWave::Tags::Psi,
                         ::Tags::dt<CurvedScalarWave::Tags::Psi>, Deriv>>;

void test_reference() {
  // Values from the unspecialized five-output CSE graph, evaluated at
  // 50-digit precision. Input columns: xp, yp, vx, vy, ax, ay, dx, dy, dz, M.
  // Output columns: Psi, inertial dt Psi, dx Psi, dy Psi, dz Psi.
  const std::array<std::array<double, 15>, 16> reference{
      {{{8., 0., 0.040000000000000001, 0.089999999999999997, -0.002, 0.001,
         0.13, -0.17000000000000001, 0.089999999999999997, 0.,
         4.3031949314616016, -0.81304466553218879, -10.32645828767615,
         13.619482836109132, -7.1715919102794423}},
       {{0., 8., -0.040000000000000001, 0.092999999999999999, -0.002, 0.001,
         0.14000000000000001, -0.16600000000000001, 0.091999999999999998, 1.,
         3.8124542159295416, -1.5699348340554813, -7.9030326940163471,
         13.484545942690973, -5.0535625020654154}},
       {{-8., 0., 0.040000000000000001, 0.096000000000000002, -0.002, 0.001,
         0.14999999999999999, -0.16200000000000001, 0.094, 1.7,
         3.6455917821840957, -0.31450616508367618, -12.495443779233332,
         8.2851150177384802, -4.4893306791207319}},
       {{0., -8., -0.040000000000000001, 0.099000000000000005, -0.002, 0.001,
         0.16, -0.158, 0.096000000000000002, 0., 4.07327911044628,
         -1.5075530953043927, -10.873429670168578, 10.827167994980908,
         -6.4878939584617106}},
       {{6., 8., 0.040000000000000001, 0.10199999999999999, -0.002, 0.001,
         0.17000000000000001, -0.15400000000000003, 0.098000000000000004, 1.,
         3.9996205326165817, -0.63425785912463895, -10.537636170630449,
         10.34291683212593, -6.2659102184773223}},
       {{-6., 8., -0.040000000000000001, 0.105, -0.002, 0.001,
         0.17999999999999999, -0.15000000000000002, 0.099999999999999992, 1.7,
         3.1050957136240034, -1.2471622004103764, -8.3392153955729302,
         8.7640694822704415, -2.9468384747612855}},
       {{-6., -8., 0.040000000000000001, 0.108, -0.002, 0.001, 0.19,
         -0.14600000000000002, 0.10199999999999999, 0., 3.8379978800690639,
         -0.46822804948493074, -10.722851458938887, 8.304591119970496,
         -5.7665274446566919}},
       {{6., -8., -0.040000000000000001, 0.111, -0.002, 0.001,
         0.20000000000000001, -0.14200000000000002, 0.104, 1.,
         3.4878815733152084, -1.2360465761526451, -9.7273767011198302,
         7.6033467723624488, -4.4341916620251585}},
       {{12., 5., 0.040000000000000001, 0.11399999999999999, -0.002, 0.001,
         0.21000000000000002, -0.13800000000000001, 0.106, 1.7,
         3.507734630309681, -0.10725929443853102, -11.173964442048728,
         4.9115731827725284, -4.5926741427458611}},
       {{-12., 5., -0.040000000000000001, 0.11699999999999999, -0.002, 0.001,
         0.22, -0.13400000000000001, 0.108, 0., 3.5662359397124299,
         -1.1281028525027557, -10.023313940766641, 6.2095522846725144,
         -4.8983969129412177}},
       {{12., -5., 0.040000000000000001, 0.12, -0.002, 0.001,
         0.23000000000000001, -0.13, 0.11, 1., 3.2660740055764483,
         -0.21472758179843784, -9.5395412228613363, 4.9791664501184263,
         -3.8455817189427326}},
       {{-12., -5., -0.040000000000000001, 0.123, -0.002, 0.001,
         0.23999999999999999, -0.126, 0.11199999999999999, 1.7,
         3.1563688339359111, -0.81937863140750478, -9.6922342561235943,
         3.4501699884507415, -3.5027599506968614}},
       {{1.0000000000000001e-09, 10., 0.040000000000000001, 0.126, -0.002,
         0.001, 0.25, -0.12200000000000001, 0.11399999999999999, 0.,
         3.3257781405000286, -0.20125321001239374, -9.188396791203207,
         4.5132057890746236, -4.1935775125624826}},
       {{10., 1.0000000000000001e-09, -0.040000000000000001, 0.129, -0.002,
         0.001, 0.26000000000000001, -0.11800000000000001, 0.11599999999999999,
         1., 3.0264904969682083, -0.75850785686428424, -8.7496359753880064,
         3.1820450859924114, -3.2334703818786323}},
       {{21., 28., 0.040000000000000001, 0.13200000000000001, -0.002, 0.001,
         0.27000000000000002, -0.11400000000000002, 0.11799999999999999, 1.7,
         3.1568997146581976, -0.10053005653848433, -8.6425180777449722,
         3.3821482068472148, -3.7132858258580468}},
       {{-21., 28., -0.040000000000000001, 0.13500000000000001, -0.002, 0.001,
         0.28000000000000003, -0.11000000000000001, 0.12, 0.,
         3.0773990668142388, -0.7750842631414766, -8.1913455586000232,
         3.3104226603130127, -3.4972984814316734}}}};
  PunctureFieldWorkspace workspace{};
  const auto tolerance = Approx::custom().epsilon(2.e-12).scale(1.);
  for (const auto& row : reference) {
    tnsr::I<DataVector, 3> coords{size_t{1}};
    for (size_t d = 0; d < 3; ++d) {
      coords.get(d)[0] = row[6 + d];
    }
    Fields result{};
    puncture_field_2(make_not_null(&result), coords,
                     tnsr::I<double, 3>{{row[0], row[1], 0.}},
                     tnsr::I<double, 3>{{row[2], row[3], 0.}},
                     tnsr::I<double, 3>{{row[4], row[5], 0.}}, row[9],
                     &workspace);
    for (size_t i = 0; i < 5; ++i) {
      CHECK(result.data()[i] == tolerance(row[10 + i]));
    }
  }
}

void test_reuse_and_derivatives() {
  PunctureFieldWorkspace workspace{};
  const tnsr::I<double, 3> position{{-6., 8., 0.}};
  const tnsr::I<double, 3> velocity{{0.03, -0.11, 0.}};
  const tnsr::I<double, 3> acceleration{{-0.002, 0.001, 0.}};
  for (const size_t size :
       {size_t{1}, size_t{7}, size_t{7}, size_t{0}, size_t{32}, size_t{1}}) {
    CAPTURE(size);
    tnsr::I<DataVector, 3> coords{size};
    for (size_t i = 0; i < size; ++i) {
      get<0>(coords)[i] = 0.3 + 0.01 * static_cast<double>(i);
      get<1>(coords)[i] = -0.2 + 0.003 * static_cast<double>(i);
      get<2>(coords)[i] = 0.1 + 0.002 * static_cast<double>(i);
    }
    Fields result{}, fresh{};
    puncture_field_2(make_not_null(&result), coords, position, velocity,
                     acceleration, 1.7, &workspace);
    puncture_field_2(make_not_null(&fresh), coords, position, velocity,
                     acceleration, 1.7);
    if (size > 0) {
      CHECK_VARIABLES_APPROX(result, fresh);
    } else {
      CHECK(result.size() == 0);
      CHECK(fresh.size() == 0);
    }
    auto migrated = serialize_and_deserialize(workspace);
    CHECK(migrated.empty());
    puncture_field_2(make_not_null(&fresh), coords, position, velocity,
                     acceleration, 1.7, &migrated);
    if (size == 0) {
      CHECK(fresh.size() == 0);
      CHECK(migrated.empty());
      continue;
    }
    CHECK_VARIABLES_APPROX(result, fresh);

    const double h = 1.e-4;
    const auto evaluate = [&](const double t, const int spatial_direction) {
      auto x = position;
      auto v = velocity;
      auto offset = coords;
      if (spatial_direction < 0) {
        for (size_t d = 0; d < 3; ++d) {
          const double displacement =
              t * velocity.get(d) + 0.5 * t * t * acceleration.get(d);
          x.get(d) += displacement;
          v.get(d) += t * acceleration.get(d);
          offset.get(d) -= displacement;
        }
      } else {
        offset.get(static_cast<size_t>(spatial_direction)) += t;
      }
      Fields shifted{};
      puncture_field_2(make_not_null(&shifted), offset, x, v, acceleration, 1.7,
                       &workspace);
      return get(get<CurvedScalarWave::Tags::Psi>(shifted));
    };
    const auto tolerance = Approx::custom().epsilon(2.e-10).scale(1.);
    for (int direction = -1; direction < 3; ++direction) {
      CAPTURE(direction);
      const DataVector finite_difference =
          (evaluate(-2. * h, direction) - 8. * evaluate(-h, direction) +
           8. * evaluate(h, direction) - evaluate(2. * h, direction)) /
          (12. * h);
      const auto& derivative =
          direction < 0
              ? get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(result))
              : get<Deriv>(result).get(static_cast<size_t>(direction));
      CHECK_ITERABLE_CUSTOM_APPROX(derivative, finite_difference, tolerance);
    }
  }
}

void test_workspace() {
  PunctureFieldWorkspace workspace{};
  CHECK(workspace.empty());
  auto& initial = workspace.buffer(3, 5);
  initial.at(0) = 2.;
  const auto* data = initial.at(0).data();
  auto& reused = workspace.buffer(3, 5);
  CHECK(reused.at(0).data() == data);
  CHECK(reused.at(0)[0] == 2.);
  auto& resized = workspace.buffer(3, 11);
  CHECK(resized.size() == 3);
  CHECK(resized.at(0).size() == 11);
  auto restored = serialize_and_deserialize(workspace);
  CHECK(restored.empty());
  CHECK(restored.buffer(2, 4).at(0).size() == 4);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.CurvedScalarWave.Worldtube.PunctureFieldOrder2",
                  "[Unit][Evolution]") {
  test_reference();
  test_reuse_and_derivatives();
  test_workspace();
}
}  // namespace CurvedScalarWave::Worldtube
