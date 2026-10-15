// Distributed under the MIT License.
// See LICENSE.txt for details.
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/SettleToConstant.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/MomentFeedback.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Framework/TestingFramework.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/Serialization/Serialize.hpp"

SPECTRE_TEST_CASE("Unit.Gh.Worldtube.MomentFeedback", "[Unit][Evolution]") {
  using gh::BoundaryConditions::detail::MomentFeedback;
  const std::string csv = "FeedbackControlTest.csv",
                    bin = "FeedbackGeometryTest.bin";
  {
    std::ofstream out(csv);
    out << std::setprecision(17);
    out << "time,E11,E22,E12,E13,E23,B11,B22,B12,B13,B23\n";
    for (size_t k = 0; k < 6; ++k)
      out << k << ',' << 2. + k << ",-0.3,0.2,0.1,0.4,0.7,0.8,0.9,1.1,1.2\n";
  }
  constexpr size_t n = 6;
  // A constant sheared reference Cholesky triad has an analytically known
  // polar rotation relative to the live identity metric: tan(theta)=s/2.
  constexpr double shear = .4;
  {
    std::ofstream out(bin, std::ios::binary);
    out.write("NPFBG001", 8);
    const uint64_t ns = 6, np = n;
    out.write(reinterpret_cast<const char*>(&ns), 8);
    out.write(reinterpret_cast<const char*>(&np), 8);
    auto real = [&](double x) {
      out.write(reinterpret_cast<const char*>(&x), 8);
    };
    for (size_t p = 0; p < n; ++p)
      real(1. / n);
    for (size_t k = 0; k < 6; ++k) {
      real(static_cast<double>(k));
      for (size_t p = 0; p < n; ++p)
        for (size_t d = 0; d < 3; ++d)
          real(.1 * p + .2 * d + (d == 0 ? .01 * k * k : 0.));
      for (size_t p = 0; p < n; ++p)
        for (size_t i = 0; i < 3; ++i)
          for (size_t j = 0; j < 3; ++j)
            real(i == j ? 1. : (i == 1 and j == 0 ? shear : 0.));
    }
  }
  tnsr::ii<DataVector, 3> metric(n, 0.);
  for (size_t d = 0; d < 3; ++d)
    metric.get(d, d) = 1.;
  tnsr::I<DataVector, 3> coords(n);
  for (size_t p = 0; p < n; ++p)
    for (size_t d = 0; d < 3; ++d)
      coords.get(d)[p] = .1 * p + .2 * d + (d == 0 ? .0625 : 0.);
  DataVector weights(n, 1. / n);
  gr::np::TidalMoments control{
      {{4.5, .7}, {-.3, .8}, {.2, .9}, {.1, 1.1}, {.4, 1.2}}};
  const double c = 1. / sqrt(1. + .25 * shear * shear), s = .5 * shear * c;
  auto rotate = [&](const gr::np::TidalMoments& q) {
    return gr::np::TidalMoments{c * c * q[0] + s * s * q[1] - 2 * c * s * q[2],
                                s * s * q[0] + c * c * q[1] + 2 * c * s * q[2],
                                c * s * (q[0] - q[1]) + (c * c - s * s) * q[2],
                                c * q[3] - s * q[4], s * q[3] + c * q[4]};
  };
  const auto live = rotate(control);
  gr::np::TidalMoments recorded{
      {{1., 2.}, {3., 4.}, {5., 6.}, {7., 8.}, {9., 10.}}},
      pulse{};
  const MomentFeedback feedback(.5, 0., 1., csv, bin);
  const auto restored = serialize_and_deserialize(feedback);
  CHECK(restored.strength(0.) == 0.);
  CHECK(restored.strength(.5) == approx(.25));
  CHECK(restored.strength(1.) == .5);
  CHECK_ITERABLE_APPROX(
      restored.imposed(2.5, live, recorded, pulse, metric, coords, weights),
      recorded);
  auto perturbed = live;
  perturbed[2] += std::complex<double>{.2, -.4};
  pulse[3] = {.1, .3};
  auto expected = recorded;
  expected[2] += std::complex<double>{.1, -.2};
  const auto rotated_pulse = rotate(pulse);
  for (size_t a = 0; a < 5; ++a)
    expected[a] += rotated_pulse[a];
  CHECK_ITERABLE_APPROX(restored.imposed(2.5, perturbed, recorded, pulse,
                                         metric, coords, weights),
                        expected);
  const MomentFeedback zero(0., 0., 1., csv, bin);
  for (size_t a = 0; a < 5; ++a)
    expected[a] = recorded[a] + rotated_pulse[a];
  CHECK_ITERABLE_APPROX(
      zero.imposed(2.5, perturbed, recorded, pulse, metric, coords, weights),
      expected);
  CHECK_THROWS_WITH(
      restored.imposed(6., live, recorded, pulse, metric, coords, weights),
      Catch::Matchers::ContainsSubstring("does not cover"));
  auto wrong_coords = coords;
  wrong_coords.get(0)[0] += .01;
  CHECK_THROWS_WITH(
      restored.imposed(2.5, live, recorded, pulse, metric, wrong_coords,
                       weights),
      Catch::Matchers::ContainsSubstring("coordinates or point ordering"));
  CHECK_THROWS_WITH(MomentFeedback(2., 0., 1., csv, bin),
                    Catch::Matchers::ContainsSubstring("Strength in [0,1]"));
  // Independent gains act on the fitted tensor after alignment. The
  // reference tide and the external pulse remain in both sectors.
  const gr::np::TidalMoments delta_reference{
      {{.2, .3}, {-.4, .5}, {.6, -.7}, {-.8, .9}, {1., -1.1}}};
  const auto delta_live = rotate(delta_reference);
  auto mixed_live = live;
  for (size_t a = 0; a < 5; ++a)
    mixed_live[a] += delta_live[a];
  for (const auto gains :
       std::array<std::array<double, 2>, 6>{{{{0., 0.}},
                                             {{.5, .5}},
                                             {{1., 1.}},
                                             {{1., 0.}},
                                             {{0., 1.}},
                                             {{.25, .75}}}}) {
    CAPTURE(gains);
    const auto selected =
        serialize_and_deserialize(TestHelpers::test_creation<MomentFeedback>(
            "Strength: [" + std::to_string(gains[0]) + ", " +
            std::to_string(gains[1]) +
            "]\n"
            "StartTime: 2\nRampDuration: 1\nMomentsFile: " +
            csv + "\nGeometryFile: " + bin));
    CHECK(selected.strengths(1.) == std::array<double, 2>{{0., 0.}});
    CHECK(selected.strengths(2.) == std::array<double, 2>{{0., 0.}});
    CHECK_ITERABLE_APPROX(
        selected.strengths(2.5),
        (std::array<double, 2>{{.5 * gains[0], .5 * gains[1]}}));
    CHECK(selected.strengths(3.) == gains);
    CHECK(selected.strengths(4.) == gains);
    auto selective_expected = recorded;
    for (size_t a = 0; a < 5; ++a) {
      selective_expected[a] +=
          std::complex<double>{.5 * gains[0] * delta_live[a].real(),
                               .5 * gains[1] * delta_live[a].imag()} +
          rotated_pulse[a];
    }
    CHECK_ITERABLE_APPROX(selected.imposed(2.5, mixed_live, recorded, pulse,
                                           metric, coords, weights),
                          selective_expected);
    if (gains[0] == gains[1]) {
      const auto scalar = TestHelpers::test_creation<MomentFeedback>(
          "Strength: " + std::to_string(gains[0]) +
          "\nStartTime: 2\nRampDuration: 1\n"
          "MomentsFile: " +
          csv + "\nGeometryFile: " + bin);
      CHECK(selected.imposed(2.5, mixed_live, recorded, pulse, metric, coords,
                             weights) == scalar.imposed(2.5, mixed_live,
                                                        recorded, pulse, metric,
                                                        coords, weights));
    }
  }
  for (const auto bad : std::array<std::array<double, 2>, 4>{
           {{{-1., 0.}},
            {{0., 2.}},
            {{std::numeric_limits<double>::quiet_NaN(), 0.}},
            {{0., std::numeric_limits<double>::infinity()}}}}) {
    CHECK_THROWS_WITH(MomentFeedback(bad, 0., 1., csv, bin),
                      Catch::Matchers::ContainsSubstring("Strength in [0,1]"));
  }
  const std::string truncated = "FeedbackTruncatedTest.bin";
  std::filesystem::copy_file(bin, truncated,
                             std::filesystem::copy_options::overwrite_existing);
  std::filesystem::resize_file(truncated, 30);
  CHECK_THROWS_WITH(MomentFeedback(.5, 0., 1., csv, truncated),
                    Catch::Matchers::ContainsSubstring("Truncated or extra"));
  // A rapidly settling map is not resolved by six-point interpolation of
  // the sparse coordinates. The version-2 exact map must still accept the
  // correct face between knots and reject an actual coordinate mismatch.
  using Translation = domain::CoordinateMaps::TimeDependent::Translation<3>;
  using Map = domain::CoordinateMap<Frame::Grid, Frame::Inertial, Translation>;
  register_classes_with_charm<Map, domain::FunctionsOfTime::SettleToConstant>();
  auto map = domain::make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
      Translation{"Shift"});
  domain::FunctionsOfTimeMap fots;
  fots["Shift"] = std::make_unique<domain::FunctionsOfTime::SettleToConstant>(
      std::array<DataVector, 3>{DataVector{0., 0., 0.}, DataVector{1., 0., 0.},
                                DataVector{0., 0., 0.}},
      0., .05);
  tnsr::I<DataVector, 3, Frame::Grid> grid(n);
  for (size_t p = 0; p < n; ++p)
    for (size_t d = 0; d < 3; ++d)
      grid.get(d)[p] = .1 * p + .2 * d;
  const std::string exact_bin = "FeedbackExactGeometryTest.bin";
  {
    std::ofstream out(exact_bin, std::ios::binary);
    out.write("NPFBG002", 8);
    const uint64_t ns = 6, np = n;
    out.write(reinterpret_cast<const char*>(&ns), 8);
    out.write(reinterpret_cast<const char*>(&np), 8);
    auto real = [&](double x) {
      out.write(reinterpret_cast<const char*>(&x), 8);
    };
    for (size_t p = 0; p < n; ++p)
      real(1. / n);
    for (size_t k = 0; k < 6; ++k) {
      real(static_cast<double>(k));
      const auto x = (*map)(grid, static_cast<double>(k), fots);
      for (size_t p = 0; p < n; ++p)
        for (size_t d = 0; d < 3; ++d)
          real(x.get(d)[p]);
      for (size_t p = 0; p < n; ++p)
        for (size_t i = 0; i < 3; ++i)
          for (size_t j = 0; j < 3; ++j)
            real(i == j ? 1. : (i == 1 and j == 0 ? shear : 0.));
    }
    const auto map_data = serialize(map), fots_data = serialize(fots);
    const uint64_t nm = map_data.size(), nf = fots_data.size();
    out.write(reinterpret_cast<const char*>(&nm), 8);
    out.write(reinterpret_cast<const char*>(&nf), 8);
    out.write(map_data.data(), static_cast<std::streamsize>(nm));
    out.write(fots_data.data(), static_cast<std::streamsize>(nf));
    for (size_t p = 0; p < n; ++p)
      for (size_t d = 0; d < 3; ++d)
        real(grid.get(d)[p]);
  }
  const std::string sparse_bin = "FeedbackSparseCoordinatesTest.bin";
  std::filesystem::copy_file(exact_bin, sparse_bin,
                             std::filesystem::copy_options::overwrite_existing);
  std::filesystem::resize_file(sparse_bin, 24 + 8 * n + 6 * (8 + 96 * n));
  {
    std::fstream out(sparse_bin,
                     std::ios::in | std::ios::out | std::ios::binary);
    out.write("NPFBG001", 8);
  }
  const MomentFeedback sparse(0., 0., 1., csv, sparse_bin);
  const gr::np::TidalMoments zero_probe{};
  CHECK_THROWS_WITH(
      sparse.imposed(.0004, live, recorded, zero_probe, metric,
                     (*map)(grid, .0004, fots), weights),
      Catch::Matchers::ContainsSubstring("coordinates or point ordering"));
  std::filesystem::remove(sparse_bin);
  const MomentFeedback exact(0., 0., 1., csv, exact_bin);
  const auto exact_restored = serialize_and_deserialize(exact);
  const auto substep_coords = (*map)(grid, .0004, fots);
  const gr::np::TidalMoments no_pulse{};
  CHECK_ITERABLE_APPROX(exact_restored.imposed(.0004, live, recorded, no_pulse,
                                               metric, substep_coords, weights),
                        recorded);
  auto bad_substep = substep_coords;
  bad_substep.get(1)[0] += 3.e-8;
  CHECK_THROWS_WITH(
      exact_restored.imposed(.0004, live, recorded, no_pulse, metric,
                             bad_substep, weights),
      Catch::Matchers::ContainsSubstring("coordinates or point ordering"));
  std::filesystem::remove(exact_bin);
  std::filesystem::remove(truncated);
  std::filesystem::remove(csv);
  std::filesystem::remove(bin);
}
