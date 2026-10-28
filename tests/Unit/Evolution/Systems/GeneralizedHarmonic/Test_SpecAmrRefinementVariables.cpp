// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <limits>

#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.RefinementVariables",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  using Kind = RefinementVariables::Kind;
  using Spectral::Basis;
  using Spectral::Quadrature;
  const auto up = amr::Flag::IncreaseResolution;
  const auto down = amr::Flag::DecreaseResolution;
  const auto hold = amr::Flag::DoNothing;

  // The layout is element-local. Different polynomial bases can coexist with
  // S2 elements under the same controller, with independent mode counts.
  const Mesh<3> polynomial{
      {{8, 10, 12}},
      {{Basis::Legendre, Basis::Chebyshev, Basis::Legendre}},
      {{Quadrature::GaussLobatto, Quadrature::Gauss, Quadrature::Gauss}}};
  const RefinementVariables cube{polynomial};
  CHECK(cube == RefinementVariables{});
  CHECK(cube.size() == 3);
  CHECK(cube.mode_counts(polynomial) == polynomial.extents().indices());
  const std::array<amr::Flag, 3> polynomial_flags{{up, hold, down}};
  CHECK(cube.logical_flags(polynomial_flags) == polynomial_flags);
  CHECK(cube.independent_flags(polynomial_flags) == polynomial_flags);
  for (size_t d = 0; d < 3; ++d) {
    CHECK(cube.kind(d) == Kind::Polynomial);
    CHECK(cube.variable_for_dimension(d) == d);
  }

  const Mesh<3> shell{{{10, 14, 27}},
                      Spectral::bases::spherical_shell<>,
                      Spectral::quadratures::spherical_shell<>};
  const RefinementVariables sphere{shell};
  CHECK(sphere.size() == 2);
  CHECK(sphere.kind(0) == Kind::Polynomial);
  CHECK(sphere.kind(1) == Kind::SphericalHarmonic);
  CHECK(sphere.variable_for_dimension(0) == 0);
  CHECK(sphere.variable_for_dimension(1) == 1);
  CHECK(sphere.variable_for_dimension(2) == 1);
  CHECK(sphere.mode_counts(shell) == (std::array<size_t, 3>{{10, 14, 0}}));
  CHECK(serialize_and_deserialize(sphere) == sphere);
  CHECK(sphere.logical_flags({{up, down, hold}}) ==
        (std::array<amr::Flag, 3>{{up, down, down}}));
  CHECK(sphere.independent_flags({{up, down, down}}) ==
        (std::array<amr::Flag, 3>{{up, down, hold}}));
  CHECK_THROWS_WITH(sphere.independent_flags({{up, up, down}}),
                    Catch::Matchers::ContainsSubstring("single angular"));
  CHECK_THROWS_WITH(sphere.logical_flags({{up, amr::Flag::Split, hold}}),
                    Catch::Matchers::ContainsSubstring("Cannot h refine"));
  CHECK_THROWS_WITH(
      RefinementVariables(
          Mesh<3>{{{8, 8, 9}},
                  {{Basis::Legendre, Basis::Legendre, Basis::Fourier}},
                  {{Quadrature::GaussLobatto, Quadrature::GaussLobatto,
                    Quadrature::Equiangular}}}),
      Catch::Matchers::ContainsSubstring("Unsupported basis combination"));

  // Feed synthetic spectra into the common controller without claiming that
  // GH spherical diagnostics are implemented. One angular history controls
  // both mesh directions, and inactive storage must never enter decisions.
  Diagnostics data{};
  data.refinement_variables = sphere;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto& group : data.spectra) {
    group = {{{1.e-5, 2.e-5, 0.5, 0., false},
              {2.e-4, 4.e-4, 0.5, 0., false},
              {nan, nan, nan, nan, false}}};
  }
  for (auto& group : data.constraints) {
    group = {{1.e-5, 2.e-3, nan}};
  }
  Settings settings{1.e-4, {{1.e-3, 1.e-3}}, true, 7, 2, -1., 3.};
  History history{};
  auto counts = sphere.mode_counts(shell);
  // First evolution call can add angular resolution, but cannot coarsen.
  CHECK(decide_p(make_not_null(&history), data, counts, 0., settings) ==
        (std::array<amr::Flag, 3>{{hold, up, up}}));
  CHECK(history.refinement_variables == sphere);
  CHECK(history.log10_target_shift[0] == 0.);
  CHECK(history.log10_target_shift[1] < 0.);
  CHECK(history.log10_target_shift[2] == 0.);
  ++counts[1];  // A realized l_max change, not a separate phi change.
  CHECK(decide_p(make_not_null(&history), data, counts, 1., settings) ==
        make_array<3>(hold));
  // The initial call-count gate has elapsed, but only the angular variable
  // has a realized-change time and remains blocked by its elapsed-time gate.
  CHECK(decide_p(make_not_null(&history), data, counts, 2., settings) ==
        (std::array<amr::Flag, 3>{{down, hold, hold}}));
  CHECK(history.last_change[1] == 0);
  CHECK(history.last_change_time[1] == 0.);
  CHECK(std::isinf(history.last_change_time[2]));
  auto restored = serialize_and_deserialize(history);
  CHECK(restored.refinement_variables == sphere);
  CHECK(restored.extents == counts);
  const auto calls = restored.calls;
  CHECK(decide_p(make_not_null(&restored), data, counts, 2., settings) ==
        history.decision);
  CHECK(restored.calls == calls);
  CHECK(decide_p(make_not_null(&restored), data, counts, 3., settings) ==
        decide_p(make_not_null(&history), data, counts, 3., settings));

  // Startup uses independent extents for recurring-mesh detection too.
  history = History{};
  history.initializing = true;
  for (auto& group : data.spectra) {
    group[0] = group[1] = Spectrum{1.e-7, 2.e-7, 0.3, 0., false};
  }
  for (auto& group : data.constraints) {
    group[0] = group[1] = 0.;
  }
  CHECK(decide_p(make_not_null(&history), data, counts, 0., settings) ==
        make_array<3>(down));
  CHECK(history.initial_meshes.size() == 1);
  --counts[1];
  CHECK(decide_p(make_not_null(&history), data, counts, 0., settings) ==
        make_array<3>(down));
  ++counts[1];
  CHECK(decide_p(make_not_null(&history), data, counts, 0., settings) ==
        make_array<3>(hold));
  CHECK(history.initial_cycle_detected);
  CHECK(history.calls == 0);
  ProjectHistory::apply(make_not_null(&history), 0);
  CHECK(history.refinement_variables == cube);
}
