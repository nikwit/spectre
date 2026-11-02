// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <string>

#include "DataStructures/ModalVector.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmrReference/ReferenceData.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Constraints.hpp"

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.SpECReference",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  // Values come from the *SpEC* diagnostic and linear-fit arithmetic, compiled
  // by GenerateReference.py. No SpEC checkout is required to run this test.
  // The agreement cases cover multiple filter counts and genuine pile-up.
  for (const auto& test : spec_amr_reference::spectra) {
    CAPTURE(test.name);
    const auto spec =
        diagnose_spec_power(DataVector{test.power}, test.filtered_modes);
    CHECK(spec.degenerate == test.degenerate);
    if (test.degenerate) {
      CHECK(spec.error == 0.);
      CHECK(spec.next_error == 0.);
    } else {
      CHECK(spec.error == approx(pow(10., test.expected[0])));
      CHECK(spec.next_error == approx(pow(10., test.expected[1])));
      CHECK(spec.convergence == approx(test.expected[2]));
      CHECK(spec.pile_up == approx(test.expected[3]));
    }
    const auto actual =
        diagnose_power(DataVector{test.power}, test.filtered_modes);
    const std::string name{test.name};
    if (name == "smooth" or name == "unfiltered" or name == "three_filtered" or
        name == "pile_up") {
      CHECK(actual.error == approx(pow(10., test.expected[0])));
      CHECK(actual.next_error == approx(pow(10., test.expected[1])));
      CHECK(actual.convergence == approx(test.expected[2]));
      CHECK(actual.pile_up == approx(test.expected[3]));
      CHECK(actual.degenerate == test.degenerate);
    } else if (name == "zero") {
      CHECK(actual.degenerate);
      CHECK(test.degenerate);
      CHECK(actual.error == 0.);
    } else {
      // Document, rather than silently "fix", the inherited differences.
      // Frozen values below were measured with the polynomial implementation
      // before the independent-variable refactor. See SpecAmrReference/README.
      if (name == "flat") {
        CHECK(test.degenerate);
        CHECK_FALSE(actual.degenerate);
        CHECK(actual.error == 1.);
      } else if (name == "parity_gaps") {
        CHECK(actual.error == approx(1.0045651363723165e-4));
        CHECK(actual.next_error == approx(1.2455767135927177e-4));
        CHECK(actual.error != approx(pow(10., test.expected[0])));
      } else if (name == "leading_zero") {
        CHECK(actual.error == approx(4.3565897335306987e-4));
        CHECK(actual.next_error == approx(1.0943258639309686e-3));
        CHECK(actual.error != approx(pow(10., test.expected[0])));
      } else if (name == "resolved_polynomial") {
        CHECK(actual.error == 0.);
        CHECK(actual.next_error == 0.);
        CHECK_FALSE(test.degenerate);
        CHECK(pow(10., test.expected[0]) == approx(.2));
      } else {
        REQUIRE(name == "roundoff_tail");
        CHECK(actual.error == 0.);
        CHECK(actual.next_error == approx(1.3283600926017064e-8));
        CHECK(pow(10., test.expected[0]) > 0.);
      }
    }
  }

  // SpEC averages over tensor components and transverse NODAL points after a
  // single-direction transform. Freeze the current full-modal averaging too:
  // for two-point LGL transverse grids it differs by exactly a factor of two.
  const Mesh<3> mesh{{{4, 2, 2}},
                     Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  DataVector polynomial_power(4, 0.);
  DataVector spec_power(4, 0.);
  for (const auto& coefficients : spec_amr_reference::coefficient_components) {
    DataVector component(mesh.number_of_grid_points());
    for (size_t stripe = 0; stripe < 4; ++stripe) {
      ModalVector modes(4);
      for (size_t k = 0; k < 4; ++k) {
        modes[k] = coefficients[stripe * 4 + k];
      }
      const auto nodal = to_nodal_coefficients(modes, mesh.slice_through(0));
      for (size_t k = 0; k < 4; ++k) {
        component[stripe * 4 + k] = nodal[k];
      }
    }
    polynomial_power +=
        square(PowerMonitors::power_monitors(component, mesh)[0]);
    spec_power += square(gh::spec_amr::polynomial_power(component, mesh)[0]);
  }
  polynomial_power = sqrt(polynomial_power / 2.);
  spec_power = sqrt(spec_power / 2.);
  for (size_t k = 0; k < 4; ++k) {
    CHECK(spec_power[k] == approx(spec_amr_reference::i1_power[k]));
    CHECK(2. * polynomial_power[k] == approx(spec_amr_reference::i1_power[k]));
  }

  // Isolate the topological grouping from the earlier projection and nodal
  // RMS. This checkout of SpEC takes another square root, even for I1. The
  // current controller uses the directional nodal RMS directly. Do not change
  // tolerance semantics as a side effect of adding a variable layout.
  for (const auto& test : spec_amr_reference::constraints) {
    CAPTURE(test.name);
    const size_t n = 3;
    Jacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial> jac(n, 0.);
    tnsr::iaa<DataVector, 3> constraint(n, 0.);
    for (size_t d = 0; d < 3; ++d) {
      jac.get(d, d) = 1.;
      constraint.get(d, 0, 0) = test.directional_rms[d];
    }
    tnsr::i<DataVector, 3, Frame::ElementLogical> normalization{}, projected{};
    DataVector buffer{};
    amr::Criteria::Constraints_detail::normalization_factor_square(
        make_not_null(&normalization), jac);
    amr::Criteria::Constraints_detail::logical_constraints(
        make_not_null(&projected), make_not_null(&buffer), constraint, jac,
        normalization);
    std::array<double, 3> rms{};
    for (size_t d = 0; d < 3; ++d) {
      rms[d] = sqrt(sum(square(projected.get(d))) / n);
      CHECK(rms[d] == approx(test.directional_rms[d]));
    }
    if (test.spec_aggregated.size() == 3) {
      const auto combined = aggregate_constraints(rms, RefinementVariables{});
      for (size_t d = 0; d < 3; ++d) {
        CHECK(combined[d] == approx(test.spec_aggregated[d]));
      }
      for (size_t d = 0; d < 3; ++d) {
        CHECK(sqrt(rms[d]) == approx(test.spec_aggregated[d]));
        CHECK(rms[d] != approx(test.spec_aggregated[d]));
      }
    } else {
      const Mesh<3> shell{
          {{5, 5, 9}},
          {{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
            Spectral::Basis::SphericalHarmonic}},
          {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::Gauss,
            Spectral::Quadrature::Equiangular}}};
      const auto combined =
          aggregate_constraints(rms, RefinementVariables{shell});
      CHECK(combined[0] == approx(test.spec_aggregated[0]));
      CHECK(combined[1] == approx(test.spec_aggregated[1]));
      CHECK(combined[2] == 0.);
      CHECK(sqrt(rms[0]) == approx(test.spec_aggregated[0]));
      CHECK(sqrt((rms[1] + rms[2]) / 2.) == approx(test.spec_aggregated[1]));
    }
  }

  // Execute the common p controller against permission expressions extracted
  // from SetAddDropLimits. Scope: one element, no neighbor pile-up, no history
  // trend, no shifted target, and no cooldown. An add permission still needs
  // LargeTrunc in ComputeGridChangesWithOrWithoutBuffer.
  for (const auto& test : spec_amr_reference::permissions) {
    CAPTURE(test.error, test.constraint_ratio, test.pile_up, test.force);
    Diagnostics data{};
    for (auto& group : data.spectra) {
      group.fill(
          Spectrum{test.error, test.next_error, .5, test.pile_up, false});
    }
    for (auto& group : data.constraints) {
      group.fill(test.constraint_ratio * 1.e-3);
    }
    History history{};
    history.calls = 3;
    history.extents = {{12, 12, 12}};
    history.constraints = data.constraints;
    Settings settings{1.e-4, {{1.e-3, 1.e-3}}, test.force, 7, 2, 0.};
    const auto expected =
        test.allow_add and test.error > settings.spectral_target
            ? amr::Flag::IncreaseResolution
        : test.allow_drop ? amr::Flag::DecreaseResolution
                          : amr::Flag::DoNothing;
    CHECK(decide_p(make_not_null(&history), data, history.extents, 0.,
                   settings) == make_array<3>(expected));
  }
}
