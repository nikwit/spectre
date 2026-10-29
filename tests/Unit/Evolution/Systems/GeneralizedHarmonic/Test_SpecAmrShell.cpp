// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Domain/Creators/SphericalShells.hpp"
#include "Domain/ElementMap.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmrReference/ReferenceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecSpectral.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Hypercube.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/None.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
using GhTags = ylm::TensorYlm::filter_detail::gh_spacetime_vars_list;
using Metric = gr::Tags::SpacetimeMetric<DataVector, 3>;
using Pi = gh::Tags::Pi<DataVector, 3>;
using Phi = gh::Tags::Phi<DataVector, 3>;
using FilterTag = Filters::Tags::SpectralFilter<3, GhTags>;
using GridJac = domain::Tags::Jacobian<3, Frame::Grid, Frame::Inertial>;
using Inverse =
    domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>;

Mesh<3> shell_mesh(const size_t nr, const size_t ell) {
  return {{{nr, ell + 1, 2 * ell + 1}},
          {{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
            Spectral::Basis::SphericalHarmonic}},
          {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::Gauss,
            Spectral::Quadrature::Equiangular}}};
}

GridJac::type identity(const size_t n) {
  GridJac::type jac(n, 0.0);
  for (size_t i = 0; i < 3; ++i) {
    jac.get(i, i) = 1.0;
  }
  return jac;
}

Filters::SphericalShell<GhTags> shell_filter(
    const std::optional<size_t> angular_power,
    const std::optional<size_t> radial_power, const bool enable = true,
    const bool volume = true) {
  return {2,
          angular_power,
          radial_power,
          enable,
          std::nullopt,
          false,
          false,
          volume ? std::optional<size_t>{1} : std::nullopt,
          1};
}

struct Metavariables {
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<amr::Criterion, tmpl::list<gh::SpecSpectral>>>;
  };
};
}  // namespace

// [[TimeOut, 30]]
SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.LocalFilter",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  // Query the real local filter and compare with its action on individual
  // modes, including a different resolution for the removal forecast.
  for (const auto basis :
       {Spectral::Basis::Legendre, Spectral::Basis::Chebyshev}) {
    const Mesh<3> mesh{{{9, 7, 6}}, basis, Spectral::Quadrature::GaussLobatto};
    const Filters::Hypercube<3, GhTags> filter{36,    true, std::nullopt, false,
                                               false, 1,    std::nullopt};
    const auto response = filter.volume_filter_response(mesh).value();
    const auto info = filtered_modes(mesh, filter);
    for (size_t d = 0; d < 3; ++d) {
      double count = 0.;
      for (size_t k = 0; k < mesh.extents(d); ++k) {
        ModalVector coefficients(mesh.number_of_grid_points(), 0.0);
        size_t stride = 1;
        for (size_t j = 0; j < d; ++j) {
          stride *= mesh.extents(j);
        }
        coefficients[k * stride] = 1.0;
        Variables<GhTags> fields(mesh.number_of_grid_points(), 0.0);
        get<Metric>(fields).get(0, 0) =
            to_nodal_coefficients(coefficients, mesh);
        filter.apply_in_volume(make_not_null(&fields), mesh, std::nullopt,
                               std::nullopt);
        const auto actual =
            to_modal_coefficients(get<Metric>(fields).get(0, 0), mesh);
        CHECK(actual[k * stride] == approx(response[d][k]));
        count += 1.0 - response[d][k];
      }
      CHECK(info.current[d] == approx(count));
      // Effective count can depend on N, even for fixed filter parameters.
      CHECK(info.after_p_down[d] > 0.);
    }
    const Filters::Hypercube<3, GhTags> disabled{
        36, false, std::nullopt, false, false, 1, 1};
    CHECK(filtered_modes(mesh, disabled).current == make_array<3>(0.));
  }
  for (const auto angular_power :
       {std::optional<size_t>{}, std::optional<size_t>{4}}) {
    const auto filter = shell_filter(angular_power, 6);
    for (const size_t ell_max : {size_t{6}, size_t{8}, size_t{6}}) {
      const auto mesh = shell_mesh(7, ell_max);
      const auto info = filtered_modes(mesh, filter);
      const auto response = filter.volume_filter_response(mesh).value();
      CHECK(response[1] == response[2]);
      CHECK(info.current[2] == 0.);
      if (not angular_power.has_value()) {
        CHECK(info.current[1] == 2.);
        CHECK(info.after_p_down[1] == 2.);
      }
      const auto& spherepack = ylm::get_spherepack_cache(ell_max);
      for (size_t ell = 0; ell <= ell_max; ++ell) {
        CAPTURE(angular_power, ell_max, ell);
        DataVector coefficients(spherepack.spectral_size(), 0.0);
        ylm::SpherepackIterator it(ell_max, ell_max);
        it.set(ell, 0);
        coefficients[it()] = 1.0;
        const auto angular = spherepack.spec_to_phys(coefficients);
        Variables<GhTags> fields(mesh.number_of_grid_points(), 0.0);
        for (size_t a = 0; a < angular.size(); ++a) {
          for (size_t r = 0; r < mesh.extents(0); ++r) {
            get<Metric>(fields).get(0, 0)[r + mesh.extents(0) * a] = angular[a];
          }
        }
        const auto jac = identity(mesh.number_of_grid_points());
        InverseJacobian<DataVector, 3, Frame::Grid, Frame::Inertial> inv(
            mesh.number_of_grid_points(), 0.0);
        for (size_t d = 0; d < 3; ++d) {
          inv.get(d, d) = 1.;
        }
        filter.apply_in_volume(make_not_null(&fields), mesh, inv, jac);
        DataVector expected(mesh.number_of_grid_points());
        for (size_t i = 0; i < expected.size(); ++i) {
          expected[i] = angular[i / mesh.extents(0)] * response[1][ell];
        }
        CHECK(max(abs(get<Metric>(fields).get(0, 0) - expected)) < 1.e-12);
      }
      CHECK(filtered_modes(mesh, shell_filter(4, 6, false)).current ==
            make_array<3>(0.));
      CHECK(filtered_modes(mesh, shell_filter(4, 6, true, false)).current ==
            make_array<3>(0.));
      const Filters::None<3, GhTags> no_filter{std::nullopt};
      CHECK(filtered_modes(mesh, no_filter).current == make_array<3>(0.));
      CHECK(filtered_modes(mesh, no_filter).after_p_down == make_array<3>(0.));
    }
  }
  const Filters::SphericalShell<GhTags> at_filter_limit{
      6, std::nullopt, std::nullopt, true, std::nullopt, false, false,
      1, std::nullopt};
  const auto limited = filtered_modes(shell_mesh(7, 6), at_filter_limit);
  CHECK(limited.current[1] == 6.);
  CHECK_FALSE(limited.can_decrease[1]);
  const Filters::SphericalShell<GhTags> disabled_at_filter_limit{
      6, std::nullopt, std::nullopt, false, std::nullopt, false, false,
      1, std::nullopt};
  const auto disabled =
      filtered_modes(shell_mesh(7, 6), disabled_at_filter_limit);
  CHECK(disabled.current == make_array<3>(0.));
  CHECK(disabled.can_decrease[1]);
  DataVector spectrum(12);
  for (size_t k = 0; k < spectrum.size(); ++k) {
    spectrum[k] = pow(10., -.4 * k) * (1. + .03 * sin(k));
  }
  CHECK(diagnose_spec_power(spectrum, 2., 1.).next_error ==
        approx(diagnose_spec_power(spectrum, 2.).error));
  CHECK(diagnose_spec_power(spectrum, 2., 3.).next_error >
        diagnose_spec_power(spectrum, 2.).next_error);
}

SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.ShellDecisions",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  const auto mesh = shell_mesh(9, 6);
  const size_t n = mesh.number_of_grid_points();
  const auto& spherepack = ylm::get_spherepack_cache(6);
  const domain::creators::SphericalShells creator{
      4., 7., 0, 9, 6, {}, domain::CoordinateMaps::Distribution::Linear};
  const auto domain = creator.create_domain();
  const ElementMap<3, Frame::Inertial> map{
      ElementId<3>{0}, domain.blocks()[0].stationary_map().get_clone()};
  const auto logical = logical_coordinates(mesh);
  const auto coords = map(logical);
  const auto inverse = map.inv_jacobian(logical);
  Metric::type metric(n, 0.);
  metric.get(0, 0) = -1.;
  for (size_t d = 1; d < 4; ++d) {
    metric.get(d, d) = 1.;
  }
  Phi::type phi(n, 0.);
  const auto filter = shell_filter(std::nullopt, 36);
  const auto info = filtered_modes(mesh, filter);
  gh::power_monitor::CartToSphereMatrices cache{};
  const gh::gauges::Harmonic gauge{};
  for (const bool resolve_radial : {true, false}) {
    ModalVector radial_modes(9);
    DataVector angular_modes(spherepack.spectral_size(), 0.);
    for (size_t k = 0; k < 9; ++k) {
      radial_modes[k] = pow(10., -(resolve_radial ? .25 : 1.5) * k);
    }
    for (size_t ell = 0; ell <= 6; ++ell) {
      ylm::SpherepackIterator it(6, 6);
      it.set(ell, 0);
      angular_modes[it()] = pow(10., -(resolve_radial ? 1.5 : .25) * ell);
    }
    const auto radial =
        to_nodal_coefficients(radial_modes, mesh.slice_through(0));
    const auto angular = spherepack.spec_to_phys(angular_modes);
    Pi::type pi(n, 0.);
    for (size_t a = 0; a < angular.size(); ++a) {
      for (size_t r = 0; r < radial.size(); ++r) {
        pi.get(0, 0)[r + radial.size() * a] = radial[r] * angular[a];
      }
    }
    const auto data = diagnose_spec(
        metric, pi, phi, mesh, inverse, identity(n), coords, 0., gauge,
        Scalar<DataVector>(n, 1.), info, make_not_null(&cache));
    DataVector expected_radial(9), expected_angular(7);
    for (size_t k = 0; k < 9; ++k) {
      expected_radial[k] = std::abs(radial_modes[k]) *
                           sqrt(sum(square(angular)) / (40. * angular.size()));
    }
    for (size_t ell = 0; ell <= 6; ++ell) {
      ylm::SpherepackIterator it(6, 6);
      it.set(ell, 0);
      expected_angular[ell] =
          std::abs(angular_modes[it()]) *
          sqrt(sum(square(radial)) /
               (radial.size() *
                static_cast<double>(
                    spec_amr_reference::angular_counts[1][ell] +
                    spec_amr_reference::angular_counts[2][ell])));
    }
    for (size_t v = 0; v < 2; ++v) {
      const auto expected =
          diagnose_spec_power(v == 0 ? expected_radial : expected_angular,
                              info.current[v], info.after_p_down[v]);
      CHECK(data.spectra[1][v].error == approx(expected.error));
      CHECK(data.spectra[1][v].next_error == approx(expected.next_error));
    }
    History history{};
    history.refinement_variables = data.refinement_variables;
    history.calls = 3;
    history.extents = data.refinement_variables.mode_counts(mesh);
    history.constraints = data.constraints;
    const Settings settings{1.e-4, {{1.e6, 1.e6}}, false, 5, 1, 0.};
    const auto flags =
        decide_p(make_not_null(&history), data, history.extents, 0., settings);
    CHECK(flags[0] == (resolve_radial ? amr::Flag::IncreaseResolution
                                      : amr::Flag::DecreaseResolution));
    CHECK(flags[1] == (resolve_radial ? amr::Flag::DecreaseResolution
                                      : amr::Flag::IncreaseResolution));
    CHECK(flags[2] == flags[1]);
  }
}

// [[TimeOut, 20]]
SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.ShellDiagnostics",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  register_factory_classes_with_charm<Metavariables>();
  const gh::gauges::Harmonic gauge{};
  const auto filter = shell_filter(std::nullopt, 36);
  gh::power_monitor::CartToSphereMatrices cache{};
  // Exercise l_max changes in both directions with the same cache, and a
  // radial extent change. Compare against fresh transforms each time.
  for (const auto mesh :
       {shell_mesh(7, 6), shell_mesh(8, 8), shell_mesh(7, 6)}) {
    const size_t n = mesh.number_of_grid_points();
    const domain::creators::SphericalShells creator{
        4.,
        7.,
        0,
        mesh.extents(0),
        mesh.extents(1) - 1,
        {},
        domain::CoordinateMaps::Distribution::Linear};
    const auto domain = creator.create_domain();
    const ElementMap<3, Frame::Inertial> map{
        ElementId<3>{0}, domain.blocks()[0].stationary_map().get_clone()};
    const auto logical = logical_coordinates(mesh);
    auto coords = map(logical);
    const auto jac = map.jacobian(logical);
    auto inv = map.inv_jacobian(logical);
    Metric::type metric(n, 0.0);
    metric.get(0, 0) = -1.;
    for (size_t d = 1; d < 4; ++d) {
      metric.get(d, d) = 1.;
    }
    Pi::type pi(n, 0.0);
    Phi::type phi(n, 0.0);
    const Scalar<DataVector> gamma2(n, 1.0);
    const auto info = filtered_modes(mesh, filter);
    const auto flat =
        diagnose_spec(metric, pi, phi, mesh, inv, identity(n), coords, 0.,
                      gauge, gamma2, info, make_not_null(&cache));
    CHECK(flat.refinement_variables.size() == 2);
    for (const auto& group : flat.constraints) {
      CHECK(max(abs(DataVector{group})) < 2.e-6);
      CHECK(group[2] == 0.);
    }
    // A known three-index constraint tests the actual projection, component
    // sum, grid-map normalization and final I1/S2 aggregation together.
    phi.get(0, 0, 0) = 1.e-3;
    phi.get(1, 0, 0) = 2.e-3;
    phi.get(2, 0, 0) = -3.e-3;
    for (const bool moving : {false, true}) {
      auto grid_jac = identity(n);
      if (moving) {
        grid_jac.get(0, 1) = .3;
        grid_jac.get(1, 2) = -.2;
        grid_jac.get(2, 0) = .4;
      }
      auto composed = jac;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t d = 0; d < 3; ++d) {
          composed.get(i, d) = 0.;
          for (size_t j = 0; j < 3; ++j) {
            composed.get(i, d) += grid_jac.get(i, j) * jac.get(j, d);
          }
        }
      }
      inv = determinant_and_inverse(composed).second;
      const auto data =
          diagnose_spec(metric, pi, phi, mesh, inv, grid_jac, coords, 0., gauge,
                        gamma2, info, make_not_null(&cache));
      std::array<double, 3> rms{};
      for (size_t d = 0; d < 3; ++d) {
        DataVector norm(n, 0.), value(n, 0.);
        for (size_t i = 0; i < 3; ++i) {
          value += composed.get(i, d) * phi.get(i, 0, 0);
          for (size_t j = 0; j < 3; ++j) {
            norm += square(grid_jac.get(i, j) * jac.get(j, d));
          }
        }
        rms[d] = sqrt(sum(square(value) / norm) / n);
      }
      CHECK_ITERABLE_APPROX(
          data.constraints[1],
          aggregate_constraints(rms, RefinementVariables{mesh}));
      gh::power_monitor::CartToSphereMatrices fresh{};
      const auto uncached =
          diagnose_spec(metric, pi, phi, mesh, inv, grid_jac, coords, 0., gauge,
                        gamma2, info, make_not_null(&fresh));
      for (size_t group = 0; group < 2; ++group) {
        for (size_t v = 0; v < 2; ++v) {
          CHECK(data.spectra[group][v].error ==
                uncached.spectra[group][v].error);
        }
      }
    }
    // The wrapper uses independent extents in its cache and enforces phi's
    // floor through the one S2 decision. Startup can raise a mesh below it.
    History history{};
    history.initializing = true;
    tnsr::I<DataVector, 3, Frame::Grid> grid_coords(n, 0.0);
    for (size_t d = 0; d < 3; ++d) {
      grid_coords.get(d) = coords.get(d);
    }
    auto box = db::create<
        tmpl::list<domain::Tags::Mesh<3>, Metric, Pi, Phi, Inverse, GridJac,
                   domain::Tags::Coordinates<3, Frame::Inertial>,
                   domain::Tags::Coordinates<3, Frame::Grid>, ::Tags::Time,
                   gh::gauges::Tags::GaugeCondition, gh::Tags::ConstraintGamma2,
                   gh::spec_amr::HistoryTag, FilterTag>>(
        mesh, metric, pi, phi, inv, identity(n), coords, grid_coords, 0.,
        std::unique_ptr<gh::gauges::GaugeCondition>(new gh::gauges::Harmonic{}),
        gamma2, history,
        FilterTag::type(new Filters::SphericalShell<GhTags>(filter)));
    Parallel::GlobalCache<Metavariables> global_cache{};
    const gh::SpecSpectral criterion{
        1.e-4,
        {{1., 1.}},
        false,
        gh::SpecSpectral::MinimumModes::type{std::array<size_t, 3>{{5, 7, 21}}},
        1,
        0.,
        0.,
        false};
    const auto flags = criterion(box, global_cache, ElementId<3>{0});
    CHECK(flags[1] == amr::Flag::IncreaseResolution);
    CHECK(flags[2] == flags[1]);
    CHECK(serialize_and_deserialize(criterion)(box, global_cache,
                                               ElementId<3>{0}) == flags);
    auto observation = make_observation_box<tmpl::list<>>(make_not_null(&box));
    CHECK(criterion.evaluate_with_history(make_not_null(&observation),
                                          global_cache,
                                          ElementId<3>{0}) == flags);
    const auto& saved = db::get<gh::spec_amr::HistoryTag>(box);
    CHECK(saved.transform_cache != nullptr);
    CHECK(saved.transform_cache->ell_max == mesh.extents(1) - 1);
    const auto restored = serialize_and_deserialize(saved);
    CHECK(restored.transform_cache == nullptr);
    CHECK(restored.refinement_variables == saved.refinement_variables);
  }
}

// A physical rescaling is relevant for q8: Phi and C_iab scale as 1/M,
// C_ia as 1/M^2. SpEC's final square root makes the AMR diagnostic scale as
// 1/sqrt(M) and 1/M respectively. There is no compensating mass factor in
// its lev targets. Exercise both masses and a nonidentity grid map.
// [[TimeOut, 20]]
SPECTRE_TEST_CASE("Unit.GeneralizedHarmonic.SpecAmr.UnequalMassConstraints",
                  "[Unit][Evolution]") {
  using namespace gh::spec_amr;
  const auto mesh = shell_mesh(7, 6);
  const size_t n = mesh.number_of_grid_points();
  const auto filter = shell_filter(std::nullopt, 36);
  const auto info = filtered_modes(mesh, filter);
  const gh::gauges::Harmonic gauge{};
  gh::power_monitor::CartToSphereMatrices cache{};
  for (const bool moving : {false, true}) {
    std::array<std::array<double, 3>, 2> reference{};
    for (const double mass : {1., .5, 8. / 9., 1. / 9.}) {
      CAPTURE(moving, mass);
      const domain::creators::SphericalShells creator{
          4. * mass, 7. * mass, 0, mesh.extents(0), mesh.extents(1) - 1,
          {}, domain::CoordinateMaps::Distribution::Linear};
      const auto domain = creator.create_domain();
      const ElementMap<3, Frame::Inertial> map{
          ElementId<3>{0}, domain.blocks()[0].stationary_map().get_clone()};
      const auto logical = logical_coordinates(mesh);
      const auto jac = map.jacobian(logical);
      auto grid_jac = identity(n);
      if (moving) {
        grid_jac.get(0, 1) = .3;
        grid_jac.get(1, 2) = -.2;
        grid_jac.get(2, 0) = .4;
      }
      auto composed = jac;
      for (size_t i = 0; i < 3; ++i) {
        for (size_t d = 0; d < 3; ++d) {
          composed.get(i, d) = 0.;
          for (size_t j = 0; j < 3; ++j) {
            composed.get(i, d) += grid_jac.get(i, j) * jac.get(j, d);
          }
        }
      }
      const auto inv = determinant_and_inverse(composed).second;
      Metric::type metric(n, 0.);
      metric.get(0, 0) = -1.;
      for (size_t d = 1; d < 4; ++d) {
        metric.get(d, d) = 1.;
      }
      Pi::type pi(n, 0.);
      Phi::type phi(n, 0.);
      phi.get(0, 0, 0) = 1.e-3 / mass;
      phi.get(1, 0, 0) = 2.e-3 / mass;
      phi.get(2, 0, 0) = -3.e-3 / mass;
      const Scalar<DataVector> gamma2(n, 1. / mass);
      const auto data = diagnose_spec(metric, pi, phi, mesh, inv, grid_jac,
                                     map(logical), 0., gauge, gamma2, info,
                                     make_not_null(&cache));
      if (mass == 1.) {
        reference = data.constraints;
      }
      for (size_t v = 0; v < 2; ++v) {
        CHECK(reference[0][v] > 0.);
        CHECK(reference[1][v] > 0.);
        CHECK(data.constraints[0][v] * mass == approx(reference[0][v]));
        CHECK(data.constraints[1][v] * sqrt(mass) == approx(reference[1][v]));
      }
    }
  }
}
