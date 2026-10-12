// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/ScalarWave/ChalmersConstraintGuard.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPoints.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"

namespace {
template <size_t Dim>
void check_polynomial() {
  const Mesh<Dim> mesh{8, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  constexpr double amplitude = 1.e-5;
  ModalVector modes(mesh.number_of_grid_points(), 0.0);
  size_t stride = 1;
  for (size_t d = 0; d < Dim; ++d) {
    modes[7 * stride] = amplitude;
    stride *= 8;
  }
  const Scalar<DataVector> psi{to_nodal_coefficients(modes, mesh)};
  for (const double scale : {0.5, 2.0}) {
    InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>
        inverse(mesh.number_of_grid_points(), 0.0);
    for (size_t d = 0; d < Dim; ++d) {
      inverse.get(d, d) = scale;
    }
    auto phi = partial_derivative(psi, mesh, inverse);
    CHECK(ScalarWave::p_coarsening_constraint_change(
              psi, phi, mesh, make_array<Dim>(false), inverse) == 0.0);
    // P_7'(1)=28, and differentiation lowers the degree. Truncating Psi
    // removes its P_7 modes while the degree-six Phi survives unchanged.
    CHECK(ScalarWave::p_coarsening_constraint_change(
              psi, phi, mesh, make_array<Dim>(true), inverse) ==
          approx(28.0 * amplitude * scale * std::sqrt(Dim)));
    auto one_axis = make_array<Dim>(false);
    one_axis[0] = true;
    CHECK(ScalarWave::p_coarsening_constraint_change(psi, phi, mesh, one_axis,
                                                     inverse) ==
          approx(28.0 * amplitude * scale));
    // A pre-existing low-mode constraint does not change the injection budget.
    phi.get(0) += 0.1;
    CHECK(ScalarWave::p_coarsening_constraint_change(psi, phi, mesh, one_axis,
                                                     inverse) ==
          approx(28.0 * amplitude * scale).epsilon(1.e-9));
    const Scalar<DataVector> constant{
        DataVector(mesh.number_of_grid_points(), 1.0)};
    phi = tnsr::i<DataVector, Dim, Frame::Inertial>(
        mesh.number_of_grid_points(), 0.1);
    CHECK(ScalarWave::p_coarsening_constraint_change(
              constant, phi, mesh, make_array<Dim>(true), inverse) < 1.e-12);
    if constexpr (Dim > 1) {
      ModalVector mixed_modes(mesh.number_of_grid_points(), 0.0);
      mixed_modes[mixed_modes.size() - 1] = amplitude;
      const Scalar<DataVector> mixed{to_nodal_coefficients(mixed_modes, mesh)};
      const auto mixed_phi = partial_derivative(mixed, mesh, inverse);
      // A single-axis reduction retains one component of Phi. Simultaneous
      // reductions discard the entire mixed field and all components of Phi,
      // preserving C=0. This checks the union of intersecting modal slabs.
      CHECK(ScalarWave::p_coarsening_constraint_change(mixed, mixed_phi, mesh,
                                                       one_axis, inverse) ==
            approx(28.0 * amplitude * scale));
      CHECK(ScalarWave::p_coarsening_constraint_change(mixed, mixed_phi, mesh,
                                                       make_array<Dim>(true),
                                                       inverse) < 1.e-12);
    }
  }
}

using Guard = ScalarWave::ChalmersWithConstraintGuard<1>;
struct Metavariables {
  using component_list = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<amr::Criterion, tmpl::list<Guard>>>;
  };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.ScalarWave.ChalmersConstraintGuard",
                  "[Unit][Evolution]") {
  check_polynomial<1>();
  check_polynomial<2>();
  check_polynomial<3>();
  // A nonlinear map with dx/dxi = 1/(1+0.2 xi). Phi now has a removable
  // P_7 coefficient too: xi P_7' contributes 7 P_7. The maximum change is
  // (1.2*28 - 0.2*7) times the Psi amplitude, attained at xi=1.
  {
    const Mesh<1> mesh{8, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
    ModalVector modes(8, 0.0);
    modes[7] = 1.e-5;
    const Scalar<DataVector> psi{to_nodal_coefficients(modes, mesh)};
    InverseJacobian<DataVector, 1, Frame::ElementLogical, Frame::Inertial>
    inverse(size_t{8}, 0.0);
    get<0, 0>(inverse) = 1.0 + 0.2 * Spectral::collocation_points(mesh);
    const auto phi = partial_derivative(psi, mesh, inverse);
    CHECK(ScalarWave::p_coarsening_constraint_change(
              psi, phi, mesh, std::array{true}, inverse) == approx(32.2e-5));
  }
  register_factory_classes_with_charm<Metavariables>();
  const auto parsed = TestHelpers::test_factory_creation<amr::Criterion, Guard>(
      "ChalmersWithConstraintGuard:\n"
      "  Chalmers:\n"
      "    VariablesToMonitor: [Psi, Pi, Phi]\n"
      "    AbsoluteTarget: 1.e-2\n"
      "    RelativeTarget: 0\n"
      "    SmoothnessThreshold: 0.5\n"
      "    PRefineOnly: false\n"
      "    CoarseningFactor: 0.1\n"
      "    HCoarsening: false\n"
      "  ConstraintChangeTarget: 1.e-4\n");
  const auto& guard = dynamic_cast<const Guard&>(*parsed);
  CHECK(parsed->type() == amr::Criteria::Type::hp);
  CHECK_FALSE(guard.needs_sibling_data());
  const Mesh<1> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  ModalVector modes(8, 0.0);
  modes[7] = 1.e-5;
  const Scalar<DataVector> psi{to_nodal_coefficients(modes, mesh)};
  const InverseJacobian<DataVector, 1, Frame::ElementLogical, Frame::Inertial>
  inverse(size_t{8}, 1.0);
  const auto phi = partial_derivative(psi, mesh, inverse);
  const auto box = db::create<tmpl::list<
      domain::Tags::Mesh<1>, ScalarWave::Tags::Psi, ScalarWave::Tags::Pi,
      ScalarWave::Tags::Phi<1, Frame::Inertial>,
      domain::Tags::InverseJacobian<1, Frame::ElementLogical,
                                    Frame::Inertial>>>(mesh, psi, psi, phi,
                                                       inverse);
  Parallel::GlobalCache<Metavariables> cache{};
  const Guard::FieldCriterion fields{
      {"Psi", "Pi", "Phi"}, 1.e-2, 0.0, 0.5, false, 0.1, false};
  CHECK(fields(box, cache, ElementId<1>{0})[0] ==
        amr::Flag::DecreaseResolution);
  CHECK(guard(box, cache, ElementId<1>{0})[0] == amr::Flag::DoNothing);
  CHECK(serialize_and_deserialize(guard)(box, cache, ElementId<1>{0})[0] ==
        amr::Flag::DoNothing);
  CHECK(Guard{fields, 1.e-3}(box, cache, ElementId<1>{0})[0] ==
        amr::Flag::DecreaseResolution);
  CHECK_THROWS_WITH((Guard{fields, 0.0}),
                    Catch::Matchers::ContainsSubstring("finite and positive"));
  const Guard::FieldCriterion joins{
      {"Psi", "Pi", "Phi"}, 1.e-2, 0.0, 0.5, false, 0.1, true};
  CHECK_THROWS_WITH((Guard{joins, 1.e-3}),
                    Catch::Matchers::ContainsSubstring("HCoarsening: false"));
}
