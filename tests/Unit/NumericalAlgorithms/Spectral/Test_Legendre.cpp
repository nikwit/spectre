// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/BasisFunctionValue.hpp"
#include "NumericalAlgorithms/Spectral/Legendre.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"

namespace {

template <size_t Dim>
tnsr::I<double, Dim, Frame::ElementLogical> make_logical_point(
    const std::array<double, Dim>& coords) {
  tnsr::I<double, Dim, Frame::ElementLogical> point{};
  for (size_t d = 0; d < Dim; ++d) {
    point.get(d) = coords[d];
  }
  return point;
}

double evaluate_by_summing_legendre_basis(const ModalVector& coefficients,
                                          const double x) {
  double sum = 0.0;
  for (size_t k = 0; k < coefficients.size(); ++k) {
    sum += coefficients[k] *
           Spectral::compute_basis_function_value<Spectral::Basis::Legendre>(
               k, x);
  }
  return sum;
}

template <size_t Dim>
void check_against_irregular(const std::array<size_t, Dim>& extents_array,
                             const std::array<double, Dim>& coords) {
  Mesh<Dim> mesh(extents_array, Spectral::Basis::Legendre,
                 Spectral::Quadrature::GaussLobatto);
  ModalVector modal_coefficients(mesh.number_of_grid_points());
  for (size_t i = 0; i < modal_coefficients.size(); ++i) {
    modal_coefficients[i] =
        -0.4 + 0.1 * static_cast<double>(i) -
        0.03 * static_cast<double>((i + 3) % 5);
  }

  const DataVector nodal_coefficients =
      to_nodal_coefficients(modal_coefficients, mesh);
  const auto logical_point = make_logical_point(coords);
  const intrp::Irregular<Dim> interpolant(mesh, logical_point);
  const DataVector nodal_value = interpolant.interpolate(nodal_coefficients);

  CHECK(Spectral::evaluate_legendre_series<Dim>(modal_coefficients, mesh,
                                                logical_point) ==
        approx(nodal_value[0]));
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Numerical.Spectral.Legendre.Series",
                  "[NumericalAlgorithms][Spectral][Unit]") {
  const std::array<double, 5> evaluation_points{{-0.95, -0.31, 0.0, 0.42, 0.99}};

  ModalVector cubic_coefficients{4};
  cubic_coefficients[0] = 1.1;
  cubic_coefficients[1] = -0.4;
  cubic_coefficients[2] = 0.25;
  cubic_coefficients[3] = -1.5;
  Mesh<1> cubic_mesh({{cubic_coefficients.size()}}, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto);
  for (const double x : evaluation_points) {
    CHECK(Spectral::evaluate_legendre_series<1>(
              cubic_coefficients, cubic_mesh, make_logical_point<1>({{x}})) ==
          approx(evaluate_by_summing_legendre_basis(cubic_coefficients, x)));
  }

  ModalVector higher_order{7};
  higher_order[0] = 0.8;
  higher_order[1] = -1.4;
  higher_order[2] = 0.2;
  higher_order[3] = 0.9;
  higher_order[4] = -0.1;
  higher_order[5] = 0.33;
  higher_order[6] = -2.1;
  Mesh<1> higher_mesh({{higher_order.size()}}, Spectral::Basis::Legendre,
                      Spectral::Quadrature::GaussLobatto);
  for (const double x : evaluation_points) {
    CHECK(Spectral::evaluate_legendre_series<1>(
              higher_order, higher_mesh, make_logical_point<1>({{x}})) ==
          approx(evaluate_by_summing_legendre_basis(higher_order, x)));
  }

  check_against_irregular<1>({{5}}, {{-0.21}});
  check_against_irregular<2>({{4, 3}}, {{0.2, -0.6}});
  check_against_irregular<3>({{3, 2, 4}}, {{0.4, -0.5, 0.1}});
}

