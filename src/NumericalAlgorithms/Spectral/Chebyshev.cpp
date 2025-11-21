// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/BasisFunctionNormalizationSquare.hpp"
#include "NumericalAlgorithms/Spectral/BasisFunctionValue.hpp"
#include "NumericalAlgorithms/Spectral/Chebyshev.hpp"
#include "NumericalAlgorithms/Spectral/Clenshaw.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPointsAndWeights.hpp"
#include "NumericalAlgorithms/Spectral/InverseWeightFunctionValues.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace Spectral {

// Algorithms to compute Chebyshev basis functions
// These functions specialize the templates declared in `Spectral.hpp`.

namespace {
template <typename T>
T compute_basis_function_value_impl(const size_t k, const T& x) {
  // Algorithm 21 in Kopriva, p. 60
  switch (k) {
    case 0:
      return make_with_value<T>(x, 1.);
    case 1:
      return x;
    default:
      // These values can be computed either through recursion
      // (implemented here), or analytically as `cos(k * acos(x))`.
      // Since the trigonometric form is expensive to compute it is useful only
      // for large k. See Kopriva, section 3.1 (p. 59) and Fig. 3.1 (p. 61) for
      // a discussion.
      T T_k_minus_2 = make_with_value<T>(x, 1.);
      T T_k_minus_1 = x;
      T T_k = make_with_value<T>(x, 0.);
      for (size_t j = 2; j <= k; j++) {
        T_k = 2 * x * T_k_minus_1 - T_k_minus_2;
        T_k_minus_2 = T_k_minus_1;
        T_k_minus_1 = T_k;
      }
      return T_k;
  }
}
}  // namespace

template <>
DataVector compute_basis_function_value<Basis::Chebyshev>(const size_t k,
                                                          const DataVector& x) {
  return compute_basis_function_value_impl(k, x);
}

template <>
double compute_basis_function_value<Basis::Chebyshev>(const size_t k,
                                                      const double& x) {
  return compute_basis_function_value_impl(k, x);
}

template <>
DataVector compute_inverse_weight_function_values<Basis::Chebyshev>(
    const DataVector& x) {
  return sqrt(1. - square(x));
}

template <>
double compute_basis_function_normalization_square<Basis::Chebyshev>(
    const size_t k) {
  if (k == 0) {
    return M_PI;
  } else {
    return M_PI_2;
  }
}

// Algorithm to compute Chebyshev-Gauss quadrature

template <>
std::pair<DataVector, DataVector>
compute_collocation_points_and_weights<Basis::Chebyshev, Quadrature::Gauss>(
    const size_t num_points) {
  // Algorithm 26 in Kopriva, p. 67
  ASSERT(num_points >= 1,
         "Chebyshev-Gauss quadrature requires at least one collocation point.");
  const size_t poly_degree = num_points - 1;
  DataVector x(num_points);
  DataVector w(num_points, M_PI / num_points);
  for (size_t j = 0; j < num_points; j++) {
    x[j] = -cos(M_PI_2 * (2. * j + 1.) / (poly_degree + 1.));
  }
  return std::make_pair(std::move(x), std::move(w));
}

// Algorithm to compute Chebyshev-Gauss-Lobatto quadrature

template <>
std::pair<DataVector, DataVector> compute_collocation_points_and_weights<
    Basis::Chebyshev, Quadrature::GaussLobatto>(const size_t num_points) {
  // Algorithm 27 in Kopriva, p. 68
  ASSERT(num_points >= 2,
         "Chebyshev-Gauss-Lobatto quadrature requires at least two collocation "
         "points.");
  const size_t poly_degree = num_points - 1;
  DataVector x(num_points);
  DataVector w(num_points, M_PI / poly_degree);
  for (size_t j = 0; j < num_points; j++) {
    x[j] = -cos(M_PI * j / poly_degree);
  }
  w[0] *= 0.5;
  w[num_points - 1] *= 0.5;
  return std::make_pair(std::move(x), std::move(w));
}

template <Basis BasisType>
Matrix spectral_indefinite_integral_matrix(size_t num_points);

template <>
Matrix spectral_indefinite_integral_matrix<Basis::Chebyshev>(
    const size_t num_points) {
  // Tridiagonal matrix that gives the indefinite integral modulo a constant
  Matrix indef_int(num_points, num_points, 0.0);
  if (LIKELY(num_points > 1)) {
    indef_int(1, 0) = 1.0;
  }
  if (LIKELY(num_points > 2)) {
    indef_int(1, 2) = -0.5;
    indef_int(num_points - 1, num_points - 2) =
        1.0 / (2.0 * (num_points - 1.0));
  }
  for (size_t i = 2; i < num_points - 1; ++i) {
    indef_int(i, i - 1) = 1.0 / (2.0 * i);
    indef_int(i, i + 1) = -1.0 / (2.0 * i);
  }

  // Matrix that ensures that BC at left of interval is 0.0
  Matrix constant(num_points, num_points, 0.0);
  double fac = 1.0;
  for (size_t i = 1; i < num_points; ++i) {
    constant(i, i) = 1.0;
    constant(0, i) = fac;
    fac = -fac;
  }
  return constant * indef_int;
}

namespace {
double evaluate_chebyshev_segment(gsl::span<const double> coefficients,
                                  const double x) {
  if (coefficients.empty()) {
    return 0.0;
  }
  double y_upper = 0.0;
  double y_lower = 0.0;
  for (size_t k = coefficients.size(); k > 1; --k) {
    const double coefficient = gsl::at(coefficients, k - 1);
    const double new_y_lower = coefficient + 2.0 * x * y_lower - y_upper;
    y_upper = y_lower;
    y_lower = new_y_lower;
  }
  return gsl::at(coefficients, 0) + x * y_lower - y_upper;
}
}  // namespace

template <size_t Dim>
double evaluate_chebyshev_series(
    const ModalVector& coefficients, const Mesh<Dim>& mesh,
    const tnsr::I<double, Dim, Frame::ElementLogical>& logical_coords) {
  ASSERT(mesh.number_of_grid_points() == coefficients.size(),
         "Mesh and coefficient sizes do not match: mesh has "
             << mesh.number_of_grid_points() << " points but coefficients hold "
             << coefficients.size() << " entries.");
  const auto& extents = mesh.extents();
  for (size_t d = 0; d < Dim; ++d) {
    ASSERT(mesh.basis(d) == Basis::Chebyshev,
           "evaluate_chebyshev_series only supports Chebyshev bases. Found "
               << mesh.basis(d) << " in dimension " << d << ".");
  }
  std::vector<double> working(coefficients.begin(), coefficients.end());
  size_t current_size = working.size();
  for (size_t dim = 0; dim < Dim; ++dim) {
    const size_t extent = extents[dim];
    ASSERT(extent > 0, "Mesh extent must be non-zero.");
    const size_t num_slices = current_size / extent;
    const double x_dim = logical_coords.get(dim);
    for (size_t slice = 0; slice < num_slices; ++slice) {
      const gsl::span<const double> slice_view(working.data() + slice * extent,
                                               extent);
      working[slice] = evaluate_chebyshev_segment(slice_view, x_dim);
    }
    current_size = num_slices;
  }
  ASSERT(current_size == 1,
         "After evaluating all dimensions there should be a single value.");
  return working.front();
}

}  // namespace Spectral

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                        \
  template double Spectral::evaluate_chebyshev_series<DIM(data)>(   \
      const ModalVector& coefficients, const Mesh<DIM(data)>& mesh, \
      const tnsr::I<double, DIM(data), Frame::ElementLogical>&      \
          logical_coords);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
