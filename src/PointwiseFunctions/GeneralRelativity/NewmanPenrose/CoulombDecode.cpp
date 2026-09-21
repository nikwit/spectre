// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/CoulombDecode.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numeric>
#include <optional>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/LinearSolver/Lapack.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
DataVector background_radial_derivative(const DataVector& radius,
                                        const double mass) {
  return sqrt(1. - 2. * mass / radius) * 3. * mass / pow<4>(radius);
}

DataVector normal_derivative_factor(const TangentBoostMember& member,
                                    const Scalar<DataVector>& rapidity) {
  return cosh(get(rapidity)) * member.radial_direction.get(0) +
         sinh(get(rapidity)) * get(member.lorentz_factor) *
             member.transverse_velocity.get(0);
}

RadiusSolve radius_from_normal_derivative(
    const DataVector& normal_derivative_of_coulomb_real,
    const DataVector& factor, const DataVector& initial_radius,
    const double mass, const size_t iterations) {
  const DataVector magnitude =
      abs(normal_derivative_of_coulomb_real) / abs(factor);
  const double floor = (coulomb_decode_turning_point_over_mass + 0.01) * mass;
  RadiusSolve solve{};
  solve.radius = initial_radius;
  for (auto& r : solve.radius) {
    r = std::max(r, floor);
  }
  for (size_t iteration = 0; iteration < iterations; ++iteration) {
    for (size_t p = 0; p < solve.radius.size(); ++p) {
      double& r = solve.radius[p];
      const double f = 1. - 2. * mass / r;
      const double b = sqrt(f) * 3. * mass / pow<4>(r);
      const double db = b * ((mass / r) / f - 4.) / r;
      const double step =
          std::clamp((b - magnitude[p]) / db, -0.2 * r, 0.2 * r);
      r = std::max(r - step, floor);
    }
  }
  const DataVector residual =
      abs(background_radial_derivative(solve.radius, mass) - magnitude) /
      magnitude;
  solve.newton_residual = max(residual);
  solve.valid = min(solve.radius) > floor * (1. + 1.e-12) and
                solve.newton_residual < 1.e-8;
  return solve;
}

std::array<Scalar<ComplexDataVector>, 10> coulomb_tide_columns(
    const TriadVector& radial_direction, const Scalar<DataVector>& radius,
    const double mass) {
  const size_t num_points = get_size(get(radius));
  const tnsr::ii<double, 3, Frame::Inertial> zero{0.};
  std::array<Scalar<ComplexDataVector>, 10> columns{};
  for (size_t column = 0; column < 10; ++column) {
    std::array<double, 5> unit{};
    gsl::at(unit, column % 5) = 1.;
    const auto direct = stf_from_components(unit);
    const ComplexMatrix tide_rest =
        column < 5 ? quadrupole_tide_tensor(direct, zero, radial_direction,
                                            radius, mass, false)
                   : quadrupole_tide_tensor(zero, direct, radial_direction,
                                            radius, mass, false);
    get(gsl::at(columns, column)) =
        ComplexDataVector(num_points, std::complex<double>{0., 0.});
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        get(gsl::at(columns, column)) += 0.5 * radial_direction.get(i) *
                                         tide_rest.get(i, j) *
                                         radial_direction.get(j);
      }
    }
  }
  return columns;
}

CoulombDecode decode_tidal_moments_from_coulomb(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const double mass,
    const Scalar<DataVector>& normal_derivative_of_coulomb_real,
    const std::optional<DataVector>& point_weights) {
  const size_t num_points = get_size(get(registration.coulomb));
  CoulombDecode result{};
  const DataVector factor =
      normal_derivative_factor(registration.member, rapidity);
  const RadiusSolve solve = radius_from_normal_derivative(
      get(normal_derivative_of_coulomb_real), factor,
      get(registration.measured_radius), mass);
  result.valid = solve.valid;
  result.newton_residual = solve.newton_residual;
  get(result.areal_radius) = solve.radius;
  const ComplexDataVector excess =
      get(registration.coulomb) + mass / cube(solve.radius);

  const auto columns = coulomb_tide_columns(registration.radial_direction,
                                            result.areal_radius, mass);

  // Real least squares for (E_a, B_a, Re/Im c_0, Re/Im c_1): the real and
  // imaginary parts of the excess are the rows, weighted, with the columns
  // normalized as a preconditioner and the normal equations solved with
  // LAPACK
  constexpr size_t n_unknown = 18;
  const size_t n_row = 2 * num_points;
  std::array<DataVector, n_unknown> design{};
  for (auto& column : design) {
    column = DataVector(n_row, 0.);
  }
  DataVector target(n_row, 0.);
  DataVector root_weight(n_row, 1.);
  for (size_t p = 0; p < num_points; ++p) {
    const double w = point_weights.has_value() ? sqrt((*point_weights)[p]) : 1.;
    root_weight[p] = w;
    root_weight[p + num_points] = w;
    for (size_t a = 0; a < 10; ++a) {
      gsl::at(design, a)[p] = std::real(get(gsl::at(columns, a))[p]);
      gsl::at(design, a)[p + num_points] =
          std::imag(get(gsl::at(columns, a))[p]);
    }
    gsl::at(design, 10)[p] = 1.;
    gsl::at(design, 11)[p + num_points] = 1.;
    for (size_t k = 0; k < 3; ++k) {
      gsl::at(design, 12 + k)[p] = registration.radial_direction.get(k)[p];
      gsl::at(design, 15 + k)[p + num_points] =
          registration.radial_direction.get(k)[p];
    }
    target[p] = std::real(excess[p]);
    target[p + num_points] = std::imag(excess[p]);
  }
  for (auto& column : design) {
    column *= root_weight;
  }
  target *= root_weight;
  std::array<double, n_unknown> norms{};
  for (size_t a = 0; a < n_unknown; ++a) {
    gsl::at(norms, a) = sqrt(
        std::inner_product(gsl::at(design, a).begin(), gsl::at(design, a).end(),
                           gsl::at(design, a).begin(), 0.));
    if (gsl::at(norms, a) <= 0.) {
      ERROR("The Coulomb decode design has an empty column " << a);
    }
  }
  Matrix normal_matrix(n_unknown, n_unknown, 0.);
  DataVector normal_rhs(n_unknown, 0.);
  for (size_t a = 0; a < n_unknown; ++a) {
    for (size_t b = 0; b < n_unknown; ++b) {
      normal_matrix(a, b) = std::inner_product(gsl::at(design, a).begin(),
                                               gsl::at(design, a).end(),
                                               gsl::at(design, b).begin(), 0.) /
                            (gsl::at(norms, a) * gsl::at(norms, b));
    }
    normal_rhs[a] =
        std::inner_product(gsl::at(design, a).begin(), gsl::at(design, a).end(),
                           target.begin(), 0.) /
        gsl::at(norms, a);
  }
  DataVector solution(n_unknown, 0.);
  const int info = lapack::general_matrix_linear_solve(
      make_not_null(&solution), normal_matrix, normal_rhs);
  if (info != 0) {
    ERROR(
        "The Coulomb decode normal equations could not be solved; LAPACK "
        "returned "
        << info);
  }
  for (size_t a = 0; a < n_unknown; ++a) {
    solution[a] /= gsl::at(norms, a);
  }
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(result.components, a) =
        std::complex<double>{solution[a], solution[a + 5]};
  }
  result.nuisance[0] = std::complex<double>{solution[10], solution[11]};
  for (size_t k = 0; k < 3; ++k) {
    gsl::at(result.nuisance, 1 + k) =
        std::complex<double>{solution[12 + k], solution[15 + k]};
  }
  double residual_squared = 0.;
  double target_squared = 0.;
  for (size_t row = 0; row < n_row; ++row) {
    double fitted = 0.;
    for (size_t a = 0; a < n_unknown; ++a) {
      fitted += gsl::at(design, a)[row] * solution[a];
    }
    residual_squared += square(fitted - target[row]);
    target_squared += square(target[row]);
  }
  result.relative_residual =
      sqrt(residual_squared) / std::max(sqrt(target_squared), 1.e-300);
  return result;
}
}  // namespace gr::np
