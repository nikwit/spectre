// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/LinearSolver/Lapack.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
std::array<double, 5> DirectPsi4Fit::electric() const {
  std::array<double, 5> result{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(result, a) = std::real(gsl::at(components, a));
  }
  return result;
}

std::array<double, 5> DirectPsi4Fit::magnetic() const {
  std::array<double, 5> result{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(result, a) = std::imag(gsl::at(components, a));
  }
  return result;
}

DirectPsi4Fit fit_psi4(
    const Scalar<ComplexDataVector>& measured_psi4,
    const std::array<Scalar<ComplexDataVector>, 5>& projected_columns,
    const std::optional<DataVector>& point_weights) {
  const size_t num_points = get_size(get(measured_psi4));
  DataVector root_weights(num_points, 1.);
  if (point_weights.has_value()) {
    root_weights = sqrt(*point_weights);
  }
  // Weighted data and design, columns normalized as a preconditioner
  const ComplexDataVector data = get(measured_psi4) * root_weights;
  std::array<ComplexDataVector, 5> design{};
  std::array<double, 5> norms{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(design, a) = get(gsl::at(projected_columns, a)) * root_weights;
    double norm_squared = 0.;
    for (size_t p = 0; p < num_points; ++p) {
      norm_squared += std::norm(gsl::at(design, a)[p]);
    }
    gsl::at(norms, a) = sqrt(norm_squared);
    if (gsl::at(norms, a) <= 0.) {
      ERROR("Psi4 design contains an empty component column " << a);
    }
  }
  // Normal equations (Z^H Z) y = Z^H d of the normalized design, written as
  // the equivalent real 10x10 system for (Re y, Im y) and solved with LAPACK
  Matrix normal_matrix(10, 10, 0.);
  DataVector normal_rhs(10, 0.);
  for (size_t a = 0; a < 5; ++a) {
    for (size_t b = 0; b < 5; ++b) {
      std::complex<double> sum{0., 0.};
      for (size_t p = 0; p < num_points; ++p) {
        sum += std::conj(gsl::at(design, a)[p]) * gsl::at(design, b)[p];
      }
      sum /= gsl::at(norms, a) * gsl::at(norms, b);
      normal_matrix(a, b) = std::real(sum);
      normal_matrix(a, b + 5) = -std::imag(sum);
      normal_matrix(a + 5, b) = std::imag(sum);
      normal_matrix(a + 5, b + 5) = std::real(sum);
    }
    std::complex<double> sum{0., 0.};
    for (size_t p = 0; p < num_points; ++p) {
      sum += std::conj(gsl::at(design, a)[p]) * data[p];
    }
    sum /= gsl::at(norms, a);
    normal_rhs[a] = std::real(sum);
    normal_rhs[a + 5] = std::imag(sum);
  }
  DataVector solution(10, 0.);
  const int info = lapack::general_matrix_linear_solve(
      make_not_null(&solution), normal_matrix, normal_rhs);
  if (info != 0) {
    ERROR("The Psi4 normal equations could not be solved; LAPACK returned "
          << info);
  }
  DirectPsi4Fit fit{};
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(fit.components, a) =
        std::complex<double>{solution[a], solution[a + 5]} / gsl::at(norms, a);
  }
  double residual_squared = 0.;
  double data_squared = 0.;
  for (size_t p = 0; p < num_points; ++p) {
    std::complex<double> fitted{0., 0.};
    for (size_t a = 0; a < 5; ++a) {
      fitted += gsl::at(design, a)[p] * gsl::at(fit.components, a);
    }
    residual_squared += std::norm(fitted - data[p]);
    data_squared += std::norm(data[p]);
  }
  fit.relative_residual =
      sqrt(residual_squared) / std::max(sqrt(data_squared), 1.e-300);
  return fit;
}

double relative_fit_residual(
    const Scalar<ComplexDataVector>& measured_psi4,
    const std::array<Scalar<ComplexDataVector>, 5>& projected_columns,
    const TidalMoments& components,
    const std::optional<DataVector>& point_weights) {
  const size_t num_points = get_size(get(measured_psi4));
  double residual_squared = 0.;
  double data_squared = 0.;
  for (size_t p = 0; p < num_points; ++p) {
    const double weight = point_weights.has_value() ? (*point_weights)[p] : 1.;
    std::complex<double> fitted{0., 0.};
    for (size_t a = 0; a < 5; ++a) {
      fitted += get(gsl::at(projected_columns, a))[p] * gsl::at(components, a);
    }
    residual_squared += weight * std::norm(fitted - get(measured_psi4)[p]);
    data_squared += weight * std::norm(get(measured_psi4)[p]);
  }
  return sqrt(residual_squared) / std::max(sqrt(data_squared), 1.e-300);
}

FrameRegistration register_frame(const WeylScalars& psi,
                                 const RealMatrix& adapted_rotation,
                                 const double mass) {
  const size_t num_points = get_size(psi.get(0));
  FrameRegistration registration{};
  registration.coulomb = coulomb_scalar(invariant_i(psi), invariant_j(psi));
  registration.rotation = solve_type_d_rotation(psi, registration.coulomb);
  registration.pulled_back =
      pull_back(psi, registration.rotation.a_bar, registration.rotation.b);
  registration.member = tangent_boost_member(registration.rotation.a_bar,
                                             registration.rotation.b);
  // r = (-M / Psi2^K)^{1/3} on the principal branch, real part
  registration.measured_radius = Scalar<DataVector>(num_points, 0.);
  for (size_t p = 0; p < num_points; ++p) {
    get(registration.measured_radius)[p] =
        std::real(std::pow(-mass / get(registration.coulomb)[p], 1. / 3.));
  }
  // Cholesky-triad components: v_i = sum_j R_ji v_j
  registration.radial_direction = TriadVector(num_points, 0.);
  registration.transverse_velocity = TriadVector(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      registration.radial_direction.get(i) +=
          adapted_rotation.get(j, i) *
          registration.member.radial_direction.get(j);
      registration.transverse_velocity.get(i) +=
          adapted_rotation.get(j, i) *
          registration.member.transverse_velocity.get(j);
    }
  }
  return registration;
}

SecondOrderEvaluation evaluate_second_order(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& adapted_rotation, const double mass,
    const std::optional<DataVector>& fit_point_weights,
    const std::optional<TidalMoments>& imposed_components) {
  SecondOrderEvaluation evaluation{};
  evaluation.direct_columns = direct_tide_scalar_columns(
      registration.radial_direction, registration.transverse_velocity, rapidity,
      registration.measured_radius, adapted_rotation, mass);
  std::array<Scalar<ComplexDataVector>, 5> projected_psi4{};
  for (size_t a = 0; a < 5; ++a) {
    get(gsl::at(projected_psi4, a)) =
        pull_back(gsl::at(evaluation.direct_columns, a),
                  registration.rotation.a_bar, registration.rotation.b)
            .get(4);
  }
  const Scalar<ComplexDataVector> pulled_back_psi4{
      registration.pulled_back.get(4)};
  if (imposed_components.has_value()) {
    evaluation.fit.components = *imposed_components;
    evaluation.fit.relative_residual =
        relative_fit_residual(pulled_back_psi4, projected_psi4,
                              *imposed_components, fit_point_weights);
  } else {
    evaluation.fit =
        fit_psi4(pulled_back_psi4, projected_psi4, fit_point_weights);
  }
  // Psi0 target: Kinnersley scalars pushed forward to the NR tetrad, plus the
  // fitted transverse tide in the NR tetrad
  const WeylScalars kinematic_nr =
      push_forward(kinnersley_scalars(registration.coulomb),
                   registration.rotation.a_bar, registration.rotation.b);
  get(evaluation.psi0_target) = kinematic_nr.get(0);
  for (size_t a = 0; a < 5; ++a) {
    get(evaluation.psi0_target) += gsl::at(evaluation.fit.components, a) *
                                   gsl::at(evaluation.direct_columns, a).get(0);
  }
  return evaluation;
}
}  // namespace gr::np
