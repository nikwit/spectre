// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"

#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
void cholesky_factor(const gsl::not_null<RealMatrix*> lower,
                     const tnsr::ii<DataVector, 3, Frame::Inertial>& metric) {
  *lower = RealMatrix(get_size(get<0, 0>(metric)), 0.);
  auto& l00 = lower->get(0, 0);
  auto& l10 = lower->get(1, 0);
  auto& l20 = lower->get(2, 0);
  auto& l11 = lower->get(1, 1);
  auto& l21 = lower->get(2, 1);
  auto& l22 = lower->get(2, 2);
  l00 = sqrt(metric.get(0, 0));
  l10 = metric.get(1, 0) / l00;
  l20 = metric.get(2, 0) / l00;
  l11 = sqrt(metric.get(1, 1) - square(l10));
  l21 = (metric.get(2, 1) - l20 * l10) / l11;
  l22 = sqrt(metric.get(2, 2) - square(l20) - square(l21));
}

RealMatrix cholesky_factor(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& metric) {
  RealMatrix lower{};
  cholesky_factor(make_not_null(&lower), metric);
  return lower;
}

RealMatrix inverse_lower_triangular(const RealMatrix& lower) {
  RealMatrix inverse(get_size(get<0, 0>(lower)), 0.);
  const auto& l00 = lower.get(0, 0);
  const auto& l10 = lower.get(1, 0);
  const auto& l20 = lower.get(2, 0);
  const auto& l11 = lower.get(1, 1);
  const auto& l21 = lower.get(2, 1);
  const auto& l22 = lower.get(2, 2);
  inverse.get(0, 0) = 1. / l00;
  inverse.get(1, 1) = 1. / l11;
  inverse.get(2, 2) = 1. / l22;
  inverse.get(1, 0) = -l10 / (l00 * l11);
  inverse.get(2, 1) = -l21 / (l11 * l22);
  inverse.get(2, 0) = (l10 * l21 - l11 * l20) / (l00 * l11 * l22);
  return inverse;
}

RealMatrix adapted_triad(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame::Inertial>& directions) {
  const size_t num_points = get_size(get<0>(directions));
  const RealMatrix cholesky = cholesky_factor(spatial_metric);
  const RealMatrix cholesky_inverse = inverse_lower_triangular(cholesky);
  RealMatrix rotation(num_points, 0.);

  // radial = normalize(L^{-1} d)
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      rotation.get(0, i) += cholesky_inverse.get(i, j) * directions.get(j);
    }
  }
  DataVector norm =
      sqrt(square(rotation.get(0, 0)) + square(rotation.get(0, 1)) +
           square(rotation.get(0, 2)));
  for (size_t i = 0; i < 3; ++i) {
    rotation.get(0, i) /= norm;
  }

  // Coordinate-vector components of the gauge angular directions, with the
  // azimuth fixed to zero at a coordinate pole
  const DataVector cylindrical_radius =
      sqrt(square(directions.get(0)) + square(directions.get(1)));
  DataVector cos_phi(num_points, 1.);
  DataVector sin_phi(num_points, 0.);
  for (size_t p = 0; p < num_points; ++p) {
    if (cylindrical_radius[p] > 1.e-14) {
      cos_phi[p] = directions.get(0)[p] / cylindrical_radius[p];
      sin_phi[p] = directions.get(1)[p] / cylindrical_radius[p];
    }
  }
  const std::array<DataVector, 3> theta_coordinate{{directions.get(2) * cos_phi,
                                                    directions.get(2) * sin_phi,
                                                    -cylindrical_radius}};
  const std::array<DataVector, 3> phi_coordinate{
      {-sin_phi, cos_phi, DataVector(num_points, 0.)}};

  // Carry to the orthonormal triad with L^T: v_i = sum_j L_ji v^j_coord
  std::array<DataVector, 3> theta{};
  std::array<DataVector, 3> phi{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(theta, i) = DataVector(num_points, 0.);
    gsl::at(phi, i) = DataVector(num_points, 0.);
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(theta, i) += cholesky.get(j, i) * gsl::at(theta_coordinate, j);
      gsl::at(phi, i) += cholesky.get(j, i) * gsl::at(phi_coordinate, j);
    }
  }

  const auto dot = [](const auto& a, const auto& b) {
    return DataVector{a[0] * b[0] + a[1] * b[1] + a[2] * b[2]};
  };
  const std::array<DataVector, 3> radial{
      {rotation.get(0, 0), rotation.get(0, 1), rotation.get(0, 2)}};

  // Gram-Schmidt: theta against radial, phi against radial and theta
  DataVector projection = dot(radial, theta);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(theta, i) -= gsl::at(radial, i) * projection;
  }
  norm = sqrt(dot(theta, theta));
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(theta, i) /= norm;
  }
  projection = dot(radial, phi);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(phi, i) -= gsl::at(radial, i) * projection;
  }
  projection = dot(theta, phi);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(phi, i) -= gsl::at(theta, i) * projection;
  }
  norm = sqrt(dot(phi, phi));
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(phi, i) /= norm;
  }

  // Orientation s x theta = phi
  const std::array<DataVector, 3> cross{
      {radial[1] * theta[2] - radial[2] * theta[1],
       radial[2] * theta[0] - radial[0] * theta[2],
       radial[0] * theta[1] - radial[1] * theta[0]}};
  const DataVector orientation = dot(cross, phi);
  for (size_t p = 0; p < num_points; ++p) {
    if (orientation[p] < 0.) {
      for (size_t i = 0; i < 3; ++i) {
        gsl::at(phi, i)[p] *= -1.;
      }
    }
  }

  for (size_t i = 0; i < 3; ++i) {
    rotation.get(1, i) = gsl::at(theta, i);
    rotation.get(2, i) = gsl::at(phi, i);
  }
  return rotation;
}

ComplexMatrix rotate_symmetric(const ComplexMatrix& tensor,
                               const RealMatrix& rotation) {
  ComplexMatrix result(get_size(get<0, 0>(tensor)),
                       std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          result.get(i, j) +=
              rotation.get(i, k) * tensor.get(k, l) * rotation.get(j, l);
        }
      }
    }
  }
  return result;
}

WeylScalars weyl_scalars_from_tidal_tensor(const ComplexMatrix& q) {
  const std::complex<double> imaginary_unit{0., 1.};
  WeylScalars psi(get_size(get<0, 0>(q)), std::complex<double>{0., 0.});
  const ComplexDataVector transverse_difference =
      0.5 * (q.get(1, 1) - q.get(2, 2));
  psi.get(0) = transverse_difference + imaginary_unit * q.get(1, 2);
  psi.get(1) = -0.5 * (q.get(0, 1) + imaginary_unit * q.get(0, 2));
  psi.get(2) = 0.5 * q.get(0, 0);
  psi.get(3) = 0.5 * (q.get(0, 1) - imaginary_unit * q.get(0, 2));
  psi.get(4) = transverse_difference - imaginary_unit * q.get(1, 2);
  return psi;
}

ComplexMatrix tidal_tensor_from_weyl_scalars(const WeylScalars& psi) {
  const std::complex<double> imaginary_unit{0., 1.};
  ComplexMatrix q(get_size(get<0>(psi)), std::complex<double>{0., 0.});
  const auto& p0 = psi.get(0);
  const auto& p1 = psi.get(1);
  const auto& p2 = psi.get(2);
  const auto& p3 = psi.get(3);
  const auto& p4 = psi.get(4);
  q.get(0, 0) = 2. * p2;
  q.get(0, 1) = p3 - p1;
  q.get(1, 0) = q.get(0, 1);
  q.get(0, 2) = imaginary_unit * (p1 + p3);
  q.get(2, 0) = q.get(0, 2);
  q.get(1, 1) = -p2 + 0.5 * (p0 + p4);
  q.get(1, 2) = 0.5 * imaginary_unit * (p4 - p0);
  q.get(2, 1) = q.get(1, 2);
  q.get(2, 2) = -p2 - 0.5 * (p0 + p4);
  return q;
}

WeylScalars weyl_scalars_from_electric_magnetic(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame::Inertial>& directions) {
  const size_t num_points = get_size(get<0>(directions));
  const std::complex<double> imaginary_unit{0., 1.};
  const RealMatrix cholesky_inverse =
      inverse_lower_triangular(cholesky_factor(spatial_metric));

  // Q_raw = L^{-1} (E + iB) L^{-T}
  ComplexMatrix q_raw(num_points, std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          q_raw.get(i, j) +=
              cholesky_inverse.get(i, k) *
              (electric.get(k, l) + imaginary_unit * magnetic.get(k, l)) *
              cholesky_inverse.get(j, l);
        }
      }
    }
  }
  // Hygiene projection: symmetrize and remove the trace
  ComplexMatrix q(num_points, std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      q.get(i, j) = 0.5 * (q_raw.get(i, j) + q_raw.get(j, i));
    }
  }
  const ComplexDataVector trace = q.get(0, 0) + q.get(1, 1) + q.get(2, 2);
  for (size_t i = 0; i < 3; ++i) {
    q.get(i, i) -= trace / 3.;
  }
  return weyl_scalars_from_tidal_tensor(
      rotate_symmetric(q, adapted_triad(spatial_metric, directions)));
}

tnsr::ii<DataVector, 3, Frame::Inertial> incoming_weyl_field(
    const Scalar<ComplexDataVector>& psi0, const RealMatrix& rotation) {
  const std::complex<double> imaginary_unit{0., 1.};
  const double one_over_sqrt2 = 1. / sqrt(2.);
  tnsr::ii<DataVector, 3, Frame::Inertial> result(get_size(get(psi0)), 0.);
  std::array<ComplexDataVector, 3> m{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(m, i) = one_over_sqrt2 *
                    (rotation.get(1, i) + imaginary_unit * rotation.get(2, i));
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      const ComplexDataVector mm = gsl::at(m, i) * gsl::at(m, j);
      result.get(i, j) =
          real(2. * (conj(get(psi0)) * mm + get(psi0) * conj(mm)));
    }
  }
  return result;
}

tnsr::ii<DataVector, 3, Frame::Inertial> orthonormal_to_coordinate_covariant(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& orthonormal,
    const RealMatrix& cholesky) {
  tnsr::ii<DataVector, 3, Frame::Inertial> result(
      get_size(get<0, 0>(orthonormal)), 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      for (size_t a = 0; a < 3; ++a) {
        for (size_t b = 0; b < 3; ++b) {
          result.get(i, j) +=
              cholesky.get(i, a) * cholesky.get(j, b) * orthonormal.get(a, b);
        }
      }
    }
  }
  return result;
}
}  // namespace gr::np
