// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"

#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gr::np {
Scalar<ComplexDataVector> coulomb_scalar(
    const Scalar<ComplexDataVector>& invariant_i,
    const Scalar<ComplexDataVector>& invariant_j) {
  return Scalar<ComplexDataVector>{-3. * get(invariant_j) / get(invariant_i)};
}

Scalar<ComplexDataVector> background_radius(
    const Scalar<ComplexDataVector>& invariant_i,
    const Scalar<ComplexDataVector>& invariant_j, const double mass) {
  const size_t num_points = get_size(get(invariant_i));
  Scalar<ComplexDataVector> radius(num_points);
  for (size_t p = 0; p < num_points; ++p) {
    get(radius)[p] = std::pow(
        mass * get(invariant_i)[p] / (3. * get(invariant_j)[p]), 1. / 3.);
  }
  return radius;
}

WeylScalars kinnersley_scalars(const Scalar<ComplexDataVector>& coulomb) {
  WeylScalars psi(get_size(get(coulomb)), std::complex<double>{0., 0.});
  psi.get(2) = get(coulomb);
  return psi;
}

WeylScalars type_d_scalars(const Scalar<ComplexDataVector>& coulomb,
                           const Scalar<ComplexDataVector>& a_bar,
                           const Scalar<ComplexDataVector>& b) {
  return type_i(type_ii(kinnersley_scalars(coulomb), b), a_bar);
}

TypeDRotation solve_type_d_rotation(const WeylScalars& psi,
                                    const Scalar<ComplexDataVector>& coulomb,
                                    const double alignment_threshold) {
  const size_t num_points = get_size(get(coulomb));
  TypeDRotation result{Scalar<ComplexDataVector>(num_points),
                       Scalar<ComplexDataVector>(num_points),
                       Scalar<ComplexDataVector>(num_points), WeylScalars{}};
  for (size_t p = 0; p < num_points; ++p) {
    const std::complex<double> ratio = psi.get(2)[p] / get(coulomb)[p];
    const std::complex<double> discriminant =
        std::sqrt(36. - 24. * (1. - ratio));
    const std::complex<double> root_plus = (-6. + discriminant) / 12.;
    const std::complex<double> root_minus = (-6. - discriminant) / 12.;
    const std::complex<double> x =
        std::abs(root_plus) <= std::abs(root_minus) ? root_plus : root_minus;
    if (std::abs(x) > 0.5) {
      ERROR(
          "The aligning root a_bar*b should be near zero; a magnitude above "
          "0.5 means the tetrad is longitudinally exchanged and the root "
          "selection is unreliable here. Got |x| = "
          << std::abs(x) << " at point " << p);
    }
    const std::complex<double> b =
        psi.get(1)[p] / (3. * (1. + 2. * x) * get(coulomb)[p]);
    if (std::abs(b) == 0. and std::abs(x) > 0.) {
      ERROR(
          "Psi1 vanishes while a_bar*b does not; the Psi1/Psi2 equation pair "
          "cannot fix the rotation at point "
          << p);
    }
    get(result.x)[p] = x;
    get(result.b)[p] = b;
    // a_bar = x / b is a ratio of two roundoff-sized numbers when the tetrad
    // is already aligned (Psi1 = Psi3 = 0, as for a hole at rest seen from
    // a radial tetrad). Below the alignment threshold, take the aligned
    // rotation instead of amplifying roundoff to O(1).
    get(result.a_bar)[p] =
        std::abs(x) <= alignment_threshold or std::abs(b) == 0.
            ? std::complex<double>{0., 0.}
            : x / b;
  }
  result.predicted_psi = type_d_scalars(coulomb, result.a_bar, result.b);
  return result;
}

WeylScalars pull_back(const WeylScalars& psi,
                      const Scalar<ComplexDataVector>& a_bar,
                      const Scalar<ComplexDataVector>& b) {
  return type_ii(type_i(psi, Scalar<ComplexDataVector>{-get(a_bar)}),
                 Scalar<ComplexDataVector>{-get(b)});
}

WeylScalars push_forward(const WeylScalars& psi,
                         const Scalar<ComplexDataVector>& a_bar,
                         const Scalar<ComplexDataVector>& b) {
  return type_i(type_ii(psi, b), a_bar);
}

Scalar<ComplexDataVector> psi0_leading(const Scalar<ComplexDataVector>& coulomb,
                                       const Scalar<ComplexDataVector>& b) {
  return Scalar<ComplexDataVector>{6. * get(b) * get(b) * get(coulomb)};
}
}  // namespace gr::np
