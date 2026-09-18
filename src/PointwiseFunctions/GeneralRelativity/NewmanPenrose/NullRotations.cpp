// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"

#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
WeylScalars type_i(const WeylScalars& psi,
                   const Scalar<ComplexDataVector>& a_bar) {
  const auto& a = get(a_bar);
  const ComplexDataVector a2 = a * a;
  const ComplexDataVector a3 = a2 * a;
  const ComplexDataVector a4 = a3 * a;
  const auto& p0 = psi.get(0);
  const auto& p1 = psi.get(1);
  const auto& p2 = psi.get(2);
  const auto& p3 = psi.get(3);
  const auto& p4 = psi.get(4);
  WeylScalars out(get_size(p0), std::complex<double>{0., 0.});
  out.get(0) = p0;
  out.get(1) = p1 + a * p0;
  out.get(2) = p2 + 2. * a * p1 + a2 * p0;
  out.get(3) = p3 + 3. * a * p2 + 3. * a2 * p1 + a3 * p0;
  out.get(4) = p4 + 4. * a * p3 + 6. * a2 * p2 + 4. * a3 * p1 + a4 * p0;
  return out;
}

WeylScalars type_ii(const WeylScalars& psi,
                    const Scalar<ComplexDataVector>& b) {
  const auto& b1 = get(b);
  const ComplexDataVector b2 = b1 * b1;
  const ComplexDataVector b3 = b2 * b1;
  const ComplexDataVector b4 = b3 * b1;
  const auto& p0 = psi.get(0);
  const auto& p1 = psi.get(1);
  const auto& p2 = psi.get(2);
  const auto& p3 = psi.get(3);
  const auto& p4 = psi.get(4);
  WeylScalars out(get_size(p0), std::complex<double>{0., 0.});
  out.get(4) = p4;
  out.get(3) = p3 + b1 * p4;
  out.get(2) = p2 + 2. * b1 * p3 + b2 * p4;
  out.get(1) = p1 + 3. * b1 * p2 + 3. * b2 * p3 + b3 * p4;
  out.get(0) = p0 + 4. * b1 * p1 + 6. * b2 * p2 + 4. * b3 * p3 + b4 * p4;
  return out;
}

WeylScalars type_iii(const WeylScalars& psi, const Scalar<DataVector>& eta,
                     const Scalar<DataVector>& chi) {
  const size_t num_points = get_size(get(eta));
  WeylScalars out(num_points, std::complex<double>{0., 0.});
  for (size_t a = 0; a < 5; ++a) {
    const double weight = 2. - static_cast<double>(a);
    for (size_t p = 0; p < num_points; ++p) {
      out.get(a)[p] =
          psi.get(a)[p] *
          std::exp(weight * std::complex<double>{get(eta)[p], get(chi)[p]});
    }
  }
  return out;
}

void invariants(const gsl::not_null<Scalar<ComplexDataVector>*> invariant_i,
                const gsl::not_null<Scalar<ComplexDataVector>*> invariant_j,
                const WeylScalars& psi) {
  const auto& p0 = psi.get(0);
  const auto& p1 = psi.get(1);
  const auto& p2 = psi.get(2);
  const auto& p3 = psi.get(3);
  const auto& p4 = psi.get(4);
  get(*invariant_i) = p0 * p4 - 4. * p1 * p3 + 3. * p2 * p2;
  get(*invariant_j) = p0 * p2 * p4 - p4 * p1 * p1 - p0 * p3 * p3 +
                      2. * p1 * p2 * p3 - p2 * p2 * p2;
}

Scalar<ComplexDataVector> invariant_i(const WeylScalars& psi) {
  Scalar<ComplexDataVector> i{};
  Scalar<ComplexDataVector> j{};
  invariants(make_not_null(&i), make_not_null(&j), psi);
  return i;
}

Scalar<ComplexDataVector> invariant_j(const WeylScalars& psi) {
  Scalar<ComplexDataVector> i{};
  Scalar<ComplexDataVector> j{};
  invariants(make_not_null(&i), make_not_null(&j), psi);
  return j;
}
}  // namespace gr::np
