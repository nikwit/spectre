// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <tuple>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube {
/// Prototype for static, unrefined spherical shells. The two coefficients
/// describe d_r(q-q_initial) = K (q-q_initial) in the t/r gauge channels.
/// The transverse gauge channels retain the zero-rate radiation condition.
struct RadialResponseGauge {
  struct RadialResponse {
    using type = std::array<double, 2>;
    static constexpr Options::String help{
        "Fixed [Kt, Kr] radial response coefficients (inverse length). "
        "Initial q and radial derivative are captured once; live GH source "
        "terms are retained. Static unrefined 3D spherical shells only."};
  };
  using options = tmpl::list<RadialResponse>;
  static constexpr Options::String help{
      "Initial-profile radial gauge response."};
  RadialResponseGauge() = default;
  explicit RadialResponseGauge(std::array<double, 2> values)
      : coefficients(values) {}
  std::array<double, 2> coefficients{};
};

struct RadialGaugeFaceData {
  tnsr::a<DataVector, 3> q{}, dr_q{}, initial_q{}, initial_dr_q{};
  tnsr::I<DataVector, 3> radial_direction{};
  std::array<size_t, 3> extents{};
  double time{0.};
  void pup(PUP::er& p) {
    p | q;
    p | dr_q;
    p | initial_q;
    p | initial_dr_q;
    p | radial_direction;
    p | extents;
    p | time;
  }
};
inline bool operator==(const RadialGaugeFaceData& a,
                       const RadialGaugeFaceData& b) {
  return std::tie(a.q, a.dr_q, a.initial_q, a.initial_dr_q, a.radial_direction,
                  a.extents,
                  a.time) == std::tie(b.q, b.dr_q, b.initial_q, b.initial_dr_q,
                                      b.radial_direction, b.extents, b.time);
}
inline bool operator!=(const RadialGaugeFaceData& a,
                       const RadialGaugeFaceData& b) {
  return not(a == b);
}

/// q is formed with the metric-dependent inner-boundary null vector at
/// every volume point, then differentiated with the native spectral operator.
/// This includes all spatial derivatives of ell, s and gamma2. It uses the
/// same operation at startup and later, so the initial correction vanishes.
inline void update_radial_gauge_face_data(
    const gsl::not_null<std::optional<RadialGaugeFaceData>*> data,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& jac,
    const Direction<3>& direction, const double time) {
  if (direction != Direction<3>::lower_xi() or
      mesh.basis(0) != Spectral::Basis::Legendre or
      mesh.quadrature(0) != Spectral::Quadrature::GaussLobatto or
      mesh.basis(1) != Spectral::Basis::SphericalHarmonic or
      mesh.basis(2) != Spectral::Basis::SphericalHarmonic) {
    ERROR(
        "RadialResponse requires the lower radial LGL face of a spherical "
        "shell.");
  }
  const size_t n = mesh.number_of_grid_points(), nr = mesh.extents(0),
               nf = n / nr;
  const std::array<size_t, 3> ext{{nr, mesh.extents(1), mesh.extents(2)}};
  if (data->has_value() and data->value().extents != ext) {
    ERROR("RadialResponse cannot change mesh resolution after initialization.");
  }
  const auto inverse = determinant_and_inverse(metric).second;
  const DataVector alpha = 1. / sqrt(-get<0, 0>(inverse));
  DataVector radius(n, 0.);
  tnsr::I<DataVector, 3> radial(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    radial.get(i) = coordinates.get(i) - center.get(i);
    radius += square(radial.get(i));
  }
  radius = sqrt(radius);
  if (min(radius) <= 0.) {
    ERROR("RadialResponse encountered the excision center.");
  }
  for (auto& c : radial) {
    c /= radius;
  }
  tnsr::A<DataVector, 3> tn(n, 0.), ell(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    tn.get(a) = -alpha * inverse.get(a, 0);
  }
  DataVector norm(n, 0.);
  tnsr::I<DataVector, 3> normal(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      normal.get(i) -=
          (inverse.get(i + 1, j + 1) + tn.get(i + 1) * tn.get(j + 1)) *
          radial.get(j);
    }
    norm -= radial.get(i) * normal.get(i);
  }
  norm = sqrt(norm);
  for (auto& c : normal) {
    c /= norm;
  }
  get<0>(ell) = get<0>(tn) / sqrt(2.);
  for (size_t i = 0; i < 3; ++i) {
    ell.get(i + 1) = (tn.get(i + 1) + normal.get(i)) / sqrt(2.);
  }
  tnsr::a<DataVector, 3> q(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      DataVector u = pi.get(a, b) - get(gamma2) * metric.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        u -= normal.get(i) * phi.get(i, a, b);
      }
      q.get(a) += ell.get(b) * u;
    }
  }
  const auto dq = partial_derivative(q, mesh, jac);
  RadialGaugeFaceData next{};
  next.q = tnsr::a<DataVector, 3>(nf, 0.);
  next.dr_q = tnsr::a<DataVector, 3>(nf, 0.);
  next.radial_direction = tnsr::I<DataVector, 3>(nf, 0.);
  next.extents = ext;
  next.time = time;
  for (size_t p = 0; p < nf; ++p) {
    const size_t v = p * nr;
    if (std::abs(radius[v] - radius[0]) > 1.e-10 * radius[0]) {
      ERROR("RadialResponse face is not a coordinate sphere.");
    }
    for (size_t i = 0; i < 3; ++i) {
      next.radial_direction.get(i)[p] = radial.get(i)[v];
    }
    for (size_t a = 0; a < 4; ++a) {
      next.q.get(a)[p] = q.get(a)[v];
      for (size_t i = 0; i < 3; ++i) {
        next.dr_q.get(a)[p] += radial.get(i)[v] * dq.get(i, a)[v];
      }
    }
  }
  if (data->has_value()) {
    if (data->value().radial_direction != next.radial_direction) {
      ERROR("RadialResponse face coordinates changed after initialization.");
    }
    next.initial_q = data->value().initial_q;
    next.initial_dr_q = data->value().initial_dr_q;
  } else {
    next.initial_q = next.q;
    next.initial_dr_q = next.dr_q;
  }
  *data = std::move(next);
}
}  // namespace gh::worldtube
