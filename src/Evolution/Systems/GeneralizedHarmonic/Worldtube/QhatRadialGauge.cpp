// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/QhatRadialGauge.hpp"

#include <cmath>
#include <tuple>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Slice.hpp"
#include "Domain/CoordinateMaps/TimeDependent/RotationMatrixHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gh::worldtube {
namespace {
Matrix rotation(const double time,
                const domain::FunctionsOfTimeMap& functions_of_time) {
  if (functions_of_time.count("Rotation") != 0) {
    return rotation_matrix<3>(time, *functions_of_time.at("Rotation"));
  }
  Matrix result(3, 3, 0.);
  for (size_t i = 0; i < 3; ++i) {
    result(i, i) = 1.;
  }
  return result;
}

struct Basis {
  tnsr::AA<DataVector, 3> inverse;
  tnsr::II<DataVector, 3> spatial_inverse;
  tnsr::A<DataVector, 3> tn, ell;
  tnsr::I<DataVector, 3> normal;
};
Basis basis(const tnsr::aa<DataVector, 3>& metric,
            const tnsr::I<DataVector, 3>& radial) {
  const size_t n = get<0, 0>(metric).size();
  Basis b{determinant_and_inverse(metric).second, tnsr::II<DataVector, 3>(n),
          tnsr::A<DataVector, 3>(n), tnsr::A<DataVector, 3>(n),
          tnsr::I<DataVector, 3>(n, 0.)};
  const DataVector alpha = 1. / sqrt(-get<0, 0>(b.inverse));
  for (size_t a = 0; a < 4; ++a) {
    b.tn.get(a) = -alpha * b.inverse.get(a, 0);
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      b.spatial_inverse.get(i, j) =
          b.inverse.get(i + 1, j + 1) + b.tn.get(i + 1) * b.tn.get(j + 1);
    }
  }
  DataVector norm(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      b.normal.get(i) -= b.spatial_inverse.get(i, j) * radial.get(j);
    }
    norm -= radial.get(i) * b.normal.get(i);
  }
  norm = sqrt(norm);
  b.ell = b.tn;
  for (size_t i = 0; i < 3; ++i) {
    b.normal.get(i) /= norm;
    b.ell.get(i + 1) += b.normal.get(i);
  }
  for (auto& component : b.ell) {
    component /= sqrt(2.);
  }
  return b;
}

// Directional metric variation with the coordinate radial direction held
// fixed. Applies to a radial derivative or to a BC's metric RHS correction.
std::pair<tnsr::I<DataVector, 3>, tnsr::A<DataVector, 3>> basis_variation(
    const Basis& b, const tnsr::aa<DataVector, 3>& dg) {
  const size_t n = get<0, 0>(dg).size();
  tnsr::I<DataVector, 3> dn(n, 0.);
  tnsr::A<DataVector, 3> dell(n, 0.);
  DataVector tt(n, 0.), nn(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t c = 0; c < 4; ++c) {
      tt += b.tn.get(a) * b.tn.get(c) * dg.get(a, c);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      nn += b.normal.get(i) * b.normal.get(j) * dg.get(i + 1, j + 1);
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    dell.get(a) = -.5 * b.tn.get(a) * tt;
    for (size_t c = 0; c < 4; ++c) {
      for (size_t d = 0; d < 4; ++d) {
        dell.get(a) -= b.inverse.get(a, c) * dg.get(c, d) * b.tn.get(d);
      }
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    dn.get(i) = .5 * b.normal.get(i) * nn;
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        dn.get(i) -= b.spatial_inverse.get(i, j) * dg.get(j + 1, k + 1) *
                     b.normal.get(k);
      }
    }
    dell.get(i + 1) += dn.get(i);
  }
  for (auto& component : dell) {
    component /= sqrt(2.);
  }
  return {std::move(dn), std::move(dell)};
}

tnsr::I<DataVector, 3> radial_direction(
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const double radius) {
  auto result = coordinates;
  DataVector norm(get<0>(coordinates).size(), 0.);
  for (size_t i = 0; i < 3; ++i) {
    result.get(i) = (coordinates.get(i) - center.get(i)) / radius;
    norm += square(result.get(i));
  }
  if (max(abs(norm - 1.)) > 1.e-9) {
    ERROR("QhatRadialResponse requires a fixed-radius coordinate sphere.");
  }
  return result;
}
}  // namespace

void QhatRadialGaugeData::pup(PUP::er& p) {
  p | initial_q;
  p | initial_dr_q;
  p | corotating_radial_direction;
  p | extents;
  p | radius;
}
bool operator==(const QhatRadialGaugeData& a, const QhatRadialGaugeData& b) {
  return std::tie(a.initial_q, a.initial_dr_q, a.corotating_radial_direction,
                  a.extents, a.radius) ==
         std::tie(b.initial_q, b.initial_dr_q, b.corotating_radial_direction,
                  b.extents, b.radius);
}
bool operator!=(const QhatRadialGaugeData& a, const QhatRadialGaugeData& b) {
  return not(a == b);
}

std::pair<tnsr::a<DataVector, 3>, tnsr::a<DataVector, 3>> qhat_radial_data(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi,
    const tnsr::I<DataVector, 3>& radial) {
  const size_t n = get<0, 0>(metric).size();
  const auto b = basis(metric, radial);
  tnsr::aa<DataVector, 3> drg(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t c = a; c < 4; ++c) {
      for (size_t i = 0; i < 3; ++i) {
        drg.get(a, c) += radial.get(i) * d_metric.get(i, a, c);
      }
    }
  }
  const auto [dn, dell] = basis_variation(b, drg);
  tnsr::a<DataVector, 3> q(n, 0.), drq(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t c = 0; c < 4; ++c) {
      DataVector w = pi.get(a, c), drw(n, 0.);
      for (size_t i = 0; i < 3; ++i) {
        w -= b.normal.get(i) * phi.get(i, a, c);
        drw += radial.get(i) * d_pi.get(i, a, c) - dn.get(i) * phi.get(i, a, c);
        for (size_t j = 0; j < 3; ++j) {
          drw -= radial.get(j) * b.normal.get(i) * d_phi.get(j, i, a, c);
        }
      }
      q.get(a) += b.ell.get(c) * w;
      drq.get(a) += b.ell.get(c) * drw + dell.get(c) * w;
    }
  }
  return {std::move(q), std::move(drq)};
}

void initialize_qhat_radial_gauge(
    const gsl::not_null<std::optional<QhatRadialGaugeData>*> data,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Direction<3>& direction, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  if (direction != Direction<3>::lower_xi() or
      mesh.basis(0) != Spectral::Basis::Legendre or
      mesh.quadrature(0) != Spectral::Quadrature::GaussLobatto or
      mesh.basis(1) != Spectral::Basis::SphericalHarmonic or
      mesh.basis(2) != Spectral::Basis::SphericalHarmonic) {
    ERROR(
        "QhatRadialResponse requires the lower radial LGL face of a spherical "
        "shell.");
  }
  const std::array<size_t, 3> extents{
      {mesh.extents(0), mesh.extents(1), mesh.extents(2)}};
  if (data->has_value()) {
    if (data->value().extents != extents) {
      ERROR(
          "QhatRadialResponse cannot change resolution after initialization.");
    }
    return;
  }
  const auto face = [&](const auto& tensor) {
    return data_on_slice(tensor, mesh.extents(), 0, 0);
  };
  const auto xyz = face(coordinates);
  double radius = 0.;
  for (size_t i = 0; i < 3; ++i) {
    radius += square(xyz.get(i)[0] - center.get(i));
  }
  radius = sqrt(radius);
  if (radius <= 0.) {
    ERROR("QhatRadialResponse encountered the excision center.");
  }
  const auto radial = radial_direction(xyz, center, radius);
  const auto [q, drq] = qhat_radial_data(
      face(metric), face(pi), face(phi),
      face(partial_derivative(metric, mesh, inverse_jacobian)),
      face(partial_derivative(pi, mesh, inverse_jacobian)),
      face(partial_derivative(phi, mesh, inverse_jacobian)), radial);
  const size_t n = get<0>(q).size();
  QhatRadialGaugeData next{q, drq, tnsr::I<DataVector, 3>(n, 0.), extents,
                           radius};
  const auto rot = rotation(time, functions_of_time);
  for (size_t i = 0; i < 3; ++i) {
    next.initial_q.get(i + 1) = 0.;
    next.initial_dr_q.get(i + 1) = 0.;
    for (size_t j = 0; j < 3; ++j) {
      next.initial_q.get(i + 1) += rot(j, i) * q.get(j + 1);
      next.initial_dr_q.get(i + 1) += rot(j, i) * drq.get(j + 1);
      next.corotating_radial_direction.get(i) += rot(j, i) * radial.get(j);
    }
  }
  *data = std::move(next);
}

tnsr::a<DataVector, 3> qhat_radial_gauge_correction(
    const QhatRadialGaugeData& initial,
    const std::array<double, 3>& coefficients,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi,
    const tnsr::I<DataVector, 3>& coordinates, const tnsr::I<double, 3>& center,
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity,
    const Scalar<DataVector>& gamma2,
    const tnsr::aa<DataVector, 3>& metric_correction, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  const size_t n = get<0>(coordinates).size();
  if (get<0>(initial.initial_q).size() != n) {
    ERROR("QhatRadialResponse face resolution changed.");
  }
  const auto radial = radial_direction(coordinates, center, initial.radius);
  const auto rot = rotation(time, functions_of_time);
  auto q0 = initial.initial_q;
  auto drq0 = initial.initial_dr_q;
  for (size_t i = 0; i < 3; ++i) {
    DataVector expected(n, 0.);
    q0.get(i + 1) = 0.;
    drq0.get(i + 1) = 0.;
    for (size_t j = 0; j < 3; ++j) {
      expected += rot(i, j) * initial.corotating_radial_direction.get(j);
      q0.get(i + 1) += rot(i, j) * initial.initial_q.get(j + 1);
      drq0.get(i + 1) += rot(i, j) * initial.initial_dr_q.get(j + 1);
    }
    if (max(abs(expected - radial.get(i))) > 1.e-9) {
      ERROR(
          "QhatRadialResponse angular labels must follow the Rotation function "
          "of time.");
    }
  }
  const auto [q, drq] =
      qhat_radial_data(metric, pi, phi, d_metric, d_pi, d_phi, radial);
  const auto b = basis(metric, radial);
  DataVector speed(n, 0.), delta_radial(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    speed -= radial.get(i) * (shift.get(i) + get(lapse) * b.normal.get(i));
    if (mesh_velocity.has_value()) {
      speed -= radial.get(i) * mesh_velocity->get(i);
    }
    delta_radial += radial.get(i) * (q.get(i + 1) - q0.get(i + 1));
  }
  tnsr::a<DataVector, 3> correction(n, 0.);
  get<0>(correction) =
      speed * (get<0>(drq) - get<0>(drq0) -
               coefficients[0] / initial.radius * (get<0>(q) - get<0>(q0)));
  for (size_t i = 0; i < 3; ++i) {
    correction.get(i + 1) =
        speed *
        (drq.get(i + 1) - drq0.get(i + 1) -
         (coefficients[2] * (q.get(i + 1) - q0.get(i + 1)) +
          (coefficients[1] - coefficients[2]) * radial.get(i) * delta_radial) /
             initial.radius);
  }
  const auto [dn, dell] = basis_variation(b, metric_correction);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t c = 0; c < 4; ++c) {
      DataVector w = pi.get(a, c);
      correction.get(a) -=
          get(gamma2) * b.ell.get(c) * metric_correction.get(a, c);
      for (size_t i = 0; i < 3; ++i) {
        w -= b.normal.get(i) * phi.get(i, a, c);
        correction.get(a) += b.ell.get(c) * dn.get(i) * phi.get(i, a, c);
      }
      correction.get(a) -= dell.get(c) * w;
    }
  }
  return correction;
}
}  // namespace gh::worldtube
