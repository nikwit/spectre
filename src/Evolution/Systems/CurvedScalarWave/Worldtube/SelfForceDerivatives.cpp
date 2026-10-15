// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SelfForceDerivatives.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "Utilities/Autodiff/Autodiff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace CurvedScalarWave::Worldtube {
#ifdef SPECTRE_AUTODIFF
namespace {
namespace ad = autodiff;
using TimeDual = ad::HigherOrderDual<3, double>;
using SpatialDual = ad::HigherOrderDual<2, double>;
template <class T>
using Vector = std::array<T, 4>;
template <class T>
using Matrix = std::array<Vector<T>, 4>;
template <class T>
using Connection = std::array<Matrix<T>, 4>;

template <class T>
struct Geometry {
  Matrix<T> metric{};
  Matrix<T> inverse{};
  Connection<T> christoffel{};
  Vector<T> contracted_christoffel{};
};

// The analytic Schwarzschild Kerr-Schild metric and its first spatial
// derivative. Evaluating with higher-order dual numbers differentiates these
// expressions further, including derivatives of the connection, without
// finite differences or separate generated third-metric-derivative code.
template <class T>
Geometry<T> geometry(const std::array<T, 3>& x, const double mass) {
  using std::sqrt;
  const T r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const T h = mass / r;
  const Vector<T> l{{T{1.}, x[0] / r, x[1] / r, x[2] / r}};
  auto l_up = l;
  l_up[0] = -1.;
  Geometry<T> result{};
  std::array<Matrix<T>, 4> d_metric{};
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      const double eta = a == b ? (a == 0 ? -1. : 1.) : 0.;
      result.metric.at(a).at(b) = eta + 2. * h * l.at(a) * l.at(b);
      result.inverse.at(a).at(b) = eta - 2. * h * l_up.at(a) * l_up.at(b);
      for (size_t i = 0; i < 3; ++i) {
        const T d_h = -h * l.at(i + 1) / r;
        const T d_l_a =
            a == 0 ? T{0.}
                   : ((a == i + 1 ? 1. : 0.) - l.at(a) * l.at(i + 1)) / r;
        const T d_l_b =
            b == 0 ? T{0.}
                   : ((b == i + 1 ? 1. : 0.) - l.at(b) * l.at(i + 1)) / r;
        d_metric.at(i + 1).at(a).at(b) =
            2. *
            (d_h * l.at(a) * l.at(b) + h * (d_l_a * l.at(b) + l.at(a) * d_l_b));
      }
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        auto& connection = result.christoffel.at(a).at(b).at(c);
        for (size_t d = 0; d < 4; ++d) {
          connection +=
              0.5 * result.inverse.at(a).at(d) *
              (d_metric.at(b).at(d).at(c) + d_metric.at(c).at(d).at(b) -
               d_metric.at(d).at(b).at(c));
        }
        result.contracted_christoffel.at(a) +=
            result.inverse.at(b).at(c) * connection;
      }
    }
  }
  return result;
}

struct FieldTimeDerivatives {
  double dt2{};
  double dt3{};
  double dt4{};
  std::array<double, 3> d_dt2{};
  std::array<double, 3> d_dt3{};
  std::array<std::array<double, 3>, 3> d2_dt2{};
};

// Evaluate the stationary wave equation with the supplied local Cauchy data.
// The dual-number seeds specify the spatial directions to differentiate.
template <class T>
T wave_rhs(const Geometry<T>& g, const T& dt_psi, const std::array<T, 3>& d_psi,
           const std::array<T, 3>& d_dt_psi,
           const tnsr::ii<double, 3>& d2_psi) {
  T result = g.contracted_christoffel[0] * dt_psi;
  for (size_t i = 0; i < 3; ++i) {
    result += -2. * g.inverse[0].at(i + 1) * d_dt_psi.at(i) +
              g.contracted_christoffel.at(i + 1) * d_psi.at(i);
    for (size_t j = 0; j < 3; ++j) {
      result -= g.inverse.at(i + 1).at(j + 1) * d2_psi.get(i, j);
    }
  }
  return result / g.inverse[0][0];
}

FieldTimeDerivatives close_wave_equation(const RegularFieldDerivatives& field,
                                         const tnsr::I<double, 3>& position,
                                         const double mass) {
  const auto dt2_at = [&](const std::array<SpatialDual, 3>& dx) -> SpatialDual {
    auto x = dx;
    SpatialDual dt_psi{field.dt_psi};
    std::array<SpatialDual, 3> d_psi{}, d_dt_psi{};
    for (size_t i = 0; i < 3; ++i) {
      x.at(i) += position.get(i);
      d_psi.at(i) = field.d_psi.get(i);
      d_dt_psi.at(i) = field.d_dt_psi.get(i);
      dt_psi += field.d_dt_psi.get(i) * dx.at(i);
      for (size_t j = 0; j < 3; ++j) {
        d_psi.at(i) += field.d2_psi.get(i, j) * dx.at(j);
        d_dt_psi.at(i) += field.d2_dt_psi.get(i, j) * dx.at(j);
        dt_psi += 0.5 * field.d2_dt_psi.get(i, j) * dx.at(i) * dx.at(j);
      }
    }
    return wave_rhs(geometry(x, mass), dt_psi, d_psi, d_dt_psi, field.d2_psi);
  };
  FieldTimeDerivatives result{};
  // As in the coordinate-map Hessian implementation, seed two spatial
  // directions and extract their mixed derivative from a second-order dual.
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      std::array<SpatialDual, 3> dx{};
      ad::seed<1>(dx.at(i), 1.);
      ad::seed<2>(dx.at(j), 1.);
      const SpatialDual dt2 = dt2_at(dx);
      result.dt2 = ad::val(dt2);
      result.d_dt2.at(i) = ad::derivative<1>(dt2);
      result.d2_dt2.at(i).at(j) = ad::derivative<2>(dt2);
      result.d2_dt2.at(j).at(i) = result.d2_dt2.at(i).at(j);
    }
  }
  // Stationarity lets us apply the same wave equation to psi_t. Only
  // psi_ttt and its gradient are needed, using the derivatives of psi_tt
  // just computed above. These mixed/time derivatives are not truncated.
  for (size_t k = 0; k < 3; ++k) {
    std::array<SpatialDual, 3> dx{};
    ad::seed<1>(dx.at(k), 1.);
    auto x = dx;
    SpatialDual dt2{result.dt2};
    std::array<SpatialDual, 3> d_dt_psi{}, d_dt2{};
    for (size_t i = 0; i < 3; ++i) {
      x.at(i) += position.get(i);
      dt2 += result.d_dt2.at(i) * dx.at(i);
      d_dt_psi.at(i) = field.d_dt_psi.get(i);
      d_dt2.at(i) = result.d_dt2.at(i);
      for (size_t j = 0; j < 3; ++j) {
        d_dt_psi.at(i) += field.d2_dt_psi.get(i, j) * dx.at(j);
        d_dt2.at(i) += result.d2_dt2.at(i).at(j) * dx.at(j);
      }
    }
    const SpatialDual dt3 =
        wave_rhs(geometry(x, mass), dt2, d_dt_psi, d_dt2, field.d2_dt_psi);
    result.dt3 = ad::val(dt3);
    result.d_dt3.at(k) = ad::derivative<1>(dt3);
  }
  std::array<double, 3> x{};
  tnsr::ii<double, 3> d2_dt2{};
  for (size_t i = 0; i < 3; ++i) {
    x.at(i) = position.get(i);
    for (size_t j = i; j < 3; ++j) {
      d2_dt2.get(i, j) = result.d2_dt2.at(i).at(j);
    }
  }
  result.dt4 = wave_rhs(geometry(x, mass), result.dt3, result.d_dt2,
                        result.d_dt3, d2_dt2);
  return result;
}

struct FieldOnWorldline {
  TimeDual psi{};
  Vector<TimeDual> gradient{};
};

// Compose the spacetime Taylor polynomial with the particle trajectory.
// Terms of total differential order four suffice for gradient derivatives
// through order three. Higher spatial derivatives of psi and psi_t are zero
// in the quadratic closure, not derivatives of psi_tt or psi_ttt.
FieldOnWorldline field_on_worldline(const RegularFieldDerivatives& field,
                                    const FieldTimeDerivatives& higher,
                                    const std::array<TimeDual, 3>& dx,
                                    const TimeDual& t) {
  const TimeDual t2 = t * t;
  const TimeDual t3 = t2 * t;
  FieldOnWorldline result{};
  result.psi = field.psi + field.dt_psi * t + 0.5 * higher.dt2 * t2 +
               higher.dt3 * t3 / 6.;
  result.gradient[0] = field.dt_psi + higher.dt2 * t + 0.5 * higher.dt3 * t2 +
                       higher.dt4 * t3 / 6.;
  for (size_t i = 0; i < 3; ++i) {
    result.psi += (field.d_psi.get(i) + field.d_dt_psi.get(i) * t +
                   0.5 * higher.d_dt2.at(i) * t2) *
                  dx.at(i);
    result.gradient[0] += (field.d_dt_psi.get(i) + higher.d_dt2.at(i) * t +
                           0.5 * higher.d_dt3.at(i) * t2) *
                          dx.at(i);
    result.gradient.at(i + 1) = field.d_psi.get(i) + field.d_dt_psi.get(i) * t +
                                0.5 * higher.d_dt2.at(i) * t2 +
                                higher.d_dt3.at(i) * t3 / 6.;
    for (size_t j = 0; j < 3; ++j) {
      result.psi += 0.5 *
                    (field.d2_psi.get(i, j) + field.d2_dt_psi.get(i, j) * t) *
                    dx.at(i) * dx.at(j);
      result.gradient[0] +=
          0.5 * (field.d2_dt_psi.get(i, j) + higher.d2_dt2.at(i).at(j) * t) *
          dx.at(i) * dx.at(j);
      result.gradient.at(i + 1) +=
          (field.d2_psi.get(i, j) + field.d2_dt_psi.get(i, j) * t +
           0.5 * higher.d2_dt2.at(i).at(j) * t2) *
          dx.at(j);
    }
  }
  return result;
}

TimeDual time_derivative(const TimeDual& f, const TimeDual& t) {
  return ad::derivative<1>(f) + ad::derivative<2>(f) * t +
         0.5 * ad::derivative<3>(f) * t * t;
}
TimeDual time_integral(const TimeDual& f, const TimeDual& t) {
  return ad::derivative<0>(f) * t + 0.5 * ad::derivative<1>(f) * t * t +
         ad::derivative<2>(f) * t * t * t / 6.;
}

struct ForceOnWorldline {
  Vector<TimeDual> u{};
  Vector<TimeDual> f{};
  std::array<TimeDual, 3> acceleration{};
  Connection<TimeDual> christoffel{};
};

ForceOnWorldline force_on_worldline(
    const RegularFieldDerivatives& field, const FieldTimeDerivatives& higher,
    const tnsr::I<double, 3>& position, const std::array<TimeDual, 3>& x,
    const std::array<TimeDual, 3>& v, const TimeDual& t, const double charge,
    const double rest_mass, const bool evolve_mass,
    const double black_hole_mass, const TimeDual& turn_on) {
  auto dx = x;
  for (size_t i = 0; i < 3; ++i) {
    dx.at(i) -= position.get(i);
  }
  const auto regular = field_on_worldline(field, higher, dx, t);
  const auto g = geometry(x, black_hole_mass);
  const Vector<TimeDual> v4{{TimeDual{1.}, v[0], v[1], v[2]}};
  TimeDual norm{0.};
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      norm -= g.metric.at(a).at(b) * v4.at(a) * v4.at(b);
    }
  }
  ASSERT(ad::derivative<0>(norm) > 0.,
         "The particle worldline must be timelike.");
  const TimeDual u0 = 1. / sqrt(norm);
  const TimeDual mass =
      evolve_mass ? rest_mass - charge * regular.psi : TimeDual{rest_mass};
  ASSERT(ad::derivative<0>(mass) > 0., "The particle mass must be positive.");
  ForceOnWorldline result{};
  result.christoffel = g.christoffel;
  for (size_t a = 0; a < 4; ++a) {
    result.u.at(a) = u0 * v4.at(a);
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      result.f.at(a) +=
          turn_on * charge / mass *
          (g.inverse.at(a).at(b) + result.u.at(a) * result.u.at(b)) *
          regular.gradient.at(b);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    result.acceleration.at(i) =
        (result.f.at(i + 1) - v.at(i) * result.f[0]) / (u0 * u0);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        result.acceleration.at(i) -= (g.christoffel.at(i + 1).at(a).at(b) -
                                      v.at(i) * g.christoffel[0].at(a).at(b)) *
                                     v4.at(a) * v4.at(b);
      }
    }
  }
  return result;
}

Vector<TimeDual> covariant_derivative(const Vector<TimeDual>& f,
                                      const ForceOnWorldline& worldline,
                                      const TimeDual& t) {
  Vector<TimeDual> result{};
  for (size_t a = 0; a < 4; ++a) {
    result.at(a) = worldline.u[0] * time_derivative(f.at(a), t);
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        result.at(a) += worldline.christoffel.at(a).at(b).at(c) *
                        worldline.u.at(b) * f.at(c);
      }
    }
  }
  return result;
}
}  // namespace
#endif  // SPECTRE_AUTODIFF

RegularFieldDerivatives reconstruct_regular_field(
    const double boundary_psi_monopole, const double boundary_dt_psi_monopole,
    const tnsr::i<double, 3>& psi_dipole,
    const tnsr::i<double, 3>& dt_psi_dipole,
    const tnsr::ii<double, 3>& psi_quadrupole,
    const tnsr::ii<double, 3>& dt_psi_quadrupole, const double psi0,
    const double total_dt_psi0, const tnsr::I<double, 3>& velocity,
    const double worldtube_radius) {
  ASSERT(worldtube_radius > 0., "The worldtube radius must be positive.");
  RegularFieldDerivatives result{psi0,           total_dt_psi0,
                                 psi_dipole,     dt_psi_dipole,
                                 psi_quadrupole, dt_psi_quadrupole};
  for (size_t i = 0; i < 3; ++i) {
    result.dt_psi -= velocity.get(i) * psi_dipole.get(i);
  }
  const double radius_squared = worldtube_radius * worldtube_radius;
  for (size_t i = 0; i < 3; ++i) {
    result.d2_psi.get(i, i) += (boundary_psi_monopole - psi0) / radius_squared;
    result.d2_dt_psi.get(i, i) +=
        (boundary_dt_psi_monopole - result.dt_psi) / radius_squared;
  }
  for (size_t i = 0; i < result.d2_psi.size(); ++i) {
    result.d2_psi[i] *= 2.;
    result.d2_dt_psi[i] *= 2.;
  }
  return result;
}

#ifdef SPECTRE_AUTODIFF
SelfForceDerivatives second_order_self_force_derivatives(
    const RegularFieldDerivatives& field, const tnsr::I<double, 3>& position,
    const tnsr::I<double, 3>& velocity, const double charge,
    const double rest_mass, const bool evolve_mass,
    const double black_hole_mass, const double time, const double turn_on_time,
    const double turn_on_interval) {
  const auto higher = close_wave_equation(field, position, black_hole_mass);
  TimeDual t{0.};
  ad::seed<1>(t, 1.);
  ad::seed<2>(t, 1.);
  ad::seed<3>(t, 1.);
  TimeDual turn_on{0.};
  if (time > turn_on_time) {
    ASSERT(turn_on_interval > 0.,
           "The self-force turn-on interval must be positive.");
    const TimeDual s = (time - turn_on_time + t) / turn_on_interval;
    turn_on = 1. - exp(-s * s * s * s);
  }
  std::array<TimeDual, 3> x{}, v{};
  for (size_t i = 0; i < 3; ++i) {
    x.at(i) = position.get(i);
    v.at(i) = velocity.get(i);
  }
  // Differentiate the local equations of motion with higher-order duals.
  // Each pass determines one more trajectory derivative; three passes give
  // everything needed for d_t^3 f at the current time. The formal time
  // variable has value zero throughout: this does not advance the orbit.
  for (size_t order = 0; order < 3; ++order) {
    const auto worldline =
        force_on_worldline(field, higher, position, x, v, t, charge, rest_mass,
                           evolve_mass, black_hole_mass, turn_on);
    for (size_t i = 0; i < 3; ++i) {
      x.at(i) = position.get(i) + time_integral(v.at(i), t);
      v.at(i) =
          velocity.get(i) + time_integral(worldline.acceleration.at(i), t);
    }
  }
  const auto worldline =
      force_on_worldline(field, higher, position, x, v, t, charge, rest_mass,
                         evolve_mass, black_hole_mass, turn_on);
  const auto j = covariant_derivative(worldline.f, worldline, t);
  const auto k = covariant_derivative(j, worldline, t);
  SelfForceDerivatives result{};
  for (size_t i = 0; i < 3; ++i) {
    result.coordinate_acceleration.get(i) =
        ad::derivative<0>(worldline.acceleration.at(i));
  }
  for (size_t a = 0; a < 4; ++a) {
    result.f.get(a) = ad::derivative<0>(worldline.f.at(a));
    result.dt_f.get(a) = ad::derivative<1>(worldline.f.at(a));
    result.Du_f.get(a) = ad::derivative<0>(j.at(a));
    result.dt_Du_f.get(a) = ad::derivative<1>(j.at(a));
    result.DuDu_f.get(a) = ad::derivative<0>(k.at(a));
    result.dt_DuDu_f.get(a) = ad::derivative<1>(k.at(a));
  }
  return result;
}
#else
SelfForceDerivatives second_order_self_force_derivatives(
    const RegularFieldDerivatives& /*field*/,
    const tnsr::I<double, 3>& /*position*/,
    const tnsr::I<double, 3>& /*velocity*/, const double /*charge*/,
    const double /*rest_mass*/, const bool /*evolve_mass*/,
    const double /*black_hole_mass*/, const double /*time*/,
    const double /*turn_on_time*/, const double /*turn_on_interval*/) {
  ERROR("Second-order self-force derivatives require SPECTRE_AUTODIFF=ON.");
}
#endif  // SPECTRE_AUTODIFF
}  // namespace CurvedScalarWave::Worldtube
