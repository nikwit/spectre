// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryEvolution.hpp"

#include <cmath>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gh::worldtube {
PhysicalTraceGeometry spherical_trace_geometry(
    const tnsr::aa<DataVector, 3>& metric,
    const tnsr::aa<DataVector, 3>& dt_metric,
    const tnsr::I<DataVector, 3>& coordinates) {
  const auto& x = get<0>(coordinates);
  const auto& y = get<1>(coordinates);
  const auto& z = get<2>(coordinates);
  const DataVector rho = sqrt(square(x) + square(y));
  const DataVector r = sqrt(square(rho) + square(z));
  if (min(rho) <= 0. or max(r) - min(r) > 1.e-10 * max(r)) {
    ERROR(
        "Physical trace prototype requires a full sphere centered at the "
        "inertial origin, without polar collocation points.");
  }
  tnsr::i<DataVector, 3> raw_normal(x.size());
  std::array<tnsr::I<DataVector, 3>, 2> seeds{
      tnsr::I<DataVector, 3>(x.size()), tnsr::I<DataVector, 3>(x.size(), 0.)};
  for (size_t i = 0; i < 3; ++i) {
    raw_normal.get(i) = -coordinates.get(i) / r;
  }
  get<0>(seeds[0]) = z * x / (r * rho);
  get<1>(seeds[0]) = z * y / (r * rho);
  get<2>(seeds[0]) = -rho / r;
  get<0>(seeds[1]) = -y / rho;
  get<1>(seeds[1]) = x / rho;
  return {metric, dt_metric, raw_normal, seeds};
}

PhysicalTraceGeometry mapped_sphere_trace_geometry(
    const tnsr::aa<DataVector, 3>& metric,
    const tnsr::aa<DataVector, 3>& dt_metric,
    const tnsr::I<DataVector, 3>& coordinates,
    const tnsr::I<DataVector, 3>& mesh_velocity, const size_t l_max) {
  const auto& sphere = ylm::get_spherepack_cache(l_max);
  const size_t n = get<0>(coordinates).size();
  if (n != sphere.physical_size()) {
    ERROR("Physical trace requires one full spherical-harmonic face.");
  }
  std::array<tnsr::I<DataVector, 3>, 2> tangents{tnsr::I<DataVector, 3>(n),
                                                 tnsr::I<DataVector, 3>(n)};
  auto dt_tangents = tangents;
  for (size_t i = 0; i < 3; ++i) {
    const auto dx = sphere.gradient(coordinates.get(i));
    const auto dv = sphere.gradient(mesh_velocity.get(i));
    for (size_t a = 0; a < 2; ++a) {
      tangents[a].get(i) = dx[a];
      dt_tangents[a].get(i) = dv[a];
    }
  }
  tnsr::i<DataVector, 3> normal(n), dt_normal(n);
  for (size_t i = 0; i < 3; ++i) {
    const size_t j = (i + 1) % 3, k = (i + 2) % 3;
    // The negative cross product points into the excised region. The
    // Pfaffian angular convention only rescales this unnormalized covector.
    normal.get(i) = -(tangents[0].get(j) * tangents[1].get(k) -
                      tangents[0].get(k) * tangents[1].get(j));
    dt_normal.get(i) = -(dt_tangents[0].get(j) * tangents[1].get(k) +
                         tangents[0].get(j) * dt_tangents[1].get(k) -
                         dt_tangents[0].get(k) * tangents[1].get(j) -
                         tangents[0].get(k) * dt_tangents[1].get(j));
  }
  return {metric, dt_metric, normal, tangents, dt_normal, dt_tangents};
}

PhysicalTraceGeometry::PhysicalTraceGeometry(
    const tnsr::aa<DataVector, 3>& metric,
    const tnsr::aa<DataVector, 3>& dt_metric,
    const tnsr::i<DataVector, 3>& normal_one_form,
    const std::array<tnsr::I<DataVector, 3>, 2>& tangent_seeds,
    const std::optional<tnsr::i<DataVector, 3>>& dt_normal_one_form,
    const std::optional<std::array<tnsr::I<DataVector, 3>, 2>>&
        dt_tangent_seeds) {
  const size_t n = get<0, 0>(metric).size();
  Scalar<DataVector> determinant(n);
  determinant_and_inverse(make_not_null(&determinant),
                          make_not_null(&inverse_metric), metric);
  tnsr::ii<DataVector, 3> spatial_metric(n);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      spatial_metric.get(i, j) = metric.get(i + 1, j + 1);
    }
  }
  determinant_and_inverse(make_not_null(&determinant),
                          make_not_null(&inverse_spatial_metric),
                          spatial_metric);
  get(lapse) = 1. / sqrt(-get<0, 0>(inverse_metric));
  get<0>(time_normal) = 1. / get(lapse);
  DataVector norm(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    shift.get(i) = square(get(lapse)) * inverse_metric.get(0, i + 1);
    time_normal.get(i + 1) = -shift.get(i) / get(lapse);
    normal_vector.get(i) = DataVector(n, 0.);
    for (size_t j = 0; j < 3; ++j) {
      normal_vector.get(i) +=
          inverse_spatial_metric.get(i, j) * normal_one_form.get(j);
    }
    norm += normal_vector.get(i) * normal_one_form.get(i);
  }
  norm = sqrt(norm);
  for (size_t i = 0; i < 3; ++i) {
    normal_vector.get(i) /= norm;
    normal_covector.get(i) = normal_one_form.get(i) / norm;
  }
  DataVector dt_gamma_ss(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      dt_gamma_ss += dt_metric.get(i + 1, j + 1) * normal_vector.get(i) *
                     normal_vector.get(j);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    dt_normal_vector.get(i) = 0.5 * normal_vector.get(i) * dt_gamma_ss;
    for (size_t j = 0; j < 3; ++j) {
      if (dt_normal_one_form.has_value()) {
        dt_normal_vector.get(i) +=
            (inverse_spatial_metric.get(i, j) -
             normal_vector.get(i) * normal_vector.get(j)) *
            dt_normal_one_form->get(j) / norm;
      }
      for (size_t k = 0; k < 3; ++k) {
        dt_normal_vector.get(i) -= inverse_spatial_metric.get(i, j) *
                                   dt_metric.get(j + 1, k + 1) *
                                   normal_vector.get(k);
      }
    }
  }

  auto e = tangent_seeds;
  std::array<tnsr::I<DataVector, 3>, 2> de{tnsr::I<DataVector, 3>(n, 0.),
                                           tnsr::I<DataVector, 3>(n, 0.)};
  if (dt_tangent_seeds.has_value()) {
    de = *dt_tangent_seeds;
  }
  for (size_t p = 0; p < 2; ++p) {
    if (p == 1) {
      DataVector dot(n, 0.), dt_dot(n, 0.);
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          dot += spatial_metric.get(i, j) * e[0].get(i) * e[1].get(j);
          dt_dot += (dt_metric.get(i + 1, j + 1) * e[0].get(i) +
                     spatial_metric.get(i, j) * de[0].get(i)) *
                    e[1].get(j);
          dt_dot += spatial_metric.get(i, j) * e[0].get(i) * de[1].get(j);
        }
      }
      for (size_t i = 0; i < 3; ++i) {
        e[1].get(i) -= dot * e[0].get(i);
        de[1].get(i) -= dt_dot * e[0].get(i) + dot * de[0].get(i);
      }
    }
    DataVector length(n, 0.), dt_length_squared(n, 0.);
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        length += spatial_metric.get(i, j) * e[p].get(i) * e[p].get(j);
        dt_length_squared +=
            dt_metric.get(i + 1, j + 1) * e[p].get(i) * e[p].get(j) +
            2. * spatial_metric.get(i, j) * de[p].get(i) * e[p].get(j);
      }
    }
    if (min(length) <= 0.) {
      ERROR("Physical trace tangent seeds must be independent.");
    }
    length = sqrt(length);
    for (size_t i = 0; i < 3; ++i) {
      de[p].get(i) = de[p].get(i) / length -
                     0.5 * e[p].get(i) * dt_length_squared / pow<3>(length);
      e[p].get(i) /= length;
    }
    DataVector normal_dot_e(n, 0.);
    for (size_t i = 0; i < 3; ++i) {
      normal_dot_e += normal_covector.get(i) * e[p].get(i);
    }
    if (max(abs(normal_dot_e)) > 1.e-10) {
      ERROR("Physical trace seeds must be tangent to the surface.");
    }
  }
  for (size_t p = 0; p < 2; ++p) {
    polarization[p] = tnsr::AA<DataVector, 3>(n, 0.);
    dt_polarization[p] = tnsr::AA<DataVector, 3>(n, 0.);
    lower_polarization[p] = tnsr::aa<DataVector, 3>(n, 0.);
  }
  const double c = 1. / sqrt(2.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      polarization[0].get(i + 1, j + 1) =
          c * (e[0].get(i) * e[0].get(j) - e[1].get(i) * e[1].get(j));
      polarization[1].get(i + 1, j + 1) =
          c * (e[0].get(i) * e[1].get(j) + e[1].get(i) * e[0].get(j));
      dt_polarization[0].get(i + 1, j + 1) =
          c * (de[0].get(i) * e[0].get(j) + e[0].get(i) * de[0].get(j) -
               de[1].get(i) * e[1].get(j) - e[1].get(i) * de[1].get(j));
      dt_polarization[1].get(i + 1, j + 1) =
          c * (de[0].get(i) * e[1].get(j) + e[0].get(i) * de[1].get(j) +
               de[1].get(i) * e[0].get(j) + e[1].get(i) * de[0].get(j));
    }
  }
  for (size_t p = 0; p < 2; ++p) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            lower_polarization[p].get(a, b) +=
                metric.get(a, i + 1) * metric.get(b, j + 1) *
                polarization[p].get(i + 1, j + 1);
          }
        }
      }
    }
  }
}

DataVector physical_trace_amplitudes(const tnsr::aa<DataVector, 3>& pi,
                                     const tnsr::iaa<DataVector, 3>& phi,
                                     const PhysicalTraceGeometry& geometry) {
  const size_t n = get<0, 0>(pi).size();
  DataVector result(2 * n, 0.);
  for (size_t p = 0; p < 2; ++p) {
    DataVector component(result.data() + p * n, n);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = 0; b < 4; ++b) {
        DataVector v = pi.get(a, b);
        for (size_t i = 0; i < 3; ++i) {
          v -= geometry.normal_vector.get(i) * phi.get(i, a, b);
        }
        component += geometry.polarization[p].get(a, b) * v;
      }
    }
  }
  return result;
}

tnsr::aa<DataVector, 3> physical_trace_tensor(
    const DataVector& amplitudes, const PhysicalTraceGeometry& geometry) {
  const size_t n = get(geometry.lapse).size();
  if (amplitudes.size() != 2 * n) {
    ERROR(
        "Physical trace state has wrong size. Initialize it before evolution; "
        "AMR of the surface trace is not supported.");
  }
  tnsr::aa<DataVector, 3> result(n, 0.);
  for (size_t p = 0; p < 2; ++p) {
    // Copy avoids a mutable alias to the evolved state.
    DataVector component(n);
    std::copy_n(amplitudes.data() + p * n, n, component.data());
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        result.get(a, b) +=
            component * geometry.lower_polarization[p].get(a, b);
      }
    }
  }
  return result;
}

DataVector physical_boundary_rhs(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi, const Scalar<DataVector>& gamma0,
    const Scalar<DataVector>& gamma2, const tnsr::a<DataVector, 3>& gauge_h,
    const tnsr::ab<DataVector, 3>& d_gauge_h,
    const tnsr::aa<DataVector, 3>& target_weyl,
    const PhysicalTraceGeometry& geometry,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity,
    const gh::gauges::DampedHarmonic* damped_harmonic,
    const tnsr::I<DataVector, 3>* coordinates) {
  using Dt = gh::TimeDerivative<gh::Solutions::all_solutions<3>, 3>;
  const size_t n = get<0, 0>(metric).size();
  Variables<Dt::temporary_tags> temp(n);
  tnsr::aa<DataVector, 3> dt_g(n), dt_pi(n), weyl(n, 0.);
  tnsr::iaa<DataVector, 3> dt_phi(n);
  const gh::gauges::Harmonic harmonic{};
  if (damped_harmonic != nullptr and coordinates == nullptr) {
    ERROR("Damped-harmonic trace evolution requires inertial coordinates.");
  }
  const gh::gauges::GaugeCondition& gauge =
      damped_harmonic == nullptr
          ? static_cast<const gh::gauges::GaugeCondition&>(harmonic)
          : *damped_harmonic;
  // Harmonic and non-rollon DampedHarmonic do not use the mesh. For an
  // explicit prescribed gauge restore H,dH below. Native DampedHarmonic uses
  // the mixed ghost fields, including their first-order metric derivatives.
  tmpl::as_pack<Dt::temporary_tags>([&]<typename... T>(tmpl::type_<T>...) {
    Dt::apply(
        make_not_null(&dt_g), make_not_null(&dt_pi), make_not_null(&dt_phi),
        make_not_null(&get<T>(temp))..., d_metric, d_pi, d_phi, metric, pi, phi,
        gamma0, Scalar<DataVector>(n, -1.), gamma2, gauge,
        Mesh<3>(2, Spectral::Basis::Legendre,
                Spectral::Quadrature::GaussLobatto),
        0.,
        coordinates == nullptr ? tnsr::I<DataVector, 3>(n, 0.) : *coordinates,
        InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>(
            n, 0.),
        mesh_velocity);
  });
  const auto& christoffel =
      get<gr::Tags::SpacetimeChristoffelFirstKind<DataVector, 3>>(temp);
  if (damped_harmonic == nullptr) {
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        DataVector source = -d_gauge_h.get(a, b) - d_gauge_h.get(b, a);
        if (a == 0) {
          source -= get(gamma0) * get(geometry.lapse) * gauge_h.get(b);
        }
        if (b == 0) {
          source -= get(gamma0) * get(geometry.lapse) * gauge_h.get(a);
        }
        for (size_t c = 0; c < 4; ++c) {
          source -= get(gamma0) * metric.get(a, b) *
                    geometry.time_normal.get(c) * gauge_h.get(c);
          for (size_t d = 0; d < 4; ++d) {
            source += 2. * christoffel.get(c, a, b) *
                      geometry.inverse_metric.get(c, d) * gauge_h.get(d);
          }
        }
        dt_pi.get(a, b) += get(geometry.lapse) * source;
      }
    }
  }
  if (mesh_velocity.has_value()) {
    // The native GH RHS includes its mesh-velocity constraint additions.
    // DG adds ordinary mesh advection separately; reproduce that here.
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t i = 0; i < 3; ++i) {
          dt_pi.get(a, b) += mesh_velocity->get(i) * d_pi.get(i, a, b);
          for (size_t j = 0; j < 3; ++j) {
            dt_phi.get(j, a, b) +=
                mesh_velocity->get(i) * d_phi.get(i, j, a, b);
          }
        }
      }
    }
  }
  const PhysicalScreen screen(metric, geometry.inverse_metric,
                              geometry.time_normal, geometry.normal_vector);
  const auto K = gh::extrinsic_curvature(geometry.time_normal, pi, phi);
  BoundaryConditions::Bjorhus::detail::add_physical_terms_to_dt_v_minus(
      make_not_null(&weyl), Scalar<DataVector>(n, 0.), geometry.normal_covector,
      geometry.normal_vector, geometry.time_normal, screen.lower, screen.mixed,
      screen.upper, geometry.inverse_spatial_metric, K, metric,
      geometry.inverse_metric, tnsr::iaa<DataVector, 3>(n, 0.),
      tnsr::aa<DataVector, 3>(n, 0.), phi, d_phi, d_pi,
      std::array<DataVector, 4>{DataVector(n, 0.), DataVector(n, 0.),
                                DataVector(n, 0.), DataVector(n, 1.)});
  DataVector speed = -get(geometry.lapse);
  for (size_t i = 0; i < 3; ++i) {
    speed -= geometry.shift.get(i) * geometry.normal_covector.get(i);
    if (mesh_velocity.has_value()) {
      speed -= mesh_velocity->get(i) * geometry.normal_covector.get(i);
    }
  }
  if (max(speed) >= 0. or min(speed + 2. * get(geometry.lapse)) <= 0.) {
    ERROR("Physical trace requires lambda_minus<0<lambda_plus.");
  }
  const auto rhs = physical_trace_rhs(
      pi, phi, dt_pi, dt_phi,
      get<gh::Tags::ThreeIndexConstraint<DataVector, 3>>(temp), weyl,
      target_weyl, gamma2, speed, geometry.normal_vector,
      geometry.dt_normal_vector, geometry.polarization,
      geometry.dt_polarization);
  DataVector result(2 * n);
  std::copy(rhs[0].begin(), rhs[0].end(), result.begin());
  std::copy(rhs[1].begin(), rhs[1].end(), result.data() + n);
  return result;
}
}  // namespace gh::worldtube
