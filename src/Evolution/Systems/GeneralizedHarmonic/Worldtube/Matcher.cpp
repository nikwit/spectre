// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::Worldtube {
namespace {

// Geometry of the data on the sphere that is independent of the model
// parameters: the outward normal, the null frame and transverse projector of
// the *numerical* metric, and gamma_2. Mirrors null_frame/decompose of the
// offline analysis (pinned_gauge_sector.py) exactly.
struct SphereFrame {
  size_t n_points;
  tnsr::AA<DataVector, 3> inverse_metric;
  DataVector lapse;
  std::array<DataVector, 3> shift;
  std::array<DataVector, 3> normal_up;  // outward radial, unit
  std::array<DataVector, 4> k_low, l_low, k_up, l_up;
  // transverse projector with mixed indices, P^c_b
  std::array<std::array<DataVector, 4>, 4> proj_ud;
  DataVector gamma2;
};

SphereFrame build_frame(const tnsr::aa<DataVector, 3>& metric,
                        const Scalar<DataVector>& gamma2,
                        const tnsr::I<DataVector, 3>& coords,
                        const std::array<double, 3>& center) {
  SphereFrame fr{};
  const size_t n = get<0, 0>(metric).size();
  fr.n_points = n;
  fr.inverse_metric = determinant_and_inverse(metric).second;
  const auto& inv = fr.inverse_metric;

  const DataVector lapse_sq = -1. / get<0, 0>(inv);
  fr.lapse = sqrt(lapse_sq);
  std::array<std::array<DataVector, 3>, 3> spatial_inv{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(fr.shift, i) = lapse_sq * inv.get(0, i + 1);
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(gsl::at(spatial_inv, i), j) =
          inv.get(i + 1, j + 1) +
          inv.get(0, i + 1) * inv.get(0, j + 1) * lapse_sq;
    }
  }

  // outward radial unit normal from the coordinate direction
  std::array<DataVector, 3> direction{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(direction, i) = coords.get(i) - gsl::at(center, i);
  }
  DataVector norm_sq(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      norm_sq += gsl::at(direction, i) * gsl::at(gsl::at(spatial_inv, i), j) *
                 gsl::at(direction, j);
    }
  }
  const DataVector norm = sqrt(norm_sq);
  std::array<DataVector, 3> normal_low{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(normal_low, i) = gsl::at(direction, i) / norm;
  }
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(fr.normal_up, i) = DataVector(n, 0.);
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(fr.normal_up, i) +=
          gsl::at(gsl::at(spatial_inv, i), j) * gsl::at(normal_low, j);
    }
  }

  // s^a: domain-outward spatial normal at an inner boundary points inward,
  // s^i = -normal_up; n_a = (-alpha, 0); k = (n - s)/sqrt2, l = (n + s)/sqrt2
  std::array<DataVector, 4> n_low{};
  n_low[0] = -fr.lapse;
  for (size_t i = 1; i < 4; ++i) {
    gsl::at(n_low, i) = DataVector(n, 0.);
  }
  std::array<DataVector, 4> s_up{};
  s_up[0] = DataVector(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(s_up, i + 1) = -gsl::at(fr.normal_up, i);
  }
  std::array<DataVector, 4> s_low{};
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(s_low, a) = DataVector(n, 0.);
    for (size_t b = 0; b < 4; ++b) {
      gsl::at(s_low, a) += metric.get(a, b) * gsl::at(s_up, b);
    }
  }
  const double root_half = 1. / std::sqrt(2.);
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(fr.k_low, a) = root_half * (gsl::at(n_low, a) - gsl::at(s_low, a));
    gsl::at(fr.l_low, a) = root_half * (gsl::at(n_low, a) + gsl::at(s_low, a));
  }
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(fr.k_up, a) = DataVector(n, 0.);
    gsl::at(fr.l_up, a) = DataVector(n, 0.);
    for (size_t b = 0; b < 4; ++b) {
      gsl::at(fr.k_up, a) += inv.get(a, b) * gsl::at(fr.k_low, b);
      gsl::at(fr.l_up, a) += inv.get(a, b) * gsl::at(fr.l_low, b);
    }
  }
  // P_ab = g_ab + k_a l_b + l_a k_b;  P^c_b = g^{ca} P_ab
  for (size_t c = 0; c < 4; ++c) {
    for (size_t b = 0; b < 4; ++b) {
      DataVector& entry = gsl::at(gsl::at(fr.proj_ud, c), b);
      entry = DataVector(n, 0.);
      for (size_t a = 0; a < 4; ++a) {
        entry += inv.get(c, a) *
                 (metric.get(a, b) +
                  gsl::at(fr.k_low, a) * gsl::at(fr.l_low, b) +
                  gsl::at(fr.l_low, a) * gsl::at(fr.k_low, b));
      }
    }
  }
  fr.gamma2 = get(gamma2);
  return fr;
}

// u^-_ab = Pi_ab + n_out^k Phi_kab - gamma_2 g_ab
tnsr::aa<DataVector, 3> u_minus_of(const tnsr::aa<DataVector, 3>& metric,
                                   const tnsr::aa<DataVector, 3>& pi,
                                   const tnsr::iaa<DataVector, 3>& phi,
                                   const SphereFrame& fr) {
  tnsr::aa<DataVector, 3> u(fr.n_points);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      u.get(a, b) = pi.get(a, b) - fr.gamma2 * metric.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        u.get(a, b) += gsl::at(fr.normal_up, k) * phi.get(k, a, b);
      }
    }
  }
  return u;
}

// gauge components: A = u(l,l), C = u(k,l), V_b = -P^c_b u_cd l^d
struct GaugeComponents {
  DataVector a, c;
  std::array<DataVector, 4> v;
};

GaugeComponents gauge_components(const tnsr::aa<DataVector, 3>& u,
                                 const SphereFrame& fr) {
  GaugeComponents out{};
  out.a = DataVector(fr.n_points, 0.);
  out.c = DataVector(fr.n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      out.a += u.get(a, b) * gsl::at(fr.l_up, a) * gsl::at(fr.l_up, b);
      out.c += u.get(a, b) * gsl::at(fr.k_up, a) * gsl::at(fr.l_up, b);
    }
  }
  for (size_t b = 0; b < 4; ++b) {
    DataVector& vb = gsl::at(out.v, b);
    vb = DataVector(fr.n_points, 0.);
    for (size_t c = 0; c < 4; ++c) {
      for (size_t d = 0; d < 4; ++d) {
        vb -= gsl::at(gsl::at(fr.proj_ud, c), b) * u.get(c, d) *
              gsl::at(fr.l_up, d);
      }
    }
  }
  return out;
}

// modal coefficients with l <= fit_l_max of the six gauge-component fields,
// concatenated
std::vector<double> gauge_modes(const GaugeComponents& gc,
                                const ylm::Spherepack& ylm_transform,
                                const std::vector<size_t>& mode_indices) {
  std::vector<double> out;
  out.reserve(6 * mode_indices.size());
  const auto append = [&out, &ylm_transform,
                       &mode_indices](const DataVector& field) {
    const DataVector spec = ylm_transform.phys_to_spec(field);
    for (const size_t idx : mode_indices) {
      out.push_back(spec[idx]);
    }
  };
  append(gc.a);
  append(gc.c);
  for (size_t b = 0; b < 4; ++b) {
    append(gsl::at(gc.v, b));
  }
  return out;
}

// x (9 free) -> p (13) with the kinematic velocity and trace pins
std::array<double, num_map_parameters> embed(
    const std::array<double, 9>& x, const MatcherConfig& config) {
  std::array<double, num_map_parameters> p{};
  p[0] = x[0];
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(p, 1 + i) = gsl::at(x, 1 + i);
    gsl::at(p, 4 + i) = (1. + x[0]) * gsl::at(config.center_velocity, i);
  }
  for (size_t i = 0; i < 5; ++i) {
    gsl::at(p, 7 + i) = gsl::at(x, 4 + i);
  }
  p[12] = 3. * config.trace_strain_pin - x[4] - x[7];
  return p;
}

std::array<double, 9> extract(const std::array<double, num_map_parameters>& p) {
  return {p[0], p[1], p[2], p[3], p[7], p[8], p[9], p[10], p[11]};
}

double norm_of(const std::vector<double>& v) {
  double out = 0.;
  for (const double x : v) {
    out += x * x;
  }
  return std::sqrt(out);
}

// solve the 9x9 normal equations by Gaussian elimination with partial
// pivoting
std::array<double, 9> solve_normal_equations(
    std::array<std::array<double, 9>, 9> a, std::array<double, 9> b) {
  for (size_t col = 0; col < 9; ++col) {
    size_t pivot = col;
    for (size_t row = col + 1; row < 9; ++row) {
      if (std::abs(gsl::at(gsl::at(a, row), col)) >
          std::abs(gsl::at(gsl::at(a, pivot), col))) {
        pivot = row;
      }
    }
    std::swap(gsl::at(a, col), gsl::at(a, pivot));
    std::swap(gsl::at(b, col), gsl::at(b, pivot));
    const double diag = gsl::at(gsl::at(a, col), col);
    for (size_t row = col + 1; row < 9; ++row) {
      const double factor = gsl::at(gsl::at(a, row), col) / diag;
      for (size_t k = col; k < 9; ++k) {
        gsl::at(gsl::at(a, row), k) -= factor * gsl::at(gsl::at(a, col), k);
      }
      gsl::at(b, row) -= factor * gsl::at(b, col);
    }
  }
  std::array<double, 9> x{};
  for (size_t row_plus_one = 9; row_plus_one > 0; --row_plus_one) {
    const size_t row = row_plus_one - 1;
    double sum = gsl::at(b, row);
    for (size_t k = row + 1; k < 9; ++k) {
      sum -= gsl::at(gsl::at(a, row), k) * gsl::at(x, k);
    }
    gsl::at(x, row) = sum / gsl::at(gsl::at(a, row), row);
  }
  return x;
}
}  // namespace

FitResult fit_map_parameters(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const std::array<double, num_map_parameters>& p_start,
    const std::array<double, num_map_parameters>& pdot_estimate) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  // modal indices with l <= fit_l_max
  std::vector<size_t> mode_indices;
  ylm::SpherepackIterator iter(ylm_transform.l_max(), ylm_transform.m_max());
  for (size_t l = 0; l <= config.fit_l_max; ++l) {
    for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
      iter.set(l, m);
      mode_indices.push_back(iter());
    }
  }

  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, config.center);
  const std::vector<double> data_modes =
      gauge_modes(gauge_components(u_minus_of(spacetime_metric, pi, phi,
                                              frame),
                                   frame),
                  ylm_transform, mode_indices);

  const auto residual_of =
      [&](const std::array<double, 9>& x) -> std::vector<double> {
    const auto p = embed(x, config);
    tnsr::aa<DataVector, 3> model_metric{};
    tnsr::aa<DataVector, 3> model_pi{};
    tnsr::iaa<DataVector, 3> model_phi{};
    gh::Solutions::affine_map_model::evolved_variables(
        make_not_null(&model_metric), make_not_null(&model_pi),
        make_not_null(&model_phi), inertial_coords, config.mass,
        config.center, p, pdot_estimate);
    std::vector<double> modes = gauge_modes(
        gauge_components(u_minus_of(model_metric, model_pi, model_phi, frame),
                         frame),
        ylm_transform, mode_indices);
    for (size_t i = 0; i < modes.size(); ++i) {
      modes[i] -= data_modes[i];
    }
    return modes;
  };

  FitResult result{};
  std::array<double, 9> x = extract(p_start);
  std::vector<double> residual = residual_of(x);
  result.residual_initial = norm_of(residual);

  constexpr size_t max_iterations = 8;
  constexpr double fd_step = 1.0e-8;
  const size_t n_res = residual.size();
  for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
    // finite-difference Jacobian
    std::array<std::vector<double>, 9> jac{};
    for (size_t a = 0; a < 9; ++a) {
      auto x_plus = x;
      gsl::at(x_plus, a) += fd_step;
      gsl::at(jac, a) = residual_of(x_plus);
      for (size_t i = 0; i < n_res; ++i) {
        gsl::at(jac, a)[i] = (gsl::at(jac, a)[i] - residual[i]) / fd_step;
      }
    }
    std::array<std::array<double, 9>, 9> jtj{};
    std::array<double, 9> jtr{};
    for (size_t a = 0; a < 9; ++a) {
      for (size_t b = a; b < 9; ++b) {
        double sum = 0.;
        for (size_t i = 0; i < n_res; ++i) {
          sum += gsl::at(jac, a)[i] * gsl::at(jac, b)[i];
        }
        gsl::at(gsl::at(jtj, a), b) = sum;
        gsl::at(gsl::at(jtj, b), a) = sum;
      }
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        sum += gsl::at(jac, a)[i] * residual[i];
      }
      gsl::at(jtr, a) = -sum;
    }
    const auto delta = solve_normal_equations(jtj, jtr);
    double max_delta = 0.;
    double max_x = 0.;
    for (size_t a = 0; a < 9; ++a) {
      gsl::at(x, a) += gsl::at(delta, a);
      max_delta = std::max(max_delta, std::abs(gsl::at(delta, a)));
      max_x = std::max(max_x, std::abs(gsl::at(x, a)));
    }
    residual = residual_of(x);
    result.iterations = iteration + 1;
    if (max_delta < std::max(1.0e-13, 1.0e-8 * max_x)) {
      break;
    }
  }
  result.residual_final = norm_of(residual);
  result.p = embed(x, config);
  return result;
}
}  // namespace gh::Worldtube
