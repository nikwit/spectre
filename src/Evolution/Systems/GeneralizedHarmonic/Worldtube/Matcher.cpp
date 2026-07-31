// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
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
// Kept spherical-harmonic modes with weights that make the least-squares
// objective rotation-invariant: Spherepack coefficients of unit-L2
// harmonics are sqrt(2) smaller for m != 0 than for m = 0, so uniform
// weighting of the raw coefficients over-weights the polar axis.
struct ModeSet {
  std::vector<size_t> indices;
  std::vector<double> weights;
};

ModeSet kept_modes(const ylm::Spherepack& ylm_transform,
                   const size_t fit_l_max) {
  ModeSet out;
  ylm::SpherepackIterator iter(ylm_transform.l_max(), ylm_transform.m_max());
  for (size_t l = 0; l <= fit_l_max; ++l) {
    for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
      iter.set(l, m);
      out.indices.push_back(iter());
      out.weights.push_back(m == 0 ? 1.0 : M_SQRT2);
    }
  }
  return out;
}

std::vector<double> gauge_modes(const GaugeComponents& gc,
                                const ylm::Spherepack& ylm_transform,
                                const ModeSet& modes) {
  std::vector<double> out;
  out.reserve(6 * modes.indices.size());
  const auto append = [&out, &ylm_transform, &modes](const DataVector& field) {
    const DataVector spec = ylm_transform.phys_to_spec(field);
    for (size_t k = 0; k < modes.indices.size(); ++k) {
      out.push_back(modes.weights[k] * spec[modes.indices[k]]);
    }
  };
  append(gc.a);
  append(gc.c);
  for (size_t b = 0; b < 4; ++b) {
    append(gsl::at(gc.v, b));
  }
  return out;
}

// x (9 or 12 free) -> (p (13), center offset q) with the kinematic velocity
// and trace pins; x[9..11] = q when the zeroth-order offset is fitted. The
// time offset q^0 is an exact zero mode of the static-in-time model and is
// never a parameter.
void embed(const gsl::not_null<std::array<double, num_map_parameters>*> p,
           const gsl::not_null<std::array<double, 3>*> center_offset,
           const std::vector<double>& x, const MatcherConfig& config) {
  (*p)[0] = x[0];
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(*p, 1 + i) = x[1 + i];
    gsl::at(*p, 4 + i) = (1. + x[0]) * gsl::at(config.center_velocity, i);
    gsl::at(*center_offset, i) = config.fit_center_offset ? x[9 + i] : 0.;
  }
  for (size_t i = 0; i < 5; ++i) {
    gsl::at(*p, 7 + i) = x[4 + i];
  }
  (*p)[12] = 3. * config.trace_strain_pin - x[4] - x[7];
}

double norm_of(const std::vector<double>& v) {
  double out = 0.;
  for (const double x : v) {
    out += x * x;
  }
  return std::sqrt(out);
}

// solve the n x n normal equations by Gaussian elimination with partial
// pivoting
std::vector<double> solve_normal_equations(std::vector<std::vector<double>> a,
                                           std::vector<double> b) {
  const size_t n = b.size();
  for (size_t col = 0; col < n; ++col) {
    size_t pivot = col;
    for (size_t row = col + 1; row < n; ++row) {
      if (std::abs(a[row][col]) > std::abs(a[pivot][col])) {
        pivot = row;
      }
    }
    std::swap(a[col], a[pivot]);
    std::swap(b[col], b[pivot]);
    const double diag = a[col][col];
    for (size_t row = col + 1; row < n; ++row) {
      const double factor = a[row][col] / diag;
      for (size_t k = col; k < n; ++k) {
        a[row][k] -= factor * a[col][k];
      }
      b[row] -= factor * b[col];
    }
  }
  std::vector<double> x(n, 0.);
  for (size_t row_plus_one = n; row_plus_one > 0; --row_plus_one) {
    const size_t row = row_plus_one - 1;
    double sum = b[row];
    for (size_t k = row + 1; k < n; ++k) {
      sum -= a[row][k] * x[k];
    }
    x[row] = sum / a[row][row];
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
    const std::array<double, 3>& center_offset_start,
    const std::array<double, num_map_parameters>& pdot_estimate) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  const ModeSet modes = kept_modes(ylm_transform, config.fit_l_max);

  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, config.center);
  const std::vector<double> data_modes =
      gauge_modes(gauge_components(u_minus_of(spacetime_metric, pi, phi,
                                              frame),
                                   frame),
                  ylm_transform, modes);

  const size_t n_free = config.fit_center_offset ? 12 : 9;
  const auto residual_of =
      [&](const std::vector<double>& x) -> std::vector<double> {
    std::array<double, num_map_parameters> p{};
    std::array<double, 3> center_offset{};
    embed(make_not_null(&p), make_not_null(&center_offset), x, config);
    std::array<double, 3> model_center = config.center;
    for (size_t i = 0; i < 3; ++i) {
      gsl::at(model_center, i) += gsl::at(center_offset, i);
    }
    tnsr::aa<DataVector, 3> model_metric{};
    tnsr::aa<DataVector, 3> model_pi{};
    tnsr::iaa<DataVector, 3> model_phi{};
    gh::Solutions::affine_map_model::evolved_variables(
        make_not_null(&model_metric), make_not_null(&model_pi),
        make_not_null(&model_phi), inertial_coords, config.mass,
        model_center, p, pdot_estimate);
    std::vector<double> modes = gauge_modes(
        gauge_components(u_minus_of(model_metric, model_pi, model_phi, frame),
                         frame),
        ylm_transform, modes);
    for (size_t i = 0; i < modes.size(); ++i) {
      modes[i] -= data_modes[i];
    }
    return modes;
  };

  FitResult result{};
  std::vector<double> x(n_free, 0.);
  x[0] = p_start[0];
  for (size_t i = 0; i < 3; ++i) {
    x[1 + i] = gsl::at(p_start, 1 + i);
  }
  for (size_t i = 0; i < 5; ++i) {
    x[4 + i] = gsl::at(p_start, 7 + i);
  }
  if (config.fit_center_offset) {
    for (size_t i = 0; i < 3; ++i) {
      x[9 + i] = gsl::at(center_offset_start, i);
    }
  }
  std::vector<double> residual = residual_of(x);
  result.residual_initial = norm_of(residual);

  constexpr size_t max_iterations = 8;
  constexpr double fd_step = 1.0e-8;
  const size_t n_res = residual.size();
  for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
    // finite-difference Jacobian
    std::vector<std::vector<double>> jac(n_free);
    for (size_t a = 0; a < n_free; ++a) {
      auto x_plus = x;
      x_plus[a] += fd_step;
      jac[a] = residual_of(x_plus);
      for (size_t i = 0; i < n_res; ++i) {
        jac[a][i] = (jac[a][i] - residual[i]) / fd_step;
      }
    }
    std::vector<std::vector<double>> jtj(n_free,
                                         std::vector<double>(n_free, 0.));
    std::vector<double> jtr(n_free, 0.);
    for (size_t a = 0; a < n_free; ++a) {
      for (size_t b = a; b < n_free; ++b) {
        double sum = 0.;
        for (size_t i = 0; i < n_res; ++i) {
          sum += jac[a][i] * jac[b][i];
        }
        jtj[a][b] = sum;
        jtj[b][a] = sum;
      }
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        sum += jac[a][i] * residual[i];
      }
      jtr[a] = -sum;
    }
    const auto delta = solve_normal_equations(std::move(jtj), std::move(jtr));
    double max_delta = 0.;
    double max_x = 0.;
    for (size_t a = 0; a < n_free; ++a) {
      x[a] += delta[a];
      max_delta = std::max(max_delta, std::abs(delta[a]));
      max_x = std::max(max_x, std::abs(x[a]));
    }
    residual = residual_of(x);
    result.iterations = iteration + 1;
    if (max_delta < std::max(1.0e-13, 1.0e-8 * max_x)) {
      break;
    }
  }
  result.residual_final = norm_of(residual);
  embed(make_not_null(&result.p), make_not_null(&result.center_offset), x,
        config);
  return result;
}

namespace {
// Shared core of the rate and acceleration fits: project a symmetric-tensor
// target onto the nine unpinned rate directions of the covariant response,
// -(g R_A g), and unfold the (rate-level) pins.
RateFitResult project_onto_rate_directions(
    const tnsr::aa<DataVector, 3>& target,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config) {
  const size_t n_points = get<0, 0>(spacetime_metric).size();
  const ModeSet modes = kept_modes(ylm_transform, config.fit_l_max);

  // weighted modes of the ten components of a symmetric rank-2 tensor
  const auto tensor_modes =
      [&ylm_transform, &modes](const tnsr::aa<DataVector, 3>& tensor) {
        std::vector<double> out;
        out.reserve(10 * modes.indices.size());
        for (size_t a = 0; a < 4; ++a) {
          for (size_t b = a; b < 4; ++b) {
            const DataVector spec =
                ylm_transform.phys_to_spec(tensor.get(a, b));
            for (size_t k = 0; k < modes.indices.size(); ++k) {
              out.push_back(modes.weights[k] * spec[modes.indices[k]]);
            }
          }
        }
        return out;
      };
  const std::vector<double> target_modes = tensor_modes(target);

  // rate-response columns for the nine unpinned rate directions: the rate
  // embedding is the derivative of the value pins, so the trace pin drops
  // out (constant) and the velocity pin loses its affine "1 + " part
  std::array<DataVector, 3> y{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(y, i) = inertial_coords.get(i) - gsl::at(config.center, i);
  }
  constexpr size_t n_free = 9;
  std::vector<std::vector<double>> columns(n_free);
  for (size_t a = 0; a < n_free; ++a) {
    std::array<double, num_map_parameters> rate_direction{};
    if (a == 0) {
      rate_direction[0] = 1.;
      for (size_t i = 0; i < 3; ++i) {
        gsl::at(rate_direction, 4 + i) = gsl::at(config.center_velocity, i);
      }
    } else if (a < 4) {
      gsl::at(rate_direction, a) = 1.;
    } else {
      gsl::at(rate_direction, 3 + a) = 1.;
      if (a == 4 or a == 7) {
        rate_direction[12] = -1.;
      }
    }
    tnsr::AA<DataVector, 3> response(n_points);
    gh::Solutions::affine_map_model::inverse_metric_combination(
        make_not_null(&response), y, config.mass, 0., rate_direction);
    tnsr::aa<DataVector, 3> column(n_points, 0.);
    for (size_t c = 0; c < 4; ++c) {
      for (size_t d = c; d < 4; ++d) {
        for (size_t e = 0; e < 4; ++e) {
          for (size_t f = 0; f < 4; ++f) {
            column.get(c, d) -= spacetime_metric.get(c, e) *
                                response.get(e, f) *
                                spacetime_metric.get(f, d);
          }
        }
      }
    }
    columns[a] = tensor_modes(column);
  }

  // linear least squares via the normal equations
  const size_t n_res = target_modes.size();
  std::vector<std::vector<double>> jtj(n_free, std::vector<double>(n_free, 0.));
  std::vector<double> jtr(n_free, 0.);
  for (size_t a = 0; a < n_free; ++a) {
    for (size_t b = a; b < n_free; ++b) {
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        sum += columns[a][i] * columns[b][i];
      }
      jtj[a][b] = sum;
      jtj[b][a] = sum;
    }
    double sum = 0.;
    for (size_t i = 0; i < n_res; ++i) {
      sum += columns[a][i] * target_modes[i];
    }
    jtr[a] = sum;
  }
  const auto x = solve_normal_equations(std::move(jtj), std::move(jtr));

  RateFitResult result{};
  result.residual_initial = norm_of(target_modes);
  std::vector<double> residual = target_modes;
  for (size_t a = 0; a < n_free; ++a) {
    for (size_t i = 0; i < n_res; ++i) {
      residual[i] -= x[a] * columns[a][i];
    }
  }
  result.residual_final = norm_of(residual);
  result.pdot[0] = x[0];
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(result.pdot, 1 + i) = x[1 + i];
    gsl::at(result.pdot, 4 + i) = x[0] * gsl::at(config.center_velocity, i);
  }
  for (size_t i = 0; i < 5; ++i) {
    gsl::at(result.pdot, 7 + i) = x[4 + i];
  }
  result.pdot[12] = -x[4] - x[7];
  return result;
}
}  // namespace

RateFitResult fit_map_parameter_rates(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  const size_t n_points = get<0, 0>(spacetime_metric).size();
  // data-side dt g = beta^k Phi_k - alpha Pi (the kinematic identity)
  const auto inverse_metric = determinant_and_inverse(spacetime_metric).second;
  const DataVector lapse_sq = -1. / get<0, 0>(inverse_metric);
  const DataVector lapse = sqrt(lapse_sq);
  std::array<DataVector, 3> shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(shift, i) = lapse_sq * inverse_metric.get(0, i + 1);
  }
  tnsr::aa<DataVector, 3> dt_metric(n_points);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_metric.get(a, b) = -lapse * pi.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        dt_metric.get(a, b) += gsl::at(shift, k) * phi.get(k, a, b);
      }
    }
  }
  return project_onto_rate_directions(dt_metric, spacetime_metric,
                                      inertial_coords, ylm_transform, config);
}

RateFitResult fit_map_parameter_accelerations(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi,
    const std::array<double, num_map_parameters>& pdot_state,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  const size_t n_points = get<0, 0>(spacetime_metric).size();
  // lapse and shift and their time derivatives from the metric and its dt:
  // dt G^{ab} = -(G dt_g G)^{ab};  alpha = (-G^tt)^{-1/2};
  // dt alpha = (alpha^3/2) dt G^tt;  beta^i = alpha^2 G^ti.
  const auto inverse_metric = determinant_and_inverse(spacetime_metric).second;
  const DataVector lapse_sq = -1. / get<0, 0>(inverse_metric);
  const DataVector lapse = sqrt(lapse_sq);
  std::array<DataVector, 3> shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(shift, i) = lapse_sq * inverse_metric.get(0, i + 1);
  }
  tnsr::AA<DataVector, 3> dt_inverse_metric(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          dt_inverse_metric.get(a, b) -= inverse_metric.get(a, c) *
                                         dt_spacetime_metric.get(c, d) *
                                         inverse_metric.get(d, b);
        }
      }
    }
  }
  const DataVector dt_lapse =
      0.5 * lapse * lapse_sq * get<0, 0>(dt_inverse_metric);
  std::array<DataVector, 3> dt_shift{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(dt_shift, i) =
        2. * lapse * dt_lapse * inverse_metric.get(0, i + 1) +
        lapse_sq * dt_inverse_metric.get(0, i + 1);
  }
  // d2t g = dt beta . Phi + beta . dt Phi - dt alpha . Pi - alpha . dt Pi
  tnsr::aa<DataVector, 3> d2t_metric(n_points);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      d2t_metric.get(a, b) =
          -dt_lapse * pi.get(a, b) - lapse * dt_pi.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        d2t_metric.get(a, b) += gsl::at(dt_shift, k) * phi.get(k, a, b) +
                                gsl::at(shift, k) * dt_phi.get(k, a, b);
      }
    }
  }
  // Subtract the pdot-quadratic Hessian term of the model's second time
  // derivative: within the amplitude-linear map, g(p) = (G_Schw^{-1} +
  // sum_A p_A R_A)^{-1} gives exactly
  //   d2t g = -(g (sum_A pddot_A R_A) g) + 2 g S g S g,
  // with S = sum_A pdot_A R_A. Moving the S-term to the data side keeps the
  // solve linear in pddot and the truncation consistent for a second-order
  // ODE (previously it was dropped as quasi-static).
  std::array<DataVector, 3> y{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(y, i) = inertial_coords.get(i) - gsl::at(config.center, i);
  }
  tnsr::AA<DataVector, 3> s_response(n_points);
  gh::Solutions::affine_map_model::inverse_metric_combination(
      make_not_null(&s_response), y, config.mass, 0., pdot_state);
  tnsr::aa<DataVector, 3> g_s_g(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          g_s_g.get(a, b) += spacetime_metric.get(a, c) *
                             s_response.get(c, d) * spacetime_metric.get(d, b);
        }
      }
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          d2t_metric.get(a, b) -= 2. * g_s_g.get(a, c) *
                                  s_response.get(c, d) *
                                  spacetime_metric.get(d, b);
        }
      }
    }
  }
  return project_onto_rate_directions(d2t_metric, spacetime_metric,
                                      inertial_coords, ylm_transform, config);
}
}  // namespace gh::Worldtube
