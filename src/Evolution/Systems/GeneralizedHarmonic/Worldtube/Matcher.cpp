// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
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
        entry += inv.get(c, a) * (metric.get(a, b) +
                                  gsl::at(fr.k_low, a) * gsl::at(fr.l_low, b) +
                                  gsl::at(fr.l_low, a) * gsl::at(fr.k_low, b));
      }
    }
  }
  fr.gamma2 = get(gamma2);
  return fr;
}

// u^\pm_ab = Pi_ab \pm n_out^k Phi_kab - gamma_2 g_ab; normal_sign +1
// gives the incoming characteristic the ghost BC sets, -1 the outgoing
// one (into the excision), which the BC never touches.
tnsr::aa<DataVector, 3> u_minus_of(const tnsr::aa<DataVector, 3>& metric,
                                   const tnsr::aa<DataVector, 3>& pi,
                                   const tnsr::iaa<DataVector, 3>& phi,
                                   const SphereFrame& fr,
                                   const double normal_sign) {
  tnsr::aa<DataVector, 3> u(fr.n_points);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      u.get(a, b) = pi.get(a, b) - fr.gamma2 * metric.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        u.get(a, b) +=
            normal_sign * gsl::at(fr.normal_up, k) * phi.get(k, a, b);
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
  std::vector<size_t> ells;
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
      out.ells.push_back(l);
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

// The free vector, in order:
//   x[0]                    qdot^0
//   x[1..3]                 beta_i
//   x[4 .. 3+n_strain]      sigma (5 pinned-trace components, or all 6)
//   next 3, if FitVelocity  qdot^i     (else pinned to v_centre)
//   next 3, if FitCenterOffset  q^i
// so 9 free at the leanest and 16 with everything on. The time offset q^0 is
// an exact zero mode of the static-in-time model and is never a parameter.
void embed(const gsl::not_null<std::array<double, num_map_parameters>*> p,
           const gsl::not_null<std::array<double, 3>*> center_offset,
           const std::vector<double>& x, const MatcherConfig& config,
           const double trace_pin,
           const std::array<double, 3>& fixed_center_offset) {
  (*p)[0] = x[0];
  const size_t n_strain = config.fit_trace_strain ? 6 : 5;
  const size_t velocity_start = 4 + n_strain;
  const size_t offset_start = velocity_start + (config.fit_velocity ? 3 : 0);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(*p, 1 + i) = x[1 + i];
    gsl::at(*p, 4 + i) = config.fit_velocity
                             ? x[velocity_start + i]
                             : gsl::at(config.center_velocity, i);
    gsl::at(*center_offset, i) = config.fit_center_offset
                                     ? x[offset_start + i]
                                     : gsl::at(fixed_center_offset, i);
  }
  for (size_t i = 0; i < n_strain; ++i) {
    gsl::at(*p, 7 + i) = x[4 + i];
  }
  if (not config.fit_trace_strain) {
    (*p)[12] = 3. * trace_pin - x[4] - x[7];
  }
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

// The squared 2-norm condition number of a design matrix is the condition
// number of its symmetric normal matrix. Find the latter's eigenvalues with a
// small Jacobi iteration; the matcher matrices have at most nine columns.
double design_matrix_condition_number(
    std::vector<std::vector<double>> normal_matrix) {
  const size_t n = normal_matrix.size();
  ASSERT(n > 0, "Cannot condition an empty design matrix.");
  for (size_t sweep = 0; sweep < 100 * n * n; ++sweep) {
    size_t p = 0;
    size_t q = 0;
    double largest_off_diagonal = 0.;
    double largest_diagonal = 0.;
    for (size_t i = 0; i < n; ++i) {
      largest_diagonal =
          std::max(largest_diagonal, std::abs(normal_matrix[i][i]));
      for (size_t j = i + 1; j < n; ++j) {
        if (std::abs(normal_matrix[i][j]) > largest_off_diagonal) {
          largest_off_diagonal = std::abs(normal_matrix[i][j]);
          p = i;
          q = j;
        }
      }
    }
    if (largest_off_diagonal <=
        10. * std::numeric_limits<double>::epsilon() * largest_diagonal) {
      break;
    }
    const double app = normal_matrix[p][p];
    const double aqq = normal_matrix[q][q];
    const double apq = normal_matrix[p][q];
    const double tau = (aqq - app) / (2. * apq);
    const double tangent =
        std::copysign(1. / (std::abs(tau) + std::sqrt(1. + tau * tau)), tau);
    const double cosine = 1. / std::sqrt(1. + tangent * tangent);
    const double sine = tangent * cosine;
    for (size_t k = 0; k < n; ++k) {
      if (k == p or k == q) {
        continue;
      }
      const double akp = normal_matrix[k][p];
      const double akq = normal_matrix[k][q];
      normal_matrix[k][p] = cosine * akp - sine * akq;
      normal_matrix[p][k] = normal_matrix[k][p];
      normal_matrix[k][q] = sine * akp + cosine * akq;
      normal_matrix[q][k] = normal_matrix[k][q];
    }
    normal_matrix[p][p] = app - tangent * apq;
    normal_matrix[q][q] = aqq + tangent * apq;
    normal_matrix[p][q] = 0.;
    normal_matrix[q][p] = 0.;
  }
  double smallest = std::numeric_limits<double>::infinity();
  double largest = 0.;
  for (size_t i = 0; i < n; ++i) {
    smallest = std::min(smallest, normal_matrix[i][i]);
    largest = std::max(largest, normal_matrix[i][i]);
  }
  if (smallest <= 0.) {
    return std::numeric_limits<double>::infinity();
  }
  return std::sqrt(largest / smallest);
}
}  // namespace

namespace detail {
WeightedModeFit fit_weighted_modes(
    const std::vector<double>& target,
    const std::vector<double>& held_out_target,
    const std::vector<std::vector<double>>& columns,
    const std::vector<size_t>& row_blocks,
    const std::array<double, 15>& block_weights) {
  ASSERT(not columns.empty(), "The weighted mode fit needs columns.");
  ASSERT(target.size() == held_out_target.size() and
             target.size() == row_blocks.size(),
         "Fitted target, held-out target, and row-block labels must have the "
         "same number of rows.");
  for (const size_t block : row_blocks) {
    ASSERT(block < block_weights.size(),
           "Every modal row must refer to a valid block weight.");
  }
  for (const auto& column : columns) {
    ASSERT(column.size() == target.size(),
           "Every fit column must have one entry per target row.");
  }
  const size_t n_free = columns.size();
  const size_t n_res = target.size();
  std::vector<std::vector<double>> jtj(n_free, std::vector<double>(n_free, 0.));
  std::vector<double> jtr(n_free, 0.);
  for (size_t a = 0; a < n_free; ++a) {
    for (size_t b = a; b < n_free; ++b) {
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        const double weight = block_weights[row_blocks[i]];
        sum += square(weight) * columns[a][i] * columns[b][i];
      }
      jtj[a][b] = sum;
      jtj[b][a] = sum;
    }
    for (size_t i = 0; i < n_res; ++i) {
      const double weight = block_weights[row_blocks[i]];
      jtr[a] += square(weight) * columns[a][i] * target[i];
    }
  }

  WeightedModeFit result{};
  result.condition_number = design_matrix_condition_number(jtj);
  ASSERT(std::isfinite(result.condition_number),
         "UPlusBlockWeights make the acceleration fit rank deficient.");
  result.coefficients = solve_normal_equations(std::move(jtj), std::move(jtr));
  result.fitted_residual = target;
  result.held_out_residual = held_out_target;
  for (size_t a = 0; a < n_free; ++a) {
    for (size_t i = 0; i < n_res; ++i) {
      result.fitted_residual[i] -= result.coefficients[a] * columns[a][i];
      result.held_out_residual[i] -= result.coefficients[a] * columns[a][i];
    }
  }
  return result;
}
}  // namespace detail

FitResult fit_map_parameters(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const std::array<double, num_map_parameters>& p_start,
    const std::array<double, 3>& center_offset_start, const double normal_sign,
    const double trace_pin) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  const ModeSet modes = kept_modes(ylm_transform, config.fit_l_max);

  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, config.center);
  const std::vector<double> data_modes = gauge_modes(
      gauge_components(
          u_minus_of(spacetime_metric, pi, phi, frame, normal_sign), frame),
      ylm_transform, modes);
  const std::vector<double> minus_data_modes = gauge_modes(
      gauge_components(
          u_minus_of(spacetime_metric, pi, phi, frame, -normal_sign), frame),
      ylm_transform, modes);

  const size_t n_strain = config.fit_trace_strain ? 6 : 5;
  const size_t n_base = 4 + n_strain + (config.fit_velocity ? 3 : 0);
  const size_t n_free = config.fit_center_offset ? n_base + 3 : n_base;
  const auto model_at =
      [&](const gsl::not_null<tnsr::aa<DataVector, 3>*> model_metric,
          const gsl::not_null<tnsr::aa<DataVector, 3>*> model_pi,
          const gsl::not_null<tnsr::iaa<DataVector, 3>*> model_phi,
          const std::vector<double>& x) {
        std::array<double, num_map_parameters> p{};
        std::array<double, 3> center_offset{};
        embed(make_not_null(&p), make_not_null(&center_offset), x, config,
              trace_pin, center_offset_start);
        std::array<double, 3> model_center = config.center;
        for (size_t i = 0; i < 3; ++i) {
          gsl::at(model_center, i) += gsl::at(center_offset, i);
        }
        gh::Solutions::affine_map_model::first_order_evolved_variables(
            model_metric, model_pi, model_phi, inertial_coords, config.mass,
            model_center, p, config.centre_advection);
      };
  const auto residual_of =
      [&](const std::vector<double>& x, const double sign,
          const std::vector<double>& data) -> std::vector<double> {
    tnsr::aa<DataVector, 3> model_metric{};
    tnsr::aa<DataVector, 3> model_pi{};
    tnsr::iaa<DataVector, 3> model_phi{};
    model_at(make_not_null(&model_metric), make_not_null(&model_pi),
             make_not_null(&model_phi), x);
    std::vector<double> model_modes = gauge_modes(
        gauge_components(
            u_minus_of(model_metric, model_pi, model_phi, frame, sign), frame),
        ylm_transform, modes);
    for (size_t i = 0; i < model_modes.size(); ++i) {
      model_modes[i] -= data[i];
    }
    return model_modes;
  };

  const size_t n_modes = modes.indices.size();
  const auto row_class = [&n_modes](const size_t row) -> size_t {
    const size_t field = row / n_modes;  // A, C, V_0..V_3
    return field == 0 ? 0 : (field == 1 ? 1 : 2);
  };
  std::vector<size_t> row_blocks(data_modes.size());
  for (size_t i = 0; i < row_blocks.size(); ++i) {
    row_blocks[i] = 5 * row_class(i) + modes.ells[i % n_modes];
  }

  FitResult result{};
  std::vector<double> x(n_free, 0.);
  x[0] = p_start[0];
  for (size_t i = 0; i < 3; ++i) {
    x[1 + i] = gsl::at(p_start, 1 + i);
  }
  for (size_t i = 0; i < n_strain; ++i) {
    x[4 + i] = gsl::at(p_start, 7 + i);
  }
  if (config.fit_velocity) {
    for (size_t i = 0; i < 3; ++i) {
      x[4 + n_strain + i] = gsl::at(p_start, 4 + i);
    }
  }
  if (config.fit_center_offset) {
    for (size_t i = 0; i < 3; ++i) {
      x[n_base + i] = gsl::at(center_offset_start, i);
    }
  }
  const std::vector<double> baseline_x(n_free, 0.);
  const std::vector<double> baseline_residual =
      residual_of(baseline_x, normal_sign, data_modes);
  const std::vector<double> minus_baseline_residual =
      residual_of(baseline_x, -normal_sign, minus_data_modes);
  std::vector<double> residual = residual_of(x, normal_sign, data_modes);
  result.residual_initial = norm_of(residual);
  result.baseline_residual = norm_of(baseline_residual);
  result.minus_baseline_residual = norm_of(minus_baseline_residual);

  constexpr size_t max_iterations = 8;
  constexpr double fd_step = 1.0e-8;
  const size_t n_res = residual.size();
  for (size_t iteration = 0; iteration < max_iterations; ++iteration) {
    // finite-difference Jacobian
    std::vector<std::vector<double>> jac(n_free);
    for (size_t a = 0; a < n_free; ++a) {
      auto x_plus = x;
      x_plus[a] += fd_step;
      jac[a] = residual_of(x_plus, normal_sign, data_modes);
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
          const double weight = config.uplus_block_weights[row_blocks[i]];
          sum += square(weight) * jac[a][i] * jac[b][i];
        }
        jtj[a][b] = sum;
        jtj[b][a] = sum;
      }
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        const double weight = config.uplus_block_weights[row_blocks[i]];
        sum += square(weight) * jac[a][i] * residual[i];
      }
      jtr[a] = -sum;
    }
    result.condition_number = design_matrix_condition_number(jtj);
    ASSERT(std::isfinite(result.condition_number),
           "UPlusBlockWeights make the value fit rank deficient.");
    const auto delta = solve_normal_equations(std::move(jtj), std::move(jtr));
    double max_delta = 0.;
    double max_x = 0.;
    for (size_t a = 0; a < n_free; ++a) {
      x[a] += delta[a];
      max_delta = std::max(max_delta, std::abs(delta[a]));
      max_x = std::max(max_x, std::abs(x[a]));
    }
    residual = residual_of(x, normal_sign, data_modes);
    result.iterations = iteration + 1;
    if (max_delta < std::max(1.0e-13, 1.0e-8 * max_x)) {
      break;
    }
  }
  result.residual_final = norm_of(residual);
  const std::vector<double> minus_residual =
      residual_of(x, -normal_sign, minus_data_modes);
  result.minus_residual_final = norm_of(minus_residual);

  for (size_t c = 0; c < 3; ++c) {
    for (size_t l = 0; l <= 4; ++l) {
      double target_squared = 0.;
      double residual_squared = 0.;
      double minus_target_squared = 0.;
      double minus_residual_squared = 0.;
      size_t count = 0;
      for (size_t i = 0; i < residual.size(); ++i) {
        if (row_class(i) == c and modes.ells[i % n_modes] == l) {
          target_squared += square(baseline_residual[i]);
          residual_squared += square(residual[i]);
          minus_target_squared += square(minus_baseline_residual[i]);
          minus_residual_squared += square(minus_residual[i]);
          ++count;
        }
      }
      const size_t block = 5 * c + l;
      gsl::at(result.block_closure, block) =
          target_squared > 0. ? std::sqrt(residual_squared / target_squared)
                              : 0.;
      gsl::at(result.block_target_rms, block) =
          count > 0 ? std::sqrt(target_squared / count) : 0.;
      gsl::at(result.block_residual_rms, block) =
          count > 0 ? std::sqrt(residual_squared / count) : 0.;
      gsl::at(result.block_minus_closure, block) =
          minus_target_squared > 0.
              ? std::sqrt(minus_residual_squared / minus_target_squared)
              : 0.;
      gsl::at(result.block_minus_target_rms, block) =
          count > 0 ? std::sqrt(minus_target_squared / count) : 0.;
      gsl::at(result.block_minus_residual_rms, block) =
          count > 0 ? std::sqrt(minus_residual_squared / count) : 0.;
    }
  }

  tnsr::aa<DataVector, 3> baseline_metric{};
  tnsr::aa<DataVector, 3> baseline_pi{};
  tnsr::iaa<DataVector, 3> baseline_phi{};
  model_at(make_not_null(&baseline_metric), make_not_null(&baseline_pi),
           make_not_null(&baseline_phi), baseline_x);
  tnsr::aa<DataVector, 3> final_metric{};
  tnsr::aa<DataVector, 3> final_pi{};
  tnsr::iaa<DataVector, 3> final_phi{};
  model_at(make_not_null(&final_metric), make_not_null(&final_pi),
           make_not_null(&final_phi), x);
  const auto tensor_difference_norm = [](const auto& lhs, const auto& rhs) {
    double norm_squared = 0.;
    for (size_t storage = 0; storage < lhs.size(); ++storage) {
      const DataVector difference = lhs[storage] - rhs[storage];
      for (size_t point = 0; point < difference.size(); ++point) {
        norm_squared += square(difference[point]);
      }
    }
    return std::sqrt(norm_squared);
  };
  result.metric_baseline_residual =
      tensor_difference_norm(baseline_metric, spacetime_metric);
  result.metric_residual_final =
      tensor_difference_norm(final_metric, spacetime_metric);
  result.phi_baseline_residual = tensor_difference_norm(baseline_phi, phi);
  result.phi_residual_final = tensor_difference_norm(final_phi, phi);
  embed(make_not_null(&result.p), make_not_null(&result.center_offset), x,
        config, trace_pin, center_offset_start);
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
  const auto tensor_modes = [&ylm_transform,
                             &modes](const tnsr::aa<DataVector, 3>& tensor) {
    std::vector<double> out;
    out.reserve(10 * modes.indices.size());
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        const DataVector spec = ylm_transform.phys_to_spec(tensor.get(a, b));
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
                                response.get(e, f) * spacetime_metric.get(f, d);
          }
        }
      }
    }
    columns[a] = tensor_modes(column);
  }

  // Per-row block labels: component class (0 = TT, 1 = Ti, 2 = ij) and l.
  // Rows are ordered (a, b) major over the 10 tensor components, mode minor.
  const size_t n_modes = modes.indices.size();
  const size_t n_res = target_modes.size();
  const auto row_class = [&n_modes](const size_t row) -> size_t {
    const size_t comp = row / n_modes;  // 0..9 in (0,0),(0,1),..,(3,3) order
    if (comp == 0) {
      return 0;  // TT
    }
    return comp < 4 ? 1 : 2;  // Ti : ij
  };
  const auto row_ell = [&modes, &n_modes](const size_t row) -> size_t {
    return modes.ells[row % n_modes];
  };

  // Solve weights on top of the isotropy weights already in the rows: the
  // spatial monopole is demoted per config (q8 block analysis: it absorbs
  // unmodeled content into trace strain and clock rate).
  std::vector<double> solve_weight(n_res, 1.);
  for (size_t i = 0; i < n_res; ++i) {
    if (row_class(i) == 2 and row_ell(i) == 0) {
      solve_weight[i] = config.spatial_monopole_weight;
    }
  }

  // weighted linear least squares via the normal equations
  const auto weighted_solve = [&](const std::vector<size_t>& free_set,
                                  const std::vector<double>& target,
                                  const std::vector<double>& weight) {
    const size_t n_sub = free_set.size();
    std::vector<std::vector<double>> jtj(n_sub, std::vector<double>(n_sub, 0.));
    std::vector<double> jtr(n_sub, 0.);
    for (size_t a = 0; a < n_sub; ++a) {
      for (size_t b = a; b < n_sub; ++b) {
        double sum = 0.;
        for (size_t i = 0; i < n_res; ++i) {
          sum += square(weight[i]) * columns[free_set[a]][i] *
                 columns[free_set[b]][i];
        }
        jtj[a][b] = sum;
        jtj[b][a] = sum;
      }
      double sum = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        sum += square(weight[i]) * columns[free_set[a]][i] * target[i];
      }
      jtr[a] = sum;
    }
    return solve_normal_equations(std::move(jtj), std::move(jtr));
  };

  std::vector<size_t> all_free(n_free);
  std::iota(all_free.begin(), all_free.end(), 0);
  const auto x = weighted_solve(all_free, target_modes, solve_weight);

  RateFitResult result{};
  result.residual_initial = norm_of(target_modes);
  std::vector<double> residual = target_modes;
  for (size_t a = 0; a < n_free; ++a) {
    for (size_t i = 0; i < n_res; ++i) {
      residual[i] -= x[a] * columns[a][i];
    }
  }
  result.residual_final = norm_of(residual);

  // Held-out closure: relative residual per (class, l) block, evaluated on
  // the isotropy-weighted rows regardless of the solve weights, so demoted
  // blocks are still measured.
  for (size_t c = 0; c < 3; ++c) {
    for (size_t l = 0; l <= 4; ++l) {
      double rr = 0.;
      double tt = 0.;
      for (size_t i = 0; i < n_res; ++i) {
        if (row_class(i) == c and row_ell(i) == l) {
          rr += square(residual[i]);
          tt += square(target_modes[i]);
        }
      }
      gsl::at(result.block_closure, 5 * c + l) =
          tt > 0. ? std::sqrt(rr / tt) : 0.;
    }
  }

  // Estimator-spread diagnostics (q4/q8 audit): re-solve restricted
  // subsystems on their audit-preferred blocks with the complementary
  // parameters frozen at the global solution, and record the largest
  // parameter shift. V1 rows (TT l=1, Ti l=0) determine the boosts; C3
  // rows (TT l=0, ij l=2) determine clock rate and strain.
  const auto restricted_spread = [&](const std::vector<size_t>& params,
                                     const auto& row_in_blocks) {
    std::vector<double> reduced = target_modes;
    for (size_t a = 0; a < n_free; ++a) {
      if (std::find(params.begin(), params.end(), a) == params.end()) {
        for (size_t i = 0; i < n_res; ++i) {
          reduced[i] -= x[a] * columns[a][i];
        }
      }
    }
    std::vector<double> mask(n_res, 0.);
    for (size_t i = 0; i < n_res; ++i) {
      mask[i] = row_in_blocks(row_class(i), row_ell(i)) ? 1. : 0.;
    }
    const auto x_sub = weighted_solve(params, reduced, mask);
    double spread = 0.;
    for (size_t a = 0; a < params.size(); ++a) {
      spread = std::max(spread, std::abs(x_sub[a] - x[params[a]]));
    }
    return spread;
  };
  result.spread_vector =
      restricted_spread({1, 2, 3}, [](const size_t c, const size_t l) {
        return (c == 0 and l == 1) or (c == 1 and l == 0);
      });
  result.spread_clock_strain =
      restricted_spread({0, 4, 5, 6, 7, 8}, [](const size_t c, const size_t l) {
        return (c == 0 and l == 0) or (c == 2 and l == 2);
      });
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
      d2t_metric.get(a, b) = -dt_lapse * pi.get(a, b) - lapse * dt_pi.get(a, b);
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
          g_s_g.get(a, b) += spacetime_metric.get(a, c) * s_response.get(c, d) *
                             spacetime_metric.get(d, b);
        }
      }
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          d2t_metric.get(a, b) -= 2. * g_s_g.get(a, c) * s_response.get(c, d) *
                                  spacetime_metric.get(d, b);
        }
      }
    }
  }
  return project_onto_rate_directions(d2t_metric, spacetime_metric,
                                      inertial_coords, ylm_transform, config);
}

namespace {
// Time derivatives of the SphereFrame quantities, from dt of the data
// metric (chain rule through the identical construction as build_frame).
struct SphereFrameDot {
  tnsr::AA<DataVector, 3> dt_inverse_metric;
  DataVector dt_lapse;
  std::array<DataVector, 3> dt_normal_up;
  std::array<DataVector, 4> dt_k_up, dt_l_up;
  std::array<std::array<DataVector, 4>, 4> dt_proj_ud;
};

SphereFrameDot build_frame_dot(const tnsr::aa<DataVector, 3>& metric,
                               const tnsr::aa<DataVector, 3>& dt_metric,
                               const SphereFrame& fr,
                               const tnsr::I<DataVector, 3>& coords,
                               const std::array<double, 3>& center) {
  SphereFrameDot out{};
  const size_t n = fr.n_points;
  const auto& inv = fr.inverse_metric;
  set_number_of_grid_points(make_not_null(&out.dt_inverse_metric), n);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      out.dt_inverse_metric.get(a, b) = 0.;
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          out.dt_inverse_metric.get(a, b) -=
              inv.get(a, c) * dt_metric.get(c, d) * inv.get(d, b);
        }
      }
    }
  }
  const auto& dt_inv = out.dt_inverse_metric;
  const DataVector lapse_sq = square(fr.lapse);
  const DataVector dt_lapse_sq = square(lapse_sq) * get<0, 0>(dt_inv);
  out.dt_lapse = 0.5 * dt_lapse_sq / fr.lapse;

  std::array<std::array<DataVector, 3>, 3> spatial_inv{};
  std::array<std::array<DataVector, 3>, 3> dt_spatial_inv{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(gsl::at(spatial_inv, i), j) =
          inv.get(i + 1, j + 1) +
          inv.get(0, i + 1) * inv.get(0, j + 1) * lapse_sq;
      gsl::at(gsl::at(dt_spatial_inv, i), j) =
          dt_inv.get(i + 1, j + 1) +
          dt_inv.get(0, i + 1) * inv.get(0, j + 1) * lapse_sq +
          inv.get(0, i + 1) * dt_inv.get(0, j + 1) * lapse_sq +
          inv.get(0, i + 1) * inv.get(0, j + 1) * dt_lapse_sq;
    }
  }

  std::array<DataVector, 3> direction{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(direction, i) = coords.get(i) - gsl::at(center, i);
  }
  DataVector norm_sq(n, 0.);
  DataVector dt_norm_sq(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      norm_sq += gsl::at(direction, i) * gsl::at(gsl::at(spatial_inv, i), j) *
                 gsl::at(direction, j);
      dt_norm_sq += gsl::at(direction, i) *
                    gsl::at(gsl::at(dt_spatial_inv, i), j) *
                    gsl::at(direction, j);
    }
  }
  const DataVector norm = sqrt(norm_sq);
  const DataVector dt_norm = 0.5 * dt_norm_sq / norm;
  std::array<DataVector, 3> normal_low{};
  std::array<DataVector, 3> dt_normal_low{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(normal_low, i) = gsl::at(direction, i) / norm;
    gsl::at(dt_normal_low, i) = -gsl::at(direction, i) * dt_norm / norm_sq;
  }
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(out.dt_normal_up, i) = DataVector(n, 0.);
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(out.dt_normal_up, i) +=
          gsl::at(gsl::at(dt_spatial_inv, i), j) * gsl::at(normal_low, j) +
          gsl::at(gsl::at(spatial_inv, i), j) * gsl::at(dt_normal_low, j);
    }
  }

  std::array<DataVector, 4> s_up{};
  std::array<DataVector, 4> dt_s_up{};
  s_up[0] = DataVector(n, 0.);
  dt_s_up[0] = DataVector(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(s_up, i + 1) = -gsl::at(fr.normal_up, i);
    gsl::at(dt_s_up, i + 1) = -gsl::at(out.dt_normal_up, i);
  }
  std::array<DataVector, 4> s_low{};
  std::array<DataVector, 4> dt_s_low{};
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(s_low, a) = DataVector(n, 0.);
    gsl::at(dt_s_low, a) = DataVector(n, 0.);
    for (size_t b = 0; b < 4; ++b) {
      gsl::at(s_low, a) += metric.get(a, b) * gsl::at(s_up, b);
      gsl::at(dt_s_low, a) += dt_metric.get(a, b) * gsl::at(s_up, b) +
                              metric.get(a, b) * gsl::at(dt_s_up, b);
    }
  }
  std::array<DataVector, 4> n_low{};
  std::array<DataVector, 4> dt_n_low{};
  n_low[0] = -fr.lapse;
  dt_n_low[0] = -out.dt_lapse;
  for (size_t i = 1; i < 4; ++i) {
    gsl::at(n_low, i) = DataVector(n, 0.);
    gsl::at(dt_n_low, i) = DataVector(n, 0.);
  }
  const double root_half = 1. / std::sqrt(2.);
  std::array<DataVector, 4> k_low{};
  std::array<DataVector, 4> l_low{};
  std::array<DataVector, 4> dt_k_low{};
  std::array<DataVector, 4> dt_l_low{};
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(k_low, a) = root_half * (gsl::at(n_low, a) - gsl::at(s_low, a));
    gsl::at(l_low, a) = root_half * (gsl::at(n_low, a) + gsl::at(s_low, a));
    gsl::at(dt_k_low, a) =
        root_half * (gsl::at(dt_n_low, a) - gsl::at(dt_s_low, a));
    gsl::at(dt_l_low, a) =
        root_half * (gsl::at(dt_n_low, a) + gsl::at(dt_s_low, a));
  }
  for (size_t a = 0; a < 4; ++a) {
    gsl::at(out.dt_k_up, a) = DataVector(n, 0.);
    gsl::at(out.dt_l_up, a) = DataVector(n, 0.);
    for (size_t b = 0; b < 4; ++b) {
      gsl::at(out.dt_k_up, a) += dt_inv.get(a, b) * gsl::at(k_low, b) +
                                 inv.get(a, b) * gsl::at(dt_k_low, b);
      gsl::at(out.dt_l_up, a) += dt_inv.get(a, b) * gsl::at(l_low, b) +
                                 inv.get(a, b) * gsl::at(dt_l_low, b);
    }
  }
  for (size_t c = 0; c < 4; ++c) {
    for (size_t b = 0; b < 4; ++b) {
      DataVector& entry = gsl::at(gsl::at(out.dt_proj_ud, c), b);
      entry = DataVector(n, 0.);
      for (size_t a = 0; a < 4; ++a) {
        const DataVector p_ab = metric.get(a, b) +
                                gsl::at(k_low, a) * gsl::at(l_low, b) +
                                gsl::at(l_low, a) * gsl::at(k_low, b);
        const DataVector dt_p_ab = dt_metric.get(a, b) +
                                   gsl::at(dt_k_low, a) * gsl::at(l_low, b) +
                                   gsl::at(k_low, a) * gsl::at(dt_l_low, b) +
                                   gsl::at(dt_l_low, a) * gsl::at(k_low, b) +
                                   gsl::at(l_low, a) * gsl::at(dt_k_low, b);
        entry += dt_inv.get(c, a) * p_ab + inv.get(c, a) * dt_p_ab;
      }
    }
  }
  return out;
}

// d/dt of the gauge components {A, C, V} of a tensor u: the dt_u part
// contracted with the static frame plus the frame-motion terms.
GaugeComponents dt_gauge_components(const tnsr::aa<DataVector, 3>& u,
                                    const tnsr::aa<DataVector, 3>& dt_u,
                                    const SphereFrame& fr,
                                    const SphereFrameDot& frd) {
  GaugeComponents out = gauge_components(dt_u, fr);
  const size_t n = fr.n_points;
  DataVector u_dtl_l(n, 0.);
  DataVector u_dtk_l(n, 0.);
  DataVector u_k_dtl(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      u_dtl_l += u.get(a, b) * gsl::at(frd.dt_l_up, a) * gsl::at(fr.l_up, b);
      u_dtk_l += u.get(a, b) * gsl::at(frd.dt_k_up, a) * gsl::at(fr.l_up, b);
      u_k_dtl += u.get(a, b) * gsl::at(fr.k_up, a) * gsl::at(frd.dt_l_up, b);
    }
  }
  out.a += 2. * u_dtl_l;
  out.c += u_dtk_l + u_k_dtl;
  for (size_t b = 0; b < 4; ++b) {
    for (size_t c = 0; c < 4; ++c) {
      for (size_t d = 0; d < 4; ++d) {
        gsl::at(out.v, b) -= gsl::at(gsl::at(frd.dt_proj_ud, c), b) *
                                 u.get(c, d) * gsl::at(fr.l_up, d) +
                             gsl::at(gsl::at(fr.proj_ud, c), b) * u.get(c, d) *
                                 gsl::at(frd.dt_l_up, d);
      }
    }
  }
  return out;
}

std::vector<double> characteristic_gauge_time_derivative_modes_with_frame(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi, const Scalar<DataVector>& dt_gamma2,
    const SphereFrame& frame, const SphereFrameDot& frame_dot,
    const ylm::Spherepack& ylm_transform, const ModeSet& modes,
    const double normal_sign) {
  const size_t n_points = get<0, 0>(spacetime_metric).size();
  const tnsr::aa<DataVector, 3> u =
      u_minus_of(spacetime_metric, pi, phi, frame, normal_sign);
  tnsr::aa<DataVector, 3> dt_u(n_points);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_u.get(a, b) = dt_pi.get(a, b) -
                       frame.gamma2 * dt_spacetime_metric.get(a, b) -
                       get(dt_gamma2) * spacetime_metric.get(a, b);
      for (size_t k = 0; k < 3; ++k) {
        dt_u.get(a, b) +=
            normal_sign *
            (gsl::at(frame.normal_up, k) * dt_phi.get(k, a, b) +
             gsl::at(frame_dot.dt_normal_up, k) * phi.get(k, a, b));
      }
    }
  }
  return gauge_modes(dt_gauge_components(u, dt_u, frame, frame_dot),
                     ylm_transform, modes);
}

// sandwich helper: (g X g)_ab for an upper-index combination X
tnsr::aa<DataVector, 3> covariant_sandwich(const tnsr::aa<DataVector, 3>& g,
                                           const tnsr::AA<DataVector, 3>& x) {
  const size_t n = get<0, 0>(g).size();
  tnsr::aa<DataVector, 3> out(n, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          out.get(a, b) += g.get(a, c) * x.get(c, d) * g.get(d, b);
        }
      }
    }
  }
  return out;
}
}  // namespace

namespace detail {
std::vector<double> characteristic_gauge_modes(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const std::array<double, 3>& center, const ylm::Spherepack& ylm_transform,
    const size_t fit_l_max, const double normal_sign) {
  ASSERT(fit_l_max <= ylm_transform.l_max(),
         "fit_l_max " << fit_l_max << " exceeds the grid l_max "
                      << ylm_transform.l_max());
  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, center);
  return gauge_modes(
      gauge_components(
          u_minus_of(spacetime_metric, pi, phi, frame, normal_sign), frame),
      ylm_transform, kept_modes(ylm_transform, fit_l_max));
}

std::vector<double> characteristic_gauge_time_derivative_modes(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& dt_gamma2,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const std::array<double, 3>& center, const ylm::Spherepack& ylm_transform,
    const size_t fit_l_max, const double normal_sign) {
  ASSERT(fit_l_max <= ylm_transform.l_max(),
         "fit_l_max " << fit_l_max << " exceeds the grid l_max "
                      << ylm_transform.l_max());
  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, center);
  const SphereFrameDot frame_dot = build_frame_dot(
      spacetime_metric, dt_spacetime_metric, frame, inertial_coords, center);
  return characteristic_gauge_time_derivative_modes_with_frame(
      spacetime_metric, pi, phi, dt_spacetime_metric, dt_pi, dt_phi, dt_gamma2,
      frame, frame_dot, ylm_transform, kept_modes(ylm_transform, fit_l_max),
      normal_sign);
}
}  // namespace detail

RateFitResult fit_map_parameter_accelerations_uplus(
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& dt_spacetime_metric,
    const tnsr::aa<DataVector, 3>& dt_pi,
    const tnsr::iaa<DataVector, 3>& dt_phi, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& dt_gamma2,
    const std::array<double, num_map_parameters>& p_state,
    const std::array<double, num_map_parameters>& pdot_state,
    const tnsr::I<DataVector, 3>& inertial_coords,
    const ylm::Spherepack& ylm_transform, const MatcherConfig& config,
    const bool evaluate_held_out_uminus) {
  ASSERT(config.fit_l_max <= ylm_transform.l_max(),
         "FitLMax " << config.fit_l_max << " exceeds the grid l_max "
                    << ylm_transform.l_max());
  const size_t n_points = get<0, 0>(spacetime_metric).size();
  const ModeSet modes = kept_modes(ylm_transform, config.fit_l_max);
  const size_t n_modes = modes.indices.size();

  const SphereFrame frame =
      build_frame(spacetime_metric, gamma2, inertial_coords, config.center);
  const SphereFrameDot frame_dot =
      build_frame_dot(spacetime_metric, dt_spacetime_metric, frame,
                      inertial_coords, config.center);

  // ---- data side: the complete fixed-sphere derivative of u^+ and of its
  // numerical null-frame projection ----
  const std::vector<double> data_modes =
      characteristic_gauge_time_derivative_modes_with_frame(
          spacetime_metric, pi, phi, dt_spacetime_metric, dt_pi, dt_phi,
          dt_gamma2, frame, frame_dot, ylm_transform, modes, -1.);
  const std::vector<double> data_minus_modes =
      evaluate_held_out_uminus
          ? characteristic_gauge_time_derivative_modes_with_frame(
                spacetime_metric, pi, phi, dt_spacetime_metric, dt_pi, dt_phi,
                dt_gamma2, frame, frame_dot, ylm_transform, modes, 1.)
          : std::vector<double>{};

  // ---- model side at the current (p, pdot), data frame throughout ----
  std::array<DataVector, 3> y{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(y, i) = inertial_coords.get(i) - gsl::at(config.center, i);
  }
  tnsr::AA<DataVector, 3> big_g(n_points);
  gh::Solutions::affine_map_model::inverse_metric_combination(
      make_not_null(&big_g), y, config.mass, 1., p_state);
  const tnsr::aa<DataVector, 3> g_m = determinant_and_inverse(big_g).second;
  tnsr::AA<DataVector, 3> s_resp(n_points);
  gh::Solutions::affine_map_model::inverse_metric_combination(
      make_not_null(&s_resp), y, config.mass, 0., pdot_state);
  const tnsr::aa<DataVector, 3> dtg_m_neg = covariant_sandwich(g_m, s_resp);
  tnsr::aa<DataVector, 3> dtg_m(n_points);
  for (size_t st = 0; st < dtg_m.size(); ++st) {
    dtg_m[st] = -dtg_m_neg[st];
  }

  // Phi_m and dt Phi_m analytically
  tnsr::iAA<DataVector, 3> dk_big_g(n_points);
  gh::Solutions::affine_map_model::
      spatial_derivative_of_inverse_metric_combination(
          make_not_null(&dk_big_g), y, config.mass, 1., p_state);
  tnsr::iAA<DataVector, 3> dk_s(n_points);
  gh::Solutions::affine_map_model::
      spatial_derivative_of_inverse_metric_combination(
          make_not_null(&dk_s), y, config.mass, 0., pdot_state);
  tnsr::iaa<DataVector, 3> phi_m(n_points, 0.);
  tnsr::iaa<DataVector, 3> dt_phi_m(n_points, 0.);
  for (size_t k = 0; k < 3; ++k) {
    // Phi_k = -(g dk_G g); build per k, completing Phi_k before its dt
    tnsr::AA<DataVector, 3> slice(n_points);
    tnsr::AA<DataVector, 3> slice_s(n_points);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        slice.get(a, b) = dk_big_g.get(k, a, b);
        slice_s.get(a, b) = dk_s.get(k, a, b);
      }
    }
    const auto phi_k = covariant_sandwich(g_m, slice);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        phi_m.get(k, a, b) = -phi_k.get(a, b);
      }
    }
    // dt Phi_k = d_k(dt g) = -(Phi_k S g + g dkS g + g S Phi_k)
    const auto g_dks_g = covariant_sandwich(g_m, slice_s);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        DataVector term(n_points, 0.);
        for (size_t c = 0; c < 4; ++c) {
          for (size_t d = 0; d < 4; ++d) {
            term += phi_m.get(k, a, c) * s_resp.get(c, d) * g_m.get(d, b) +
                    g_m.get(a, c) * s_resp.get(c, d) * phi_m.get(k, d, b);
          }
        }
        dt_phi_m.get(k, a, b) = -(term + g_dks_g.get(a, b));
      }
    }
  }

  // model lapse/shift and their rates (dt G^-1_m = S exactly)
  const DataVector lapse_sq_m = -1. / get<0, 0>(big_g);
  const DataVector lapse_m = sqrt(lapse_sq_m);
  const DataVector dt_lapse_m = 0.5 * lapse_m * lapse_sq_m * get<0, 0>(s_resp);
  std::array<DataVector, 3> shift_m{};
  std::array<DataVector, 3> dt_shift_m{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(shift_m, i) = lapse_sq_m * big_g.get(0, i + 1);
    gsl::at(dt_shift_m, i) = 2. * lapse_m * dt_lapse_m * big_g.get(0, i + 1) +
                             lapse_sq_m * s_resp.get(0, i + 1);
  }

  // Pi_m and the pddot-independent part of dt Pi_m (the pddot part is the
  // column response below); d2t g at pddot = 0 is the Hessian term only
  tnsr::aa<DataVector, 3> pi_m(n_points);
  tnsr::aa<DataVector, 3> dt_pi_m0(n_points);
  const auto gsg = covariant_sandwich(g_m, s_resp);  // = -dt g_m
  tnsr::aa<DataVector, 3> d2tg_m0(n_points, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          d2tg_m0.get(a, b) +=
              2. * gsg.get(a, c) * s_resp.get(c, d) * g_m.get(d, b);
        }
      }
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      DataVector beta_phi(n_points, 0.);
      DataVector dt_beta_phi(n_points, 0.);
      for (size_t k = 0; k < 3; ++k) {
        beta_phi += gsl::at(shift_m, k) * phi_m.get(k, a, b);
        dt_beta_phi += gsl::at(dt_shift_m, k) * phi_m.get(k, a, b) +
                       gsl::at(shift_m, k) * dt_phi_m.get(k, a, b);
      }
      pi_m.get(a, b) = (beta_phi - dtg_m.get(a, b)) / lapse_m;
      dt_pi_m0.get(a, b) = (dt_beta_phi - d2tg_m0.get(a, b)) / lapse_m -
                           pi_m.get(a, b) * dt_lapse_m / lapse_m;
    }
  }

  if (config.centre_advection) {
    // The formulas above differentiate the affine response at a fixed model
    // center.  When qdot^i is nonzero the model also translates through the
    // fixed numerical sphere.  Evaluate that pddot-free transport derivative
    // by a symmetric time shift of the same production model used by the value
    // fit.  This includes, consistently, -qdot^k Phi_k in dt g and its spatial
    // derivative/Hessian contributions in dt Phi and dt Pi without introducing
    // a separate finite-difference approximation to the numerical data.
    constexpr double model_time_step = 1.e-3;
    std::array<double, num_map_parameters> p_plus = p_state;
    std::array<double, num_map_parameters> p_minus = p_state;
    std::array<double, 3> center_plus = config.center;
    std::array<double, 3> center_minus = config.center;
    for (size_t a = 0; a < num_map_parameters; ++a) {
      gsl::at(p_plus, a) += model_time_step * gsl::at(pdot_state, a);
      gsl::at(p_minus, a) -= model_time_step * gsl::at(pdot_state, a);
    }
    for (size_t i = 0; i < 3; ++i) {
      const double velocity = gsl::at(p_state, 4 + i);
      const double acceleration = gsl::at(pdot_state, 4 + i);
      gsl::at(center_plus, i) += model_time_step * velocity +
                                 0.5 * square(model_time_step) * acceleration;
      gsl::at(center_minus, i) += -model_time_step * velocity +
                                  0.5 * square(model_time_step) * acceleration;
    }

    tnsr::aa<DataVector, 3> model_g{};
    tnsr::aa<DataVector, 3> model_pi{};
    tnsr::iaa<DataVector, 3> model_phi{};
    gh::Solutions::affine_map_model::evolved_variables(
        make_not_null(&model_g), make_not_null(&model_pi),
        make_not_null(&model_phi), inertial_coords, config.mass, config.center,
        p_state, pdot_state, true);
    pi_m = std::move(model_pi);
    phi_m = std::move(model_phi);

    tnsr::aa<DataVector, 3> g_plus{};
    tnsr::aa<DataVector, 3> pi_plus{};
    tnsr::iaa<DataVector, 3> phi_plus{};
    tnsr::aa<DataVector, 3> g_minus{};
    tnsr::aa<DataVector, 3> pi_minus{};
    tnsr::iaa<DataVector, 3> phi_minus{};
    gh::Solutions::affine_map_model::evolved_variables(
        make_not_null(&g_plus), make_not_null(&pi_plus),
        make_not_null(&phi_plus), inertial_coords, config.mass, center_plus,
        p_plus, pdot_state, true);
    gh::Solutions::affine_map_model::evolved_variables(
        make_not_null(&g_minus), make_not_null(&pi_minus),
        make_not_null(&phi_minus), inertial_coords, config.mass, center_minus,
        p_minus, pdot_state, true);
    for (size_t st = 0; st < dtg_m.size(); ++st) {
      dtg_m[st] = (g_plus[st] - g_minus[st]) / (2. * model_time_step);
      dt_pi_m0[st] = (pi_plus[st] - pi_minus[st]) / (2. * model_time_step);
    }
    for (size_t st = 0; st < dt_phi_m.size(); ++st) {
      dt_phi_m[st] = (phi_plus[st] - phi_minus[st]) / (2. * model_time_step);
    }
  }

  const auto assemble_u_and_dt = [&](const tnsr::aa<DataVector, 3>& g_t,
                                     const tnsr::aa<DataVector, 3>& pi_t,
                                     const tnsr::iaa<DataVector, 3>& phi_t,
                                     const tnsr::aa<DataVector, 3>& dtg_t,
                                     const tnsr::aa<DataVector, 3>& dtpi_t,
                                     const tnsr::iaa<DataVector, 3>& dtphi_t,
                                     const double normal_sign) {
    const tnsr::aa<DataVector, 3> u =
        u_minus_of(g_t, pi_t, phi_t, frame, normal_sign);
    tnsr::aa<DataVector, 3> dt_u(n_points);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        dt_u.get(a, b) = dtpi_t.get(a, b) - frame.gamma2 * dtg_t.get(a, b) -
                         get(dt_gamma2) * g_t.get(a, b);
        for (size_t k = 0; k < 3; ++k) {
          dt_u.get(a, b) +=
              normal_sign * gsl::at(frame.normal_up, k) * dtphi_t.get(k, a, b) +
              normal_sign * gsl::at(frame_dot.dt_normal_up, k) *
                  phi_t.get(k, a, b);
        }
      }
    }
    return gauge_modes(dt_gauge_components(u, dt_u, frame, frame_dot),
                       ylm_transform, modes);
  };
  const std::vector<double> base_modes =
      assemble_u_and_dt(g_m, pi_m, phi_m, dtg_m, dt_pi_m0, dt_phi_m, -1.);
  const std::vector<double> base_minus_modes =
      evaluate_held_out_uminus
          ? assemble_u_and_dt(g_m, pi_m, phi_m, dtg_m, dt_pi_m0, dt_phi_m, 1.)
          : std::vector<double>{};

  // target for the linear solve: data minus the pddot-free model part
  std::vector<double> target(data_modes.size());
  for (size_t i = 0; i < target.size(); ++i) {
    target[i] = data_modes[i] - base_modes[i];
  }
  std::vector<double> minus_target = target;
  if (evaluate_held_out_uminus) {
    for (size_t i = 0; i < minus_target.size(); ++i) {
      minus_target[i] = data_minus_modes[i] - base_minus_modes[i];
    }
  }

  // ---- columns: pddot enters only through dt Pi_m -> -d2tg/alpha, so the
  // column tensor is (g_m R(E_A) g_m)/alpha_m contracted like dt_u (pure
  // dt-part, no frame-motion terms) ----
  constexpr size_t n_free = 9;
  std::vector<std::vector<double>> columns(n_free);
  for (size_t a = 0; a < n_free; ++a) {
    std::array<double, num_map_parameters> dir{};
    if (a == 0) {
      dir[0] = 1.;
      for (size_t i = 0; i < 3; ++i) {
        gsl::at(dir, 4 + i) = gsl::at(config.center_velocity, i);
      }
    } else if (a < 4) {
      gsl::at(dir, a) = 1.;
    } else {
      gsl::at(dir, 3 + a) = 1.;
      if (a == 4 or a == 7) {
        dir[12] = -1.;
      }
    }
    tnsr::AA<DataVector, 3> r_dir(n_points);
    gh::Solutions::affine_map_model::inverse_metric_combination(
        make_not_null(&r_dir), y, config.mass, 0., dir);
    auto col_tensor = covariant_sandwich(g_m, r_dir);
    for (size_t st = 0; st < col_tensor.size(); ++st) {
      col_tensor[st] /= lapse_m;
    }
    columns[a] =
        gauge_modes(gauge_components(col_tensor, frame), ylm_transform, modes);
  }

  // Block-weighted least squares. The Spherepack isotropy weights are already
  // inside the modal rows. Apply the fitted acceleration to D_T u^- as a
  // genuinely held-out test: those rows never enter the solve.
  const size_t n_res = target.size();
  const auto row_class = [&n_modes](const size_t row) -> size_t {
    const size_t field = row / n_modes;  // A, C, V_0..V_3
    return field == 0 ? 0 : (field == 1 ? 1 : 2);
  };
  std::vector<size_t> row_blocks(n_res);
  for (size_t i = 0; i < n_res; ++i) {
    row_blocks[i] = 5 * row_class(i) + modes.ells[i % n_modes];
  }
  const detail::WeightedModeFit fit = detail::fit_weighted_modes(
      target, minus_target, columns, row_blocks, config.uplus_block_weights);
  const std::vector<double>& x = fit.coefficients;

  RateFitResult result{};
  result.residual_initial = norm_of(target);
  result.residual_final = norm_of(fit.fitted_residual);
  if (evaluate_held_out_uminus) {
    result.minus_residual_initial = norm_of(minus_target);
    result.minus_residual_final = norm_of(fit.held_out_residual);
  }
  result.condition_number = fit.condition_number;

  // Relative closure and absolute RMS per {A, C, V} x l block. Keep both:
  // relative closure is ill-defined as a quality indicator when a quiet-hole
  // target is at roundoff, while absolute RMS alone hides scale differences.
  for (size_t c = 0; c < 3; ++c) {
    for (size_t l = 0; l <= 4; ++l) {
      double plus_residual_squared = 0.;
      double plus_target_squared = 0.;
      double minus_residual_squared = 0.;
      double minus_target_squared = 0.;
      size_t count = 0;
      for (size_t i = 0; i < n_res; ++i) {
        if (row_class(i) == c and modes.ells[i % n_modes] == l) {
          plus_residual_squared += square(fit.fitted_residual[i]);
          plus_target_squared += square(target[i]);
          if (evaluate_held_out_uminus) {
            minus_residual_squared += square(fit.held_out_residual[i]);
            minus_target_squared += square(minus_target[i]);
          }
          ++count;
        }
      }
      const size_t block = 5 * c + l;
      gsl::at(result.block_closure, block) =
          plus_target_squared > 0.
              ? std::sqrt(plus_residual_squared / plus_target_squared)
              : 0.;
      gsl::at(result.block_target_rms, block) =
          count > 0 ? std::sqrt(plus_target_squared / count) : 0.;
      gsl::at(result.block_residual_rms, block) =
          count > 0 ? std::sqrt(plus_residual_squared / count) : 0.;
      gsl::at(result.block_minus_closure, block) =
          minus_target_squared > 0.
              ? std::sqrt(minus_residual_squared / minus_target_squared)
              : 0.;
      gsl::at(result.block_minus_target_rms, block) =
          count > 0 ? std::sqrt(minus_target_squared / count) : 0.;
      gsl::at(result.block_minus_residual_rms, block) =
          count > 0 ? std::sqrt(minus_residual_squared / count) : 0.;
    }
  }

  result.pdot[0] = x[0];
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(result.pdot, 1 + i) = x[1 + i];
    gsl::at(result.pdot, 4 + i) = x[0] * gsl::at(config.center_velocity, i);
  }
  for (size_t i = 0; i < 5; ++i) {
    gsl::at(result.pdot, 7 + i) = x[4 + i];
  }
  result.pdot[12] = -x[4] - x[7];
  for (const double acceleration : result.pdot) {
    result.parameter_derivative_norm += square(acceleration);
  }
  result.parameter_derivative_norm =
      std::sqrt(result.parameter_derivative_norm);
  return result;
}
}  // namespace gh::Worldtube
