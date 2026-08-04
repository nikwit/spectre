// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "PointwiseFunctions/SpecialRelativity/LorentzBoostMatrix.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace {
constexpr double mass = 1.3;
const std::array<double, 3> centre{{0.2, -0.1, 0.4}};
const std::array<double, 13> zero_p{};

// A handful of points comfortably outside the horizon (harmonic r = M).
tnsr::I<DataVector, 3> sample_points() {
  const DataVector x{2.5, -3.1, 0.7, 4.2, -1.9};
  const DataVector y{-1.4, 2.2, 3.6, -0.8, 2.7};
  const DataVector z{0.9, 1.1, -2.8, 3.3, -3.9};
  tnsr::I<DataVector, 3> out(x.size());
  get<0>(out) = x + centre[0];
  get<1>(out) = y + centre[1];
  get<2>(out) = z + centre[2];
  return out;
}

// The rest-frame image of the lab points, computed independently of the
// production code: xbar^abar = Lambda(-v)^abar_b (t, x - centre)^b.
tnsr::I<DataVector, 3> rest_frame_points(const tnsr::I<DataVector, 3>& x,
                                         const double time,
                                         const std::array<double, 3>& v) {
  tnsr::I<double, 3, Frame::NoFrame> minus_v{};
  for (size_t i = 0; i < 3; ++i) {
    minus_v.get(i) = -gsl::at(v, i);
  }
  const auto matrix = sr::lorentz_boost_matrix(minus_v);
  tnsr::I<DataVector, 3> out(get<0>(x).size());
  for (size_t i = 0; i < 3; ++i) {
    out.get(i) = matrix.get(i + 1, 0) * time;
    for (size_t j = 0; j < 3; ++j) {
      out.get(i) += matrix.get(i + 1, j + 1) * (x.get(j) - gsl::at(centre, j));
    }
  }
  return out;
}
}  // namespace

// The headline invariant: boost the solution, then undo the boost by hand,
// and the static solution must come back. The inverse transformation here is
// written from lorentz_boost_matrix directly rather than by calling the
// production routine with -v, so a sign or index error in the boost cannot
// cancel itself.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.InverseBoost",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const size_t n_points = get<0>(x).size();
  for (const double time : {0., 1.7, -3.4}) {
    for (const std::array<double, 3>& velocity :
         {std::array<double, 3>{{0.3, 0., 0.}},
          std::array<double, 3>{{0., -0.25, 0.}},
          std::array<double, 3>{{0.1, 0.2, -0.15}}}) {
      tnsr::aa<DataVector, 3> boosted_metric{};
      tnsr::aa<DataVector, 3> boosted_pi{};
      tnsr::iaa<DataVector, 3> boosted_phi{};
      gh::Solutions::affine_map_model::boosted_evolved_variables(
          make_not_null(&boosted_metric), make_not_null(&boosted_pi),
          make_not_null(&boosted_phi), x, time, mass, centre, zero_p, zero_p,
          velocity);

      // Undo the boost: with N^a_abar = dx^a / dxbar^abar = Lambda(+v),
      // gbar_abar bbar = N^a_abar N^b_bbar g_ab.
      tnsr::I<double, 3, Frame::NoFrame> plus_v{};
      for (size_t i = 0; i < 3; ++i) {
        plus_v.get(i) = gsl::at(velocity, i);
      }
      const auto inverse_matrix = sr::lorentz_boost_matrix(plus_v);
      tnsr::aa<DataVector, 3> unboosted_metric(n_points, 0.);
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = a; b < 4; ++b) {
          for (size_t c = 0; c < 4; ++c) {
            for (size_t d = 0; d < 4; ++d) {
              unboosted_metric.get(a, b) += inverse_matrix.get(c, a) *
                                            inverse_matrix.get(d, b) *
                                            boosted_metric.get(c, d);
            }
          }
        }
      }

      // ... and compare with the static solution at the rest-frame points.
      const auto rest_points = rest_frame_points(x, time, velocity);
      tnsr::aa<DataVector, 3> expected_metric{};
      tnsr::aa<DataVector, 3> expected_pi{};
      tnsr::iaa<DataVector, 3> expected_phi{};
      gh::Solutions::affine_map_model::evolved_variables(
          make_not_null(&expected_metric), make_not_null(&expected_pi),
          make_not_null(&expected_phi), rest_points, mass, {{0., 0., 0.}},
          zero_p, zero_p);

      Approx custom = Approx::custom().epsilon(1.e-12).scale(1.);
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = a; b < 4; ++b) {
          CHECK_ITERABLE_CUSTOM_APPROX(unboosted_metric.get(a, b),
                                       expected_metric.get(a, b), custom);
        }
      }
      // The boost must actually make the solution time dependent. Note
      // that Pi itself does not vanish in the static frame: this harmonic
      // form has a nonzero shift (G^{0i} = 4 M^2 n^i / r^2), and
      // Pi = (beta^k Phi_k - dt g)/alpha. So compare dt g, recovered from
      // each triple, which is zero in the rest frame and not in the lab.
      const auto dt_metric_of = [](const tnsr::aa<DataVector, 3>& g,
                                   const tnsr::aa<DataVector, 3>& local_pi,
                                   const tnsr::iaa<DataVector, 3>& local_phi) {
        const auto inv = determinant_and_inverse(g).second;
        const DataVector lapse = 1. / sqrt(-get<0, 0>(inv));
        DataVector out = -lapse * get<0, 0>(local_pi);
        for (size_t k = 0; k < 3; ++k) {
          out += -inv.get(0, k + 1) / get<0, 0>(inv) * local_phi.get(k, 0, 0);
        }
        return out;
      };
      CHECK(max(abs(dt_metric_of(expected_metric, expected_pi, expected_phi))) <
            1.e-13);
      CHECK(max(abs(dt_metric_of(boosted_metric, boosted_pi, boosted_phi))) >
            1.e-4);
    }
  }
}

// Zero velocity must reduce to the unboosted routine exactly, so that
// turning the option on for a hole at rest cannot perturb a run.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.ZeroVelocity",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const std::array<double, 13> p{{1.e-3, 2.e-3, -1.e-3, 5.e-4, 0., 0., 0.,
                                  1.e-3, -2.e-4, 3.e-4, 5.e-4, 1.e-4, -1.e-3}};
  const std::array<double, 13> pdot{
      {1.e-5, 0., 2.e-5, 0., 0., 0., 0., 3.e-5, 0., 0., -1.e-5, 0., 2.e-5}};
  tnsr::aa<DataVector, 3> metric_a{};
  tnsr::aa<DataVector, 3> pi_a{};
  tnsr::iaa<DataVector, 3> phi_a{};
  gh::Solutions::affine_map_model::boosted_evolved_variables(
      make_not_null(&metric_a), make_not_null(&pi_a), make_not_null(&phi_a), x,
      2.5, mass, centre, p, pdot, {{0., 0., 0.}});
  tnsr::aa<DataVector, 3> metric_b{};
  tnsr::aa<DataVector, 3> pi_b{};
  tnsr::iaa<DataVector, 3> phi_b{};
  gh::Solutions::affine_map_model::evolved_variables(
      make_not_null(&metric_b), make_not_null(&pi_b), make_not_null(&phi_b), x,
      mass, centre, p, pdot);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      CHECK_ITERABLE_APPROX(metric_a.get(a, b), metric_b.get(a, b));
      CHECK_ITERABLE_APPROX(pi_a.get(a, b), pi_b.get(a, b));
      for (size_t k = 0; k < 3; ++k) {
        CHECK_ITERABLE_APPROX(phi_a.get(k, a, b), phi_b.get(k, a, b));
      }
    }
  }
}

// The Dhesi slow-time first-order model must be exactly affine in every
// retained coefficient after epsilon is set to one.  Any even-in-p remainder
// here would be an uncontrolled O(epsilon^2) resummation.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.StrictFirstOrderLinearity",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const std::array<double, 13> p{{2.e-3, -1.e-3, 3.e-3, -2.e-3, 5.e-4, -7.e-4,
                                  9.e-4, 4.e-4, -3.e-4, 2.e-4, -6.e-4, 8.e-4,
                                  2.e-4}};
  std::array<double, 13> minus_p{};
  for (size_t a = 0; a < p.size(); ++a) {
    gsl::at(minus_p, a) = -gsl::at(p, a);
  }
  const auto variables = [&x](const std::array<double, 13>& parameters) {
    tnsr::aa<DataVector, 3> metric{};
    tnsr::aa<DataVector, 3> local_pi{};
    tnsr::iaa<DataVector, 3> local_phi{};
    gh::Solutions::affine_map_model::first_order_evolved_variables(
        make_not_null(&metric), make_not_null(&local_pi),
        make_not_null(&local_phi), x, mass, centre, parameters);
    return std::make_tuple(metric, local_pi, local_phi);
  };
  const auto [background_metric, background_pi, background_phi] =
      variables(zero_p);
  const auto [plus_metric, plus_pi, plus_phi] = variables(p);
  const auto [minus_metric, minus_pi, minus_phi] = variables(minus_p);
  Approx custom = Approx::custom().epsilon(2.e-12).scale(1.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      const DataVector metric_sum =
          plus_metric.get(a, b) + minus_metric.get(a, b);
      const DataVector twice_background_metric =
          2. * background_metric.get(a, b);
      CHECK_ITERABLE_CUSTOM_APPROX(metric_sum, twice_background_metric, custom);
      const DataVector pi_sum = plus_pi.get(a, b) + minus_pi.get(a, b);
      const DataVector twice_background_pi = 2. * background_pi.get(a, b);
      CHECK_ITERABLE_CUSTOM_APPROX(pi_sum, twice_background_pi, custom);
      for (size_t k = 0; k < 3; ++k) {
        const DataVector phi_sum =
            plus_phi.get(k, a, b) + minus_phi.get(k, a, b);
        const DataVector twice_background_phi =
            2. * background_phi.get(k, a, b);
        CHECK_ITERABLE_CUSTOM_APPROX(phi_sum, twice_background_phi, custom);
      }
    }
  }
}

// At t=0 an exact constant Lorentz boost has first-order coefficients
// beta_z=qdot_z=v.  The strict model must differ from it by O(v^2), including
// Pi and Phi, rather than by an O(v) derivative-ordering error.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.StrictFirstOrderBoostScaling",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const auto error = [&x](const double velocity) {
    std::array<double, 13> p{};
    p[3] = velocity;
    p[6] = velocity;
    tnsr::aa<DataVector, 3> first_order_metric{};
    tnsr::aa<DataVector, 3> first_order_pi{};
    tnsr::iaa<DataVector, 3> first_order_phi{};
    gh::Solutions::affine_map_model::first_order_evolved_variables(
        make_not_null(&first_order_metric), make_not_null(&first_order_pi),
        make_not_null(&first_order_phi), x, mass, centre, p);
    tnsr::aa<DataVector, 3> exact_metric{};
    tnsr::aa<DataVector, 3> exact_pi{};
    tnsr::iaa<DataVector, 3> exact_phi{};
    gh::Solutions::affine_map_model::boosted_evolved_variables(
        make_not_null(&exact_metric), make_not_null(&exact_pi),
        make_not_null(&exact_phi), x, 0., mass, centre, zero_p, zero_p,
        {{0., 0., velocity}});
    double squared_error = 0.;
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        squared_error +=
            sum(square(first_order_metric.get(a, b) - exact_metric.get(a, b)));
        squared_error +=
            sum(square(first_order_pi.get(a, b) - exact_pi.get(a, b)));
        for (size_t k = 0; k < 3; ++k) {
          squared_error += sum(
              square(first_order_phi.get(k, a, b) - exact_phi.get(k, a, b)));
        }
      }
    }
    return std::sqrt(squared_error);
  };
  const double error_at_v = error(1.e-3);
  const double error_at_half_v = error(5.e-4);
  CAPTURE(error_at_v, error_at_half_v, error_at_v / error_at_half_v);
  CHECK(error_at_v > 0.);
  CHECK(error_at_v / error_at_half_v == Approx::custom().epsilon(5.e-3)(4.));
}

// Pi and Phi are assembled from transformed derivatives, which is where a
// sign error would hide. Check them against finite differences of the
// boosted metric -- the one place a finite difference is the right tool,
// since it is the independent standard rather than the production path.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.BoostedDerivatives",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const size_t n_points = get<0>(x).size();
  const std::array<double, 3> velocity{{0.12, -0.07, 0.2}};
  const double time = 0.9;
  const double step = 1.e-6;

  const auto metric_at = [&x, &velocity](const double t, const size_t direction,
                                         const double shift) {
    auto shifted = x;
    if (direction < 3) {
      shifted.get(direction) += shift;
    }
    tnsr::aa<DataVector, 3> metric{};
    tnsr::aa<DataVector, 3> pi{};
    tnsr::iaa<DataVector, 3> phi{};
    gh::Solutions::affine_map_model::boosted_evolved_variables(
        make_not_null(&metric), make_not_null(&pi), make_not_null(&phi),
        shifted, t, mass, centre, zero_p, zero_p, velocity);
    return metric;
  };

  tnsr::aa<DataVector, 3> metric{};
  tnsr::aa<DataVector, 3> pi{};
  tnsr::iaa<DataVector, 3> phi{};
  gh::Solutions::affine_map_model::boosted_evolved_variables(
      make_not_null(&metric), make_not_null(&pi), make_not_null(&phi), x, time,
      mass, centre, zero_p, zero_p, velocity);
  const auto inverse_metric = determinant_and_inverse(metric).second;
  const DataVector lapse = 1. / sqrt(-get<0, 0>(inverse_metric));

  Approx custom = Approx::custom().epsilon(1.e-6).scale(1.e-3);
  // Phi_kab = d_k g_ab
  for (size_t k = 0; k < 3; ++k) {
    const auto plus = metric_at(time, k, step);
    const auto minus = metric_at(time, k, -step);
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        const DataVector fd = (plus.get(a, b) - minus.get(a, b)) / (2. * step);
        CHECK_ITERABLE_CUSTOM_APPROX(phi.get(k, a, b), fd, custom);
      }
    }
  }
  // Pi_ab = (beta^k Phi_kab - d_t g_ab) / alpha
  const auto later = metric_at(time + step, 3, 0.);
  const auto earlier = metric_at(time - step, 3, 0.);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      const DataVector dt_metric =
          (later.get(a, b) - earlier.get(a, b)) / (2. * step);
      DataVector expected(n_points, 0.);
      for (size_t k = 0; k < 3; ++k) {
        const DataVector shift_k =
            -inverse_metric.get(0, k + 1) / get<0, 0>(inverse_metric);
        expected += shift_k * phi.get(k, a, b);
      }
      expected = (expected - dt_metric) / lapse;
      CHECK_ITERABLE_CUSTOM_APPROX(pi.get(a, b), expected, custom);
    }
  }
}

// Is a boost inside the span of the thirteen response columns? Project the
// O(v) part of the exact boosted metric at t = 0 (where the secular
// displacement -v t is absent, so only the tensor-transformation content
// remains) onto the columns and report both the coefficients and what is
// left over. If the leftover is at round-off the model can represent a
// boost and any online failure is the fit's; if not, it cannot.
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Gr."
    "AffineMappedHarmonicSchwarzschild.BoostInSpan",
    "[PointwiseFunctions][Unit]") {
  const auto x = sample_points();
  const size_t n_points = get<0>(x).size();
  const double v = 1.e-4;
  const double eps = 1.e-6;

  const auto model = [&x](const std::array<double, 13>& p) {
    tnsr::aa<DataVector, 3> g{};
    tnsr::aa<DataVector, 3> local_pi{};
    tnsr::iaa<DataVector, 3> local_phi{};
    gh::Solutions::affine_map_model::first_order_evolved_variables(
        make_not_null(&g), make_not_null(&local_pi), make_not_null(&local_phi),
        x, mass, centre, p);
    return g;
  };

  const auto base = model(zero_p);
  tnsr::aa<DataVector, 3> boosted{};
  tnsr::aa<DataVector, 3> boosted_pi{};
  tnsr::iaa<DataVector, 3> boosted_phi{};
  gh::Solutions::affine_map_model::boosted_evolved_variables(
      make_not_null(&boosted), make_not_null(&boosted_pi),
      make_not_null(&boosted_phi), x, 0., mass, centre, zero_p, zero_p,
      {{0., 0., v}});

  // flatten (10 components x n_points)
  const size_t n_rows = 10 * n_points;
  std::vector<double> target(n_rows, 0.);
  std::vector<std::vector<double>> columns(13, std::vector<double>(n_rows, 0.));
  size_t row = 0;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      for (size_t k = 0; k < n_points; ++k) {
        target[row + k] = (boosted.get(a, b)[k] - base.get(a, b)[k]) / v;
      }
      row += n_points;
    }
  }
  for (size_t col = 0; col < 13; ++col) {
    std::array<double, 13> p{};
    gsl::at(p, col) = eps;
    const auto perturbed = model(p);
    row = 0;
    for (size_t a = 0; a < 4; ++a) {
      for (size_t b = a; b < 4; ++b) {
        for (size_t k = 0; k < n_points; ++k) {
          columns[col][row + k] =
              (perturbed.get(a, b)[k] - base.get(a, b)[k]) / eps;
        }
        row += n_points;
      }
    }
  }

  // normal equations, Gaussian elimination with partial pivoting
  std::vector<std::vector<double>> ata(13, std::vector<double>(14, 0.));
  for (size_t i = 0; i < 13; ++i) {
    for (size_t j = 0; j < 13; ++j) {
      double sum = 0.;
      for (size_t r = 0; r < n_rows; ++r) {
        sum += columns[i][r] * columns[j][r];
      }
      ata[i][j] = sum;
    }
    double rhs = 0.;
    for (size_t r = 0; r < n_rows; ++r) {
      rhs += columns[i][r] * target[r];
    }
    ata[i][13] = rhs;
  }
  for (size_t i = 0; i < 13; ++i) {
    size_t pivot = i;
    for (size_t r = i + 1; r < 13; ++r) {
      if (std::abs(ata[r][i]) > std::abs(ata[pivot][i])) {
        pivot = r;
      }
    }
    std::swap(ata[i], ata[pivot]);
    for (size_t r = 0; r < 13; ++r) {
      if (r == i or std::abs(ata[i][i]) < 1.e-300) {
        continue;
      }
      const double factor = ata[r][i] / ata[i][i];
      for (size_t col = i; col < 14; ++col) {
        ata[r][col] -= factor * ata[i][col];
      }
    }
  }
  std::array<double, 13> coefficients{};
  for (size_t i = 0; i < 13; ++i) {
    gsl::at(coefficients, i) =
        std::abs(ata[i][i]) < 1.e-300 ? 0. : ata[i][13] / ata[i][i];
  }

  double target_norm = 0.;
  double leftover = 0.;
  for (size_t r = 0; r < n_rows; ++r) {
    double fit = 0.;
    for (size_t col = 0; col < 13; ++col) {
      fit += gsl::at(coefficients, col) * columns[col][r];
    }
    leftover += square(target[r] - fit);
    target_norm += square(target[r]);
  }
  const double relative = std::sqrt(leftover / target_norm);

  const std::array<std::string, 13> names{
      {"qdot0", "b_x", "b_y", "b_z", "qdot_x", "qdot_y", "qdot_z", "s_xx",
       "s_xy", "s_xz", "s_yy", "s_yz", "s_zz"}};
  for (size_t i = 0; i < 13; ++i) {
    INFO("coefficient of " + gsl::at(names, i) +
         " per unit v: " + std::to_string(gsl::at(coefficients, i)));
    CHECK(true);
  }
  // A boost is exactly beta_i = qdot^i = +v_i and nothing else.
  // to O(v), the accuracy of the finite difference in v that built target
  Approx unit = Approx::custom().epsilon(2. * v).scale(1.);
  CHECK(unit(gsl::at(coefficients, 3)) == 1.);  // b_z
  CHECK(unit(gsl::at(coefficients, 6)) == 1.);  // qdot_z
  for (const size_t other :
       {0_st, 1_st, 2_st, 4_st, 5_st, 7_st, 8_st, 9_st, 10_st, 11_st, 12_st}) {
    INFO("spurious content in " + gsl::at(names, other));
    CHECK(std::abs(gsl::at(coefficients, other)) < 20. * v);
  }
  // What is left over is the O(v) remainder of the finite difference in v,
  // not a gap in the model space, so it must scale linearly with v. Assert
  // that rather than a fixed floor: the model DOES contain a boost.
  INFO("relative leftover after projecting the boost onto the columns: " +
       std::to_string(relative));
  CHECK(relative < 2. * v);
  CHECK(relative > 0.1 * v);
}
