// Distributed under the MIT License.
// See LICENSE.txt for details.

// Profiles Shape::coords_frame_velocity_jacobian at production sizes and
// decomposes the cost into its Spherepack building blocks, plus a
// cached-basis-matrix (GEMM) evaluation of the same interpolations to bound
// the payoff of caching per-element interpolation data.
//
// Usage: ShapeJacobianBenchmark <num_points> <l_max> <iters> <warmup>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Shape.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/Gsl.hpp"

extern "C" void CkRegisterMainModule(void) {}

namespace {
using clk = std::chrono::steady_clock;

template <typename F>
double time_min_us(F&& fn, const size_t reps, const size_t iters) {
  double best = 1.0e30;
  for (size_t r = 0; r < reps; ++r) {
    const auto t0 = clk::now();
    for (size_t k = 0; k < iters; ++k) {
      fn();
    }
    const double dt =
        std::chrono::duration<double, std::micro>(clk::now() - t0).count() /
        static_cast<double>(iters);
    best = std::min(best, dt);
  }
  return best;
}
}  // namespace

int main(const int argc, const char* const* const argv) {
  const size_t num_points =
      argc > 1 ? static_cast<size_t>(std::strtoull(argv[1], nullptr, 10))
               : 1000;
  const size_t l_max =
      argc > 2 ? static_cast<size_t>(std::strtoull(argv[2], nullptr, 10)) : 8;
  const size_t iterations =
      argc > 3 ? static_cast<size_t>(std::strtoull(argv[3], nullptr, 10))
               : 10000;
  const size_t warmup =
      argc > 4 ? static_cast<size_t>(std::strtoull(argv[4], nullptr, 10)) : 200;

  // --- Shape map setup (BBH-like: sphere transition, small coefficients) ---
  const std::array<double, 3> center{0.0, 0.0, 0.0};
  const double r_min = 1.5;
  const double r_max = 10.0;
  domain::CoordinateMaps::TimeDependent::Shape shape_map{
      center, 1.0e-14,
      std::make_unique<domain::CoordinateMaps::ShapeMapTransitionFunctions::
                           SphereTransition>(r_min, r_max),
      "Shape", std::nullopt};

  const ylm::SpherepackIterator iter{l_max, l_max};
  DataVector coefs{iter.spherepack_array_size(), 0.0};
  // Deterministic smooth coefficients, decaying in l; keep the deformation
  // small so the map is well behaved over the annulus.
  {
    ylm::SpherepackIterator it{l_max, l_max};
    for (size_t l = 0; l <= l_max; ++l) {
      for (int m = -static_cast<int>(l); m <= static_cast<int>(l); ++m) {
        it.set(l, m);
        coefs[it()] =
            0.02 *
            std::cos(1.7 * static_cast<double>(l) +
                     0.9 * static_cast<double>(m)) /
            ((1.0 + static_cast<double>(l)) * (1.0 + static_cast<double>(l)));
      }
    }
  }
  DataVector dt_coefs{iter.spherepack_array_size(), 0.0};
  for (size_t i = 0; i < dt_coefs.size(); ++i) {
    dt_coefs[i] = 0.01 * coefs[i];
  }
  domain::CoordinateMaps::TimeDependent::Shape::FunctionsOfTimeMap
      functions_of_time{};
  functions_of_time["Shape"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0,
          std::array<DataVector, 3>{coefs, dt_coefs,
                                    DataVector{coefs.size(), 0.0}},
          100.0);
  const double time = 1.3;

  // --- source points spread over the annulus ---
  std::array<DataVector, 3> source_coords{
      DataVector{num_points}, DataVector{num_points}, DataVector{num_points}};
  for (size_t p = 0; p < num_points; ++p) {
    const double frac =
        static_cast<double>(p) / static_cast<double>(num_points);
    const double r = 2.0 + 6.0 * frac;
    const double theta =
        0.1 + 2.9 * std::fmod(1.618 * static_cast<double>(p), 1.0);
    const double phi = 6.28 * std::fmod(0.7548 * static_cast<double>(p), 1.0);
    source_coords[0][p] = r * std::sin(theta) * std::cos(phi);
    source_coords[1][p] = r * std::sin(theta) * std::sin(phi);
    source_coords[2][p] = r * std::cos(theta);
  }

  // --- 1. the full production call ---
  std::array<DataVector, 3> target_coords = source_coords;
  std::array<DataVector, 3> frame_vel{
      DataVector{num_points}, DataVector{num_points}, DataVector{num_points}};
  tnsr::Ij<DataVector, 3, Frame::NoFrame> jac{num_points};
  const auto full_call = [&]() {
    target_coords = source_coords;
    shape_map.coords_frame_velocity_jacobian(
        make_not_null(&target_coords), make_not_null(&frame_vel),
        make_not_null(&jac), time, functions_of_time);
  };
  for (size_t i = 0; i < warmup; ++i) {
    full_call();
  }
  const double t_full = time_min_us(full_call, 5, iterations);

  // --- 2. Spherepack building blocks at the same shapes ---
  // The jacobian path extends the gradient to l_max + 1, so measure both.
  for (const size_t l : {l_max, l_max + 1}) {
    const ylm::Spherepack ylm_srf{l, l};
    // target angles of the source points
    std::array<DataVector, 2> theta_phis{DataVector{num_points},
                                         DataVector{num_points}};
    for (size_t p = 0; p < num_points; ++p) {
      const double x = source_coords[0][p];
      const double y = source_coords[1][p];
      const double z = source_coords[2][p];
      const double r = std::sqrt(x * x + y * y + z * z);
      theta_phis[0][p] = std::acos(z / r);
      theta_phis[1][p] = std::atan2(y, x);
    }
    DataVector spec_coefs{ylm_srf.spectral_size(), 0.0};
    for (size_t i = 0; i < std::min(coefs.size(), spec_coefs.size()); ++i) {
      spec_coefs[i] = coefs[i];
    }
    DataVector phys{ylm_srf.physical_size(), 0.0};
    ylm_srf.spec_to_phys(make_not_null(phys.data()), spec_coefs.data());
    DataVector result{num_points};

    const auto info = ylm_srf.set_up_interpolation_info(theta_phis);
    const double t_setup = time_min_us(
        [&]() { (void)ylm_srf.set_up_interpolation_info(theta_phis); }, 5,
        iterations / 4 + 1);
    const double t_eval = time_min_us(
        [&]() {
          ylm_srf.interpolate_from_coefs(make_not_null(&result), spec_coefs,
                                         info);
        },
        5, iterations);
    const double t_interp = time_min_us(
        [&]() {
          ylm_srf.interpolate(make_not_null(&result), phys.data(), info);
        },
        5, iterations);
    const double t_p2s = time_min_us(
        [&]() {
          ylm_srf.phys_to_spec(make_not_null(spec_coefs.data()), phys.data());
        },
        5, iterations);
    const double t_grad =
        time_min_us([&]() { (void)ylm_srf.gradient_from_coefs(spec_coefs); }, 5,
                    iterations / 4 + 1);

    std::cout << "l_max " << l << ":\n"
              << "  set_up_interpolation_info: " << t_setup << " us\n"
              << "  interpolate_from_coefs   : " << t_eval << " us\n"
              << "  interpolate (p2s + eval) : " << t_interp << " us\n"
              << "  phys_to_spec             : " << t_p2s << " us\n"
              << "  gradient_from_coefs      : " << t_grad << " us\n";

    // --- 3. cached-basis-matrix evaluation of 5 fields (lever-2 bound) ---
    // Build B column-by-column with unit coefficient vectors (setup cost is
    // irrelevant: in the cached scheme it happens once per element), then
    // time the steady-state: one dgemm evaluating 5 fields at once.
    if (l == l_max) {
      const size_t k_coefs = ylm_srf.spectral_size();
      std::vector<double> basis(num_points * k_coefs);
      DataVector unit{k_coefs, 0.0};
      DataVector column{num_points};
      for (size_t k = 0; k < k_coefs; ++k) {
        unit = 0.0;
        unit[k] = 1.0;
        ylm_srf.interpolate_from_coefs(make_not_null(&column), unit, info);
        for (size_t p = 0; p < num_points; ++p) {
          basis[p + num_points * k] = column[p];
        }
      }
      std::vector<double> field_coefs(k_coefs * 5);
      for (size_t f = 0; f < 5; ++f) {
        for (size_t k = 0; k < k_coefs; ++k) {
          field_coefs[k + k_coefs * f] =
              spec_coefs[k] * (1.0 + 0.1 * static_cast<double>(f));
        }
      }
      std::vector<double> out(num_points * 5);
      const double t_gemm = time_min_us(
          [&]() {
            dgemm_<true>('N', 'N', num_points, 5, k_coefs, 1.0, basis.data(),
                         num_points, field_coefs.data(), k_coefs, 0.0,
                         out.data(), num_points);
          },
          5, iterations);
      // correctness of the B-based evaluation vs Clenshaw
      ylm_srf.interpolate_from_coefs(make_not_null(&result), spec_coefs, info);
      double max_diff = 0.0;
      for (size_t p = 0; p < num_points; ++p) {
        max_diff = std::max(max_diff, std::abs(out[p] - result[p]));
      }
      std::cout << "  basis matrix B           : " << num_points << " x "
                << k_coefs << " (" << num_points * k_coefs * 8 / 1024
                << " KB)\n"
                << "  gemm eval of 5 fields    : " << t_gemm << " us"
                << "  (max |B*c - clenshaw| = " << max_diff << ")\n";
    }
  }

  std::cout << "num_points: " << num_points << "\n"
            << "l_max: " << l_max << "\n"
            << "full coords_frame_velocity_jacobian: " << t_full << " us\n"
            << "checksum: " << jac.get(0, 0)[0] + jac.get(1, 2)[num_points / 2]
            << "\n";
  return 0;
}
