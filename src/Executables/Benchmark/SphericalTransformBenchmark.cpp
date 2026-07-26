// Distributed under the MIT License.
// See LICENSE.txt for details.

// Benchmark harness for the SPHEREPACK-replacement study: scaling of the
// production spherical-harmonic operations over l_max = 8..100 and A/B
// comparisons against candidate re-implementations.
//
// Usage: SphericalTransformBenchmark <l_max> <n_r> <n_comp> <iters> <warmup>
//                                    <mode>
//   mode 0 = SPHEREPACK baseline: phys_to_spec / spec_to_phys / gradient
//            (each via the *_all_offsets production entry points, looped over
//            n_comp components with radial stride n_r), the filter round trip
//            (analysis + diagonal modal scaling + synthesis), construction
//            time of a ylm::Spherepack, and its constant workspace size.
//   mode 1 = libsharp on the identical grid and memory layout: map2alm,
//            alm2map, alm2map_deriv1 with ntrans = n_comp*n_r simultaneous
//            transforms. Validates against SPHEREPACK first: synthesis round
//            trip to roundoff, gradient to roundoff for fields band-limited
//            to l <= l_max-1, and prints the expected O(1e-2) top-degree
//            disagreement caused by SPHEREPACK's (l = l_max, m even) gradient
//            truncation.
//   mode 2 = fused per-m theta matrices applied with dgemm + SPHEREPACK's
//            hrfft phi FFTs (the high-l variant of the mode-7 prototype in
//            PartialDerivativesBenchmark): analysis/synthesis pair and the
//            gradient, both validated against SPHEREPACK. Reports the
//            per-l_max matrix table memory.
//   mode 3 = arbitrary-point evaluation (the Shape-map / horizon path):
//            Clenshaw interpolate_from_coefs for K=1000 points x 5 fields
//            (DataVector batch, the Shape jacobian pattern), per-point
//            double calls (the Shape::inverse pattern), and the equivalent
//            cached-dense-matrix dgemm.
//
// All timings are min-of-<iters> wall-clock per full operation (all
// components), printed as "  <op>: <us> us/callset  <us>/comp  <ns>/pt".

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Spherepack.hpp"

#include <sharp_almhelpers.h>
#include <sharp_geomhelpers.h>
#include <sharp_lowlevel.h>

extern "C" void CkRegisterMainModule(void) {}

// Batched real FFTs from the vendored SPHEREPACK (FFTPACK halfcomplex
// layout, batch index fastest). Not declared in Utilities/Spherepack.hpp.
extern "C" {
void hrffti_(const int&, double*);
void hrfftf_(const int&, const int&, double*, const int&, const double*,
             double*);
void hrfftb_(const int&, const int&, double*, const int&, const double*,
             double*);
}

namespace {

// --- timing -----------------------------------------------------------

template <typename F>
double time_min_us(const size_t iters, const size_t warmup, F&& f) {
  for (size_t i = 0; i < warmup; ++i) {
    f();
  }
  double best = std::numeric_limits<double>::max();
  for (size_t i = 0; i < iters; ++i) {
    const auto start = std::chrono::steady_clock::now();
    f();
    const std::chrono::duration<double, std::micro> dt =
        std::chrono::steady_clock::now() - start;
    best = std::min(best, dt.count());
  }
  return best;
}

struct Sizes {
  size_t l_max;
  size_t n_r;
  size_t n_comp;
  size_t nth;
  size_t nph;
  size_t npts;       // n_r * nth * nph, one component
  size_t spec_size;  // SPHEREPACK padded coefficient buffer, one offset
};

void print_op(const std::string& name, const double us_per_callset,
              const Sizes& s) {
  const double per_comp = us_per_callset / static_cast<double>(s.n_comp);
  const double ns_per_pt =
      1.0e3 * us_per_callset / static_cast<double>(s.n_comp * s.npts);
  std::cout << "  " << name << ": " << us_per_callset << " us/callset  "
            << per_comp << " us/comp  " << ns_per_pt << " ns/pt\n";
}

double max_abs(const double* const d, const size_t n) {
  double m = 0.0;
  for (size_t i = 0; i < n; ++i) {
    m = std::max(m, std::abs(d[i]));
  }
  return m;
}

double max_rel_diff(const double* const a, const double* const b,
                    const size_t n) {
  const double scale = std::max(max_abs(a, n), max_abs(b, n));
  if (scale == 0.0) {
    return 0.0;
  }
  double m = 0.0;
  for (size_t i = 0; i < n; ++i) {
    m = std::max(m, std::abs(a[i] - b[i]));
  }
  return m / scale;
}

// Band-limited random data: fill the valid SPHEREPACK coefficients with
// uniform random values (optionally zeroing l > l_cut) and synthesize.
void random_coefs(const gsl::not_null<std::vector<double>*> coefs,
                  const size_t l_max, const size_t stride, const size_t l_cut,
                  std::mt19937* const gen) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::fill(coefs->begin(), coefs->end(), 0.0);
  for (ylm::SpherepackIterator it(l_max, l_max); it; ++it) {
    if (it.l() > l_cut) {
      continue;
    }
    for (size_t s = 0; s < stride; ++s) {
      (*coefs)[it() * stride + s] = dist(*gen);
    }
  }
}

// --- mode 0: SPHEREPACK baseline ---------------------------------------

// Constant workspace of one ylm::Spherepack (formulas from
// SpherepackHelper.cpp).
size_t spherepack_const_storage_bytes(const size_t l_max, const size_t m_max) {
  const auto n_theta = static_cast<long>(l_max + 1);
  const auto n_phi = static_cast<long>(2 * m_max + 1);
  const auto l1 = static_cast<long>(m_max + 1);
  const long l2 = (n_theta + 1) / 2;
  const long scalar_work = n_phi + 15 + n_theta * (3 * (l1 + l2) - 2) +
                           (l1 - 1) * (l2 * (2 * n_theta - l1) - 3 * l1) / 2;
  const long vector_work =
      n_theta * l2 * (n_theta + 1) + n_phi + 15 + 2 * n_theta;
  const long quad_work = n_theta * n_phi;
  const long interp_work = n_theta * l1 - l1 * (l1 - 1) / 2;
  const long total = 2 * scalar_work + vector_work + quad_work +
                     4 * (static_cast<long>(l_max) + 1) + 2 * n_phi +
                     2 * interp_work + (l1 - 1 + 1) +
                     interp_work /* index array, size_t */;
  return static_cast<size_t>(total) * 8;
}

void run_spherepack_mode(const Sizes& s, const size_t iters,
                         const size_t warmup) {
  const ylm::Spherepack ylm(s.l_max, s.l_max);
  std::mt19937 gen(314159);

  std::vector<double> coefs(s.spec_size * s.n_r * s.n_comp);
  std::vector<double> u(s.npts * s.n_comp);
  std::vector<double> u2(s.npts * s.n_comp);
  std::vector<double> dth(s.npts * s.n_comp);
  std::vector<double> dph(s.npts * s.n_comp);
  for (size_t c = 0; c < s.n_comp; ++c) {
    std::vector<double> ctmp(s.spec_size * s.n_r);
    random_coefs(&ctmp, s.l_max, s.n_r, s.l_max, &gen);
    ylm.spec_to_phys_all_offsets(make_not_null(u.data() + c * s.npts),
                                 make_not_null(ctmp.data()), s.n_r);
  }

  // Diagonal modal filter factors (exponential in l), applied via the
  // iterator like the production tensor-Ylm filter's modal stage (the real
  // filter multiplies a sparse tensor-Ylm matrix; this is a lower bound on
  // that stage but the round trip is transform-dominated).
  std::vector<double> filter_factor(s.l_max + 1);
  for (size_t l = 0; l <= s.l_max; ++l) {
    filter_factor[l] = std::exp(
        -36.0 *
        std::pow(static_cast<double>(l) / static_cast<double>(s.l_max), 32));
  }

  const double t_analysis = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      ylm.phys_to_spec_all_offsets(
          make_not_null(coefs.data() + c * s.spec_size * s.n_r),
          make_not_null(u.data() + c * s.npts), s.n_r);
    }
  });
  print_op("phys_to_spec", t_analysis, s);

  const double t_synthesis = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      ylm.spec_to_phys_all_offsets(
          make_not_null(u2.data() + c * s.npts),
          make_not_null(coefs.data() + c * s.spec_size * s.n_r), s.n_r);
    }
  });
  print_op("spec_to_phys", t_synthesis, s);

  const double t_gradient = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      const std::array<double*, 2> df{dth.data() + c * s.npts,
                                      dph.data() + c * s.npts};
      ylm.gradient_all_offsets(df, make_not_null(u.data() + c * s.npts),
                               s.n_r);
    }
  });
  print_op("gradient", t_gradient, s);

  const double t_filter = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      double* const cc = coefs.data() + c * s.spec_size * s.n_r;
      ylm.phys_to_spec_all_offsets(make_not_null(cc),
                                   make_not_null(u.data() + c * s.npts),
                                   s.n_r);
      for (ylm::SpherepackIterator it(s.l_max, s.l_max); it; ++it) {
        const double f = filter_factor[it.l()];
        for (size_t r = 0; r < s.n_r; ++r) {
          cc[it() * s.n_r + r] *= f;
        }
      }
      ylm.spec_to_phys_all_offsets(make_not_null(u2.data() + c * s.npts),
                                   make_not_null(cc), s.n_r);
    }
  });
  print_op("filter_roundtrip", t_filter, s);

  const double t_construct = time_min_us(3, 1, [&]() {
    const ylm::Spherepack fresh(s.l_max, s.l_max);
    (void)fresh;
  });
  std::cout << "  construction: " << t_construct << " us\n";
  std::cout << "  const_workspace: "
            << spherepack_const_storage_bytes(s.l_max, s.l_max) << " bytes\n";
}

// --- mode 1: libsharp ---------------------------------------------------

struct SharpPlan {
  sharp_geom_info* geom = nullptr;
  sharp_alm_info* alm_info = nullptr;
  size_t nalm{};

  // Geometry over the SPHEREPACK Gauss grid in the production memory
  // layout: index = i_r + n_r*(i_theta + n_theta*i_phi), so for a fixed
  // radial offset the theta stride is n_r and the phi stride is n_r*n_theta.
  SharpPlan(const ylm::Spherepack& ylm, const size_t n_r) {
    const size_t nth = ylm.l_max() + 1;
    const size_t nph = 2 * ylm.m_max() + 1;
    std::vector<int> nph_arr(nth, static_cast<int>(nph));
    std::vector<ptrdiff_t> ofs(nth);
    std::vector<int> stride(nth, static_cast<int>(n_r * nth));
    std::vector<double> phi0(nth, 0.0);
    std::vector<double> theta(nth);
    std::vector<double> wgt(nth);
    for (size_t i = 0; i < nth; ++i) {
      ofs[i] = static_cast<ptrdiff_t>(i * n_r);
      theta[i] = ylm.theta_points()[i];
      // SPHEREPACK stores (2 pi / n_phi) * gauss_weight, which is exactly
      // libsharp's ring-weight convention for map2alm.
      wgt[i] = ylm.integration_weights()[i];
    }
    sharp_make_geom_info(static_cast<int>(nth), nph_arr.data(), ofs.data(),
                         stride.data(), phi0.data(), theta.data(), wgt.data(),
                         &geom);
    sharp_make_triangular_alm_info(static_cast<int>(ylm.l_max()),
                                   static_cast<int>(ylm.m_max()), 1, &alm_info);
    nalm = (ylm.l_max() + 1) * (ylm.l_max() + 2) / 2;
  }
  SharpPlan(const SharpPlan&) = delete;
  SharpPlan& operator=(const SharpPlan&) = delete;
  ~SharpPlan() {
    sharp_destroy_geom_info(geom);
    sharp_destroy_alm_info(alm_info);
  }
};

void validate_sharp(const ylm::Spherepack& ylm, std::mt19937* const gen) {
  const size_t l_max = ylm.l_max();
  const size_t nth = l_max + 1;
  const size_t nph = 2 * l_max + 1;
  const size_t npts = nth * nph;
  const SharpPlan plan(ylm, 1);

  std::vector<double> coefs(ylm.spectral_size());
  std::vector<std::complex<double>> alm(plan.nalm);
  std::vector<double> g2(npts);
  std::vector<double> gth(npts);
  std::vector<double> gph(npts);

  const auto run_case = [&](const size_t l_cut, const char* const label) {
    random_coefs(&coefs, l_max, 1, l_cut, gen);
    const DataVector g_dv =
        ylm.spec_to_phys(DataVector(coefs.data(), coefs.size()));
    void* alm_ptr[1] = {alm.data()};
    void* map_ptr[1] = {const_cast<double*>(g_dv.data())};
    sharp_execute(SHARP_MAP2ALM, 0, alm_ptr, map_ptr, plan.geom, plan.alm_info,
                  1, SHARP_DP, nullptr, nullptr);
    void* map2_ptr[1] = {g2.data()};
    sharp_execute(SHARP_ALM2MAP, 0, alm_ptr, map2_ptr, plan.geom, plan.alm_info,
                  1, SHARP_DP, nullptr, nullptr);
    std::cout << "  validate " << label << " roundtrip max_rel_diff = "
              << max_rel_diff(g_dv.data(), g2.data(), npts) << "\n";

    void* dmap_ptr[2] = {gth.data(), gph.data()};
    sharp_execute(SHARP_ALM2MAP_DERIV1, 1, alm_ptr, dmap_ptr, plan.geom,
                  plan.alm_info, 1, SHARP_DP, nullptr, nullptr);
    const auto grad = ylm.gradient(g_dv);
    std::cout << "  validate " << label << " gradient max_rel_diff = ("
              << max_rel_diff(grad.get(0).data(), gth.data(), npts) << ", "
              << max_rel_diff(grad.get(1).data(), gph.data(), npts) << ")\n";
  };
  run_case(l_max - 1, "bandlimited(l<=l_max-1)");
  // Full band limit: the gradient comparison shows SPHEREPACK's silent
  // truncation of the (l = l_max, m even) modes; libsharp computes the
  // exact band-limited gradient.
  run_case(l_max, "fullband(l<=l_max)  [expected O(1) mismatch]");
}

void run_sharp_mode(const Sizes& s, const size_t iters, const size_t warmup) {
  const ylm::Spherepack ylm(s.l_max, s.l_max);
  std::mt19937 gen(314159);
  validate_sharp(ylm, &gen);

  const SharpPlan plan(ylm, s.n_r);
  const size_t ntrans = s.n_comp * s.n_r;

  std::vector<double> u(s.npts * s.n_comp);
  for (size_t c = 0; c < s.n_comp; ++c) {
    std::vector<double> ctmp(s.spec_size * s.n_r);
    random_coefs(&ctmp, s.l_max, s.n_r, s.l_max, &gen);
    ylm.spec_to_phys_all_offsets(make_not_null(u.data() + c * s.npts),
                                 make_not_null(ctmp.data()), s.n_r);
  }
  std::vector<std::complex<double>> alm(plan.nalm * ntrans);
  std::vector<double> du(2 * s.npts * s.n_comp);

  std::vector<void*> alm_ptrs(ntrans);
  std::vector<void*> map_ptrs(ntrans);
  std::vector<void*> dmap_ptrs(2 * ntrans);
  for (size_t c = 0; c < s.n_comp; ++c) {
    for (size_t r = 0; r < s.n_r; ++r) {
      const size_t j = c * s.n_r + r;
      alm_ptrs[j] = alm.data() + j * plan.nalm;
      map_ptrs[j] = u.data() + c * s.npts + r;
      dmap_ptrs[2 * j] = du.data() + 2 * c * s.npts + r;
      dmap_ptrs[2 * j + 1] = du.data() + (2 * c + 1) * s.npts + r;
    }
  }

  // libsharp caps a single call at SHARP_MAXTRANS = 100 simultaneous
  // transforms; chunk the batch.
  const auto chunked_execute = [&](const sharp_jobtype type, const int spin,
                                   void** alm_p, void** map_p,
                                   const size_t ptrs_per_trans) {
    constexpr size_t max_trans = 100;
    for (size_t j0 = 0; j0 < ntrans; j0 += max_trans) {
      const size_t n = std::min(max_trans, ntrans - j0);
      sharp_execute(type, spin, alm_p + j0, map_p + j0 * ptrs_per_trans,
                    plan.geom, plan.alm_info, static_cast<int>(n), SHARP_DP,
                    nullptr, nullptr);
    }
  };

  const double t_analysis = time_min_us(iters, warmup, [&]() {
    chunked_execute(SHARP_MAP2ALM, 0, alm_ptrs.data(), map_ptrs.data(), 1);
  });
  print_op("sharp_map2alm", t_analysis, s);

  const double t_synthesis = time_min_us(iters, warmup, [&]() {
    chunked_execute(SHARP_ALM2MAP, 0, alm_ptrs.data(), map_ptrs.data(), 1);
  });
  print_op("sharp_alm2map", t_synthesis, s);

  const double t_deriv = time_min_us(iters, warmup, [&]() {
    chunked_execute(SHARP_ALM2MAP_DERIV1, 1, alm_ptrs.data(), dmap_ptrs.data(),
                    2);
  });
  print_op("sharp_alm2map_deriv1", t_deriv, s);

  // The production gradient path starts from nodal values, so the fair
  // comparison to SPHEREPACK's `gradient` is map2alm + alm2map_deriv1;
  // the filter round trip is map2alm + alm2map.
  print_op("sharp_gradient(map2alm+deriv1)", t_analysis + t_deriv, s);
  print_op("sharp_filter(map2alm+alm2map)", t_analysis + t_synthesis, s);
}

// --- mode 2: fused per-m theta matrices via dgemm + hrfft ----------------

class FusedMatrixTransform {
 public:
  FusedMatrixTransform(const size_t l_max, const size_t n_r)
      : l_(l_max),
        nth_(l_max + 1),
        nph_(2 * l_max + 1),
        nr_(n_r),
        wsave_(2 * nph_ + 15),
        fft_work_(nr_ * nth_ * nph_),
        g_(nr_ * nth_ * nph_),
        h_(nr_ * nth_ * nph_),
        // Analysis A[m] is (nth x (l+1-m)); synthesis S[m] is ((l+1-m) x
        // nth); gradient D_v[m], D_w[m] are (nth x nth), stored row-major.
        a_mats_(nth_ * total_coef_cols()),
        s_mats_(nth_ * total_coef_cols()),
        dv_mats_((l_ + 1) * nth_ * nth_),
        dw_mats_((l_ + 1) * nth_ * nth_),
        coef_((l_ + 1) * (l_ + 1) * nr_ * 2) {
    const int nph_int = static_cast<int>(nph_);
    hrffti_(nph_int, wsave_.data());
    build_matrices();
  }

  size_t table_bytes() const {
    return 8 * (a_mats_.size() + s_mats_.size() + dv_mats_.size() +
                dw_mats_.size());
  }
  size_t transform_table_bytes() const {
    return 8 * (a_mats_.size() + s_mats_.size());
  }

  // Analysis + synthesis round trip on one component slab (in -> out); used
  // both for validation and as the filter-cost analog.
  void roundtrip(double* const out, const double* const in) {
    analysis(in);
    synthesis(out);
  }

  void analysis(const double* const in) {
    const int batch = static_cast<int>(nr_ * nth_);
    const int nph_int = static_cast<int>(nph_);
    std::copy(in, in + nr_ * nth_ * nph_, g_.data());
    hrfftf_(batch, nph_int, g_.data(), batch, wsave_.data(), fft_work_.data());
    // m = 0: single (cos) column; m >= 1: (re, im) columns.
    size_t coef_offset = 0;
    for (size_t m = 0; m <= l_; ++m) {
      const size_t ncols = l_ + 1 - m;
      const double* const amat = a_mats_.data() + a_offset(m);
      const size_t ncopies = (m == 0 ? 1 : 2);
      for (size_t p = 0; p < ncopies; ++p) {
        const double* const gm =
            g_.data() + (m == 0 ? 0 : (2 * m - 1 + p)) * nr_ * nth_;
        dgemm_('N', 'N', nr_, ncols, nth_, 1.0, gm, nr_, amat, nth_, 0.0,
               coef_.data() + coef_offset, nr_);
        coef_offset += nr_ * ncols;
      }
    }
  }

  void synthesis(double* const out) {
    const int batch = static_cast<int>(nr_ * nth_);
    const int nph_int = static_cast<int>(nph_);
    size_t coef_offset = 0;
    for (size_t m = 0; m <= l_; ++m) {
      const size_t ncols = l_ + 1 - m;
      const double* const smat = s_mats_.data() + a_offset(m);
      const size_t ncopies = (m == 0 ? 1 : 2);
      for (size_t p = 0; p < ncopies; ++p) {
        double* const hm =
            h_.data() + (m == 0 ? 0 : (2 * m - 1 + p)) * nr_ * nth_;
        // S stored as (nth x ncols) row-major-transposed: we need
        // h(r, i) = sum_n coef(r, n) * Pbar_n(theta_i); use the same table
        // as analysis without the weights, so multiply by S^T via 'T'.
        dgemm_('N', 'T', nr_, nth_, ncols, 1.0, coef_.data() + coef_offset, nr_,
               smat, nth_, 0.0, hm, nr_);
        coef_offset += nr_ * ncols;
      }
    }
    hrfftb_(batch, nph_int, h_.data(), batch, wsave_.data(), fft_work_.data());
    std::copy(h_.begin(), h_.end(), out);
  }

  void gradient(double* const du_theta, double* const du_phi,
                const double* const in) {
    const int batch = static_cast<int>(nr_ * nth_);
    const int nph_int = static_cast<int>(nph_);
    const size_t npts = nr_ * nth_ * nph_;
    std::copy(in, in + npts, g_.data());
    hrfftf_(batch, nph_int, g_.data(), batch, wsave_.data(), fft_work_.data());

    // m = 0.
    apply_theta(du_theta, g_.data(), dv_mats_.data(), 1.0);
    std::fill(du_phi, du_phi + nr_ * nth_, 0.0);
    for (size_t m = 1; m <= l_; ++m) {
      const size_t col = (2 * m - 1) * nr_ * nth_;
      const double* const dv = dv_mats_.data() + m * nth_ * nth_;
      const double* const dw = dw_mats_.data() + m * nth_ * nth_;
      apply_theta(du_theta + col, g_.data() + col, dv, 1.0);
      apply_theta(du_theta + col + nr_ * nth_, g_.data() + col + nr_ * nth_, dv,
                  1.0);
      // w_re = -D_w g_im ; w_im = +D_w g_re
      apply_theta(du_phi + col, g_.data() + col + nr_ * nth_, dw, -1.0);
      apply_theta(du_phi + col + nr_ * nth_, g_.data() + col, dw, 1.0);
    }
    hrfftb_(batch, nph_int, du_theta, batch, wsave_.data(), fft_work_.data());
    hrfftb_(batch, nph_int, du_phi, batch, wsave_.data(), fft_work_.data());
  }

 private:
  size_t total_coef_cols() const {
    // sum_m (l+1-m) for m = 0..l
    return (l_ + 1) * (l_ + 2) / 2;
  }
  size_t a_offset(const size_t m) const {
    // Offset of A[m] in units of doubles: nth * sum_{k<m} (l+1-k).
    return nth_ * (m * (l_ + 1) - m * (m - 1) / 2);
  }

  // out(r, i) = sign * sum_j D(i, j) g(r, j): D row-major (i, j) is
  // column-major (j, i), i.e. already D^T, so plain 'N','N' dgemm.
  void apply_theta(double* const out, const double* const in,
                   const double* const dmat, const double sign) {
    dgemm_('N', 'N', nr_, nth_, nth_, sign, in, nr_, dmat, nth_, 0.0, out, nr_);
  }

  void build_matrices() {
    std::vector<double> theta(nth_ + 1);
    std::vector<double> wts(nth_ + 1);
    std::vector<double> work(nth_ + 1);
    int err = 0;
    gaqd_(static_cast<int>(nth_), theta.data(), wts.data(), work.data(),
          static_cast<int>(work.size()), make_not_null(&err));
    if (err != 0) {
      std::cerr << "gaqd error " << err << "\n";
      std::exit(1);
    }

    std::vector<double> pbar(nth_ * (l_ + 1));
    std::vector<double> dpbar(nth_ * (l_ + 1));
    std::vector<double> x(nth_), sx(nth_);
    for (size_t j = 0; j < nth_; ++j) {
      x[j] = std::cos(theta[j]);
      sx[j] = std::sin(theta[j]);
    }

    for (size_t m = 0; m <= l_; ++m) {
      for (size_t j = 0; j < nth_; ++j) {
        double pmm = M_SQRT1_2;
        for (size_t k = 1; k <= m; ++k) {
          pmm *= sx[j] * std::sqrt((2.0 * k + 1.0) / (2.0 * k));
        }
        pbar[j + m * nth_] = pmm;
      }
      if (m + 1 <= l_) {
        for (size_t j = 0; j < nth_; ++j) {
          pbar[j + (m + 1) * nth_] =
              std::sqrt(2.0 * m + 3.0) * x[j] * pbar[j + m * nth_];
        }
      }
      for (size_t n = m + 2; n <= l_; ++n) {
        const double nn = static_cast<double>(n);
        const double mm = static_cast<double>(m);
        const double alpha =
            std::sqrt((4.0 * nn * nn - 1.0) / (nn * nn - mm * mm));
        const double beta = std::sqrt(
            ((2.0 * nn + 1.0) / (2.0 * nn - 3.0)) *
            (((nn - 1.0) * (nn - 1.0) - mm * mm) / (nn * nn - mm * mm)));
        for (size_t j = 0; j < nth_; ++j) {
          pbar[j + n * nth_] = alpha * x[j] * pbar[j + (n - 1) * nth_] -
                               beta * pbar[j + (n - 2) * nth_];
        }
      }
      for (size_t n = m; n <= l_; ++n) {
        const double nn = static_cast<double>(n);
        const double mm = static_cast<double>(m);
        const double e_n = n == m
                               ? 0.0
                               : std::sqrt((nn * nn - mm * mm) *
                                           (2.0 * nn + 1.0) / (2.0 * nn - 1.0));
        for (size_t j = 0; j < nth_; ++j) {
          const double prev = n == m ? 0.0 : pbar[j + (n - 1) * nth_];
          dpbar[j + n * nth_] =
              (nn * x[j] * pbar[j + n * nth_] - e_n * prev) / sx[j];
        }
      }

      const double scale = 1.0 / static_cast<double>(nph_);
      // Analysis table A[m](j, n) = scale * w_j * Pbar_n(theta_j), stored
      // column-major (nth x ncols); synthesis table S[m](j, n) =
      // Pbar_n(theta_j) in the same layout (used transposed).
      double* const amat = a_mats_.data() + a_offset(m);
      double* const smat = s_mats_.data() + a_offset(m);
      for (size_t n = m; n <= l_; ++n) {
        for (size_t j = 0; j < nth_; ++j) {
          amat[j + (n - m) * nth_] = scale * wts[j] * pbar[j + n * nth_];
          smat[j + (n - m) * nth_] = pbar[j + n * nth_];
        }
      }

      // Fused gradient matrices, replicating SPHEREPACK's top-degree
      // truncation (vhsgs1 drops (l = l_max, m even) from the gradient).
      double* const dv = dv_mats_.data() + m * nth_ * nth_;
      double* const dw = dw_mats_.data() + m * nth_ * nth_;
      for (size_t i = 0; i < nth_; ++i) {
        for (size_t j = 0; j < nth_; ++j) {
          double acc_v = 0.0;
          double acc_w = 0.0;
          for (size_t n = m; n <= l_; ++n) {
            if (n == l_ and m % 2 == 0) {
              continue;
            }
            acc_v += dpbar[i + n * nth_] * pbar[j + n * nth_];
            acc_w += pbar[i + n * nth_] * pbar[j + n * nth_];
          }
          dv[i * nth_ + j] = scale * acc_v * wts[j];
          dw[i * nth_ + j] =
              scale * static_cast<double>(m) / sx[i] * acc_w * wts[j];
        }
      }
    }
  }

  size_t l_, nth_, nph_, nr_;
  std::vector<double> wsave_, fft_work_, g_, h_;
  std::vector<double> a_mats_, s_mats_, dv_mats_, dw_mats_;
  std::vector<double> coef_;
};

void run_fused_mode(const Sizes& s, const size_t iters, const size_t warmup) {
  const ylm::Spherepack ylm(s.l_max, s.l_max);
  std::mt19937 gen(314159);
  FusedMatrixTransform fused(s.l_max, s.n_r);

  std::vector<double> u(s.npts * s.n_comp);
  for (size_t c = 0; c < s.n_comp; ++c) {
    std::vector<double> ctmp(s.spec_size * s.n_r);
    random_coefs(&ctmp, s.l_max, s.n_r, s.l_max, &gen);
    ylm.spec_to_phys_all_offsets(make_not_null(u.data() + c * s.npts),
                                 make_not_null(ctmp.data()), s.n_r);
  }
  std::vector<double> u2(s.npts * s.n_comp);
  std::vector<double> dth(s.npts * s.n_comp);
  std::vector<double> dph(s.npts * s.n_comp);
  std::vector<double> dth_ref(s.npts);
  std::vector<double> dph_ref(s.npts);

  // Validation: round trip reproduces the band-limited input; gradient
  // matches SPHEREPACK (including its top-degree truncation).
  fused.roundtrip(u2.data(), u.data());
  std::cout << "  validate roundtrip max_rel_diff = "
            << max_rel_diff(u.data(), u2.data(), s.npts) << "\n";
  fused.gradient(dth.data(), dph.data(), u.data());
  {
    const std::array<double*, 2> df{dth_ref.data(), dph_ref.data()};
    ylm.gradient_all_offsets(df, make_not_null(u.data()), s.n_r);
  }
  std::cout << "  validate gradient max_rel_diff = ("
            << max_rel_diff(dth.data(), dth_ref.data(), s.npts) << ", "
            << max_rel_diff(dph.data(), dph_ref.data(), s.npts) << ")\n";

  const double t_analysis = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      fused.analysis(u.data() + c * s.npts);
    }
  });
  print_op("fused_analysis", t_analysis, s);

  const double t_synthesis = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      fused.synthesis(u2.data() + c * s.npts);
    }
  });
  print_op("fused_synthesis", t_synthesis, s);

  const double t_roundtrip = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      fused.roundtrip(u2.data() + c * s.npts, u.data() + c * s.npts);
    }
  });
  print_op("fused_filter_roundtrip", t_roundtrip, s);

  const double t_gradient = time_min_us(iters, warmup, [&]() {
    for (size_t c = 0; c < s.n_comp; ++c) {
      fused.gradient(dth.data() + c * s.npts, dph.data() + c * s.npts,
                     u.data() + c * s.npts);
    }
  });
  print_op("fused_gradient", t_gradient, s);

  std::cout << "  fused tables: " << fused.table_bytes()
            << " bytes total (transform pair " << fused.transform_table_bytes()
            << ")\n";
}

// --- mode 3: arbitrary-point evaluation (Shape / horizon path) -----------

void run_interp_mode(const Sizes& s, const size_t iters, const size_t warmup) {
  const ylm::Spherepack ylm(s.l_max, s.l_max);
  std::mt19937 gen(314159);
  constexpr size_t num_points = 1000;
  constexpr size_t num_fields = 5;

  std::uniform_real_distribution<double> theta_dist(0.05, M_PI - 0.05);
  std::uniform_real_distribution<double> phi_dist(0.0, 2.0 * M_PI);
  std::array<DataVector, 2> points{DataVector(num_points),
                                   DataVector(num_points)};
  for (size_t i = 0; i < num_points; ++i) {
    points[0][i] = theta_dist(gen);
    points[1][i] = phi_dist(gen);
  }

  std::vector<std::vector<double>> field_coefs(num_fields);
  for (auto& coefs : field_coefs) {
    coefs.resize(ylm.spectral_size());
    random_coefs(&coefs, s.l_max, 1, s.l_max, &gen);
  }

  const double t_setup = time_min_us(iters, warmup, [&]() {
    const auto info = ylm.set_up_interpolation_info(points);
    (void)info;
  });
  std::cout << "  interp_setup(" << num_points << " pts): " << t_setup
            << " us\n";

  const auto info = ylm.set_up_interpolation_info(points);
  DataVector result(num_points);
  const double t_clenshaw = time_min_us(iters, warmup, [&]() {
    for (size_t f = 0; f < num_fields; ++f) {
      ylm.interpolate_from_coefs(
          make_not_null(&result),
          DataVector(field_coefs[f].data(), field_coefs[f].size()), info);
    }
  });
  std::cout << "  interp_clenshaw(" << num_points << " pts x " << num_fields
            << " fields): " << t_clenshaw << " us\n";

  // Per-point calls with T = double: the Shape::inverse /
  // block_logical_coordinates pattern (one target point per call, setup +
  // Clenshaw every time).
  const double t_per_point = time_min_us(iters, warmup, [&]() {
    double acc = 0.0;
    for (size_t i = 0; i < num_points; ++i) {
      const std::array<double, 2> pt{points[0][i], points[1][i]};
      acc += ylm.interpolate_from_coefs<double>(
          DataVector(field_coefs[0].data(), field_coefs[0].size()), pt);
    }
    if (acc == 0.12345) {
      std::cout << "";  // defeat dead-code elimination
    }
  });
  std::cout << "  interp_per_point(" << num_points
            << " single-point calls, 1 field): " << t_per_point << " us\n";

  // Cached dense-matrix alternative: result = B * C with B (num_points x
  // n_modes) built once per (grid, l_max) and C the (n_modes x num_fields)
  // coefficients. Timing only; values are irrelevant to the cost.
  const size_t n_modes = (s.l_max + 1) * (s.l_max + 1);
  std::vector<double> basis_matrix(num_points * n_modes);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto& v : basis_matrix) {
    v = dist(gen);
  }
  std::vector<double> coef_matrix(n_modes * num_fields);
  for (auto& v : coef_matrix) {
    v = dist(gen);
  }
  std::vector<double> out(num_points * num_fields);
  const double t_matrix = time_min_us(iters, warmup, [&]() {
    dgemm_('N', 'N', num_points, num_fields, n_modes, 1.0, basis_matrix.data(),
           num_points, coef_matrix.data(), n_modes, 0.0, out.data(),
           num_points);
  });
  std::cout << "  interp_matrix(" << num_points << " pts x " << num_fields
            << " fields dgemm): " << t_matrix
            << " us  (B = " << 8 * num_points * n_modes << " bytes)\n";
}

}  // namespace

int main(const int argc, const char* const* const argv) {
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0]
              << " <l_max> <n_r> <n_comp> <iters> <warmup> <mode>\n";
    return 1;
  }
  const size_t l_max = std::stoul(argv[1]);
  const size_t n_r = std::stoul(argv[2]);
  const size_t n_comp = std::stoul(argv[3]);
  const size_t iters = std::stoul(argv[4]);
  const size_t warmup = std::stoul(argv[5]);
  const size_t mode = std::stoul(argv[6]);

  Sizes s{};
  s.l_max = l_max;
  s.n_r = n_r;
  s.n_comp = n_comp;
  s.nth = l_max + 1;
  s.nph = 2 * l_max + 1;
  s.npts = n_r * s.nth * s.nph;
  s.spec_size = ylm::Spherepack::spectral_size(l_max, l_max);

  std::cout << "l_max=" << l_max << " n_r=" << n_r << " n_comp=" << n_comp
            << " grid=" << s.nth << "x" << s.nph << " npts/comp=" << s.npts
            << " mode=" << mode << "\n";

  switch (mode) {
    case 0:
      run_spherepack_mode(s, iters, warmup);
      break;
    case 1:
      run_sharp_mode(s, iters, warmup);
      break;
    case 2:
      run_fused_mode(s, iters, warmup);
      break;
    case 3:
      run_interp_mode(s, iters, warmup);
      break;
    default:
      std::cerr << "Unknown mode " << mode << "\n";
      return 1;
  }
  return 0;
}
