// Distributed under the MIT License.
// See LICENSE.txt for details.

// Benchmark harness for the spectral derivative kernels at production sizes
// (GH-like Variables with 50 components).
//
// Usage: PartialDerivativesBenchmark <pts_per_dim> <iters> <warmup> <mode>
//   mode 0 = logical_partial_derivatives
//   mode 1 = partial_derivatives (logical derivatives + Jacobian contraction)
//   mode 2 = exponential filter (apply_matrices along all dimensions)
//   mode 3 = PROTOTYPE: transpose-free fixed-size derivative kernels with the
//            Jacobian contraction fused per component (32 KB working set).
//            Prints max relative difference against mode 1 before timing.
//   mode 4 = logical_partial_derivatives on a SphericalHarmonic-basis shell
//            mesh. Pass extents as "n_r x n_theta x n_phi" with
//            n_theta = l_max+1 and n_phi = 2*l_max+1, e.g. 12x13x25.
//   mode 5 = partial_derivatives on the SphericalHarmonic shell mesh
//   mode 6 = component decomposition of the SphericalHarmonic derivative:
//            radial dgemm, SPHEREPACK analysis (shags), gradient synthesis
//            (gradgs), and the combined gradient, each timed separately.
//   mode 7 = PROTOTYPE: fused per-m angular derivative. The SPHEREPACK chain
//            (scalar analysis -> coefficient scaling -> vector synthesis)
//            collapses, per Fourier mode m, into two small theta-matrices
//            D_v[m] = dPbar * W * Pbar^T and D_w[m] = m csc(theta) Pbar * W *
//            Pbar^T applied to the (re, im) columns, batched over the
//            contiguous radial index. Phi transforms stay hrfftf/hrfftb.
//            Prints max relative difference against the SPHEREPACK gradient
//            before timing.

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/DifferentiationMatrix.hpp"
#include "NumericalAlgorithms/Spectral/Filtering.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Spherepack.hpp"
#include "Utilities/TMPL.hpp"

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
constexpr size_t Dim = 3;

struct Metric : db::SimpleTag {
  using type = tnsr::aa<DataVector, Dim, Frame::Inertial>;
};
struct Pi : db::SimpleTag {
  using type = tnsr::aa<DataVector, Dim, Frame::Inertial>;
};
struct Phi : db::SimpleTag {
  using type = tnsr::iaa<DataVector, Dim, Frame::Inertial>;
};

using vars_tags = tmpl::list<Metric, Pi, Phi>;
using deriv_tags = db::wrap_tags_in<::Tags::deriv, vars_tags, tmpl::size_t<Dim>,
                                    Frame::Inertial>;

// --- Prototype kernels -------------------------------------------------
//
// All three kernels differentiate one component (N^3 points, x fastest) with
// the (column-major, unpadded) differentiation matrix `dmat`, writing plain
// unit-stride FMA loops that auto-vectorize on NEON here and AVX-512 on the
// x86 targets. No transposes, no BLAS packing.

// out(i,s) = sum_k D(i,k) u(k,s) for each of the N^2 contiguous x-lines s.
// Vectorized over the output line held in registers.
template <size_t N>
void dxi_kernel(double* const out, const double* const uc,
                const double* const dmat) {
  for (size_t s = 0; s < N * N; ++s) {
    const double* const line = uc + s * N;
    double* const o = out + s * N;
    std::array<double, N> acc{};
    for (size_t k = 0; k < N; ++k) {
      const double uk = line[k];
      for (size_t i = 0; i < N; ++i) {
        acc[i] += dmat[i + N * k] * uk;
      }
    }
    for (size_t i = 0; i < N; ++i) {
      o[i] = acc[i];
    }
  }
}

// out(x,i,z) = sum_k D(i,k) u(x,k,z). Vectorized over the contiguous x-line.
template <size_t N>
void deta_kernel(double* const out, const double* const uc,
                 const double* const dmat) {
  for (size_t z = 0; z < N; ++z) {
    const double* const slab = uc + z * N * N;
    double* const oslab = out + z * N * N;
    for (size_t i = 0; i < N; ++i) {
      double* const o = oslab + i * N;
      std::array<double, N> acc{};
      for (size_t k = 0; k < N; ++k) {
        const double d = dmat[i + N * k];
        const double* const line = slab + k * N;
        for (size_t x = 0; x < N; ++x) {
          acc[x] += d * line[x];
        }
      }
      for (size_t x = 0; x < N; ++x) {
        o[x] = acc[x];
      }
    }
  }
}

// out(xy,i) = sum_k D(i,k) u(xy,k). Register-blocked over xy-chunks of N so
// each output value is accumulated in registers and stored exactly once
// (a single read-modify-write pass over the output plane per k was ~1.6x
// slower).
template <size_t N>
void dzeta_kernel(double* const out, const double* const uc,
                  const double* const dmat) {
  for (size_t chunk = 0; chunk < N * N; chunk += N) {
    const double* const ublock = uc + chunk;
    for (size_t i = 0; i < N; ++i) {
      std::array<double, N> acc{};
      for (size_t k = 0; k < N; ++k) {
        const double d = dmat[i + N * k];
        const double* const plane = ublock + k * N * N;
        for (size_t x = 0; x < N; ++x) {
          acc[x] += d * plane[x];
        }
      }
      double* const o = out + i * N * N + chunk;
      for (size_t x = 0; x < N; ++x) {
        o[x] = acc[x];
      }
    }
  }
}

// Component-at-a-time derivative with fused Jacobian contraction. The three
// logical-derivative buffers (3 x 8 KB at N=10) plus the component itself stay
// L1-resident; the logical derivatives are never written to main memory.
// Output layout matches partial_derivatives: [component][deriv index][points].
// jac[d * Dim + i] points to inverse_jacobian.get(d, i).
template <size_t N>
void prototype_apply(double* const du, const double* const u,
                     const size_t number_of_components,
                     const double* const dmat,
                     const std::array<const double*, Dim * Dim>& jac,
                     double* const buffer) {
  constexpr size_t npts = N * N * N;
  double* const b0 = buffer;
  double* const b1 = buffer + npts;
  double* const b2 = buffer + 2 * npts;
  for (size_t c = 0; c < number_of_components; ++c) {
    const double* const uc = u + c * npts;
    dxi_kernel<N>(b0, uc, dmat);
    deta_kernel<N>(b1, uc, dmat);
    dzeta_kernel<N>(b2, uc, dmat);
    for (size_t i = 0; i < Dim; ++i) {
      double* const o = du + (c * Dim + i) * npts;
      const double* const j0 = jac[0 * Dim + i];
      const double* const j1 = jac[1 * Dim + i];
      const double* const j2 = jac[2 * Dim + i];
      for (size_t p = 0; p < npts; ++p) {
        o[p] = j0[p] * b0[p] + j1[p] * b1[p] + j2[p] * b2[p];
      }
    }
  }
}

// --- Spherical-harmonic angular-derivative prototype (mode 7) ------------
//
// Grid: (n_r, n_theta, n_phi) with r fastest, theta Gauss points, phi
// equiangular; n_theta = l_max+1, n_phi = 2*l_max+1. Per component:
//   1. forward real FFT along phi (batch n_r*n_theta, hrfftf)
//   2. for each Fourier mode m, apply the fused theta-matrices to the
//      (re, im) columns:  v = D_v[m] g  (d/dtheta) and
//      w_re = -D_w[m] g_im, w_im = +D_w[m] g_re  (csc(theta) d/dphi)
//   3. backward real FFTs of v and w (hrfftb)
// The matrices absorb the quadrature weights, the FFT round-trip factor
// 1/n_phi, and the factor m of the phi derivative. Normalization of the
// Legendre functions cancels in the product dPbar * W * Pbar^T, so the
// tables are built with standard fully-normalized recurrences and validated
// against SPHEREPACK's gradient to roundoff.
class SphericalGradientPrototype {
 public:
  SphericalGradientPrototype(const size_t l_max, const size_t n_r)
      : l_(l_max),
        nth_(l_max + 1),
        nph_(2 * l_max + 1),
        nr_(n_r),
        wsave_(2 * nph_ + 15),
        fft_work_(nr_ * nth_ * nph_),
        g_(nr_ * nth_ * nph_),
        dv_mats_((l_ + 1) * nth_ * nth_),
        dw_mats_((l_ + 1) * nth_ * nth_) {
    const int nph_int = static_cast<int>(nph_);
    hrffti_(nph_int, wsave_.data());
    // hrfft only has fast kernels for radices 2, 3, 4, 5; a PRIME n_phi
    // makes it fall through to a single O(n^2) generic-radix pass, which the
    // batch-vectorized folded real-DFT kernels beat. For composite
    // non-smooth sizes (21, 33, 49, ...) hrfft's mixed radix still wins.
    bool is_prime = nph_ > 1;
    for (size_t p = 2; p * p <= nph_; ++p) {
      if (nph_ % p == 0) {
        is_prime = false;
        break;
      }
    }
    use_dft_kernels_ = is_prime and (nph_ - 1) / 2 <= 32;
    if (use_dft_kernels_) {
      const size_t h = (nph_ - 1) / 2;
      dft_cos_.resize(h * h);
      dft_sin_.resize(h * h);
      const double wphi = 2.0 * M_PI / static_cast<double>(nph_);
      for (size_t k = 1; k <= h; ++k) {
        for (size_t i = 1; i <= h; ++i) {
          dft_cos_[(k - 1) * h + (i - 1)] =
              std::cos(wphi * static_cast<double>(k * i));
          dft_sin_[(k - 1) * h + (i - 1)] =
              std::sin(wphi * static_cast<double>(k * i));
        }
      }
      vhat_.resize(nr_ * nth_ * nph_);
      what_.resize(nr_ * nth_ * nph_);
    }
    build_matrices();
  }

  // du_theta, du_phi, u: pointers to one component's (n_r*n_theta*n_phi)
  // slab.
  void apply(double* const du_theta, double* const du_phi,
             const double* const u) {
    const int batch = static_cast<int>(nr_ * nth_);
    const int nph_int = static_cast<int>(nph_);
    const size_t npts = nr_ * nth_ * nph_;
    if (use_dft_kernels_) {
      forward_dft(g_.data(), u);
    } else {
      std::copy(u, u + npts, g_.data());
      hrfftf_(batch, nph_int, g_.data(), batch, wsave_.data(),
              fft_work_.data());
    }

    double* const vdst = use_dft_kernels_ ? vhat_.data() : du_theta;
    double* const wdst = use_dft_kernels_ ? what_.data() : du_phi;
    // m = 0: single column j=0. dv output column, dw output is zero.
    apply_theta_matrix(vdst, g_.data(), dv_mat(0));
    std::fill(wdst, wdst + nr_ * nth_, 0.0);
    for (size_t m = 1; m <= l_; ++m) {
      double* const v_re = vdst + (2 * m - 1) * nr_ * nth_;
      double* const v_im = vdst + (2 * m) * nr_ * nth_;
      double* const w_re = wdst + (2 * m - 1) * nr_ * nth_;
      double* const w_im = wdst + (2 * m) * nr_ * nth_;
      const double* const g_re = g_.data() + (2 * m - 1) * nr_ * nth_;
      const double* const g_im = g_.data() + (2 * m) * nr_ * nth_;
      apply_theta_matrix(v_re, g_re, dv_mat(m));
      apply_theta_matrix(v_im, g_im, dv_mat(m));
      apply_theta_matrix_negated(w_re, g_im, dw_mat(m));
      apply_theta_matrix(w_im, g_re, dw_mat(m));
    }
    if (use_dft_kernels_) {
      backward_dft(du_theta, vhat_.data());
      backward_dft(du_phi, what_.data());
    } else {
      hrfftb_(batch, nph_int, du_theta, batch, wsave_.data(),
              fft_work_.data());
      hrfftb_(batch, nph_int, du_phi, batch, wsave_.data(), fft_work_.data());
    }
  }

 private:
  const double* dv_mat(const size_t m) const {
    return dv_mats_.data() + m * nth_ * nth_;
  }
  const double* dw_mat(const size_t m) const {
    return dw_mats_.data() + m * nth_ * nth_;
  }

  // out(ig, i) = Sign * sum_j D(i, j) in(ig, j), batch ig contiguous
  // (size NR = n_r, compile time). D is (nth x nth) row-major. Output rows
  // are accumulated in registers and stored once (two rows at a time to
  // reuse the loaded input lines).
  template <size_t NR, int Sign>
  static void theta_kernel(double* const out, const double* const in,
                           const double* const dmat, const size_t nth) {
    size_t i = 0;
    // Two-row blocking halves the input-line loads but doubles the register
    // accumulators; past NR ~ 12 it spills and single-row wins.
    constexpr bool two_rows = NR <= 12;
    for (; two_rows and i + 2 <= nth; i += 2) {
      std::array<double, NR> acc0{};
      std::array<double, NR> acc1{};
      const double* const row0 = dmat + i * nth;
      const double* const row1 = dmat + (i + 1) * nth;
      for (size_t j = 0; j < nth; ++j) {
        const double d0 = row0[j];
        const double d1 = row1[j];
        const double* const line = in + j * NR;
        for (size_t ig = 0; ig < NR; ++ig) {
          acc0[ig] += d0 * line[ig];
          acc1[ig] += d1 * line[ig];
        }
      }
      double* const o0 = out + i * NR;
      double* const o1 = out + (i + 1) * NR;
      for (size_t ig = 0; ig < NR; ++ig) {
        o0[ig] = Sign * acc0[ig];
        o1[ig] = Sign * acc1[ig];
      }
    }
    for (; i < nth; ++i) {
      std::array<double, NR> acc{};
      const double* const row = dmat + i * nth;
      for (size_t j = 0; j < nth; ++j) {
        const double d = row[j];
        const double* const line = in + j * NR;
        for (size_t ig = 0; ig < NR; ++ig) {
          acc[ig] += d * line[ig];
        }
      }
      double* const o = out + i * NR;
      for (size_t ig = 0; ig < NR; ++ig) {
        o[ig] = Sign * acc[ig];
      }
    }
  }

  template <int Sign>
  void apply_theta_matrix_impl(double* const out, const double* const in,
                               const double* const dmat) const {
    switch (nr_) {
#define SPH_PROTO_CASE(NR)                            \
  case (NR):                                          \
    theta_kernel<(NR), Sign>(out, in, dmat, nth_);    \
    break
      SPH_PROTO_CASE(2);
      SPH_PROTO_CASE(3);
      SPH_PROTO_CASE(4);
      SPH_PROTO_CASE(5);
      SPH_PROTO_CASE(6);
      SPH_PROTO_CASE(7);
      SPH_PROTO_CASE(8);
      SPH_PROTO_CASE(9);
      SPH_PROTO_CASE(10);
      SPH_PROTO_CASE(11);
      SPH_PROTO_CASE(12);
      SPH_PROTO_CASE(13);
      SPH_PROTO_CASE(14);
      SPH_PROTO_CASE(15);
      SPH_PROTO_CASE(16);
      SPH_PROTO_CASE(17);
      SPH_PROTO_CASE(18);
      SPH_PROTO_CASE(19);
      SPH_PROTO_CASE(20);
#undef SPH_PROTO_CASE
      default:
        std::cerr << "Prototype needs n_r in 2..20, got " << nr_ << "\n";
        std::exit(1);
    }
  }

  void apply_theta_matrix(double* const out, const double* const in,
                          const double* const dmat) const {
    apply_theta_matrix_impl<1>(out, in, dmat);
  }
  void apply_theta_matrix_negated(double* const out, const double* const in,
                                  const double* const dmat) const {
    apply_theta_matrix_impl<-1>(out, in, dmat);
  }

  // Folded real DFT along phi, hrfftf-compatible halfcomplex layout, n odd.
  // Data is (batch, n) column-major, batch = n_r*n_theta contiguous. With
  // u_i = x_i + x_{n-i}, v_i = x_i - x_{n-i} (i = 1..H, H = (n-1)/2):
  //   out(:,0)    = x_0 + sum_i u_i
  //   out(:,2k-1) = x_0 + sum_i u_i cos(2 pi k i / n)
  //   out(:,2k)   =     - sum_i v_i sin(2 pi k i / n)
  // Processed in strips of W batch lanes; output (re, im) pairs accumulate
  // in registers over i.
  static constexpr size_t dft_strip = 8;

  // Strip worker with a compile-time lane count so the accumulators live in
  // registers (a runtime lane bound blocks register promotion entirely).
  template <size_t LANES>
  void forward_dft_strip(double* const out, const double* const in,
                         const size_t b0) const {
    const size_t h = (nph_ - 1) / 2;
    const size_t batch = nr_ * nth_;
    std::array<double, 32 * LANES> u_fold;
    std::array<double, 32 * LANES> v_fold;
    const double* const x0 = in + b0;
    for (size_t i = 1; i <= h; ++i) {
      const double* const xi = in + i * batch + b0;
      const double* const xni = in + (nph_ - i) * batch + b0;
      double* const uf = u_fold.data() + (i - 1) * LANES;
      double* const vf = v_fold.data() + (i - 1) * LANES;
      for (size_t lane = 0; lane < LANES; ++lane) {
        uf[lane] = xi[lane] + xni[lane];
        vf[lane] = xi[lane] - xni[lane];
      }
    }
    {
      std::array<double, LANES> acc{};
      for (size_t i = 0; i < h; ++i) {
        const double* const uf = u_fold.data() + i * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          acc[lane] += uf[lane];
        }
      }
      for (size_t lane = 0; lane < LANES; ++lane) {
        out[b0 + lane] = x0[lane] + acc[lane];
      }
    }
    size_t k = 1;
    for (; k + 1 <= h; k += 2) {
      std::array<double, LANES> re0{};
      std::array<double, LANES> im0{};
      std::array<double, LANES> re1{};
      std::array<double, LANES> im1{};
      const double* const c0 = dft_cos_.data() + (k - 1) * h;
      const double* const s0 = dft_sin_.data() + (k - 1) * h;
      const double* const c1 = dft_cos_.data() + k * h;
      const double* const s1 = dft_sin_.data() + k * h;
      for (size_t i = 0; i < h; ++i) {
        const double* const uf = u_fold.data() + i * LANES;
        const double* const vf = v_fold.data() + i * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          re0[lane] += c0[i] * uf[lane];
          im0[lane] -= s0[i] * vf[lane];
          re1[lane] += c1[i] * uf[lane];
          im1[lane] -= s1[i] * vf[lane];
        }
      }
      double* const ore0 = out + (2 * k - 1) * batch + b0;
      double* const oim0 = out + (2 * k) * batch + b0;
      double* const ore1 = out + (2 * k + 1) * batch + b0;
      double* const oim1 = out + (2 * k + 2) * batch + b0;
      for (size_t lane = 0; lane < LANES; ++lane) {
        ore0[lane] = x0[lane] + re0[lane];
        oim0[lane] = im0[lane];
        ore1[lane] = x0[lane] + re1[lane];
        oim1[lane] = im1[lane];
      }
    }
    for (; k <= h; ++k) {
      std::array<double, LANES> re{};
      std::array<double, LANES> im{};
      const double* const ck = dft_cos_.data() + (k - 1) * h;
      const double* const sk = dft_sin_.data() + (k - 1) * h;
      for (size_t i = 0; i < h; ++i) {
        const double* const uf = u_fold.data() + i * LANES;
        const double* const vf = v_fold.data() + i * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          re[lane] += ck[i] * uf[lane];
          im[lane] -= sk[i] * vf[lane];
        }
      }
      double* const ore = out + (2 * k - 1) * batch + b0;
      double* const oim = out + (2 * k) * batch + b0;
      for (size_t lane = 0; lane < LANES; ++lane) {
        ore[lane] = x0[lane] + re[lane];
        oim[lane] = im[lane];
      }
    }
  }

  void forward_dft(double* const out, const double* const in) const {
    const size_t batch = nr_ * nth_;
    constexpr size_t W = dft_strip;
    size_t b0 = 0;
    for (; b0 + W <= batch; b0 += W) {
      forward_dft_strip<W>(out, in, b0);
    }
    for (; b0 < batch; ++b0) {
      forward_dft_strip<1>(out, in, b0);
    }
  }

  // Inverse of hrfftf's layout (hrfftb-compatible, unnormalized):
  //   x_0     = c_0 + 2 sum_k cre_k
  //   x_i     = P_i - Q_i,  x_{n-i} = P_i + Q_i
  //   P_i = c_0 + sum_k 2 cos(2 pi k i/n) cre_k, Q_i = sum_k 2 sin(...) cim_k
  template <size_t LANES>
  void backward_dft_strip(double* const out, const double* const in,
                          const size_t b0) const {
    const size_t h = (nph_ - 1) / 2;
    const size_t batch = nr_ * nth_;
    std::array<double, 32 * LANES> cre;
    std::array<double, 32 * LANES> cim;
    const double* const c0 = in + b0;
    for (size_t k = 1; k <= h; ++k) {
      const double* const pre = in + (2 * k - 1) * batch + b0;
      const double* const pim = in + (2 * k) * batch + b0;
      double* const re = cre.data() + (k - 1) * LANES;
      double* const im = cim.data() + (k - 1) * LANES;
      for (size_t lane = 0; lane < LANES; ++lane) {
        re[lane] = pre[lane];
        im[lane] = pim[lane];
      }
    }
    {
      std::array<double, LANES> acc{};
      for (size_t k = 0; k < h; ++k) {
        const double* const re = cre.data() + k * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          acc[lane] += re[lane];
        }
      }
      for (size_t lane = 0; lane < LANES; ++lane) {
        out[b0 + lane] = c0[lane] + 2.0 * acc[lane];
      }
    }
    size_t i = 1;
    for (; i + 1 <= h; i += 2) {
      std::array<double, LANES> p0{};
      std::array<double, LANES> q0{};
      std::array<double, LANES> p1{};
      std::array<double, LANES> q1{};
      for (size_t k = 1; k <= h; ++k) {
        const double ca = 2.0 * dft_cos_[(k - 1) * h + (i - 1)];
        const double sa = 2.0 * dft_sin_[(k - 1) * h + (i - 1)];
        const double cb = 2.0 * dft_cos_[(k - 1) * h + i];
        const double sb = 2.0 * dft_sin_[(k - 1) * h + i];
        const double* const re = cre.data() + (k - 1) * LANES;
        const double* const im = cim.data() + (k - 1) * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          p0[lane] += ca * re[lane];
          q0[lane] += sa * im[lane];
          p1[lane] += cb * re[lane];
          q1[lane] += sb * im[lane];
        }
      }
      double* const oi0 = out + i * batch + b0;
      double* const oni0 = out + (nph_ - i) * batch + b0;
      double* const oi1 = out + (i + 1) * batch + b0;
      double* const oni1 = out + (nph_ - i - 1) * batch + b0;
      for (size_t lane = 0; lane < LANES; ++lane) {
        oi0[lane] = c0[lane] + p0[lane] - q0[lane];
        oni0[lane] = c0[lane] + p0[lane] + q0[lane];
        oi1[lane] = c0[lane] + p1[lane] - q1[lane];
        oni1[lane] = c0[lane] + p1[lane] + q1[lane];
      }
    }
    for (; i <= h; ++i) {
      std::array<double, LANES> p{};
      std::array<double, LANES> q{};
      for (size_t k = 1; k <= h; ++k) {
        const double c = 2.0 * dft_cos_[(k - 1) * h + (i - 1)];
        const double s = 2.0 * dft_sin_[(k - 1) * h + (i - 1)];
        const double* const re = cre.data() + (k - 1) * LANES;
        const double* const im = cim.data() + (k - 1) * LANES;
        for (size_t lane = 0; lane < LANES; ++lane) {
          p[lane] += c * re[lane];
          q[lane] += s * im[lane];
        }
      }
      double* const oi = out + i * batch + b0;
      double* const oni = out + (nph_ - i) * batch + b0;
      for (size_t lane = 0; lane < LANES; ++lane) {
        oi[lane] = c0[lane] + p[lane] - q[lane];
        oni[lane] = c0[lane] + p[lane] + q[lane];
      }
    }
  }

  void backward_dft(double* const out, const double* const in) const {
    const size_t batch = nr_ * nth_;
    constexpr size_t W = dft_strip;
    size_t b0 = 0;
    for (; b0 + W <= batch; b0 += W) {
      backward_dft_strip<W>(out, in, b0);
    }
    for (; b0 < batch; ++b0) {
      backward_dft_strip<1>(out, in, b0);
    }
  }

  void build_matrices() {
    // Gauss nodes (colatitudes) and weights from SPHEREPACK's gaqd.
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

    // Fully normalized associated Legendre tables at the Gauss nodes:
    // pbar(j, n) and dpbar/dtheta(j, n) for n = m..l, built per m.
    std::vector<double> pbar(nth_ * (l_ + 1));
    std::vector<double> dpbar(nth_ * (l_ + 1));
    std::vector<double> x(nth_), sx(nth_);
    for (size_t j = 0; j < nth_; ++j) {
      x[j] = std::cos(theta[j]);
      sx[j] = std::sin(theta[j]);
    }

    for (size_t m = 0; m <= l_; ++m) {
      // pmm
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
      // dpbar/dtheta = (n x pbar_n - e_n pbar_{n-1}) / sin(theta),
      // e_n = sqrt((n^2 - m^2)(2n+1)/(2n-1)); pbar_{m-1} term absent at n=m.
      for (size_t n = m; n <= l_; ++n) {
        const double nn = static_cast<double>(n);
        const double mm = static_cast<double>(m);
        const double e_n =
            n == m ? 0.0
                   : std::sqrt((nn * nn - mm * mm) * (2.0 * nn + 1.0) /
                               (2.0 * nn - 1.0));
        for (size_t j = 0; j < nth_; ++j) {
          const double prev = n == m ? 0.0 : pbar[j + (n - 1) * nth_];
          dpbar[j + n * nth_] =
              (nn * x[j] * pbar[j + n * nth_] - e_n * prev) / sx[j];
        }
      }

      // Fused matrices. Scale: FFT round trip contributes n_phi, so 1/n_phi;
      // analysis of the halfcomplex (re, im) columns contributes a further
      // factor 2/n_phi... both absorbed uniformly (validated numerically).
      const double scale = 1.0 / static_cast<double>(nph_);
      double* const dv = dv_mats_.data() + m * nth_ * nth_;
      double* const dw = dw_mats_.data() + m * nth_ * nth_;
      for (size_t i = 0; i < nth_; ++i) {
        for (size_t j = 0; j < nth_; ++j) {
          double acc_v = 0.0;
          double acc_w = 0.0;
          for (size_t n = m; n <= l_; ++n) {
            // SPHEREPACK's vector synthesis (vhsgs1 ndo1/ndo2 bounds) drops
            // the top-degree modes (n = l_max, m even) from the gradient;
            // replicate that truncation to match its output exactly.
            if (n == l_ and m % 2 == 0) {
              continue;
            }
            acc_v += dpbar[i + n * nth_] * pbar[j + n * nth_];
            acc_w += pbar[i + n * nth_] * pbar[j + n * nth_];
          }
          dv[i * nth_ + j] = scale * acc_v * wts[j];
          dw[i * nth_ + j] = scale * static_cast<double>(m) / sx[i] *
                             acc_w * wts[j];
        }
      }
    }
  }

  size_t l_, nth_, nph_, nr_;
  std::vector<double> wsave_, fft_work_, g_;
  std::vector<double> dv_mats_, dw_mats_;
  bool use_dft_kernels_ = false;
  std::vector<double> dft_cos_{}, dft_sin_{};
  std::vector<double> vhat_{}, what_{};
};

DataVector field(const DataVector& x, const DataVector& y, const DataVector& z,
                 const double offset, const double ax, const double ay,
                 const double az) {
  return offset + ax * x + ay * square(y) + az * sin(z + 0.25 * x);
}

struct BenchmarkState {
  explicit BenchmarkState(const Mesh<3>& mesh_in)
      : mesh(mesh_in),
        number_of_points(mesh.number_of_grid_points()),
        u(number_of_points),
        du(number_of_points),
        filtered(number_of_points),
        du_prototype(Variables<vars_tags>::number_of_independent_components *
                         Dim * number_of_points,
                     0.0),
        prototype_buffer(Dim * number_of_points, 0.0) {
    for (auto& logical_deriv : logical_derivs) {
      logical_deriv.initialize(number_of_points);
    }
    const auto logical_coords = logical_coordinates(mesh);
    const auto& x = get<0>(logical_coords);
    const auto& y = get<1>(logical_coords);
    const auto& z = get<2>(logical_coords);

    const size_t number_of_components =
        Variables<vars_tags>::number_of_independent_components;
    for (size_t c = 0; c < number_of_components; ++c) {
      DataVector component_view{u.data() + c * number_of_points,
                                number_of_points};
      component_view = field(x, y, z, 1.0 + 0.03 * static_cast<double>(c),
                             0.01 * static_cast<double>(c + 1),
                             -0.008 * static_cast<double>(c + 2),
                             0.006 * static_cast<double>(c + 3));
    }

    // Smooth, diagonally dominant, non-identity inverse Jacobian so that
    // index mistakes in the prototype cannot hide.
    for (size_t d = 0; d < Dim; ++d) {
      for (size_t i = 0; i < Dim; ++i) {
        inverse_jacobian.get(d, i) =
            (d == i ? 1.0 : 0.0) +
            0.05 * sin(x + 0.3 * y - 0.2 * z +
                       0.4 * static_cast<double>(d + Dim * i)) +
            0.03 * static_cast<double>(d) * 0.1;
        jac_pointers[d * Dim + i] = inverse_jacobian.get(d, i).data();
      }
    }

    const Matrix& diff_matrix =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t n0 = mesh.extents(0);
    diff_matrix_flat.resize(n0 * n0);
    for (size_t i = 0; i < n0; ++i) {
      for (size_t k = 0; k < n0; ++k) {
        diff_matrix_flat[i + n0 * k] = diff_matrix(i, k);
      }
    }

    if (mesh.basis(1) == Spectral::Basis::Legendre) {
      for (size_t d = 0; d < Dim; ++d) {
        gsl::at(filter_matrices, d) = Spectral::filtering::exponential_filter(
            mesh.slice_through(d), 36.0, 32);
      }
    } else if (mesh.basis(1) == Spectral::Basis::SphericalHarmonic) {
      const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
      spectral_coefs.resize(number_of_components * ylm.spectral_size() *
                            mesh.extents(0));
      sph_proto = std::make_unique<SphericalGradientPrototype>(
          mesh.extents(1) - 1, mesh.extents(0));
      sph_proto_du = DataVector(2 * number_of_components * number_of_points,
                                0.0);
    }
  }

  void apply_logical() {
    logical_partial_derivatives(make_not_null(&logical_derivs), u, mesh);
  }

  void apply_full() {
    partial_derivatives(make_not_null(&du), u, mesh, inverse_jacobian);
  }

  void apply_filter() {
    apply_matrices(make_not_null(&filtered), filter_matrices, u,
                   mesh.extents());
  }

  // --- SphericalHarmonic-mesh decomposition pieces (mode 6) ---------------
  // The production path (spherical_apply) is: radial dgemm over all
  // components, then per component gradient_all_offsets = shags analysis +
  // gradgs synthesis. These run each piece in isolation.

  void sh_radial() {
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t deriv_size =
        Variables<vars_tags>::number_of_independent_components *
        number_of_points;
    partial_derivatives_detail::apply_matrix_in_first_dim(
        logical_derivs[0].data(), u.data(), differentiation_matrix_xi,
        deriv_size);
  }

  void sh_analysis() {
    const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
    const size_t n_r = mesh.extents(0);
    const size_t spectral_stride_size = ylm.spectral_size() * n_r;
    for (size_t n = 0;
         n < Variables<vars_tags>::number_of_independent_components; ++n) {
      ylm.phys_to_spec_all_offsets(
          make_not_null(spectral_coefs.data() + n * spectral_stride_size),
          make_not_null(u.data() + n * number_of_points), n_r);
    }
  }

  void sh_synthesis() {
    const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
    const size_t n_r = mesh.extents(0);
    const size_t spectral_stride_size = ylm.spectral_size() * n_r;
    for (size_t n = 0;
         n < Variables<vars_tags>::number_of_independent_components; ++n) {
      const size_t offset = n * number_of_points;
      const auto du_ang = std::array{logical_derivs[1].data() + offset,
                                     logical_derivs[2].data() + offset};
      ylm.gradient_from_coefs_all_offsets(
          du_ang,
          make_not_null(spectral_coefs.data() + n * spectral_stride_size),
          n_r);
    }
  }

  void sh_gradient() {
    const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
    const size_t n_r = mesh.extents(0);
    for (size_t n = 0;
         n < Variables<vars_tags>::number_of_independent_components; ++n) {
      const size_t offset = n * number_of_points;
      const auto du_ang = std::array{logical_derivs[1].data() + offset,
                                     logical_derivs[2].data() + offset};
      ylm.gradient_all_offsets(du_ang, make_not_null(u.data() + offset), n_r);
    }
  }

  // Mode 7: radial dgemm (as in spherical_apply) + prototype fused angular
  // gradient for every component.
  void apply_sph_prototype() {
    sh_radial();
    for (size_t n = 0;
         n < Variables<vars_tags>::number_of_independent_components; ++n) {
      const size_t offset = n * number_of_points;
      sph_proto->apply(sph_proto_du.data() + 2 * offset,
                       sph_proto_du.data() + 2 * offset + number_of_points,
                       u.data() + offset);
    }
  }

  // Max relative difference of the prototype angular derivatives against the
  // SPHEREPACK gradient (normalized by the max magnitude, since the
  // derivative passes through zero).
  // With SPH_DIAG=1, sweep single-(l,m) inputs (same radial profile at all
  // shells) and report the per-mode disagreement, to localize convention
  // errors.
  void sph_mode_sweep() {
    const size_t l_max = mesh.extents(1) - 1;
    const auto& ylm = ylm::get_spherepack_cache(l_max);
    const size_t n_r = mesh.extents(0);
    DataVector coefs(ylm.spectral_size(), 0.0);
    DataVector phys(ylm.physical_size());
    for (size_t l = 0; l <= l_max; ++l) {
      for (size_t m = 0; m <= l; ++m) {
        for (const auto part : {ylm::SpherepackIterator::CoefficientArray::a,
                                ylm::SpherepackIterator::CoefficientArray::b}) {
          if (m == 0 and
              part == ylm::SpherepackIterator::CoefficientArray::b) {
            continue;
          }
          ylm::SpherepackIterator iter(l_max, l_max);
          coefs = 0.0;
          coefs[iter.set(l, m, part)()] = 1.0;
          ylm.spec_to_phys(make_not_null(phys.data()),
                           make_not_null(coefs.data()));
          // replicate at all radial points for component 0
          for (size_t ang = 0; ang < ylm.physical_size(); ++ang) {
            for (size_t r = 0; r < n_r; ++r) {
              u.data()[ang * n_r + r] = phys[ang];
            }
          }
          const auto du_ang =
              std::array{logical_derivs[1].data(), logical_derivs[2].data()};
          ylm.gradient_all_offsets(du_ang, make_not_null(u.data()), n_r);
          sph_proto->apply(sph_proto_du.data(),
                           sph_proto_du.data() + number_of_points, u.data());
          double dth = 0.0;
          double dph = 0.0;
          for (size_t p = 0; p < number_of_points; ++p) {
            dth = std::max(dth,
                           std::abs(logical_derivs[1].data()[p] -
                                    sph_proto_du[p]));
            dph = std::max(dph, std::abs(logical_derivs[2].data()[p] -
                                         sph_proto_du[number_of_points + p]));
          }
          if (dth > 1.0e-10 or dph > 1.0e-10) {
            std::cout << "l=" << l << " m=" << m << " part="
                      << (part == ylm::SpherepackIterator::CoefficientArray::a
                              ? "a"
                              : "b")
                      << "  theta diff " << dth << "  phi diff " << dph
                      << "\n";
          }
        }
      }
    }
  }

  double sph_prototype_max_rel_diff() {
    if (std::getenv("SPH_DIAG") != nullptr) {
      sph_mode_sweep();
    }
    apply_logical();
    apply_sph_prototype();
    double max_mag = 0.0;
    for (size_t c = 0;
         c < Variables<vars_tags>::number_of_independent_components; ++c) {
      for (size_t d = 0; d < 2; ++d) {
        const double* const ref =
            logical_derivs[d + 1].data() + c * number_of_points;
        for (size_t p = 0; p < number_of_points; ++p) {
          max_mag = std::max(max_mag, std::abs(ref[p]));
        }
      }
    }
    std::array<double, 2> max_diff{{0.0, 0.0}};
    for (size_t c = 0;
         c < Variables<vars_tags>::number_of_independent_components; ++c) {
      for (size_t d = 0; d < 2; ++d) {
        const double* const ref =
            logical_derivs[d + 1].data() + c * number_of_points;
        const double* const mine = sph_proto_du.data() +
                                   2 * c * number_of_points +
                                   d * number_of_points;
        for (size_t p = 0; p < number_of_points; ++p) {
          gsl::at(max_diff, d) =
              std::max(gsl::at(max_diff, d), std::abs(ref[p] - mine[p]));
        }
      }
    }
    std::cout << "  theta-deriv diff: " << max_diff[0] / max_mag
              << "  phi-deriv diff: " << max_diff[1] / max_mag << "\n";
    return std::max(max_diff[0], max_diff[1]) / max_mag;
  }

  void apply_prototype() {
    const size_t number_of_components =
        Variables<vars_tags>::number_of_independent_components;
    switch (mesh.extents(0)) {
      case 8:
        prototype_apply<8>(du_prototype.data(), u.data(), number_of_components,
                           diff_matrix_flat.data(), jac_pointers,
                           prototype_buffer.data());
        break;
      case 9:
        prototype_apply<9>(du_prototype.data(), u.data(), number_of_components,
                           diff_matrix_flat.data(), jac_pointers,
                           prototype_buffer.data());
        break;
      case 10:
        prototype_apply<10>(du_prototype.data(), u.data(), number_of_components,
                            diff_matrix_flat.data(), jac_pointers,
                            prototype_buffer.data());
        break;
      case 11:
        prototype_apply<11>(du_prototype.data(), u.data(), number_of_components,
                            diff_matrix_flat.data(), jac_pointers,
                            prototype_buffer.data());
        break;
      case 12:
        prototype_apply<12>(du_prototype.data(), u.data(), number_of_components,
                            diff_matrix_flat.data(), jac_pointers,
                            prototype_buffer.data());
        break;
      default:
        std::cerr << "Prototype only instantiated for 8..12 points per "
                     "dimension, got "
                  << mesh.extents(0) << "\n";
        std::exit(1);
    }
  }

  // Max relative difference between the prototype and partial_derivatives.
  double prototype_max_rel_diff() {
    apply_full();
    apply_prototype();
    double max_rel_diff = 0.0;
    for (size_t j = 0; j < du.size(); ++j) {
      const double rel_diff = std::abs(du.data()[j] - du_prototype[j]) /
                              (std::abs(du.data()[j]) + 1.0e-100);
      max_rel_diff = std::max(max_rel_diff, rel_diff);
    }
    return max_rel_diff;
  }

  double checksum(const size_t mode) const {
    switch (mode) {
      case 0:
      case 4:
        return logical_derivs[0].data()[7] + logical_derivs[1].data()[13] +
               logical_derivs[2].data()[21];
      case 1:
      case 5:
        return du.data()[7] + du.data()[13] + du.data()[21];
      case 2:
        return filtered.data()[7] + filtered.data()[13] + filtered.data()[21];
      case 7:
        return sph_proto_du[7] + sph_proto_du[13] + sph_proto_du[21];
      default:
        return du_prototype[7] + du_prototype[13] + du_prototype[21];
    }
  }

  Mesh<Dim> mesh;
  size_t number_of_points;
  Variables<vars_tags> u;
  InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian{};
  std::array<Variables<vars_tags>, Dim> logical_derivs{};
  Variables<deriv_tags> du;
  Variables<vars_tags> filtered;
  std::array<Matrix, Dim> filter_matrices{};
  std::vector<double> diff_matrix_flat{};
  std::array<const double*, Dim * Dim> jac_pointers{};
  DataVector du_prototype;
  DataVector prototype_buffer;
  std::vector<double> spectral_coefs{};
  std::unique_ptr<SphericalGradientPrototype> sph_proto{};
  DataVector sph_proto_du{};
};
}  // namespace

int main(const int argc, const char* const* const argv) {
  // argv[1] is either a single points-per-dimension or "N0xN1xN2"
  std::array<size_t, 3> extents{10, 10, 10};
  if (argc > 1) {
    const std::string extents_arg{argv[1]};
    if (extents_arg.find('x') != std::string::npos) {
      size_t pos = 0;
      for (size_t d = 0; d < 3; ++d) {
        extents[d] = std::strtoull(extents_arg.c_str() + pos, nullptr, 10);
        pos = extents_arg.find('x', pos) + 1;
      }
    } else {
      extents.fill(std::strtoull(argv[1], nullptr, 10));
    }
  }
  const size_t iterations =
      argc > 2 ? static_cast<size_t>(std::strtoull(argv[2], nullptr, 10))
               : 10000;
  const size_t warmup_iterations =
      argc > 3 ? static_cast<size_t>(std::strtoull(argv[3], nullptr, 10)) : 200;
  const size_t mode =
      argc > 4 ? static_cast<size_t>(std::strtoull(argv[4], nullptr, 10)) : 1;

  const bool spherical = mode >= 4;
  if (spherical and extents[2] != 2 * extents[1] - 1) {
    std::cerr << "Spherical modes need n_phi = 2*n_theta - 1 (l_max = "
                 "n_theta - 1), got "
              << extents[1] << "x" << extents[2] << "\n";
    return 1;
  }
  const Mesh<3> mesh =
      spherical
          ? Mesh<3>{extents,
                    {{Spectral::Basis::Legendre,
                      Spectral::Basis::SphericalHarmonic,
                      Spectral::Basis::SphericalHarmonic}},
                    {{Spectral::Quadrature::GaussLobatto,
                      Spectral::Quadrature::Gauss,
                      Spectral::Quadrature::Equiangular}}}
          : Mesh<3>{extents, Spectral::Basis::Legendre,
                    Spectral::Quadrature::GaussLobatto};

  BenchmarkState state{mesh};

  if (mode == 6) {
    // Time each piece of the spherical derivative separately.
    const auto time_piece = [iterations, warmup_iterations](
                                const char* const name, auto&& piece) {
      for (size_t i = 0; i < warmup_iterations; ++i) {
        piece();
      }
      const auto start = std::chrono::steady_clock::now();
      for (size_t i = 0; i < iterations; ++i) {
        piece();
      }
      const std::chrono::duration<double> elapsed =
          std::chrono::steady_clock::now() - start;
      std::cout << name << ": "
                << 1.0e6 * elapsed.count() / static_cast<double>(iterations)
                << " us\n";
    };
    time_piece("radial_dgemm", [&state]() { state.sh_radial(); });
    time_piece("shags_analysis", [&state]() { state.sh_analysis(); });
    time_piece("gradgs_synthesis", [&state]() { state.sh_synthesis(); });
    time_piece("gradient_combined", [&state]() { state.sh_gradient(); });
    time_piece("logical_all", [&state]() { state.apply_logical(); });
    time_piece("full_partial_derivatives", [&state]() { state.apply_full(); });
    return 0;
  }

  if (mode == 3) {
    std::cout << "prototype_max_rel_diff: " << state.prototype_max_rel_diff()
              << "\n";
  }
  if (mode == 7) {
    std::cout << "sph_prototype_max_rel_diff: "
              << state.sph_prototype_max_rel_diff() << "\n";
  }

  const auto run_once = [&state, mode]() {
    switch (mode) {
      case 0:
      case 4:
        state.apply_logical();
        break;
      case 1:
      case 5:
        state.apply_full();
        break;
      case 2:
        state.apply_filter();
        break;
      case 7:
        state.apply_sph_prototype();
        break;
      default:
        state.apply_prototype();
        break;
    }
  };

  for (size_t i = 0; i < warmup_iterations; ++i) {
    run_once();
  }

  const auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < iterations; ++i) {
    run_once();
  }
  const auto stop = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = stop - start;

  std::cout << "extents: " << extents[0] << "x" << extents[1] << "x"
            << extents[2] << "\n"
            << "number_of_points: " << state.number_of_points << "\n"
            << "iterations: " << iterations << "\n"
            << "mode: " << mode << "\n"
            << "seconds: " << elapsed.count() << "\n"
            << "microseconds_per_call: "
            << 1.0e6 * elapsed.count() / static_cast<double>(iterations) << "\n"
            << "checksum: " << state.checksum(mode) << "\n";
  return 0;
}
