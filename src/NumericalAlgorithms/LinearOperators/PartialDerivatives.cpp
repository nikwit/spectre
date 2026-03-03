// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"

#include <array>
#include <blaze/math/DynamicMatrix.h>
#include <cstddef>
#include <functional>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Transpose.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/Spectral/DifferentiationMatrix.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/StdArrayHelpers.hpp"

namespace partial_derivatives_detail {
void apply_matrix_in_first_dim(double* result, const double* const input,
                               const Matrix& matrix, const size_t size,
                               const bool add_to_result) {
  dgemm_<true>(
      'N', 'N',
      matrix.rows(),              // rows of matrix and result
      size / matrix.columns(),    // columns of result and input
      matrix.columns(),           // columns of matrix and rows of input
      1.0,                        // overall multiplier
      matrix.data(),              // matrix
      matrix.spacing(),           // rows of matrix including padding
      input,                      // input
      matrix.columns(),           // rows of input
      add_to_result ? 1.0 : 0.0,  // overwrite output with result or add to it
      result,                     // result
      matrix.rows());             // rows of result
}
void apply_matrix_in_first_dim(std::complex<double>* result,
                               const std::complex<double>* const input,
                               const Matrix& matrix, const size_t size,
                               const bool add_to_result) {
  // BLAS zgemm operates on complex matrices, so we need to copy the real matrix
  // to a complex matrix with zero imaginary part before calling zgemm.
  // Possible performance optimization: avoid the copy here by storing the
  // complex matrix in a static cache. We probably only want to add this to
  // Spectral.hpp once profiling shows that it becomes necessary.
  const blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>
      matrix_complex{matrix};
  zgemm_<true>('N', 'N',
               matrix.rows(),            // rows of matrix and result
               size / matrix.columns(),  // columns of result and input
               matrix.columns(),         // columns of matrix and rows of input
               std::complex{1.0, 0.0},   // overall multiplier
               matrix_complex.data(),    // matrix
               matrix.spacing(),         // rows of matrix including padding
               input,                    // input
               matrix.columns(),         // rows of input
               std::complex{add_to_result ? 1.0 : 0.0,
                            0.0},  // overwrite output with result or add to it
               result,             // result
               matrix.rows());     // rows of result
  // This implementation is ~1.35x slower than the implementation above (based
  // on the "Partial derivatives complex" benchmark in
  // Test_PartialDerivatives.cpp).
  //   DataVector buffer(size * 2);
  //   raw_transpose(make_not_null(reinterpret_cast<double*>(result)),
  //                 reinterpret_cast<const double*>(input), 2, size);
  //   apply_matrix_in_first_dim(buffer.data(),
  //                             reinterpret_cast<const double*>(result),
  //                             matrix, size * 2, add_to_result);
  //   raw_transpose(make_not_null(reinterpret_cast<double*>(result)),
  //                 buffer.data(), size, 2);
}
#ifdef SPECTRE_KOKKOS
template <size_t DerivDim, size_t Dim, bool AddToResult>
void apply_matrix_in_dim(
    Kokkos::View<double**> result,        // [total_num_points, num_components]
    const Kokkos::View<double**>& input,  // [total_num_points, num_components]
    const MatrixViewRO& matrix,           // [num_points_this_dim^2]
    const Mesh<Dim>& mesh,
    const std::array<const double*, Dim * Dim>& inv_jacobian,
    const std::array<size_t, Dim * Dim>& inv_jacobian_strides) {
  const size_t num_components = input.extent(1);
  const std::array<size_t, Dim> extents = mesh.extents().indices();
  const std::array<size_t, Dim - 1> extents_transverse =
      all_but_specified_element_of(extents, DerivDim);
  const size_t total_num_points = mesh.extents().product();
  const size_t num_points_this_dim = extents[DerivDim];
  const size_t num_points_transverse = total_num_points / num_points_this_dim;
  ASSERT(input.extent(0) == total_num_points,
         "Input has " << input.extent(0) << " points, but mesh has "
                      << total_num_points << " points.");
  ASSERT(matrix.extent(0) == matrix.extent(1), "Matrix must be square, but has "
                                                   << matrix.extent(0) << "x"
                                                   << matrix.extent(1));
  ASSERT(matrix.extent(0) == num_points_this_dim,
         "Matrix has " << matrix.extent(0) << " rows, but mesh has "
                       << num_points_this_dim << " points in dim " << DerivDim
                       << '.');

  // Precompute strides for column-major flattening
  std::array<size_t, Dim> strides;
  {
    size_t s = 1;
    for (size_t d = 0; d < Dim; ++d) {
      strides[d] = s;
      s *= extents[d];
    }
  }

  // Precompute storage indices for the inverse Jacobian. Also store as pointers
  // because capturing the full `inv_jacobian` by value and indexing it in
  // nested Kokkos lambdas has a surprisingly large overhead.
  using JacobianStructure = Tensor_detail::Structure<
      Symmetry<1, 2>,
      Tensor_detail::TensorIndexType<Dim, UpLo::Up, Frame::NoFrame,
                                     IndexType::Spatial>,
      Tensor_detail::TensorIndexType<Dim, UpLo::Lo, Frame::NoFrame,
                                     IndexType::Spatial>>;
  std::array<const double*, Dim> inv_jac_components{};
  std::array<size_t, Dim> inv_jac_component_strides{};
  for (size_t d = 0; d < Dim; ++d) {
    const size_t storage_index =
        JacobianStructure::get_storage_index(DerivDim, d);
    gsl::at(inv_jac_components, d) = inv_jacobian[storage_index];
    gsl::at(inv_jac_component_strides, d) = inv_jacobian_strides[storage_index];
  }

  // Parallelization strategy on GPU hardware:
  //
  // We currently sweep over the tensor data Dim times, once for each logical
  // dimension, and assemble the result by contracting with the Jacobian:
  //   \partial_i u = J^\hat{j}_i \partial_\hat{j} u
  // (calling this function Dim times, the first time with AddToResult=false and
  // then later with AddToResult=true to do the sum over $\hat{j}$). In each
  // sweep we parallelize over all independent "stripes" of tensor data in the
  // derivative dimension (typically O(10^3) stripes of length N~6-20 per DG
  // element) by assigning a team of threads (thread block) to each stripe. We
  // load the stripe data into shared team memory, then apply the
  // differentiation matrix to the stripe by parallelizing over the N outputs /
  // rows of the matrix with the thread team (TeamThreadRange) and parallelizing
  // over the N inputs / columns of the matrix with a vectorized sum
  // (ThreadVectorRange). The number of threads in a team (team size) and the
  // number of vector lanes (vector length) are configurable via the TeamPolicy.
  // Currently we just let Kokkos pick the team size automatically and set the
  // vector length to 1 because that gave the best result in a preliminary
  // benchmark with 32 random tensor components in a 3D element with 20^3 grid
  // points (18.5 us on an A100).
  //
  // Notes on performance and possible improvements:
  //
  // - We currently parallelize a single element over GPU threads. It would be
  //   better to batch elements together to saturate the GPU with more work.
  //   This would just extend the `num_components` over the tensor data in all
  //   batched elements to apply the logical differentiation, with the caveat
  //   that the Jacobian data has to be strided over the different elements.
  // - Currently the amount of work per team is very small because it works only
  //   on a single stripe of tensor data with length N~6-20. It could be better
  //   to have each team work on multiple stripes at once, each thread in the
  //   team working on one stripe, and to batch these threads/stripes into
  //   vector lanes so that one warp completes multiple stripes at once (doing
  //   the full matrix multiply in one vector lane). This approach probably
  //   needs to batch elements together or launch kernels for multiple elements
  //   at once to saturate the GPU with enough teams (thread blocks).
  // - The data layout of the input and output data is very important and hasn't
  //   been tested very thoroughly yet, in particular in situations where memory
  //   bandwidth is the bottleneck. It is also important to balance the optimal
  //   data layout for differentiation with the optimal layout for pointwise RHS
  //   evaluations, which is likely different (differentiation prefers locality
  //   of grid points and pointwise evaluations prefer locality of tensor
  //   components). Also, large systems of PDEs like GH might prefer locality of
  //   tensor components whereas small systems like scalar wave might prefer
  //   locality of grid points.

  const size_t num_teams = num_components * num_points_transverse;
  const size_t shmem_size = num_points_this_dim * (1 + Dim) * sizeof(double);
  using TeamPolicy = Kokkos::TeamPolicy<>;
  TeamPolicy team_policy(num_teams, Kokkos::AUTO);
  team_policy = team_policy.set_scratch_size(0, Kokkos::PerTeam(shmem_size));

  Kokkos::parallel_for(
      "apply_matrix_in_dim", team_policy,
      KOKKOS_LAMBDA(const TeamPolicy::member_type& team_member) {
        const size_t team_rank = team_member.league_rank();
        const size_t component_index = team_rank / num_points_transverse;
        size_t remaining = team_rank % num_points_transverse;
        // Decode stripe index in the transverse dims
        std::array<size_t, Dim - 1> stripe_index{};
        for (size_t t = 0; t < Dim - 1; ++t) {
          stripe_index[t] = remaining % extents_transverse[t];
          remaining /= extents_transverse[t];
        }
        // Map stripe to an offset and stride into the column-major data layout
        size_t base_offset = 0;
        {
          size_t t = 0;
          for (size_t d = 0; d < Dim; ++d) {
            if (d == DerivDim) {
              continue;
            }
            base_offset += stripe_index[t++] * strides[d];
          }
        }
        const size_t stride = strides[DerivDim];

        // Preload the input data for this component and stripe
        double* shmem = (double*)team_member.team_shmem().get_shmem(shmem_size);
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> input_stripe(
            shmem, num_points_this_dim);
        Kokkos::View<double* [Dim], Kokkos::MemoryUnmanaged> inv_jac_stripe(
            shmem + num_points_this_dim, num_points_this_dim);

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, num_points_this_dim),
            [=](int j) {
              const size_t point_index =
                  base_offset + static_cast<size_t>(j) * stride;
              input_stripe[j] = input(point_index, component_index);
              for (size_t d = 0; d < Dim; ++d) {
                inv_jac_stripe(j, d) =
                    inv_jac_components[d][point_index *
                                          inv_jac_component_strides[d]];
              }
            });

        team_member.team_barrier();

        // Apply the differentiation matrix and the Jacobian
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, num_points_this_dim),
            [=](int i) {
              // Dense matrix-vector product along DerivDim
              double sum = 0.0;
              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(team_member, num_points_this_dim),
                  [=](int j, double& local_sum) {
                    local_sum += matrix(i, j) * input_stripe[j];
                  },
                  sum);
              // Contract with Jacobian and write out to result
              const size_t point_index =
                  base_offset + static_cast<size_t>(i) * stride;
              for (size_t d = 0; d < Dim; ++d) {
                const double contracted = inv_jac_stripe(i, d) * sum;
                const size_t deriv_component_index = component_index * Dim + d;
                (void)result;  // capture `result` for `if constexpr`
                if constexpr (AddToResult) {
                  result(point_index, deriv_component_index) += contracted;
                } else {
                  result(point_index, deriv_component_index) = contracted;
                }
              }
            });
      });
}

template <size_t DerivDim, size_t Dim, bool AddToResult>
void apply_matrix_in_dim_batched(
    Kokkos::View<double**> result,        // [total_num_points, num_components]
    const Kokkos::View<double**>& input,  // [total_num_points, num_components]
    const MatrixViewRO& matrix,           // [num_points_this_dim^2]
    const Mesh<Dim>& mesh,
    const Kokkos::View<double***>&
        inverse_jacobian) {  // [num_elements, points_per_element, Dim*Dim]
  const size_t num_components = input.extent(1);
  const std::array<size_t, Dim> extents = mesh.extents().indices();
  const std::array<size_t, Dim - 1> extents_transverse =
      all_but_specified_element_of(extents, DerivDim);
  const size_t points_per_element = mesh.extents().product();
  const size_t total_num_points = input.extent(0);
  const size_t num_points_this_dim = extents[DerivDim];
  const size_t num_points_transverse = points_per_element / num_points_this_dim;
  ASSERT(points_per_element > 0, "Points per element must be positive.");
  ASSERT(total_num_points == result.extent(0),
         "Input and result point extents do not match.");
  ASSERT(total_num_points % points_per_element == 0,
         "Input has " << total_num_points
                      << " points, not divisible by points-per-element "
                      << points_per_element << ".");
  const size_t num_elements = total_num_points / points_per_element;
  ASSERT(inverse_jacobian.extent(0) == num_elements,
         "Inverse Jacobian has " << inverse_jacobian.extent(0)
                                 << " elements, expected " << num_elements
                                 << ".");
  ASSERT(inverse_jacobian.extent(1) == points_per_element,
         "Inverse Jacobian has " << inverse_jacobian.extent(1)
                                 << " points per element, expected "
                                 << points_per_element << ".");
  ASSERT(inverse_jacobian.extent(2) == Dim * Dim,
         "Inverse Jacobian has " << inverse_jacobian.extent(2)
                                 << " components, expected " << Dim * Dim
                                 << ".");
  ASSERT(matrix.extent(0) == matrix.extent(1), "Matrix must be square, but has "
                                                   << matrix.extent(0) << "x"
                                                   << matrix.extent(1));
  ASSERT(matrix.extent(0) == num_points_this_dim,
         "Matrix has " << matrix.extent(0) << " rows, but mesh has "
                       << num_points_this_dim << " points in dim " << DerivDim
                       << '.');

  // Precompute strides for column-major flattening
  std::array<size_t, Dim> strides;
  {
    size_t s = 1;
    for (size_t d = 0; d < Dim; ++d) {
      strides[d] = s;
      s *= extents[d];
    }
  }

  const size_t num_teams =
      num_elements * num_components * num_points_transverse;
  const size_t shmem_size = num_points_this_dim * (1 + Dim) * sizeof(double);
  using TeamPolicy = Kokkos::TeamPolicy<>;
  TeamPolicy team_policy(num_teams, Kokkos::AUTO);
  team_policy = team_policy.set_scratch_size(0, Kokkos::PerTeam(shmem_size));

  Kokkos::parallel_for(
      "apply_matrix_in_dim_batched", team_policy,
      KOKKOS_LAMBDA(const TeamPolicy::member_type& team_member) {
        const size_t team_rank = team_member.league_rank();
        const size_t teams_per_element = num_components * num_points_transverse;
        const size_t element_index = team_rank / teams_per_element;
        const size_t element_offset = element_index * points_per_element;
        const size_t element_rank = team_rank % teams_per_element;
        const size_t component_index = element_rank / num_points_transverse;
        size_t remaining = element_rank % num_points_transverse;

        // Decode stripe index in the transverse dims
        std::array<size_t, Dim - 1> stripe_index{};
        for (size_t t = 0; t < Dim - 1; ++t) {
          stripe_index[t] = remaining % extents_transverse[t];
          remaining /= extents_transverse[t];
        }
        // Map stripe to an offset and stride into the element-local data layout
        size_t base_offset = 0;
        {
          size_t t = 0;
          for (size_t d = 0; d < Dim; ++d) {
            if (d == DerivDim) {
              continue;
            }
            base_offset += stripe_index[t++] * strides[d];
          }
        }
        const size_t stride = strides[DerivDim];

        // Preload the input data for this component and stripe
        double* shmem = (double*)team_member.team_shmem().get_shmem(shmem_size);
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> input_stripe(
            shmem, num_points_this_dim);
        Kokkos::View<double* [Dim], Kokkos::MemoryUnmanaged> inv_jac_stripe(
            shmem + num_points_this_dim, num_points_this_dim);

        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, num_points_this_dim),
            [=](int j) {
              const size_t local_point_index =
                  base_offset + static_cast<size_t>(j) * stride;
              const size_t point_index = element_offset + local_point_index;
              input_stripe[j] = input(point_index, component_index);
              for (size_t d = 0; d < Dim; ++d) {
                // Packed inverse Jacobian ordering:
                // [element, point, logical_dim * Dim + inertial_dim]
                inv_jac_stripe(j, d) = inverse_jacobian(
                    element_index, local_point_index, DerivDim * Dim + d);
              }
            });

        team_member.team_barrier();

        // Apply the differentiation matrix and the Jacobian
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team_member, num_points_this_dim),
            [=](int i) {
              // Dense matrix-vector product along DerivDim
              double sum = 0.0;
              Kokkos::parallel_reduce(
                  Kokkos::ThreadVectorRange(team_member, num_points_this_dim),
                  [=](int j, double& local_sum) {
                    local_sum += matrix(i, j) * input_stripe[j];
                  },
                  sum);
              // Contract with Jacobian and write out to result
              const size_t local_point_index =
                  base_offset + static_cast<size_t>(i) * stride;
              const size_t point_index = element_offset + local_point_index;
              for (size_t d = 0; d < Dim; ++d) {
                const double contracted = inv_jac_stripe(i, d) * sum;
                const size_t deriv_component_index = component_index * Dim + d;
                (void)result;  // capture `result` for `if constexpr`
                if constexpr (AddToResult) {
                  result(point_index, deriv_component_index) += contracted;
                } else {
                  result(point_index, deriv_component_index) = contracted;
                }
              }
            });
      });
}

template <size_t Dim>
inline void apply_diff_matrices_fused_batched(
    Kokkos::View<double**> result,        // [total_points, num_components*3]
    const Kokkos::View<double**>& input,  // [total_points, num_components]
    const MatrixViewRO& D0,               // [n0,n0]
    const MatrixViewRO& D1,               // [n1,n1]
    const MatrixViewRO& D2,               // [n2,n2]
    const Mesh<Dim>& mesh,
    const Kokkos::View<double***>& inverse_jacobian) {  // [nelems, ppe, 9]
  static_assert(Dim == 3, "Assumes Dim==3.");

  const auto ext = mesh.extents().indices();
  const int n0 = static_cast<int>(ext[0]);
  const int n1 = static_cast<int>(ext[1]);
  const int n2 = static_cast<int>(ext[2]);
  const int ppe = static_cast<int>(mesh.extents().product());

  const int total_points = static_cast<int>(input.extent(0));
  const int num_components = static_cast<int>(input.extent(1));
  if (ppe == 0 || total_points == 0 || num_components == 0)
    return;

  ASSERT(static_cast<size_t>(total_points) % static_cast<size_t>(ppe) == 0,
         "total_points must be divisible by points_per_element");
  const int nelems = total_points / ppe;

  using exec_space = typename Kokkos::View<double**>::execution_space;
  using team_policy = Kokkos::TeamPolicy<exec_space>;
  using member_type = team_policy::member_type;

  // Tuneable: 4 works well on A100; 5 may also be good for 35 comps.
  // Keep small to control registers.
  constexpr int tileC = 2;

  // Shared holds u(lp, tileC) for one element. Layout: [tileC][ppe] contiguous
  // in lp.
  const size_t sh_bytes =
      static_cast<size_t>(tileC) * static_cast<size_t>(ppe) * sizeof(double);

  const int ntiles = (num_components + tileC - 1) / tileC;

  // league over (element, tile)
  team_policy pol(nelems * ntiles, 256);
  Kokkos::parallel_for(
      "apply_diff_matrices_fused_tiled_components_batched",
      pol.set_scratch_size(0, Kokkos::PerTeam(sh_bytes)),
      KOKKOS_LAMBDA(const member_type& team) {
        const int league = team.league_rank();
        const int e = league / ntiles;
        const int tile = league - e * ntiles;
        const int c0 = tile * tileC;
        const int cN =
            (c0 + tileC <= num_components) ? tileC : (num_components - c0);

        const size_t elem_base =
            static_cast<size_t>(e) * static_cast<size_t>(ppe);

        double* sh = (double*)team.team_shmem().get_shmem(sh_bytes);

        // load only this tile's components into shared
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ppe), [&](int lp) {
          const size_t pidx = elem_base + static_cast<size_t>(lp);
#pragma unroll
          for (int t = 0; t < tileC; ++t) {
            if (t < cN)
              sh[t * ppe + lp] = input(pidx, c0 + t);
          }
        });
        team.team_barrier();

        // Compute per point; each thread does tileC components worth of work
        // (ILP).
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ppe), [&](int lp) {
          int tmp = lp;
          const int i0 = tmp % n0;
          tmp /= n0;
          const int i1 = tmp % n1;
          const int i2 = tmp / n1;

          // Prefetch invJ row block (LayoutLeft: e,lp contiguous; idx is
          // strided but only 9)
          double j[9];
#pragma unroll
          for (int k = 0; k < 9; ++k) {
            j[k] = inverse_jacobian(e, lp, k);
          }

          double sum0[tileC], sum1[tileC], sum2[tileC];
#pragma unroll
          for (int t = 0; t < tileC; ++t) {
            sum0[t] = 0.0;
            sum1[t] = 0.0;
            sum2[t] = 0.0;
          }

          // Dim 0
          {
            const int base = (i1 + n1 * i2) * n0;  // in lp indexing
#pragma unroll
            for (int k0 = 0; k0 < 16; ++k0) {
              if (k0 >= n0)
                break;
              const int src_lp = base + k0;
              const double a = D0(i0, k0);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum0[t] += a * sh[t * ppe + src_lp];
              }
            }
            for (int k0 = 16; k0 < n0; ++k0) {
              const int src_lp = base + k0;
              const double a = D0(i0, k0);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum0[t] += a * sh[t * ppe + src_lp];
              }
            }
          }

          // Dim 1
          {
#pragma unroll
            for (int k1 = 0; k1 < 16; ++k1) {
              if (k1 >= n1)
                break;
              const int src_lp = i0 + n0 * (k1 + n1 * i2);
              const double a = D1(i1, k1);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum1[t] += a * sh[t * ppe + src_lp];
              }
            }
            for (int k1 = 16; k1 < n1; ++k1) {
              const int src_lp = i0 + n0 * (k1 + n1 * i2);
              const double a = D1(i1, k1);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum1[t] += a * sh[t * ppe + src_lp];
              }
            }
          }

          // Dim 2
          {
#pragma unroll
            for (int k2 = 0; k2 < 16; ++k2) {
              if (k2 >= n2)
                break;
              const int src_lp = i0 + n0 * (i1 + n1 * k2);
              const double a = D2(i2, k2);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum2[t] += a * sh[t * ppe + src_lp];
              }
            }
            for (int k2 = 16; k2 < n2; ++k2) {
              const int src_lp = i0 + n0 * (i1 + n1 * k2);
              const double a = D2(i2, k2);
#pragma unroll
              for (int t = 0; t < tileC; ++t) {
                if (t < cN)
                  sum2[t] += a * sh[t * ppe + src_lp];
              }
            }
          }

          const size_t pidx = elem_base + static_cast<size_t>(lp);

#pragma unroll
          for (int t = 0; t < tileC; ++t) {
            if (t >= cN)
              break;

            const double out0 =
                j[0] * sum0[t] + j[3] * sum1[t] + j[6] * sum2[t];
            const double out1 =
                j[1] * sum0[t] + j[4] * sum1[t] + j[7] * sum2[t];
            const double out2 =
                j[2] * sum0[t] + j[5] * sum1[t] + j[8] * sum2[t];

            const int col0 = (c0 + t) * 3 + 0;
            const int col1 = (c0 + t) * 3 + 1;
            const int col2 = (c0 + t) * 3 + 2;

            result(pidx, col0) = out0;
            result(pidx, col1) = out1;
            result(pidx, col2) = out2;
          }
        });

        team.team_barrier();
      });
}

// generate instantations
#define INSTANTIATE_APPLY_MATRIX_IN_DIM(DerivDim, Dim, AddToResult)       \
  template void apply_matrix_in_dim<DerivDim, Dim, AddToResult>(          \
      Kokkos::View<double**> result, const Kokkos::View<double**>& input, \
      const MatrixViewRO& matrix, const Mesh<Dim>& mesh,                  \
      const std::array<const double*, Dim * Dim>& inv_jacobian,           \
      const std::array<size_t, Dim * Dim>& inv_jacobian_strides);
INSTANTIATE_APPLY_MATRIX_IN_DIM(0, 1, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM(0, 2, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM(0, 3, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM(1, 2, true)
INSTANTIATE_APPLY_MATRIX_IN_DIM(1, 3, true)
INSTANTIATE_APPLY_MATRIX_IN_DIM(2, 3, true)
#undef INSTANTIATE_APPLY_MATRIX_IN_DIM
#define INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(DerivDim, Dim, AddToResult) \
  template void apply_matrix_in_dim_batched<DerivDim, Dim, AddToResult>(    \
      Kokkos::View<double**> result, const Kokkos::View<double**>& input,   \
      const MatrixViewRO& matrix, const Mesh<Dim>& mesh,                    \
      const Kokkos::View<double***>& inverse_jacobian);
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(0, 1, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(0, 2, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(0, 3, false)
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(1, 2, true)
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(1, 3, true)
INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED(2, 3, true)
#undef INSTANTIATE_APPLY_MATRIX_IN_DIM_BATCHED
#define INSTANTIATE_APPLY_DIFF_MATRICES_FUSED_BATCHED(Dim)                \
  template void                                                           \
  partial_derivatives_detail::apply_diff_matrices_fused_batched<Dim>(     \
      Kokkos::View<double**> result, const Kokkos::View<double**>& input, \
      const MatrixViewRO& matrix_dim_0, const MatrixViewRO& matrix_dim_1, \
      const MatrixViewRO& matrix_dim_2, const Mesh<Dim>& mesh,            \
      const Kokkos::View<double***>& inverse_jacobian);

INSTANTIATE_APPLY_DIFF_MATRICES_FUSED_BATCHED(3)

#undef INSTANTIATE_APPLY_DIFF_MATRICES_FUSED_BATCHED
#endif  // SPECTRE_KOKKOS
}  // namespace partial_derivatives_detail

template <typename DataType, typename SymmList, typename IndexList, size_t Dim>
void logical_partial_derivative(
    const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>*>
        logical_derivative_of_u,
    const gsl::not_null<gsl::span<typename DataType::value_type>*> buffer,
    const Tensor<DataType, SymmList, IndexList>& u, const Mesh<Dim>& mesh) {
  static_assert(
      Dim > 0 and Dim < 4,
      "logical_partial_derivative is only implemented for 1, 2, and 3d");
  const size_t num_grid_points = mesh.number_of_grid_points();
  ASSERT(buffer->size() >= num_grid_points,
         "The buffer in logical_partial_derivative must be at least of size "
             << num_grid_points << " but is of size " << buffer->size());

  set_number_of_grid_points(logical_derivative_of_u,
                            mesh.number_of_grid_points());
  if (Dim == 3 and mesh.basis(1) == Spectral::Basis::SphericalHarmonic) {
    if constexpr (std::is_same_v<typename DataType::value_type, double>) {
      const Matrix& differentiation_matrix_xi =
          Spectral::differentiation_matrix(mesh.slice_through(0));
      const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
      for (size_t storage_index = 0; storage_index < u.size();
           ++storage_index) {
        const auto u_tensor_index = u.get_tensor_index(storage_index);
        partial_derivatives_detail::apply_matrix_in_first_dim(
            // NOLINTNEXTLINE(readability-redundant-smartptr-get)
            logical_derivative_of_u->get(prepend(u_tensor_index, 0_st)).data(),
            u[storage_index].data(), differentiation_matrix_xi,
            num_grid_points);
        const auto du = std::array{
            // NOLINTNEXTLINE(readability-redundant-smartptr-get)
            logical_derivative_of_u->get(prepend(u_tensor_index, 1_st)).data(),
            // NOLINTNEXTLINE(readability-redundant-smartptr-get)
            logical_derivative_of_u->get(prepend(u_tensor_index, 2_st)).data()};
        ylm.gradient_all_offsets(du, make_not_null(u[storage_index].data()),
                                 mesh.extents(0));
      }
    } else {
      ERROR(
          "Support for complex numbers with spherical harmonics is not yet "
          "implemented for logical_partial_derivative.");
    }
  } else {
    const Matrix empty_matrix{};
    std::array<std::reference_wrapper<const Matrix>, Dim> diff_matrices{
        make_array<Dim, std::reference_wrapper<const Matrix>>(empty_matrix)};
    for (size_t d = 0; d < Dim; ++d) {
      gsl::at(diff_matrices, d) =
          std::cref(Spectral::differentiation_matrix(mesh.slice_through(d)));
    }

    // It would be possible to check if the memory is contiguous and then
    // differentiate all components at once. Note that the buffer in that case
    // would also need to be the size of all components.
    for (size_t storage_index = 0; storage_index < u.size(); ++storage_index) {
      const auto u_tensor_index = u.get_tensor_index(storage_index);
      const auto xi_deriv_tensor_index = prepend(u_tensor_index, 0_st);
      partial_derivatives_detail::apply_matrix_in_first_dim(
          // NOLINTNEXTLINE(readability-redundant-smartptr-get)
          logical_derivative_of_u->get(xi_deriv_tensor_index).data(),
          u[storage_index].data(), diff_matrices[0].get(), num_grid_points);
      for (size_t i = 1; i < Dim; ++i) {
        const auto deriv_tensor_index = prepend(u_tensor_index, i);
        DataType& deriv_component =
            logical_derivative_of_u->get(deriv_tensor_index);
        size_t chunk_size =
            diff_matrices[0].get().rows() *
            (i == 1 ? 1 : gsl::at(diff_matrices, 1).get().rows());
        raw_transpose(make_not_null(deriv_component.data()),
                      u[storage_index].data(), chunk_size,
                      num_grid_points / chunk_size);
        partial_derivatives_detail::apply_matrix_in_first_dim(
            buffer->data(), deriv_component.data(),
            gsl::at(diff_matrices, i).get(), num_grid_points);
        chunk_size =
            i == 1 ? (Dim == 2 ? gsl::at(diff_matrices, 1).get().rows()
                               : gsl::at(diff_matrices, 1).get().rows() *
                                     gsl::at(diff_matrices, 2).get().rows())
                   : gsl::at(diff_matrices, 2).get().rows();
        raw_transpose(make_not_null(deriv_component.data()), buffer->data(),
                      chunk_size, num_grid_points / chunk_size);
      }
    }
  }
}

template <typename DataType, typename SymmList, typename IndexList, size_t Dim>
void logical_partial_derivative(
    gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>*>
        logical_derivative_of_u,
    const Tensor<DataType, SymmList, IndexList>& u, const Mesh<Dim>& mesh) {
  using ValueType = typename DataType::value_type;  // double or complex<double>
  std::vector<ValueType> buffer(mesh.number_of_grid_points());
  gsl::span<ValueType> buffer_view{buffer.data(), buffer.size()};
  logical_partial_derivative(logical_derivative_of_u,
                             make_not_null(&buffer_view), u, mesh);
}

template <typename DataType, typename SymmList, typename IndexList, size_t Dim>
auto logical_partial_derivative(const Tensor<DataType, SymmList, IndexList>& u,
                                const Mesh<Dim>& mesh)
    -> TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical> {
  TensorMetafunctions::prepend_spatial_index<
      Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
      Frame::ElementLogical>
      result{mesh.number_of_grid_points()};
  logical_partial_derivative(make_not_null(&result), u, mesh);
  return result;
}

template <typename DataType, typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
void partial_derivative(
    const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo, DerivativeFrame>*>
        du,
    const TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
        Frame::ElementLogical>& logical_partial_derivative_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  for (size_t storage_index = 0;
       storage_index < Tensor<DataType, SymmList, IndexList>::size();
       ++storage_index) {
    const auto u_multi_index =
        Tensor<DataType, SymmList,
               IndexList>::structure::get_canonical_tensor_index(storage_index);
    for (size_t i = 0; i < Dim; i++) {
      const auto du_multi_index = prepend(u_multi_index, i);
      du->get(du_multi_index) =
          inverse_jacobian.get(0, i) *
          logical_partial_derivative_of_u.get(prepend(u_multi_index, 0_st));
      for (size_t j = 1; j < Dim; j++) {
        du->get(du_multi_index) +=
            inverse_jacobian.get(j, i) *
            logical_partial_derivative_of_u.get(prepend(u_multi_index, j));
      }
    }
  }
}

template <typename DataType, typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
void partial_derivative(
    const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo, DerivativeFrame>*>
        du,
    const Tensor<DataType, SymmList, IndexList>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  TensorMetafunctions::prepend_spatial_index<
      Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo,
      Frame::ElementLogical>
      logical_partial_derivative_of_u{mesh.number_of_grid_points()};
  logical_partial_derivative(make_not_null(&logical_partial_derivative_of_u), u,
                             mesh);
  partial_derivative<DataType, SymmList, IndexList>(
      du, logical_partial_derivative_of_u, inverse_jacobian);
}

template <typename DataType, typename SymmList, typename IndexList, size_t Dim,
          typename DerivativeFrame>
auto partial_derivative(
    const Tensor<DataType, SymmList, IndexList>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian)
    -> TensorMetafunctions::prepend_spatial_index<
        Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo, DerivativeFrame> {
  TensorMetafunctions::prepend_spatial_index<
      Tensor<DataType, SymmList, IndexList>, Dim, UpLo::Lo, DerivativeFrame>
      result{mesh.number_of_grid_points()};
  partial_derivative(make_not_null(&result), u, mesh, inverse_jacobian);
  return result;
}

#define GET_DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(1, data)
#define GET_FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)
#define GET_TENSOR(data) BOOST_PP_TUPLE_ELEM(3, data)

#define INSTANTIATION(r, data)                                                 \
  template void logical_partial_derivative(                                    \
      gsl::not_null<TensorMetafunctions::prepend_spatial_index<                \
                        GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),     \
                        GET_FRAME(data)>,                                      \
                    GET_DIM(data), UpLo::Lo, Frame::ElementLogical>* >         \
          logical_derivative_of_u,                                             \
      gsl::not_null<gsl::span<typename GET_DTYPE(data)::value_type>*> buffer,  \
      const GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),                 \
      GET_FRAME(data) > &u, const Mesh<GET_DIM(data)>& mesh);                  \
  template void logical_partial_derivative(                                    \
      gsl::not_null<TensorMetafunctions::prepend_spatial_index<                \
                        GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),     \
                        GET_FRAME(data)>,                                      \
                    GET_DIM(data), UpLo::Lo, Frame::ElementLogical>* >         \
          logical_derivative_of_u,                                             \
      const GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),                 \
      GET_FRAME(data) > &u, const Mesh<GET_DIM(data)>& mesh);                  \
  template TensorMetafunctions::prepend_spatial_index<                         \
      GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>,     \
      GET_DIM(data), UpLo::Lo,                                                 \
      Frame::ElementLogical >                                                  \
          logical_partial_derivative(const GET_TENSOR(data) < GET_DTYPE(data), \
                                     GET_DIM(data), GET_FRAME(data) > &u,      \
                                     const Mesh<GET_DIM(data)>& mesh);         \
  template void                                                                \
      partial_derivative<GET_DTYPE(data), GET_TENSOR(data) < GET_DTYPE(data),  \
                         GET_DIM(data), GET_FRAME(data)>::symmetry,            \
      GET_TENSOR(                                                              \
          data)<GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>::index_list > \
          (const gsl::not_null<TensorMetafunctions::prepend_spatial_index<     \
                                   GET_TENSOR(data) < GET_DTYPE(data),         \
                                   GET_DIM(data), GET_FRAME(data)>,            \
                               GET_DIM(data), UpLo::Lo, GET_FRAME(data)>* >    \
               du,                                                             \
           const TensorMetafunctions::prepend_spatial_index<                   \
               GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),              \
               GET_FRAME(data)>,                                               \
           GET_DIM(data), UpLo::Lo,                                            \
           Frame::ElementLogical > &logical_partial_derivative_of_u,           \
           const InverseJacobian<DataVector, GET_DIM(data),                    \
                                 Frame::ElementLogical, GET_FRAME(data)>       \
               & inverse_jacobian);                                            \
  template void partial_derivative(                                            \
      const gsl::not_null<TensorMetafunctions::prepend_spatial_index<          \
                              GET_TENSOR(data) < GET_DTYPE(data),              \
                              GET_DIM(data), GET_FRAME(data)>,                 \
                          GET_DIM(data), UpLo::Lo, GET_FRAME(data)>* > du,     \
      const GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),                 \
      GET_FRAME(data) > &u, const Mesh<GET_DIM(data)>& mesh,                   \
      const InverseJacobian<DataVector, GET_DIM(data), Frame::ElementLogical,  \
                            GET_FRAME(data)>& inverse_jacobian);               \
  template TensorMetafunctions::prepend_spatial_index<                         \
      GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>,     \
      GET_DIM(data), UpLo::Lo,                                                 \
      GET_FRAME(data) >                                                        \
          partial_derivative(                                                  \
              const GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),         \
              GET_FRAME(data) > &u, const Mesh<GET_DIM(data)>& mesh,           \
              const InverseJacobian<DataVector, GET_DIM(data),                 \
                                    Frame::ElementLogical, GET_FRAME(data)>&   \
                  inverse_jacobian);

GENERATE_INSTANTIATIONS(INSTANTIATION, (DataVector, ComplexDataVector),
                        (1, 2, 3),
                        (Frame::Grid, Frame::Distorted, Frame::Inertial),
                        (tnsr::a, tnsr::A, tnsr::i, tnsr::I, tnsr::ab, tnsr::Ab,
                         tnsr::aB, tnsr::AB, tnsr::ij, tnsr::iJ, tnsr::Ij,
                         tnsr::IJ, tnsr::iA, tnsr::ia, tnsr::aa, tnsr::AA,
                         tnsr::ii, tnsr::II, tnsr::ijj, tnsr::Ijj, tnsr::iaa))

#undef INSTANTIATION

// Some additional mixed-dimension instantiations
template TensorMetafunctions::prepend_spatial_index<
    tnsr::aa<ComplexDataVector, 3, Frame::Inertial>, 2, UpLo::Lo,
    Frame::Inertial>
partial_derivative(const tnsr::aa<ComplexDataVector, 3, Frame::Inertial>& u,
                   const Mesh<2>& mesh,
                   const InverseJacobian<DataVector, 2, Frame::ElementLogical,
                                         Frame::Inertial>& inverse_jacobian);

#define INSTANTIATION(r, data)                                                 \
  template void logical_partial_derivative(                                    \
      gsl::not_null<TensorMetafunctions::prepend_spatial_index<                \
          Scalar<GET_DTYPE(data)>, GET_DIM(data), UpLo::Lo,                    \
          Frame::ElementLogical>*>                                             \
          logical_derivative_of_u,                                             \
      gsl::not_null<gsl::span<typename GET_DTYPE(data)::value_type>*> buffer,  \
      const Scalar<GET_DTYPE(data)>& u, const Mesh<GET_DIM(data)>& mesh);      \
  template void logical_partial_derivative(                                    \
      gsl::not_null<TensorMetafunctions::prepend_spatial_index<                \
          Scalar<GET_DTYPE(data)>, GET_DIM(data), UpLo::Lo,                    \
          Frame::ElementLogical>*>                                             \
          logical_derivative_of_u,                                             \
      const Scalar<GET_DTYPE(data)>& u, const Mesh<GET_DIM(data)>& mesh);      \
  template TensorMetafunctions::prepend_spatial_index<                         \
      Scalar<GET_DTYPE(data)>, GET_DIM(data), UpLo::Lo, Frame::ElementLogical> \
  logical_partial_derivative(const Scalar<GET_DTYPE(data)>& u,                 \
                             const Mesh<GET_DIM(data)>& mesh);

GENERATE_INSTANTIATIONS(INSTANTIATION, (DataVector, ComplexDataVector),
                        (1, 2, 3))

#undef INSTANTIATION

#define INSTANTIATE_JACOBIANS(r, data)                                       \
  template TensorMetafunctions::prepend_spatial_index<                       \
      GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),                     \
      Frame::ElementLogical, GET_FRAME(data)>,                               \
      GET_DIM(data), UpLo::Lo,                                               \
      GET_FRAME(data) >                                                      \
          partial_derivative(                                                \
              const GET_TENSOR(data) < GET_DTYPE(data), GET_DIM(data),       \
              Frame::ElementLogical, GET_FRAME(data) > &u,                   \
              const Mesh<GET_DIM(data)>& mesh,                               \
              const InverseJacobian<DataVector, GET_DIM(data),               \
                                    Frame::ElementLogical, GET_FRAME(data)>& \
                  inverse_jacobian);

GENERATE_INSTANTIATIONS(INSTANTIATE_JACOBIANS, (DataVector), (1, 2, 3),
                        (Frame::Inertial), (InverseJacobian))

#undef INSTANTIATE_JACOBIANS

#define INSTANTIATE_SCALAR(r, data)                                            \
  template void partial_derivative<GET_DTYPE(data), Symmetry<>, index_list<>>( \
      gsl::not_null<tnsr::i<GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>*> \
          du,                                                                  \
      const tnsr::i<GET_DTYPE(data), GET_DIM(data), Frame::ElementLogical>&    \
          logical_partial_derivative_of_u,                                     \
      const InverseJacobian<DataVector, GET_DIM(data), Frame::ElementLogical,  \
                            GET_FRAME(data)>& inverse_jacobian);               \
  template void partial_derivative(                                            \
      gsl::not_null<tnsr::i<GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>*> \
          du,                                                                  \
      const Scalar<GET_DTYPE(data)>& u, const Mesh<GET_DIM(data)>& mesh,       \
      const InverseJacobian<DataVector, GET_DIM(data), Frame::ElementLogical,  \
                            GET_FRAME(data)>& inverse_jacobian);               \
  template tnsr::i<GET_DTYPE(data), GET_DIM(data), GET_FRAME(data)>            \
  partial_derivative(                                                          \
      const Scalar<GET_DTYPE(data)>& u, const Mesh<GET_DIM(data)>& mesh,       \
      const InverseJacobian<DataVector, GET_DIM(data), Frame::ElementLogical,  \
                            GET_FRAME(data)>& inverse_jacobian);

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALAR, (DataVector, ComplexDataVector),
                        (1, 2, 3),
                        (Frame::Grid, Frame::Distorted, Frame::Inertial))

#undef INSTANTIATE_SCALAR
#undef GET_FRAME
#undef GET_DIM
#undef GET_TENSOR
