// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TaggedTuple.hpp"

#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime_api.h>
#endif

namespace ScalarWave::Actions {
namespace detail {

inline void check_cuda_error_and_clear(const char* const context) {
#ifdef KOKKOS_ENABLE_CUDA
  const int expected_device = Kokkos::Cuda().cuda_device();
  int current_device = -1;
  auto err = cudaGetDevice(&current_device);
  ASSERT(err == cudaSuccess,
         "cudaGetDevice failed at " << context << ": "
                                     << cudaGetErrorString(err));
  if (current_device != expected_device) {
    err = cudaSetDevice(expected_device);
    ASSERT(err == cudaSuccess,
           "cudaSetDevice(" << expected_device << ") failed at " << context
                            << ": " << cudaGetErrorString(err));
    err = cudaGetDevice(&current_device);
    ASSERT(err == cudaSuccess && current_device == expected_device,
           "Failed to switch CUDA device at " << context
                                              << ". expected=" << expected_device
                                              << " current=" << current_device);
  }

  Kokkos::fence(context);
  err = cudaGetLastError();
  ASSERT(err == cudaSuccess,
         "CUDA error at " << context << ": " << cudaGetErrorString(err));
#else
  (void)context;
#endif
}

}  // namespace detail

struct CheckCudaOnEvolveEntry {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& /*box*/,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    detail::check_cuda_error_and_clear("ScalarWaveKokkosEvolveEntry");
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct CheckCudaAfterInitialization {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& /*box*/,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    detail::check_cuda_error_and_clear("ScalarWaveKokkosAfterInitialization");
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct CheckCudaOnInitializeTimeStepperHistoryEntry {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& /*box*/,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    detail::check_cuda_error_and_clear(
        "ScalarWaveKokkosInitializeTimeStepperHistoryEntry");
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct CheckCudaOnInitializeTimeStepperHistoryExit {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& /*box*/,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    detail::check_cuda_error_and_clear(
        "ScalarWaveKokkosInitializeTimeStepperHistoryExit");
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace ScalarWave::Actions
