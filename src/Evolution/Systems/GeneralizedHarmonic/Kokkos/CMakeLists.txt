# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY GeneralizedHarmonicKokkos)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.cpp
  Batched/ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched.cpp
  Batched/ComputeInternalBoundaryTermsBatched.cpp
  Batched/ComputeTimeDerivativeBatched.cpp
  Batched/FilterKokkosBatched.cpp
  Batched/HostElementView.cpp
  Batched/InitializeBatchedData.cpp
  Batched/LiftInternalBoundaryTermsBatched.cpp
  Batched/UpdateConstraintDampingBatched.cpp
  Batched/UpdateFrameFieldsBatched.cpp
  Batched/ObserveFieldsBatchedEvent.cpp
  Batched/ObserveTimeStepBatchedEvent.cpp
  Batched/PackageLocalFacesBatched.cpp
  ComputeTimeDerivativeKokkos.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp
  Batched/Actions.hpp
  Batched/BoundaryCorrectionHelpers.hpp
  Batched/ApplyExternalBoundaryCorrectionsToTimeDerivativeBatched.hpp
  Batched/ComputeInternalBoundaryTermsBatched.hpp
  Batched/ComputeTimeDerivativeBatched.hpp
  Batched/DriverState.hpp
  Batched/FilterKokkosBatched.hpp
  Batched/FindApparentHorizonBatchedEvent.hpp
  Batched/HostElementView.hpp
  Batched/InitializeBatchedData.hpp
  Batched/InitializeBoundaryBatchMetadata.hpp
  Batched/LiftInternalBoundaryTermsBatched.hpp
  Batched/ObserveFieldsBatchedEvent.hpp
  Batched/ObserveTimeStepBatchedEvent.hpp
  Batched/PackageLocalFacesBatched.hpp
  Batched/UpdateConstraintDampingBatched.hpp
  Batched/UpdateFrameFieldsBatched.hpp
  CleanHistoryKokkos.hpp
  ComputeTimeDerivativeKokkos.hpp
  FilterKokkos.hpp
  InitializeKokkosBoundaryCommunication.hpp
  InitializeKokkosTags.hpp
  InitializeKokkosTimeStepperState.hpp
  KokkosBoundaryCommunication.hpp
  KokkosTimeStepperTags.hpp
  RecordTimeStepperDataKokkos.hpp
  SyncKokkosToHost.hpp
  UpdateUKokkos.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Domain
  ErrorHandling
  GeneralizedHarmonic
  LinearOperators
  Parallel
  Serialization
  Time
  TimeKokkos
  Utilities
  )
