# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarWaveKokkos)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.cpp
  ComputeTimeDerivativeKokkos.cpp
  RecordTimeStepperDataKokkos3D.cpp
  UpdateUKokkos3D.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp
  CleanHistoryKokkos.hpp
  ComputeTimeDerivativeKokkos.hpp
  CudaDiagnostics.hpp
  InitializeKokkosBoundaryCommunication.hpp
  InitializeKokkosTags.hpp
  InitializeKokkosTimeStepperState.hpp
  KokkosBoundaryCommunication.hpp
  KokkosTimeStepperTags.hpp
  RecordTimeStepperDataKokkos.hpp
  RecordTimeStepperDataKokkos3D.hpp
  SendBoundaryCorrectionsKokkos.hpp
  SyncKokkosToHost.hpp
  UpdateUKokkos.hpp
  UpdateUKokkos3D.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Domain
  ErrorHandling
  EventsAndTriggers
  Options
  Parallel
  ScalarWave
  Serialization
  Time
  Utilities
  PRIVATE
  LinearOperators
  )
