# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY GeneralizedHarmonicKokkos)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.cpp
  ComputeTimeDerivativeKokkos.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ApplyBoundaryCorrectionsToTimeDerivativeKokkos.hpp
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
  Parallel
  Serialization
  Time
  TimeKokkos
  Utilities
  )
