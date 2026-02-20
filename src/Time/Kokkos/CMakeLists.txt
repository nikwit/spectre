# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY TimeKokkos)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  DebugChecks.cpp
  KokkosModule.cpp
  RecordTimeStepperData.cpp
  UpdateU.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  CleanHistory.hpp
  DebugChecks.hpp
  RecordTimeStepperData.hpp
  UpdateU.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  ErrorHandling
  Time
  Utilities
  )
