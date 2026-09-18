# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY NewmanPenrose)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  NullRotations.cpp
  Psi4Fit.cpp
  RestFrame.cpp
  Tetrad.cpp
  TidalResponse.cpp
  TypeD.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  NullRotations.hpp
  Psi4Fit.hpp
  RestFrame.hpp
  Tetrad.hpp
  TidalResponse.hpp
  TypeD.hpp
  Types.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Utilities
  PRIVATE
  ErrorHandling
  LinearSolver
  )
