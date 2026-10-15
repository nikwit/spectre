# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY NewmanPenrose)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  CoulombDecode.cpp
  GeometricTide.cpp
  GeometricTideFixedPoint.cpp
  GeometricTideTransport.cpp
  LongitudinalTide.cpp
  NullRotations.cpp
  Psi4Fit.cpp
  RestFrame.cpp
  Tetrad.cpp
  ThirdOrderTide.cpp
  TidalResponse.cpp
  TypeD.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  CoulombDecode.hpp
  GeometricTide.hpp
  GeometricTideFixedPoint.hpp
  GeometricTideTransport.hpp
  LongitudinalTide.hpp
  NullRotations.hpp
  Psi4Fit.hpp
  RestFrame.hpp
  Tetrad.hpp
  ThirdOrderTide.hpp
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
  GSL::gsl
  ErrorHandling
  LinearSolver
  )
