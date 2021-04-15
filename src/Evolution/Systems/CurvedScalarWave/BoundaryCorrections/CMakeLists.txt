# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  CurvedScalarWave
  PRIVATE
  RegisterDerived.cpp
  UpwindPenalty.cpp
  )

spectre_target_headers(
  CurvedScalarWave
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  BoundaryCorrection.hpp
  Factory.hpp
  RegisterDerived.hpp
  UpwindPenalty.hpp
  )
