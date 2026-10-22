# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarWaveWorldtube)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  AccelerationTermsOrder0.cpp
  AccelerationTermsOrder1.cpp
  AccelerationTermsOrder2.cpp
  KerrSchildDerivatives.cpp
  Tags.cpp
  PunctureField.cpp
  PunctureFieldOrder0.cpp
  PunctureFieldOrder1.cpp
  PunctureFieldOrder2.cpp
  RadiusFunctions.cpp
  SelfForce.cpp
  SelfForceDerivatives.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Inboxes.hpp
  KerrSchildDerivatives.hpp
  SingletonChare.hpp
  Tags.hpp
  PunctureField.hpp
  PunctureFieldWorkspace.hpp
  RadiusFunctions.hpp
  SelfForce.hpp
  SelfForceDerivatives.hpp
  System.hpp
  Worldtube.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  Autodiff
  PUBLIC
  DataStructures
  Domain
  GeneralRelativity
  GeneralRelativitySolutions
  LinearOperators
  Options
  Parallel
  Time
  Utilities
  )

add_subdirectory(SingletonActions)
add_subdirectory(ElementActions)
add_subdirectory(Instantiations)
add_subdirectory(Triggers)
add_subdirectory(InitialData)
