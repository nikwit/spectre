# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ScalarWaveWorldtube)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  AccelerationTermsOrder0.cpp
  AccelerationTermsOrder1.cpp
  KerrSchildDerivatives.cpp
  Tags.cpp
  PunctureField.cpp
  PunctureFieldOrder0.cpp
  PunctureFieldOrder1.cpp
  PunctureFieldOrder2.cpp
  PunctureFieldOrder2Part0.cpp
  PunctureFieldOrder2Part1.cpp
  PunctureFieldOrder2Part2.cpp
  PunctureFieldOrder2Part3.cpp
  PunctureFieldOrder2Part4.cpp
  PunctureFieldOrder2Part5.cpp
  PunctureFieldOrder2Part6.cpp
  PunctureFieldOrder2Part7.cpp
  PunctureFieldOrder2Part8.cpp
  PunctureFieldOrder2Part9.cpp
  PunctureFieldOrder2Part10.cpp
  PunctureFieldOrder2Part11.cpp
  RadiusFunctions.cpp
  SelfForce.cpp
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
  PunctureFieldOrder2Impl.hpp
  RadiusFunctions.hpp
  SelfForce.hpp
  System.hpp
  Worldtube.hpp
  )

target_link_libraries(
  ${LIBRARY}
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
