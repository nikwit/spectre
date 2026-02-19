# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_TimeKokkos")

set(LIBRARY_SOURCES
  Test_KokkosTimeStepperActions.cpp
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DataStructures
  Evolution
  ScalarWave
  ScalarWaveKokkos
  Time
  TimeKokkos
  Utilities
  )
