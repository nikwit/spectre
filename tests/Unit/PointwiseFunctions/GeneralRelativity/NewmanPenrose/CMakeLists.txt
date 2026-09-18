# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_NewmanPenrose")

set(LIBRARY_SOURCES
  Test_NullRotations.cpp
  Test_Psi4Fit.cpp
  Test_RestFrame.cpp
  Test_Tetrad.cpp
  Test_TidalResponse.cpp
  Test_TypeD.cpp
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DataStructures
  DataStructuresHelpers
  NewmanPenrose
  Utilities
  )
