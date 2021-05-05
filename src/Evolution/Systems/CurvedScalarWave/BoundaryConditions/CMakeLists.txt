# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  BoundaryCondition.cpp
  ConstraintPreserving.cpp
  ConstraintPreservingBaylissTurkel.cpp
  #DirichletAnalytic.cpp
  Freezing.cpp
  RegisterDerivedWithCharm.cpp
  #SphericalRadiation.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  BoundaryCondition.hpp
  ConstraintPreserving.hpp
  ConstraintPreservingBaylissTurkel.hpp
  #DirichletAnalytic.hpp
  Factory.hpp
  Freezing.hpp
  RegisterDerivedWithCharm.hpp
  #SphericalRadiation.hpp
  )
