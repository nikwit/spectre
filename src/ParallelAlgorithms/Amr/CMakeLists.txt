# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ParallelAmr)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  MeshPlan.cpp
  Tags.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  MeshPlan.hpp
  CriterionData.hpp
  RedistributeElements.hpp
  Tags.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Domain
  DomainStructure
  AmrPolicies
  INTERFACE
  Amr
  AmrCriteria
  AmrPolicies
  AmrProjectors
  Initialization
  Logging
  Printf
  Serialization
  Utilities
  )

add_subdirectory(Actions)
add_subdirectory(Criteria)
add_subdirectory(Events)
add_subdirectory(Policies)
add_subdirectory(Projectors)
add_subdirectory(Protocols)
