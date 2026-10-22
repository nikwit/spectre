# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  CalculateGeodesicPunctureField.cpp
  InitializeConstraintGammas.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  CalculateGeodesicPunctureField.hpp
  InitializeConstraintGammas.hpp
  InitializeCurrentIteration.hpp
  IteratePunctureField.hpp
  ReceiveWorldtubeData.hpp
  SendToWorldtube.hpp
)
