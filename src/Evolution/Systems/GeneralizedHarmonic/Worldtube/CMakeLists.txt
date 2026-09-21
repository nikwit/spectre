# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Events/ObserveWorldtubeMatching.cpp
  KretschmannFaceData.cpp
  Matching.cpp
  WeylCurvature.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Actions/UpdateKretschmannFaceData.hpp
  Events/ObserveWorldtubeMatching.hpp
  KretschmannFaceData.hpp
  Matching.hpp
  WeylCurvature.hpp
  )
