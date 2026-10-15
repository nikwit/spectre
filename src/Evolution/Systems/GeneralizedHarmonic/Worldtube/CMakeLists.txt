# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ComputeRecordingData.cpp
  Events/ObserveWorldtubeMatching.cpp
  FaceReplay.cpp
  KretschmannFaceData.cpp
  Matching.cpp
  WeylCurvature.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ComputeExcisionRecordingData.hpp
  ComputeRecordingData.hpp
  Actions/UpdateKretschmannFaceData.hpp
  Events/ObserveWorldtubeMatching.hpp
  FaceReplay.hpp
  KretschmannFaceData.hpp
  RadialGauge.hpp
  Matching.hpp
  WeylCurvature.hpp
  )
