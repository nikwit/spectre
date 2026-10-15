# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Actions/EvolvePhysicalTrace.cpp
  ComputeRecordingData.cpp
  Events/ObserveWorldtubeMatching.cpp
  FaceReplay.cpp
  KretschmannFaceData.cpp
  Matching.cpp
  PhysicalBoundaryEvolution.cpp
  PhysicalBoundaryTrace.cpp
  QhatRadialGauge.cpp
  SphereReplay.cpp
  WeylCurvature.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Actions/EvolvePhysicalTrace.hpp
  Actions/UpdateKretschmannFaceData.hpp
  ComputeExcisionRecordingData.hpp
  ComputeRecordingData.hpp
  Events/ObserveWorldtubeMatching.hpp
  FaceReplay.hpp
  KretschmannFaceData.hpp
  Matching.hpp
  PhysicalBoundaryEvolution.hpp
  PhysicalBoundaryTrace.hpp
  QhatRadialGauge.hpp
  RadialGauge.hpp
  SphereReplay.hpp
  WeylCurvature.hpp
  )
