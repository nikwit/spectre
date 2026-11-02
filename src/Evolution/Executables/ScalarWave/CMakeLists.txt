# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  Actions
  AmrEvents
  Charmxx::main
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  Events
  EventsAndDenseTriggers
  EventsAndTriggers
  Evolution
  Informer
  InitialDataUtilities
  LinearOperators
  MathFunctions
  ScalarWave
  Observer
  Options
  Parallel
  ParallelAmr
  PhaseControl
  Serialization
  Time
  Utilities
  WaveEquationSolutions
  )

function(add_scalar_wave_executable DIM AMR_IN_PHASE)
  set(EXECUTABLE "EvolveScalarWave${DIM}D")
  if(AMR_IN_PHASE)
    string(APPEND EXECUTABLE "AmrPhase")
  endif()
  add_spectre_executable(
    ${EXECUTABLE}
    EXCLUDE_FROM_ALL
    EvolveScalarWave.cpp
    )
  target_compile_definitions(
    ${EXECUTABLE}
    PRIVATE
    DIM=${DIM}
    AMR_IN_PHASE=${AMR_IN_PHASE}
    )
  target_link_libraries(${EXECUTABLE} PRIVATE ${LIBS_TO_LINK})
endfunction(add_scalar_wave_executable)

add_scalar_wave_executable(1 false)
add_scalar_wave_executable(2 false)
add_scalar_wave_executable(3 false)
# Exercise synchronized p refinement with ordinary testing criteria, using
# the same AMR phases as the GH SpecSpectral controller.
add_scalar_wave_executable(3 true)
