# Distributed under the MIT License.
# See LICENSE.txt for details.

# Since benchmarking is only interesting in release mode the executable isn't
# added for Debug builds. Charm++'s main function is overridden with the main
# from the Google Benchmark library. The executable is not added to the `all` make
# target since it is only interesting in specific circumstances.
if("${GoogleBenchmark_FOUND}" AND NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(executable Benchmark)

  add_spectre_executable(
    ${executable}
    EXCLUDE_FROM_ALL
    Benchmark.cpp
    )

  # Add specific libraries needed for the benchmark you are interested in.
  target_link_libraries(
    ${executable}
    PRIVATE
    CoordinateMaps
    Domain
    Informer
    GoogleBenchmark
    LinearOperators
    Spectral
    )
endif()

if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(executable GhTimeDerivativeBenchmark)

  add_spectre_executable(
    ${executable}
    EXCLUDE_FROM_ALL
    GhTimeDerivativeBenchmark.cpp
    )

  target_link_libraries(
    ${executable}
    PRIVATE
    DataStructures
    GeneralizedHarmonic
    GeneralRelativity
    GeneralRelativitySolutions
    Spectral
    Utilities
    )

  set(executable PartialDerivativesBenchmark)

  add_spectre_executable(
    ${executable}
    EXCLUDE_FROM_ALL
    PartialDerivativesBenchmark.cpp
    )

  target_link_libraries(
    ${executable}
    PRIVATE
    DataStructures
    LinearOperators
    Spectral
    Utilities
    )

  set(executable SphericalTransformBenchmark)

  add_spectre_executable(
    ${executable}
    EXCLUDE_FROM_ALL
    SphericalTransformBenchmark.cpp
    )

  target_link_libraries(
    ${executable}
    PRIVATE
    DataStructures
    Libsharp
    SphericalHarmonics
    Utilities
    )

  set(executable ShapeJacobianBenchmark)

  add_spectre_executable(
    ${executable}
    EXCLUDE_FROM_ALL
    ShapeJacobianBenchmark.cpp
    )

  target_link_libraries(
    ${executable}
    PRIVATE
    DataStructures
    Domain
    FunctionsOfTime
    SphericalHarmonics
    Utilities
    )
endif()
