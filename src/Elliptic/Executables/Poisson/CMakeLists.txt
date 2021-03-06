# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  Convergence
  CoordinateMaps
  DiscontinuousGalerkin
  DomainCreators
  Elliptic
  EllipticDg
  Events
  Informer
  IO
  LinearOperators
  MathFunctions
  Options
  Parallel
  ParallelLinearSolver
  Poisson
  PoissonSolutions
  Utilities
  )

function(add_poisson_executable DIM)
  add_spectre_parallel_executable(
    "SolvePoisson${DIM}D"
    SolvePoisson
    Elliptic/Executables/Poisson
    "Metavariables<${DIM}>"
    "${LIBS_TO_LINK}"
    )
endfunction(add_poisson_executable)

add_poisson_executable(1)
add_poisson_executable(2)
add_poisson_executable(3)
