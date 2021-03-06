# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY FunctionsOfTime)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  PiecewisePolynomial.cpp
  ReadSpecThirdOrderPiecewisePolynomial.cpp
  RegisterDerivedWithCharm.cpp
  SettleToConstant.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  FunctionOfTime.hpp
  OptionTags.hpp
  PiecewisePolynomial.hpp
  ReadSpecThirdOrderPiecewisePolynomial.hpp
  RegisterDerivedWithCharm.hpp
  SettleToConstant.hpp
  Tags.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Options
  Utilities
  PRIVATE
  ErrorHandling
  IO
  )
