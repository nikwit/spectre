// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/TimeDependent/WorldtubeExpansion.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Domain/CoordinateMaps/TestMapHelpers.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TypeTraits.hpp"

namespace domain {
namespace {
void test() {
  INFO("Map");
  static constexpr size_t deriv_order = 2;
  const double initial_time = -0.5;
  double t = -0.4;
  const double dt = 0.6;
  const std::array<DataVector, deriv_order + 1> init_func{
      {{1.1}, {0.0}, {0.0}}};
  using Polynomial = domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>;
  using FoftPtr = std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>;
  std::unordered_map<std::string, FoftPtr> f_of_t_list{};
  f_of_t_list["Expansion"] = std::make_unique<Polynomial>(0.1, init_func, 20.);
  const FoftPtr& expansion_a_base = f_of_t_list.at("Expansion");
  const double inner_radius = 1.5;
  const double outer_radius = 5.;
  const CoordinateMaps::TimeDependent::WorldtubeExpansion expansion_map(
      inner_radius, outer_radius, "Expansion");
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.TimeDependent.WorldtubeExpansion",
                  "[Domain][Unit]") {
  test();
}
}  // namespace domain
