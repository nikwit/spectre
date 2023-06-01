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
  const double radial_shift = -0.1;
  const std::array<DataVector, deriv_order + 1> init_func{
      {{radial_shift}, {0.0}, {0.0}}};
  using Polynomial = domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>;
  using FoftPtr = std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>;
  std::unordered_map<std::string, FoftPtr> f_of_t_list{};
  f_of_t_list["Expansion"] = std::make_unique<Polynomial>(0.1, init_func, 20.);
  const FoftPtr& expansion_a_base = f_of_t_list.at("Expansion");
  const double inner_radius = 1.5;
  const double outer_radius = 5.;
  const CoordinateMaps::TimeDependent::WorldtubeExpansion expansion_map(
      inner_radius, outer_radius, "Expansion");
  const double time = 1.1;

  const std::array<double, 3> test_point_inner{{0.1, 0.2, 0.3}};
  const auto mapped_test_point_inner =
      expansion_map(test_point_inner, time, f_of_t_list);
  CHECK_ITERABLE_APPROX(test_point_inner, mapped_test_point_inner);

  const std::array<double, 3> test_point_outer{{5., 0.2, 0.3}};
  const double outer_point_radius = magnitude(test_point_outer);
  const auto mapped_outer_point =
      expansion_map(test_point_outer, time, f_of_t_list);
  const double mapped_radius_outer = magnitude(mapped_outer_point);
  CHECK(mapped_radius_outer - radial_shift == approx(outer_point_radius));

  const std::array<double, 3> test_point_trans{{2., 0.2, 0.3}};
  const double trans_point_radius = magnitude(test_point_trans);
  const auto mapped_point_trans =
      expansion_map(test_point_trans, time, f_of_t_list);
  const double mapped_radius_trans = magnitude(mapped_point_trans);
  CHECK(mapped_radius_trans <= trans_point_radius);
  CHECK(mapped_radius_trans >= trans_point_radius + radial_shift);
  test_inverse_map(expansion_map, test_point_inner, 4., f_of_t_list);
  test_inverse_map(expansion_map, test_point_outer, 4., f_of_t_list);
  test_inverse_map(expansion_map, test_point_trans, 4., f_of_t_list);
  test_jacobian(expansion_map, test_point_inner, 4., f_of_t_list);
  test_jacobian(expansion_map, test_point_outer, 4., f_of_t_list);
  test_jacobian(expansion_map, test_point_trans, 4., f_of_t_list);
  test_inverse_map(expansion_map, test_point_trans, 4., f_of_t_list);
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.TimeDependent.WorldtubeExpansion",
                  "[Domain][Unit]") {
  test();
}
}  // namespace domain
