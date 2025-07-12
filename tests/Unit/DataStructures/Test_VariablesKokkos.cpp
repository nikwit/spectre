// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/VariablesKokkos.hpp"
#include "Helpers/DataStructures/TestTags.hpp"

SPECTRE_TEST_CASE("Unit.DataStructures.VariablesKokkos",
                  "[DataStructures][Unit]") {
  using VectorType = Kokkos::View<double*>;
  const size_t num_points = 5;
  using VectorTag = TestHelpers::Tags::Vector<VectorType>;
  using ScalarTag = TestHelpers::Tags::Scalar<VectorType>;
  Variables<tmpl::list<VectorTag, ScalarTag>> vars{num_points};
  CHECK(vars.number_of_grid_points() == num_points);
  CHECK(vars.size() == num_points * 4);
  Kokkos::parallel_for(
      "fill_vars", num_points, KOKKOS_LAMBDA(const size_t i) {
        auto vec = get<VectorTag>(vars);
        for (size_t d = 0; d < 3; ++d) {
          vec.get(d)[i] = static_cast<double>(d + i);
        }
        auto scal = get<ScalarTag>(vars);
        get(scal)[i] = static_cast<double>(i);
      });
  const auto vars_host = vars.create_mirror_view_and_copy(Kokkos::HostSpace{});
  for (size_t i = 0; i < num_points; ++i) {
    for (size_t d = 0; d < 3; ++d) {
      CHECK(
          get<::Tags::MirrorView<VectorTag, Kokkos::HostSpace>>(vars_host).get(
              d)[i] == static_cast<double>(d + i));
    }
    CHECK(get(get<::Tags::MirrorView<ScalarTag, Kokkos::HostSpace>>(
              vars_host))[i] == static_cast<double>(i));
  }
}
