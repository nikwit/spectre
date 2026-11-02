// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestHelpers.hpp"
#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <optional>

#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ApplyTensorYlmFilter.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.tpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"

// Changing lmax must rebuild the angular matrices even when filter parameters
// are unchanged. Exercise both production volume and face paths, including a
// return to a previously used resolution, with the fused radial multiply.
SPECTRE_TEST_CASE("Unit.Gh.SphericalFilter.PRefinement", "[Unit][Evolution]") {
  using Tags = ylm::TensorYlm::filter_detail::gh_spacetime_vars_list;
  using Filter = Filters::SphericalShell<Tags>;
  const auto make_filter = []() {
    return Filter{4,    std::nullopt, 420,          true,        std::nullopt,
                  true, false,        std::nullopt, std::nullopt};
  };
  auto reused_volume = make_filter();
  auto reused_face = make_filter();
  for (const auto& extents : std::array<std::array<size_t, 2>, 5>{
           {{12, 6}, {14, 8}, {10, 5}, {12, 6}, {14, 6}}}) {
    const auto nr = extents[0];
    const auto lmax = extents[1];
    CAPTURE(nr, lmax);
    const Mesh<3> mesh{
        {nr, lmax + 1, 2 * lmax + 1},
        {Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
         Spectral::Basis::SphericalHarmonic},
        {Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::Gauss,
         Spectral::Quadrature::Equiangular}};
    for (const bool volume : {true, false}) {
      CAPTURE(volume);
      const size_t points = volume ? mesh.number_of_grid_points()
                                   : mesh.slice_away(0).number_of_grid_points();
      Variables<Tags> actual(points);
      for (size_t i = 0; i < actual.size(); ++i) {
        actual.data()[i] = sin(0.17 * static_cast<double>(i)) +
                           cos(0.031 * static_cast<double>(i));
      }
      auto expected = actual;
      InverseJacobian<DataVector, 3, Frame::Grid, Frame::Inertial> inverse(
          points, 0.0);
      Jacobian<DataVector, 3, Frame::Grid, Frame::Inertial> jacobian(points,
                                                                     0.0);
      for (size_t i = 0; i < 3; ++i) {
        inverse.get(i, i) = 1.0;
        jacobian.get(i, i) = 1.0;
      }
      const auto fresh = make_filter();
      if (volume) {
        reused_volume.apply_in_volume(make_not_null(&actual), mesh, inverse,
                                      jacobian);
        fresh.apply_in_volume(make_not_null(&expected), mesh, inverse,
                              jacobian);
      } else {
        reused_face.apply_on_boundary(make_not_null(&actual),
                                      mesh.slice_away(0), inverse, jacobian);
        fresh.apply_on_boundary(make_not_null(&expected), mesh.slice_away(0),
                                inverse, jacobian);
      }
      CHECK_VARIABLES_APPROX(actual, expected);
    }
  }
}
