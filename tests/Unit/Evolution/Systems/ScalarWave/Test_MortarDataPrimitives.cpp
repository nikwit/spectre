// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/MortarDataPrimitives.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

namespace {

constexpr size_t Dim = 3;
using package_field_tags =
    typename ScalarWave::BoundaryCorrections::UpwindPenalty<
        Dim>::dg_package_field_tags;
using device_package_field_tags =
    db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;

void fill_with_linear_data(
    const gsl::not_null<Variables<package_field_tags>*> variables) {
  for (size_t i = 0; i < variables->size(); ++i) {
    variables->data()[i] = 0.1 + static_cast<double>(i);
  }
}

void test_orient_each_component_device_matches_host() {
  const Mesh<2> mortar_mesh{
      {{4, 3}}, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto};
  Variables<package_field_tags> host_vars{mortar_mesh.number_of_grid_points()};
  fill_with_linear_data(make_not_null(&host_vars));

  const size_t sliced_dim = 0;
  const OrientationMap<3> orientation{std::array{Direction<3>::upper_xi(),
                                                 Direction<3>::lower_eta(),
                                                 Direction<3>::upper_zeta()}};

  DataVector host_vars_raw{host_vars.size()};
  std::copy_n(host_vars.data(), host_vars.size(), host_vars_raw.data());
  const DataVector expected_oriented_raw = orient_variables_on_slice(
      host_vars_raw, mortar_mesh.extents(), sliced_dim, orientation);

  std::vector<double> grid_point_indices(mortar_mesh.number_of_grid_points());
  std::iota(grid_point_indices.begin(), grid_point_indices.end(), 0.0);
  const std::vector<double> oriented_indices = orient_variables_on_slice(
      grid_point_indices, mortar_mesh.extents(), sliced_dim, orientation);

  Kokkos::View<size_t*> oriented_source_index(
      "OrientedMortarGridPointSourceIndex", oriented_indices.size());
  auto host_oriented_source_index =
      Kokkos::create_mirror_view(oriented_source_index);
  for (size_t i = 0; i < oriented_indices.size(); ++i) {
    host_oriented_source_index(i) = static_cast<size_t>(oriented_indices[i]);
  }
  Kokkos::deep_copy(oriented_source_index, host_oriented_source_index);

  const auto device_vars = copy_to_device(host_vars);
  Variables<device_package_field_tags> oriented_device_vars{
      mortar_mesh.number_of_grid_points()};
  ScalarWave::Kokkos::orient_each_component_device(
      make_not_null(&oriented_device_vars), device_vars, oriented_source_index);

  Variables<package_field_tags> oriented_host_vars{
      mortar_mesh.number_of_grid_points()};
  copy_to_host(make_not_null(&oriented_host_vars), oriented_device_vars);
  DataVector oriented_host_raw{oriented_host_vars.size()};
  std::copy_n(oriented_host_vars.data(), oriented_host_vars.size(),
              oriented_host_raw.data());

  CHECK_ITERABLE_APPROX(oriented_host_raw, expected_oriented_raw);
}

void test_project_to_and_from_mortar_device_match_host() {
  const Mesh<2> face_mesh{
      {{3, 4}}, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto};
  const Mesh<2> mortar_mesh{
      {{4, 4}}, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto};
  const std::array<Spectral::SegmentSize, 2> mortar_size{
      Spectral::SegmentSize::UpperHalf, Spectral::SegmentSize::Full};

  Variables<package_field_tags> face_host_vars{
      face_mesh.number_of_grid_points()};
  fill_with_linear_data(make_not_null(&face_host_vars));

  const Variables<package_field_tags> projected_to_mortar_host =
      ::dg::project_to_mortar(face_host_vars, face_mesh, mortar_mesh,
                              mortar_size);

  const auto face_device_vars = copy_to_device(face_host_vars);
  Variables<device_package_field_tags> projected_to_mortar_device{
      mortar_mesh.number_of_grid_points()};
  ScalarWave::Kokkos::project_to_mortar_device(
      make_not_null(&projected_to_mortar_device), face_device_vars, face_mesh,
      mortar_mesh, mortar_size);
  Variables<package_field_tags> projected_to_mortar_host_from_device{
      mortar_mesh.number_of_grid_points()};
  copy_to_host(make_not_null(&projected_to_mortar_host_from_device),
               projected_to_mortar_device);
  CHECK_ITERABLE_APPROX(
      gsl::make_span(projected_to_mortar_host_from_device.data(),
                     projected_to_mortar_host_from_device.size()),
      gsl::make_span(projected_to_mortar_host.data(),
                     projected_to_mortar_host.size()));

  const Variables<package_field_tags> projected_from_mortar_host =
      ::dg::project_from_mortar(projected_to_mortar_host, face_mesh,
                                mortar_mesh, mortar_size);
  const auto projected_to_mortar_device_input =
      copy_to_device(projected_to_mortar_host);
  Variables<device_package_field_tags> projected_from_mortar_device{
      face_mesh.number_of_grid_points()};
  ScalarWave::Kokkos::project_from_mortar_device(
      make_not_null(&projected_from_mortar_device),
      projected_to_mortar_device_input, face_mesh, mortar_mesh, mortar_size);
  Variables<package_field_tags> projected_from_mortar_host_from_device{
      face_mesh.number_of_grid_points()};
  copy_to_host(make_not_null(&projected_from_mortar_host_from_device),
               projected_from_mortar_device);
  CHECK_ITERABLE_APPROX(
      gsl::make_span(projected_from_mortar_host_from_device.data(),
                     projected_from_mortar_host_from_device.size()),
      gsl::make_span(projected_from_mortar_host.data(),
                     projected_from_mortar_host.size()));
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.ScalarWave.MortarDataPrimitives",
                  "[Unit][Evolution]") {
  test_orient_each_component_device_matches_host();
  test_project_to_and_from_mortar_device_match_host();
}
