// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "IO/Exporter/ModalSpacetimeInterpolator.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace spectre::Exporter {

namespace {

void write_test_volume_data(const std::string& h5_file_name,
                            const std::vector<double>& times,
                            const Domain<3>& domain, const Mesh<3>& mesh,
                            const std::vector<ElementId<3>>& element_ids) {
  const auto logical_coords = logical_coordinates(mesh);

  h5::H5File<h5::AccessType::ReadWrite> h5_file{h5_file_name};
  auto& volume_file = h5_file.insert<h5::VolumeData>("/VolumeData", 0);

  size_t obs_id = 0;
  for (const double time : times) {
    std::vector<ElementVolumeData> element_volume_data{};
    element_volume_data.reserve(element_ids.size());

    for (const auto& element_id : element_ids) {
      const ElementMap<3, Frame::Inertial> element_map{
          element_id, domain.blocks()[element_id.block_id()]};
      const auto inertial_coords = element_map(logical_coords);

      DataVector psi = get<0>(inertial_coords);
      psi += 2.0 * get<1>(inertial_coords);
      psi += 3.0 * get<2>(inertial_coords);
      psi += cos(time);

      DataVector phi = 0.5 * get<0>(inertial_coords);
      phi -= get<2>(inertial_coords);
      phi += 2.0 * sin(time) + 1.0;

      element_volume_data.push_back(
          ElementVolumeData{element_id,
                            {TensorComponent{"Psi", std::move(psi)},
                             TensorComponent{"Phi", std::move(phi)}},
                            mesh});
    }

    volume_file.write_volume_data(obs_id, time, element_volume_data,
                                  serialize(domain));
    ++obs_id;
  }
}

std::array<double, 2> expected_values(
    const tnsr::I<double, 3, Frame::Inertial>& x, const double time) {
  const double psi = cos(time) + x.get(0) + 2.0 * x.get(1) + 3.0 * x.get(2);
  const double phi = (2.0 * sin(time) + 1.0) + 0.5 * x.get(0) - x.get(2);
  return {psi, phi};
}

}  // namespace

SPECTRE_TEST_CASE("Unit.IO.Exporter.ModalSpacetimeInterpolator", "[Unit]") {
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  const std::string h5_file_name_1{
      "Unit.IO.Exporter.ModalSpacetimeInterpolator.1.h5"};
  const std::string h5_file_name_2{
      "Unit.IO.Exporter.ModalSpacetimeInterpolator.2.h5"};
  if (file_system::check_if_file_exists(h5_file_name_1)) {
    file_system::rm(h5_file_name_1, true);
  }
  if (file_system::check_if_file_exists(h5_file_name_2)) {
    file_system::rm(h5_file_name_2, true);
  }

  const domain::creators::Brick domain_creator{
      {{0.0, 0.0, 0.0}},
      {{1.0, 1.0, 1.0}},
      {{1, 0, 0}},
      {{3, 3, 3}}};
  const auto domain = domain_creator.create_domain();
  const auto all_element_ids =
      initial_element_ids(domain_creator.initial_refinement_levels());
  const Mesh<3> mesh{4, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  std::vector<ElementId<3>> element_ids_file_1{};
  std::vector<ElementId<3>> element_ids_file_2{};
  element_ids_file_1.reserve(all_element_ids.size() / 2 + 1);
  element_ids_file_2.reserve(all_element_ids.size() / 2 + 1);
  for (size_t i = 0; i < all_element_ids.size(); ++i) {
    if (i % 2 == 0) {
      element_ids_file_1.push_back(all_element_ids[i]);
    } else {
      element_ids_file_2.push_back(all_element_ids[i]);
    }
  }

  const double final_time = 4.0;
  std::vector<double> times{};
  const double step_size = 0.01;
  for (double time = 0.; time <= final_time; time += step_size) {
    times.push_back(time);
  }
  write_test_volume_data(h5_file_name_1, times, domain, mesh,
                         element_ids_file_1);
  write_test_volume_data(h5_file_name_2, times, domain, mesh,
                         element_ids_file_2);

  const tnsr::I<double, 3, Frame::Inertial> target_point{{0.2, 0.3, 0.4}};
  const std::vector<double> query_times{0.0, 0.12345, 1.7654, 2.75345};
  const std::vector<double> tolerances{1.0e-4, 1.0e-8, 1.0e-12};
  const std::vector<std::string> volume_files{h5_file_name_1, h5_file_name_2};

  for (const double tolerance : tolerances) {
    CAPTURE(tolerance);
    ModalSpacetimeInterpolator<3, Frame::Inertial> interpolator(
        volume_files, "VolumeData", {"Psi", "Phi"}, tolerance);
    std::vector<double> result{};
    for (const double time : query_times) {
      CAPTURE(time);
      const auto expected = expected_values(target_point, time);
      interpolator.interpolate_to_point(make_not_null(&result), target_point,
                                        time);
      CHECK(result.size() == expected.size());
      CHECK(result[0] == approx(expected[0]).epsilon(tolerance));
      CHECK(result[1] == approx(expected[1]).epsilon(tolerance));
      CHECK_FALSE(std::abs(result[0] - expected[0]) >
                  tolerance * std::max(std::abs(expected[0]), 1.0));
      CHECK_FALSE(std::abs(result[1] - expected[1]) >
                  tolerance * std::max(std::abs(expected[1]), 1.0));
    }
  }

  if (file_system::check_if_file_exists(h5_file_name_1)) {
    file_system::rm(h5_file_name_1, true);
  }
  if (file_system::check_if_file_exists(h5_file_name_2)) {
    file_system::rm(h5_file_name_2, true);
  }
}

}  // namespace spectre::Exporter
