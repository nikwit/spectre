// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "IO/Exporter/ModalSpacetimeInterpolator.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel//Printf/Printf.hpp"
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

struct ValidationElement {
  ElementId<3> element_id;
  Mesh<3> mesh;
  DataVector nodal_values;
};

void validate_against_reference_data() {
  const std::string volume_file_path{
      "/Users/niko/caltech/simulations/spacetime-interpolator/BbhVolume0.h5"};
  const std::string volume_file_path2{
      "/Users/niko/caltech/simulations/spacetime-interpolator/BbhVolume1.h5"};

  h5::H5File<h5::AccessType::ReadOnly> h5_file(volume_file_path);
  const auto& sparse_volume = h5_file.get<h5::VolumeData>("/VolumeData");

  auto sparse_observation_ids = sparse_volume.list_observation_ids();
  std::ranges::sort(sparse_observation_ids,
                    [&sparse_volume](const size_t lhs, const size_t rhs) {
                      return sparse_volume.get_observation_value(lhs) <
                             sparse_volume.get_observation_value(rhs);
                    });

  const auto serialized_functions =
      sparse_volume.get_functions_of_time(sparse_observation_ids.back());
  domain::FunctionsOfTimeMap functions_of_time{};
  if (serialized_functions.has_value()) {
    functions_of_time =
        deserialize<domain::FunctionsOfTimeMap>(serialized_functions->data());
  }
  ModalSpacetimeInterpolator<3, Frame::Inertial> interpolator(
      "serialized_interpolator.h5", "/Interpolator");
  /*ModalSpacetimeInterpolator<3, Frame::Inertial> interpolator(
      std::vector<std::string>{volume_file_path, volume_file_path2},
      std::vector<std::string>{"VerySparseModal", "SparseModal", "FullModal"},
      {"Lapse"}, 1.0e-8);
  interpolator.write_to_h5("serialized_interpolator.h5", "/Interpolator");*/
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> logical_dist(-1.0, 1.0);
  h5_file.close_current_object();

  const auto& validation_volume = h5_file.get<h5::VolumeData>("/VolumeData");
  auto validation_observation_ids = validation_volume.list_observation_ids();

  const double cutoff_time_front = 0.0;
  const double cutoff_time_back = 1500.;
  validation_observation_ids.erase(
      std::remove_if(validation_observation_ids.begin(),
                     validation_observation_ids.end(),
                     [cutoff_time_front, cutoff_time_back,
                      &validation_volume](const auto& id) {
                       return (validation_volume.get_observation_value(id) <
                               cutoff_time_front) ||
                              (validation_volume.get_observation_value(id) >
                               cutoff_time_back);
                     }),
      validation_observation_ids.end());

  const auto serialized_domain =
      validation_volume.get_domain(validation_observation_ids.back());
  REQUIRE(serialized_domain.has_value());
  const Domain<3> file_domain =
      deserialize<Domain<3>>(serialized_domain->data());
  for (const size_t validation_obs_id : validation_observation_ids) {
    const double validation_time =
        validation_volume.get_observation_value(validation_obs_id);
    const auto grid_names = validation_volume.get_grid_names(validation_obs_id);
    const auto extents = validation_volume.get_extents(validation_obs_id);
    const auto bases = validation_volume.get_bases(validation_obs_id);
    const auto quadratures =
        validation_volume.get_quadratures(validation_obs_id);
    const auto tensor_data =
        validation_volume.get_tensor_component(validation_obs_id, "Lapse");
    // REQUIRE(std::holds_alternative<DataVector>(tensor_data.data));
    const auto& validation_data = std::get<DataVector>(tensor_data.data);
    const auto reference_element_id = ElementId<3>("[B0,(L2I0,L2I3,L2I3)]");
    std::vector<ValidationElement> validation_elements{};
    validation_elements.reserve(grid_names.size());
    for (size_t grid_index = 0; grid_index < grid_names.size(); ++grid_index) {
      const ElementId<3> element_id(grid_names[grid_index]);
      if (element_id.block_id() != 0) {
        continue;
      }
      std::array<size_t, 3> extent_array{};
      std::array<Spectral::Basis, 3> basis_array{};
      std::array<Spectral::Quadrature, 3> quadrature_array{};
      for (size_t dim = 0; dim < 3; ++dim) {
        extent_array[dim] = extents[grid_index][dim];
        basis_array[dim] = bases[grid_index][dim];
        quadrature_array[dim] = quadratures[grid_index][dim];
      }
      const auto [offset, length] = h5::offset_and_length_for_grid(
          grid_names[grid_index], grid_names, extents);
      DataVector nodal_values(length);
      for (size_t i = 0; i < length; ++i) {
        nodal_values[i] = validation_data[offset + i];
      }
      validation_elements.push_back(ValidationElement{
          element_id, Mesh<3>{extent_array, basis_array, quadrature_array},
          std::move(nodal_values)});
    }

    REQUIRE_FALSE(validation_elements.empty());
    std::uniform_int_distribution<size_t> element_dist(
        0, validation_elements.size() - 1);
    const double relative_tolerance = 1.0e-4;

    for (size_t sample_index = 0; sample_index < 2; ++sample_index) {
      const ValidationElement& element_data =
          validation_elements[element_dist(generator)];
      tnsr::I<double, 3, Frame::ElementLogical> logical_point{
          {logical_dist(generator), logical_dist(generator),
           logical_dist(generator)}};
      const ElementMap<3, Frame::Inertial> element_map(
          element_data.element_id,
          file_domain.blocks()[element_data.element_id.block_id()]);
      const auto inertial_point =
          element_map(logical_point, validation_time, functions_of_time);

      std::vector<double> interpolated_values{};
      interpolator.interpolate_to_point(make_not_null(&interpolated_values),
                                        inertial_point, validation_time);
      REQUIRE(interpolated_values.size() == 1);

      const intrp::Irregular<3> irregular(element_data.mesh, logical_point);
      const gsl::span<const double> nodal_values(
          element_data.nodal_values.data(), element_data.nodal_values.size());
      double validation_value = 0.0;
      gsl::span<double> output_span(&validation_value, 1);
      irregular.interpolate(make_not_null(&output_span), nodal_values);
      Parallel::printf("Relative error: %e\n",
                       std::abs(interpolated_values[0] - validation_value) /
                           std::abs(validation_value));
      CHECK(interpolated_values[0] ==
            approx(validation_value).epsilon(relative_tolerance));
    }
  }
}

}  // namespace

SPECTRE_TEST_CASE("Unit.IO.Exporter.ModalSpacetimeInterpolator", "[Unit]") {
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  /*const std::string h5_file_name_1{
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
      {{1, 2, 1}},
      {{4, 4, 4}}};
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
  }*/

  validate_against_reference_data();
}

}  // namespace spectre::Exporter
