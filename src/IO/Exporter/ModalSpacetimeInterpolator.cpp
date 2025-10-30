// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#include "IO/Exporter/ModalSpacetimeInterpolator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/ModalVector.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "IO/Exporter/Exporter.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel//Printf/Printf.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Overloader.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace spectre::Exporter {

namespace {

std::vector<std::string> resolve_filenames(
    const std::variant<std::vector<std::string>, std::string>&
        volume_files_or_glob) {
  std::vector<std::string> filenames =
      std::visit(Overloader{[](const std::vector<std::string>& volume_files) {
                              return volume_files;
                            },
                            [](const std::string& volume_files_glob) {
                              return file_system::glob(volume_files_glob);
                            }},
                 volume_files_or_glob);
  if (filenames.empty()) {
    ERROR("No volume files found. Specify at least one volume file.");
  }
  return filenames;
}

std::vector<std::pair<size_t, double>> load_observation_ids(
    const std::vector<std::string>& filenames,
    const std::string& subfile_name) {
  std::vector<std::pair<size_t, double>> obs_ids_and_times;
  const h5::H5File<h5::AccessType::ReadOnly> first_h5file(filenames.front());
  const auto& first_volfile = first_h5file.get<h5::VolumeData>(subfile_name);
  const auto all_observation_ids = first_volfile.list_observation_ids();
  if (all_observation_ids.empty()) {
    ERROR("No observation IDs found in the volume data files.");
  }
  obs_ids_and_times.reserve(all_observation_ids.size());
  for (const size_t obs_id : all_observation_ids) {
    const double obs_value = first_volfile.get_observation_value(obs_id);
    obs_ids_and_times.emplace_back(obs_id, obs_value);
  }
  first_h5file.close();
  std::ranges::sort(obs_ids_and_times,
                    [](const std::pair<size_t, double>& lhs,
                       const std::pair<size_t, double>& rhs) {
                      return lhs.second < rhs.second;
                    });

  const double cutoff_time = 500.0;
  obs_ids_and_times.erase(
      std::remove_if(obs_ids_and_times.begin(), obs_ids_and_times.end(),
                     [cutoff_time](const auto& id_and_time) {
                       return id_and_time.second < cutoff_time;
                     }),
      obs_ids_and_times.end());

  // Check that all other files have the same observation ids
  for (size_t file_index = 1; file_index < filenames.size(); ++file_index) {
    const h5::H5File<h5::AccessType::ReadOnly> other_h5file(
        filenames[file_index]);
    const auto& other_volfile = other_h5file.get<h5::VolumeData>(subfile_name);
    const auto other_observation_ids = other_volfile.list_observation_ids();
    if (other_observation_ids.size() != all_observation_ids.size()) {
      ERROR("Mismatched number of observation IDs between volume data files.");
    }
    for (size_t id_index = 0; id_index < all_observation_ids.size();
         ++id_index) {
      if (other_observation_ids[id_index] != all_observation_ids[id_index]) {
        ERROR("Mismatched observation IDs between volume data files.");
      }
    }
    other_h5file.close();
  }
  return obs_ids_and_times;
}

template <size_t Dim>
auto load_grids(const h5::VolumeData& volfile, const size_t obs_id) {
  const auto grid_names = volfile.get_grid_names(obs_id);
  const auto all_extents = volfile.get_extents(obs_id);
  const auto all_bases = volfile.get_bases(obs_id);
  const auto all_quadratures = volfile.get_quadratures(obs_id);
  std::vector<ElementId<Dim>> element_ids{};
  std::unordered_map<ElementId<Dim>, Mesh<Dim>> meshes{};
  element_ids.reserve(grid_names.size());
  for (const auto& grid_name : grid_names) {
    const ElementId<Dim> element_id(grid_name);
    element_ids.push_back(element_id);
    meshes[element_id] = h5::mesh_for_grid<Dim>(
        grid_name, grid_names, all_extents, all_bases, all_quadratures);
  }
  return std::make_pair(std::move(element_ids), std::move(meshes));
}

auto interpolant_for_mode(
    const std::vector<std::pair<size_t, double>>& obs_ids_and_times,
    const std::vector<double>& values, double relative_error) {
  const size_t num_observations = obs_ids_and_times.size();
  std::vector<size_t> selected_indices{{0, num_observations / 3,
                                        2 * num_observations / 3,
                                        num_observations - 1}};
  auto build_interpolant = [&obs_ids_and_times, &selected_indices, &values]() {
    std::vector<double> selected_times;
    std::vector<double> selected_values;
    selected_times.reserve(selected_indices.size());
    selected_values.reserve(selected_indices.size());
    for (const size_t index : selected_indices) {
      selected_times.push_back(obs_ids_and_times[index].second);
      selected_values.push_back(values[index]);
    }
    auto selected_times_copy = selected_times;
    auto selected_values_copy = selected_values;
    return boost::math::interpolators::pchip<std::vector<double>>(
        std::move(selected_times_copy), std::move(selected_values_copy));
  };
  auto compute_max_residual = [obs_ids_and_times, &values,
                               num_observations](const auto& interpolant) {
    double current_max_error = 0.0;
    size_t max_index = 0;
    for (size_t obs_index = 0; obs_index < num_observations; ++obs_index) {
      const double prediction =
          interpolant(obs_ids_and_times[obs_index].second);
      const double actual = values[obs_index];
      const double abs_error = std::abs(prediction - actual);
      if (abs_error > current_max_error) {
        current_max_error = abs_error;
        max_index = obs_index;
      }
    }
    return std::make_pair(current_max_error, max_index);
  };
  auto interpolant = build_interpolant();
  double max_error = std::numeric_limits<double>::infinity();

  while (true) {
    const auto [current_max_error, max_index] =
        compute_max_residual(interpolant);
    max_error = current_max_error;
    if (std::find(selected_indices.begin(), selected_indices.end(),
                  max_index) != selected_indices.end()) {
      // max_index is already included
      ERROR(
          "The interpolator tried to add an observation index that is "
          "already included. This indicates a logic error or numerical "
          "instability.");
    }
    if (max_error <= relative_error) {
      break;
    }
    selected_indices.push_back(max_index);
    std::sort(selected_indices.begin(), selected_indices.end());
    interpolant = build_interpolant();
    if (selected_indices.size() == num_observations - 1) {
      break;
    }
  }
  return std::make_tuple(interpolant, max_error, selected_indices.size());
}

}  // namespace

template <size_t Dim, typename Frame>
ModalSpacetimeInterpolator<Dim, Frame>::ModalSpacetimeInterpolator(
    std::variant<std::vector<std::string>, std::string> volume_files_or_glob,
    std::string subfile_name, std::vector<std::string> tensor_components,
    const double relative_error, const size_t max_interpolation_order)
    : volume_files_or_glob_(std::move(volume_files_or_glob)),
      tensor_components_(std::move(tensor_components)),
      relative_error_(relative_error),
      max_interpolation_order_(max_interpolation_order) {
  if (relative_error_ < 0.0) {
    ERROR("Relative error tolerance must be non-negative but is "
          << relative_error_ << ".");
  }
  const auto filenames = resolve_filenames(volume_files_or_glob_);
  const auto& obs_ids_and_times = load_observation_ids(filenames, subfile_name);
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();
  const h5::H5File<h5::AccessType::ReadOnly> first_h5file(filenames.front());
  const auto& first_volfile = first_h5file.get<h5::VolumeData>(subfile_name);
  const size_t reference_obs_id = obs_ids_and_times.back().first;
  auto serialized_domain = first_volfile.get_domain(reference_obs_id);
  domain_ = deserialize<Domain<Dim>>(serialized_domain->data());
  auto serialized_fots = first_volfile.get_functions_of_time(reference_obs_id);
  if (serialized_fots.has_value()) {
    functions_of_time_ =
        deserialize<domain::FunctionsOfTimeMap>(serialized_fots->data());
  } else {
    functions_of_time_.clear();
  }
  gather_element_metadata(filenames, subfile_name, reference_obs_id);
  build_interpolators(filenames, subfile_name, obs_ids_and_times);
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::gather_element_metadata(
    const std::vector<std::string>& filenames, const std::string& subfile_name,
    const size_t reference_obs_id) {
  // we do a first pass over all files to gather metadata about which elements
  // exist in which files, and what their meshes are. This avoids loading all
  // data into memory at once.
  for (size_t file_index = 0; file_index < filenames.size(); ++file_index) {
    h5::H5File<h5::AccessType::ReadOnly> h5file(filenames[file_index]);
    const auto& volfile = h5file.get<h5::VolumeData>(subfile_name);
    const auto [element_ids, meshes] =
        load_grids<Dim>(volfile, reference_obs_id);
    for (const auto& element_id : element_ids) {
      const ElementMetadata element_metadata{meshes.at(element_id), file_index};
      element_metadata_.emplace(element_id, element_metadata);
      element_search_trees_[element_id.block_id()].insert(element_id);
    }
  }
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::build_interpolators(
    const std::vector<std::string>& filenames, const std::string& subfile_name,
    const std::vector<std::pair<size_t, double>>& obs_ids_and_times) {
  // Build a time interpolant for every tensor component / grid point pair by
  // streaming the observations from disk element-by-element.
  const size_t num_observations = obs_ids_and_times.size();
  const size_t num_components = tensor_components_.size();

  const ElementId<Dim> reference_id("[B0,(L2I0,L2I3,L2I3)]");
  for (const auto& [_, tree] : element_search_trees_) {
    for (const auto& element_id : tree) {
      if (element_id != reference_id) {
        continue;
      }
      const auto& metadata = element_metadata_.at(element_id);
      ElementInterpolator element_interpolator{};
      element_interpolator.mesh = metadata.mesh;
      element_interpolator.component_interpolators.resize(num_components);
      const size_t num_grid_points =
          element_interpolator.mesh.number_of_grid_points();

      for (size_t component_index = 0; component_index < num_components;
           ++component_index) {
        auto& component_interpolator =
            element_interpolator.component_interpolators[component_index];
        const auto per_mode_values = load_component_time_series(
            metadata, element_id, component_index, filenames, subfile_name,
            obs_ids_and_times);
        component_interpolator.modal_interpolants.reserve(num_grid_points);
        for (size_t mode_index = 0; mode_index < num_grid_points;
             ++mode_index) {
          const auto& values = per_mode_values[mode_index];
          const auto [interpolant, max_error, num_abscissae] =
              interpolant_for_mode(obs_ids_and_times, values, relative_error_);
          if (max_error > relative_error_) {
            Parallel::printf(
                "For element %s, component %s, mode %zu, could not achieve "
                "the requested relative error tolerance %.3e; achieved "
                "maximum error %.3e using all %zu observation points.\n",
                get_output(element_id).c_str(),
                tensor_components_[component_index].c_str(), mode_index,
                relative_error_, max_error, num_observations);
          } else {
            Parallel::printf(
                "For element %s, component %s, mode %zu, achieved "
                "maximum error %.3e using %zu observation points.\n",
                get_output(element_id).c_str(),
                tensor_components_[component_index].c_str(), mode_index,
                max_error, num_abscissae);
          }
          component_interpolator.modal_interpolants.push_back(interpolant);
        }
      }
      Parallel::printf(
          "Constructed interpolator for element %s with %zu tensor components "
          "and %zu grid points per component.\n",
          get_output(element_id).c_str(), num_components, num_grid_points);
      interpolators_.emplace(element_id, std::move(element_interpolator));
    }
  }
}

template <size_t Dim, typename Frame>
std::vector<std::vector<double>>
ModalSpacetimeInterpolator<Dim, Frame>::load_component_time_series(
    const ElementMetadata& metadata, const ElementId<Dim>& element_id,
    const size_t component_index, const std::vector<std::string>& filenames,
    const std::string& subfile_name,
    const std::vector<std::pair<size_t, double>>& obs_ids_and_times) const {
  const size_t num_grid_points = metadata.mesh.number_of_grid_points();
  const size_t file_index = metadata.file_index;
  const size_t num_observations = obs_ids_and_times.size();
  // Every row is a modal coefficient, every column corresponds to an
  // observation.
  std::vector<std::vector<double>> per_mode_values(
      num_grid_points, std::vector<double>(num_observations, 0.0));

  const std::string element_name = get_output(element_id);
  h5::H5File<h5::AccessType::ReadOnly> element_file(filenames[file_index]);
  const auto& element_volfile = element_file.get<h5::VolumeData>(subfile_name);
  DataVector nodal_values(num_grid_points, 0.0);
  for (size_t obs_index = 0; obs_index < num_observations; ++obs_index) {
    const size_t obs_id = obs_ids_and_times[obs_index].first;
    const auto grid_names = element_volfile.get_grid_names(obs_id);
    ASSERT(
        std::find(grid_names.begin(), grid_names.end(), element_name) !=
            grid_names.end(),
        "Element "
            << element_id << " is not present in file index " << file_index
            << " for observation " << obs_id
            << ". Each element is expected to reside in the same volume file "
               "for all observations.");
    const auto extents = element_volfile.get_extents(obs_id);
    const auto [offset, length] =
        h5::offset_and_length_for_grid(element_name, grid_names, extents);
    const auto component_data =
        element_volfile
            .get_tensor_component(obs_id, tensor_components_[component_index])
            .data;
    if (std::holds_alternative<DataVector>(component_data)) {
      const auto& data = std::get<DataVector>(component_data);
      const double* element_data =
          data.data() + static_cast<std::ptrdiff_t>(offset);
      for (size_t i = 0; i < num_grid_points; ++i) {
        nodal_values[i] = element_data[i];
      }
    } else {
      const auto& data = std::get<std::vector<float>>(component_data);
      for (size_t i = 0; i < num_grid_points; ++i) {
        nodal_values[i] = static_cast<double>(data[offset + i]);
      }
    }
    const ModalVector modal_coefficients =
        to_modal_coefficients<Dim>(nodal_values, metadata.mesh);
    for (size_t mode = 0; mode < num_grid_points; ++mode) {
      per_mode_values[mode][obs_index] = modal_coefficients[mode];
    }
  }
  return per_mode_values;
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::interpolate_to_point(
    const gsl::not_null<std::vector<double>*> result,
    const tnsr::I<double, Dim, Frame>& target_point, const double time,
    const std::optional<gsl::not_null<std::vector<size_t>*>> block_order)
    const {
  if (UNLIKELY(interpolators_.empty())) {
    ERROR("ModalSpacetimeInterpolator has not been initialized.");
  }
  if (UNLIKELY(time < time_bounds_[0] || time > time_bounds_[1])) {
    ERROR("Requested time " << time
                            << " lies outside the available data interval "
                            << time_bounds_ << ".");
  }

  const auto block_logical_coords = block_logical_coordinates_single_point(
      target_point, domain_, time, functions_of_time_, block_order);
  if (UNLIKELY(not block_logical_coords.has_value())) {
    ERROR("Point is not in any block:\n" << target_point);
  }

  const auto& block_pair = block_logical_coords.value();

  const auto element_coords =
      element_logical_coordinates(block_pair, element_search_trees_);
  if (UNLIKELY(not element_coords.has_value())) {
    ERROR("Failed to determine element logical coordinates for point "
          << target_point << " at time " << time << ".");
  }

  const auto& element_id = element_coords->first;
  const auto& logical_coords = element_coords->second;

  const auto element_it = interpolators_.find(element_id);
  if (UNLIKELY(element_it == interpolators_.end())) {
    ERROR("No interpolator data found for element " << element_id << ".");
  }

  const auto& element_interpolator = element_it->second;
  const Mesh<Dim>& mesh = element_interpolator.mesh;
  const auto& component_interpolators =
      element_interpolator.component_interpolators;
  const size_t num_components = component_interpolators.size();

  if (UNLIKELY(num_components != tensor_components_.size())) {
    ERROR("Inconsistent number of tensor components stored in interpolator.");
  }

  const intrp::Irregular<Dim> spatial_interpolant(mesh, logical_coords);
  const size_t num_grid_points = mesh.number_of_grid_points();
  ModalVector modal_values(num_grid_points);
  DataVector nodal_values(num_grid_points);

  result->resize(num_components);
  for (size_t component_index = 0; component_index < num_components;
       ++component_index) {
    const auto& component_interpolator =
        component_interpolators[component_index];
    if (UNLIKELY(component_interpolator.modal_interpolants.size() !=
                 num_grid_points)) {
      ERROR("Stored modal interpolants do not match mesh size for element "
            << element_id << ".");
    }
    for (size_t point = 0; point < num_grid_points; ++point) {
      modal_values[point] =
          component_interpolator.modal_interpolants[point](time);
    }
    to_nodal_coefficients<Dim>(make_not_null(&nodal_values), modal_values,
                               mesh);
    const gsl::span<const double> input_span(nodal_values.data(),
                                             nodal_values.size());
    gsl::span<double> output_span(&(*result)[component_index], 1);
    spatial_interpolant.interpolate(make_not_null(&output_span), input_span);
  }
}

// Explicit instantiations

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data) \
  template class ModalSpacetimeInterpolator<DIM(data), FRAME(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (Frame::Inertial))

#undef INSTANTIATE
#undef DIM
#undef FRAME

}  // namespace spectre::Exporter
