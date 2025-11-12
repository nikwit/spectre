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
#include "IO/H5/CheckH5.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "IO/H5/Wrappers.hpp"
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
#include "Utilities/PrettyType.hpp"
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

  const double cutoff_time = 0.0;
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

struct ModeInterpolationResult {
  boost::math::interpolators::pchip<std::vector<double>> interpolant;
  double max_error = std::numeric_limits<double>::quiet_NaN();
  size_t num_abscissae = 0;
  std::vector<double> nodes{};
  std::vector<double> values{};
};

ModeInterpolationResult interpolant_for_mode(
    const std::vector<std::pair<size_t, double>>& obs_ids_and_times,
    const std::vector<double>& values, double absolute_error) {
  const size_t num_observations = obs_ids_and_times.size();
  std::vector<size_t> selected_indices{{0, num_observations / 3,
                                        2 * num_observations / 3,
                                        num_observations - 1}};

  std::vector<size_t> largest_residual_indices(num_observations);
  std::iota(largest_residual_indices.begin(), largest_residual_indices.end(),
            0);
  auto build_interpolant = [&obs_ids_and_times,
                            &values](std::vector<size_t>& current_indices) {
    std::vector<double> selected_times;
    std::vector<double> selected_values;
    selected_times.reserve(current_indices.size());
    selected_values.reserve(current_indices.size());
    for (const size_t index : current_indices) {
      selected_times.push_back(obs_ids_and_times[index].second);
      selected_values.push_back(values[index]);
    }
    auto selected_times_copy = selected_times;
    auto selected_values_copy = selected_values;
    return boost::math::interpolators::pchip<std::vector<double>>(
        std::move(selected_times_copy), std::move(selected_values_copy));
  };
  auto interpolant = build_interpolant(selected_indices);
  double max_error = std::numeric_limits<double>::max();
  while (true) {
    std::vector<double> residuals(num_observations);
    std::transform(
        obs_ids_and_times.begin(), obs_ids_and_times.end(), values.begin(),
        residuals.begin(),
        [&interpolant](const auto& id_and_time, const double actual_value) {
          const double predicted_value = interpolant(id_and_time.second);
          return std::abs(predicted_value - actual_value);
        });
    const size_t batch_size = 1;
    std::partial_sort(largest_residual_indices.begin(),
                      largest_residual_indices.begin() + batch_size,
                      largest_residual_indices.end(),
                      [&residuals](const size_t lhs, const size_t rhs) {
                        return residuals[lhs] > residuals[rhs];
                      });
    const size_t max_index = largest_residual_indices.front();
    max_error = residuals[max_index];
    if (max_error <= absolute_error) {
      // achieved the desired accuracy
      break;
    }
    for (size_t i = 0; i < batch_size; ++i) {
      if (std::find(selected_indices.begin(), selected_indices.end(),
                    largest_residual_indices[i]) == selected_indices.end()) {
        selected_indices.push_back(largest_residual_indices[i]);
      }
    }
    std::ranges::sort(selected_indices);
    if (std::adjacent_find(selected_indices.begin(), selected_indices.end()) !=
        selected_indices.end()) {
      // duplicate indices, cannot build a valid interpolant
      ERROR("Failed to build interpolant: duplicate indices selected.");
    }
    interpolant = build_interpolant(selected_indices);
    if (selected_indices.size() == num_observations - 1) {
      // all points have been selected, the max error is still too large
      break;
    }
  }
  std::vector<double> selected_times;
  std::vector<double> selected_values;
  selected_times.reserve(selected_indices.size());
  selected_values.reserve(selected_indices.size());
  for (const size_t index : selected_indices) {
    selected_times.push_back(obs_ids_and_times[index].second);
    selected_values.push_back(values[index]);
  }
  return ModeInterpolationResult{
      std::move(interpolant), max_error, selected_indices.size(),
      std::move(selected_times), std::move(selected_values)};
}

}  // namespace

template <size_t Dim, typename Frame>
ModalSpacetimeInterpolator<Dim, Frame>::ModalSpacetimeInterpolator(
    std::variant<std::vector<std::string>, std::string> volume_files_or_glob,
    std::vector<std::string> subfile_names,
    std::vector<std::string> tensor_components, const double absolute_error)
    : volume_files_or_glob_(std::move(volume_files_or_glob)),
      tensor_components_(std::move(tensor_components)),
      absolute_error_(absolute_error) {
  if (absolute_error < 0.0) {
    ERROR("Absolute error tolerance must be non-negative but is "
          << absolute_error_ << ".");
  }
  const auto filenames = resolve_filenames(volume_files_or_glob_);
  const std::string last_subfile_name = subfile_names.back();
  const auto& obs_ids_and_times =
      load_observation_ids(filenames, last_subfile_name);
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();
  const h5::H5File<h5::AccessType::ReadOnly> first_h5file(filenames.front());
  const auto& first_volfile =
      first_h5file.get<h5::VolumeData>(last_subfile_name);
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
  gather_element_metadata(filenames, last_subfile_name, reference_obs_id);
  build_interpolators(filenames, subfile_names);
}

template <size_t Dim, typename Frame>
ModalSpacetimeInterpolator<Dim, Frame>::ModalSpacetimeInterpolator(
    const std::string& h5_filename, const std::string& group_path)
    : volume_files_or_glob_(std::vector<std::string>{}) {
  hid_t file_id =
      H5Fopen(h5_filename.c_str(), H5F_ACC_RDONLY, h5::h5p_default());
  CHECK_H5(file_id, "Failed to open H5 file '" << h5_filename << "'.");

  h5::detail::OpenGroup root_group(file_id, group_path,
                                   h5::AccessType::ReadOnly);
  const size_t stored_dim =
      h5::read_value_attribute<size_t>(root_group.id(), "Dimension");
  if (stored_dim != Dim) {
    ERROR("Stored dimension "
          << stored_dim << " does not match interpolator dimension " << Dim
          << ".");
  }
  const std::string stored_frame =
      h5::read_value_attribute<std::string>(root_group.id(), "Frame");
  if (stored_frame != pretty_type::name<Frame>()) {
    ERROR("Stored frame '" << stored_frame
                           << "' does not match interpolator frame '"
                           << pretty_type::name<Frame>() << "'.");
  }

  absolute_error_ =
      h5::read_value_attribute<double>(root_group.id(), "AbsoluteError");

  const auto stored_time_bounds =
      h5::read_rank1_attribute<double>(root_group.id(), "TimeBounds");
  time_bounds_[0] = stored_time_bounds[0];
  time_bounds_[1] = stored_time_bounds[1];

  tensor_components_ = h5::read_rank1_attribute<std::string>(
      root_group.id(), "TensorComponents");

  const auto serialized_domain =
      h5::read_data<1, std::vector<char>>(root_group.id(), "Domain");
  domain_ = deserialize<Domain<Dim>>(serialized_domain.data());

  if (h5::contains_dataset_or_group(root_group.id(), "", "FunctionsOfTime")) {
    const auto serialized_fots =
        h5::read_data<1, std::vector<char>>(root_group.id(), "FunctionsOfTime");
    if (serialized_fots.empty()) {
      functions_of_time_.clear();
    } else {
      functions_of_time_ =
          deserialize<domain::FunctionsOfTimeMap>(serialized_fots.data());
    }
  } else {
    functions_of_time_.clear();
  }

  element_search_trees_.clear();
  element_data_.clear();

  if (not h5::contains_dataset_or_group(root_group.id(), "", "Elements")) {
    ERROR("Interpolator file '" << h5_filename
                                << "' does not contain an 'Elements' group.");
  }

  h5::detail::OpenGroup elements_group(root_group.id(), "Elements",
                                       h5::AccessType::ReadOnly);
  const auto element_names = h5::get_group_names(elements_group.id(), "");
  for (const auto& element_name : element_names) {
    h5::detail::OpenGroup element_group(elements_group.id(), element_name,
                                        h5::AccessType::ReadOnly);
    const Index<Dim> extents =
        h5::read_extents<Dim>(element_group.id(), "Extents");
    const auto basis_ints =
        h5::read_rank1_attribute<int>(element_group.id(), "Basis");
    const auto quadrature_ints =
        h5::read_rank1_attribute<int>(element_group.id(), "Quadrature");
    if (basis_ints.size() != Dim or quadrature_ints.size() != Dim) {
      ERROR("Stored basis or quadrature has incorrect size for element "
            << element_name << ".");
    }
    std::array<size_t, Dim> extents_array{};
    std::array<Spectral::Basis, Dim> bases{};
    std::array<Spectral::Quadrature, Dim> quadratures{};
    for (size_t d = 0; d < Dim; ++d) {
      extents_array[d] = extents[d];
      bases[d] = static_cast<Spectral::Basis>(basis_ints[d]);
      quadratures[d] = static_cast<Spectral::Quadrature>(quadrature_ints[d]);
    }
    Mesh<Dim> mesh(extents_array, bases, quadratures);
    ElementData element_data{};
    element_data.mesh = mesh;
    element_data.file_index = 0;
    element_data.component_interpolators.resize(tensor_components_.size());
    for (auto& component : element_data.component_interpolators) {
      component.modal_interpolants.resize(mesh.number_of_grid_points());
    }

    for (size_t component_index = 0;
         component_index < tensor_components_.size(); ++component_index) {
      const std::string component_group_name =
          "Component_" + std::to_string(component_index);
      if (not h5::contains_dataset_or_group(element_group.id(), "",
                                            component_group_name)) {
        ERROR("Missing component group '" << component_group_name
                                          << "' for element " << element_name
                                          << ".");
      }
      h5::detail::OpenGroup component_group(
          element_group.id(), component_group_name, h5::AccessType::ReadOnly);
      const auto offsets = h5::read_data<1, std::vector<size_t>>(
          component_group.id(), "Offsets");
      const size_t num_modes = mesh.number_of_grid_points();
      if (offsets.size() != num_modes + 1) {
        ERROR("Offsets array for element "
              << element_name << " component "
              << tensor_components_[component_index] << " has size "
              << offsets.size() << ", expected " << num_modes + 1 << ".");
      }

      std::vector<double> flat_times{};
      if (h5::contains_dataset_or_group(component_group.id(), "", "Times")) {
        flat_times = h5::read_data<1, std::vector<double>>(component_group.id(),
                                                           "Times");
      }
      std::vector<double> flat_values{};
      if (h5::contains_dataset_or_group(component_group.id(), "", "Values")) {
        flat_values = h5::read_data<1, std::vector<double>>(
            component_group.id(), "Values");
      }
      if (flat_times.size() != flat_values.size()) {
        ERROR("Times and values arrays have different sizes ("
              << flat_times.size() << " vs " << flat_values.size()
              << ") for element " << element_name << ", component "
              << tensor_components_[component_index] << ".");
      }

      const auto has_data =
          h5::read_data<1, std::vector<int>>(component_group.id(), "HasData");
      if (has_data.size() != num_modes) {
        ERROR("HasData array for element "
              << element_name << " component "
              << tensor_components_[component_index] << " has size "
              << has_data.size() << ", expected " << num_modes << ".");
      }

      for (size_t mode_index = 0; mode_index < num_modes; ++mode_index) {
        const size_t begin = offsets[mode_index];
        const size_t end = offsets[mode_index + 1];
        auto& mode_data = element_data.component_interpolators[component_index]
                              .modal_interpolants[mode_index];
        if (has_data[mode_index] == 0 or end <= begin) {
          mode_data.interpolant.reset();
          mode_data.times.clear();
          mode_data.values.clear();
          continue;
        }
        if (end > flat_times.size()) {
          ERROR("Offset entry " << end << " exceeds stored data length "
                                << flat_times.size() << " for element "
                                << element_name << ", component "
                                << tensor_components_[component_index] << ".");
        }
        mode_data.times.assign(flat_times.begin() + begin,
                               flat_times.begin() + end);
        mode_data.values.assign(flat_values.begin() + begin,
                                flat_values.begin() + end);
        mode_data.interpolant =
            boost::math::interpolators::pchip<std::vector<double>>(
                std::vector<double>(mode_data.times),
                std::vector<double>(mode_data.values));
      }
    }

    const ElementId<Dim> element_id(element_name);
    element_search_trees_[element_id.block_id()].insert(element_id);
    element_data_.emplace(element_id, std::move(element_data));
  }
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
      const auto& mesh = meshes.at(element_id);
      const size_t number_of_grid_points = mesh.number_of_grid_points();
      ComponentInterpolator component_interpolator{};
      component_interpolator.modal_interpolants.resize(number_of_grid_points);
      const std::vector<ComponentInterpolator> component_interpolators(
          tensor_components_.size(), component_interpolator);
      ElementData element_data{mesh, file_index, component_interpolators};
      element_data_.emplace(element_id, std::move(element_data));
      element_search_trees_[element_id.block_id()].insert(element_id);
    }
  }
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::build_interpolators(
    const std::vector<std::string>& filenames,
    const std::vector<std::string>& subfile_names) {
  // Build a time interpolant for every tensor component / grid point pair by
  // streaming the observations from disk element-by-element.
  const size_t num_components = tensor_components_.size();
  const ElementId<Dim> reference_id("[B0,(L2I0,L2I3,L2I3)]");
  for (size_t i = 0; i < subfile_names.size(); ++i) {
    const std::string& subfile_name = subfile_names[i];
    const auto obs_ids_and_times =
        load_observation_ids(filenames, subfile_name);
    const double lower_time = obs_ids_and_times.front().second;
    const double upper_time = obs_ids_and_times.back().second;
    if (not std::isfinite(time_bounds_[0]) or lower_time < time_bounds_[0]) {
      time_bounds_[0] = lower_time;
    }
    if (not std::isfinite(time_bounds_[1]) or upper_time > time_bounds_[1]) {
      time_bounds_[1] = upper_time;
    }

    for (const auto& [block_id, tree] : element_search_trees_) {
      if (block_id != 0) {
        continue;
      }
      for (const auto& element_id : tree) {
        // if (element_id != reference_id) {
        //   continue;
        // }
        auto& element_data = element_data_.at(element_id);
        const auto& total_mesh = element_data.mesh;
        auto& component_interpolators = element_data.component_interpolators;

        for (size_t component_index = 0; component_index < num_components;
             ++component_index) {
          auto& component_interpolator =
              component_interpolators[component_index];
          const auto [per_mode_values, extents] = load_component_time_series(
              element_data.file_index, element_id, component_index, filenames,
              subfile_name, obs_ids_and_times);
          const size_t num_grid_points = per_mode_values.size();
          ASSERT(num_grid_points == extents.product(),
                 "Number of grid points does not match extents.");
          for (size_t mode_index = 0; mode_index < num_grid_points;
               ++mode_index) {
            const auto full_index = expanded_index<Dim>(mode_index, extents);
            const auto collapsed_total_index =
                collapsed_index<Dim>(full_index, total_mesh.extents());
            auto& mode_data = component_interpolator.modal_interpolants.at(
                collapsed_total_index);
            if (mode_data.interpolant.has_value()) {
              continue;
            }
            const auto& values = per_mode_values[mode_index];
            auto interpolation_result = interpolant_for_mode(
                obs_ids_and_times, values, absolute_error_);
            const auto full_index_str = get_output(full_index);
            if (interpolation_result.max_error > absolute_error_) {
              Parallel::printf(
                  "For element %s, component %s, mode %s, could not achieve "
                  "the requested relative error tolerance %.3e; achieved "
                  "maximum error %.3e using %zu observation points.\n",
                  get_output(element_id).c_str(),
                  tensor_components_[component_index].c_str(),
                  full_index_str.c_str(), absolute_error_,
                  interpolation_result.max_error,
                  interpolation_result.num_abscissae);
            } else {
              Parallel::printf(
                  "For element %s, component %s, mode %s, achieved "
                  "maximum error %.3e using %zu observation points.\n",
                  get_output(element_id).c_str(),
                  tensor_components_[component_index].c_str(),
                  full_index_str.c_str(), interpolation_result.max_error,
                  interpolation_result.num_abscissae);
            }
            mode_data.times = std::move(interpolation_result.nodes);
            mode_data.values = std::move(interpolation_result.values);
            mode_data.interpolant = std::move(interpolation_result.interpolant);
          }
        }
        Parallel::printf(
            "Constructed interpolator for element %s with %zu tensor "
            "components.\n",
            get_output(element_id).c_str(), num_components);
      }
    }
  }
}

template <size_t Dim, typename Frame>
std::pair<std::vector<std::vector<double>>, Index<Dim>>
ModalSpacetimeInterpolator<Dim, Frame>::load_component_time_series(
    const size_t file_index, const ElementId<Dim>& element_id,
    const size_t component_index, const std::vector<std::string>& filenames,
    const std::string& subfile_name,
    const std::vector<std::pair<size_t, double>>& obs_ids_and_times) const {
  const size_t num_observations = obs_ids_and_times.size();
  // Every row is a modal coefficient, every column corresponds to an
  // observation.

  const std::string element_name = get_output(element_id);
  h5::H5File<h5::AccessType::ReadOnly> element_file(filenames[file_index]);
  const auto& element_volfile = element_file.get<h5::VolumeData>(subfile_name);
  const auto first_obs_id = obs_ids_and_times.front().first;
  const auto all_grid_names_first_obs =
      element_volfile.get_grid_names(first_obs_id);
  const auto all_extents_first_obs = element_volfile.get_extents(first_obs_id);
  const auto all_bases_first_obs = element_volfile.get_bases(first_obs_id);
  const auto all_quadratures_first_obs =
      element_volfile.get_quadratures(first_obs_id);
  const auto mesh_first_obs = h5::mesh_for_grid<Dim>(
      element_name, all_grid_names_first_obs, all_extents_first_obs,
      all_bases_first_obs, all_quadratures_first_obs);
  const size_t num_grid_points = mesh_first_obs.number_of_grid_points();
  std::vector<std::vector<double>> per_mode_values(
      num_grid_points, std::vector<double>(num_observations, 0.0));
  for (size_t obs_index = 0; obs_index < num_observations; ++obs_index) {
    const size_t obs_id = obs_ids_and_times[obs_index].first;
    const auto all_grid_names = element_volfile.get_grid_names(obs_id);
    const auto all_extents = element_volfile.get_extents(obs_id);
    const auto all_bases = element_volfile.get_bases(obs_id);
    const auto all_quadratures = element_volfile.get_quadratures(obs_id);
    const auto mesh = h5::mesh_for_grid<Dim>(
        element_name, all_grid_names, all_extents, all_bases, all_quadratures);
    const auto [offset, length] = h5::offset_and_length_for_grid(
        element_name, all_grid_names, all_extents);
    ASSERT(
        std::find(all_grid_names.begin(), all_grid_names.end(), element_name) !=
            all_grid_names.end(),
        "Element "
            << element_id << " is not present in file index " << file_index
            << " for observation " << obs_id
            << ". Each element is expected to reside in the same volume file "
               "for all observations.");

    if (mesh != mesh_first_obs) {
      ERROR("Element " << element_id
                       << " has inconsistent mesh between observations. AMR is "
                          "not yet supported");
    }
    const auto component_data =
        element_volfile
            .get_tensor_component(obs_id, tensor_components_[component_index])
            .data;
    DataVector modal_data(num_grid_points, 0.0);
    if (std::holds_alternative<DataVector>(component_data)) {
      const auto& data = std::get<DataVector>(component_data);
      const double* element_data =
          data.data() + static_cast<std::ptrdiff_t>(offset);
      for (size_t i = 0; i < num_grid_points; ++i) {
        modal_data[i] = element_data[i];
      }
    } else {
      const auto& data = std::get<std::vector<float>>(component_data);
      for (size_t i = 0; i < num_grid_points; ++i) {
        modal_data[i] = static_cast<double>(data[offset + i]);
      }
    }
    for (size_t mode = 0; mode < num_grid_points; ++mode) {
      per_mode_values[mode][obs_index] = modal_data[mode];
    }
  }
  return std::make_pair(per_mode_values, mesh_first_obs.extents());
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::write_to_h5(
    const std::string& h5_filename, const std::string& group_path) const {
  if (UNLIKELY(element_data_.empty())) {
    ERROR(
        "Cannot write ModalSpacetimeInterpolator to H5 file because no "
        "element data is available.");
  }
  if (UNLIKELY(not std::isfinite(time_bounds_[0]) or
               not std::isfinite(time_bounds_[1]))) {
    ERROR(
        "Cannot write ModalSpacetimeInterpolator to H5 file because the time "
        "bounds are not finite: "
        << time_bounds_ << ".");
  }

  hid_t file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC,
                            h5::h5p_default(), h5::h5p_default());
  CHECK_H5(file_id, "Failed to create H5 file '" << h5_filename << "'.");

  h5::detail::OpenGroup root_group(file_id, group_path,
                                   h5::AccessType::ReadWrite);

  h5::write_to_attribute(root_group.id(), "Dimension",
                         static_cast<size_t>(Dim));
  h5::write_to_attribute(root_group.id(), "Frame", pretty_type::name<Frame>());
  h5::write_to_attribute(root_group.id(), "AbsoluteError", absolute_error_);
  h5::write_to_attribute(root_group.id(), "TimeBounds",
                         std::vector<double>{time_bounds_[0], time_bounds_[1]});
  h5::write_to_attribute(root_group.id(), "TensorComponents",
                         tensor_components_);

  const auto serialized_domain = serialize<Domain<Dim>>(domain_);
  ASSERT(not serialized_domain.empty(),
         "Domain serialization produced an empty buffer.");
  h5::write_data(root_group.id(), serialized_domain, {serialized_domain.size()},
                 "Domain");

  const auto serialized_fots = serialize(functions_of_time_);
  if (not serialized_fots.empty()) {
    h5::write_data(root_group.id(), serialized_fots, {serialized_fots.size()},
                   "FunctionsOfTime");
  }

  h5::detail::OpenGroup elements_group(root_group.id(), "Elements",
                                       h5::AccessType::ReadWrite);
  for (const auto& [element_id, data] : element_data_) {
    const std::string element_name = get_output(element_id);
    h5::detail::OpenGroup element_group(elements_group.id(), element_name,
                                        h5::AccessType::ReadWrite);
    h5::write_extents<Dim>(element_group.id(), data.mesh.extents());
    const auto bases = data.mesh.basis();
    const auto quadratures = data.mesh.quadrature();
    std::vector<int> basis_ints(Dim);
    std::vector<int> quadrature_ints(Dim);
    for (size_t d = 0; d < Dim; ++d) {
      basis_ints[d] = static_cast<int>(bases[d]);
      quadrature_ints[d] = static_cast<int>(quadratures[d]);
    }
    h5::write_to_attribute(element_group.id(), "Basis", basis_ints);
    h5::write_to_attribute(element_group.id(), "Quadrature", quadrature_ints);

    for (size_t component_index = 0;
         component_index < tensor_components_.size(); ++component_index) {
      const std::string component_group_name =
          "Component_" + std::to_string(component_index);
      h5::detail::OpenGroup component_group(
          element_group.id(), component_group_name, h5::AccessType::ReadWrite);
      const auto& component_interpolator =
          data.component_interpolators[component_index];
      const size_t num_modes = component_interpolator.modal_interpolants.size();
      if (num_modes == 0) {
        continue;
      }

      std::vector<size_t> offsets;
      offsets.reserve(num_modes + 1);
      offsets.push_back(0);
      std::vector<double> times_flat{};
      std::vector<double> values_flat{};
      std::vector<int> has_data(num_modes, 0);

      for (size_t mode_index = 0; mode_index < num_modes; ++mode_index) {
        const auto& mode_data =
            component_interpolator.modal_interpolants[mode_index];
        if (mode_data.times.size() != mode_data.values.size()) {
          ERROR("Times and values have different sizes for element "
                << element_name << ", component "
                << tensor_components_[component_index] << ", mode "
                << mode_index << ".");
        }
        if (mode_data.interpolant.has_value() and not mode_data.times.empty()) {
          has_data[mode_index] = 1;
          times_flat.insert(times_flat.end(), mode_data.times.begin(),
                            mode_data.times.end());
          values_flat.insert(values_flat.end(), mode_data.values.begin(),
                             mode_data.values.end());
        }
        offsets.push_back(times_flat.size());
      }

      if (not times_flat.empty()) {
        h5::write_data(component_group.id(), times_flat, {times_flat.size()},
                       "Times");
        h5::write_data(component_group.id(), values_flat, {values_flat.size()},
                       "Values");
      }
      h5::write_data(component_group.id(), offsets, {offsets.size()},
                     "Offsets");
      h5::write_data(component_group.id(), has_data, {has_data.size()},
                     "HasData");
    }
  }
}

template <size_t Dim, typename Frame>
void ModalSpacetimeInterpolator<Dim, Frame>::interpolate_to_point(
    const gsl::not_null<std::vector<double>*> result,
    const tnsr::I<double, Dim, Frame>& target_point, const double time,
    const std::optional<gsl::not_null<std::vector<size_t>*>> block_order)
    const {
  if (UNLIKELY(element_data_.empty())) {
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

  if (not element_data_.contains(element_id)) {
    ERROR("No interpolator data found for element " << element_id << ".");
  }

  const auto& element_data = element_data_.at(element_id);
  const auto& component_interpolators = element_data.component_interpolators;
  const Mesh<Dim>& mesh = element_data.mesh;
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
      const auto& mode_data = component_interpolator.modal_interpolants[point];
      if (UNLIKELY(not mode_data.interpolant.has_value())) {
        ERROR("Missing modal interpolant for component "
              << tensor_components_[component_index] << " in element "
              << element_id << " at grid point " << point << ".");
      }
      modal_values[point] = mode_data.interpolant.value()(time);
      if (point == 0) {
        Parallel::printf(
            "Time interpolation for component %s, element %s, point %zu: value "
            "%.16e\n",
            tensor_components_[component_index].c_str(),
            get_output(element_id).c_str(), point, modal_values[point]);
      }
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

GENERATE_INSTANTIATIONS(INSTANTIATE, (3), (Frame::Inertial))

#undef INSTANTIATE
#undef DIM
#undef FRAME

}  // namespace spectre::Exporter
