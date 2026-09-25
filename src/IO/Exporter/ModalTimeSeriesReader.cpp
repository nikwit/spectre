// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "IO/Exporter/ModalTimeSeriesReader.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Overloader.hpp"

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
    ERROR_NO_TRACE("No volume files found. Specify at least one volume file.");
  }
  return filenames;
}

template <size_t Dim>
void enforce_dimension(const h5::VolumeData& volfile,
                       const std::string& filename) {
  if (volfile.get_dimension() != Dim) {
    ERROR_NO_TRACE("Mismatched dimensions: expected "
                   << Dim << "D volume data, but got "
                   << volfile.get_dimension() << "D in file '" << filename
                   << "'.");
  }
}

std::vector<std::pair<size_t, double>> observation_ids_and_times(
    const h5::VolumeData& volfile, const std::optional<double> start_time,
    const std::optional<double> end_time, const Verbosity verbosity) {
  const auto start = std::chrono::steady_clock::now();
  if (verbosity >= Verbosity::Verbose) {
    Parallel::printf(
        "Listing and sorting all observations before filtering "
        "by the requested time interval...\n");
  }
  const auto observation_ids = volfile.list_observation_ids();
  std::vector<std::pair<size_t, double>> observations{};
  observations.reserve(observation_ids.size());
  for (const size_t observation_id : observation_ids) {
    observations.emplace_back(observation_id,
                              volfile.get_observation_value(observation_id));
  }
  if (verbosity >= Verbosity::Verbose) {
    Parallel::printf(
        "Listed %zu observations in %.2f s; filtering observation times "
        "in memory...\n",
        observations.size(),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count());
  }
  std::vector<std::pair<size_t, double>> result{};
  result.reserve(observations.size());
  for (const auto& [observation_id, observation_time] : observations) {
    if ((start_time.has_value() and observation_time < start_time.value()) or
        (end_time.has_value() and observation_time > end_time.value())) {
      continue;
    }
    result.emplace_back(observation_id, observation_time);
  }
  if (verbosity >= Verbosity::Verbose) {
    Parallel::printf(
        "Selected %zu of %zu observations; metadata scan took %.2f s.\n",
        result.size(), observations.size(),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count());
  }
  return result;
}

template <size_t Dim>
void enforce_legendre_basis(const Mesh<Dim>& mesh, const std::string& context) {
  for (size_t d = 0; d < Dim; ++d) {
    if (gsl::at(mesh.basis(), d) != Spectral::Basis::Legendre) {
      ERROR_NO_TRACE("Only the Legendre basis is supported, but found "
                     << gsl::at(mesh.basis(), d) << " in dimension " << d
                     << " for " << context << ".");
    }
  }
}

// Data of one observation in one file that is needed to locate and interpret
// each element's subset of the tensor component datasets
struct ObservationCache {
  std::unordered_map<std::string, size_t> grid_index_by_name{};
  // Offset of each grid into the contiguous datasets, plus a total size at
  // the end
  std::vector<size_t> grid_offsets{};
  std::vector<std::vector<size_t>> extents{};
  std::vector<std::vector<Spectral::Basis>> bases{};
  std::vector<std::vector<Spectral::Quadrature>> quadratures{};
};

ObservationCache make_observation_cache(const h5::VolumeData& volfile,
                                        const size_t obs_id,
                                        const std::string& filename) {
  ObservationCache cache{};
  const auto grid_names = volfile.get_grid_names(obs_id);
  cache.extents = volfile.get_extents(obs_id);
  cache.bases = volfile.get_bases(obs_id);
  cache.quadratures = volfile.get_quadratures(obs_id);
  cache.grid_index_by_name.reserve(grid_names.size());
  cache.grid_offsets.reserve(grid_names.size() + 1);
  cache.grid_offsets.push_back(0);
  for (size_t grid_index = 0; grid_index < grid_names.size(); ++grid_index) {
    if (not cache.grid_index_by_name.emplace(grid_names[grid_index], grid_index)
                .second) {
      ERROR_NO_TRACE("Grid '" << grid_names[grid_index]
                              << "' appears multiple times at observation ID "
                              << obs_id << " in file '" << filename << "'.");
    }
    size_t num_points = 1;
    for (const size_t extent : cache.extents[grid_index]) {
      num_points *= extent;
    }
    cache.grid_offsets.push_back(cache.grid_offsets.back() + num_points);
  }
  return cache;
}

// Read a whole tensor component dataset and check that it holds exactly the
// grid points of the file's elements, so that copying an element's subset at
// its cached offset stays in bounds.
TensorComponent read_tensor_component(const h5::VolumeData& volfile,
                                      const size_t obs_id,
                                      const std::string& component,
                                      const size_t expected_size,
                                      const std::string& filename) {
  auto tensor_component = volfile.get_tensor_component(obs_id, component);
  const size_t size = std::visit([](const auto& data) { return data.size(); },
                                 tensor_component.data);
  if (size != expected_size) {
    ERROR_NO_TRACE("Tensor component '"
                   << component << "' at observation ID " << obs_id
                   << " in file '" << filename << "' has " << size
                   << " values, but the elements in this file have "
                   << expected_size << " grid points in total.");
  }
  return tensor_component;
}

template <size_t Dim>
Mesh<Dim> mesh_from_cache(const ObservationCache& cache,
                          const size_t grid_index) {
  std::array<size_t, Dim> extents{};
  std::array<Spectral::Basis, Dim> bases{};
  std::array<Spectral::Quadrature, Dim> quadratures{};
  for (size_t d = 0; d < Dim; ++d) {
    gsl::at(extents, d) = cache.extents[grid_index][d];
    gsl::at(bases, d) = cache.bases[grid_index][d];
    gsl::at(quadratures, d) = cache.quadratures[grid_index][d];
  }
  return Mesh<Dim>{extents, bases, quadratures};
}

}  // namespace

// Keep only validated offsets, rather than copies of every observation's
// element names and meshes. These remain cached across all component passes.
template <size_t Dim>
struct ModalTimeSeriesReader<Dim>::FileCache {
  FileCache(const std::string& filename, const std::string& subfile_name,
            const std::vector<std::pair<size_t, double>>& obs_ids_and_times,
            const std::vector<std::pair<ElementId<Dim>, Mesh<Dim>>>& elements,
            const Verbosity verbosity)
      : h5file(filename), volfile(&h5file.get<h5::VolumeData>(subfile_name)) {
    std::vector<std::string> element_names{};
    for (size_t i = 0; i < elements.size(); ++i) {
      element_index_by_id.emplace(elements[i].first, i);
      element_names.push_back(get_output(elements[i].first));
      num_points += elements[i].second.number_of_grid_points();
    }
    offsets.reserve(obs_ids_and_times.size());
    const auto start = std::chrono::steady_clock::now();
    for (const auto& [obs_id, obs_time] : obs_ids_and_times) {
      const auto observation =
          make_observation_cache(*volfile, obs_id, filename);
      if (observation.grid_index_by_name.size() != elements.size()) {
        ERROR_NO_TRACE("File '"
                       << filename << "' contains "
                       << observation.grid_index_by_name.size()
                       << " elements at observation ID " << obs_id << " but "
                       << elements.size()
                       << " elements at the last observation. Adaptive mesh "
                          "refinement and element migration between files "
                          "(e.g. by load balancing) are not supported.");
      }
      std::vector<size_t> observation_offsets{};
      observation_offsets.reserve(elements.size());
      for (size_t i = 0; i < elements.size(); ++i) {
        const auto grid_it =
            observation.grid_index_by_name.find(element_names[i]);
        if (grid_it == observation.grid_index_by_name.end()) {
          ERROR_NO_TRACE("Element "
                         << element_names[i] << " is missing in file '"
                         << filename << "' at observation ID " << obs_id
                         << ". Each element must reside in the same volume "
                            "file for all observations, so element migration "
                            "between files (e.g. by load balancing) is not "
                            "supported.");
        }
        const auto mesh = mesh_from_cache<Dim>(observation, grid_it->second);
        if (mesh != elements[i].second) {
          ERROR_NO_TRACE("Element "
                         << element_names[i] << " in file '" << filename
                         << "' has mesh " << mesh << " at observation ID "
                         << obs_id << " but mesh " << elements[i].second
                         << " at the last observation. Mesh changes between "
                            "observations (e.g. by adaptive mesh refinement) "
                            "are not supported.");
        }
        observation_offsets.push_back(
            observation.grid_offsets[grid_it->second]);
      }
      offsets.push_back(std::move(observation_offsets));
      if (verbosity >= Verbosity::Verbose and
          (offsets.size() % 100 == 0 or
           offsets.size() == obs_ids_and_times.size())) {
        Parallel::printf(
            "Cached mesh metadata for %zu/%zu observations in %.2f s.\n",
            offsets.size(), obs_ids_and_times.size(),
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          start)
                .count());
      }
    }
  }

  h5::H5File<h5::AccessType::ReadOnly> h5file;
  const h5::VolumeData* volfile;
  size_t num_points = 0;
  std::unordered_map<ElementId<Dim>, size_t> element_index_by_id{};
  // Indexed by observation, then by the canonical element order for this file.
  std::vector<std::vector<size_t>> offsets{};
};

template <size_t Dim>
ModalTimeSeriesReader<Dim>::ModalTimeSeriesReader(
    const std::variant<std::vector<std::string>, std::string>&
        volume_files_or_glob,
    std::string subfile_name, std::vector<std::string> tensor_components,
    const std::optional<double> start_time,
    const std::optional<double> end_time, const Verbosity verbosity,
    const size_t observation_batch_size)
    : filenames_(resolve_filenames(volume_files_or_glob)),
      subfile_name_(std::move(subfile_name)),
      tensor_components_(std::move(tensor_components)),
      verbosity_(verbosity),
      observation_batch_size_(observation_batch_size) {
  if (observation_batch_size_ == 0) {
    ERROR_NO_TRACE("The observation batch size must be positive.");
  }
  if (verbosity_ >= Verbosity::Verbose) {
    Parallel::printf(
        "Reader settings for subfile %s: observation batch size %zu.\n",
        subfile_name_.c_str(), observation_batch_size_);
  }
  // Load the observation IDs and times from the first file, restrict them to
  // the requested time interval, and check they are uniformly spaced
  {
    if (verbosity_ >= Verbosity::Verbose) {
      Parallel::printf(
          "Reading observation metadata for subfile %s, file %s "
          "(1/%zu)...\n",
          subfile_name_.c_str(), filenames_.front().c_str(), filenames_.size());
    }
    const h5::H5File<h5::AccessType::ReadOnly> first_h5file(filenames_.front());
    const auto& volfile = first_h5file.get<h5::VolumeData>(subfile_name_);
    enforce_dimension<Dim>(volfile, filenames_.front());
    obs_ids_and_times_ =
        observation_ids_and_times(volfile, start_time, end_time, verbosity_);
  }
  if (obs_ids_and_times_.size() < 2) {
    ERROR_NO_TRACE("At least 2 observations are required, but found "
                   << obs_ids_and_times_.size() << " in subfile '"
                   << subfile_name_ << "'"
                   << (start_time.has_value() or end_time.has_value()
                           ? " after restricting to the requested time "
                             "interval."
                           : "."));
  }
  const double first_time = obs_ids_and_times_.front().second;
  time_step_ = obs_ids_and_times_[1].second - first_time;
  const double relative_epsilon =
      100.0 * std::numeric_limits<double>::epsilon();
  if (time_step_ <= 0.0) {
    ERROR_NO_TRACE("Observation times in subfile '"
                   << subfile_name_
                   << "' must be strictly increasing, but the first time step "
                      "is "
                   << time_step_
                   << ". Select a strictly increasing interval with "
                      "start_time and end_time or preprocess the volume "
                      "files.");
  }
  for (size_t i = 2; i < obs_ids_and_times_.size(); ++i) {
    const double actual_time = obs_ids_and_times_[i].second;
    const double expected_time =
        first_time + static_cast<double>(i) * time_step_;
    const double scale =
        std::max({1.0, std::abs(expected_time), std::abs(actual_time)});
    if (not equal_within_roundoff(actual_time, expected_time, relative_epsilon,
                                  scale)) {
      ERROR_NO_TRACE(
          "Observation times in subfile '"
          << subfile_name_
          << "' must be uniformly spaced, but observation index " << i
          << " has time " << actual_time << " where " << expected_time
          << " was expected from the first time " << first_time
          << " and the time step " << time_step_
          << ". Select a uniformly spaced interval with start_time and "
             "end_time or preprocess the volume files.");
    }
  }

  // Check that all other files contain the same observations
  for (size_t file_index = 1; file_index < filenames_.size(); ++file_index) {
    if (verbosity_ >= Verbosity::Verbose) {
      Parallel::printf(
          "Reading observation metadata for subfile %s, file %s "
          "(%zu/%zu)...\n",
          subfile_name_.c_str(), filenames_[file_index].c_str(), file_index + 1,
          filenames_.size());
    }
    const h5::H5File<h5::AccessType::ReadOnly> h5file(filenames_[file_index]);
    const auto& volfile = h5file.get<h5::VolumeData>(subfile_name_);
    enforce_dimension<Dim>(volfile, filenames_[file_index]);
    const auto other_obs_ids_and_times =
        observation_ids_and_times(volfile, start_time, end_time, verbosity_);
    if (other_obs_ids_and_times.size() != obs_ids_and_times_.size()) {
      ERROR_NO_TRACE("File '"
                     << filenames_[file_index] << "' contains "
                     << other_obs_ids_and_times.size()
                     << " observations in the requested time interval, but "
                     << "file '" << filenames_.front() << "' contains "
                     << obs_ids_and_times_.size()
                     << ". All volume files must contain the same "
                        "observations in the requested time interval.");
    }
    for (size_t obs_index = 0; obs_index < obs_ids_and_times_.size();
         ++obs_index) {
      const auto& [obs_id, obs_time] = obs_ids_and_times_[obs_index];
      const auto& [other_obs_id, other_time] =
          other_obs_ids_and_times[obs_index];
      if (other_obs_id != obs_id) {
        ERROR_NO_TRACE("Mismatched observation ID at observation index "
                       << obs_index << ": expected " << obs_id << " from file '"
                       << filenames_.front() << "' but found " << other_obs_id
                       << " in file '" << filenames_[file_index]
                       << "'. All volume files must contain the same "
                          "observations in the requested time interval.");
      }
      if (other_time != obs_time) {
        ERROR_NO_TRACE("Mismatched observation value for observation ID "
                       << obs_id << ": expected " << obs_time << " from file '"
                       << filenames_.front() << "' but found " << other_time
                       << " in file '" << filenames_[file_index] << "'.");
      }
    }
  }

  // Gather the elements and their meshes from all files at the last
  // observation. Consistency with the other observations is checked when
  // populating the per-file metadata cache.
  elements_by_file_.resize(filenames_.size());
  file_caches_.resize(filenames_.size());
  const size_t reference_obs_id = obs_ids_and_times_.back().first;
  for (size_t file_index = 0; file_index < filenames_.size(); ++file_index) {
    if (verbosity_ >= Verbosity::Verbose) {
      Parallel::printf("Reading element meshes for subfile %s, file %s...\n",
                       subfile_name_.c_str(), filenames_[file_index].c_str());
    }
    const h5::H5File<h5::AccessType::ReadOnly> h5file(filenames_[file_index]);
    const auto& volfile = h5file.get<h5::VolumeData>(subfile_name_);
    const auto grid_names = volfile.get_grid_names(reference_obs_id);
    const auto all_extents = volfile.get_extents(reference_obs_id);
    const auto all_bases = volfile.get_bases(reference_obs_id);
    const auto all_quadratures = volfile.get_quadratures(reference_obs_id);
    for (const auto& grid_name : grid_names) {
      const ElementId<Dim> element_id(grid_name);
      const auto mesh = h5::mesh_for_grid<Dim>(
          grid_name, grid_names, all_extents, all_bases, all_quadratures);
      enforce_legendre_basis(mesh, "element " + grid_name + " in file '" +
                                       filenames_[file_index] + "'");
      if (const auto existing_it = element_info_.find(element_id);
          existing_it != element_info_.end()) {
        ERROR_NO_TRACE("Element "
                       << grid_name << " in file '" << filenames_[file_index]
                       << "' already exists in file '"
                       << filenames_[existing_it->second.second]
                       << "'. Each element must reside in exactly one volume "
                          "file.");
      }
      element_info_.emplace(element_id, std::make_pair(mesh, file_index));
      elements_.emplace_back(element_id, mesh);
      elements_by_file_[file_index].emplace_back(element_id, mesh);
    }
  }
}

template <size_t Dim>
ModalTimeSeriesReader<Dim>::ModalTimeSeriesReader(ModalTimeSeriesReader&&) =
    default;
template <size_t Dim>
ModalTimeSeriesReader<Dim>& ModalTimeSeriesReader<Dim>::operator=(
    ModalTimeSeriesReader&&) = default;
template <size_t Dim>
ModalTimeSeriesReader<Dim>::~ModalTimeSeriesReader() = default;

template <size_t Dim>
typename ModalTimeSeriesReader<Dim>::Series
ModalTimeSeriesReader<Dim>::modal_time_series(
    const ElementId<Dim>& element_id) {
  const auto info_it = element_info_.find(element_id);
  if (info_it == element_info_.end()) {
    ERROR_NO_TRACE("Element "
                   << element_id
                   << " does not exist in the volume files. The available "
                      "elements are listed by `elements()`.");
  }
  const auto& [mesh, file_index] = info_it->second;
  auto& cache = file_cache(file_index);
  const auto& volfile = *cache.volfile;
  const size_t element_index = cache.element_index_by_id.at(element_id);
  const size_t num_observations = obs_ids_and_times_.size();
  const size_t num_components = tensor_components_.size();
  const size_t num_points = mesh.number_of_grid_points();
  Series series(num_components,
                std::vector<std::vector<double>>(
                    num_points, std::vector<double>(num_observations)));
  DataVector nodal_data(num_points);
  ModalVector modal_data(num_points);
  for (size_t obs_index = 0; obs_index < num_observations; ++obs_index) {
    const size_t obs_id = obs_ids_and_times_[obs_index].first;
    const size_t offset = cache.offsets[obs_index][element_index];
    for (size_t component_index = 0; component_index < num_components;
         ++component_index) {
      const auto tensor_component = read_tensor_component(
          volfile, obs_id, tensor_components_[component_index],
          cache.num_points, filenames_[file_index]);
      std::visit(
          [&nodal_data, offset, num_points](const auto& data) {
            std::copy_n(data.data() + offset, num_points, nodal_data.begin());
          },
          tensor_component.data);
      to_modal_coefficients(make_not_null(&modal_data), nodal_data, mesh);
      auto& component_series = series[component_index];
      for (size_t mode = 0; mode < num_points; ++mode) {
        component_series[mode][obs_index] = modal_data[mode];
      }
    }
  }
  return series;
}

template <size_t Dim>
typename ModalTimeSeriesReader<Dim>::FileCache&
ModalTimeSeriesReader<Dim>::file_cache(const size_t file_index) {
  if (file_index >= filenames_.size()) {
    ERROR_NO_TRACE("Volume file index " << file_index << " is out of range for "
                                        << filenames_.size() << " files.");
  }
  auto& cache = file_caches_[file_index];
  if (cache == nullptr) {
    if (verbosity_ >= Verbosity::Verbose) {
      Parallel::printf(
          "Caching mesh metadata for %zu observations in subfile "
          "%s, file %s...\n",
          obs_ids_and_times_.size(), subfile_name_.c_str(),
          filenames_[file_index].c_str());
    }
    cache = std::make_unique<FileCache>(
        filenames_[file_index], subfile_name_, obs_ids_and_times_,
        elements_by_file_[file_index], verbosity_);
  }
  return *cache;
}

template <size_t Dim>
typename ModalTimeSeriesReader<Dim>::FileSeries
ModalTimeSeriesReader<Dim>::component_modal_time_series(
    const size_t file_index, const size_t component_index) {
  if (component_index >= tensor_components_.size()) {
    ERROR_NO_TRACE("Tensor component index "
                   << component_index << " is out of range for "
                   << tensor_components_.size() << " components.");
  }
  auto& cache = file_cache(file_index);
  const auto& elements = elements_by_file_[file_index];
  const auto& component = tensor_components_[component_index];
  const size_t num_observations = obs_ids_and_times_.size();
  if (verbosity_ >= Verbosity::Verbose) {
    Parallel::printf(
        "Reading component %s (%zu/%zu), subfile %s, file %s (%zu/%zu): "
        "%zu elements, %zu observations, %.3f GiB modal buffer.\n",
        component.c_str(), component_index + 1, tensor_components_.size(),
        subfile_name_.c_str(), filenames_[file_index].c_str(), file_index + 1,
        filenames_.size(), elements.size(), num_observations,
        static_cast<double>(cache.num_points) *
            static_cast<double>(num_observations) * sizeof(double) /
            (1024.0 * 1024.0 * 1024.0));
  }
  const auto start = std::chrono::steady_clock::now();
  FileSeries result{};
  result.reserve(elements.size());
  for (const auto& [element_id, mesh] : elements) {
    result.emplace_back(element_id,
                        ComponentSeries(mesh.number_of_grid_points(),
                                        std::vector<double>(num_observations)));
  }
  if (elements.empty()) {
    return result;
  }
  const size_t batch_size = std::min(observation_batch_size_, num_observations);
  // The staging buffer is observation-major. Each transformed element is
  // copied contiguously into it; a batch is then transposed into short,
  // contiguous stretches of the mode histories.
  DataVector staged_modes(batch_size * cache.num_points);
  DataVector nodal_data{};
  ModalVector modal_data{};
  using Clock = std::chrono::steady_clock;
  const auto elapsed = [](const auto since) {
    return std::chrono::duration<double>(Clock::now() - since).count();
  };
  const double allocation_seconds = elapsed(start);
  double read_seconds = 0.0;
  double transform_seconds = 0.0;
  double history_seconds = 0.0;
  if (verbosity_ >= Verbosity::Verbose) {
    Parallel::printf(
        "Allocated histories in %.2f s; staging %zu observations "
        "(%.3f GiB).\n",
        allocation_seconds, batch_size,
        static_cast<double>(staged_modes.size()) * sizeof(double) /
            (1024.0 * 1024.0 * 1024.0));
  }
  for (size_t batch_start = 0; batch_start < num_observations;
       batch_start += batch_size) {
    const size_t batch_count =
        std::min(batch_size, num_observations - batch_start);
    for (size_t j = 0; j < batch_count; ++j) {
      const size_t obs_index = batch_start + j;
      const auto read_start = Clock::now();
      // Read the entire concatenated component in one H5 read. This timer
      // includes HDF5 metadata access, allocation, and decompression.
      const auto tensor_component = read_tensor_component(
          *cache.volfile, obs_ids_and_times_[obs_index].first, component,
          cache.num_points, filenames_[file_index]);
      read_seconds += elapsed(read_start);
      const auto transform_start = Clock::now();
      size_t mode_offset = j * cache.num_points;
      for (size_t element_index = 0; element_index < elements.size();
           ++element_index) {
        const auto& mesh = elements[element_index].second;
        const size_t num_points = mesh.number_of_grid_points();
        const size_t offset = cache.offsets[obs_index][element_index];
        nodal_data.destructive_resize(num_points);
        modal_data.destructive_resize(num_points);
        std::visit(
            [&nodal_data, offset, num_points](const auto& data) {
              std::copy_n(data.data() + offset, num_points, nodal_data.begin());
            },
            tensor_component.data);
        to_modal_coefficients(make_not_null(&modal_data), nodal_data, mesh);
        std::copy_n(modal_data.data(), num_points,
                    staged_modes.data() + mode_offset);
        mode_offset += num_points;
      }
      // Includes nodal copies and contiguous writes to the staging buffer.
      transform_seconds += elapsed(transform_start);
    }
    const auto history_start = Clock::now();
    size_t mode_offset = 0;
    for (auto& [element_id, series] : result) {
      for (auto& history : series) {
        for (size_t j = 0; j < batch_count; ++j) {
          history[batch_start + j] =
              staged_modes[j * cache.num_points + mode_offset];
        }
        ++mode_offset;
      }
    }
    history_seconds += elapsed(history_start);
    const size_t completed = batch_start + batch_count;
    if (verbosity_ >= Verbosity::Verbose and
        (batch_start == 0 or completed / 100 != batch_start / 100 or
         completed == num_observations)) {
      Parallel::printf(
          "Read and transformed %zu/%zu observations for component %s in "
          "%.2f s (HDF5 %.2f s, transforms/staging %.2f s, history copies "
          "%.2f s).\n",
          completed, num_observations, component.c_str(), elapsed(start),
          read_seconds, transform_seconds, history_seconds);
    }
  }
  if (verbosity_ >= Verbosity::Verbose) {
    Parallel::printf(
        "Reader timings for component %s, subfile %s, file %zu/%zu: "
        "allocation %.2f s, HDF5 %.2f s, transforms/staging %.2f s, "
        "history copies %.2f s, total %.2f s.\n",
        component.c_str(), subfile_name_.c_str(), file_index + 1,
        filenames_.size(), allocation_seconds, read_seconds, transform_seconds,
        history_seconds, elapsed(start));
  }
  return result;
}

// Explicit instantiations

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) template class ModalTimeSeriesReader<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM

}  // namespace spectre::Exporter
