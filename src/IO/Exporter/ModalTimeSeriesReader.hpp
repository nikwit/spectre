// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "Domain/Structure/ElementId.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace spectre::Exporter {

/*!
 * \brief Reads time series of modal coefficients from volume data files.
 *
 * This class reads tensor components written by multiple observations of a
 * volume data subfile, transforms the nodal data of each element to modal
 * coefficients, and returns their time series. Use
 * `component_modal_time_series()` to read one component for all elements in
 * one file, reading each observation's component dataset only once. The
 * returned buffer holds one component's modal history for that file. The
 * per-element `modal_time_series()` API uses less memory but reads every
 * observation of a tensor component from disk once per element residing in
 * the file.
 *
 * Request the time series with `modal_time_series()`, e.g.:
 *
 * \snippet Test_ModalTimeSeriesReader.cpp modal_time_series_reader_example
 *
 * The observation times must be uniformly spaced. Multiple files can be used
 * when the elements of each observation are distributed across files, such as
 * files written by different nodes. All files must contain the same
 * observations in the requested time interval, and each element must reside
 * in the same file with the same mesh for all observations. Files from
 * different simulation segments must be joined first, with overlapping
 * observations removed. Adaptive mesh refinement and element migration
 * between files (e.g. by load balancing) are not supported yet and raise
 * errors.
 *
 * Only the Legendre basis is supported for now.
 *
 * \note When observing fields on a smaller mesh, use the `ProjectToMesh`
 * option of the ObserveFields event to ensure the data is truncated cleanly
 * in modal space. Do not use `InterpolateToMesh` as this creates a
 * catastrophic aliasing error.
 */
template <size_t Dim>
class ModalTimeSeriesReader {
 public:
  /// One component's modal histories, indexed by mode, then observation.
  /// Modes are ordered by the collapsed index of the element's mesh.
  using ComponentSeries = std::vector<std::vector<double>>;
  /// A single element's histories, indexed by component, mode, observation.
  using Series = std::vector<ComponentSeries>;
  /// One component's modal histories for all elements in a file.
  using FileSeries = std::vector<std::pair<ElementId<Dim>, ComponentSeries>>;

  /*!
   * \brief Construct from one or more volume files.
   *
   * Reads and validates metadata from all files. The volume data itself is
   * only read when calling either time-series accessor.
   *
   * \param volume_files_or_glob A list of volume H5 files, or a glob string
   *     that resolves to volume files. The files can distribute elements of
   *     the same observations across nodes, but files from different
   *     simulation segments must be joined first.
   * \param subfile_name The name of the volume data subfile in the H5 files.
   * \param tensor_components Tensor component names to read. Each component
   *     must exist in every file at every observation.
   * \param start_time Optional lower bound to restrict observations.
   * \param end_time Optional upper bound to restrict observations.
   * \param verbosity Controls metadata and file-cache progress output.
   * \param observation_batch_size Number of observations staged for contiguous
   *     history writes. Uses this many additional component-sized buffers.
   */
  ModalTimeSeriesReader(const std::variant<std::vector<std::string>,
                                           std::string>& volume_files_or_glob,
                        std::string subfile_name,
                        std::vector<std::string> tensor_components,
                        std::optional<double> start_time = std::nullopt,
                        std::optional<double> end_time = std::nullopt,
                        Verbosity verbosity = Verbosity::Silent,
                        size_t observation_batch_size = 16);

  ModalTimeSeriesReader(ModalTimeSeriesReader&&);
  ModalTimeSeriesReader& operator=(ModalTimeSeriesReader&&);
  ModalTimeSeriesReader(const ModalTimeSeriesReader&) = delete;
  ModalTimeSeriesReader& operator=(const ModalTimeSeriesReader&) = delete;
  ~ModalTimeSeriesReader();

  /// Time of the first observation (after filtering by `start_time` and
  /// `end_time`)
  double start_time() const { return obs_ids_and_times_.front().second; }

  /// Uniform spacing between observation times
  double time_step() const { return time_step_; }

  /// Number of observations (after filtering by `start_time` and `end_time`)
  size_t num_observations() const { return obs_ids_and_times_.size(); }

  /// The tensor components that will be read
  const std::vector<std::string>& tensor_components() const {
    return tensor_components_;
  }

  /*!
   * \brief All elements in the volume files and their meshes.
   *
   * The elements are grouped by the volume file they reside in.
   */
  const std::vector<std::pair<ElementId<Dim>, Mesh<Dim>>>& elements() const {
    return elements_;
  }

  /*!
   * \brief The time series of all modal coefficients of all tensor
   * components of the given element.
   *
   * Metadata is cached between calls to either accessor. Reading the volume
   * data is not optimized: each tensor component dataset is read as a whole
   * and only this element's subset is kept. Prefer
   * `component_modal_time_series()` to read all elements in a file.
   */
  Series modal_time_series(const ElementId<Dim>& element_id);

  /// Number of volume files, in the order supplied to the constructor (or
  /// resolved by the glob).
  size_t num_files() const { return filenames_.size(); }

  /*!
   * \brief Read one component's modal histories for all elements in a file.
   *
   * Indices refer to the input files and `tensor_components()`. Each selected
   * observation's entire component dataset is read once, transformed element
   * by element, and scattered into histories indexed by mode and observation.
   * Elements are returned in their order in `elements()`. Mesh metadata and
   * offsets are validated once per file and reused across component passes.
   */
  FileSeries component_modal_time_series(size_t file_index,
                                         size_t component_index);

 private:
  struct FileCache;
  FileCache& file_cache(size_t file_index);

  std::vector<std::string> filenames_;
  std::string subfile_name_;
  std::vector<std::string> tensor_components_;
  std::vector<std::pair<size_t, double>> obs_ids_and_times_;
  double time_step_{};
  Verbosity verbosity_ = Verbosity::Silent;
  size_t observation_batch_size_ = 16;
  /// Elements grouped by file, in the order they appear in the files
  std::vector<std::pair<ElementId<Dim>, Mesh<Dim>>> elements_;
  /// Mesh and file index of each element for fast lookup
  std::unordered_map<ElementId<Dim>, std::pair<Mesh<Dim>, size_t>>
      element_info_;
  std::vector<std::vector<std::pair<ElementId<Dim>, Mesh<Dim>>>>
      elements_by_file_;
  /// Validated offsets for every accessed file, retained across components.
  std::vector<std::unique_ptr<FileCache>> file_caches_;
};

}  // namespace spectre::Exporter
