// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/ElementSearchTree.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/Interpolation/UniformCardinalBSpline.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

namespace spectre::Exporter {

/*!
 * \brief Interpolate tensor components in space and time using modal data.
 *
 * This class constructs a time interpolant for every retained modal
 * coefficient of every tensor component on each element. At evaluation time,
 * it locates the element containing the target point, evaluates the modal
 * time interpolants, and evaluates the resulting Legendre series at the
 * element-logical point.
 *
 * The input consists of one or more volume subfiles in priority order,
 * typically from a coarse spatial mesh observed frequently to a fine spatial
 * mesh observed less frequently. The first subfile containing a mode owns that
 * mode. Its time series is either compressed into an
 * `intrp::UniformCardinalBSpline` or dropped if its amplitude does not exceed
 * the inferred error tolerance. Later subfiles supply modes not present on
 * earlier meshes. Subfile meshes must be nested and the last subfile defines
 * the final spatial mesh.
 *
 * The error tolerance for each element and tensor component is estimated from
 * the lowest mode in the first subfile. This is a per-mode compression
 * heuristic. It does not bound the total reconstructed field error and does
 * not detect temporal undersampling in the original volume data.
 *
 * Construction reads one element's modal time series at a time, so the raw
 * volume data is not loaded into memory as a whole. The compressed
 * interpolants accumulated for all elements remain in memory.
 *
 * The observation times in each subfile must be uniformly spaced, but the
 * subfiles may have different time steps and nonaligned observation times.
 * The available time interval is their intersection. Files from different
 * simulation segments must be joined first, with overlapping observations
 * removed. Changes to the element topology or per-element mesh within one
 * time series, and element migration between files within one time series,
 * are not supported yet. Elements may reside in different files in different
 * subfiles.
 *
 * Moving domains use the global functions of time stored with the final
 * subfile. All functions required by the domain must cover the full available
 * time interval. Only the Legendre basis is supported for now.
 *
 * \note When observing fields on a smaller mesh, use the `ProjectToMesh`
 * option of the ObserveFields event to ensure the data is truncated cleanly
 * in modal space. Do not use `InterpolateToMesh` as this creates a
 * catastrophic aliasing error.
 *
 * \note Requires Boost 1.81 or newer; see
 * `intrp::UniformCardinalBSpline` for details.
 */
template <size_t Dim, typename Frame = ::Frame::Inertial>
class ModalSpacetimeInterpolator {
 public:
  /*!
   * \brief Construct from one or more volume files and nested volume
   * subfiles.
   *
   * \param volume_files_or_glob A list of volume H5 files, or a glob string
   *     that resolves to volume files. The files may distribute the elements
   *     of each observation across nodes.
   * \param subfiles_in_priority_order Volume subfile names ordered from the
   *     preferred source of low modes to the source of the final spatial
   *     mesh. The meshes must have component-wise nondecreasing extents.
   * \param tensor_components Tensor component names to interpolate. The
   *     returned values have this order.
   * \param start_time Optional lower bound used to select observations from
   *     every subfile.
   * \param end_time Optional upper bound used to select observations from
   *     every subfile.
   * \param verbosity Controls diagnostic output during construction.
   */
  ModalSpacetimeInterpolator(
      const std::variant<std::vector<std::string>, std::string>&
          volume_files_or_glob,
      std::vector<std::string> subfiles_in_priority_order,
      std::vector<std::string> tensor_components,
      std::optional<double> start_time = std::nullopt,
      std::optional<double> end_time = std::nullopt,
      Verbosity verbosity = Verbosity::Quiet);

  ModalSpacetimeInterpolator(ModalSpacetimeInterpolator&&) = default;
  ModalSpacetimeInterpolator& operator=(ModalSpacetimeInterpolator&&) = default;
  ModalSpacetimeInterpolator(const ModalSpacetimeInterpolator&) = delete;
  ModalSpacetimeInterpolator& operator=(const ModalSpacetimeInterpolator&) =
      delete;
  ~ModalSpacetimeInterpolator() = default;

  /*!
   * \brief Interpolate the tensor components at a spacetime point.
   *
   * \param result Output buffer, resized to the number of tensor components.
   * \param target_point Coordinates of the target point in `Frame`.
   * \param time Time at which to evaluate the interpolant.
   * \param block_order Optional block-ordering hint that speeds up repeated
   *     block searches.
   */
  void interpolate_to_point(gsl::not_null<std::vector<double>*> result,
                            const tnsr::I<double, Dim, Frame>& target_point,
                            double time,
                            std::optional<gsl::not_null<std::vector<size_t>*>>
                                block_order = std::nullopt) const;

  /// Tensor components in the order returned by `interpolate_to_point()`
  const std::vector<std::string>& tensor_components() const {
    return tensor_components_;
  }

  /// Inclusive time interval shared by all input subfiles
  const std::array<double, 2>& time_bounds() const { return time_bounds_; }

 private:
  using ModeInterpolant = std::optional<intrp::UniformCardinalBSpline>;

  struct ElementData {
    Mesh<Dim> mesh{};
    // Indexed by [component][mode on the final mesh]
    std::vector<std::vector<ModeInterpolant>> interpolants{};
  };

  std::vector<std::string> tensor_components_{};
  std::array<double, 2> time_bounds_{};
  Domain<Dim> domain_{};
  domain::FunctionsOfTimeMap functions_of_time_{};
  std::map<size_t, domain::ElementSearchTree<Dim>> element_search_trees_{};
  std::unordered_map<ElementId<Dim>, ElementData> element_data_{};
};

}  // namespace spectre::Exporter
