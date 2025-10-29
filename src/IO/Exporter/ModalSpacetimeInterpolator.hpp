// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <boost/math/interpolators/barycentric_rational.hpp>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/ElementSearchTree.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

namespace spectre::Exporter {

/*!
 * \brief Interpolate tensor components in both space and time using nodal data.
 *
 * \details For each grid point of an element we adaptively construct a
 * barycentric rational interpolant over the observation times. The
 * interpolation adds observation slices until a requested relative error
 * tolerance is satisfied for that grid point. At evaluation we locate the
 * element containing the target point, evaluate the nodal interpolants at the
 * requested time, and interpolate within the element using an
 * `intrp::Irregular` interpolant.
 */
template <size_t Dim, typename Frame = ::Frame::Inertial>
class ModalSpacetimeInterpolator {
 public:
  ModalSpacetimeInterpolator() = default;

  ModalSpacetimeInterpolator(
      std::variant<std::vector<std::string>, std::string> volume_files_or_glob,
      std::string subfile_name, std::vector<std::string> tensor_components,
      double relative_error, size_t max_interpolation_order = 8);

  void interpolate_to_point(gsl::not_null<std::vector<double>*> result,
                            const tnsr::I<double, Dim, Frame>& target_point,
                            double time,
                            std::optional<gsl::not_null<std::vector<size_t>*>>
                                block_order = std::nullopt) const;

 private:
  struct ElementMetadata {
    Mesh<Dim> mesh{};
    size_t offset{};
    size_t length{};
    size_t file_index{};
  };

  struct ComponentInterpolator {
    std::vector<boost::math::interpolators::barycentric_rational<double>>
        modal_interpolants;
  };

  struct ElementInterpolator {
    Mesh<Dim> mesh{};
    std::vector<ComponentInterpolator> component_interpolators{};
  };

  void load_observation_ids(const std::vector<std::string>& filenames);
  void gather_element_metadata(const std::vector<std::string>& filenames);
  void build_interpolators(const std::vector<std::string>& filenames);
  std::vector<std::vector<double>> load_component_time_series(
      const ElementMetadata& metadata, const ElementId<Dim>& element_id,
      size_t component_index, const std::vector<std::string>& filenames) const;

  std::variant<std::vector<std::string>, std::string> volume_files_or_glob_;
  std::string subfile_name_;
  std::vector<std::string> tensor_components_;

  std::vector<std::pair<size_t, double>> obs_ids_and_times_;
  double relative_error_ = 0.0;
  size_t max_interpolation_order_ = 12;
  std::array<double, 2> time_bounds_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};

  Domain<Dim> domain_{};
  domain::FunctionsOfTimeMap functions_of_time_{};
  std::map<size_t, domain::ElementSearchTree<Dim>> element_search_trees_;
  std::unordered_map<ElementId<Dim>, ElementMetadata> element_metadata_;
  std::vector<std::unordered_set<size_t>> file_observation_ids_;
  std::unordered_map<ElementId<Dim>, ElementInterpolator> interpolators_;
};

}  // namespace spectre::Exporter
