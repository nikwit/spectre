// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include <boost/math/interpolators/pchip.hpp>

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
      std::vector<std::string> subfile_names,
      std::vector<std::string> tensor_components, double absolute_error);

  ModalSpacetimeInterpolator(const std::string& h5_filename,
                             const std::string& group_path);

  void write_to_h5(const std::string& h5_filename,
                   const std::string& group_path) const;

  void interpolate_to_point(gsl::not_null<std::vector<double>*> result,
                            const tnsr::I<double, Dim, Frame>& target_point,
                            double time,
                            std::optional<gsl::not_null<std::vector<size_t>*>>
                                block_order = std::nullopt) const;

 private:
  struct ModeInterpolator {
    std::optional<boost::math::interpolators::pchip<std::vector<double>>>
        interpolant{};
    std::vector<double> times{};
    std::vector<double> values{};
  };

  struct ComponentInterpolator {
    std::vector<ModeInterpolator> modal_interpolants{};
  };

  struct ElementData {
    Mesh<Dim> mesh{};
    size_t file_index{};
    std::vector<ComponentInterpolator> component_interpolators{};
  };

  void gather_element_metadata(const std::vector<std::string>& filenames,
                               const std::string& subfile_name,
                               size_t reference_obs_id);
  void build_interpolators(const std::vector<std::string>& filenames,
                           const std::vector<std::string>& subfile_names);
  std::pair<std::vector<std::vector<double>>, Index<Dim>>
  load_component_time_series(
      size_t file_index, const ElementId<Dim>& element_id,
      size_t component_index, const std::vector<std::string>& filenames,
      const std::string& subfile_name,
      const std::vector<std::pair<size_t, double>>& obs_id_and_times) const;

  std::variant<std::vector<std::string>, std::string> volume_files_or_glob_;
  std::vector<std::string> tensor_components_;

  double absolute_error_ = 0.0;
  std::array<double, 2> time_bounds_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};

  Domain<Dim> domain_{};
  domain::FunctionsOfTimeMap functions_of_time_{};
  std::map<size_t, domain::ElementSearchTree<Dim>> element_search_trees_;
  std::unordered_map<ElementId<Dim>, ElementData> element_data_;
};

}  // namespace spectre::Exporter
