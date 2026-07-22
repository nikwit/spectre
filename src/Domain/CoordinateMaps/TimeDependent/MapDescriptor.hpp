// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace domain::CoordinateMaps::TimeDependent {

/// Host-side metadata for a RotScaleTrans coordinate map.
struct RotScaleTransMapDescriptor {
  enum class Region { Inner, Transition, Outer };

  std::optional<std::pair<std::string, std::string>> scale_names{};
  std::optional<std::string> rotation_name{};
  std::optional<std::string> translation_name{};
  double inner_radius{};
  double outer_radius{};
  Region region{Region::Inner};
};

/// Host-side metadata for a Skew coordinate map.
struct SkewMapDescriptor {
  std::string function_of_time_name{};
  std::array<double, 3> center{};
  double outer_radius{};
};

/// Host-side metadata for the spherical Shape-map transition.
struct SphereTransitionDescriptor {
  double inner_radius{};
  double outer_radius{};
  double radial_slope{};
  double radial_intercept{};
  bool interior{};
};

/// Host-side metadata for the wedge Shape-map transition.
struct WedgeTransitionDescriptor {
  std::array<double, 3> inner_center{};
  double inner_radius{};
  double inner_sphericity{};
  std::array<double, 3> outer_center{};
  double outer_radius{};
  double outer_sphericity{};
  int axis{};
  bool reverse{};
};

using ShapeMapTransitionDescriptor =
    std::variant<SphereTransitionDescriptor, WedgeTransitionDescriptor>;

/// Host-side metadata for a Shape coordinate map.
struct ShapeMapDescriptor {
  std::array<double, 3> center{};
  double truncation_limit{};
  ShapeMapTransitionDescriptor transition{};
  std::string shape_function_of_time_name{};
  std::optional<std::string> size_function_of_time_name{};
};

using MapDescriptor = std::variant<RotScaleTransMapDescriptor,
                                   SkewMapDescriptor, ShapeMapDescriptor>;

}  // namespace domain::CoordinateMaps::TimeDependent
