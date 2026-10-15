// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"

#include <array>
#include <cstddef>
#include <pup.h>

#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Options/Options.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/StdHelpers.hpp"

namespace intrp::OptionHolders {
namespace {
struct SphereVisitor {
  const Options::Context& context;

  std::set<double> operator()(const double radius) {
    positive_radius(radius);
    return std::set<double>{radius};
  }

  std::set<double> operator()(const std::vector<double>& radii) {
    std::set<double> result;
    for (const double radius : radii) {
      if (result.count(radius) != 0) {
        using ::operator<<;
        PARSE_ERROR(context,
                    "Cannot insert radius "
                        << radius
                        << " into radii for Sphere interpolation target. It "
                           "already exists. Existing radii are "
                        << result);
      }
      positive_radius(radius);
      result.emplace(radius);
    }
    return result;
  }

 private:
  void positive_radius(const double radius) {
    if (radius <= 0) {
      PARSE_ERROR(context, "Radius must be positive, not " << radius);
    }
  }
};
}  // namespace

Sphere::Sphere(const size_t l_max_in, const std::array<double, 3> center_in,
               const typename Radius::type& radius_in,
               const ylm::AngularOrdering angular_ordering_in,
               const Options::Context& context)
    : l_max(l_max_in),
      center(center_in),
      radii(std::visit(SphereVisitor{context}, radius_in)),
      angular_ordering(angular_ordering_in) {}

void Sphere::pup(PUP::er& p) {
  p | l_max;
  p | center;
  p | radii;
  p | angular_ordering;
}

bool operator==(const Sphere& lhs, const Sphere& rhs) {
  return lhs.l_max == rhs.l_max and lhs.center == rhs.center and
         lhs.radii == rhs.radii and
         lhs.angular_ordering == rhs.angular_ordering;
}

bool operator!=(const Sphere& lhs, const Sphere& rhs) {
  return not(lhs == rhs);
}

}  // namespace intrp::OptionHolders

std::array<double, 3> intrp::sphere_detail::excision_center(
    const std::array<double, 3>& grid_center, const Domain<3>& domain,
    const double time, const domain::FunctionsOfTimeMap& functions_of_time) {
  for (const auto& [name, excision] : domain.excision_spheres()) {
    (void)name;
    if (std::array<double, 3>{{get<0>(excision.center()),
                               get<1>(excision.center()),
                               get<2>(excision.center())}} != grid_center) {
      continue;
    }
    if (excision.is_time_dependent()) {
      const auto mapped = excision.moving_mesh_grid_to_inertial_map()(
          excision.center(), time, functions_of_time);
      return {{get<0>(mapped), get<1>(mapped), get<2>(mapped)}};
    }
    return grid_center;
  }
  ERROR(
      "The moving sphere Center must equal a domain excision center in grid "
      "coordinates. Use the same input values or YAML anchors as the domain.");
}
