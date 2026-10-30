// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Amr/MeshPlan.hpp"

#include <algorithm>
#include <unordered_set>
#include <vector>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Structure/IsValidDgMesh.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace amr::planning {
namespace {
bool is_shell(const Mesh<3>& mesh) {
  return mesh.basis(1) == Spectral::Basis::SphericalHarmonic and
         mesh.basis(2) == Spectral::Basis::SphericalHarmonic;
}
size_t retained_l(const Request<3>& r) {
  if (not r.filtered_angular_modes.has_value() or
      *r.filtered_angular_modes >= r.proposed.new_mesh.extents(1)) {
    ERROR("Unknown or invalid angular filter at " << r.element.id());
  }
  return r.proposed.new_mesh.extents(1) - 1 - *r.filtered_angular_modes;
}
}  // namespace

Plan<3> make_plan(const Requests<3>& requests, const Domain<3>& domain,
                  const Policies& policies) {
  if (requests.empty()) {
    ERROR("Cannot plan an empty AMR snapshot");
  }
  if (not policies.object_shell_floors().empty() and
      policies.isotropy() != Isotropy::Anisotropic) {
    ERROR("ObjectShellFloors currently requires anisotropic AMR policies");
  }
  Plan<3> result;
  for (const auto& [id, r] : requests) {
    if (id != r.element.id() or id.grid_index() != 0) {
      ERROR("Invalid element in AMR snapshot: " << id);
    }
    for (const auto flag : r.proposed.flags) {
      if (flag == Flag::Split or flag == Flag::Join or
          flag == Flag::Undefined) {
        ERROR(
            "Central mesh planning currently requires defined p-only "
            "decisions on a fixed topology: "
            << id);
      }
    }
    if (r.current_mesh.basis() != r.proposed.new_mesh.basis() or
        r.current_mesh.quadrature() != r.proposed.new_mesh.quadrature()) {
      ERROR("AMR planner cannot change basis or quadrature: " << id);
    }
    for (const auto& [direction, neighbors] : r.element.neighbors()) {
      (void)direction;
      for (const auto& neighbor : neighbors.ids()) {
        if (requests.count(neighbor) != 1) {
          ERROR("Missing neighbor " << neighbor << " in AMR snapshot");
        }
      }
    }
    result.emplace(id, Decision<3>{r.proposed, std::nullopt});
  }

  std::unordered_set<std::string> seen_spheres;
  std::unordered_set<ElementId<3>> assigned_shells;
  for (const auto& sphere_name : policies.object_shell_floors()) {
    if (not seen_spheres.insert(sphere_name).second or
        domain.excision_spheres().count(sphere_name) != 1) {
      ERROR("Unknown or repeated ObjectShellFloors sphere: " << sphere_name);
    }
    const auto& sphere = domain.excision_spheres().at(sphere_name);
    std::optional<ElementId<3>> anchor;
    for (const auto& [id, r] : requests) {
      const auto direction = sphere.abutting_direction(id);
      if (direction.has_value()) {
        if (anchor.has_value() or not is_shell(r.proposed.new_mesh) or
            direction->dimension() != 0) {
          ERROR("ObjectShellFloors requires one full-S2 inner element for "
                << sphere_name);
        }
        anchor = id;
      }
    }
    if (not anchor.has_value()) {
      ERROR("Missing inner element for ObjectShellFloors sphere "
            << sphere_name);
    }
    const size_t floor = retained_l(requests.at(*anchor));
    std::optional<ElementId<3>> previous;
    auto current = *anchor;
    while (true) {
      const auto& r = requests.at(current);
      if (not is_shell(r.proposed.new_mesh)) {
        break;
      }
      if (not assigned_shells.insert(current).second) {
        ERROR("Cyclic or overlapping ObjectShellFloors chains at " << current);
      }
      const size_t local_l = retained_l(r);
      if (local_l < floor) {
        auto extents = r.proposed.new_mesh.extents().indices();
        const size_t target_l = floor + *r.filtered_angular_modes;
        if (target_l + 1 > policies.limits().maximum_resolution(
                               r.proposed.new_mesh.basis(1),
                               r.proposed.new_mesh.quadrature(1))) {
          ERROR("ObjectShellFloors exceeds maximum resolution at " << current);
        }
        extents[1] = target_l + 1;
        extents[2] = 2 * target_l + 1;
        const Mesh<3> target{extents, r.proposed.new_mesh.basis(),
                             r.proposed.new_mesh.quadrature()};
        result.at(current) = {p_decision(r.current_mesh, target), sphere_name};
      }
      std::vector<ElementId<3>> next;
      for (const auto& [direction, neighbors] : r.element.neighbors()) {
        if (direction.dimension() != 0) {
          ERROR("Unexpected angular neighbor of full-S2 element " << current);
        }
        for (const auto& neighbor : neighbors.ids()) {
          if (not previous.has_value() or neighbor != *previous) {
            next.push_back(neighbor);
          }
        }
      }
      // A shell-to-wedge transition terminates the chain. A branching chain
      // of full spheres is unsupported and must not silently omit a shell.
      if (next.size() != 1) {
        for (const auto& id : next) {
          if (is_shell(requests.at(id).proposed.new_mesh)) {
            ERROR("Branching ObjectShellFloors chain at " << current);
          }
        }
        break;
      }
      previous = current;
      current = next.front();
    }
  }

  for (const auto& [id, decision] : result) {
    const auto& mesh = decision.info.new_mesh;
    if (not domain::is_valid_dg_mesh(mesh, requests.at(id).element)) {
      ERROR("Invalid planned DG mesh for " << id);
    }
    for (size_t d = 0; d < 3; ++d) {
      if (mesh.extents(d) > policies.limits().maximum_resolution(
                                mesh.basis(d), mesh.quadrature(d))) {
        ERROR("Planned mesh exceeds maximum resolution at "
              << id
              << "; ObjectShellFloors cannot be satisfied within the limits");
      }
    }
  }
  return result;
}
}  // namespace amr::planning
