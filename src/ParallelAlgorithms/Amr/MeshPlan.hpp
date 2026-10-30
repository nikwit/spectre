// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Info.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "ParallelAlgorithms/Amr/Policies/Policies.hpp"

namespace amr::planning {
/*!
 * Central p-mesh policy planning is enabled under Amr -> Policies with:
 * \code{.yaml}
 * ObjectShellFloors: [ExcisionSphereA, ExcisionSphereB]
 * \endcode
 * Omitting the option (or using []) keeps the existing distributed path.
 * Requires anisotropic policies and full spherical shells at the selected
 * inner boundaries. Schedule EvaluateAmrCriteria/AdjustDomain and use criteria
 * that run in that phase (e.g. SpecSpectral); local RefineMesh changes are
 * rejected while this policy is active. Each object's proposed inner retained
 * cutoff is a floor on its radial shell chain, ending at the shell-to-wedge
 * transition. This also works for a worldtube represented by excision metadata.
 * It does not impose a horizon-finder floor or a maximum difference between
 * neighbors.
 *
 * Planning currently rejects h changes and multigrid. Exact target meshes are
 * delivered together with neighbor targets before projection starts. Policy
 * conflicts fail before projection; maximum-resolution limits are never
 * silently exceeded. The plan is rebuilt from fresh proposals each cycle.
 */

/// Immutable local proposal. No evolution fields are collected by the planner.
/// The topology operation remains in Info::flags; new_mesh is the exact target,
/// including p changes larger than a single mode. The first planner supports
/// fixed topology only and rejects Split/Join before delivering any plan.
template <size_t Dim>
struct Request {
  Element<Dim> element;
  Mesh<Dim> current_mesh;
  Info<Dim> proposed;
  // Number of exactly zeroed top angular bands, not exponential attenuation.
  // nullopt means the filter response is unknown (an error on selected shells).
  std::optional<size_t> filtered_angular_modes;
  void pup(PUP::er& p) {
    p | element;
    p | current_mesh;
    p | proposed;
    p | filtered_angular_modes;
  }
};

template <size_t Dim>
struct Decision {
  Info<Dim> info;
  // Empty when no floor changed the local proposal. Suitable for diagnostics
  // without printing one line per element on every AMR cycle.
  std::optional<std::string> floor_source;
  void pup(PUP::er& p) {
    p | info;
    p | floor_source;
  }
};

template <size_t Dim>
using Requests = std::unordered_map<ElementId<Dim>, Request<Dim>>;
template <size_t Dim>
using Plan = std::unordered_map<ElementId<Dim>, Decision<Dim>>;

/// A pure calculation on a complete snapshot. Starting at each named excision
/// sphere, follow its radial full-S2 shell chain up to the non-S2 transition.
/// Raise proposed retained cutoffs to that of the anchor's proposed mesh.
/// Outer high resolution never raises the anchor. Floors are recomputed each
/// cycle, so they may decrease between cycles. Reject missing/ambiguous chains,
/// conflicting floors, h requests, and targets above hard resolution limits.
Plan<3> make_plan(const Requests<3>& requests, const Domain<3>& domain,
                  const Policies& policies);

/// Convert an exact target mesh into direction flags for existing projectors.
/// These flags describe the sign of a change, not its magnitude.
template <size_t Dim>
Info<Dim> p_decision(const Mesh<Dim>& current, const Mesh<Dim>& target) {
  std::array<Flag, Dim> flags{};
  for (size_t d = 0; d < Dim; ++d) {
    flags[d] = target.extents(d) > current.extents(d) ? Flag::IncreaseResolution
               : target.extents(d) < current.extents(d)
                   ? Flag::DecreaseResolution
                   : Flag::DoNothing;
  }
  return {flags, target};
}
}  // namespace amr::planning
