// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <vector>

#include "Options/Options.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Amr/Policies/Isotropy.hpp"
#include "ParallelAlgorithms/Amr/Policies/Limits.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace amr {
/// \brief A set of runtime policies controlling adaptive mesh refinement
class Policies {
 public:
  /// The isotropy of AMR
  struct Isotropy {
    using type = amr::Isotropy;
    static constexpr Options::String help = {
        "Isotropy of adaptive mesh refinement (whether or not each dimension "
        "can be refined independently)."};
  };

  /// The limits on refinement level and resolution for AMR
  struct Limits {
    using type = amr::Limits;
    static constexpr Options::String help = {
        "Limits on refinement level and resolution for adaptive mesh "
        "refinement."};
  };

  /// Whether or not to enforce 2:1 balance in the direction normal to element
  /// faces
  struct EnforceTwoToOneBalanceInNormalDirection {
    using type = bool;
    static constexpr Options::String help = {
        "Whether or not to enforce 2:1 balance in the direction normal to "
        "element faces."};
  };

  /// Whether or not to allow coarsening or only refining
  struct AllowCoarsening {
    using type = bool;
    static constexpr Options::String help = {
        "Whether or not to allow coarsening or only refining."};
  };

  /// Central p-mesh planning, anchored at the named excision spheres.
  struct ObjectShellFloors {
    using type = std::vector<std::string>;
    static constexpr Options::String help =
        "Excision sphere names whose innermost retained angular resolution "
        "is a floor on their connected spherical shells. For a binary use "
        "[ExcisionSphereA, ExcisionSphereB]. Empty disables central planning. "
        "Currently supports p changes on a fixed topology only.";
  };
  using legacy_options =
      tmpl::list<Isotropy, Limits, EnforceTwoToOneBalanceInNormalDirection,
                 AllowCoarsening>;
  using options = tmpl::push_back<legacy_options, ObjectShellFloors>;

  static constexpr Options::String help = {
      "Policies controlling adaptive mesh refinement."};

  Policies() = default;

  Policies(amr::Isotropy isotropy, const amr::Limits& limits,
           bool enforce_two_to_one_balance_in_normal_direction,
           bool allow_coarsening,
           std::vector<std::string> object_shell_floors = {});

  amr::Isotropy isotropy() const { return isotropy_; }

  amr::Limits limits() const { return limits_; }

  bool enforce_two_to_one_balance_in_normal_direction() const {
    return enforce_two_to_one_balance_in_normal_direction_;
  }

  bool allow_coarsening() const { return allow_coarsening_; }

  const std::vector<std::string>& object_shell_floors() const {
    return object_shell_floors_;
  }

  void pup(PUP::er& p);

 private:
  amr::Isotropy isotropy_{amr::Isotropy::Anisotropic};
  amr::Limits limits_{};
  bool enforce_two_to_one_balance_in_normal_direction_{true};
  bool allow_coarsening_{true};
  std::vector<std::string> object_shell_floors_{};
};

std::ostream& operator<<(std::ostream& os, const Policies& policies);

bool operator==(const Policies& lhs, const Policies& rhs);

bool operator!=(const Policies& lhs, const Policies& rhs);
}  // namespace amr

// Older inputs omit ObjectShellFloors. Preserve that syntax with an empty list.
template <>
struct Options::create_from_yaml<amr::Policies> {
  template <typename Metavariables>
  static amr::Policies create(const Options::Option& option) {
    return create<void>(option);
  }
};
template <>
amr::Policies Options::create_from_yaml<amr::Policies>::create<void>(
    const Options::Option& option);
