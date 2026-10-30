// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Amr/Policies/Policies.hpp"

#include <pup.h>
#include <utility>
#include <yaml-cpp/yaml.h>
#include "Options/ParseOptions.hpp"
#include "Utilities/StdHelpers.hpp"

namespace amr {
Policies::Policies(const amr::Isotropy isotropy, const amr::Limits& limits,
                   const bool enforce_two_to_one_balance_in_normal_direction,
                   const bool allow_coarsening,
                   std::vector<std::string> object_shell_floors)
    : isotropy_(isotropy),
      limits_(limits),
      enforce_two_to_one_balance_in_normal_direction_(
          enforce_two_to_one_balance_in_normal_direction),
      allow_coarsening_(allow_coarsening),
      object_shell_floors_(std::move(object_shell_floors)) {}

void Policies::pup(PUP::er& p) {
  p | isotropy_;
  p | limits_;
  p | enforce_two_to_one_balance_in_normal_direction_;
  p | allow_coarsening_;
  p | object_shell_floors_;
}

std::ostream& operator<<(std::ostream& os, const Policies& policies) {
  using ::operator<<;
  os << "Isotropy: " << policies.isotropy() << "\n";
  os << "EnforceTwoToOneBalanceInNormalDirection: "
     << policies.enforce_two_to_one_balance_in_normal_direction() << "\n";
  os << "AllowCoarsening: " << policies.allow_coarsening() << "\n";
  os << "ObjectShellFloors: " << policies.object_shell_floors() << "\n";
  os << "Limits: " << policies.limits() << "\n";
  return os;
}

bool operator==(const Policies& lhs, const Policies& rhs) {
  return lhs.isotropy() == rhs.isotropy() and lhs.limits() == rhs.limits() and
         lhs.enforce_two_to_one_balance_in_normal_direction() ==
             rhs.enforce_two_to_one_balance_in_normal_direction() and
         lhs.allow_coarsening() == rhs.allow_coarsening() and
         lhs.object_shell_floors() == rhs.object_shell_floors();
}

bool operator!=(const Policies& lhs, const Policies& rhs) {
  return not(lhs == rhs);
}

}  // namespace amr

template <>
amr::Policies Options::create_from_yaml<amr::Policies>::create<void>(
    const Options::Option& option) {
  auto node = YAML::Clone(option.node());
  if (node.IsMap() and not node["ObjectShellFloors"]) {
    node["ObjectShellFloors"] = std::vector<std::string>{};
  }
  Options::Option complete{option.context()};
  complete.set_node(std::move(node));
  Options::Parser<amr::Policies::options> parser{amr::Policies::help};
  parser.parse(complete);
  return parser.apply<amr::Policies::options>(
      [](const auto&... args) { return amr::Policies{args...}; });
}
