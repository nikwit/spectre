// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <unordered_set>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Block.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/Topology.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "ParallelAlgorithms/Amr/MeshPlan.hpp"
#include "Utilities/MakeArray.hpp"

namespace {
ElementId<3> id(const size_t object, const size_t shell) {
  return ElementId<3>{object, std::array{SegmentId{2, shell}, SegmentId{0, 0},
                                         SegmentId{0, 0}}};
}
Mesh<3> mesh(const size_t l, const size_t nr = 9) {
  return Mesh<3>{
      std::array{nr, l + 1, 2 * l + 1},
      std::array{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
                 Spectral::Basis::SphericalHarmonic},
      std::array{Spectral::Quadrature::GaussLobatto,
                 Spectral::Quadrature::Gauss,
                 Spectral::Quadrature::Equiangular}};
}
Domain<3> make_domain() {
  // The planner needs topology and excision metadata, not coordinate maps.
  std::unordered_map<std::string, ExcisionSphere<3>> spheres;
  for (size_t object = 0; object < 2; ++object) {
    spheres.emplace(object == 0 ? "ExcisionSphereA" : "ExcisionSphereB",
                    ExcisionSphere<3>{1.,
                                      tnsr::I<double, 3, Frame::Grid>{0.},
                                      {{object, Direction<3>::lower_xi()}}});
  }
  return Domain<3>{std::vector<Block<3>>{}, std::move(spheres)};
}
amr::planning::Requests<3> requests() {
  amr::planning::Requests<3> result;
  for (size_t object = 0; object < 2; ++object) {
    for (size_t shell = 0; shell < 4; ++shell) {
      DirectionMap<3, Neighbors<3>> ns;
      if (shell > 0) {
        ns.emplace(Direction<3>::lower_xi(),
                   Neighbors<3>{std::unordered_set{id(object, shell - 1)},
                                OrientationMap<3>::create_aligned()});
      }
      if (shell < 3) {
        ns.emplace(Direction<3>::upper_xi(),
                   Neighbors<3>{std::unordered_set{id(object, shell + 1)},
                                OrientationMap<3>::create_aligned()});
      }
      const auto local_mesh = mesh(shell == 0 ? (object == 0 ? 15 : 12) : 14);
      result.emplace(id(object, shell),
                     amr::planning::Request<3>{
                         Element<3>{id(object, shell), std::move(ns),
                                    std::array{domain::Topology::I1,
                                               domain::Topology::S2Colatitude,
                                               domain::Topology::S2Longitude}},
                         local_mesh,
                         amr::planning::p_decision(local_mesh, local_mesh), 4});
    }
  }
  return result;
}
amr::Policies policies(const std::vector<std::string>& names = {
                           "ExcisionSphereA", "ExcisionSphereB"}) {
  return {amr::Isotropy::Anisotropic, amr::Limits{}, false, true, names};
}
void propose(amr::planning::Request<3>& r, const size_t l,
             const size_t nr = 9) {
  r.proposed = amr::planning::p_decision(r.current_mesh, mesh(l, nr));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Amr.MeshPlan",
                  "[ParallelAlgorithms][Unit]") {
  const auto d = make_domain();
  auto r = requests();
  // The observed 15/14 mismatch is corrected on every outer shell, without
  // coupling the other object's smaller floor or changing radial proposals.
  propose(r.at(id(0, 2)), 14, 10);
  propose(r.at(id(0, 3)), 23);
  auto plan = amr::planning::make_plan(r, d, policies());
  CHECK(plan.at(id(0, 0)).info.new_mesh == mesh(15));
  CHECK(plan.at(id(0, 1)).info.new_mesh == mesh(15));
  CHECK(plan.at(id(0, 2)).info.new_mesh == mesh(15, 10));
  CHECK(plan.at(id(0, 3)).info.new_mesh == mesh(23));
  CHECK(plan.at(id(1, 0)).info.new_mesh == mesh(12));
  CHECK(plan.at(id(1, 1)).info.new_mesh == mesh(14));
  CHECK(plan.at(id(0, 1)).floor_source == "ExcisionSphereA");
  CHECK(not plan.at(id(0, 0)).floor_source.has_value());

  // Same-cycle anchor refinement and a preexisting multi-mode deficit.
  propose(r.at(id(0, 0)), 16);
  propose(r.at(id(0, 1)), 12);
  plan = amr::planning::make_plan(r, d, policies());
  CHECK(plan.at(id(0, 1)).info.new_mesh == mesh(16));
  CHECK(plan.at(id(0, 2)).info.new_mesh == mesh(16, 10));
  CHECK(plan.at(id(0, 0)).info.new_mesh == mesh(16));

  // Compare retained bands, including differing filter widths. B's anchor
  // retains 12-2=10, so B's 8-mode-filtered shell must reach L=18.
  r.at(id(1, 0)).filtered_angular_modes = 2;
  r.at(id(1, 1)).filtered_angular_modes = 8;
  plan = amr::planning::make_plan(r, d, policies());
  CHECK(plan.at(id(1, 1)).info.new_mesh == mesh(18));

  // A later cycle can lower the floor. Decisions use proposals, not a cached
  // historical maximum or the previous (16) inner mesh.
  for (auto& [element_id, request] : r) {
    request.current_mesh = plan.at(element_id).info.new_mesh;
  }
  for (size_t shell = 0; shell < 3; ++shell) {
    propose(r.at(id(0, shell)), 15);
  }
  plan = amr::planning::make_plan(r, d, policies());
  CHECK(plan.at(id(0, 1)).info.new_mesh == mesh(15));
  CHECK(plan.at(id(0, 1)).info.flags[1] == amr::Flag::DecreaseResolution);
  // Order of policies and arrival/insertion order cannot change the result.
  const auto reversed = amr::planning::make_plan(
      r, d, policies({"ExcisionSphereB", "ExcisionSphereA"}));
  for (const auto& [element_id, decision] : plan) {
    CHECK(decision.info == reversed.at(element_id).info);
  }
  const auto copied = serialize_and_deserialize(r);
  CHECK(amr::planning::make_plan(copied, d, policies()).at(id(0, 1)).info ==
        plan.at(id(0, 1)).info);

  // Disabling the policy leaves proposed meshes unchanged.
  const auto disabled = amr::planning::make_plan(r, d, policies({}));
  for (const auto& [element_id, request] : r) {
    CHECK(disabled.at(element_id).info == request.proposed);
  }
  CHECK_THROWS_WITH(amr::planning::make_plan(r, d, policies({"missing"})),
                    Catch::Matchers::ContainsSubstring("Unknown or repeated"));
  CHECK_THROWS_WITH(amr::planning::make_plan(
                        r, d, policies({"ExcisionSphereA", "ExcisionSphereA"})),
                    Catch::Matchers::ContainsSubstring("Unknown or repeated"));
  auto invalid = r;
  invalid.erase(id(0, 1));
  CHECK_THROWS_WITH(amr::planning::make_plan(invalid, d, policies()),
                    Catch::Matchers::ContainsSubstring("Missing neighbor"));
  invalid = r;
  invalid.at(id(0, 0)).proposed.flags[0] = amr::Flag::Split;
  CHECK_THROWS_WITH(amr::planning::make_plan(invalid, d, policies()),
                    Catch::Matchers::ContainsSubstring("p-only"));
  invalid = r;
  invalid.at(id(0, 0)).filtered_angular_modes.reset();
  CHECK_THROWS_WITH(
      amr::planning::make_plan(invalid, d, policies()),
      Catch::Matchers::ContainsSubstring("Unknown or invalid angular filter"));
  // Stop at the transition to nonspherical wedges; their p requests stay local.
  auto transition = requests();
  const ElementId<3> wedge0{2};
  const ElementId<3> wedge1{3};
  auto shell_neighbors = transition.at(id(0, 3)).element.neighbors();
  shell_neighbors.emplace(
      Direction<3>::upper_xi(),
      Neighbors<3>{std::unordered_set{wedge0, wedge1},
                   {{2, OrientationMap<3>::create_aligned()},
                    {3, OrientationMap<3>::create_aligned()}},
                   false});
  const auto topology = transition.at(id(0, 3)).element.topologies();
  transition.at(id(0, 3)).element =
      Element<3>{id(0, 3), shell_neighbors, topology};
  const Mesh<3> wedge_mesh{5, Spectral::Basis::Legendre,
                           Spectral::Quadrature::GaussLobatto};
  for (const auto& wedge : {wedge0, wedge1}) {
    const Element<3> element{
        wedge,
        {{Direction<3>::lower_xi(),
          Neighbors<3>{std::unordered_set{id(0, 3)},
                       OrientationMap<3>::create_aligned()}}}};
    transition.emplace(
        wedge,
        amr::planning::Request<3>{
            element, wedge_mesh,
            amr::planning::p_decision(wedge_mesh, wedge_mesh), std::nullopt});
  }
  const auto transition_plan =
      amr::planning::make_plan(transition, d, policies());
  CHECK(transition_plan.at(wedge0).info.new_mesh == wedge_mesh);
  CHECK(transition_plan.at(wedge1).info.new_mesh == wedge_mesh);
  CHECK(transition_plan.at(id(0, 3)).info.new_mesh == mesh(15));
  auto cancel = requests();
  cancel.at(id(0, 1)).current_mesh = mesh(15);
  propose(cancel.at(id(0, 1)), 14);
  CHECK(amr::planning::make_plan(cancel, d, policies())
            .at(id(0, 1))
            .info.flags[1] == amr::Flag::DoNothing);
  propose(r.at(id(0, 3)), 15);
  const amr::Policies limited{
      amr::Isotropy::Anisotropic,
      amr::Limits{std::nullopt, std::nullopt, std::nullopt,
                  std::array<size_t, 2>{2, 17}, false},
      false,
      true,
      {"ExcisionSphereA", "ExcisionSphereB"}};
  CHECK_THROWS_WITH(amr::planning::make_plan(r, d, limited),
                    Catch::Matchers::ContainsSubstring("maximum resolution"));
}
