// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/FilterKokkosBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Tags/Filter.hpp"
#include "NumericalAlgorithms/LinearOperators/ExponentialFilter.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Actions/FilterAction.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {

constexpr size_t dim = 3;
using gh_system = gh::System<dim>;
using filter_type = Filters::Exponential<0>;

struct FilterKokkosBatchedMetavars {
  using system = gh_system;
  using component_list = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<dim>, Filters::Tags::Filter<filter_type>>;
};

Domain<dim> make_domain() {
  using identity_map = domain::CoordinateMaps::Identity<dim>;
  using map =
      domain::CoordinateMap<Frame::BlockLogical, Frame::Inertial, identity_map>;
  register_classes_with_charm(tmpl::list<map>{});

  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::BlockLogical, Frame::Inertial, dim>>>
      maps{};
  maps.emplace_back(std::make_unique<map>(identity_map{}));

  return Domain<dim>{
      std::move(maps),
      {},
      {"Block0"},
      std::unordered_map<std::string, std::unordered_set<std::string>>{}};
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.FilterKokkosBatched",
    "[Unit][Evolution]") {
  using spacetime_metric_tag =
      gr::Tags::SpacetimeMetric<DataVector, dim, Frame::Inertial>;
  using pi_tag = gh::Tags::Pi<DataVector, dim, Frame::Inertial>;
  using phi_tag = gh::Tags::Phi<DataVector, dim, Frame::Inertial>;
  using batched_filter_action = gh::Actions::FilterKokkosBatched<filter_type>;
  using host_filter_action =
      dg::Actions::Filter<filter_type,
                          typename gh_system::variables_tag::tags_list>;

  const Mesh<dim> mesh{std::array<size_t, dim>{{4, 3, 5}},
                       Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  const size_t num_points = mesh.number_of_grid_points();
  const Element<dim> element{ElementId<dim>{0}, {}};

  typename gh_system::variables_tag::type host_vars{num_points, 0.0};
  auto& spacetime_metric = get<spacetime_metric_tag>(host_vars);
  auto& pi = get<pi_tag>(host_vars);
  auto& phi = get<phi_tag>(host_vars);
  for (size_t c = 0; c < spacetime_metric.size(); ++c) {
    for (size_t s = 0; s < num_points; ++s) {
      spacetime_metric[c][s] =
          0.1 * static_cast<double>(c) + static_cast<double>(s);
      pi[c][s] = 10.0 + 0.2 * static_cast<double>(c) + static_cast<double>(s);
    }
  }
  for (size_t c = 0; c < phi.size(); ++c) {
    for (size_t s = 0; s < num_points; ++s) {
      phi[c][s] = 20.0 + 0.3 * static_cast<double>(c) + static_cast<double>(s);
    }
  }

  const filter_type filter{18.0, 8, true, std::nullopt};

  auto host_box = db::create<
      db::AddSimpleTags<domain::Tags::Mesh<dim>, domain::Tags::Element<dim>,
                        typename gh_system::variables_tag>>(mesh, element,
                                                            host_vars);
  tuples::TaggedTuple<> inboxes{};
  Parallel::GlobalCache<FilterKokkosBatchedMetavars> cache{
      typename Parallel::GlobalCache<
          FilterKokkosBatchedMetavars>::ConstTagsTuple{make_domain(), filter}};

  const auto result_host = host_filter_action::apply(
      host_box, inboxes, cache, size_t{0}, tmpl::list<>{},
      static_cast<const int*>(nullptr));
  CHECK(std::get<0>(result_host) == Parallel::AlgorithmExecution::Continue);
  CHECK_FALSE(std::get<1>(result_host).has_value());

  evolution::Kokkos::PackedTopology<gh_system> packed_topology{};
  packed_topology.local_element_ids.push_back(element.id());
  packed_topology.local_elements.push_back(element);
  packed_topology.element_index_by_id.insert_or_assign(element.id(), 0);
  packed_topology.element_extents_host.push_back(mesh.extents().indices());
  packed_topology.element_point_offsets_host = {0, num_points};
  packed_topology.total_points = num_points;
  packed_topology.points_per_element = num_points;
  packed_topology.uniform_extents_host = mesh.extents().indices();
  packed_topology.uniform_basis_host = {mesh.basis(0), mesh.basis(1),
                                        mesh.basis(2)};
  packed_topology.uniform_quadrature_host = {
      mesh.quadrature(0), mesh.quadrature(1), mesh.quadrature(2)};

  evolution::Kokkos::PackedEvolutionState<gh_system> packed_evolution_state{};
  packed_evolution_state.device_variables = copy_to_device(host_vars);
  evolution::Kokkos::PackedBoundaryScratch<gh_system> packed_boundary_scratch{};

  batched_filter_action::apply(make_not_null(&packed_evolution_state),
                               make_not_null(&packed_boundary_scratch),
                               packed_topology, filter);

  typename gh_system::variables_tag::type filtered_batched{num_points, 0.0};
  copy_to_host(make_not_null(&filtered_batched),
               packed_evolution_state.device_variables);
  const auto& expected_host =
      db::get<typename gh_system::variables_tag>(host_box);

  CHECK_ITERABLE_APPROX(get<spacetime_metric_tag>(filtered_batched),
                        get<spacetime_metric_tag>(expected_host));
  CHECK_ITERABLE_APPROX(get<pi_tag>(filtered_batched),
                        get<pi_tag>(expected_host));
  CHECK_ITERABLE_APPROX(get<phi_tag>(filtered_batched),
                        get<phi_tag>(expected_host));
}
