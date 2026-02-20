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
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/FilterKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
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

struct FilterKokkosMetavars {
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

  return Domain<dim>{std::move(maps), {}, {"Block0"},
                     std::unordered_map<std::string,
                                        std::unordered_set<std::string>>{}};
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.GeneralizedHarmonic.FilterKokkos",
                  "[Unit][Evolution]") {
  using spacetime_metric_tag =
      gr::Tags::SpacetimeMetric<DataVector, dim, Frame::Inertial>;
  using pi_tag = gh::Tags::Pi<DataVector, dim, Frame::Inertial>;
  using phi_tag = gh::Tags::Phi<DataVector, dim, Frame::Inertial>;
  using device_variables_tag = gh::KokkosTags::DeviceVariables<gh_system>;
  using device_filter_action = gh::Actions::FilterKokkos<filter_type>;
  using host_filter_action =
      dg::Actions::Filter<filter_type,
                          typename gh_system::variables_tag::tags_list>;

  const Mesh<dim> mesh{5, Spectral::Basis::Legendre,
                       Spectral::Quadrature::GaussLobatto};
  const Element<dim> element{ElementId<dim>{0}, {}};

  typename gh_system::variables_tag::type host_vars{
      mesh.number_of_grid_points(), 0.0};
  auto& spacetime_metric = get<spacetime_metric_tag>(host_vars);
  auto& pi = get<pi_tag>(host_vars);
  auto& phi = get<phi_tag>(host_vars);
  for (size_t c = 0; c < spacetime_metric.size(); ++c) {
    for (size_t s = 0; s < mesh.number_of_grid_points(); ++s) {
      spacetime_metric[c][s] = 0.1 * static_cast<double>(c) +
                               static_cast<double>(s);
      pi[c][s] = 10.0 + 0.2 * static_cast<double>(c) + static_cast<double>(s);
    }
  }
  for (size_t c = 0; c < phi.size(); ++c) {
    for (size_t s = 0; s < mesh.number_of_grid_points(); ++s) {
      phi[c][s] = 20.0 + 0.3 * static_cast<double>(c) + static_cast<double>(s);
    }
  }

  const filter_type filter{18.0, 8, true, std::nullopt};

  auto box =
      db::create<db::AddSimpleTags<domain::Tags::Mesh<dim>,
                                   domain::Tags::Element<dim>,
                                   device_variables_tag>>(mesh, element,
                                                          copy_to_device(
                                                              host_vars));
  auto host_box = db::create<db::AddSimpleTags<
      domain::Tags::Mesh<dim>, domain::Tags::Element<dim>,
      typename gh_system::variables_tag>>(mesh, element, host_vars);
  tuples::TaggedTuple<> inboxes{};

  Parallel::GlobalCache<FilterKokkosMetavars> cache{
      typename Parallel::GlobalCache<FilterKokkosMetavars>::ConstTagsTuple{
          make_domain(), filter}};

  const auto result_host = host_filter_action::apply(
      host_box, inboxes, cache, size_t{0}, tmpl::list<>{},
      static_cast<const int*>(nullptr));
  CHECK(std::get<0>(result_host) == Parallel::AlgorithmExecution::Continue);
  CHECK_FALSE(std::get<1>(result_host).has_value());

  const auto result_device = device_filter_action::apply(
      box, inboxes, cache, size_t{0}, tmpl::list<>{},
      static_cast<const int*>(nullptr));
  CHECK(std::get<0>(result_device) == Parallel::AlgorithmExecution::Continue);
  CHECK_FALSE(std::get<1>(result_device).has_value());

  typename gh_system::variables_tag::type filtered_host{
      mesh.number_of_grid_points(), 0.0};
  copy_to_host(make_not_null(&filtered_host),
               db::get<device_variables_tag>(box));
  const auto& expected_host =
      db::get<typename gh_system::variables_tag>(host_box);

  CHECK_ITERABLE_APPROX(get<spacetime_metric_tag>(filtered_host),
                        get<spacetime_metric_tag>(expected_host));
  CHECK_ITERABLE_APPROX(get<pi_tag>(filtered_host), get<pi_tag>(expected_host));
  CHECK_ITERABLE_APPROX(get<phi_tag>(filtered_host),
                        get<phi_tag>(expected_host));
}
