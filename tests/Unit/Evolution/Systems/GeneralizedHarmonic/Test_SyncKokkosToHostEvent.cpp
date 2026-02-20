// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct SyncEventMetavars {
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<tmpl::pair<Event, tmpl::list<>>>;
  };
  using component_list = tmpl::list<>;
};
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.GeneralizedHarmonic.SyncKokkosToHostEvent",
    "[Unit][Evolution]") {
  constexpr size_t dim = 1;
  using system = gh::System<dim>;
  using event = gh::Events::SyncKokkosToHost<system>;
  using spacetime_metric_tag =
      gr::Tags::SpacetimeMetric<DataVector, dim, Frame::Inertial>;
  using pi_tag = gh::Tags::Pi<DataVector, dim, Frame::Inertial>;
  using phi_tag = gh::Tags::Phi<DataVector, dim, Frame::Inertial>;

  constexpr size_t num_points = 4;
  typename system::variables_tag::type host_vars{num_points, 0.0};

  auto& spacetime_metric = get<spacetime_metric_tag>(host_vars);
  auto& pi = get<pi_tag>(host_vars);
  auto& phi = get<phi_tag>(host_vars);
  for (size_t c = 0; c < spacetime_metric.size(); ++c) {
    for (size_t s = 0; s < num_points; ++s) {
      spacetime_metric[c][s] = 10.0 * static_cast<double>(c) + s;
      pi[c][s] = 100.0 + 10.0 * static_cast<double>(c) + s;
    }
  }
  for (size_t c = 0; c < phi.size(); ++c) {
    for (size_t s = 0; s < num_points; ++s) {
      phi[c][s] = 200.0 + 10.0 * static_cast<double>(c) + s;
    }
  }

  const auto device_vars = copy_to_device(host_vars);

  typename system::variables_tag::type host_copy_target{num_points, -1.0};

  event sync_event{};
  Parallel::GlobalCache<SyncEventMetavars> cache{
      typename Parallel::GlobalCache<SyncEventMetavars>::ConstTagsTuple{}};
  const Event::ObservationValue observation_value{"Time", 0.0};
  sync_event(make_not_null(&host_copy_target), device_vars, cache, size_t{0},
             static_cast<const void*>(nullptr), observation_value);

  CHECK_ITERABLE_APPROX(get<spacetime_metric_tag>(host_copy_target),
                        get<spacetime_metric_tag>(host_vars));
  CHECK_ITERABLE_APPROX(get<pi_tag>(host_copy_target), get<pi_tag>(host_vars));
  CHECK_ITERABLE_APPROX(get<phi_tag>(host_copy_target),
                        get<phi_tag>(host_vars));
}
