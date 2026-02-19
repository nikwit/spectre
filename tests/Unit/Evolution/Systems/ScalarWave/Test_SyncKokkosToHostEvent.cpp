// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
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

SPECTRE_TEST_CASE("Unit.Evolution.Systems.ScalarWave.SyncKokkosToHostEvent",
                  "[Unit][Evolution]") {
  using system = ScalarWave::System<1>;
  using event = ScalarWave::Events::SyncKokkosToHost<system>;

  constexpr size_t num_points = 4;
  typename system::variables_tag::type host_vars{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_vars)) = DataVector{{1.0, 2.0, 3.0, 4.0}};
  get(get<ScalarWave::Tags::Pi>(host_vars)) = DataVector{{5.0, 6.0, 7.0, 8.0}};
  get<ScalarWave::Tags::Phi<1>>(host_vars).get(0) =
      DataVector{{9.0, 10.0, 11.0, 12.0}};

  auto device_vars = copy_to_device(host_vars);

  typename system::variables_tag::type host_copy_target{num_points, 0.0};
  get(get<ScalarWave::Tags::Psi>(host_copy_target)) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};
  get(get<ScalarWave::Tags::Pi>(host_copy_target)) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};
  get<ScalarWave::Tags::Phi<1>>(host_copy_target).get(0) =
      DataVector{{-1.0, -1.0, -1.0, -1.0}};

  event sync_event{};
  Parallel::GlobalCache<SyncEventMetavars> cache{
      typename Parallel::GlobalCache<SyncEventMetavars>::ConstTagsTuple{}};
  const Event::ObservationValue observation_value{"Time", 0.0};
  sync_event(make_not_null(&host_copy_target), device_vars, cache, size_t{0},
             static_cast<const void*>(nullptr), observation_value);

  CHECK_ITERABLE_APPROX(get(get<ScalarWave::Tags::Psi>(host_copy_target)),
                        get(get<ScalarWave::Tags::Psi>(host_vars)));
  CHECK_ITERABLE_APPROX(get(get<ScalarWave::Tags::Pi>(host_copy_target)),
                        get(get<ScalarWave::Tags::Pi>(host_vars)));
  CHECK_ITERABLE_APPROX(get<ScalarWave::Tags::Phi<1>>(host_copy_target),
                        get<ScalarWave::Tags::Phi<1>>(host_vars));
}
