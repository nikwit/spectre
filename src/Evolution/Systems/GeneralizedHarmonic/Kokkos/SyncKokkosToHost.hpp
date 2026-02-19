// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Events {

// Event to synchronize device-resident evolved variables back to host
// variables so standard observation events can run unchanged.
template <typename System>
class SyncKokkosToHost : public Event {
 public:
  /// \cond
  explicit SyncKokkosToHost(CkMigrateMessage* /*msg*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(SyncKokkosToHost);  // NOLINT
  /// \endcond

  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "Synchronize GH Kokkos device evolved variables to host memory for "
      "observation."};

  SyncKokkosToHost() = default;

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<typename System::variables_tag>;
  using argument_tags = tmpl::list<gh::KokkosTags::DeviceVariables<System>>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const gsl::not_null<typename System::variables_tag::type*> host_vars,
      const typename gh::KokkosTags::DeviceVariables<System>::type&
          device_vars,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& /*observation_value*/) const {
    copy_to_host(host_vars, device_vars);
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override { Event::pup(p); }
};

template <typename System>
PUP::able::PUP_ID SyncKokkosToHost<System>::my_PUP_ID = 0;  // NOLINT

}  // namespace gh::Events
