// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryEvolution.hpp"
#include "ParallelAlgorithms/Amr/Protocols/AmrMetavariables.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "PointwiseFunctions/ConstraintDamping/DampingFunction.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/LtsMode.hpp"
#include "Time/Tags/Time.hpp"

namespace gh::worldtube {
using PhysicalTraceBoundaryConditions = std::vector<DirectionMap<
    3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>;

/// Full-sphere prototype; all spatial differentiation is performed on
/// Cartesian components of the reconstructed mixed ghost, before projection
/// into the dyad. There is no scalar-harmonic differentiation of spin-2 B.
DataVector physical_trace_on_sphere(
    bool initialize, const DataVector& amplitudes,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const tnsr::I<DataVector, 3>& coordinates,
    const tnsr::I<DataVector, 3, Frame::Grid>& grid_coordinates,
    const Domain<3>& domain, const Element<3>& element,
    const PhysicalTraceBoundaryConditions& boundary_conditions, double time,
    LtsMode lts_mode,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity,
    const gh::gauges::GaugeCondition& gauge,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma0,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma1,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma2,
    const domain::FunctionsOfTimeMap& functions_of_time);

using PhysicalTraceArguments = tmpl::list<
    gr::Tags::SpacetimeMetric<DataVector, 3>, gh::Tags::Pi<DataVector, 3>,
    gh::Tags::Phi<DataVector, 3>, domain::Tags::Mesh<3>,
    domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>,
    domain::Tags::Coordinates<3, Frame::Inertial>,
    domain::Tags::Coordinates<3, Frame::Grid>, domain::Tags::Domain<3>,
    domain::Tags::Element<3>, domain::Tags::ExternalBoundaryConditions<3>,
    ::Tags::Time, ::Tags::LtsMode, domain::Tags::MeshVelocity<3>,
    gh::gauges::Tags::GaugeCondition,
    gh::Tags::DampingFunctionGamma0<3, Frame::Grid>,
    gh::Tags::DampingFunctionGamma1<3, Frame::Grid>,
    gh::Tags::DampingFunctionGamma2<3, Frame::Grid>,
    domain::Tags::FunctionsOfTime>;

namespace Initialization {
struct PhysicalTrace {
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<Tags::PhysicalTrace>;
  using compute_tags = tmpl::list<>;
  using return_tags = simple_tags;
  using argument_tags = tmpl::list<>;
  static void apply(const gsl::not_null<DataVector*> trace) { *trace = {}; }
};
}  // namespace Initialization

/// Called once after GH initial data and gauge initialization, before
/// self-start.
struct InitializePhysicalTrace {
  using return_tags = tmpl::list<Tags::PhysicalTrace>;
  using argument_tags = PhysicalTraceArguments;
  template <typename... Args>
  static void apply(const gsl::not_null<DataVector*> trace,
                    const Args&... args) {
    *trace = physical_trace_on_sphere(true, DataVector{}, args...);
  }
};

/// Called once per RHS evaluation, before recording the SAME TimeStepId for
/// both volume and surface state. Does not mutate B or any persistent cache.
struct ComputePhysicalTraceRhs {
  using return_tags = tmpl::list<::Tags::dt<Tags::PhysicalTrace>>;
  using argument_tags =
      tmpl::push_front<PhysicalTraceArguments, Tags::PhysicalTrace>;
  template <typename... Args>
  static void apply(const gsl::not_null<DataVector*> dt_trace,
                    const DataVector& trace, const Args&... args) {
    *dt_trace = physical_trace_on_sphere(false, trace, args...);
  }
};

/// Keep empty auxiliary histories consistent for ordinary single-hole runs.
/// Refinement of an active surface trace is deliberately rejected until a
/// spin-2 surface-state and history projection has been implemented.
template <typename BulkVariablesTag>
struct ProjectPhysicalTrace : tt::ConformsTo<amr::protocols::Projector> {
  using return_tags =
      tmpl::list<Tags::PhysicalTrace, ::Tags::dt<Tags::PhysicalTrace>,
                 ::Tags::HistoryEvolvedVariables<Tags::PhysicalTrace>>;
  using argument_tags =
      tmpl::list<::Tags::HistoryEvolvedVariables<BulkVariablesTag>>;
  static void check_empty(const DataVector& value) {
    if (value.size() != 0) {
      ERROR("AMR of an evolved physical boundary trace is not implemented.");
    }
  }
  static void check_old(const std::pair<Mesh<3>, Element<3>>& /*old*/) {}
  template <typename... T>
  static void check_old(const tuples::TaggedTuple<T...>& parent) {
    check_empty(get<Tags::PhysicalTrace>(parent));
  }
  template <typename K, typename V, typename... T>
  static void check_old(const std::unordered_map<K, V, T...>& children) {
    for (const auto& [id, child] : children) {
      (void)id;
      check_old(child);
    }
  }
  template <typename Old>
  static void apply(
      const gsl::not_null<DataVector*> trace,
      const gsl::not_null<DataVector*> dt_trace,
      const gsl::not_null<TimeSteppers::History<DataVector>*> history,
      const typename ::Tags::HistoryEvolvedVariables<BulkVariablesTag>::type&
          bulk,
      const Old& old) {
    check_empty(*trace);
    check_old(old);
    *trace = DataVector{};
    *dt_trace = DataVector{};
    *history = TimeSteppers::History<DataVector>(bulk.integration_order());
    for (const auto& entry : bulk) {
      history->insert(entry.time_step_id, DataVector{}, DataVector{});
    }
    for (const auto& entry : bulk.substeps()) {
      history->insert(entry.time_step_id, DataVector{}, DataVector{});
    }
  }
};
}  // namespace gh::worldtube
