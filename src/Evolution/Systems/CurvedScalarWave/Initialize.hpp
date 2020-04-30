// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>  // IWYU pragma: keep
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Norms.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Initialization/DiscontinuousGalerkin.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
#include "Evolution/Systems/CurvedScalarWave/Constraints.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits.hpp"

namespace CurvedScalarWave {
namespace Actions {
/// \ingroup InitializationGroup
/// \brief Initialize items adding constraint damping parameters of the
/// CurvedScalarWave system
///
/// DataBox changes:
/// - Adds:
///   * `CurvedScalarWave::Tags::ConstraintGamma1`
///   * `CurvedScalarWave::Tags::ConstraintGamma2`
/// - Removes: nothing
/// - Modifies: nothing
///
/// \note This action relies on the `SetupDataBox` aggregated initialization
/// mechanism, so `Actions::SetupDataBox` must be present in the
/// `Initialization` phase action list prior to this action.
struct InitializeConstraintDampingGammas {
  using compute_tags =
      tmpl::list<Tags::ConstraintGamma1Compute, Tags::ConstraintGamma2Compute>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    return std::make_tuple(std::move(box));
  }
};

/// \ingroup InitializationGroup
/// \brief Initialize items related to constraints of the CurvedScalarWave
/// system
///
/// We add background spacetime geometry variables to the evolution databox.
///
/// DataBox changes:
/// - Adds:
///   * `CurvedScalarWave::Tags::ConstraintGamma1`
///   * `CurvedScalarWave::Tags::ConstraintGamma2`
///   * `gr:Tags::SpatialChristoffelFirstKind`
///   * `gr:Tags::SpatialChristoffelSecondKind`
///   * `gr:Tags::TraceSpatialChristoffelSecondKind`
///   * `gr:Tags::TraceExtrinsicCurvature`
///   * `CurvedScalarWave::System::spacetime_variables_tag`
/// - Removes: nothing
/// - Modifies: nothing
///
/// \note This action relies on the `SetupDataBox` aggregated initialization
/// mechanism, so `Actions::SetupDataBox` must be present in the
/// `Initialization` phase action list prior to this action.
template <size_t Dim>
struct InitializeGrVars {
  using compute_tags =
      tmpl::list<gr::Tags::SpatialChristoffelFirstKindCompute<
                     Dim, Frame::Inertial, DataVector>,
                 gr::Tags::SpatialChristoffelSecondKindCompute<
                     Dim, Frame::Inertial, DataVector>,
                 gr::Tags::TraceSpatialChristoffelSecondKindCompute<
                     Dim, Frame::Inertial, DataVector>,
                 GeneralizedHarmonic::Tags::TraceExtrinsicCurvatureCompute<
                     Dim, Frame::Inertial>>;

  using simple_tags = tmpl::list<
      typename CurvedScalarWave::System<Dim>::spacetime_variables_tag>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const double initial_time = db::get<Initialization::Tags::InitialTime>(box);
    const size_t num_grid_points =
        db::get<domain::Tags::Mesh<Dim>>(box).number_of_grid_points();
    const auto inertial_coords =
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);

    using GrVars =
        typename Metavariables::system::spacetime_variables_tag::type;

    // Set initial data from analytic solution
    GrVars gr_vars{num_grid_points};
    gr_vars.assign_subset(evolution::initial_data(
        Parallel::get<::Tags::AnalyticSolutionOrData>(cache), inertial_coords,
        initial_time, typename GrVars::tags_list{}));

    Initialization::mutate_assign<simple_tags>(make_not_null(&box),
                                               std::move(gr_vars));
    return std::make_tuple(std::move(box));
  }
};

/// \ingroup InitializationGroup
/// \brief Initialize items related to constraints of the CurvedScalarWave
/// system
///
/// We add constraints to the evolution databox.
///
/// DataBox changes:
/// - Adds:
///   * `CurvedScalarWave::Tags::OneIndexConstraint<Dim>`
///   * `CurvedScalarWave::Tags::TwoIndexConstraint<Dim>`
///   * `::Tags::PointwiseL2Norm<
///                  CurvedScalarWave::Tags::OneIndexConstraint<Dim>>`
///   * `::Tags::PointwiseL2Norm<
///                  CurvedScalarWave::Tags::TwoIndexConstraint<Dim>>`
/// - Removes: nothing
/// - Modifies: nothing
///
/// \note This action relies on the `SetupDataBox` aggregated initialization
/// mechanism, so `Actions::SetupDataBox` must be present in the
/// `Initialization` phase action list prior to this action.
template <size_t Dim>
struct InitializeConstraints {
  using compute_tags =
      tmpl::list<Tags::OneIndexConstraintCompute<Dim>,
                 Tags::TwoIndexConstraintCompute<Dim>,
                 // following tags added to observe constraints
                 ::Tags::PointwiseL2NormCompute<Tags::OneIndexConstraint<Dim>>,
                 ::Tags::PointwiseL2NormCompute<Tags::TwoIndexConstraint<Dim>>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    return std::make_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace CurvedScalarWave
