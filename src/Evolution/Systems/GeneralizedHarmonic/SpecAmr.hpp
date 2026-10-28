// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SphericalShellPowerMonitor.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Filter.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::gauges {
class GaugeCondition;
}

namespace gh::spec_amr {
/// The independent spectral extents of one element, as in SpEC's
/// IndependentGridExtents. A spherical shell has one radial and one angular
/// extent, even though its Mesh has three dimensions.
class RefinementVariables {
 public:
  enum class Kind { Polynomial, SphericalHarmonic };

  /// Default to three independent polynomial orders (the existing controller).
  RefinementVariables() = default;
  explicit RefinementVariables(const Mesh<3>& mesh);

  size_t size() const { return number_of_variables_; }
  Kind kind(size_t variable) const;
  size_t variable_for_dimension(size_t dimension) const;

  /// Number of spectral modes per independent variable: N for a polynomial
  /// and l_max + 1 for S2. Unused storage entries are zero.
  std::array<size_t, 3> mode_counts(const Mesh<3>& mesh) const;
  /// Expand independent decisions to mesh directions. For S2 the same angular
  /// decision controls both theta and phi; h changes in S2 are forbidden.
  std::array<amr::Flag, 3> logical_flags(
      const std::array<amr::Flag, 3>& independent_flags) const;
  /// Inverse mapping for an existing request. Coupled directions must agree.
  std::array<amr::Flag, 3> independent_flags(
      const std::array<amr::Flag, 3>& logical_flags) const;

  void pup(PUP::er& p);
  friend bool operator==(const RefinementVariables&,
                         const RefinementVariables&) = default;

 private:
  size_t number_of_variables_{3};
  std::array<Kind, 3> kinds_{
      {Kind::Polynomial, Kind::Polynomial, Kind::Polynomial}};
  std::array<size_t, 3> variable_for_dimension_{{0, 1, 2}};
};

/// Relative high-mode amplitude, including the prospective filtered p-down
/// spectrum. Pile-up follows SpEC's multi-fit diagnostic, excluding roundoff
/// and the highest filtered mode. Mode indices are retained across parity gaps.
struct Spectrum {
  double error{0.0};
  double next_error{0.0};
  double convergence{0.0};
  double pile_up{0.0};
  bool degenerate{true};
};
Spectrum diagnose_power(const DataVector& power, size_t filtered_modes = 1);

/// Source-faithful SpEC diagnostic: compact exact zeros, retain their original
/// indices for fits, and round the effective filter count. The optional count
/// describes the same local filter evaluated at one lower independent extent.
/// The legacy diagnose_power remains available to the hybrid controller.
Spectrum diagnose_spec_power(
    const DataVector& power, double filtered_modes,
    std::optional<double> after_p_down_filtered_modes = std::nullopt);

struct FilteredModes {
  std::array<double, 3> current{};
  std::array<double, 3> after_p_down{};
  std::array<bool, 3> can_decrease{{true, true, true}};
};
/// Sum 1-response per independent topology, using the actual element filter
/// both at the current mesh and at each prospective one-mode-lower mesh.
FilteredModes filtered_modes(
    const Mesh<3>& mesh,
    const Filters::Filter<
        3, ylm::TensorYlm::filter_detail::gh_spacetime_vars_list>& filter);

/// One-direction polynomial transforms, RMS over transverse nodal points.
std::array<DataVector, 3> polynomial_power(const DataVector& component,
                                           const Mesh<3>& mesh);
/// SpEC's final topological aggregation (including its outer square root).
std::array<double, 3> aggregate_constraints(
    const std::array<double, 3>& directional_rms,
    const RefinementVariables& variables);

struct Diagnostics {
  RefinementVariables refinement_variables{};
  // Metric and combined Pi/Phi, respectively. diagnose_spec averages transverse
  // nodal points on I1 and tensor-harmonic coefficients on S2. The legacy
  // diagnose uses full tensor-product modal averaging.
  // Arrays are indexed by independent refinement variable, not mesh dimension.
  // Only refinement_variables.size() entries are active. In particular, an S2
  // spectrum must appear once, not once each for theta and phi.
  std::array<std::array<Spectrum, 3>, 2> spectra{};
  // Projected C_ia and C_iab. diagnose_spec uses SpEC's topological aggregation
  // of directional nodal RMS; the legacy diagnose returns directional RMS.
  // Each independent remaining tensor component is counted once.
  std::array<std::array<double, 3>, 2> constraints{};
  FilteredModes filter{};
};
Diagnostics diagnose(const tnsr::aa<DataVector, 3>& metric,
                     const tnsr::aa<DataVector, 3>& pi,
                     const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
                     const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                                           Frame::Inertial>& inverse_jacobian,
                     const tnsr::I<DataVector, 3>& coordinates, double time,
                     const gauges::GaugeCondition& gauge,
                     const Scalar<DataVector>& gamma2);

/// SpEC diagnostics for polynomial and I1 x S2 elements. Angular spectra use
/// tensor harmonics in the grid frame; radial spectra use the inertial tensor
/// components. The transient transform cache is bounded to the current l_max.
/// Projected constraints follow NormOfProjectedDerivConstraint, including its
/// separate element/grid and grid/inertial Jacobian normalization.
Diagnostics diagnose_spec(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Jacobian<DataVector, 3, Frame::Grid, Frame::Inertial>& grid_jacobian,
    const tnsr::I<DataVector, 3>& coordinates, double time,
    const gauges::GaugeCondition& gauge, const Scalar<DataVector>& gamma2,
    const FilteredModes& filter,
    gsl::not_null<power_monitor::CartToSphereMatrices*> transform_cache);

/// Serialized per-element history. Reset on h changes; retained on p changes
/// so the next call can detect the *realized* mesh change. Observations do not
/// update it. Target adjustment retains constraint history across p changes,
/// as in SpEC, including in the prospective-removal forecast.
struct History {
  RefinementVariables refinement_variables{};
  // Transient, per-element cache. Sharing avoids copying sparse matrices when
  // evaluating on a history snapshot. It carries no numerical history and is
  // rebuilt after migration; global criteria must not cache element meshes.
  std::shared_ptr<power_monitor::CartToSphereMatrices> transform_cache{};
  // Startup iterations do not advance evolution history or its cooldowns.
  bool initializing{false};
  bool initial_cycle_detected{false};
  std::vector<std::array<size_t, 3>> initial_meshes{};
  // Cached grid-frame radius-plus-time restriction, retained across p changes.
  double radius_release_time{std::numeric_limits<double>::signaling_NaN()};
  size_t calls{0};
  double time{-std::numeric_limits<double>::infinity()};
  // Extents and all directional histories use independent-variable indices.
  // The third entry is unused on a spherical shell. decision is expanded into
  // logical mesh directions only at the controller boundary.
  std::array<size_t, 3> extents{};
  std::array<size_t, 3> last_change{};
  std::array<double, 3> last_change_time{
      {-std::numeric_limits<double>::infinity(),
       -std::numeric_limits<double>::infinity(),
       -std::numeric_limits<double>::infinity()}};
  // Spatial base target, sampled once and retained across p changes.
  double spectral_target{std::numeric_limits<double>::signaling_NaN()};
  std::array<std::array<double, 3>, 2> constraints{};
  std::array<double, 3> log10_target_shift{};
  std::array<amr::Flag, 3> decision{};
  void pup(PUP::er& p) {
    if (p.isUnpacking()) {
      transform_cache.reset();
    }
    p | refinement_variables;
    p | initializing;
    p | initial_cycle_detected;
    p | initial_meshes;
    p | radius_release_time;
    p | calls;
    p | time;
    p | extents;
    p | last_change;
    p | last_change_time;
    p | spectral_target;
    p | constraints;
    p | log10_target_shift;
    p | decision;
  }
};

struct Settings {
  double spectral_target{1.e-4};
  std::array<double, 2> constraint_targets{{1.e-3, 1.e-3}};
  bool force_pile_up{false};
  size_t minimum_modes{5};
  size_t cooldown{2};
  // Zero disables stored target adjustment. The hybrid retains its legacy
  // addition gate in that case; decide_p always requires LargeTrunc.
  // Negative values enable SpEC's local target adjustment (inspiral uses -1).
  double adjustment_limit{0.0};
  double minimum_time_between_changes{0.0};
};

/// Compute the next local log10 target shifts without changing history.
/// Use the most restrictive constraint update and clamp to [limit, 0].
/// Nondecaying or degenerate spectra provide no usable adjustment step.
std::array<double, 3> updated_target_shifts(const History& history,
                                            const Diagnostics& diagnostics,
                                            const Settings& settings);

/// Replace p-coarsening and limit p-addition using SpEC's tolerance and
/// pile-up rules. Preserve Chalmers h-refinement. Existing prospective-parent
/// checks still decide joins; both constraint targets must additionally pass.
/// With target adjustment enabled, p addition and prospective removal use the
/// same shifted spectral target. Chalmers still requests refinement and selects
/// h versus p. Ordinary permission budgets retain the unshifted target, as in
/// SpEC. This omits SpEC's neighbor smoothing, adaptive filters and h
/// controller.
std::array<amr::Flag, 3> decide(gsl::not_null<History*> history,
                                const Diagnostics& diagnostics,
                                const std::array<size_t, 3>& extents,
                                double time,
                                std::array<amr::Flag, 3> field_decision,
                                const Settings& settings);

/// Standalone SpEC-style p control. Both addition and removal are decided
/// from the grouped spectra and projected constraints, without an upstream
/// refinement request. Even with adjustment_limit == 0, addition requires
/// spectral error above the target. No h changes are requested. When
/// History::initializing is set, bypass evolution cooldowns and permit removal
/// on the first call. Detect recurring meshes, and retain only the target shift
/// and startup mesh history. extents contains independent spectral mode counts
/// (RefinementVariables::mode_counts), not the raw mesh extents. Changes remain
/// one spectral mode per iteration. The returned flags use logical dimensions.
std::array<amr::Flag, 3> decide_p(gsl::not_null<History*> history,
                                  const Diagnostics& diagnostics,
                                  const std::array<size_t, 3>& extents,
                                  double time, const Settings& settings);

struct HistoryTag : db::SimpleTag {
  using type = History;
};
struct InitializeHistory {
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<HistoryTag>;
  using compute_tags = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<>;
  static void apply() {}
};
struct ProjectHistory : tt::ConformsTo<amr::protocols::Projector> {
  using return_tags = tmpl::list<HistoryTag>;
  using argument_tags = tmpl::list<>;
  static void apply(gsl::not_null<History*>,
                    const std::pair<Mesh<3>, Element<3>>&) {}
  template <typename ParentOrChildren>
  static void apply(const gsl::not_null<History*> history,
                    const ParentOrChildren&) {
    *history = History{};
  }
};
}  // namespace gh::spec_amr
