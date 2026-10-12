// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <pup.h>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Structure/Element.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::gauges {
class GaugeCondition;
}

namespace gh::spec_amr {
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

struct Diagnostics {
  // Metric and combined Pi/Phi, respectively. RMS over independent tensor
  // components and transverse modal indices (not directional Legendre L2).
  std::array<std::array<Spectrum, 3>, 2> spectra{};
  // Nodal RMS of directionally projected C_ia and C_iab. As in SpEC's
  // component sum, each independent tensor component is counted once.
  std::array<std::array<double, 3>, 2> constraints{};
};
Diagnostics diagnose(const tnsr::aa<DataVector, 3>& metric,
                     const tnsr::aa<DataVector, 3>& pi,
                     const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
                     const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                                           Frame::Inertial>& inverse_jacobian,
                     const tnsr::I<DataVector, 3>& coordinates, double time,
                     const gauges::GaugeCondition& gauge,
                     const Scalar<DataVector>& gamma2);

/// Serialized per-element history. Reset on h changes; retained on p changes
/// so the next call can detect the *realized* mesh change. Observations do not
/// update it. Target adjustment retains constraint history across p changes,
/// as in SpEC; the prospective-removal forecast uses only an unchanged mesh.
struct History {
  size_t calls{0};
  double time{-std::numeric_limits<double>::infinity()};
  std::array<size_t, 3> extents{};
  std::array<size_t, 3> last_change{};
  std::array<std::array<double, 3>, 2> constraints{};
  std::array<double, 3> log10_target_shift{};
  std::array<amr::Flag, 3> decision{};
  void pup(PUP::er& p) {
    p | calls;
    p | time;
    p | extents;
    p | last_change;
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
  // Zero preserves the original coarsening experiment. Negative values enable
  // SpEC's local, history-dependent target adjustment (inspiral uses -1).
  double adjustment_limit{0.0};
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
