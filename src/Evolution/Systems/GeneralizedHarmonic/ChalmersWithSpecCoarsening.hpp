// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

#include "Evolution/Systems/GeneralizedHarmonic/ChalmersConstraintFeedback.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"

namespace gh {
/// Experimental SpEC-style p coarsening around the existing GH hp criterion.
/// Its relative tensor targets are independent of Chalmers' componentwise
/// absolute/relative refinement budgets. Optional local target adjustment
/// applies the same stored factor to both budgets, replacing instantaneous GH
/// feedback. Chalmers still requests refinement and selects h versus p.
class ChalmersWithSpecCoarsening : public amr::Criterion {
 public:
  struct RefinementCriterion {
    using type = ChalmersWithConstraintFeedback<3>;
    static constexpr Options::String help = "Existing GH hp refinement policy.";
  };
  struct RelativeTruncationTarget {
    using type = double;
    static constexpr Options::String help = "Relative tensor spectral target.";
  };
  struct ConstraintTargets {
    using type = std::array<double, 2>;
    static constexpr Options::String help =
        "Projected nodal RMS targets for C_ia and C_iab, respectively.";
  };
  struct EnableCoarsening {
    using type = bool;
    static constexpr Options::String help =
        "Apply the new policy; false only records diagnostics.";
  };
  struct ForceDropPileUpModes {
    using type = bool;
    static constexpr Options::String help =
        "Permit SpEC's pile-up override with threshold 3 and add threshold "
        "1.5.";
  };
  struct MinimumModes {
    using type = size_t;
    static type lower_bound() { return 7; }
    static constexpr Options::String help =
        "Minimum extent after p-coarsening: three modes per parity plus the "
        "filtered mode requires at least seven.";
  };
  struct CallsBetweenChanges {
    using type = size_t;
    static type lower_bound() { return 1; }
    static constexpr Options::String help =
        "AMR-call cooldown after p changes.";
  };
  struct StartTime {
    using type = double;
    static constexpr Options::String help =
        "Apply the policy from this time, allowing identical early controls.";
  };
  struct Verbose {
    using type = bool;
    static constexpr Options::String help = "Record pre-AMR diagnostics.";
  };
  struct TruncationErrorAdjustmentLimit {
    using type = double;
    static type upper_bound() { return 0.0; }
    static constexpr Options::String help =
        "Minimum stored log10 target shift. Zero retains the legacy policy; "
        "negative values replace instantaneous feedback with SpEC's local "
        "controller. Inspiral uses -1. No neighbor smoothing is applied.";
  };
  using options =
      tmpl::list<RefinementCriterion, RelativeTruncationTarget,
                 ConstraintTargets, EnableCoarsening, ForceDropPileUpModes,
                 MinimumModes, CallsBetweenChanges, StartTime, Verbose,
                 TruncationErrorAdjustmentLimit>;
  static constexpr Options::String help =
      "SpEC-style tolerance and pile-up coarsening on filtered Legendre cubes.";
  ChalmersWithSpecCoarsening() = default;
  ChalmersWithSpecCoarsening(ChalmersWithConstraintFeedback<3> refinement,
                             const double spectral_target,
                             const std::array<double, 2>& constraint_targets,
                             const bool enabled, const bool force,
                             const size_t minimum_modes, const size_t cooldown,
                             const double start_time, const bool verbose,
                             const double adjustment_limit = 0.0,
                             const Options::Context& context = {})
      : refinement_(std::move(refinement)),
        settings_{spectral_target, constraint_targets, force,
                  minimum_modes,   cooldown,           adjustment_limit},
        enabled_(enabled),
        start_time_(start_time),
        verbose_(verbose) {
    if (not std::isfinite(spectral_target) or spectral_target <= 0.0 or
        not std::isfinite(start_time) or minimum_modes < 7 or cooldown < 1) {
      PARSE_ERROR(context,
                  "Invalid spectral target, start time, minimum modes "
                  "or cooldown for SpEC coarsening.");
    }
    for (const double target : constraint_targets) {
      if (not std::isfinite(target) or target <= 0.0) {
        PARSE_ERROR(context, "Constraint targets must be finite and positive.");
      }
    }
    const double minimum_log = log10(std::numeric_limits<double>::min());
    if (not std::isfinite(adjustment_limit) or adjustment_limit > 0.0 or
        adjustment_limit < minimum_log or
        log10(spectral_target) + adjustment_limit < minimum_log) {
      PARSE_ERROR(context,
                  "TruncationErrorAdjustmentLimit must be finite and "
                  "nonpositive and leave a representable positive target.");
    }
  }
  explicit ChalmersWithSpecCoarsening(CkMigrateMessage* msg)
      : amr::Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ChalmersWithSpecCoarsening);  // NOLINT
  amr::Criteria::Type type() override { return amr::Criteria::Type::hp; }
  std::string observation_name() override {
    return "ChalmersWithSpecCoarsening";
  }
  bool needs_sibling_data() const override {
    return refinement_.needs_sibling_data();
  }
  template <typename DbTags>
  amr::CriterionTensorData pack_sibling_data(
      const db::DataBox<DbTags>& box) const {
    return refinement_.pack_sibling_data(box);
  }
  using compute_tags_for_observation_box = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, 3> operator()(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<3>& id) const {
    auto history = db::get<spec_amr::HistoryTag>(box);
    return evaluate_impl(make_not_null(&history), box, cache, id, false);
  }

  template <typename ComputeTags, typename DbType, typename Metavariables>
  std::array<amr::Flag, 3> evaluate_with_history(
      const gsl::not_null<ObservationBox<ComputeTags, DbType>*> box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<3>& id) const {
    auto history = db::get<spec_amr::HistoryTag>(box->databox());
    const auto flags =
        evaluate_impl(make_not_null(&history), box->databox(), cache, id, true);
    db::mutate<spec_amr::HistoryTag>(
        [&history](const auto h) { *h = std::move(history); },
        make_not_null(&box->databox()));
    box->reset();
    return flags;
  }
  void pup(PUP::er& p) override {
    amr::Criterion::pup(p);
    p | refinement_;
    p | settings_.spectral_target;
    p | settings_.constraint_targets;
    p | settings_.force_pile_up;
    p | settings_.minimum_modes;
    p | settings_.cooldown;
    p | settings_.adjustment_limit;
    p | enabled_;
    p | start_time_;
    p | verbose_;
  }

 private:
  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, 3> evaluate_impl(
      const gsl::not_null<spec_amr::History*> history,
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& cache, const ElementId<3>& id,
      const bool log) const {
    const auto& mesh = db::get<domain::Tags::Mesh<3>>(box);
    for (size_t d = 0; d < 3; ++d) {
      if (mesh.extents(d) < 7 or mesh.basis(d) != Spectral::Basis::Legendre) {
        ERROR("The parity refinement wrapper requires Legendre extents >= 7 in "
              << id);
      }
    }
    const double time = db::get<::Tags::Time>(box);
    if (history->calls > 0 and time == history->time and
        mesh.extents().indices() == history->extents) {
      return history->decision;
    }
    const auto diagnostics = spec_amr::diagnose(
        db::get<gr::Tags::SpacetimeMetric<DataVector, 3>>(box),
        db::get<gh::Tags::Pi<DataVector, 3>>(box),
        db::get<gh::Tags::Phi<DataVector, 3>>(box), mesh,
        db::get<domain::Tags::InverseJacobian<3, Frame::ElementLogical,
                                              Frame::Inertial>>(box),
        db::get<domain::Tags::Coordinates<3, Frame::Inertial>>(box), time,
        db::get<gauges::Tags::GaugeCondition>(box),
        db::get<gh::Tags::ConstraintGamma2>(box));
    const bool active = enabled_ and time >= start_time_;
    auto settings = settings_;
    if (not active) {
      // Do not accumulate a target shift before the comparison switches on.
      settings.adjustment_limit = 0.0;
    }
    const auto shifts =
        spec_amr::updated_target_shifts(*history, diagnostics, settings);
    const auto before = [&]() {
      if (settings.adjustment_limit == 0.0) {
        return refinement_(box, cache, id);
      }
      std::array<double, 3> factors{};
      for (size_t d = 0; d < 3; ++d) {
        factors[d] = pow(10.0, shifts[d]);
      }
      return refinement_.evaluate_with_target_factors(box, cache, id, factors);
    }();
    auto result = spec_amr::decide(
        history, diagnostics, mesh.extents().indices(), time, before, settings);
    if (not active) {
      result = before;
      history->decision = result;
    }
    if (verbose_ and log) {
      for (size_t d = 0; d < 3; ++d) {
        const auto& g = diagnostics.spectra[0][d];
        const auto& k = diagnostics.spectra[1][d];
        Parallel::printf(
            "SpecCoarsening time=%.16e element=%s direction=%zu extent=%zu "
            "metric=(%.16e,%.16e,%.16e,%.16e) "
            "derivative=(%.16e,%.16e,%.16e,%.16e) "
            "constraints=(%.16e,%.16e) before=%s after=%s "
            "target_shift=%.16e effective_target=%.16e\n",
            time, get_output(id), d, mesh.extents(d), g.error, g.next_error,
            g.convergence, g.pile_up, k.error, k.next_error, k.convergence,
            k.pile_up, diagnostics.constraints[0][d],
            diagnostics.constraints[1][d], get_output(before[d]),
            get_output(result[d]), history->log10_target_shift[d],
            settings.spectral_target *
                pow(10.0, history->log10_target_shift[d]));
      }
    }
    return result;
  }
  ChalmersWithConstraintFeedback<3> refinement_{};
  spec_amr::Settings settings_{};
  bool enabled_{false};
  double start_time_{0.0};
  bool verbose_{false};
};
}  // namespace gh
