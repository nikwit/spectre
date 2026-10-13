// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/ParseError.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/MakeArray.hpp"

namespace gh {
/*!
 * \brief SpEC-style spectral p refinement and coarsening for GH.
 *
 * Both decisions use relative, grouped metric and Pi/Phi spectra and projected
 * two- and three-index constraints. The local stored target adjustment,
 * prospective-removal checks, pile-up override and cooldown are shared with
 * the experimental hybrid controller. No Chalmers estimate or h decision is
 * involved. Neighbor smoothing and adaptive filtering are not implemented.
 *
 * Requires tensor-product Legendre meshes and a steep filter that removes
 * only the highest mode in each direction. That mode is excluded from the
 * diagnostics. StartTime holds all decisions and history updates until the
 * specified time; this global gate is not a radius-dependent junk gate.
 * InitializeHistory and ProjectHistory must be installed in the executable.
 * Schedule EvaluateAmrCriteria/AdjustDomain phases; this criterion opts out of
 * the RefineMesh event so that changes update neighbor meshes as well.
 */
class SpecSpectral : public amr::Criterion {
 public:
  struct RelativeTruncationTarget {
    using type = double;
    static constexpr Options::String help = "Relative tensor spectral target.";
  };
  struct ConstraintTargets {
    using type = std::array<double, 2>;
    static constexpr Options::String help =
        "Projected nodal RMS targets for C_ia and C_iab, respectively.";
  };
  struct ForceDropPileUpModes {
    using type = bool;
    static constexpr Options::String help =
        "Apply SpEC's exceptional pile-up removal and p-addition hysteresis.";
  };
  struct MinimumModes {
    using type = size_t;
    static type lower_bound() { return 5; }
    static constexpr Options::String help =
        "Minimum points per direction, including the filtered highest mode.";
  };
  struct CallsBetweenChanges {
    using type = size_t;
    static type lower_bound() { return 1; }
    static constexpr Options::String help =
        "AMR-call cooldown after p changes.";
  };
  struct TruncationErrorAdjustmentLimit {
    using type = double;
    static type upper_bound() { return 0.0; }
    static constexpr Options::String help =
        "Minimum stored log10 target shift. Inspiral uses -1; zero disables "
        "target adjustment. No neighbor smoothing is applied.";
  };
  struct StartTime {
    using type = double;
    static constexpr Options::String help =
        "Hold all p decisions and history updates before this time.";
  };
  struct Verbose {
    using type = bool;
    static constexpr Options::String help = "Record pre-AMR diagnostics.";
  };
  using options =
      tmpl::list<RelativeTruncationTarget, ConstraintTargets,
                 ForceDropPileUpModes, MinimumModes, CallsBetweenChanges,
                 TruncationErrorAdjustmentLimit, StartTime, Verbose>;
  static constexpr Options::String help =
      "SpEC-style GH p refinement and coarsening on filtered Legendre meshes.";

  SpecSpectral() = default;
  SpecSpectral(double spectral_target,
               const std::array<double, 2>& constraint_targets, bool force,
               size_t minimum_modes, size_t cooldown, double adjustment_limit,
               double start_time, bool verbose,
               const Options::Context& context = {});
  explicit SpecSpectral(CkMigrateMessage* msg) : amr::Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(SpecSpectral);  // NOLINT
  amr::Criteria::Type type() override { return amr::Criteria::Type::p; }
  bool p_refine_in_domain_adjustment_phase() const override { return true; }
  std::string observation_name() override { return "SpecSpectral"; }
  using compute_tags_for_observation_box = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, 3> operator()(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<3>& id) const {
    auto history = db::get<spec_amr::HistoryTag>(box);
    return evaluate_impl(make_not_null(&history), box, id, false);
  }

  template <typename ComputeTags, typename DbType, typename Metavariables>
  std::array<amr::Flag, 3> evaluate_with_history(
      const gsl::not_null<ObservationBox<ComputeTags, DbType>*> box,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<3>& id) const {
    auto history = db::get<spec_amr::HistoryTag>(box->databox());
    const auto flags =
        evaluate_impl(make_not_null(&history), box->databox(), id, true);
    db::mutate<spec_amr::HistoryTag>(
        [&history](const auto h) { *h = std::move(history); },
        make_not_null(&box->databox()));
    box->reset();
    return flags;
  }
  void pup(PUP::er& p) override;

 private:
  template <typename DbTags>
  std::array<amr::Flag, 3> evaluate_impl(
      const gsl::not_null<spec_amr::History*> history,
      const db::DataBox<DbTags>& box, const ElementId<3>& id,
      const bool log) const {
    const double time = db::get<::Tags::Time>(box);
    if (time < start_time_) {
      return make_array<3>(amr::Flag::DoNothing);
    }
    const auto& mesh = db::get<domain::Tags::Mesh<3>>(box);
    for (size_t d = 0; d < 3; ++d) {
      if (mesh.extents(d) < 5 or mesh.basis(d) != Spectral::Basis::Legendre) {
        ERROR("SpecSpectral requires Legendre extents >= 5 in " << id);
      }
    }
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
    const auto result = spec_amr::decide_p(
        history, diagnostics, mesh.extents().indices(), time, settings_);
    if (verbose_ and log) {
      for (size_t d = 0; d < 3; ++d) {
        const auto& g = diagnostics.spectra[0][d];
        const auto& k = diagnostics.spectra[1][d];
        Parallel::printf(
            "SpecSpectral time=%.16e element=%s direction=%zu extent=%zu "
            "metric=(%.16e,%.16e,%.16e,%.16e) "
            "derivative=(%.16e,%.16e,%.16e,%.16e) "
            "constraints=(%.16e,%.16e) after=%s "
            "target_shift=%.16e effective_target=%.16e\n",
            time, get_output(id), d, mesh.extents(d), g.error, g.next_error,
            g.convergence, g.pile_up, k.error, k.next_error, k.convergence,
            k.pile_up, diagnostics.constraints[0][d],
            diagnostics.constraints[1][d], get_output(result[d]),
            history->log10_target_shift[d],
            settings_.spectral_target *
                pow(10.0, history->log10_target_shift[d]));
      }
    }
    return result;
  }
  spec_amr::Settings settings_{};
  double start_time_{0.0};
  bool verbose_{false};
};
}  // namespace gh
