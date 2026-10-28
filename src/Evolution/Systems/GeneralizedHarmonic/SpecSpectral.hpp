// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Tag.hpp"
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
 * Supports polynomial meshes and polynomial x S2 shells. Diagnostics use
 * the actual element-local volume filter at the current and prospective lower
 * resolution. S2 has one coupled angular decision. Constraint targets use
 * SpEC's topological aggregation of projected nodal norms. Scalar StartTime
 * holds decisions and history updates until the specified time. Its
 * RadiusDependentStart alternative additionally freezes p changes in selected
 * outer blocks until max(grid-frame radius)+RadiusPlusTime. The radial gate
 * retains diagnostic/history updates and is bypassed during initial p
 * adjustment, as is the global start-time gate. InitializeHistory and
 * ProjectHistory must be installed in the executable. Schedule
 * EvaluateAmrCriteria/AdjustDomain phases; this criterion opts out of the
 * RefineMesh event so that changes update neighbor meshes as well.
 */
class SpecSpectral : public amr::Criterion {
 public:
  struct GaussianTarget {
    struct Center {
      using type = std::array<double, 3>;
      static constexpr Options::String help = "Center in grid coordinates.";
    };
    struct Width {
      using type = double;
      static constexpr Options::String help = "Gaussian e-folding radius.";
    };
    struct Amplitude {
      using type = double;
      static constexpr Options::String help = "Nonnegative blending weight.";
    };
    struct Target {
      using type = double;
      static constexpr Options::String help = "Relative spectral target.";
    };
    using options = tmpl::list<Center, Width, Amplitude, Target>;
    static constexpr Options::String help =
        "Gaussian contribution to the log10 spectral target.";
    GaussianTarget() = default;
    GaussianTarget(std::array<double, 3> center_in, double width_in,
                   double amplitude_in, double target_in,
                   const Options::Context& context = {});
    void pup(PUP::er& p);
    std::array<double, 3> center{};
    double width{1.0};
    double amplitude{1.0};
    double target{1.e-4};
  };
  struct SpatialTarget {
    struct Background {
      using type = double;
      static constexpr Options::String help = "Constant target of weight one.";
    };
    struct Gaussians {
      using type = std::vector<GaussianTarget>;
      static constexpr Options::String help = "Grid-frame Gaussian profiles.";
    };
    using options = tmpl::list<Background, Gaussians>;
    static constexpr Options::String help =
        "SpEC log10 target blend, using each profile's maximum over the "
        "element's initial grid points. Frozen across p changes.";
    SpatialTarget();
    SpatialTarget(double background_in,
                  std::vector<GaussianTarget> gaussians_in,
                  const Options::Context& context = {});
    double operator()(
        const tnsr::I<DataVector, 3, Frame::Grid>& coordinates) const;
    void pup(PUP::er& p);
    double background{1.e-4};
    std::vector<GaussianTarget> gaussians{};
  };
  struct RelativeTruncationTarget {
    using type = std::variant<double, SpatialTarget>;
    static constexpr Options::String help =
        "Constant relative spectral target or a spatial log10 target blend.";
  };
  struct ConstraintTargets {
    using type = std::array<double, 2>;
    static constexpr Options::String help =
        "SpEC projected constraint targets for C_ia and C_iab, respectively.";
  };
  struct ForceDropPileUpModes {
    using type = bool;
    static constexpr Options::String help =
        "Apply SpEC's exceptional pile-up removal and p-addition hysteresis.";
  };
  struct MinimumModes {
    using type = std::variant<size_t, std::array<size_t, 3>,
                              std::vector<std::array<size_t, 3>>>;
    static constexpr Options::String help =
        "Minimum points including filtered modes: a scalar, three "
        "logical-direction values, or one triple per block in domain order. "
        "Each must be >= 5. On S2, theta and phi bounds give a single "
        "l_max floor. Initialize the mesh at or above these bounds.";
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
  struct RadiusDependentStart {
    struct Time {
      using type = double;
      static constexpr Options::String help = "Global AMR start time.";
    };
    struct Blocks {
      using type = std::vector<std::string>;
      static constexpr Options::String help =
          "Outer block names or groups subject to the junk timing restriction.";
    };
    struct RadiusPlusTime {
      using type = double;
      static constexpr Options::String help =
          "Release p changes at max(initial grid-frame radius) plus this time. "
          "For a start at zero, SpEC uses max(200,3*separation) for SKS/SSphKS "
          "or 2*separation for SHK bound-orbit data.";
    };
    using options = tmpl::list<Time, Blocks, RadiusPlusTime>;
    static constexpr Options::String help =
        "Global start time plus a radius-dependent restriction in outer "
        "blocks.";
    RadiusDependentStart();
    RadiusDependentStart(double time_in, std::vector<std::string> blocks_in,
                         double radius_plus_time_in,
                         const Options::Context& context = {});
    double release_time(
        const Domain<3>& domain, size_t block_id,
        const tnsr::I<DataVector, 3, Frame::Grid>& coordinates) const;
    void pup(PUP::er& p);
    double time{0.0};
    std::vector<std::string> blocks{};
    double radius_plus_time{0.0};
  };
  struct StartTime {
    using type = std::variant<double, RadiusDependentStart>;
    static constexpr Options::String help =
        "Global start time, or Time/Blocks/RadiusPlusTime for a junk gate. "
        "Startup p adjustment bypasses both restrictions.";
  };
  struct Verbose {
    using type = bool;
    static constexpr Options::String help = "Record pre-AMR diagnostics.";
  };
  struct TimeBetweenChanges {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help =
        "Minimum elapsed evolution time after a realized p change. Zero "
        "disables this restriction; the call-count cooldown still applies.";
  };
  using options =
      tmpl::list<RelativeTruncationTarget, ConstraintTargets,
                 ForceDropPileUpModes, MinimumModes, CallsBetweenChanges,
                 TruncationErrorAdjustmentLimit, StartTime, Verbose,
                 TimeBetweenChanges>;
  static constexpr Options::String help =
      "SpEC-style GH p refinement on polynomial and spherical-harmonic meshes.";

  SpecSpectral() = default;
  SpecSpectral(RelativeTruncationTarget::type spectral_target,
               const std::array<double, 2>& constraint_targets, bool force,
               MinimumModes::type minimum_modes, size_t cooldown,
               double adjustment_limit, StartTime::type start_time,
               bool verbose, double time_between_changes = 0.0,
               const Options::Context& context = {});
  SpecSpectral(RelativeTruncationTarget::type spectral_target,
               const std::array<double, 2>& constraint_targets, bool force,
               size_t minimum_modes, size_t cooldown, double adjustment_limit,
               StartTime::type start_time, bool verbose,
               double time_between_changes = 0.0,
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
    if (not history->initializing and time < start_time_) {
      return make_array<3>(amr::Flag::DoNothing);
    }
    const auto& mesh = db::get<domain::Tags::Mesh<3>>(box);
    const spec_amr::RefinementVariables variables{mesh};
    const auto modes = variables.mode_counts(mesh);
    auto minimum = minimum_modes(id.block_id());
    if (variables.size() == 2) {
      minimum[1] = std::max(minimum[1], (minimum[2] + 2) / 2);
      minimum[2] = 0;
    }
    for (size_t v = 0; v < variables.size(); ++v) {
      if (modes[v] < 5) {
        ERROR("SpecSpectral requires independent extents >= 5 in " << id);
      }
      if (not history->initializing and modes[v] < minimum[v]) {
        ERROR("SpecSpectral mesh is below MinimumModes in "
              << id << " refinement variable " << v
              << ". Initialize at or above the floor.");
      }
    }
    if (not history->initializing and history->calls > 0 and
        time == history->time and modes == history->extents) {
      return history->decision;
    }
    spec_amr::Diagnostics diagnostics{};
    using FilterTag = Filters::Tags::SpectralFilter<
        3, ylm::TensorYlm::filter_detail::gh_spacetime_vars_list>;
    using GridJacobian =
        domain::Tags::Jacobian<3, Frame::Grid, Frame::Inertial>;
    if constexpr (db::tag_is_retrievable_v<FilterTag, db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<GridJacobian, db::DataBox<DbTags>>) {
      power_monitor::CartToSphereMatrices unused_polynomial_cache{};
      if (variables.size() == 2 and history->transform_cache == nullptr) {
        history->transform_cache =
            std::make_shared<power_monitor::CartToSphereMatrices>();
      }
      diagnostics = spec_amr::diagnose_spec(
          db::get<gr::Tags::SpacetimeMetric<DataVector, 3>>(box),
          db::get<gh::Tags::Pi<DataVector, 3>>(box),
          db::get<gh::Tags::Phi<DataVector, 3>>(box), mesh,
          db::get<domain::Tags::InverseJacobian<3, Frame::ElementLogical,
                                                Frame::Inertial>>(box),
          db::get<GridJacobian>(box),
          db::get<domain::Tags::Coordinates<3, Frame::Inertial>>(box), time,
          db::get<gauges::Tags::GaugeCondition>(box),
          db::get<gh::Tags::ConstraintGamma2>(box),
          spec_amr::filtered_modes(mesh, db::get<FilterTag>(box)),
          make_not_null(variables.size() == 2 ? history->transform_cache.get()
                                              : &unused_polynomial_cache));
    } else {
      ERROR(
          "SpecSpectral requires the element's SpectralFilter and "
          "Grid-to-Inertial Jacobian tags.");
    }
    auto local_settings = settings_;
    if (spatial_target_.has_value()) {
      if (not std::isfinite(history->spectral_target)) {
        history->spectral_target = (*spatial_target_)(
            db::get<domain::Tags::Coordinates<3, Frame::Grid>>(box));
      }
      local_settings.spectral_target = history->spectral_target;
    }
    auto result =
        spec_amr::decide_p(history, diagnostics,
                           diagnostics.refinement_variables.mode_counts(mesh),
                           time, local_settings);
    if (radius_dependent_start_.has_value()) {
      if (not std::isfinite(history->radius_release_time)) {
        if constexpr (db::tag_is_retrievable_v<domain::Tags::Domain<3>,
                                               db::DataBox<DbTags>>) {
          history->radius_release_time = radius_dependent_start_->release_time(
              db::get<domain::Tags::Domain<3>>(box), id.block_id(),
              db::get<domain::Tags::Coordinates<3, Frame::Grid>>(box));
        } else {
          ERROR("Radius-dependent SpecSpectral start requires the Domain tag.");
        }
      }
      if (not history->initializing and time < history->radius_release_time) {
        result = make_array<3>(amr::Flag::DoNothing);
      }
    }
    // Resolution restrictions also take priority over forced pile-up removal.
    // Store the restricted result so repeated evaluation remains idempotent.
    auto independent_result = variables.independent_flags(result);
    for (size_t v = 0; v < variables.size(); ++v) {
      if (history->initializing and modes[v] < minimum[v]) {
        independent_result[v] = amr::Flag::IncreaseResolution;
      }
      if (independent_result[v] == amr::Flag::DecreaseResolution and
          (modes[v] <= minimum[v] or not diagnostics.filter.can_decrease[v])) {
        independent_result[v] = amr::Flag::DoNothing;
      }
    }
    result = variables.logical_flags(independent_result);
    history->decision = result;
    if (verbose_ and log) {
      for (size_t d = 0; d < variables.size(); ++d) {
        const auto& g = diagnostics.spectra[0][d];
        const auto& k = diagnostics.spectra[1][d];
        Parallel::printf(
            "SpecSpectral time=%.16e element=%s direction=%zu extent=%zu "
            "metric=(%.16e,%.16e,%.16e,%.16e) "
            "derivative=(%.16e,%.16e,%.16e,%.16e) "
            "constraints=(%.16e,%.16e) after=%s "
            "target_shift=%.16e effective_target=%.16e "
            "filtered_modes=%.16e p_down_filtered_modes=%.16e\n",
            time, get_output(id), d, modes[d], g.error, g.next_error,
            g.convergence, g.pile_up, k.error, k.next_error, k.convergence,
            k.pile_up, diagnostics.constraints[0][d],
            diagnostics.constraints[1][d], get_output(independent_result[d]),
            history->log10_target_shift[d],
            local_settings.spectral_target *
                pow(10.0, history->log10_target_shift[d]),
            diagnostics.filter.current[d], diagnostics.filter.after_p_down[d]);
      }
    }
    return result;
  }
  std::array<size_t, 3> minimum_modes(size_t block_id) const;
  spec_amr::Settings settings_{};
  std::optional<SpatialTarget> spatial_target_{};
  MinimumModes::type minimum_modes_{size_t{5}};
  double start_time_{0.0};
  std::optional<RadiusDependentStart> radius_dependent_start_{};
  bool verbose_{false};
};
}  // namespace gh
