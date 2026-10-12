// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <pup.h>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "DataStructures/DataBox/ValidateSelection.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/Chalmers.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Criterion.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Type.hpp"
#include "ParallelAlgorithms/Amr/CriterionData.hpp"
#include "ParallelAlgorithms/Amr/Policies/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

namespace amr::Criteria {
/*!
 * \brief Experimental hp criterion with optional hp-coarsening based on modal
 * decay.
 *
 * Uses the estimator in `::Chalmers::estimate`. Above the absolute plus
 * relative target, chooses p when sigma exceeds SmoothnessThreshold and h
 * otherwise. Split takes precedence if different monitored components disagree.
 * PRefineOnly provides a control experiment using the same estimator and
 * adaptation phase. An above-target fit with too
 * few meaningful modes stops the run with a diagnostic rather than silently
 * choosing h or p.
 * Components with sqrt(2^Dim) max|u| <= AbsoluteTarget are ignored before
 * fitting, so a nondecaying fit to a negligible far tail cannot trigger
 * unbounded refinement. This bounds the reference L2 norm of the sampled
 * Legendre interpolant, not the error in an unsampled analytic function.
 * With UseDirectionalEnergy in multiple dimensions, uses
 * ::Chalmers::estimate_directional_energy to fit
 * the transverse L2 norm of the directional modal coefficients. Both the
 * error and smoothness come from this same sequence. Directions with sampled
 * nonconstant reference L2 norm below the absolute plus relative target are
 * ignored before fitting. These are experimental extensions of the stripe
 * rule; the 1D behavior is unchanged.
 *
 * A positive CoarseningFactor enables p reduction when every monitored
 * component's candidate error is below factor * target / Dim. The candidate
 * includes the discarded modal norm and a fit of the retained spectrum.
 * Division by Dim budgets simultaneous directional reductions conservatively.
 * A factor strictly below one provides hysteresis. Zero retains refinement-only
 * behavior. HCoarsening additionally enables one-axis joins of exact face
 * siblings after checking projection loss and the prospective parent's
 * spectrum. PrepareRefinementCriteria exchanges the required fields, rotating
 * the candidate axis each AMR round. A join retains the maximum child extents
 * and does not simultaneously decrease polynomial order. Sibling agreement,
 * refinement limits and 2:1 balance are enforced by the existing AMR actions.
 *
 * IgnoreHighestMode registers the name ChalmersIgnoringHighestMode. Its fits
 * omit the highest mode in each direction, appropriate for a steep filter
 * that zeros that mode. Coarsening also budgets the next mode that filtering
 * would discard at the lower resolution. Other filter transfer functions
 * are not modeled by this approximation.
 *
 * SeparateParity registers ChalmersParity[IgnoringHighestMode] and uses
 * independent even/odd fits in refinement and both coarsening checks. Requires
 * directional energy in multiple dimensions. The legacy names retain their
 * original estimator for control experiments.
 *
 * The joint decision runs in EvaluateAmrCriteria / AdjustDomain, not
 * RefineMesh. Only tensor-product Legendre meshes with real-valued fields are
 * supported.
 */
template <size_t Dim, typename TensorTags, bool IgnoreHighestMode = false,
          bool SeparateParity = false>
class Chalmers : public Criterion {
 public:
  struct VariablesToMonitor {
    using type = std::vector<std::string>;
    static size_t lower_bound_on_size() { return 1; }
    static constexpr Options::String help = "Tensors whose modal decay to fit.";
  };
  struct AbsoluteTarget {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help = "Absolute spectral error target.";
  };
  struct RelativeTarget {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help = "Error target relative to max |u|.";
  };
  struct SmoothnessThreshold {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help =
        "Choose p when the natural-log decay rate exceeds this positive value.";
  };
  struct PRefineOnly {
    using type = bool;
    static constexpr Options::String help =
        "Always choose p above target, for a comparison with the same "
        "estimator.";
  };
  struct UseDirectionalEnergy {
    using type = bool;
    static constexpr Options::String help =
        "Fit the directional L2 spectrum. False retains individual stripe "
        "fits for control experiments. Available only in 2D and 3D.";
  };
  struct CoarseningFactor {
    using type = double;
    static type lower_bound() { return 0.0; }
    static type upper_bound() { return 1.0; }
    static constexpr Options::String help =
        "P-coarsen if the candidate error is below factor * target / Dim. "
        "Use a factor below one for hysteresis; zero disables coarsening.";
  };
  struct HCoarsening {
    using type = bool;
    static constexpr Options::String help =
        "Check the projected parent before joining exact face siblings. "
        "Requires positive CoarseningFactor, hp mode, directional energy, "
        "and anisotropic AMR policies in multiple dimensions.";
  };
  using options =
      tmpl::append<tmpl::list<VariablesToMonitor, AbsoluteTarget,
                              RelativeTarget, SmoothnessThreshold, PRefineOnly>,
                   tmpl::conditional_t<Dim == 1, tmpl::list<>,
                                       tmpl::list<UseDirectionalEnergy>>,
                   tmpl::list<CoarseningFactor, HCoarsening>>;
  static constexpr Options::String help =
      "Chalmers hp experiment with optional hp-coarsening on Legendre "
      "elements.";

  Chalmers() = default;
  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const Options::Context& context = {})
      : Chalmers(std::move(variables), absolute_target, relative_target,
                 smoothness_threshold, p_refine_only, true, 0.0, context) {}

  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const bool use_directional_energy,
           const Options::Context& context = {})
      : Chalmers(std::move(variables), absolute_target, relative_target,
                 smoothness_threshold, p_refine_only, use_directional_energy,
                 0.0, context) {}

  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const double coarsening_factor,
           const Options::Context& context = {})
      : Chalmers(std::move(variables), absolute_target, relative_target,
                 smoothness_threshold, p_refine_only, true, coarsening_factor,
                 context) {}

  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const bool use_directional_energy,
           const double coarsening_factor, const Options::Context& context = {})
      : Chalmers(std::move(variables), absolute_target, relative_target,
                 smoothness_threshold, p_refine_only, use_directional_energy,
                 coarsening_factor, false, context) {}

  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const double coarsening_factor,
           const bool h_coarsening, const Options::Context& context = {})
      : Chalmers(std::move(variables), absolute_target, relative_target,
                 smoothness_threshold, p_refine_only, true, coarsening_factor,
                 h_coarsening, context) {}

  Chalmers(std::vector<std::string> variables, const double absolute_target,
           const double relative_target, const double smoothness_threshold,
           const bool p_refine_only, const bool use_directional_energy,
           const double coarsening_factor, const bool h_coarsening,
           const Options::Context& context = {})
      : variables_(std::move(variables)),
        absolute_target_(absolute_target),
        relative_target_(relative_target),
        smoothness_threshold_(smoothness_threshold),
        p_refine_only_(p_refine_only),
        use_directional_energy_(use_directional_energy),
        coarsening_factor_(coarsening_factor),
        h_coarsening_(h_coarsening) {
    db::validate_selection<TensorTags>(variables_, context);
    if (SeparateParity and Dim > 1 and not use_directional_energy_) {
      PARSE_ERROR(context, "Separate-parity fits require directional energy.");
    }
    if (h_coarsening_ and (coarsening_factor_ <= 0.0 or p_refine_only_ or
                           (Dim > 1 and not use_directional_energy_))) {
      PARSE_ERROR(context,
                  "HCoarsening requires positive CoarseningFactor, hp mode, "
                  "and directional-energy estimates.");
    }
    if (not std::isfinite(coarsening_factor_) or coarsening_factor_ < 0.0 or
        coarsening_factor_ >= 1.0) {
      PARSE_ERROR(context, "CoarseningFactor must be finite and in [0, 1).");
    }
    if (not std::isfinite(smoothness_threshold_) or
        smoothness_threshold_ <= 0.0) {
      PARSE_ERROR(context, "SmoothnessThreshold must be finite and positive.");
    }
    if (not std::isfinite(absolute_target_) or absolute_target_ < 0.0 or
        not std::isfinite(relative_target_) or relative_target_ < 0.0 or
        (absolute_target_ == 0.0 and relative_target_ == 0.0)) {
      PARSE_ERROR(context,
                  "Set finite nonnegative targets, at least one positive.");
    }
  }

  explicit Chalmers(CkMigrateMessage* msg) : Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Chalmers);  // NOLINT

  static std::string name() {
    if constexpr (SeparateParity) {
      return IgnoreHighestMode ? "ChalmersParityIgnoringHighestMode"
                               : "ChalmersParity";
    }
    return IgnoreHighestMode ? "ChalmersIgnoringHighestMode" : "Chalmers";
  }
  Type type() override { return Type::hp; }
  std::string observation_name() override { return name(); }
  bool needs_sibling_data() const override { return h_coarsening_; }

  template <typename DbTags>
  CriterionTensorData pack_sibling_data(const db::DataBox<DbTags>& box) const {
    CriterionTensorData result;
    tmpl::for_each<TensorTags>([&](const auto tag_v) {
      using tag = tmpl::type_from<std::decay_t<decltype(tag_v)>>;
      const auto name = db::tag_name<tag>();
      if (alg::found(variables_, name)) {
        auto& components = result[name];
        for (const auto& component : db::get<tag>(box)) {
          components.push_back(component);
        }
      }
    });
    return result;
  }
  using compute_tags_for_observation_box = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTags, typename Metavariables>
  std::array<Flag, Dim> operator()(const db::DataBox<DbTags>& box,
                                   Parallel::GlobalCache<Metavariables>& cache,
                                   const ElementId<Dim>& element_id) const {
    return evaluate_with_target_factors(box, cache, element_id,
                                        make_array<Dim>(1.0));
  }

  /// Apply external directional tolerance factors in (0, 1]. The smallest
  /// factor sets conservative common fit-negligibility and join budgets;
  /// directional refinement and p-coarsening use their own factors.
  template <typename DbTags, typename Metavariables>
  std::array<Flag, Dim> evaluate_with_target_factors(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& element_id,
      const std::array<double, Dim>& target_factors) const {
    for (const double factor : target_factors) {
      ASSERT(std::isfinite(factor) and factor > 0.0 and factor <= 1.0,
             "Chalmers target factors must lie in (0, 1]");
    }
    const double minimum_factor =
        *std::min_element(target_factors.begin(), target_factors.end());
    auto result = make_array<Dim>(
        coarsening_factor_ > 0.0 ? Flag::DecreaseResolution : Flag::DoNothing);
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    tmpl::for_each<TensorTags>([&](const auto tag_v) {
      using tag = tmpl::type_from<std::decay_t<decltype(tag_v)>>;
      if (not alg::found(variables_, db::tag_name<tag>())) {
        return;
      }
      for (const auto& component : db::get<tag>(box)) {
        const double amplitude = max(abs(component));
        const double base_target =
            absolute_target_ + relative_target_ * amplitude;
        const double target = minimum_factor * base_target;
        const double coarsening_budget = coarsening_factor_ * target / Dim;
        const bool negligible =
            std::sqrt(static_cast<double>(size_t{1} << Dim)) * amplitude <=
            minimum_factor * absolute_target_;
        if (negligible and coarsening_factor_ == 0.0) {
          continue;
        }
        const auto coarsening = coarsening_factor_ > 0.0
                                    ? ::Chalmers::p_coarsening_error(
                                          component, mesh, coarsening_budget,
                                          IgnoreHighestMode, SeparateParity)
                                    : std::array<double, Dim>{};
        const auto estimates = [&]() {
          if (negligible) {
            return std::array<::Chalmers::DirectionEstimate, Dim>{};
          }
          if constexpr (SeparateParity) {
            return ::Chalmers::estimate_directional_energy(
                component, mesh, target, IgnoreHighestMode, true);
          } else if constexpr (Dim == 1) {
            return ::Chalmers::estimate(component, mesh, 0.0,
                                        IgnoreHighestMode);
          } else {
            return use_directional_energy_
                       ? ::Chalmers::estimate_directional_energy(
                             component, mesh, target, IgnoreHighestMode)
                       : ::Chalmers::estimate(component, mesh,
                                              minimum_factor * absolute_target_,
                                              IgnoreHighestMode);
          }
        }();
        for (size_t d = 0; d < Dim; ++d) {
          const auto& fit = gsl::at(estimates, d);
          const double directional_target = target_factors[d] * base_target;
          if (fit.error <= directional_target) {
            const auto decision =
                coarsening_factor_ > 0.0 and
                        gsl::at(coarsening, d) <=
                            coarsening_factor_ * directional_target / Dim
                    ? Flag::DecreaseResolution
                    : Flag::DoNothing;
            gsl::at(result, d) = std::max(gsl::at(result, d), decision);
            continue;
          }
          if (not std::isfinite(fit.sigma)) {
            ERROR("Chalmers fit is undetermined in element "
                  << element_id << ", field " << db::tag_name<tag>()
                  << ", direction " << d
                  << ". Increase the initial resolution; no h/p fallback was "
                     "selected.");
          }
          const auto decision =
              p_refine_only_ or fit.sigma > smoothness_threshold_
                  ? Flag::IncreaseResolution
                  : Flag::Split;
          gsl::at(result, d) = std::max(gsl::at(result, d), decision);
        }
      }
    });
    if (h_coarsening_ and not alg::any_of(result, [](const Flag flag) {
          return flag == Flag::Split or flag == Flag::IncreaseResolution;
        })) {
      if constexpr (db::tag_is_retrievable_v<
                        amr::Tags::CriterionNeighborData<Dim>,
                        db::DataBox<DbTags>>) {
        const auto& state = db::get<amr::Tags::CriterionNeighborData<Dim>>(box);
        const auto& policies = db::get<amr::Tags::Policies>(box);
        if (not state.sibling.has_value() or state.round == 0 or
            not policies.allow_coarsening() or
            (Dim > 1 and policies.isotropy() != amr::Isotropy::Anisotropic)) {
          return result;
        }
        const size_t axis = (state.round - 1) % Dim;
        if (element_id.segment_id(axis).refinement_level() <=
            policies.limits().minimum_refinement_level()) {
          return result;
        }
        const auto& sibling = state.sibling->second;
        ASSERT(sibling.round == state.round, "Stale Chalmers sibling data");
        const bool is_lower = element_id.segment_id(axis).index() % 2 == 0;
        bool accept = true;
        tmpl::for_each<TensorTags>([&](const auto tag_v) {
          using tag = tmpl::type_from<std::decay_t<decltype(tag_v)>>;
          if (not alg::found(variables_, db::tag_name<tag>())) {
            return;
          }
          const auto& other = sibling.tensors.at(db::tag_name<tag>());
          const auto& tensor = db::get<tag>(box);
          ASSERT(other.size() == tensor.size(), "Sibling tensor sizes differ");
          size_t i = 0;
          for (const auto& component : tensor) {
            const auto candidate = ::Chalmers::h_coarsening_error(
                is_lower ? component : other[i], is_lower ? mesh : sibling.mesh,
                is_lower ? other[i] : component, is_lower ? sibling.mesh : mesh,
                axis, coarsening_factor_ * minimum_factor * absolute_target_,
                IgnoreHighestMode, SeparateParity);
            const double target = absolute_target_ +
                                  relative_target_ * candidate.parent_amplitude;
            accept = accept and candidate.error <= coarsening_factor_ *
                                                       minimum_factor * target;
            ++i;
          }
        });
        if (accept) {
          result.fill(Flag::DoNothing);
          result[axis] = Flag::Join;
        }
      } else {
        ERROR("HCoarsening requires prepared sibling data in the AMR DataBox.");
      }
    }
    return result;
  }

  void pup(PUP::er& p) override {
    Criterion::pup(p);
    p | variables_;
    p | absolute_target_;
    p | relative_target_;
    p | smoothness_threshold_;
    p | p_refine_only_;
    p | use_directional_energy_;
    p | coarsening_factor_;
    p | h_coarsening_;
  }

 private:
  std::vector<std::string> variables_{};
  double absolute_target_{0.0};
  double relative_target_{0.0};
  double smoothness_threshold_{1.0};
  bool p_refine_only_{false};
  bool use_directional_energy_{true};
  double coarsening_factor_{0.0};
  bool h_coarsening_{false};
};

template <size_t Dim, typename TensorTags, bool IgnoreHighestMode,
          bool SeparateParity>
PUP::able::PUP_ID
    Chalmers<Dim, TensorTags, IgnoreHighestMode, SeparateParity>::my_PUP_ID =
        0;  // NOLINT
}  // namespace amr::Criteria
