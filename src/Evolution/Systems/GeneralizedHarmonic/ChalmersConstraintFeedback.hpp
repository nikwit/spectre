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

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/ParseError.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Chalmers.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/GetOutput.hpp"

namespace gh {
/*!
 * \brief Nodal RMS of C_iab = d_i g_ab - Phi_iab along logical directions.
 *
 * The derivative index is contracted with each normalized column of dx/dxi.
 * The remaining indices are summed with symmetric-tensor multiplicities.
 * This is a Euclidean component norm, not a spacetime invariant or a
 * physical-volume integral. It uses the current discrete fields, before any
 * proposed projection or filtering. Unlike the existing Constraints criterion,
 * the sum includes both copies of each off-diagonal spacetime component.
 */
template <size_t Dim>
std::array<double, Dim> chalmers_directional_constraints(
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian) {
  const auto derivative = partial_derivative(metric, mesh, inverse_jacobian);
  const auto jacobian = determinant_and_inverse(inverse_jacobian).second;
  std::array<double, Dim> result{};
  const size_t points = mesh.number_of_grid_points();
  for (size_t d = 0; d < Dim; ++d) {
    DataVector norm_squared(points, 0.0);
    for (size_t i = 0; i < Dim; ++i) {
      norm_squared += square(jacobian.get(i, d));
    }
    DataVector constraint_squared(points, 0.0);
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        DataVector projected(points, 0.0);
        for (size_t i = 0; i < Dim; ++i) {
          projected +=
              jacobian.get(i, d) * (derivative.get(i, a, b) - phi.get(i, a, b));
        }
        constraint_squared += (a == b ? 1.0 : 2.0) * square(projected);
      }
    }
    result[d] = std::sqrt(sum(constraint_squared / norm_squared) /
                          static_cast<double>(points));
  }
  return result;
}

/*!
 * \brief Experimental GH constraint feedback around parity-aware Chalmers AMR.
 *
 * Each direction uses target factor max(MinimumTargetFactor, min(1, T/C_d)).
 * MinimumTargetFactor=1 disables tightening. This first experiment is an
 * instantaneous, bounded adjustment, not SpEC's history-dependent controller.
 * If EnableCoarseningVeto, C_d > T vetoes p reduction in direction d. A join
 * is vetoed if any direction is above target, because parent projection and
 * filtering can affect every derivative direction. Both siblings evaluate
 * their own constraints and must agree to join through the usual AMR actions.
 * Refinement flags are preserved. No prospective-parent constraint guarantee
 * is made. The field criterion's projection/filter-loss checks still apply.
 */
template <size_t Dim>
class ChalmersWithConstraintFeedback : public amr::Criterion {
 public:
  using FieldCriterion = amr::Criteria::Chalmers<
      Dim, typename gh::System<Dim>::variables_tag::tags_list, true, true>;
  struct Chalmers {
    using type = FieldCriterion;
    static constexpr Options::String help =
        "Parity-aware filtered hp criterion.";
  };
  struct ConstraintTarget {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help =
        "Target directional nodal RMS C_iab.";
  };
  struct MinimumTargetFactor {
    using type = double;
    static type lower_bound() { return 0.0; }
    static type upper_bound() { return 1.0; }
    static constexpr Options::String help =
        "Lower bound on spectral target factor; 1 disables tightening.";
  };
  struct EnableCoarseningVeto {
    using type = bool;
    static constexpr Options::String help =
        "Veto coarsening where current constraints exceed the target.";
  };
  struct Verbose {
    using type = bool;
    static constexpr Options::String help =
        "Log directional constraints, factors and flags at each evaluation.";
  };
  using options = tmpl::list<Chalmers, ConstraintTarget, MinimumTargetFactor,
                             EnableCoarseningVeto, Verbose>;
  static constexpr Options::String help =
      "Experimental GH constraint feedback for Chalmers hp AMR.";

  ChalmersWithConstraintFeedback() = default;
  ChalmersWithConstraintFeedback(FieldCriterion criterion, const double target,
                                 const double minimum_factor, const bool veto,
                                 const bool verbose,
                                 const Options::Context& context = {})
      : criterion_(std::move(criterion)),
        target_(target),
        minimum_factor_(minimum_factor),
        veto_(veto),
        verbose_(verbose) {
    if (not std::isfinite(target_) or target_ <= 0.0) {
      PARSE_ERROR(context, "ConstraintTarget must be finite and positive.");
    }
    if (not std::isfinite(minimum_factor_) or minimum_factor_ <= 0.0 or
        minimum_factor_ > 1.0) {
      PARSE_ERROR(context, "MinimumTargetFactor must lie in (0, 1].");
    }
  }
  explicit ChalmersWithConstraintFeedback(CkMigrateMessage* msg)
      : amr::Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ChalmersWithConstraintFeedback);  // NOLINT
  static std::string name() { return "ChalmersWithConstraintFeedback"; }
  amr::Criteria::Type type() override { return amr::Criteria::Type::hp; }
  std::string observation_name() override { return name(); }
  bool needs_sibling_data() const override {
    return criterion_.needs_sibling_data();
  }
  template <typename DbTags>
  amr::CriterionTensorData pack_sibling_data(
      const db::DataBox<DbTags>& box) const {
    return criterion_.pack_sibling_data(box);
  }
  using compute_tags_for_observation_box = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::DataBox>;

  /// Evaluate the field criterion with externally supplied target factors.
  /// This replaces this wrapper's instantaneous constraint feedback and veto;
  /// the external controller is responsible for its constraint checks.
  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, Dim> evaluate_with_target_factors(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id,
      const std::array<double, Dim>& factors) const {
    return criterion_.evaluate_with_target_factors(box, cache, element_id,
                                                   factors);
  }

  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, Dim> operator()(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id) const {
    const auto constraints = chalmers_directional_constraints(
        db::get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(box),
        db::get<gh::Tags::Phi<DataVector, Dim>>(box),
        db::get<domain::Tags::Mesh<Dim>>(box),
        db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                              Frame::Inertial>>(box));
    std::array<double, Dim> factors{};
    bool any_excess = false;
    for (size_t d = 0; d < Dim; ++d) {
      if (not std::isfinite(constraints[d])) {
        ERROR("Nonfinite Chalmers constraint in " << element_id);
      }
      any_excess |= constraints[d] > target_;
      factors[d] = constraints[d] <= target_
                       ? 1.0
                       : std::max(minimum_factor_, target_ / constraints[d]);
    }
    auto flags = criterion_.evaluate_with_target_factors(box, cache, element_id,
                                                         factors);
    const auto before_veto = flags;
    if (veto_) {
      for (size_t d = 0; d < Dim; ++d) {
        if ((flags[d] == amr::Flag::DecreaseResolution and
             constraints[d] > target_) or
            (flags[d] == amr::Flag::Join and any_excess)) {
          flags[d] = amr::Flag::DoNothing;
        }
      }
    }
    if (verbose_) {
      Parallel::printf(
          "ChalmersFeedback time=%.16e element=%s constraints=%s factors=%s "
          "before=%s after=%s\n",
          db::get<::Tags::Time>(box), get_output(element_id),
          get_output(constraints), get_output(factors), get_output(before_veto),
          get_output(flags));
    }
    return flags;
  }
  void pup(PUP::er& p) override {
    amr::Criterion::pup(p);
    p | criterion_;
    p | target_;
    p | minimum_factor_;
    p | veto_;
    p | verbose_;
  }

 private:
  FieldCriterion criterion_{};
  double target_{0.0};
  double minimum_factor_{1.0};
  bool veto_{false};
  bool verbose_{false};
};

template <size_t Dim>
PUP::able::PUP_ID ChalmersWithConstraintFeedback<Dim>::my_PUP_ID = 0;  // NOLINT
}  // namespace gh
