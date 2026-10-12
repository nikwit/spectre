// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <pup.h>
#include <string>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/ParseError.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Chalmers.hpp"
#include "Utilities/Algorithm.hpp"

namespace ScalarWave {
/*!
 * \brief Maximum nodal change in the one-index constraint from p coarsening.
 *
 * Removes the last mode in every selected direction, exactly as the native
 * p projector does. Evaluates C_old - C_new = d_i(discarded Psi) - discarded
 * Phi_i on the original grid using its physical inverse Jacobian. This tests
 * all simultaneous reductions together, including mixed modal terms.
 *
 * The result is the maximum nodal Euclidean norm of this change. It is an
 * experimental per-event guard, not a bound on the full accumulated constraint
 * or its continuum maximum between nodes. Existing constraint violation is
 * not repaired or included in the budget. Pi does not enter this constraint.
 */
template <size_t Dim>
double p_coarsening_constraint_change(
    const Scalar<DataVector>& psi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi, const Mesh<Dim>& mesh,
    const std::array<bool, Dim>& decrease,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian) {
  const auto discarded = [&mesh, &decrease](const DataVector& values) {
    auto modes = to_modal_coefficients(values, mesh);
    for (size_t i = 0; i < modes.size(); ++i) {
      size_t index = i;
      bool remove = false;
      for (size_t d = 0; d < Dim; ++d) {
        remove = remove or (decrease[d] and
                            index % mesh.extents(d) == mesh.extents(d) - 1);
        index /= mesh.extents(d);
      }
      if (not remove) {
        modes[i] = 0.0;
      }
    }
    return to_nodal_coefficients(modes, mesh);
  };
  const auto derivative = partial_derivative(
      Scalar<DataVector>{discarded(get(psi))}, mesh, inverse_jacobian);
  DataVector norm_squared(mesh.number_of_grid_points(), 0.0);
  for (size_t i = 0; i < Dim; ++i) {
    const DataVector change = derivative.get(i) - discarded(phi.get(i));
    norm_squared += change * change;
  }
  return std::sqrt(max(norm_squared));
}

/*!
 * \brief Experimental Chalmers rule with a scalar-wave p-coarsening veto.
 *
 * Applies the Chalmers field decision first. If its simultaneous p reductions
 * would inject too much one-index constraint, veto all those reductions while
 * preserving refinement flags. Sibling joining must be disabled: this small
 * prototype validates only p projection, without changing the generic AMR
 * criterion or pretending to certify a prospective h parent.
 */
template <size_t Dim>
class ChalmersWithConstraintGuard : public amr::Criterion {
 public:
  using FieldCriterion = amr::Criteria::Chalmers<
      Dim, tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<Dim, Frame::Inertial>>>;
  struct Chalmers {
    using type = FieldCriterion;
    static constexpr Options::String help = "Underlying field criterion.";
  };
  struct ConstraintChangeTarget {
    using type = double;
    static type lower_bound() { return 0.0; }
    static constexpr Options::String help =
        "Maximum nodal Euclidean change in d_i Psi - Phi_i per p reduction.";
  };
  using options = tmpl::list<Chalmers, ConstraintChangeTarget>;
  static constexpr Options::String help =
      "Veto Chalmers p coarsening that injects excessive scalar-wave "
      "constraint.";

  ChalmersWithConstraintGuard() = default;
  ChalmersWithConstraintGuard(FieldCriterion criterion, const double target,
                              const Options::Context& context = {})
      : criterion_(std::move(criterion)), target_(target) {
    if (not std::isfinite(target_) or target_ <= 0.0) {
      PARSE_ERROR(context,
                  "ConstraintChangeTarget must be finite and positive.");
    }
    if (criterion_.needs_sibling_data()) {
      PARSE_ERROR(context,
                  "The p-constraint prototype requires HCoarsening: false.");
    }
  }
  explicit ChalmersWithConstraintGuard(CkMigrateMessage* msg)
      : amr::Criterion(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ChalmersWithConstraintGuard);  // NOLINT

  static std::string name() { return "ChalmersWithConstraintGuard"; }
  amr::Criteria::Type type() override { return amr::Criteria::Type::hp; }
  std::string observation_name() override { return name(); }
  using compute_tags_for_observation_box = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTags, typename Metavariables>
  std::array<amr::Flag, Dim> operator()(
      const db::DataBox<DbTags>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id) const {
    auto flags = criterion_(box, cache, element_id);
    std::array<bool, Dim> decrease{};
    for (size_t d = 0; d < Dim; ++d) {
      decrease[d] = flags[d] == amr::Flag::DecreaseResolution;
    }
    if (alg::any_of(decrease, [](const bool value) { return value; })) {
      const double change = p_coarsening_constraint_change(
          db::get<Tags::Psi>(box),
          db::get<Tags::Phi<Dim, Frame::Inertial>>(box),
          db::get<domain::Tags::Mesh<Dim>>(box), decrease,
          db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                                Frame::Inertial>>(box));
      if (not std::isfinite(change) or change > target_) {
        for (size_t d = 0; d < Dim; ++d) {
          if (decrease[d]) {
            flags[d] = amr::Flag::DoNothing;
          }
        }
      }
    }
    return flags;
  }

  void pup(PUP::er& p) override {
    amr::Criterion::pup(p);
    p | criterion_;
    p | target_;
  }

 private:
  FieldCriterion criterion_{};
  double target_{0.0};
};

template <size_t Dim>
PUP::able::PUP_ID ChalmersWithConstraintGuard<Dim>::my_PUP_ID = 0;  // NOLINT
}  // namespace ScalarWave
