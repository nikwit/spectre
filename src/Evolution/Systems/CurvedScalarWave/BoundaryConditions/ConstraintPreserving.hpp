// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace CurvedScalarWave::BoundaryConditions {
template <size_t Dim>
class ConstraintPreserving final : public BoundaryCondition<Dim> {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Constraint-preserving spherical radiation boundary conditions setting "
      "the time derivatives of Psi, Phi, and Pi to avoid incoming constraint "
      "violations as according to Holst (2004)"};
  ConstraintPreserving() = default;
  /// \cond
  ConstraintPreserving(ConstraintPreserving&&) noexcept = default;
  ConstraintPreserving& operator=(ConstraintPreserving&&) noexcept = default;
  ConstraintPreserving(const ConstraintPreserving&) = default;
  ConstraintPreserving& operator=(const ConstraintPreserving&) = default;
  /// \endcond
  ~ConstraintPreserving() override = default;

  explicit ConstraintPreserving(CkMigrateMessage* msg) noexcept;

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, ConstraintPreserving);

  auto get_clone() const noexcept -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::TimeDerivative;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<Pi, Phi<Dim>, Psi>;
  using dg_interior_temporary_tags = tmpl::list<
      domain::Tags::Coordinates<Dim, Frame::Inertial>, Tags::ConstraintGamma1,
      Tags::ConstraintGamma2, gr::Tags::Lapse<DataVector>,
      ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<Dim>,
                    Frame::Inertial>,
      gr::Tags::Shift<Dim, Frame::Inertial, DataVector>,
      ::Tags::deriv<gr::Tags::Shift<Dim, Frame::Inertial, DataVector>,
                    tmpl::size_t<Dim>, Frame::Inertial>,
      gr::Tags::InverseSpatialMetric<Dim, Frame::Inertial, DataVector>>;
  using dg_interior_dt_vars_tags =
      tmpl::list<::Tags::dt<Pi>, ::Tags::dt<Phi<Dim>>, ::Tags::dt<Psi>>;
  using dg_interior_deriv_vars_tags =
      tmpl::list<::Tags::deriv<Pi, tmpl::size_t<Dim>, Frame::Inertial>,
                 ::Tags::deriv<Phi<Dim>, tmpl::size_t<Dim>, Frame::Inertial>>;
  using dg_gridless_tags = tmpl::list<>;

  std::optional<std::string> dg_time_derivative(
      gsl::not_null<Scalar<DataVector>*> dt_pi_correction,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
          dt_phi_correction,
      gsl::not_null<Scalar<DataVector>*> dt_psi_correction,
      const std::optional<tnsr::I<DataVector, Dim>>& face_mesh_velocity,
      const tnsr::i<DataVector, Dim>& normal_covector,
      const tnsr::I<DataVector, Dim>& normal_vector,
      const Scalar<DataVector>& pi, const tnsr::i<DataVector, Dim>& phi,
      const Scalar<DataVector>& psi,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const Scalar<DataVector>& lapse, const tnsr::i<DataVector, Dim>& d_lapse,
      const tnsr::I<DataVector, Dim>& shift,
      const tnsr::iJ<DataVector, Dim>& d_shift,
      const tnsr::II<DataVector, Dim>& inverse_spatial_metric,
      const Scalar<DataVector>& dt_pi, const tnsr::i<DataVector, Dim>& dt_phi,
      const Scalar<DataVector>& dt_psi, const tnsr::i<DataVector, Dim>& d_pi,
      const tnsr::ij<DataVector, Dim>& d_phi) const noexcept;
};
}  // namespace CurvedScalarWave::BoundaryConditions
