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
class Outflowing final : public BoundaryCondition<Dim> {
 public:
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Outflowing boundary conditions. These check that all characteristic "
      "fields are outflowing and throw and error otherwise."};
  Outflowing() = default;
  /// \cond
  Outflowing(Outflowing&&) noexcept = default;
  Outflowing& operator=(Outflowing&&) noexcept = default;
  Outflowing(const Outflowing&) = default;
  Outflowing& operator=(const Outflowing&) = default;
  /// \endcond
  ~Outflowing() override = default;

  explicit Outflowing(CkMigrateMessage* msg) noexcept;

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, Outflowing);

  auto get_clone() const noexcept -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::TimeDerivative;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<Pi, Phi<Dim>, Psi>;
  using dg_interior_temporary_tags =
      tmpl::list<Tags::ConstraintGamma1, gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<Dim, Frame::Inertial, DataVector>>;
  using dg_interior_dt_vars_tags = tmpl::list<>;
  using dg_interior_deriv_vars_tags = tmpl::list<>;
  using dg_gridless_tags = tmpl::list<>;

  std::optional<std::string> dg_time_derivative(
      gsl::not_null<Scalar<DataVector>*> dt_pi_correction,
      gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
      dt_phi_correction,
      gsl::not_null<Scalar<DataVector>*> dt_psi_correction,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      face_mesh_velocity,
      const tnsr::i<DataVector, Dim>& normal_covector,
      const Scalar<DataVector>& pi, const tnsr::i<DataVector, Dim>& phi,
      const Scalar<DataVector>& psi,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim>& shift) const noexcept;
};
}  // namespace CurvedScalarWave::BoundaryConditions
