// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/DiscontinuousGalerkin/TimeDerivativeDecisions.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl

class DataVector;
/// \endcond

namespace ScalarWave {
/*!
 * \brief Compute the time derivatives for scalar wave system
 */
template <size_t Dim>
struct TimeDerivative {
  using temporary_tags = tmpl::list<Tags::ConstraintGamma2>;
  using argument_tags =
      tmpl::list<Tags::Pi, Tags::Phi<Dim>, Tags::ConstraintGamma2>;

  template <typename DataType>
  KOKKOS_FUNCTION static evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      gsl::not_null<Scalar<DataType>*> dt_psi,
      gsl::not_null<Scalar<DataType>*> dt_pi,
      gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*> dt_phi,

      gsl::not_null<Scalar<DataType>*> result_gamma2,

      // Partial derivative arguments. Listed in the system struct as
      // gradient_variables.
      const tnsr::i<DataType, Dim, Frame::Inertial>& d_psi,
      const tnsr::i<DataType, Dim, Frame::Inertial>& d_pi,
      const tnsr::ij<DataType, Dim, Frame::Inertial>& d_phi,

      // Terms list in argument_tags above
      const Scalar<DataType>& pi,
      const tnsr::i<DataType, Dim, Frame::Inertial>& phi,
      const Scalar<DataType>& gamma2);
};

template <size_t Dim>
template <typename DataType>
KOKKOS_FUNCTION evolution::dg::TimeDerivativeDecisions<Dim>
TimeDerivative<Dim>::apply(
    const gsl::not_null<Scalar<DataType>*> dt_psi,
    const gsl::not_null<Scalar<DataType>*> dt_pi,
    const gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*> dt_phi,

    const gsl::not_null<Scalar<DataType>*> result_gamma2,

    const tnsr::i<DataType, Dim, Frame::Inertial>& d_psi,
    const tnsr::i<DataType, Dim, Frame::Inertial>& d_pi,
    const tnsr::ij<DataType, Dim, Frame::Inertial>& d_phi,

    const Scalar<DataType>& pi,
    const tnsr::i<DataType, Dim, Frame::Inertial>& phi,
    const Scalar<DataType>& gamma2) {
  // The constraint damping parameter gamma2 is needed for boundary corrections,
  // which means we need it as a temporary tag in order to project it to the
  // boundary. We prevent slicing/projecting directly from the volume to prevent
  // people from adding many compute tags to the DataBox, instead preferring
  // quantities be computed inside the TimeDerivative/Flux/Source structs. This
  // keeps related code together and makes figuring out where something is
  // computed a lot easier.
  *result_gamma2 = gamma2;

  get(*dt_psi) = -get(pi);
  get(*dt_pi) = -get<0, 0>(d_phi);
  for (size_t d = 1; d < Dim; ++d) {
    get(*dt_pi) -= d_phi.get(d, d);
  }
  for (size_t d = 0; d < Dim; ++d) {
    dt_phi->get(d) = -d_pi.get(d) + get(gamma2) * (d_psi.get(d) - phi.get(d));
  }
  return {true};
}
}  // namespace ScalarWave
