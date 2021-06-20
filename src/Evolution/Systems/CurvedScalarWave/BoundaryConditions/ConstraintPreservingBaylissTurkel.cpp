// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/ConstraintPreservingBaylissTurkel.hpp"

#include <cstddef>
#include <memory>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace CurvedScalarWave::BoundaryConditions {

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
ConstraintPreservingBaylissTurkel<Dim>::get_clone() const noexcept {
  return std::make_unique<ConstraintPreservingBaylissTurkel>(*this);
}

template <size_t Dim>
void ConstraintPreservingBaylissTurkel<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
}

template <size_t Dim>
ConstraintPreservingBaylissTurkel<Dim>::ConstraintPreservingBaylissTurkel(
    CkMigrateMessage* const msg) noexcept
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::optional<std::string>
ConstraintPreservingBaylissTurkel<Dim>::dg_time_derivative(
    const gsl::not_null<Scalar<DataVector>*> dt_pi_correction,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        dt_phi_correction,
    const gsl::not_null<Scalar<DataVector>*> dt_psi_correction,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim>& normal_covector,
    const tnsr::I<DataVector, Dim>& normal_vector, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim>& phi, const Scalar<DataVector>& psi,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse, const tnsr::i<DataVector, Dim>& d_lapse,
    const tnsr::I<DataVector, Dim>& shift,
    const tnsr::iJ<DataVector, Dim>& d_shift, const Scalar<DataVector>& dt_pi,
    const tnsr::i<DataVector, Dim>& dt_phi, const Scalar<DataVector>& dt_psi,
    const tnsr::i<DataVector, Dim>& d_pi, const tnsr::i<DataVector, Dim>& d_psi,
    const tnsr::ij<DataVector, Dim>& d_phi) const noexcept {
  const DataVector inv_radius = 1. / get(magnitude(coords));

  /*
  auto dt_char_fields = characteristic_fields(
      gamma2, inverse_spatial_metric, dt_psi, dt_pi, dt_phi, normal_covector);

  // set them to the negative value, so the correction starts from 0
  dt_char_fields *= -1.;

  auto& dt_VPsi = get<Tags::VPsi>(dt_char_fields);
  auto& dt_VZero = get<Tags::VZero<Dim>>(dt_char_fields);
  auto& dt_VPlus = get<Tags::VPlus>(dt_char_fields);
  auto& dt_VMinus = get<Tags::VMinus>(dt_char_fields);

  get(dt_VPsi) += get(dot_product(shift, phi)) - get(lapse) * get(pi);

  auto tmp = make_with_value<tnsr::i<DataVector, Dim, Frame::Inertial>>(pi, 0.);
  for (size_t j = 0; j < Dim; ++j) {
    tmp.get(j) = -get(lapse) * d_pi.get(j) - get(pi) * d_lapse.get(j);
    for (size_t k = 0; k < Dim; ++k) {
      tmp.get(j) +=
          d_shift.get(j, k) * phi.get(k) + shift.get(k) * d_phi.get(j, k);
    }
  }

  for (size_t i = 0; i < Dim; ++i) {
    dt_VZero.get(i) += tmp.get(i);
    for (size_t j = 0; j < Dim; ++j) {
      dt_VZero.get(i) -=
          normal_covector.get(i) *
          raise_or_lower_index(normal_covector, inverse_spatial_metric).get(j) *
          tmp.get(j);
    }
  }
  get(dt_VMinus) -= get(gamma2) * get(dt_VPsi);  // * inv_radius;
*/

  auto char_speeds =
      characteristic_speeds(gamma1, lapse, shift, normal_covector);

  if (face_mesh_velocity.has_value()) {
    const auto face_speed = dot_product(normal_covector, *face_mesh_velocity);
    for (auto& char_speed : char_speeds) {
      char_speed -= get(face_speed);
    }
  }

  get(*dt_psi_correction) =
      get<0>(normal_vector) * (get<0>(d_psi) - get<0>(phi));
  for (size_t i = 1; i < Dim; ++i) {
    get(*dt_psi_correction) +=
        normal_vector.get(i) * (d_psi.get(i) - phi.get(i));
  }

  // Compute dt Phi 2-index constraint correction
  for (size_t i = 0; i < Dim; ++i) {
    dt_phi_correction->get(i) =
        get<0>(normal_vector) * (d_phi.get(0, i) - d_phi.get(i, 0));
    for (size_t j = 1; j < Dim; ++j) {
      dt_phi_correction->get(i) +=
          normal_vector.get(j) * (d_phi.get(j, i) - d_phi.get(i, j));
    }
  }

  for (size_t i = 0; i < get(*dt_psi_correction).size(); ++i) {
    get(*dt_psi_correction)[i] *= std::min(0., gsl::at(char_speeds[0], i));
    for (size_t j = 0; j < Dim; ++j) {
      dt_phi_correction->get(j)[i] *=
          0.5 * std::min(0., gsl::at(char_speeds[1], i));
    }
  }

  get(*dt_pi_correction) =
      -get(dt_pi) - 2.0 * inv_radius * (2.0 * get(pi) - inv_radius * get(psi));
  for (size_t i = 0; i < Dim; ++i) {
    get(*dt_pi_correction) +=
        normal_vector.get(i) *
        (dt_phi.get(i) - d_pi.get(i) + 4.0 * inv_radius * phi.get(i));
    for (size_t j = 0; j < Dim; ++j) {
      get(*dt_pi_correction) +=
          normal_vector.get(i) * normal_vector.get(j) * d_phi.get(i, j);
    }
  }
  get(*dt_pi_correction) += get(gamma2) * get(*dt_psi_correction);

  return {};
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID ConstraintPreservingBaylissTurkel<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class ConstraintPreservingBaylissTurkel<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryConditions
