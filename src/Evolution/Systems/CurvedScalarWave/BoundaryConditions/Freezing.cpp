// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Freezing.hpp"

#include <cstddef>
#include <memory>
#include <optional>
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
Freezing<Dim>::get_clone() const noexcept {
  return std::make_unique<Freezing>(*this);
}

template <size_t Dim>
void Freezing<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
}

template <size_t Dim>
Freezing<Dim>::Freezing(CkMigrateMessage* const msg) noexcept
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::optional<std::string> Freezing<Dim>::dg_time_derivative(
    gsl::not_null<Scalar<DataVector>*> dt_pi_correction,
    gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*> dt_phi_correction,
    gsl::not_null<Scalar<DataVector>*> dt_psi_correction,
    const std::optional<tnsr::I<DataVector, Dim>>& face_mesh_velocity,
    const tnsr::i<DataVector, Dim>& normal_covector,
    const Scalar<DataVector>& pi, const tnsr::i<DataVector, Dim>& phi,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim>& shift,
    const tnsr::II<DataVector, Dim>& inverse_spatial_metric,
    const Scalar<DataVector>& dt_pi, const tnsr::i<DataVector, Dim>& dt_phi,
    const Scalar<DataVector>& dt_psi) const noexcept {
  auto dt_char_fields = characteristic_fields(
      gamma2, inverse_spatial_metric, dt_psi, dt_pi, dt_phi, normal_covector);

  // set them to the negative value
  dt_char_fields *= -1.;

  auto& dt_VPsi = get<Tags::VPsi>(dt_char_fields);
  auto& dt_VZero = get<Tags::VZero<Dim>>(dt_char_fields);
  auto& dt_VPlus = get<Tags::VPlus>(dt_char_fields);
  auto& dt_VMinus = get<Tags::VMinus>(dt_char_fields);

  get(dt_VMinus) -= get(gamma2) * get(dt_psi);

  const auto char_speeds =
      characteristic_speeds(gamma1, lapse, shift, normal_covector);

  const auto set_outgoing_field_corrections_to_zero =
      [](DataVector& char_field, const DataVector& char_speed) {
        for (size_t i = 0; i < char_field.size(); ++i) {
          if (char_speed[i] > 0.) {
            char_field[i] = 0.;
          }
        }
      };

  set_outgoing_field_corrections_to_zero(get(dt_VPsi), char_speeds[0]);
  for (size_t i = 0; i < Dim; ++i) {
    set_outgoing_field_corrections_to_zero(dt_VZero.get(i), char_speeds[1]);
  }
  set_outgoing_field_corrections_to_zero(get(dt_VPlus), char_speeds[2]);
  set_outgoing_field_corrections_to_zero(get(dt_VMinus), char_speeds[3]);

  // Convert them to desired values on dt<U>
  const auto dt_evolved_fields_corrections =
      evolved_fields_from_characteristic_fields(
          gamma2, dt_VPsi, dt_VZero, dt_VPlus, dt_VMinus, normal_covector);

  // change signature of `evolved_fields_from_characteristics` to return tensors
  // individually by reference to avoid this extra allocation and copy
  *dt_pi_correction = get<Pi>(dt_evolved_fields_corrections);
  *dt_psi_correction = get<Psi>(dt_evolved_fields_corrections);
  *dt_phi_correction = get<Phi<Dim>>(dt_evolved_fields_corrections);
  return {};
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID Freezing<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) template class Freezing<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryConditions
