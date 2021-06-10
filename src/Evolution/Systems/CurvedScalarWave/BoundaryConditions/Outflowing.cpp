// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/Outflowing.hpp"

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
Outflowing<Dim>::get_clone() const noexcept {
  return std::make_unique<Outflowing>(*this);
}

template <size_t Dim>
void Outflowing<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
}

template <size_t Dim>
Outflowing<Dim>::Outflowing(CkMigrateMessage* const msg) noexcept
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::optional<std::string> Outflowing<Dim>::dg_time_derivative(
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
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim>& shift) const noexcept {
  auto char_speeds =
      characteristic_speeds(gamma1, lapse, shift, normal_covector);

  if (face_mesh_velocity.has_value()) {
    const auto face_speed = dot_product(normal_covector, *face_mesh_velocity);
    for (auto& char_speed : char_speeds) {
      char_speed -= get(face_speed);
    }
  }
  for (const auto& char_speed : char_speeds) {
    for (const auto& char_speed_point : char_speed) {
      ASSERT(char_speed_point > 0.,
             "characteristic speed should be outflowing but detected negative "
             "speed "
                 << char_speed_point);
    }
  }
  get(*dt_pi_correction) = 0.;
  get(*dt_psi_correction) = 0.;
  for (auto& val : *dt_phi_correction) {
    val = 0.;
  }
  return {};
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID Outflowing<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) template class Outflowing<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryConditions
