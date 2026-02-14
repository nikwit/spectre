// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma nv_diag_suppress 20014

#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenaltyImpl.tpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::BoundaryCorrections {

template <size_t Dim>
UpwindPenalty<Dim>::UpwindPenalty(CkMigrateMessage* msg)
    : BoundaryCorrection<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<BoundaryCorrection<Dim>> UpwindPenalty<Dim>::get_clone() const {
  return std::make_unique<UpwindPenalty>(*this);
}

template <size_t Dim>
void UpwindPenalty<Dim>::pup(PUP::er& p) {
  BoundaryCorrection<Dim>::pup(p);
}

template <size_t Dim>
double UpwindPenalty<Dim>::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_psi,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_v_zero,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_plus,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_v_minus,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_minus,
    const gsl::not_null<Scalar<DataVector>*> packaged_char_speed_gamma2_v_psi,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        packaged_char_speeds,

    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,

    const Scalar<DataVector>& constraint_gamma2,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
  const auto* normal_dot_mesh_velocity_ptr =
      normal_dot_mesh_velocity.has_value() ? &(*normal_dot_mesh_velocity)
                                           : nullptr;
  return detail::dg_package_data_impl(
      packaged_char_speed_v_psi, packaged_char_speed_v_zero,
      packaged_char_speed_v_plus, packaged_char_speed_v_minus,
      packaged_char_speed_n_times_v_plus, packaged_char_speed_n_times_v_minus,
      packaged_char_speed_gamma2_v_psi, packaged_char_speeds, psi, pi, phi,
      constraint_gamma2, normal_covector, normal_dot_mesh_velocity_ptr);
}

template <size_t Dim>
void UpwindPenalty<Dim>::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
    const gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        phi_boundary_correction,

    const Scalar<DataVector>& char_speed_v_psi_int,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& char_speed_v_zero_int,
    const Scalar<DataVector>& char_speed_v_plus_int,
    const Scalar<DataVector>& char_speed_v_minus_int,
    const tnsr::i<DataVector, Dim, Frame::Inertial>&
        char_speed_normal_times_v_plus_int,
    const tnsr::i<DataVector, Dim, Frame::Inertial>&
        char_speed_normal_times_v_minus_int,
    const Scalar<DataVector>& char_speed_constraint_gamma2_v_psi_int,
    const tnsr::i<DataVector, 3, Frame::Inertial>& char_speeds_int,

    const Scalar<DataVector>& char_speed_v_psi_ext,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& char_speed_v_zero_ext,
    const Scalar<DataVector>& char_speed_v_plus_ext,
    const Scalar<DataVector>& char_speed_v_minus_ext,
    const tnsr::i<DataVector, Dim, Frame::Inertial>&
        char_speed_minus_normal_times_v_plus_ext,
    const tnsr::i<DataVector, Dim, Frame::Inertial>&
        char_speed_minus_normal_times_v_minus_ext,
    const Scalar<DataVector>& char_speed_constraint_gamma2_v_psi_ext,
    const tnsr::i<DataVector, 3, Frame::Inertial>& char_speeds_ext,
    dg::Formulation /*dg_formulation*/) const {
  detail::dg_boundary_terms_impl(
      psi_boundary_correction, pi_boundary_correction, phi_boundary_correction,
      char_speed_v_psi_int, char_speed_v_zero_int, char_speed_v_plus_int,
      char_speed_v_minus_int, char_speed_normal_times_v_plus_int,
      char_speed_normal_times_v_minus_int, char_speed_constraint_gamma2_v_psi_int,
      char_speeds_int, char_speed_v_psi_ext, char_speed_v_zero_ext,
      char_speed_v_plus_ext, char_speed_v_minus_ext,
      char_speed_minus_normal_times_v_plus_ext,
      char_speed_minus_normal_times_v_minus_ext,
      char_speed_constraint_gamma2_v_psi_ext, char_speeds_ext);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID UpwindPenalty<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(_, data) template class UpwindPenalty<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM


}  // namespace ScalarWave::BoundaryCorrections
