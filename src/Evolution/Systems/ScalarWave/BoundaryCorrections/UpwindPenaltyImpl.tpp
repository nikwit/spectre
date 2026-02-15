// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once
#pragma nv_diag_suppress 20014


#include <cstddef>
#include <type_traits>

#include "DataStructures/Blaze/StepFunction.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::BoundaryCorrections::detail {

template <size_t Dim, typename DataType>
KOKKOS_FUNCTION double dg_package_data_impl(
    const gsl::not_null<Scalar<DataType>*> packaged_char_speed_v_psi,
    const gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*>
        packaged_char_speed_v_zero,
    const gsl::not_null<Scalar<DataType>*> packaged_char_speed_v_plus,
    const gsl::not_null<Scalar<DataType>*> packaged_char_speed_v_minus,
    const gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*>
        packaged_char_speed_n_times_v_minus,
    const gsl::not_null<Scalar<DataType>*> packaged_char_speed_gamma2_v_psi,
    const gsl::not_null<tnsr::i<DataType, 3, Frame::Inertial>*>
        packaged_char_speeds,
    const Scalar<DataType>& psi, const Scalar<DataType>& pi,
    const tnsr::i<DataType, Dim, Frame::Inertial>& phi,
    const Scalar<DataType>& constraint_gamma2,
    const tnsr::i<DataType, Dim, Frame::Inertial>& normal_covector,
    const Scalar<DataType>* const normal_dot_mesh_velocity) {
  if (normal_dot_mesh_velocity != nullptr) {
    get<0>(*packaged_char_speeds) = -get(*normal_dot_mesh_velocity);
    get<1>(*packaged_char_speeds) = 1.0 - get(*normal_dot_mesh_velocity);
    get<2>(*packaged_char_speeds) = -1.0 - get(*normal_dot_mesh_velocity);
  } else {
    get<0>(*packaged_char_speeds) = 0.0;
    get<1>(*packaged_char_speeds) = 1.0;
    get<2>(*packaged_char_speeds) = -1.0;
  }

  get(*packaged_char_speed_gamma2_v_psi) = get(constraint_gamma2) * get(psi);

  DataType normal_dot_phi = normal_covector.get(0) * phi.get(0);
  for (size_t i = 1; i < Dim; ++i) {
    normal_dot_phi += normal_covector.get(i) * phi.get(i);
  }

  for (size_t i = 0; i < Dim; ++i) {
    packaged_char_speed_v_zero->get(i) =
        get<0>(*packaged_char_speeds) *
        (phi.get(i) - normal_covector.get(i) * normal_dot_phi);
  }

  get(*packaged_char_speed_v_plus) =
      get<1>(*packaged_char_speeds) *
      (get(pi) + normal_dot_phi - get(*packaged_char_speed_gamma2_v_psi));
  get(*packaged_char_speed_v_minus) =
      get<2>(*packaged_char_speeds) *
      (get(pi) - normal_dot_phi - get(*packaged_char_speed_gamma2_v_psi));

  for (size_t d = 0; d < Dim; ++d) {
    packaged_char_speed_n_times_v_plus->get(d) =
        get(*packaged_char_speed_v_plus) * normal_covector.get(d);
    packaged_char_speed_n_times_v_minus->get(d) =
        get(*packaged_char_speed_v_minus) * normal_covector.get(d);
  }

  get(*packaged_char_speed_v_psi) = get<0>(*packaged_char_speeds) * get(psi);
  get(*packaged_char_speed_gamma2_v_psi) *= get<0>(*packaged_char_speeds);

  if constexpr (std::is_same_v<DataType, double>) {
    const double a = get<0>(*packaged_char_speeds);
    const double b = get<1>(*packaged_char_speeds);
    const double c = get<2>(*packaged_char_speeds);
    return a > b ? (a > c ? a : c) : (b > c ? b : c);
  } else {
    return max(max(get<0>(*packaged_char_speeds), get<1>(*packaged_char_speeds),
                   get<2>(*packaged_char_speeds)));
  }
}

template <size_t Dim, typename DataType>
KOKKOS_FUNCTION void dg_boundary_terms_impl(
    const gsl::not_null<Scalar<DataType>*> psi_boundary_correction,
    const gsl::not_null<Scalar<DataType>*> pi_boundary_correction,
    const gsl::not_null<tnsr::i<DataType, Dim, Frame::Inertial>*>
        phi_boundary_correction,
    const Scalar<DataType>& char_speed_v_psi_int,
    const tnsr::i<DataType, Dim, Frame::Inertial>& char_speed_v_zero_int,
    const Scalar<DataType>& char_speed_v_plus_int,
    const Scalar<DataType>& char_speed_v_minus_int,
    const tnsr::i<DataType, Dim, Frame::Inertial>&
        char_speed_normal_times_v_plus_int,
    const tnsr::i<DataType, Dim, Frame::Inertial>&
        char_speed_normal_times_v_minus_int,
    const Scalar<DataType>& char_speed_constraint_gamma2_v_psi_int,
    const tnsr::i<DataType, 3, Frame::Inertial>& char_speeds_int,
    const Scalar<DataType>& char_speed_v_psi_ext,
    const tnsr::i<DataType, Dim, Frame::Inertial>& char_speed_v_zero_ext,
    const Scalar<DataType>& char_speed_v_plus_ext,
    const Scalar<DataType>& char_speed_v_minus_ext,
    const tnsr::i<DataType, Dim, Frame::Inertial>&
        char_speed_minus_normal_times_v_plus_ext,
    const tnsr::i<DataType, Dim, Frame::Inertial>&
        char_speed_minus_normal_times_v_minus_ext,
    const Scalar<DataType>& char_speed_constraint_gamma2_v_psi_ext,
    const tnsr::i<DataType, 3, Frame::Inertial>& char_speeds_ext) {
  DataType weighted_lambda_psi_int{};
  DataType weighted_lambda_psi_ext{};
  DataType weighted_lambda_zero_int{};
  DataType weighted_lambda_zero_ext{};
  DataType weighted_lambda_plus_int{};
  DataType weighted_lambda_plus_ext{};
  DataType weighted_lambda_minus_int{};
  DataType weighted_lambda_minus_ext{};
  if constexpr (std::is_same_v<DataType, double>) {
    const auto heaviside = [](const double x) {
      return x < 0.0 ? 0.0 : 1.0;
    };
    weighted_lambda_psi_int = heaviside(-char_speeds_int[0]);
    weighted_lambda_psi_ext = -heaviside(char_speeds_ext[0]);
    weighted_lambda_zero_int = heaviside(-char_speeds_int[0]);
    weighted_lambda_zero_ext = -heaviside(char_speeds_ext[0]);
    weighted_lambda_plus_int = heaviside(-char_speeds_int[1]);
    weighted_lambda_plus_ext = -heaviside(char_speeds_ext[1]);
    weighted_lambda_minus_int = heaviside(-char_speeds_int[2]);
    weighted_lambda_minus_ext = -heaviside(char_speeds_ext[2]);
  } else {
    weighted_lambda_psi_int = step_function(-char_speeds_int[0]);
    weighted_lambda_psi_ext = -step_function(char_speeds_ext[0]);
    weighted_lambda_zero_int = step_function(-char_speeds_int[0]);
    weighted_lambda_zero_ext = -step_function(char_speeds_ext[0]);
    weighted_lambda_plus_int = step_function(-char_speeds_int[1]);
    weighted_lambda_plus_ext = -step_function(char_speeds_ext[1]);
    weighted_lambda_minus_int = step_function(-char_speeds_int[2]);
    weighted_lambda_minus_ext = -step_function(char_speeds_ext[2]);
  }

  get(*psi_boundary_correction) =
      weighted_lambda_psi_ext * get(char_speed_v_psi_ext) -
      weighted_lambda_psi_int * get(char_speed_v_psi_int);

  get(*pi_boundary_correction) =
      0.5 * (weighted_lambda_plus_ext * get(char_speed_v_plus_ext) +
             weighted_lambda_minus_ext * get(char_speed_v_minus_ext)) +
      weighted_lambda_psi_ext * get(char_speed_constraint_gamma2_v_psi_ext) -
      0.5 * (weighted_lambda_plus_int * get(char_speed_v_plus_int) +
             weighted_lambda_minus_int * get(char_speed_v_minus_int)) -
      weighted_lambda_psi_int * get(char_speed_constraint_gamma2_v_psi_int);

  for (size_t d = 0; d < Dim; ++d) {
    phi_boundary_correction->get(d) =
        0.5 * (weighted_lambda_plus_ext *
                   char_speed_minus_normal_times_v_plus_ext.get(d) -
               weighted_lambda_minus_ext *
                   char_speed_minus_normal_times_v_minus_ext.get(d)) +
        weighted_lambda_zero_ext * char_speed_v_zero_ext.get(d) -
        0.5 * (weighted_lambda_plus_int *
                   char_speed_normal_times_v_plus_int.get(d) -
               weighted_lambda_minus_int *
                   char_speed_normal_times_v_minus_int.get(d)) -
        weighted_lambda_zero_int * char_speed_v_zero_int.get(d);
  }
}

}  // namespace ScalarWave::BoundaryCorrections::detail
