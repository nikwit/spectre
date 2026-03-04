// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/UpwindPenalty.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/ApplyTensorYlmFilter.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::BoundaryCorrections {
template <size_t Dim>
UpwindPenalty<Dim>::UpwindPenalty(CkMigrateMessage* msg)
    : BoundaryCorrection(msg) {}

template <size_t Dim>
std::unique_ptr<evolution::BoundaryCorrection> UpwindPenalty<Dim>::get_clone()
    const {
  return std::make_unique<UpwindPenalty>(*this);
}

template <size_t Dim>
void UpwindPenalty<Dim>::pup(PUP::er& p) {
  BoundaryCorrection::pup(p);
}

template <size_t Dim>
double UpwindPenalty<Dim>::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> packaged_v_psi,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_v_zero,
    const gsl::not_null<Scalar<DataVector>*> packaged_v_plus,
    const gsl::not_null<Scalar<DataVector>*> packaged_v_minus,
    const gsl::not_null<Scalar<DataVector>*> packaged_gamma2,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        packaged_interface_unit_normal,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        packaged_interface_unit_normal_vector,
    const gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*>
        packaged_char_speeds,

    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& phi,

    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
    const Scalar<DataVector>& constraint_gamma1,
    const Scalar<DataVector>& constraint_gamma2,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& interface_unit_normal,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        interface_unit_normal_vector,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
  *packaged_gamma2 = constraint_gamma2;
  *packaged_interface_unit_normal = interface_unit_normal;
  *packaged_interface_unit_normal_vector = interface_unit_normal_vector;
  characteristic_fields(packaged_v_psi, packaged_v_zero, packaged_v_plus,
                        packaged_v_minus, constraint_gamma2, psi, pi, phi,
                        interface_unit_normal, interface_unit_normal_vector);
  characteristic_speeds(packaged_char_speeds, constraint_gamma1, lapse, shift,
                        interface_unit_normal);
  if (normal_dot_mesh_velocity.has_value()) {
    get<0>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<1>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<2>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
    get<3>(*packaged_char_speeds) -= get(*normal_dot_mesh_velocity);
  }

  return max(max(get<0>(*packaged_char_speeds), get<1>(*packaged_char_speeds),
                 get<2>(*packaged_char_speeds)));
}

template <size_t Dim>
void UpwindPenalty<Dim>::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> psi_boundary_correction,
    const gsl::not_null<Scalar<DataVector>*> pi_boundary_correction,
    const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*>
        phi_boundary_correction,

    const Scalar<DataVector>& v_psi_int,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& v_zero_int,
    const Scalar<DataVector>& v_plus_int, const Scalar<DataVector>& v_minus_int,
    const Scalar<DataVector>& gamma2_int,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& interface_unit_normal_int,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        interface_unit_normal_vector_int,
    const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_int,

    const Scalar<DataVector>& v_psi_ext,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& v_zero_ext,
    const Scalar<DataVector>& v_plus_ext, const Scalar<DataVector>& v_minus_ext,
    const Scalar<DataVector>& gamma2_ext,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& interface_unit_normal_ext,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        interface_unit_normal_vector_ext,
    const tnsr::a<DataVector, 3, Frame::Inertial>& char_speeds_ext,
    dg::Formulation dg_formulation, const Mesh<Dim>& volume_mesh,
    const ylm::TensorYlm::CurvedScalarWaveTensorYlmFilter& tensor_ylm_filter)
    const {
  (void)dg_formulation;
  (void)interface_unit_normal_vector_int;
  const Scalar<DataVector>* v_psi_ext_ptr = &v_psi_ext;
  const tnsr::i<DataVector, Dim, Frame::Inertial>* v_zero_ext_ptr = &v_zero_ext;
  const Scalar<DataVector>* v_plus_ext_ptr = &v_plus_ext;
  const Scalar<DataVector>* v_minus_ext_ptr = &v_minus_ext;
  std::optional<
      Variables<ylm::TensorYlm::filter_detail::sw_vars_list<Frame::Inertial>>>
      filtered_ext_vars{};
  std::optional<Variables<
      tmpl::list<Tags::VPsi, Tags::VZero<Dim>, Tags::VPlus, Tags::VMinus>>>
      filtered_ext_char_vars{};

  if constexpr (Dim == 3) {
    if (volume_mesh.basis(1) == Spectral::Basis::SphericalHarmonic and
        volume_mesh.basis(2) == Spectral::Basis::SphericalHarmonic) {
      const Mesh<2> angular_mesh = volume_mesh.slice_away(0);
      const size_t num_pts = get(v_psi_ext).size();
      filtered_ext_vars.emplace(num_pts);
      ASSERT(angular_mesh.number_of_grid_points() == num_pts,
             "Expected number of points in the angular mesh to match the "
             "number of points in the volume mesh, but got "
                 << angular_mesh.number_of_grid_points() << " and " << num_pts);

      evolved_fields_from_characteristic_fields(
          make_not_null(&*filtered_ext_vars), gamma2_ext, v_psi_ext, v_zero_ext,
          v_plus_ext, v_minus_ext, interface_unit_normal_ext);

      InverseJacobian<DataVector, 3, Frame::Grid, Frame::Inertial>
          jac_grid_to_inertial{num_pts};
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          jac_grid_to_inertial.get(i, j) = (i == j) ? 1.0 : 0.0;
        }
      }
      tensor_ylm_filter(make_not_null(&*filtered_ext_vars), angular_mesh,
                        jac_grid_to_inertial);

      filtered_ext_char_vars.emplace(num_pts);
      characteristic_fields(
          make_not_null(&*filtered_ext_char_vars), gamma2_ext,
          get<Tags::Psi>(*filtered_ext_vars), get<Tags::Pi>(*filtered_ext_vars),
          get<Tags::Phi<Dim>>(*filtered_ext_vars), interface_unit_normal_ext,
          interface_unit_normal_vector_ext);

      v_psi_ext_ptr = &get<Tags::VPsi>(*filtered_ext_char_vars);
      v_zero_ext_ptr = &get<Tags::VZero<Dim>>(*filtered_ext_char_vars);
      v_plus_ext_ptr = &get<Tags::VPlus>(*filtered_ext_char_vars);
      v_minus_ext_ptr = &get<Tags::VMinus>(*filtered_ext_char_vars);
    } else {
      ERROR("not at spherical harmonic boundary?");
    }
  } else {
    (void)volume_mesh;
    (void)tensor_ylm_filter;
    (void)gamma2_ext;
    (void)interface_unit_normal_ext;
    (void)interface_unit_normal_vector_ext;
  }

  const auto& v_psi_ext_used = *v_psi_ext_ptr;
  const auto& v_zero_ext_used = *v_zero_ext_ptr;
  const auto& v_plus_ext_used = *v_plus_ext_ptr;
  const auto& v_minus_ext_used = *v_minus_ext_ptr;
  // The implementations assumes that the external unit normal vector is
  // exactly negative the internal unit normal vector and that the internal
  // gamma2 is equal to the external gamma2. For Gauss quadrature this will not
  // be exactly true due to interpolation error but that should not matter.
  get(*psi_boundary_correction) = -step_function(get<0>(char_speeds_ext)) *
                                      get<0>(char_speeds_ext) *
                                      get(v_psi_ext_used) -
                                  step_function(-get<0>(char_speeds_int)) *
                                      get<0>(char_speeds_int) * get(v_psi_int);

  auto& temp_1 = get<Dim - 1>(*phi_boundary_correction);
  temp_1 = -step_function(get<2>(char_speeds_ext)) * get<2>(char_speeds_ext) *
               get(v_plus_ext_used) -
           step_function(-get<3>(char_speeds_int)) * get<3>(char_speeds_int) *
               get(v_minus_int);

  // in 2+ dimensions the calculation is done without any memory allocations
  DataVector temp_2{};
  if constexpr (Dim > 1) {
    temp_2.set_data_ref(make_not_null(&get<Dim - 2>(*phi_boundary_correction)));
  }
  temp_2 = step_function(-get<2>(char_speeds_int)) * get<2>(char_speeds_int) *
               get(v_plus_int) +
           step_function(get<3>(char_speeds_ext)) * get<3>(char_speeds_ext) *
               get(v_minus_ext_used);
  get(*pi_boundary_correction) =
      0.5 * (temp_1 - temp_2) + get(gamma2_int) * get(*psi_boundary_correction);

  temp_1 = -0.5 * (temp_1 + temp_2);
  for (size_t i = 0; i < Dim; ++i) {
    phi_boundary_correction->get(i) =
        temp_1 * interface_unit_normal_int.get(i) -
        step_function(get<1>(char_speeds_ext)) * get<1>(char_speeds_ext) *
            v_zero_ext_used.get(i) -
        step_function(-get<1>(char_speeds_int)) * get<1>(char_speeds_int) *
            v_zero_int.get(i);
  }
}

template <size_t Dim>
bool operator==(const UpwindPenalty<Dim>& /*lhs*/,
                const UpwindPenalty<Dim>& /*rhs*/) {
  return true;
}

template <size_t Dim>
bool operator!=(const UpwindPenalty<Dim>& lhs, const UpwindPenalty<Dim>& rhs) {
  return not(lhs == rhs);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID UpwindPenalty<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(_, data)                                       \
  template class UpwindPenalty<DIM(data)>;                           \
  template bool operator==(const UpwindPenalty<DIM(data)>& /*lhs*/,  \
                           const UpwindPenalty<DIM(data)>& /*rhs*/); \
  template bool operator!=(const UpwindPenalty<DIM(data)>& /*lhs*/,  \
                           const UpwindPenalty<DIM(data)>& /*rhs*/);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryCorrections
