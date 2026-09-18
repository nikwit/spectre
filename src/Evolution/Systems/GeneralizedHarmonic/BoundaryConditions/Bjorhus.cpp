// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"

#include <algorithm>
#include <optional>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gh::BoundaryConditions {

namespace detail {
ConstraintPreservingBjorhusType
convert_constraint_preserving_bjorhus_type_from_yaml(
    const Options::Option& options) {
  const auto type_read = options.parse_as<std::string>();
  if (type_read == "ConstraintPreserving") {
    return ConstraintPreservingBjorhusType::ConstraintPreserving;
  } else if (type_read == "ConstraintPreservingPhysical") {
    return ConstraintPreservingBjorhusType::ConstraintPreservingPhysical;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to "
              "ConstraintPreservingBjorhusType::Type. Must "
              "be one of ConstraintPreserving or ConstraintPreservingPhysical");
}
}  // namespace detail

template <size_t Dim>
ConstraintPreservingBjorhus<Dim>::ConstraintPreservingBjorhus(
    const detail::ConstraintPreservingBjorhusType type,
    std::optional<std::unique_ptr<::MathFunction<1, Frame::Inertial>>>
        incoming_wave_profile)
    : type_(type),
      incoming_wave_profile_(incoming_wave_profile.has_value()
                                 ? std::move(incoming_wave_profile.value())
                                 : nullptr) {
  if constexpr (Dim < 3) {
    if (incoming_wave_profile_ != nullptr) {
      ERROR(
          "IncomingWaveProfile can only be used for "
          "ConstraintPreservingBjorhus in 3 spatial dimensions.");
    }
  }
}

template <size_t Dim>
ConstraintPreservingBjorhus<Dim>::ConstraintPreservingBjorhus(
    const ConstraintPreservingBjorhus& rhs)
    : BoundaryCondition<Dim>{dynamic_cast<const BoundaryCondition<Dim>&>(rhs)},
      type_(rhs.type_),
      incoming_wave_profile_(rhs.incoming_wave_profile_ != nullptr
                                 ? rhs.incoming_wave_profile_->get_clone()
                                 : nullptr) {}

template <size_t Dim>
ConstraintPreservingBjorhus<Dim>& ConstraintPreservingBjorhus<Dim>::operator=(
    const ConstraintPreservingBjorhus& rhs) {
  if (&rhs == this) {
    return *this;
  }
  type_ = rhs.type_;
  incoming_wave_profile_ = rhs.incoming_wave_profile_ != nullptr
                               ? rhs.incoming_wave_profile_->get_clone()
                               : nullptr;
  return *this;
}

template <size_t Dim>
ConstraintPreservingBjorhus<Dim>::ConstraintPreservingBjorhus(
    CkMigrateMessage* const msg)
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
ConstraintPreservingBjorhus<Dim>::get_clone() const {
  return std::make_unique<ConstraintPreservingBjorhus>(*this);
}

template <size_t Dim>
void ConstraintPreservingBjorhus<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
  p | type_;
  p | incoming_wave_profile_;
}

template <size_t Dim>
std::optional<std::string> ConstraintPreservingBjorhus<Dim>::dg_time_derivative(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        dt_spacetime_metric_correction,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        dt_pi_correction,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        dt_phi_correction,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& /*normal_vector*/,
    // c.f. dg_interior_evolved_variables_tags
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    // c.f. dg_interior_temporary_tags
    const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
    const tnsr::AA<DataVector, Dim, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, Dim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, Dim, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, Dim, Frame::Inertial>&
        spacetime_deriv_gauge_source,
    // c.f. dg_interior_dt_vars_tags
    const tnsr::aa<DataVector, Dim, Frame::Inertial>&
        logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& logical_dt_phi,
    // c.f. dg_interior_deriv_vars_tags
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi,
    const double time) const {
  Bjorhus::IntermediateVariables<Dim> vars{get_size(get<0>(normal_covector))};
  Bjorhus::compute_intermediate_variables(
      make_not_null(&vars), face_mesh_velocity, normal_covector,
      spacetime_metric, pi, phi, gamma1, gamma2, lapse, shift,
      inverse_spacetime_metric, spacetime_unit_normal_vector,
      three_index_constraint, gauge_source, spacetime_deriv_gauge_source,
      logical_dt_spacetime_metric, logical_dt_pi, logical_dt_phi,
      d_spacetime_metric, d_pi, d_phi);

  // If no point on the boundary has any incoming characteristic, return here
  if (Bjorhus::min_characteristic_speed(vars.char_speeds) >= 0.) {
    std::fill(dt_spacetime_metric_correction->begin(),
              dt_spacetime_metric_correction->end(), 0.);
    std::fill(dt_pi_correction->begin(), dt_pi_correction->end(), 0.);
    std::fill(dt_phi_correction->begin(), dt_phi_correction->end(), 0.);
    return {};
  }

  auto bc_dt_v_psi =
      make_with_value<tnsr::aa<DataVector, Dim, Frame::Inertial>>(get(gamma2),
                                                                  0.);
  auto bc_dt_v_zero =
      make_with_value<tnsr::iaa<DataVector, Dim, Frame::Inertial>>(get(gamma2),
                                                                   0.);
  auto bc_dt_v_plus =
      make_with_value<tnsr::aa<DataVector, Dim, Frame::Inertial>>(get(gamma2),
                                                                  0.);
  auto bc_dt_v_minus =
      make_with_value<tnsr::aa<DataVector, Dim, Frame::Inertial>>(get(gamma2),
                                                                  0.);

  Bjorhus::constraint_preserving_corrections_dt_v_psi(
      make_not_null(&bc_dt_v_psi), vars.unit_interface_normal_vector,
      three_index_constraint, vars.char_speeds);

  Bjorhus::constraint_preserving_corrections_dt_v_zero(
      make_not_null(&bc_dt_v_zero), vars.unit_interface_normal_vector,
      vars.four_index_constraint, vars.char_speeds);

  // In order to set dt<V+> = 0, the correction term returned here must be
  // b_correction = -1*existing(dt<V+>), such that
  // final(dt<V+>) = existing(dt<V+>) + b_correction = 0
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = a; b <= Dim; ++b) {
      bc_dt_v_plus.get(a, b) = -vars.char_projected_rhs_dt_v_plus.get(a, b);
    }
  }

  if (type_ == detail::ConstraintPreservingBjorhusType::ConstraintPreserving) {
    Bjorhus::constraint_preserving_gauge_corrections_dt_v_minus(
        make_not_null(&bc_dt_v_minus), gamma2, coords,
        vars.incoming_null_one_form, vars.outgoing_null_one_form,
        vars.incoming_null_vector, vars.outgoing_null_vector,
        vars.projection_ab, vars.projection_Ab, vars.projection_AB,
        vars.char_projected_rhs_dt_v_psi, vars.char_projected_rhs_dt_v_minus,
        vars.constraint_char_zero_plus, vars.constraint_char_zero_minus,
        vars.char_speeds);
  } else if (type_ == detail::ConstraintPreservingBjorhusType::
                          ConstraintPreservingPhysical) {
    Bjorhus::constraint_preserving_gauge_physical_corrections_dt_v_minus(
        make_not_null(&bc_dt_v_minus), gamma2, coords, time, normal_covector,
        vars.unit_interface_normal_vector, spacetime_unit_normal_vector,
        vars.incoming_null_one_form, vars.outgoing_null_one_form,
        vars.incoming_null_vector, vars.outgoing_null_vector,
        vars.projection_ab, vars.projection_Ab, vars.projection_AB,
        vars.inverse_spatial_metric, vars.extrinsic_curvature, spacetime_metric,
        inverse_spacetime_metric, three_index_constraint,
        vars.char_projected_rhs_dt_v_psi, vars.char_projected_rhs_dt_v_minus,
        vars.constraint_char_zero_plus, vars.constraint_char_zero_minus, phi,
        d_phi, d_pi, vars.char_speeds, incoming_wave_profile_.get());
  } else {
    ERROR(
        "Failed to set dtVMinus. Input option must be one of "
        "ConstraintPreserving or ConstraintPreservingPhysical");
  }

  Bjorhus::project_corrections_onto_evolved_variables(
      dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
      make_not_null(&bc_dt_v_psi), make_not_null(&bc_dt_v_zero),
      make_not_null(&bc_dt_v_plus), make_not_null(&bc_dt_v_minus),
      vars.char_speeds, gamma2, normal_covector);

  if (face_mesh_velocity.has_value()) {
    const auto radial_mesh_velocity =
        get(dot_product(normal_covector, *face_mesh_velocity));
    // we use 1e-10 instead of 0 below to allow for purely tangentially
    // moving grids, eg a rotating sphere, with some leeway for
    // floating-point errors.
    if (max(radial_mesh_velocity) > 1.e-10) {
      return {
          "We found the radial mesh velocity points in the direction "
          "of the outward normal, i.e. we possibly have an expanding "
          "domain. Its unclear if proper boundary conditions are "
          "imposed in this case."};
    }
  }

  return {};
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID ConstraintPreservingBjorhus<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class ConstraintPreservingBjorhus<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace gh::BoundaryConditions
