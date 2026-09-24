// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/AveragedUpwindPenalty.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Element.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/AlgebraicGauge.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "Options/ParseError.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gh::BoundaryConditions {
namespace detail {
SectorImposition
convert_sector_imposition_from_yaml(const Options::Option &options) {
  const auto read = options.parse_as<std::string>();
  if (read == "Bjorhus") {
    return SectorImposition::Bjorhus;
  } else if (read == "Frozen") {
    return SectorImposition::Frozen;
  } else if (read == "SommerfeldAbsorbing") {
    return SectorImposition::SommerfeldAbsorbing;
  } else if (read == "SommerfeldOutgoing") {
    return SectorImposition::SommerfeldOutgoing;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to a sector imposition. Must be "
              "one of Bjorhus, Frozen, SommerfeldAbsorbing or "
              "SommerfeldOutgoing. The two Sommerfeld conditions are only "
              "available for the gauge sector.");
}

bool is_sommerfeld(const SectorImposition imposition) {
  return imposition == SectorImposition::SommerfeldAbsorbing or
         imposition == SectorImposition::SommerfeldOutgoing;
}

double sommerfeld_one_over_r_sign(const SectorImposition imposition) {
  ASSERT(is_sommerfeld(imposition),
         "sommerfeld_one_over_r_sign called on a non-Sommerfeld imposition.");
  // The coefficient is gamma2 - c/r. At an inner boundary the domain's outward
  // normal points into the hole, so d_n = -d_r, and the ingoing condition
  // (d_t - d_r - 1/r) f = 0 carries the opposite 1/r sign to the outgoing
  // outer-boundary form (d_t + d_r + 1/r) f = 0.
  return imposition == SectorImposition::SommerfeldAbsorbing ? -1.0 : 1.0;
}

std::ostream &operator<<(std::ostream &os, const SectorImposition imposition) {
  switch (imposition) {
  case SectorImposition::Bjorhus:
    return os << "Bjorhus";
  case SectorImposition::Frozen:
    return os << "Frozen";
  case SectorImposition::SommerfeldAbsorbing:
    return os << "SommerfeldAbsorbing";
  case SectorImposition::SommerfeldOutgoing:
    return os << "SommerfeldOutgoing";
  case SectorImposition::Algebraic:
    return os << "Algebraic";
  case SectorImposition::OutgoingDriven:
    return os << "OutgoingDriven";
  case SectorImposition::SchwarzschildReference:
    return os << "SchwarzschildReference";
  case SectorImposition::RadialResponse:
    return os << "RadialResponse";
  default:
    ERROR("Unknown SectorImposition");
  }
}

PerFieldConstraintSectors::PerFieldConstraintSectors(
    const SectorImposition v_psi_in, const SectorImposition v_zero_in,
    const SectorImposition v_minus_in)
    : v_psi(v_psi_in), v_zero(v_zero_in), v_minus(v_minus_in) {}

template <size_t Dim>
tnsr::I<double, Dim, Frame::Inertial> excision_sphere_center(
    const std::unordered_map<std::string, ExcisionSphere<Dim>>
        &excision_spheres,
    const ElementId<Dim> &element_id, const double time,
    const domain::FunctionsOfTimeMap &functions_of_time) {
  for (const auto &[name, excision_sphere] : excision_spheres) {
    (void)name;
    if (not excision_sphere.abutting_direction(element_id).has_value()) {
      continue;
    }
    if (excision_sphere.is_time_dependent()) {
      return excision_sphere.moving_mesh_grid_to_inertial_map()(
          excision_sphere.center(), time, functions_of_time);
    }
    tnsr::I<double, Dim, Frame::Inertial> center{};
    for (size_t i = 0; i < Dim; ++i) {
      center.get(i) = excision_sphere.center().get(i);
    }
    return center;
  }
  ERROR("The element "
        << element_id
        << " abuts none of the domain's excision spheres, but the worldtube "
           "boundary condition's Sommerfeld gauge sector needs the center of "
           "the excision sphere it is applied on. The domain has "
        << excision_spheres.size() << " excision sphere(s).");
}

tnsr::ii<DataVector, 3, Frame::Inertial> type_d_incoming_mode(
    const tnsr::ii<DataVector, 3, Frame::Inertial> &electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial> &magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial> &spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial> &unit_normal_covector) {
  // Lapse, shift and face data are not consumed by the type-D model
  const Scalar<DataVector> unused_lapse{};
  const tnsr::I<DataVector, 3, Frame::Inertial> unused_shift{};
  return worldtube::evaluate_matching(worldtube::PhysicalModel::TypeD,
                                      std::nullopt, electric, magnetic,
                                      spatial_metric, unit_normal_covector,
                                      unused_lapse, unused_shift, nullptr)
      .incoming_mode;
}
} // namespace detail

namespace {
// Resolve the aggregate-or-per-field option into one imposition per field.
detail::PerFieldConstraintSectors resolve_constraint_sectors(
    const std::variant<detail::SectorImposition,
                       detail::PerFieldConstraintSectors> &sector) {
  if (std::holds_alternative<detail::SectorImposition>(sector)) {
    const auto imposition = std::get<detail::SectorImposition>(sector);
    return {imposition, imposition, imposition};
  }
  return std::get<detail::PerFieldConstraintSectors>(sector);
}
} // namespace

template <size_t Dim>
WorldtubeTypeD<Dim>::WorldtubeTypeD(
    const std::variant<detail::SectorImposition,
                       detail::PerFieldConstraintSectors>
        constraint_preserving_sector,
    const detail::SectorImposition physical_sector,
    const std::variant<detail::SectorImposition, detail::AlgebraicGauge,
                       detail::OutgoingDrivenGauge,
                       detail::SchwarzschildReferenceGauge, worldtube::ReferenceReplayGauge,
                   worldtube::RadialResponseGauge>
        gauge_sector,
    const detail::PhysicalModel physical_model,
    const std::optional<double> mass,
    const std::optional<double> moment_relaxation_time,
    const Options::Context &context)
    : constraint_v_psi_(
          resolve_constraint_sectors(constraint_preserving_sector).v_psi),
      constraint_v_zero_(
          resolve_constraint_sectors(constraint_preserving_sector).v_zero),
      constraint_v_minus_(
          resolve_constraint_sectors(constraint_preserving_sector).v_minus),
      physical_sector_(physical_sector),
      gauge_sector_(
          std::holds_alternative<worldtube::RadialResponseGauge>(gauge_sector)
              ? detail::SectorImposition::RadialResponse
              : std::holds_alternative<worldtube::ReferenceReplayGauge>(gauge_sector)
              ? detail::SectorImposition::Algebraic
              : std::holds_alternative<detail::AlgebraicGauge>(gauge_sector)
              ? detail::SectorImposition::Algebraic
              : (std::holds_alternative<detail::OutgoingDrivenGauge>(
                     gauge_sector)
                     ? detail::SectorImposition::OutgoingDriven
                     : (std::holds_alternative<
                            detail::SchwarzschildReferenceGauge>(gauge_sector)
                            ? detail::SectorImposition::SchwarzschildReference
                            : std::get<detail::SectorImposition>(
                                  gauge_sector)))),
      algebraic_gauge_rate_(
          std::holds_alternative<detail::AlgebraicGauge>(gauge_sector)
              ? std::get<detail::AlgebraicGauge>(gauge_sector).rate
              : 0.),
      outgoing_gauge_rate_(
          std::holds_alternative<detail::OutgoingDrivenGauge>(gauge_sector)
              ? std::get<detail::OutgoingDrivenGauge>(gauge_sector).rate
              : 0.),
      radial_response_(std::holds_alternative<worldtube::RadialResponseGauge>(gauge_sector)
          ? std::optional{std::get<worldtube::RadialResponseGauge>(gauge_sector).coefficients}
          : std::nullopt),
      schwarzschild_reference_(
          std::holds_alternative<detail::SchwarzschildReferenceGauge>(
              gauge_sector)
              ? std::optional{std::get<detail::SchwarzschildReferenceGauge>(
                                  gauge_sector)
                                  .parameters}
              : std::nullopt),
      physical_model_(physical_model), mass_(mass),
      moment_relaxation_time_(moment_relaxation_time) {
  if (std::holds_alternative<worldtube::ReferenceReplayGauge>(gauge_sector)) {
    face_replay_ = std::get<worldtube::ReferenceReplayGauge>(gauge_sector).parameters;
    schwarzschild_reference_ = face_replay_->reference;
    gauge_sector_ = schwarzschild_reference_.has_value() ? detail::SectorImposition::SchwarzschildReference : detail::SectorImposition::Algebraic;
    if constexpr (Dim != 3) { PARSE_ERROR(context,"Face replay requires three dimensions."); }
  }
  if constexpr (Dim != 3) {
    if (schwarzschild_reference_.has_value()) {
      PARSE_ERROR(context, "SchwarzschildReference requires three dimensions.");
    }
  }
  if (radial_response_.has_value()) {
    if constexpr (Dim != 3) { PARSE_ERROR(context, "RadialResponse requires three dimensions."); }
    for (const auto k : *radial_response_) {
      if (not std::isfinite(k)) { PARSE_ERROR(context, "RadialResponse coefficients must be finite."); }
    }
  }
  if (not std::isfinite(outgoing_gauge_rate_) or outgoing_gauge_rate_ < 0.) {
    PARSE_ERROR(context,
                "OutgoingDriven gauge rate must be finite and nonnegative.");
  }
  if (not std::isfinite(algebraic_gauge_rate_) or algebraic_gauge_rate_ < 0.) {
    PARSE_ERROR(context,
                "Algebraic gauge rate must be finite and nonnegative.");
  }
  if (gauge_sector_ == detail::SectorImposition::Bjorhus) {
    PARSE_ERROR(
        context,
        "GaugeSector: Bjorhus is not implemented. A relaxation towards "
        "a model needs the time derivative of the model's u^-, which "
        "no model supplies. Use {Algebraic: rate} for homogeneous gauge "
        "data, Frozen for a zero projected RHS, "
        "or SommerfeldAbsorbing / SommerfeldOutgoing for a model-free "
        "time-derivative condition.");
  }
  const auto reject_sommerfeld = [&context](
                                     const detail::SectorImposition imposition,
                                     const std::string &option_name) {
    if (imposition == detail::SectorImposition::Algebraic or
        imposition == detail::SectorImposition::OutgoingDriven or
        imposition == detail::SectorImposition::SchwarzschildReference or
        imposition == detail::SectorImposition::RadialResponse) {
      PARSE_ERROR(context,
                  option_name
                      << ": Algebraic applies only to the gauge sector.");
    }
    if (detail::is_sommerfeld(imposition)) {
      PARSE_ERROR(context,
                  option_name
                      << ": the Sommerfeld conditions apply only to the gauge "
                         "sector, which is the sector whose characteristic "
                         "they describe.");
    }
  };
  reject_sommerfeld(constraint_v_psi_, "ConstraintPreservingSector VPsi");
  reject_sommerfeld(constraint_v_zero_, "ConstraintPreservingSector VZero");
  reject_sommerfeld(constraint_v_minus_, "ConstraintPreservingSector VMinus");
  reject_sommerfeld(physical_sector_, "PhysicalSector");
  if (physical_model_ != detail::PhysicalModel::None and
      physical_sector_ != detail::SectorImposition::Bjorhus) {
    PARSE_ERROR(context,
                "PhysicalModel: "
                    << physical_model_
                    << " enters through the Bjorhus term of the physical "
                       "sector, but PhysicalSector is "
                    << physical_sector_ << ". Use PhysicalSector: Bjorhus.");
  }
  if (worldtube::is_order_two(physical_model_)) {
    if (not mass_.has_value()) {
      PARSE_ERROR(context, "PhysicalModel: "
                               << physical_model_
                               << " needs the mass of the excised hole. Set "
                                  "Mass.");
    }
    if (*mass_ <= 0.) {
      PARSE_ERROR(context, "Mass must be positive, not " << *mass_);
    }
  } else if (mass_.has_value()) {
    PARSE_ERROR(context,
                "Mass is only used by the order-two models Quadrupole and "
                "QuadrupoleCoulomb, but PhysicalModel is "
                    << physical_model_ << ". Set Mass: None.");
  }
  if (moment_relaxation_time_.has_value()) {
    if (not worldtube::is_order_two(physical_model_)) {
      PARSE_ERROR(context,
                  "MomentRelaxationTime is only used by the order-two models "
                  "Quadrupole and QuadrupoleCoulomb, but PhysicalModel is "
                      << physical_model_
                      << ". Set MomentRelaxationTime: None.");
    }
    if (*moment_relaxation_time_ <= 0.) {
      PARSE_ERROR(context, "MomentRelaxationTime must be positive, not "
                               << *moment_relaxation_time_);
    }
  }
  if constexpr (Dim != 3) {
    if (physical_model_ != detail::PhysicalModel::None) {
      PARSE_ERROR(context,
                  "PhysicalModel: " << physical_model_
                                    << " is only implemented in 3 dimensions.");
    }
  }
}

template <size_t Dim>
WorldtubeTypeD<Dim>::WorldtubeTypeD(CkMigrateMessage *const msg)
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubeTypeD<Dim>::get_clone() const {
  return std::make_unique<WorldtubeTypeD>(*this);
}

template <size_t Dim> void WorldtubeTypeD<Dim>::pup(PUP::er &p) {
  BoundaryCondition<Dim>::pup(p);
  p | constraint_v_psi_;
  p | constraint_v_zero_;
  p | constraint_v_minus_;
  p | physical_sector_;
  p | gauge_sector_;
  p | algebraic_gauge_rate_;
  p | outgoing_gauge_rate_;
  p | radial_response_;
  p | schwarzschild_reference_;
  p | physical_model_;
  p | mass_;
  p | moment_relaxation_time_;
  p | face_replay_;
}

template <size_t Dim>
std::optional<std::string> WorldtubeTypeD<Dim>::dg_time_derivative(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial> *>
        dt_spacetime_metric_correction,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial> *>
        dt_pi_correction,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial> *>
        dt_phi_correction,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>
        &face_mesh_velocity,
    const tnsr::i<DataVector, Dim, Frame::Inertial> &normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial> & /*normal_vector*/,
    const tnsr::aa<DataVector, Dim, Frame::Inertial> &spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial> &pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial> &phi,
    const tnsr::I<DataVector, Dim, Frame::Inertial> &coords,
    const Scalar<DataVector> &gamma1, const Scalar<DataVector> &gamma2,
    const Scalar<DataVector> &lapse,
    const tnsr::I<DataVector, Dim, Frame::Inertial> &shift,
    const tnsr::AA<DataVector, Dim, Frame::Inertial> &inverse_spacetime_metric,
    const tnsr::A<DataVector, Dim, Frame::Inertial>
        &spacetime_unit_normal_vector,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial> &three_index_constraint,
    const tnsr::a<DataVector, Dim, Frame::Inertial> &gauge_source,
    const tnsr::ab<DataVector, Dim, Frame::Inertial>
        &spacetime_deriv_gauge_source,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>
        &logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial> &logical_dt_pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial> &logical_dt_phi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial> &d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial> &d_pi,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial> &d_phi,
    const double time, const Domain<Dim> &domain, const Element<Dim> &element,
    const domain::FunctionsOfTimeMap &functions_of_time,
    const worldtube::KretschmannFaceData<Dim> &face_data,
    const ConstraintDamping::DampingFunction<Dim, Frame::Grid> &damping_gamma2)
    const {
  if (gauge_sector_ == detail::SectorImposition::Algebraic or
      gauge_sector_ == detail::SectorImposition::OutgoingDriven or
      gauge_sector_ == detail::SectorImposition::SchwarzschildReference or
      gauge_sector_ == detail::SectorImposition::RadialResponse) {
    if (domain.blocks()[element.id().block_id()].is_time_dependent()) {
      return "Algebraic gauge currently requires a static mesh.";
    }
    if (dynamic_cast<const ConstraintDamping::Constant<Dim, Frame::Grid> *>(
            &damping_gamma2) == nullptr and
        dynamic_cast<
            const ConstraintDamping::GaussianPlusConstant<Dim, Frame::Grid> *>(
            &damping_gamma2) == nullptr) {
      return "Algebraic gauge requires time-independent gamma2 (Constant or "
             "GaussianPlusConstant).";
    }
  }
  const size_t num_points = get_size(get<0>(normal_covector));
  if (gauge_sector_ == detail::SectorImposition::OutgoingDriven and
      (not face_data.initial_gauge_difference.has_value() or
       get<0>(*face_data.initial_gauge_difference).size() != num_points)) {
    return "OutgoingDriven requires initial gauge face data at fixed "
           "resolution.";
  }
  if (radial_response_.has_value() and
      (not face_data.radial_gauge.has_value() or
       get<0>(face_data.radial_gauge->q).size() != num_points or
       face_data.radial_gauge->time != time)) {
    return "RadialResponse requires current fixed-resolution radial gauge face data.";
  }
  Bjorhus::IntermediateVariables<Dim> vars{num_points};
  Bjorhus::compute_intermediate_variables(
      make_not_null(&vars), face_mesh_velocity, normal_covector,
      spacetime_metric, pi, phi, gamma1, gamma2, lapse, shift,
      inverse_spacetime_metric, spacetime_unit_normal_vector,
      three_index_constraint, gauge_source, spacetime_deriv_gauge_source,
      logical_dt_spacetime_metric, logical_dt_pi, logical_dt_phi,
      d_spacetime_metric, d_pi, d_phi);

  if (face_replay_.has_value() and face_replay_->record and
      Bjorhus::min_characteristic_speed(vars.char_speeds) < 0.) {
    return "The recording donor must have pure outflow at its excision boundary.";
  }
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

  // Each sector is imposed independently. The three sector projections of u^-
  // partition the identity (see BjorhusImpl), so the total correction to
  // dt v^- splits by sector:
  //
  //   Bjorhus -> that sector's condition terms alone. The `add_*_terms`
  //              helpers add P_X(char_projected_rhs) + condition, and the
  //              P_X(char_projected_rhs) piece is exactly what unfreezes the
  //              sector, so subtracting it leaves the condition.
  //   Frozen  -> -P_X(char_projected_rhs), i.e. dt u^-|X = 0.
  //
  // With Bjorhus for the constraint-preserving and physical sectors,
  // SommerfeldOutgoing for the gauge sector and no physical model, this
  // reproduces ConstraintPreservingBjorhus with Type:
  // ConstraintPreservingPhysical on an excision centered at the origin.
  const DataVector minus_one(num_points, -1.0);

  if (constraint_v_psi_ == detail::SectorImposition::Bjorhus) {
    Bjorhus::constraint_preserving_corrections_dt_v_psi(
        make_not_null(&bc_dt_v_psi), vars.unit_interface_normal_vector,
        three_index_constraint, vars.char_speeds);
  } else {
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = a; b <= Dim; ++b) {
        bc_dt_v_psi.get(a, b) = -vars.char_projected_rhs_dt_v_psi.get(a, b);
      }
    }
  }

  if (constraint_v_zero_ == detail::SectorImposition::Bjorhus) {
    Bjorhus::constraint_preserving_corrections_dt_v_zero(
        make_not_null(&bc_dt_v_zero), vars.unit_interface_normal_vector,
        vars.four_index_constraint, vars.char_speeds);
  } else {
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = a; b <= Dim; ++b) {
        for (size_t i = 0; i < Dim; ++i) {
          bc_dt_v_zero.get(i, a, b) =
              -vars.char_projected_rhs_dt_v_zero.get(i, a, b);
        }
      }
    }
  }

  // v^+ is always frozen: dt<V+> = 0 requires the correction
  // -existing(dt<V+>).
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = a; b <= Dim; ++b) {
      bc_dt_v_plus.get(a, b) = -vars.char_projected_rhs_dt_v_plus.get(a, b);
    }
  }

  // Constraint-preserving sector of v^-
  if (constraint_v_minus_ == detail::SectorImposition::Bjorhus) {
    Bjorhus::detail::add_constraint_dependent_terms_to_dt_v_minus(
        make_not_null(&bc_dt_v_minus), vars.outgoing_null_one_form,
        vars.incoming_null_vector, vars.outgoing_null_vector,
        vars.projection_ab, vars.projection_Ab, vars.projection_AB,
        vars.constraint_char_zero_plus, vars.constraint_char_zero_minus,
        vars.char_projected_rhs_dt_v_minus, vars.char_speeds);
  }
  Bjorhus::detail::add_constraint_sector_projection(
      make_not_null(&bc_dt_v_minus), minus_one, vars.outgoing_null_one_form,
      vars.incoming_null_vector, vars.projection_ab, vars.projection_Ab,
      vars.projection_AB, vars.char_projected_rhs_dt_v_minus);

  // Physical sector of v^-
  if (physical_sector_ == detail::SectorImposition::Bjorhus) {
    // The incoming Weyl mode the model supplies, if any
    std::optional<tnsr::ii<DataVector, Dim, Frame::Inertial>> incoming_mode{};
    if (physical_model_ != detail::PhysicalModel::None) {
      if constexpr (Dim == 3) {
        tnsr::ii<DataVector, Dim, Frame::Inertial> spatial_metric(num_points);
        for (size_t i = 0; i < Dim; ++i) {
          for (size_t j = i; j < Dim; ++j) {
            spatial_metric.get(i, j) = spacetime_metric.get(i + 1, j + 1);
          }
        }
        tnsr::ii<DataVector, Dim, Frame::Inertial> electric{};
        tnsr::ii<DataVector, Dim, Frame::Inertial> magnetic{};
        worldtube::weyl_electric_magnetic(
            make_not_null(&electric), make_not_null(&magnetic), phi, d_phi,
            d_pi, spacetime_unit_normal_vector, spatial_metric,
            vars.inverse_spatial_metric, vars.extrinsic_curvature,
            inverse_spacetime_metric);
        incoming_mode =
            worldtube::evaluate_matching(
                physical_model_, mass_, electric, magnetic, spatial_metric,
                normal_covector, lapse, shift, &face_data,
                // The relaxed moments maintained by UpdateKretschmannFaceData,
                // if the face has them; the instantaneous fit otherwise
                face_data.filtered_moments)
                .incoming_mode;
      } else {
        (void)face_data;
        ERROR("PhysicalModel: " << physical_model_
                                << " is only implemented in 3 dimensions.");
      }
    }
    Bjorhus::detail::add_physical_terms_to_dt_v_minus(
        make_not_null(&bc_dt_v_minus), gamma2, normal_covector,
        vars.unit_interface_normal_vector, spacetime_unit_normal_vector,
        vars.projection_ab, vars.projection_Ab, vars.projection_AB,
        vars.inverse_spatial_metric, vars.extrinsic_curvature, spacetime_metric,
        inverse_spacetime_metric, three_index_constraint,
        vars.char_projected_rhs_dt_v_minus, phi, d_phi, d_pi, vars.char_speeds,
        std::numeric_limits<double>::signaling_NaN(), nullptr,
        Bjorhus::default_incoming_wave_components,
        incoming_mode.has_value() ? &*incoming_mode : nullptr);
  }
  Bjorhus::detail::add_physical_sector_projection(
      make_not_null(&bc_dt_v_minus), minus_one, vars.projection_ab,
      vars.projection_Ab, vars.projection_AB,
      vars.char_projected_rhs_dt_v_minus);

  // Gauge sector of v^-
  if (gauge_sector_ == detail::SectorImposition::Algebraic or
      gauge_sector_ == detail::SectorImposition::OutgoingDriven or
      gauge_sector_ == detail::SectorImposition::SchwarzschildReference or
      gauge_sector_ == detail::SectorImposition::RadialResponse) {
    auto full_dt_metric = vars.dt_spacetime_metric;
    // A correction of the zero-speed metric characteristic, if selected,
    // also changes the time derivatives of the characteristic basis.
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = a; b <= Dim; ++b) {
        for (size_t p = 0; p < num_points; ++p) {
          if (vars.char_speeds[0][p] <= 0.) {
            full_dt_metric.get(a, b)[p] += bc_dt_v_psi.get(a, b)[p];
          }
        }
      }
    }
    std::optional<tnsr::iaa<DataVector, Dim>> reference_phi{};
    if constexpr (Dim == 3) {
      if (schwarzschild_reference_.has_value()) {
        const auto center = detail::excision_sphere_center(
            domain.excision_spheres(), element.id(), time, functions_of_time);
        reference_phi = detail::schwarzschild_reference_fields(
                            *schwarzschild_reference_, coords, center)
                            .second;
      }
    }
    const auto target = detail::algebraic_gauge_rhs(
        spacetime_metric, pi, phi, gamma2, vars.inverse_spatial_metric,
        inverse_spacetime_metric, spacetime_unit_normal_vector,
        vars.unit_interface_normal_vector, vars.incoming_null_one_form,
        vars.outgoing_null_vector, full_dt_metric, algebraic_gauge_rate_,
        gauge_sector_ == detail::SectorImposition::OutgoingDriven
            ? &*face_data.initial_gauge_difference
            : nullptr,
        outgoing_gauge_rate_,
        reference_phi.has_value() ? &*reference_phi : nullptr,
        schwarzschild_reference_.has_value()
            ? 1. / schwarzschild_reference_->values[8]
            : 0.);
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = a; b <= Dim; ++b) {
        bc_dt_v_minus.get(a, b) += target.get(a, b);
      }
    }
  }
  if (detail::is_sommerfeld(gauge_sector_)) {
    // 1/r is the distance from the center of the excision sphere, i.e. the
    // local worldtube radius, not the distance from the origin of the
    // inertial coordinates that the outer-boundary condition uses.
    const tnsr::I<double, Dim, Frame::Inertial> center =
        detail::excision_sphere_center(domain.excision_spheres(), element.id(),
                                       time, functions_of_time);
    DataVector radius_squared(num_points, 0.);
    for (size_t i = 0; i < Dim; ++i) {
      radius_squared += square(coords.get(i) - center.get(i));
    }
    const DataVector gauge_coefficient =
        get(gamma2) - detail::sommerfeld_one_over_r_sign(gauge_sector_) /
                          sqrt(radius_squared);
    Bjorhus::detail::add_gauge_sector_terms_to_dt_v_minus(
        make_not_null(&bc_dt_v_minus), gauge_coefficient,
        vars.incoming_null_one_form, vars.outgoing_null_one_form,
        vars.incoming_null_vector, vars.outgoing_null_vector,
        vars.projection_Ab, vars.char_projected_rhs_dt_v_psi);
  }
  Bjorhus::detail::add_gauge_sector_projection(
      make_not_null(&bc_dt_v_minus), minus_one, vars.incoming_null_one_form,
      vars.outgoing_null_one_form, vars.incoming_null_vector,
      vars.outgoing_null_vector, vars.projection_Ab,
      vars.char_projected_rhs_dt_v_minus);

  if (radial_response_.has_value()) {
    if constexpr (Dim == 3) {
      // The preceding algebraic branch gives the existing full-basis
      // zero-rate radiation condition. Replace only its t/r gauge
      // contractions; its tangential contractions stay unchanged.
      const auto& rd = *face_data.radial_gauge;
      DataVector speed(num_points, 0.), delta_r(num_points, 0.);
      DataVector dr_r(num_points, 0.), old_r(num_points, 0.);
      tnsr::a<DataVector, 3> old_q(num_points, 0.);
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = 0; b < 4; ++b) {
          old_q.get(a) += vars.outgoing_null_vector.get(b) * bc_dt_v_minus.get(a, b);
        }
      }
      for (size_t i = 0; i < 3; ++i) {
        speed -= rd.radial_direction.get(i) *
                 (shift.get(i) + get(lapse) * vars.unit_interface_normal_vector.get(i));
        delta_r += rd.radial_direction.get(i) * (rd.q.get(i+1) - rd.initial_q.get(i+1));
        dr_r += rd.radial_direction.get(i) * (rd.dr_q.get(i+1) - rd.initial_dr_q.get(i+1));
        old_r += rd.radial_direction.get(i) * old_q.get(i+1);
      }
      tnsr::a<DataVector, 3> extra(num_points, 0.);
      get<0>(extra) = speed * (get<0>(rd.dr_q) - get<0>(rd.initial_dr_q) -
                              (*radial_response_)[0] * (get<0>(rd.q) - get<0>(rd.initial_q))) - get<0>(old_q);
      for (size_t i = 0; i < 3; ++i) {
        extra.get(i+1) = rd.radial_direction.get(i) *
                        (speed * (dr_r - (*radial_response_)[1] * delta_r) - old_r);
      }
      DataVector ell_extra(num_points, 0.);
      for (size_t a = 0; a < 4; ++a) { ell_extra += vars.outgoing_null_vector.get(a) * extra.get(a); }
      for (size_t a = 0; a < 4; ++a) {
        for (size_t b = a; b < 4; ++b) {
          bc_dt_v_minus.get(a,b) -= vars.incoming_null_one_form.get(a) * extra.get(b) +
              vars.incoming_null_one_form.get(b) * extra.get(a) +
              vars.incoming_null_one_form.get(a) * vars.incoming_null_one_form.get(b) * ell_extra;
        }
      }
    }
  }

  Bjorhus::project_corrections_onto_evolved_variables(
      dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
      make_not_null(&bc_dt_v_psi), make_not_null(&bc_dt_v_zero),
      make_not_null(&bc_dt_v_plus), make_not_null(&bc_dt_v_minus),
      vars.char_speeds, gamma2, normal_covector);

  if (face_replay_.has_value() and not face_replay_->record and face_replay_->model_sector != "All") {
    if constexpr (Dim == 3) {
      if (face_mesh_velocity.has_value() or not face_data.replay.has_value())
        return "Face replay requires static-mesh donor data refreshed before this RHS.";
      for(size_t p=0;p<num_points;++p) {
        if(vars.char_speeds[0][p]<0. or vars.char_speeds[1][p]<0. or vars.char_speeds[2][p]<0. or vars.char_speeds[3][p]>=0.)
          return "Diagnostic face replay requires only u-minus to be incoming.";
      }
      const auto& donor=*face_data.replay;
      tnsr::ii<DataVector,3> spatial(num_points,0.);
      for(size_t i=0;i<3;++i)for(size_t j=i;j<3;++j)spatial.get(i,j)=donor.metric.get(i+1,j+1);
      const auto inv_donor=determinant_and_inverse(spatial).second;
      DataVector mag(num_points,0.),donor_mag(num_points,0.);
      for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j){
        mag+=donor.raw_normal.get(i)*vars.inverse_spatial_metric.get(i,j)*donor.raw_normal.get(j);
        donor_mag+=donor.raw_normal.get(i)*inv_donor.get(i,j)*donor.raw_normal.get(j);
      }
      mag=sqrt(mag);donor_mag=sqrt(donor_mag);
      tnsr::i<DataVector,3> donor_normal(num_points,0.);
      for(size_t i=0;i<3;++i)donor_normal.get(i)=-donor.raw_normal.get(i)/donor_mag;
      tnsr::aa<DataVector,3> rg(num_points,0.),rp(num_points,0.);
      tnsr::iaa<DataVector,3> rphi(num_points,0.);
      const tnsr::I<DataVector,3> velocity(num_points,0.);
      gh::BoundaryCorrections::AveragedUpwindPenalty<3>{}.dg_boundary_terms(
          make_not_null(&rg),make_not_null(&rp),make_not_null(&rphi),spacetime_metric,pi,phi,gamma1,gamma2,normal_covector,velocity,
          donor.metric,donor.pi,donor.phi,gamma1,gamma2,donor_normal,velocity,dg::Formulation::StrongInertial);
      const double nr=static_cast<double>(donor.radial_points);
      const DataVector lift=-.5*nr*(nr-1.)*mag;
      for(auto& c:rg)c*=lift;for(auto& c:rp)c*=lift;for(auto& c:rphi)c*=lift;
      // Replace the selected live u-minus projections. Keep all other
      // characteristic parts of the reference numerical interface unchanged.
      auto difference=*dt_pi_correction;
      for(size_t a=0;a<4;++a)for(size_t b=a;b<4;++b){
        difference.get(a,b)-=rp.get(a,b)+get(gamma2)*((*dt_spacetime_metric_correction).get(a,b)-rg.get(a,b));
        for(size_t i=0;i<3;++i)difference.get(a,b)-=vars.unit_interface_normal_vector.get(i)*((*dt_phi_correction).get(i,a,b)-rphi.get(i,a,b));
      }
      tnsr::aa<DataVector,3> selected(num_points,0.);const DataVector one(num_points,1.);
      if(face_replay_->model_sector=="Gauge")
        Bjorhus::detail::add_gauge_sector_projection(make_not_null(&selected),one,vars.incoming_null_one_form,vars.outgoing_null_one_form,vars.incoming_null_vector,vars.outgoing_null_vector,vars.projection_Ab,difference);
      if(face_replay_->model_sector=="CP" or face_replay_->model_sector=="CPAndPhysical")
        Bjorhus::detail::add_constraint_sector_projection(make_not_null(&selected),one,vars.outgoing_null_one_form,vars.incoming_null_vector,vars.projection_ab,vars.projection_Ab,vars.projection_AB,difference);
      if(face_replay_->model_sector=="Physical" or face_replay_->model_sector=="CPAndPhysical")
        Bjorhus::detail::add_physical_sector_projection(make_not_null(&selected),one,vars.projection_ab,vars.projection_Ab,vars.projection_AB,difference);
      for(size_t a=0;a<4;++a)for(size_t b=a;b<4;++b){
        rp.get(a,b)+=.5*selected.get(a,b);
        for(size_t i=0;i<3;++i)rphi.get(i,a,b)-=.5*normal_covector.get(i)*selected.get(a,b);
      }
      *dt_spacetime_metric_correction=std::move(rg);*dt_pi_correction=std::move(rp);*dt_phi_correction=std::move(rphi);
    }
  }
  // No veto of a mesh velocity along the outward normal here: at an inner
  // boundary that normal points into the excision, so an excision tracking a
  // moving hole has such a velocity over half of the sphere by construction.
  // The mesh velocity is accounted for in the inertial time derivatives and
  // in the characteristic speeds.
  return {};
}

template <size_t Dim>
bool operator==(const WorldtubeTypeD<Dim> &lhs,
                const WorldtubeTypeD<Dim> &rhs) {
  return lhs.face_replay() == rhs.face_replay() and lhs.constraint_v_psi() == rhs.constraint_v_psi() and
         lhs.constraint_v_zero() == rhs.constraint_v_zero() and
         lhs.constraint_v_minus() == rhs.constraint_v_minus() and
         lhs.physical_sector() == rhs.physical_sector() and
         lhs.gauge_sector() == rhs.gauge_sector() and
         lhs.algebraic_gauge_rate() == rhs.algebraic_gauge_rate() and
         lhs.outgoing_gauge_rate() == rhs.outgoing_gauge_rate() and
         lhs.radial_response() == rhs.radial_response() and
         lhs.schwarzschild_reference() == rhs.schwarzschild_reference() and
         lhs.physical_model() == rhs.physical_model() and
         lhs.mass() == rhs.mass() and
         lhs.moment_relaxation_time() == rhs.moment_relaxation_time();
}

template <size_t Dim>
bool operator!=(const WorldtubeTypeD<Dim> &lhs,
                const WorldtubeTypeD<Dim> &rhs) {
  return not(lhs == rhs);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID WorldtubeTypeD<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                                 \
  template class WorldtubeTypeD<DIM(data)>;                                    \
  template bool operator==(const WorldtubeTypeD<DIM(data)> &lhs,               \
                           const WorldtubeTypeD<DIM(data)> &rhs);              \
  template bool operator!=(const WorldtubeTypeD<DIM(data)> &lhs,               \
                           const WorldtubeTypeD<DIM(data)> &rhs);              \
  template tnsr::I<double, DIM(data), Frame::Inertial>                         \
  detail::excision_sphere_center(                                              \
      const std::unordered_map<std::string, ExcisionSphere<DIM(data)>>         \
          &excision_spheres,                                                   \
      const ElementId<DIM(data)> &element_id, double time,                     \
      const domain::FunctionsOfTimeMap &functions_of_time);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
} // namespace gh::BoundaryConditions
