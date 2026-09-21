// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"

#include <algorithm>
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
SectorImposition convert_sector_imposition_from_yaml(
    const Options::Option& options) {
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

std::ostream& operator<<(std::ostream& os, const SectorImposition imposition) {
  switch (imposition) {
    case SectorImposition::Bjorhus:
      return os << "Bjorhus";
    case SectorImposition::Frozen:
      return os << "Frozen";
    case SectorImposition::SommerfeldAbsorbing:
      return os << "SommerfeldAbsorbing";
    case SectorImposition::SommerfeldOutgoing:
      return os << "SommerfeldOutgoing";
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
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const ElementId<Dim>& element_id, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  for (const auto& [name, excision_sphere] : excision_spheres) {
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
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial>& unit_normal_covector) {
  // Lapse, shift and face data are not consumed by the type-D model
  const Scalar<DataVector> unused_lapse{};
  const tnsr::I<DataVector, 3, Frame::Inertial> unused_shift{};
  return worldtube::evaluate_matching(worldtube::PhysicalModel::TypeD,
                                      std::nullopt, electric, magnetic,
                                      spatial_metric, unit_normal_covector,
                                      unused_lapse, unused_shift, nullptr)
      .incoming_mode;
}
}  // namespace detail

namespace {
// Resolve the aggregate-or-per-field option into one imposition per field.
detail::PerFieldConstraintSectors resolve_constraint_sectors(
    const std::variant<detail::SectorImposition,
                       detail::PerFieldConstraintSectors>& sector) {
  if (std::holds_alternative<detail::SectorImposition>(sector)) {
    const auto imposition = std::get<detail::SectorImposition>(sector);
    return {imposition, imposition, imposition};
  }
  return std::get<detail::PerFieldConstraintSectors>(sector);
}
}  // namespace

template <size_t Dim>
WorldtubeTypeD<Dim>::WorldtubeTypeD(
    const std::variant<detail::SectorImposition,
                       detail::PerFieldConstraintSectors>
        constraint_preserving_sector,
    const detail::SectorImposition physical_sector,
    const detail::SectorImposition gauge_sector,
    const detail::PhysicalModel physical_model,
    const std::optional<double> mass, const Options::Context& context)
    : constraint_v_psi_(
          resolve_constraint_sectors(constraint_preserving_sector).v_psi),
      constraint_v_zero_(
          resolve_constraint_sectors(constraint_preserving_sector).v_zero),
      constraint_v_minus_(
          resolve_constraint_sectors(constraint_preserving_sector).v_minus),
      physical_sector_(physical_sector),
      gauge_sector_(gauge_sector),
      physical_model_(physical_model),
      mass_(mass) {
  if (gauge_sector_ == detail::SectorImposition::Bjorhus) {
    PARSE_ERROR(context,
                "GaugeSector: Bjorhus is not implemented. A relaxation towards "
                "a model needs the time derivative of the model's u^-, which "
                "no model supplies. Use Frozen for the no-condition control, "
                "or SommerfeldAbsorbing / SommerfeldOutgoing for a model-free "
                "time-derivative condition.");
  }
  const auto reject_sommerfeld = [&context](
                                     const detail::SectorImposition imposition,
                                     const std::string& option_name) {
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
  if (physical_model_ == detail::PhysicalModel::Quadrupole) {
    if (not mass_.has_value()) {
      PARSE_ERROR(context,
                  "PhysicalModel: Quadrupole needs the mass of the excised "
                  "hole. Set Mass.");
    }
    if (*mass_ <= 0.) {
      PARSE_ERROR(context, "Mass must be positive, not " << *mass_);
    }
  } else if (mass_.has_value()) {
    PARSE_ERROR(context,
                "Mass is only used by PhysicalModel: Quadrupole, but "
                "PhysicalModel is "
                    << physical_model_ << ". Set Mass: None.");
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
WorldtubeTypeD<Dim>::WorldtubeTypeD(CkMigrateMessage* const msg)
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubeTypeD<Dim>::get_clone() const {
  return std::make_unique<WorldtubeTypeD>(*this);
}

template <size_t Dim>
void WorldtubeTypeD<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
  p | constraint_v_psi_;
  p | constraint_v_zero_;
  p | constraint_v_minus_;
  p | physical_sector_;
  p | gauge_sector_;
  p | physical_model_;
  p | mass_;
}

template <size_t Dim>
std::optional<std::string> WorldtubeTypeD<Dim>::dg_time_derivative(
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
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
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
    const tnsr::aa<DataVector, Dim, Frame::Inertial>&
        logical_dt_spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& logical_dt_pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& logical_dt_phi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi,
    const double time, const Domain<Dim>& domain, const Element<Dim>& element,
    const domain::FunctionsOfTimeMap& functions_of_time,
    const worldtube::KretschmannFaceData<Dim>& face_data) const {
  const size_t num_points = get_size(get<0>(normal_covector));
  Bjorhus::IntermediateVariables<Dim> vars{num_points};
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
                normal_covector, lapse, shift, &face_data)
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
        incoming_mode.has_value() ? &*incoming_mode : nullptr);
  }
  Bjorhus::detail::add_physical_sector_projection(
      make_not_null(&bc_dt_v_minus), minus_one, vars.projection_ab,
      vars.projection_Ab, vars.projection_AB,
      vars.char_projected_rhs_dt_v_minus);

  // Gauge sector of v^-
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

  Bjorhus::project_corrections_onto_evolved_variables(
      dt_spacetime_metric_correction, dt_pi_correction, dt_phi_correction,
      make_not_null(&bc_dt_v_psi), make_not_null(&bc_dt_v_zero),
      make_not_null(&bc_dt_v_plus), make_not_null(&bc_dt_v_minus),
      vars.char_speeds, gamma2, normal_covector);

  // No veto of a mesh velocity along the outward normal here: at an inner
  // boundary that normal points into the excision, so an excision tracking a
  // moving hole has such a velocity over half of the sphere by construction.
  // The mesh velocity is accounted for in the inertial time derivatives and
  // in the characteristic speeds.
  return {};
}

template <size_t Dim>
bool operator==(const WorldtubeTypeD<Dim>& lhs,
                const WorldtubeTypeD<Dim>& rhs) {
  return lhs.constraint_v_psi() == rhs.constraint_v_psi() and
         lhs.constraint_v_zero() == rhs.constraint_v_zero() and
         lhs.constraint_v_minus() == rhs.constraint_v_minus() and
         lhs.physical_sector() == rhs.physical_sector() and
         lhs.gauge_sector() == rhs.gauge_sector() and
         lhs.physical_model() == rhs.physical_model() and
         lhs.mass() == rhs.mass();
}

template <size_t Dim>
bool operator!=(const WorldtubeTypeD<Dim>& lhs,
                const WorldtubeTypeD<Dim>& rhs) {
  return not(lhs == rhs);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID WorldtubeTypeD<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                          \
  template class WorldtubeTypeD<DIM(data)>;                             \
  template bool operator==(const WorldtubeTypeD<DIM(data)>& lhs,        \
                           const WorldtubeTypeD<DIM(data)>& rhs);       \
  template bool operator!=(const WorldtubeTypeD<DIM(data)>& lhs,        \
                           const WorldtubeTypeD<DIM(data)>& rhs);       \
  template tnsr::I<double, DIM(data), Frame::Inertial>                  \
  detail::excision_sphere_center(                                       \
      const std::unordered_map<std::string, ExcisionSphere<DIM(data)>>& \
          excision_spheres,                                             \
      const ElementId<DIM(data)>& element_id, double time,              \
      const domain::FunctionsOfTimeMap& functions_of_time);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace gh::BoundaryConditions
