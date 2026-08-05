// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"

#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Characteristics.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/InterfaceNullNormal.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::BoundaryConditions {
namespace {
double min_characteristic_speed(const std::array<DataVector, 4>& char_speeds) {
  std::array<double, 4> min_speeds{{min(char_speeds[0]), min(char_speeds[1]),
                                    min(char_speeds[2]), min(char_speeds[3])}};
  return *std::min_element(min_speeds.begin(), min_speeds.end());
}
template <typename T>
void set_bc_corr_zero_when_char_speed_is_positive(
    const gsl::not_null<T*> dt_v_corr, const DataVector& char_speed_u) {
  for (DataVector& component : *dt_v_corr) {
    for (size_t i = 0; i < component.size(); ++i) {
      if (char_speed_u[i] > 0.) {
        component[i] = 0.;
      }
    }
  }
}
}  // namespace

namespace detail {
WorldtubeTypeDType convert_worldtube_type_d_type_from_yaml(
    const Options::Option& options) {
  const auto type_read = options.parse_as<std::string>();
  if (type_read == "ConstraintPreserving") {
    return WorldtubeTypeDType::ConstraintPreserving;
  } else if (type_read == "ConstraintPreservingPhysical") {
    return WorldtubeTypeDType::ConstraintPreservingPhysical;
  } else if (type_read == "ConstraintPreservingPhysicalFrozenGauge") {
    return WorldtubeTypeDType::ConstraintPreservingPhysicalFrozenGauge;
  } else if (type_read == "ConstraintPreservingPhysicalAnalyticGhostGauge") {
    return WorldtubeTypeDType::ConstraintPreservingPhysicalAnalyticGhostGauge;
  } else if (type_read == "ConstraintPreservingPhysicalSommerfeldGhostGauge") {
    return WorldtubeTypeDType::ConstraintPreservingPhysicalSommerfeldGhostGauge;
  } else if (type_read == "ConstraintPreservingPhysicalOnlineGhostGauge") {
    return WorldtubeTypeDType::ConstraintPreservingPhysicalOnlineGhostGauge;
  } else if (type_read == "ConstraintPreservingPhysicalGhostGauge") {
    return WorldtubeTypeDType::ConstraintPreservingPhysicalGhostGauge;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to "
              "WorldtubeTypeDType::Type. Must "
              "be one of ConstraintPreserving, ConstraintPreservingPhysical, "
              "ConstraintPreservingPhysicalFrozenGauge, "
              "ConstraintPreservingPhysicalAnalyticGhostGauge, "
              "ConstraintPreservingPhysicalSommerfeldGhostGauge, "
              "ConstraintPreservingPhysicalOnlineGhostGauge or "
              "ConstraintPreservingPhysicalGhostGauge");
}
}  // namespace detail

template <size_t Dim>
WorldtubeTypeD<Dim>::WorldtubeTypeD(
    const detail::WorldtubeTypeDType type,
    std::optional<std::unique_ptr<evolution::initial_data::InitialData>>
        analytic_gauge_prescription,
    const double gauge_relaxation_rate, const Options::Context& context)
    : type_(type),
      analytic_gauge_prescription_(analytic_gauge_prescription.has_value()
                                       ? std::move(*analytic_gauge_prescription)
                                       : nullptr),
      gauge_relaxation_rate_(gauge_relaxation_rate) {
  if ((type_ == detail::WorldtubeTypeDType::
                    ConstraintPreservingPhysicalAnalyticGhostGauge or
       type_ == detail::WorldtubeTypeDType::
                    ConstraintPreservingPhysicalSommerfeldGhostGauge) and
      analytic_gauge_prescription_ == nullptr) {
    PARSE_ERROR(context,
                "This Type requires an AnalyticGaugePrescription, but None "
                "was given.");
  }
  if (gauge_relaxation_rate_ < 0.) {
    PARSE_ERROR(context,
                "GaugeRelaxationRate must be non-negative, but got "
                    << gauge_relaxation_rate_
                    << ". A negative rate drives the gauge sector away from "
                       "the model value.");
  }
}

template <size_t Dim>
WorldtubeTypeD<Dim>::WorldtubeTypeD(const WorldtubeTypeD& rhs)
    : BoundaryCondition<Dim>(rhs),
      type_(rhs.type_),
      analytic_gauge_prescription_(
          rhs.analytic_gauge_prescription_ == nullptr
              ? nullptr
              : rhs.analytic_gauge_prescription_->get_clone()),
      gauge_relaxation_rate_(rhs.gauge_relaxation_rate_) {}

template <size_t Dim>
WorldtubeTypeD<Dim>& WorldtubeTypeD<Dim>::operator=(const WorldtubeTypeD& rhs) {
  if (&rhs == this) {
    return *this;
  }
  BoundaryCondition<Dim>::operator=(rhs);
  type_ = rhs.type_;
  analytic_gauge_prescription_ =
      rhs.analytic_gauge_prescription_ == nullptr
          ? nullptr
          : rhs.analytic_gauge_prescription_->get_clone();
  gauge_relaxation_rate_ = rhs.gauge_relaxation_rate_;
  return *this;
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
  p | type_;
  p | analytic_gauge_prescription_;
  p | gauge_relaxation_rate_;
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
    // c.f. dg_gridless_tags
    const double time,
    const std::optional<gh::Worldtube::MatcherConfig>& matcher_config,
    const gh::Worldtube::MapParameterData& map_parameters) const {
  TempBuffer<tmpl::list<::Tags::TempI<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempiaa<1, Dim, Frame::Inertial, DataVector>,
                        ::Tags::TempII<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempii<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempa<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempa<1, Dim, Frame::Inertial, DataVector>,
                        ::Tags::TempA<1, Dim, Frame::Inertial, DataVector>,
                        ::Tags::TempA<2, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::TempAb<0, Dim, Frame::Inertial, DataVector>,
                        ::Tags::TempAA<1, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<1, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempiaa<2, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<2, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<3, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempa<2, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempa<3, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<4, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempiaa<3, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<5, Dim, Frame::Inertial, DataVector>,
                        ::Tags::Tempaa<6, Dim, Frame::Inertial, DataVector>,
                        gr::Tags::SpacetimeNormalOneForm<DataVector, Dim>,
                        // inertial time derivatives
                        ::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>,
                        ::Tags::dt<Tags::Pi<DataVector, Dim>>,
                        ::Tags::dt<Tags::Phi<DataVector, Dim>>>>
      local_buffer(get_size(get<0>(normal_covector)), 0.);
  get<0>(get<gr::Tags::SpacetimeNormalOneForm<DataVector, Dim>>(local_buffer)) =
      -get(lapse);
  const tnsr::a<DataVector, Dim, Frame::Inertial>&
      spacetime_unit_normal_one_form =
          get<gr::Tags::SpacetimeNormalOneForm<DataVector, Dim>>(local_buffer);
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_spacetime_metric;
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_pi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> dt_phi;
  if (face_mesh_velocity.has_value()) {
    for (size_t storage_index = 0; storage_index < dt_pi.size();
         ++storage_index) {
      dt_spacetime_metric[storage_index].set_data_ref(make_not_null(
          &get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              local_buffer)[storage_index]));
      dt_pi[storage_index].set_data_ref(
          make_not_null(&get<::Tags::dt<Tags::Pi<DataVector, Dim>>>(
              local_buffer)[storage_index]));
    }
    for (size_t storage_index = 0; storage_index < dt_phi.size();
         ++storage_index) {
      dt_phi[storage_index].set_data_ref(
          make_not_null(&get<::Tags::dt<Tags::Phi<DataVector, Dim>>>(
              local_buffer)[storage_index]));
    }
    // Compute inertial time derivative by subtracting mesh velocity from
    // logical time derivative.
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        dt_spacetime_metric.get(a, b) = logical_dt_spacetime_metric.get(a, b);
        dt_pi.get(a, b) = logical_dt_pi.get(a, b);
        for (size_t d = 0; d < Dim; ++d) {
          dt_spacetime_metric.get(a, b) -=
              face_mesh_velocity->get(d) * d_spacetime_metric.get(d, a, b);
          dt_pi.get(a, b) -= face_mesh_velocity->get(d) * d_pi.get(d, a, b);
        }
        for (size_t i = 0; i < Dim; ++i) {
          dt_phi.get(i, a, b) = logical_dt_phi.get(i, a, b);
          for (size_t d = 0; d < Dim; ++d) {
            dt_phi.get(i, a, b) -=
                face_mesh_velocity->get(d) * d_phi.get(d, i, a, b);
          }
        }
      }
    }
  } else {
    for (size_t storage_index = 0; storage_index < dt_pi.size();
         ++storage_index) {
      dt_spacetime_metric[storage_index].set_data_ref(
          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
          make_not_null(&const_cast<DataVector&>(
              logical_dt_spacetime_metric[storage_index])));
      dt_pi[storage_index].set_data_ref(make_not_null(
          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
          &const_cast<DataVector&>(logical_dt_pi[storage_index])));
    }
    for (size_t storage_index = 0; storage_index < dt_phi.size();
         ++storage_index) {
      dt_phi[storage_index].set_data_ref(make_not_null(
          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
          &const_cast<DataVector&>(logical_dt_phi[storage_index])));
    }
  }

  auto& unit_interface_normal_vector =
      get<::Tags::TempI<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& four_index_constraint =
      get<::Tags::Tempiaa<1, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& inverse_spatial_metric =
      get<::Tags::TempII<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& extrinsic_curvature =
      get<::Tags::Tempii<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& incoming_null_one_form =
      get<::Tags::Tempa<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& outgoing_null_one_form =
      get<::Tags::Tempa<1, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& incoming_null_vector =
      get<::Tags::TempA<1, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& outgoing_null_vector =
      get<::Tags::TempA<2, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& projection_ab =
      get<::Tags::Tempaa<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& projection_Ab =
      get<::Tags::TempAb<0, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& projection_AB =
      get<::Tags::TempAA<1, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& char_projected_rhs_dt_v_psi =
      get<::Tags::Tempaa<1, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& char_projected_rhs_dt_v_zero =
      get<::Tags::Tempiaa<2, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& char_projected_rhs_dt_v_plus =
      get<::Tags::Tempaa<2, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& char_projected_rhs_dt_v_minus =
      get<::Tags::Tempaa<3, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& constraint_char_zero_plus =
      get<::Tags::Tempa<2, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& constraint_char_zero_minus =
      get<::Tags::Tempa<3, Dim, Frame::Inertial, DataVector>>(local_buffer);

  typename Tags::CharacteristicSpeeds<DataVector, Dim>::type char_speeds;

  auto& bc_dt_v_psi =
      get<::Tags::Tempaa<4, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& bc_dt_v_zero =
      get<::Tags::Tempiaa<3, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& bc_dt_v_plus =
      get<::Tags::Tempaa<5, Dim, Frame::Inertial, DataVector>>(local_buffer);
  auto& bc_dt_v_minus =
      get<::Tags::Tempaa<6, Dim, Frame::Inertial, DataVector>>(local_buffer);

  compute_intermediate_vars(
      make_not_null(&unit_interface_normal_vector),
      make_not_null(&four_index_constraint),
      make_not_null(&inverse_spatial_metric),
      make_not_null(&extrinsic_curvature),
      make_not_null(&incoming_null_one_form),
      make_not_null(&outgoing_null_one_form),
      make_not_null(&incoming_null_vector),
      make_not_null(&outgoing_null_vector), make_not_null(&projection_ab),
      make_not_null(&projection_Ab), make_not_null(&projection_AB),
      make_not_null(&char_projected_rhs_dt_v_psi),
      make_not_null(&char_projected_rhs_dt_v_zero),
      make_not_null(&char_projected_rhs_dt_v_plus),
      make_not_null(&char_projected_rhs_dt_v_minus),
      make_not_null(&constraint_char_zero_plus),
      make_not_null(&constraint_char_zero_minus), make_not_null(&char_speeds),
      face_mesh_velocity, normal_covector, pi, phi, spacetime_metric, coords,
      gamma1, gamma2, lapse, shift, inverse_spacetime_metric,
      spacetime_unit_normal_vector, spacetime_unit_normal_one_form,
      three_index_constraint, gauge_source, spacetime_deriv_gauge_source, dt_pi,
      dt_phi, dt_spacetime_metric, d_pi, d_phi, d_spacetime_metric);

  // Account for moving mesh: char speeds -> cher speeds - n_i v^i_g
  if (face_mesh_velocity.has_value()) {
    const auto radial_mesh_velocity =
        get(dot_product(normal_covector, *face_mesh_velocity));
    for (size_t a = 0; a < 4; ++a) {
      char_speeds.at(a) -= radial_mesh_velocity;
    }
  }

  // If no point on the boundary has any incoming characteristic, return here
  if (min_characteristic_speed(char_speeds) >= 0.) {
    std::fill(dt_spacetime_metric_correction->begin(),
              dt_spacetime_metric_correction->end(), 0.);
    std::fill(dt_pi_correction->begin(), dt_pi_correction->end(), 0.);
    std::fill(dt_phi_correction->begin(), dt_phi_correction->end(), 0.);
    return {};
  }

  Bjorhus::constraint_preserving_corrections_dt_v_psi(
      make_not_null(&bc_dt_v_psi), unit_interface_normal_vector,
      three_index_constraint, char_speeds);

  Bjorhus::constraint_preserving_corrections_dt_v_zero(
      make_not_null(&bc_dt_v_zero), unit_interface_normal_vector,
      four_index_constraint, char_speeds);

  // In order to set dt<V+> = 0, the correction term returned here must be
  // b_correction = -1*existing(dt<V+>), such that
  // final(dt<V+>) = existing(dt<V+>) + b_correction = 0
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = a; b <= Dim; ++b) {
      bc_dt_v_plus.get(a, b) = -char_projected_rhs_dt_v_plus.get(a, b);
      // bc_dt_v_minus.get(a, b) = -char_projected_rhs_dt_v_minus.get(a, b);
      // bc_dt_v_psi.get(a, b) = -char_projected_rhs_dt_v_psi.get(a, b);
      for (size_t c = 0; c <= Dim; ++c) {
        // bc_dt_v_zero.get(a, b, c) = -char_projected_rhs_dt_v_zero.get(a, b,
        // c);
      }
    }
  }

  if (type_ == detail::WorldtubeTypeDType::ConstraintPreserving) {
    ERROR("wrong path");
    Bjorhus::constraint_preserving_gauge_corrections_dt_v_minus(
        make_not_null(&bc_dt_v_minus), gamma2, coords, incoming_null_one_form,
        outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
        projection_ab, projection_Ab, projection_AB,
        char_projected_rhs_dt_v_psi, char_projected_rhs_dt_v_minus,
        constraint_char_zero_plus, constraint_char_zero_minus, char_speeds);
  } else if (type_ ==
                 detail::WorldtubeTypeDType::ConstraintPreservingPhysical or
             type_ == detail::WorldtubeTypeDType::
                          ConstraintPreservingPhysicalFrozenGauge or
             type_ == detail::WorldtubeTypeDType::
                          ConstraintPreservingPhysicalAnalyticGhostGauge or
             type_ == detail::WorldtubeTypeDType::
                          ConstraintPreservingPhysicalSommerfeldGhostGauge or
             type_ == detail::WorldtubeTypeDType::
                          ConstraintPreservingPhysicalOnlineGhostGauge or
             type_ == detail::WorldtubeTypeDType::
                          ConstraintPreservingPhysicalGhostGauge) {
    // AnalyticGhostGauge leaves the gauge sector frozen here and adds only the
    // relaxation below; SommerfeldGhostGauge keeps the Sommerfeld radiation
    // term and adds the relaxation on top of it.
    const auto gauge_sector_condition =
        (type_ == detail::WorldtubeTypeDType::ConstraintPreservingPhysical or
         type_ == detail::WorldtubeTypeDType::
                      ConstraintPreservingPhysicalSommerfeldGhostGauge)
            ? Bjorhus::GaugeSectorCondition::Sommerfeld
            : Bjorhus::GaugeSectorCondition::Frozen;
    Bjorhus::
        constraint_preserving_gauge_physical_corrections_dt_v_minus_worldtube(
            make_not_null(&bc_dt_v_minus), gauge_sector_condition, gamma2,
            coords, normal_covector, unit_interface_normal_vector,
            spacetime_unit_normal_vector, incoming_null_one_form,
            outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
            projection_ab, projection_Ab, projection_AB, inverse_spatial_metric,
            extrinsic_curvature, spacetime_metric, inverse_spacetime_metric,
            three_index_constraint, char_projected_rhs_dt_v_psi,
            char_projected_rhs_dt_v_minus, constraint_char_zero_plus,
            constraint_char_zero_minus, phi, d_phi, d_pi, char_speeds);
  } else {
    ERROR(
        "Failed to set dtVMinus. Input option must be one of "
        "ConstraintPreserving, ConstraintPreservingPhysical, "
        "ConstraintPreservingPhysicalFrozenGauge or "
        "ConstraintPreservingPhysicalAnalyticGhostGauge");
  }

  if (type_ ==
      detail::WorldtubeTypeDType::ConstraintPreservingPhysicalGhostGauge) {
    // The gauge sector is imposed weakly through the ghost data and the
    // upwind penalty (see dg_ghost), so it must receive no Bjorhus
    // correction at all: the corrections above start every sector frozen
    // (bc_dt_v_minus = -char_projected_rhs), so add back the gauge
    // projection of the projected volume RHS. The helper adds
    // -kappa * P_gauge(source); with kappa = -1 and source =
    // char_projected_rhs this is exactly + P_gauge(char_projected_rhs).
    const DataVector minus_one_kappa(get_size(get(gamma2)), -1.0);
    Bjorhus::detail::add_gauge_sector_terms_to_dt_v_minus(
        make_not_null(&bc_dt_v_minus), minus_one_kappa, incoming_null_one_form,
        outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
        projection_Ab, char_projected_rhs_dt_v_minus);
  }

  if (type_ == detail::WorldtubeTypeDType::
                   ConstraintPreservingPhysicalAnalyticGhostGauge or
      type_ == detail::WorldtubeTypeDType::
                   ConstraintPreservingPhysicalSommerfeldGhostGauge or
      type_ == detail::WorldtubeTypeDType::
                   ConstraintPreservingPhysicalOnlineGhostGauge) {
    // Relax the gauge sector towards the model value,
    //   dt u^-_ab|gauge = -kappa (u^-_ab - u^-_model,ab)|gauge.
    // The gauge sector was left frozen above, i.e. the correction there is
    // -dt u^-_ab|gauge, so the term added here is the whole gauge condition.
    // Only u^-_model is needed, never its time derivative: that is the point of
    // the ghost form, since a Bjorhus form would need a second time derivative
    // of the model.
    using evolved_vars_tags = typename System<Dim>::variables_tag::tags_list;
    using EvolvedVars = tuples::tagged_tuple_from_typelist<evolved_vars_tags>;
    std::optional<EvolvedVars> model{};
    if (type_ == detail::WorldtubeTypeDType::
                     ConstraintPreservingPhysicalOnlineGhostGauge) {
      // The closed loop: the model comes from the online matcher's latest
      // fit.  In the algebraic value mode this is the strict slow-time
      // first-order model, so p_(1) is held fixed between fits; extrapolating
      // epsilon*p_(1) with D_t p_(1) would add O(epsilon^2) content. The
      // zeroth-order center is nevertheless advanced with its O(epsilon)
      // fitted velocity. The ODE modes retain their experimental
      // rate-resummed behavior. While no fit exists yet the gauge sector
      // simply stays frozen.
      if constexpr (Dim == 3) {
        if (not matcher_config.has_value()) {
          ERROR(
              "ConstraintPreservingPhysicalOnlineGhostGauge requires the "
              "WorldtubeMatcher option to be active, but it is None.");
        }
        if (map_parameters.valid) {
          const bool first_order_value_mode =
              not matcher_config->rate_ode and
              not matcher_config->second_order_ode and
              not matcher_config->stepper_ode;
          std::array<double, gh::Worldtube::num_map_parameters> p =
              map_parameters.p;
          const double dt_extrapolate = time - map_parameters.last_fit_time;
          if (not first_order_value_mode) {
            for (size_t a = 0; a < gh::Worldtube::num_map_parameters; ++a) {
              gsl::at(p, a) += dt_extrapolate * gsl::at(map_parameters.pdot, a);
            }
          }
          const std::array<double, 3> model_center =
              gh::Worldtube::detail::model_center(*matcher_config,
                                                  map_parameters, time);
          model.emplace();
          if (first_order_value_mode) {
            gh::Solutions::affine_map_model::
                first_order_boosted_evolved_variables(
                    make_not_null(
                        &get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(
                            *model)),
                    make_not_null(&get<Tags::Pi<DataVector, Dim>>(*model)),
                    make_not_null(&get<Tags::Phi<DataVector, Dim>>(*model)),
                    coords, matcher_config->mass, model_center, p,
                    map_parameters.bulk_velocity,
                    matcher_config->centre_advection);
          } else {
            gh::Solutions::affine_map_model::evolved_variables(
                make_not_null(
                    &get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(*model)),
                make_not_null(&get<Tags::Pi<DataVector, Dim>>(*model)),
                make_not_null(&get<Tags::Phi<DataVector, Dim>>(*model)), coords,
                matcher_config->mass, model_center, p, map_parameters.pdot,
                matcher_config->centre_advection);
          }
        }
      } else {
        ERROR(
            "ConstraintPreservingPhysicalOnlineGhostGauge is only "
            "implemented in 3 dimensions.");
      }
    } else {
      // all_solutions plus the affine-map model; keep in sync with the
      // executable's initial_data_list or the dispatch ERRORs at runtime
      using gauge_prescriptions = tmpl::conditional_t<
          Dim == 3,
          tmpl::push_back<gh::Solutions::all_solutions<Dim>,
                          gh::Solutions::AffineMappedHarmonicSchwarzschild>,
          gh::Solutions::all_solutions<Dim>>;
      model = call_with_dynamic_type<EvolvedVars, gauge_prescriptions>(
          analytic_gauge_prescription_.get(),
          [&coords, &time](const auto* const solution_or_data) {
            if constexpr (is_analytic_solution_v<
                              std::decay_t<decltype(*solution_or_data)>>) {
              return solution_or_data->variables(coords, time,
                                                 evolved_vars_tags{});
            } else {
              (void)time;
              return solution_or_data->variables(coords, evolved_vars_tags{});
            }
          });
    }

    if (model.has_value()) {
      const auto v_minus_numerical = get<Tags::VMinus<DataVector, Dim>>(
          characteristic_fields(gamma2, inverse_spatial_metric,
                                spacetime_metric, pi, phi, normal_covector));
      const auto v_minus_model =
          get<Tags::VMinus<DataVector, Dim>>(characteristic_fields(
              gamma2, inverse_spatial_metric,
              get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(*model),
              get<Tags::Pi<DataVector, Dim>>(*model),
              get<Tags::Phi<DataVector, Dim>>(*model), normal_covector));

      auto delta_v_minus = v_minus_numerical;
      for (size_t a = 0; a <= Dim; ++a) {
        for (size_t b = a; b <= Dim; ++b) {
          delta_v_minus.get(a, b) -= v_minus_model.get(a, b);
        }
      }
      const DataVector kappa(get_size(get(gamma2)), gauge_relaxation_rate_);
      Bjorhus::detail::add_gauge_sector_terms_to_dt_v_minus(
          make_not_null(&bc_dt_v_minus), kappa, incoming_null_one_form,
          outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
          projection_Ab, delta_v_minus);
    }
  }

  // Only add corrections at grid points where the char speeds are negative
  set_bc_corr_zero_when_char_speed_is_positive(make_not_null(&bc_dt_v_psi),
                                               char_speeds[0]);
  set_bc_corr_zero_when_char_speed_is_positive(make_not_null(&bc_dt_v_zero),
                                               char_speeds[1]);
  set_bc_corr_zero_when_char_speed_is_positive(make_not_null(&bc_dt_v_plus),
                                               char_speeds[2]);
  set_bc_corr_zero_when_char_speed_is_positive(make_not_null(&bc_dt_v_minus),
                                               char_speeds[3]);

  // The boundary conditions here are imposed as corrections to the projections
  // of the right-hand-sides of the GH evolution equations (i.e. using Bjorhus'
  // method), and are written down in Eq. (63) - (65) of Lindblom et al (2005).
  // Now that we have calculated those corrections, we project them back as
  // corrections to dt<evolved variables>
  auto dt_evolved_vars = evolved_fields_from_characteristic_fields(
      gamma2, bc_dt_v_psi, bc_dt_v_zero, bc_dt_v_plus, bc_dt_v_minus,
      normal_covector);

  *dt_pi_correction = get<Tags::Pi<DataVector, Dim>>(dt_evolved_vars);
  *dt_phi_correction = get<Tags::Phi<DataVector, Dim>>(dt_evolved_vars);
  *dt_spacetime_metric_correction =
      get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(dt_evolved_vars);

  // Note: the outer-boundary class vetoes any radial mesh velocity here, on the
  // grounds that an outward-moving boundary means an expanding domain whose
  // characteristic analysis is unclear. That veto is inappropriate at a
  // worldtube, and was inherited rather than re-derived.
  //
  // At an inner boundary `normal_covector` points into the excision, so an
  // excision that tracks a moving hole has normal . v_mesh > 0 over half the
  // sphere by construction — exactly the situation the worldtube scheme needs
  // for a binary, where the worldtube follows the smaller hole. The veto
  // therefore rejects the intended use.
  //
  // It is also redundant. The mesh velocity is already accounted for in the two
  // places where it matters: the advective terms subtracted from dt of the
  // evolved variables above, and the characteristic speeds, which have
  // normal . v_mesh removed before any correction is applied or zeroed (compare
  // DemandOutgoingCharSpeeds, which tests those corrected speeds rather than
  // the mesh velocity itself). So no check is needed here.

  return {};
}

template <size_t Dim>
std::optional<std::string> WorldtubeTypeD<Dim>::dg_ghost(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        spacetime_metric_ghost,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*> pi_ghost,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*> phi_ghost,
    const gsl::not_null<Scalar<DataVector>*> gamma1_ghost,
    const gsl::not_null<Scalar<DataVector>*> gamma2_ghost,
    const gsl::not_null<Scalar<DataVector>*> lapse_ghost,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> shift_ghost,
    const gsl::not_null<tnsr::II<DataVector, Dim, Frame::Inertial>*>
        inv_spatial_metric_ghost,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*face_mesh_velocity*/,
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
    const tnsr::iaa<DataVector, Dim,
                    Frame::Inertial>& /*three_index_constraint*/,
    const tnsr::a<DataVector, Dim, Frame::Inertial>& /*gauge_source*/,
    const tnsr::ab<DataVector, Dim, Frame::Inertial>&
    /*spacetime_deriv_gauge_source*/,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>&
    /*logical_dt_spacetime_metric*/,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& /*logical_dt_pi*/,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& /*logical_dt_phi*/,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& /*d_spacetime_metric*/,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& /*d_pi*/,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& /*d_phi*/,
    const double time,
    const std::optional<gh::Worldtube::MatcherConfig>& matcher_config,
    const gh::Worldtube::MapParameterData& map_parameters) const {
  // Interior copies first: for every Type except
  // ConstraintPreservingPhysicalGhostGauge the ghost state equals the
  // interior, the upwind flux is consistent, and the penalty contributes
  // exactly zero.
  *spacetime_metric_ghost = spacetime_metric;
  *pi_ghost = pi;
  *phi_ghost = phi;
  *gamma1_ghost = gamma1;
  *gamma2_ghost = gamma2;
  *lapse_ghost = lapse;
  *shift_ghost = shift;
  const DataVector one_over_lapse_sqrd = 1.0 / (get(lapse) * get(lapse));
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      inv_spatial_metric_ghost->get(i, j) =
          inverse_spacetime_metric.get(1 + i, 1 + j) +
          shift.get(i) * shift.get(j) * one_over_lapse_sqrd;
    }
  }

  if (type_ !=
      detail::WorldtubeTypeDType::ConstraintPreservingPhysicalGhostGauge) {
    return {};
  }

  if constexpr (Dim == 3) {
    if (not matcher_config.has_value()) {
      ERROR(
          "ConstraintPreservingPhysicalGhostGauge requires the "
          "WorldtubeMatcher option to be active, but it is None.");
    }
    if (not map_parameters.valid) {
      // no fit yet: ghost = interior, no driving
      return {};
    }

    // Model fields from the online matcher. The algebraic mode is the strict
    // slow-time first-order system: hold p_(1) fixed within the fit interval,
    // omit D_t p_(1), and kinematically advance the zeroth-order center. The
    // derivative ODE modes retain the old resummed extrapolation as
    // higher-order experiments.
    const bool first_order_value_mode = not matcher_config->rate_ode and
                                        not matcher_config->second_order_ode and
                                        not matcher_config->stepper_ode;
    std::array<double, gh::Worldtube::num_map_parameters> p = map_parameters.p;
    const double dt_extrapolate = time - map_parameters.last_fit_time;
    if (not first_order_value_mode) {
      for (size_t a = 0; a < gh::Worldtube::num_map_parameters; ++a) {
        gsl::at(p, a) += dt_extrapolate * gsl::at(map_parameters.pdot, a);
      }
    }
    const std::array<double, 3> model_center =
        gh::Worldtube::detail::model_center(*matcher_config, map_parameters,
                                            time);
    tnsr::aa<DataVector, Dim, Frame::Inertial> model_metric{};
    tnsr::aa<DataVector, Dim, Frame::Inertial> model_pi{};
    tnsr::iaa<DataVector, Dim, Frame::Inertial> model_phi{};
    if (first_order_value_mode) {
      gh::Solutions::affine_map_model::first_order_boosted_evolved_variables(
          make_not_null(&model_metric), make_not_null(&model_pi),
          make_not_null(&model_phi), coords, matcher_config->mass, model_center,
          p, map_parameters.bulk_velocity, matcher_config->centre_advection);
    } else {
      gh::Solutions::affine_map_model::evolved_variables(
          make_not_null(&model_metric), make_not_null(&model_pi),
          make_not_null(&model_phi), coords, matcher_config->mass, model_center,
          p, map_parameters.pdot, matcher_config->centre_advection);
    }

    // characteristic fields of the interior and of the model, in the same
    // (interior) frame
    const auto char_fields_interior =
        characteristic_fields(gamma2, *inv_spatial_metric_ghost,
                              spacetime_metric, pi, phi, normal_covector);
    const auto v_minus_model = get<Tags::VMinus<DataVector, Dim>>(
        characteristic_fields(gamma2, *inv_spatial_metric_ghost, model_metric,
                              model_pi, model_phi, normal_covector));
    const auto& v_minus_interior =
        get<Tags::VMinus<DataVector, Dim>>(char_fields_interior);

    // null frame pieces for the gauge projector
    const size_t n_points = get(lapse).size();
    tnsr::a<DataVector, Dim, Frame::Inertial> normal_one_form(n_points, 0.);
    get<0>(normal_one_form) = -get(lapse);
    tnsr::I<DataVector, Dim, Frame::Inertial> unit_interface_normal_vector(
        n_points);
    raise_or_lower_index(make_not_null(&unit_interface_normal_vector),
                         normal_covector, *inv_spatial_metric_ghost);
    tnsr::a<DataVector, Dim, Frame::Inertial> incoming_null_one_form(n_points);
    tnsr::a<DataVector, Dim, Frame::Inertial> outgoing_null_one_form(n_points);
    tnsr::A<DataVector, Dim, Frame::Inertial> incoming_null_vector(n_points);
    tnsr::A<DataVector, Dim, Frame::Inertial> outgoing_null_vector(n_points);
    gr::interface_null_normal(make_not_null(&incoming_null_one_form),
                              normal_one_form, normal_covector, shift, -1.);
    gr::interface_null_normal(make_not_null(&outgoing_null_one_form),
                              normal_one_form, normal_covector, shift, 1.);
    gr::interface_null_normal(make_not_null(&incoming_null_vector),
                              spacetime_unit_normal_vector,
                              unit_interface_normal_vector, -1.);
    gr::interface_null_normal(make_not_null(&outgoing_null_vector),
                              spacetime_unit_normal_vector,
                              unit_interface_normal_vector, 1.);
    tnsr::Ab<DataVector, Dim, Frame::Inertial> projection_Ab(n_points);
    gr::transverse_projection_operator(
        make_not_null(&projection_Ab), spacetime_unit_normal_vector,
        normal_one_form, unit_interface_normal_vector, normal_covector, shift);

    // v^-_ghost = v^-_interior + P_gauge(v^-_model - v^-_interior), using
    // the gauge-sector helper: it adds -kappa * P_gauge(delta), so with
    // kappa = 1 and delta = (v^-_interior - v^-_model) it adds exactly the
    // required replacement term.
    auto v_minus_ghost = v_minus_interior;
    auto delta_v_minus = v_minus_interior;
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = a; b <= Dim; ++b) {
        delta_v_minus.get(a, b) -= v_minus_model.get(a, b);
      }
    }
    const DataVector unit_kappa(n_points, 1.0);
    Bjorhus::detail::add_gauge_sector_terms_to_dt_v_minus(
        make_not_null(&v_minus_ghost), unit_kappa, incoming_null_one_form,
        outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
        projection_Ab, delta_v_minus);

    // reassemble the ghost evolved fields; only Pi and the normal part of
    // Phi change (the metric is the v_psi characteristic)
    const auto ghost_evolved = evolved_fields_from_characteristic_fields(
        gamma2,
        get<Tags::VSpacetimeMetric<DataVector, Dim>>(char_fields_interior),
        get<Tags::VZero<DataVector, Dim>>(char_fields_interior),
        get<Tags::VPlus<DataVector, Dim>>(char_fields_interior), v_minus_ghost,
        normal_covector);
    *spacetime_metric_ghost =
        get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(ghost_evolved);
    *pi_ghost = get<Tags::Pi<DataVector, Dim>>(ghost_evolved);
    *phi_ghost = get<Tags::Phi<DataVector, Dim>>(ghost_evolved);
    return {};
  } else {
    ERROR(
        "ConstraintPreservingPhysicalGhostGauge is only implemented in 3 "
        "dimensions.");
  }
}

template <size_t Dim>
void WorldtubeTypeD<Dim>::compute_intermediate_vars(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        unit_interface_normal_vector,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        four_index_constraint,
    const gsl::not_null<tnsr::II<DataVector, Dim, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<tnsr::ii<DataVector, Dim, Frame::Inertial>*>
        extrinsic_curvature,
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
        incoming_null_one_form,
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
        outgoing_null_one_form,
    const gsl::not_null<tnsr::A<DataVector, Dim, Frame::Inertial>*>
        incoming_null_vector,
    const gsl::not_null<tnsr::A<DataVector, Dim, Frame::Inertial>*>
        outgoing_null_vector,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        projection_ab,
    const gsl::not_null<tnsr::Ab<DataVector, Dim, Frame::Inertial>*>
        projection_Ab,
    const gsl::not_null<tnsr::AA<DataVector, Dim, Frame::Inertial>*>
        projection_AB,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        char_projected_rhs_dt_v_psi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        char_projected_rhs_dt_v_zero,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        char_projected_rhs_dt_v_plus,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        char_projected_rhs_dt_v_minus,
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
        constraint_char_zero_plus,
    const gsl::not_null<tnsr::a<DataVector, Dim, Frame::Inertial>*>
        constraint_char_zero_minus,
    const gsl::not_null<std::array<DataVector, 4>*> char_speeds,

    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& /*coords*/,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& shift,
    const tnsr::AA<DataVector, Dim, Frame::Inertial>& inverse_spacetime_metric,
    const tnsr::A<DataVector, Dim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::a<DataVector, Dim, Frame::Inertial>&
        spacetime_unit_normal_one_form,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& three_index_constraint,
    const tnsr::a<DataVector, Dim, Frame::Inertial>& gauge_source,
    const tnsr::ab<DataVector, Dim, Frame::Inertial>&
        spacetime_deriv_gauge_source,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& dt_pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& dt_phi,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& dt_spacetime_metric,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& /* d_spacetime_metric */)
    const {
  TempBuffer<tmpl::list<::Tags::TempScalar<0, DataVector>,
                        ::Tags::Tempia<0, Dim, Frame::Inertial, DataVector>>>
      local_buffer(get_size(get<0>(normal_covector)), 0.);
  auto& one_over_lapse_sqrd =
      get(get<::Tags::TempScalar<0, DataVector>>(local_buffer));
  auto& two_index_constraint =
      get<::Tags::Tempia<0, Dim, Frame::Inertial, DataVector>>(local_buffer);

  one_over_lapse_sqrd = 1.0 / (get(lapse) * get(lapse));
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      inverse_spatial_metric->get(i, j) =
          inverse_spacetime_metric.get(1 + i, 1 + j) +
          (shift.get(i) * shift.get(j) * one_over_lapse_sqrd);
    }
  }

  raise_or_lower_index(unit_interface_normal_vector, normal_covector,
                       *inverse_spatial_metric);
  gh::extrinsic_curvature(extrinsic_curvature, spacetime_unit_normal_vector, pi,
                          phi);

  if (LIKELY(Dim == 3)) {
    gh::four_index_constraint(four_index_constraint, d_phi);
  } else if (UNLIKELY(Dim == 2)) {
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = 0; b <= Dim; ++b) {
        four_index_constraint->get(0, a, b) =
            d_phi.get(0, 1, a, b) - d_phi.get(1, 0, a, b);
        four_index_constraint->get(1, a, b) =
            -four_index_constraint->get(0, a, b);
      }
    }
  } else {
    std::fill(four_index_constraint->begin(), four_index_constraint->end(), 0.);
  }

  gr::interface_null_normal(incoming_null_one_form,
                            spacetime_unit_normal_one_form, normal_covector,
                            shift, -1.);
  gr::interface_null_normal(outgoing_null_one_form,
                            spacetime_unit_normal_one_form, normal_covector,
                            shift, 1.);
  gr::interface_null_normal(incoming_null_vector, spacetime_unit_normal_vector,
                            *unit_interface_normal_vector, -1.);
  gr::interface_null_normal(outgoing_null_vector, spacetime_unit_normal_vector,
                            *unit_interface_normal_vector, 1.);

  gr::transverse_projection_operator(projection_ab, spacetime_metric,
                                     spacetime_unit_normal_one_form,
                                     normal_covector, shift);
  gr::transverse_projection_operator(
      projection_Ab, spacetime_unit_normal_vector,
      spacetime_unit_normal_one_form, *unit_interface_normal_vector,
      normal_covector, shift);
  gr::transverse_projection_operator(projection_AB, inverse_spacetime_metric,
                                     spacetime_unit_normal_vector,
                                     *unit_interface_normal_vector);

  const auto dt_char_fields = characteristic_fields(
      gamma2, *inverse_spatial_metric, dt_spacetime_metric, dt_pi, dt_phi,
      normal_covector);
  *char_projected_rhs_dt_v_psi =
      get<Tags::VSpacetimeMetric<DataVector, Dim>>(dt_char_fields);
  *char_projected_rhs_dt_v_zero =
      get<Tags::VZero<DataVector, Dim>>(dt_char_fields);
  *char_projected_rhs_dt_v_plus =
      get<Tags::VPlus<DataVector, Dim>>(dt_char_fields);
  *char_projected_rhs_dt_v_minus =
      get<Tags::VMinus<DataVector, Dim>>(dt_char_fields);

  // c^{\hat{0}-}_a = F_a + n^k C_{ka}
  gh::two_index_constraint(
      make_not_null(&two_index_constraint), spacetime_deriv_gauge_source,
      spacetime_unit_normal_one_form, spacetime_unit_normal_vector,
      *inverse_spatial_metric, inverse_spacetime_metric, pi, phi, d_pi, d_phi,
      gamma2, three_index_constraint);
  f_constraint(constraint_char_zero_plus, gauge_source,
               spacetime_deriv_gauge_source, spacetime_unit_normal_one_form,
               spacetime_unit_normal_vector, *inverse_spatial_metric,
               inverse_spacetime_metric, pi, phi, d_pi, d_phi, gamma2,
               three_index_constraint);
  for (size_t a = 0; a < Dim + 1; ++a) {
    constraint_char_zero_minus->get(a) = constraint_char_zero_plus->get(a);
    for (size_t i = 0; i < Dim; ++i) {
      constraint_char_zero_plus->get(a) -=
          unit_interface_normal_vector->get(i) * two_index_constraint.get(i, a);
      constraint_char_zero_minus->get(a) +=
          unit_interface_normal_vector->get(i) * two_index_constraint.get(i, a);
    }
  }
  characteristic_speeds(char_speeds, gamma1, lapse, shift, normal_covector,
                        face_mesh_velocity);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID WorldtubeTypeD<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) template class WorldtubeTypeD<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace gh::BoundaryConditions
