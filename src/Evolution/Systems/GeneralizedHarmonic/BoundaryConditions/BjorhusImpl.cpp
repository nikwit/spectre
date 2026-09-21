// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Characteristics.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/CovariantDerivOfExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/InterfaceNullNormal.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylPropagating.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gh::BoundaryConditions::Bjorhus {
template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_psi(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_psi,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const std::array<DataType, 4>& char_speeds) {
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      bc_dt_v_psi->get(a, b) = char_speeds[0] *
                               unit_interface_normal_vector.get(0) *
                               three_index_constraint.get(0, a, b);
      for (size_t i = 1; i < VolumeDim; ++i) {
        bc_dt_v_psi->get(a, b) += char_speeds[0] *
                                  unit_interface_normal_vector.get(i) *
                                  three_index_constraint.get(i, a, b);
      }
    }
  }
}

template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_zero(
    const gsl::not_null<tnsr::iaa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_zero,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        four_index_constraint,
    const std::array<DataType, 4>& char_speeds) {
  set_number_of_grid_points(bc_dt_v_zero, unit_interface_normal_vector);
  std::fill(bc_dt_v_zero->begin(), bc_dt_v_zero->end(), 0.);

  if (VolumeDim == 3) {
    for (size_t a = 0; a <= VolumeDim; ++a) {
      for (size_t b = a; b <= VolumeDim; ++b) {
        // Lets say this term is T2_{iab} := - n_l \beta^l n^j C_{jiab}.
        // But we store D_{iab} = LeviCivita^{ijk} dphi_{jkab},
        // and C_{ijab} = LeviCivita^{kij} D_{kab}
        // where D is `four_index_constraint`.
        // therefore, T2_{iab} =  char_speed<VZero> n^j C_{jiab}
        // (since char_speed<VZero> = - n_l \beta^l), and therefore:
        // T2_{iab} = char_speed<VZero> n^j LeviCivita^{ikj} D_{kab}.
        // Let LeviCivitaIterator be indexed by
        // it[0] <--> i,
        // it[1] <--> j,
        // it[2] <--> k, then
        // T2_{it[0], ab} += char_speed<VZero> n^it[2] it.sign() D_{it[1], ab};
        for (LeviCivitaIterator<VolumeDim> it; it; ++it) {
          bc_dt_v_zero->get(it[0], a, b) +=
              it.sign() * char_speeds[1] *
              unit_interface_normal_vector.get(it[2]) *
              four_index_constraint.get(it[1], a, b);
        }
      }
    }
  } else if (VolumeDim == 2) {
    for (size_t a = 0; a <= VolumeDim; ++a) {
      for (size_t b = a; b <= VolumeDim; ++b) {
        // Lets say this term is T2_{kab} := - n_l \beta^l n^j C_{jkab}.
        // In 2+1 spacetime, we store the four index constraint to
        // be D_{1ab} = C_{12ab}, C_{2ab} = C_{21ab}. Therefore,
        // T_{kab} = -n_l \beta^l (n^1 C_{1kab} + n^2 C_{2kab}), i.e.
        // T_{1ab} = -n_l \beta^l n^2 D_{2ab}, T_{2ab} = -n_l \beta^l n^1
        // D_{1ab}.
        bc_dt_v_zero->get(0, a, b) +=
            char_speeds[1] * (unit_interface_normal_vector.get(1) *
                              four_index_constraint.get(1, a, b));
        bc_dt_v_zero->get(1, a, b) +=
            char_speeds[1] * (unit_interface_normal_vector.get(0) *
                              four_index_constraint.get(0, a, b));
      }
    }
  }
}

namespace detail {
template <size_t VolumeDim, typename DataType>
void add_constraint_sector_projection(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source) {
  // The tensor structure multiplying dt v^-_{cd} in T^C_{ab} of Eq. (64) of
  // Lindblom et al (2005).
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      for (size_t c = 0; c <= VolumeDim; ++c) {
        for (size_t d = 0; d <= VolumeDim; ++d) {
          result->get(a, b) +=
              0.5 * scalar_coefficient *
              (2. * incoming_null_vector.get(c) * incoming_null_vector.get(d) *
                   outgoing_null_one_form.get(a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(c) * projection_Ab.get(d, a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(c) * projection_Ab.get(d, b) *
                   outgoing_null_one_form.get(a) -
               incoming_null_vector.get(d) * projection_Ab.get(c, a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(d) * projection_Ab.get(c, b) *
                   outgoing_null_one_form.get(a) +
               projection_AB.get(c, d) * projection_ab.get(a, b)) *
              source.get(c, d);
        }
      }
    }
  }
}

template <size_t VolumeDim, typename DataType>
void add_physical_sector_projection(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source) {
  // The transverse-traceless projection with respect to the two-surface.
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      for (size_t c = 0; c <= VolumeDim; ++c) {
        for (size_t d = 0; d <= VolumeDim; ++d) {
          result->get(a, b) +=
              scalar_coefficient *
              (projection_Ab.get(c, a) * projection_Ab.get(d, b) -
               0.5 * projection_ab.get(a, b) * projection_AB.get(c, d)) *
              source.get(c, d);
        }
      }
    }
  }
}

template <size_t VolumeDim, typename DataType>
void add_gauge_sector_projection(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*> result,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source) {
  // Minus the tensor structure multiplying dt v^g_{cd} in T^G_{ab} of Eq. (64)
  // of Lindblom et al (2005): with this sign the three sector projections
  // partition the identity.
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      for (size_t c = 0; c <= VolumeDim; ++c) {
        for (size_t d = 0; d <= VolumeDim; ++d) {
          result->get(a, b) -=
              scalar_coefficient *
              (incoming_null_one_form.get(a) * projection_Ab.get(c, b) *
                   outgoing_null_vector.get(d) +
               incoming_null_one_form.get(b) * projection_Ab.get(c, a) *
                   outgoing_null_vector.get(d) -
               (incoming_null_one_form.get(a) * outgoing_null_one_form.get(b) *
                    incoming_null_vector.get(c) * outgoing_null_vector.get(d) +
                incoming_null_one_form.get(b) * outgoing_null_one_form.get(a) *
                    incoming_null_vector.get(c) * outgoing_null_vector.get(d) +
                incoming_null_one_form.get(a) * incoming_null_one_form.get(b) *
                    outgoing_null_vector.get(c) *
                    outgoing_null_vector.get(d))) *
              source.get(c, d);
        }
      }
    }
  }
}

template <size_t VolumeDim, typename DataType>
void add_gauge_sector_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const DataType& scalar_coefficient,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& source) {
  const DataType minus_coefficient = -scalar_coefficient;
  add_gauge_sector_projection(bc_dt_v_minus, minus_coefficient,
                              incoming_null_one_form, outgoing_null_one_form,
                              incoming_null_vector, outgoing_null_vector,
                              projection_Ab, source);
}

template <size_t VolumeDim, typename DataType>
void add_gauge_sommerfeld_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi) {
  // gauge_bc_coeff below is hard-coded here to its default value in SpEC
  constexpr double gauge_bc_coeff = 1.;

  DataType inertial_radius_or_scalar_factor(get_size(get<0>(inertial_coords)),
                                            0.);
  for (size_t i = 0; i < VolumeDim; ++i) {
    inertial_radius_or_scalar_factor += square(inertial_coords.get(i));
  }
  inertial_radius_or_scalar_factor =
      get(gamma2) - (gauge_bc_coeff / sqrt(inertial_radius_or_scalar_factor));
  add_gauge_sector_terms_to_dt_v_minus(
      bc_dt_v_minus, inertial_radius_or_scalar_factor, incoming_null_one_form,
      outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
      projection_Ab, char_projected_rhs_dt_v_psi);
}

template <size_t VolumeDim, typename DataType>
void add_constraint_dependent_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const std::array<DataType, 4>& char_speeds) {
  constexpr double mu = 0.;  // hard-coded value from SpEC Bbh input file Mu = 0
  const double one_by_sqrt_2 = 1. / sqrt(2.);

  // Add corrections c.f. Eq (64) of gr-qc/0512093
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      for (size_t c = 0; c <= VolumeDim; ++c) {
        for (size_t d = 0; d <= VolumeDim; ++d) {
          bc_dt_v_minus->get(a, b) +=
              0.5 *
              (2. * incoming_null_vector.get(c) * incoming_null_vector.get(d) *
                   outgoing_null_one_form.get(a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(c) * projection_Ab.get(d, a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(c) * projection_Ab.get(d, b) *
                   outgoing_null_one_form.get(a) -
               incoming_null_vector.get(d) * projection_Ab.get(c, a) *
                   outgoing_null_one_form.get(b) -
               incoming_null_vector.get(d) * projection_Ab.get(c, b) *
                   outgoing_null_one_form.get(a) +
               projection_AB.get(c, d) * projection_ab.get(a, b)) *
              char_projected_rhs_dt_v_minus.get(c, d);
        }
      }
      if constexpr (mu == 0.) {
        for (size_t c = 0; c <= VolumeDim; ++c) {
          bc_dt_v_minus->get(a, b) +=
              one_by_sqrt_2 * char_speeds[3] *
              constraint_char_zero_minus.get(c) *
              (outgoing_null_one_form.get(a) * outgoing_null_one_form.get(b) *
                   incoming_null_vector.get(c) +
               projection_ab.get(a, b) * outgoing_null_vector.get(c) -
               projection_Ab.get(c, b) * outgoing_null_one_form.get(a) -
               projection_Ab.get(c, a) * outgoing_null_one_form.get(b));
        }
      } else {
        for (size_t c = 0; c <= VolumeDim; ++c) {
          bc_dt_v_minus->get(a, b) +=
              one_by_sqrt_2 * char_speeds[3] *
              (constraint_char_zero_minus.get(c) -
               mu * constraint_char_zero_plus.get(c)) *
              (outgoing_null_one_form.get(a) * outgoing_null_one_form.get(b) *
                   incoming_null_vector.get(c) +
               projection_ab.get(a, b) * outgoing_null_vector.get(c) -
               projection_Ab.get(c, b) * outgoing_null_one_form.get(a) -
               projection_Ab.get(c, a) * outgoing_null_one_form.get(b));
        }
      }
    }
  }
}

template <size_t VolumeDim, typename DataType>
void add_physical_terms_to_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::i<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::II<DataType, VolumeDim, Frame::Inertial>&
        inverse_spatial_metric,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>& extrinsic_curvature,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& spacetime_metric,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>&
        inverse_spacetime_metric,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& phi,
    const tnsr::ijaa<DataType, VolumeDim, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& d_pi,
    const std::array<DataType, 4>& char_speeds, const double time,
    const MathFunction<1, Frame::Inertial>* const incoming_wave_profile,
    const std::array<double, 6>& incoming_wave_components,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>* const
        incoming_weyl_propagating_minus) {
  // hard-coded value from SpEC Bbh input file Mu = MuPhys = 0
  constexpr double mu_phys = 0.;
  constexpr bool adjust_phys_using_c4 = true;
  constexpr bool gamma2_in_phys = true;

  // In what follows, we follow Kidder, Scheel & Teukolsky (2001)
  // https://arxiv.org/pdf/gr-qc/0105031.pdf.
  TempBuffer<tmpl::list<::Tags::Tempaa<0, VolumeDim, Frame::Inertial, DataType>,
                        ::Tags::Tempaa<1, VolumeDim, Frame::Inertial, DataType>,
                        ::Tags::Tempaa<2, VolumeDim, Frame::Inertial, DataType>,
                        ::Tags::TempScalar<0, DataType>>>
      u3_buffer(get_size(get<0>(unit_interface_normal_vector)), 0.);
  auto& U3p =
      get<::Tags::Tempaa<0, VolumeDim, Frame::Inertial, DataType>>(u3_buffer);
  auto& U3m =
      get<::Tags::Tempaa<1, VolumeDim, Frame::Inertial, DataType>>(u3_buffer);
  auto& injected_wave =
      get<::Tags::Tempaa<2, VolumeDim, Frame::Inertial, DataType>>(u3_buffer);

  {
    TempBuffer<
        tmpl::list<::Tags::Tempijj<0, VolumeDim, Frame::Inertial, DataType>,
                   // cov deriv of Kij
                   ::Tags::Tempijj<1, VolumeDim, Frame::Inertial, DataType>,
                   // spatial Ricci
                   ::Tags::Tempii<0, VolumeDim, Frame::Inertial, DataType>,
                   // spatial projection operators P_ij, P^ij, and P^i_j
                   ::Tags::TempII<0, VolumeDim, Frame::Inertial, DataType>,
                   ::Tags::Tempii<1, VolumeDim, Frame::Inertial, DataType>,
                   ::Tags::TempIj<0, VolumeDim, Frame::Inertial, DataType>,
                   // weyl propagating modes
                   ::Tags::Tempii<2, VolumeDim, Frame::Inertial, DataType>,
                   ::Tags::Tempii<3, VolumeDim, Frame::Inertial, DataType>,
                   ::Tags::Tempii<4, VolumeDim, Frame::Inertial, DataType>>>
        local_buffer(get_size(get<0>(unit_interface_normal_vector)), 0.);

    auto& spatial_phi =
        get<::Tags::Tempijj<0, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& cov_deriv_ex_curv =
        get<::Tags::Tempijj<1, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& ricci_3 =
        get<::Tags::Tempii<0, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& spatial_projection_IJ =
        get<::Tags::TempII<0, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& spatial_projection_ij =
        get<::Tags::Tempii<1, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& spatial_projection_Ij =
        get<::Tags::TempIj<0, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto weyl_prop_minus =
        get<::Tags::Tempii<2, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);
    auto& spatial_metric =
        get<::Tags::Tempii<3, VolumeDim, Frame::Inertial, DataType>>(
            local_buffer);

    // D_(k,i,j) = (1/2) \partial_k g_(ij) and its derivative
    for (size_t i = 0; i < VolumeDim; ++i) {
      for (size_t j = i; j < VolumeDim; ++j) {
        for (size_t k = 0; k < VolumeDim; ++k) {
          spatial_phi.get(k, i, j) = phi.get(k, i + 1, j + 1);
        }
      }
    }

    // Compute covariant deriv of extrinsic curvature
    gh::covariant_deriv_of_extrinsic_curvature(
        make_not_null(&cov_deriv_ex_curv), extrinsic_curvature,
        spacetime_unit_normal_vector,
        raise_or_lower_first_index(gr::christoffel_first_kind(spatial_phi),
                                   inverse_spatial_metric),
        inverse_spacetime_metric, phi, d_pi, d_phi);

    // Compute spatial Ricci tensor
    gh::spatial_ricci_tensor(make_not_null(&ricci_3), phi, d_phi,
                             inverse_spatial_metric);

    if (adjust_phys_using_c4) {
      // This adds 4-index constraint terms to 3Ricci so as to cancel
      // out normal derivatives from the final expression for U8.
      // It is much easier to add them here than to recalculate U8
      // from scratch.

      // Add some 4-index constraint terms to 3Ricci.
      for (size_t i = 0; i < VolumeDim; ++i) {
        for (size_t j = i; j < VolumeDim; ++j) {
          for (size_t k = 0; k < VolumeDim; ++k) {
            for (size_t l = 0; l < VolumeDim; ++l) {
              ricci_3.get(i, j) += 0.25 * inverse_spatial_metric.get(k, l) *
                                   (d_phi.get(i, k, 1 + l, 1 + j) -
                                    d_phi.get(k, i, 1 + l, 1 + j) +
                                    d_phi.get(j, k, 1 + l, 1 + i) -
                                    d_phi.get(k, j, 1 + l, 1 + i));
            }
          }
        }
      }

      // Add more 4-index constraint terms to 3Ricci
      // These compensate for some of the cov_deriv_ex_curv terms.
      for (size_t i = 0; i < VolumeDim; ++i) {
        for (size_t j = i; j < VolumeDim; ++j) {
          for (size_t a = 0; a <= VolumeDim; ++a) {
            for (size_t k = 0; k < VolumeDim; ++k) {
              ricci_3.get(i, j) +=
                  0.5 * unit_interface_normal_vector.get(k) *
                  spacetime_unit_normal_vector.get(a) *
                  (d_phi.get(i, k, j + 1, a) - d_phi.get(k, i, j + 1, a) +
                   d_phi.get(j, k, i + 1, a) - d_phi.get(k, j, i + 1, a));
            }
          }
        }
      }
    }

    // Make spatial projection operators
    for (size_t j = 0; j < VolumeDim; ++j) {
      for (size_t k = j; k < VolumeDim; ++k) {
        spatial_metric.get(j, k) = spacetime_metric.get(1 + j, 1 + k);
      }
    }
    gr::transverse_projection_operator(make_not_null(&spatial_projection_IJ),
                                       inverse_spatial_metric,
                                       unit_interface_normal_vector);
    gr::transverse_projection_operator(make_not_null(&spatial_projection_ij),
                                       spatial_metric,
                                       unit_interface_normal_one_form);
    gr::transverse_projection_operator(make_not_null(&spatial_projection_Ij),
                                       unit_interface_normal_vector,
                                       unit_interface_normal_one_form);

    // Weyl propagating mode
    gr::weyl_propagating(make_not_null(&weyl_prop_minus), ricci_3,
                         extrinsic_curvature, inverse_spatial_metric,
                         cov_deriv_ex_curv, unit_interface_normal_vector,
                         spatial_projection_IJ, spatial_projection_ij,
                         spatial_projection_Ij, -1);
    if (incoming_weyl_propagating_minus != nullptr) {
      // Drive the incoming mode towards the supplied value instead of zero.
      for (size_t i = 0; i < VolumeDim; ++i) {
        for (size_t j = i; j < VolumeDim; ++j) {
          weyl_prop_minus.get(i, j) -=
              incoming_weyl_propagating_minus->get(i, j);
        }
      }
    }

    if constexpr (mu_phys == 0.) {
      // No need to compute U3p or weyl_prop_plus in this case
      for (size_t a = 0; a <= VolumeDim; ++a) {
        for (size_t b = a; b <= VolumeDim; ++b) {
          for (size_t i = 0; i < VolumeDim; ++i) {
            for (size_t j = 0; j < VolumeDim; ++j) {
              U3m.get(a, b) += 2. * projection_Ab.get(i + 1, a) *
                               projection_Ab.get(j + 1, b) *
                               weyl_prop_minus.get(i, j);
            }
          }
        }
      }
    } else {
      auto& weyl_prop_plus =
          get<::Tags::Tempii<3, VolumeDim, Frame::Inertial, DataType>>(
              local_buffer);
      gr::weyl_propagating(make_not_null(&weyl_prop_plus), ricci_3,
                           extrinsic_curvature, inverse_spatial_metric,
                           cov_deriv_ex_curv, unit_interface_normal_vector,
                           spatial_projection_IJ, spatial_projection_ij,
                           spatial_projection_Ij, 1);

      for (size_t a = 0; a <= VolumeDim; ++a) {
        for (size_t b = a; b <= VolumeDim; ++b) {
          for (size_t i = 0; i < VolumeDim; ++i) {
            for (size_t j = 0; j < VolumeDim; ++j) {
              U3p.get(a, b) += 2. * projection_Ab.get(i + 1, a) *
                               projection_Ab.get(j + 1, b) *
                               weyl_prop_plus.get(i, j);
              U3m.get(a, b) += 2. * projection_Ab.get(i + 1, a) *
                               projection_Ab.get(j + 1, b) *
                               weyl_prop_minus.get(i, j);
            }
          }
        }
      }
    }
  }

  if (incoming_wave_profile != nullptr) {
    if constexpr (VolumeDim == 3) {
      // The configured profile is the strain g(t) of the incoming wave at the
      // boundary. The slot below is the incoming Weyl mode 2 U^{8-}, which
      // for an incoming plane wave is -2 d^2 h/dt^2, so the injected tensor
      // is -2 g''(t) h_ij: the strain, its rate and the curvature all follow
      // g and return to zero once a pulse has passed.
      const double injected_wave_profile_value =
          -2. * incoming_wave_profile->second_deriv(time);
      // The spatial block of the injected wave, in the storage order
      // (xx, xy, xz, yy, yz, zz). The transverse-traceless projection applied
      // below keeps only the part that is transverse and trace free with
      // respect to the boundary normal, so the components supplied here need
      // be neither.
      size_t component = 0;
      for (size_t i = 0; i < VolumeDim; ++i) {
        for (size_t j = i; j < VolumeDim; ++j) {
          injected_wave.get(i + 1, j + 1) =
              gsl::at(incoming_wave_components, component) *
              injected_wave_profile_value;
          ++component;
        }
      }
    } else {
      ERROR("IncomingWaveProfile can only be used in 3 spatial dimensions.");
    }
  }

  // Add physical boundary corrections
  if (gamma2_in_phys) {
    auto& normal_dot_three_index_constraint_gamma2 =
        get(get<::Tags::TempScalar<0, DataType>>(u3_buffer));
    for (size_t a = 0; a <= VolumeDim; ++a) {
      for (size_t b = a; b <= VolumeDim; ++b) {
        for (size_t c = 0; c <= VolumeDim; ++c) {
          for (size_t d = 0; d <= VolumeDim; ++d) {
            normal_dot_three_index_constraint_gamma2 =
                get<0>(unit_interface_normal_vector) *
                three_index_constraint.get(0, c, d);
            for (size_t i = 1; i < VolumeDim; ++i) {
              normal_dot_three_index_constraint_gamma2 +=
                  unit_interface_normal_vector.get(i) *
                  three_index_constraint.get(i, c, d);
            }
            normal_dot_three_index_constraint_gamma2 *= get(gamma2);

            if constexpr (mu_phys == 0.) {
              bc_dt_v_minus->get(a, b) +=
                  (projection_Ab.get(c, a) * projection_Ab.get(d, b) -
                   0.5 * projection_ab.get(a, b) * projection_AB.get(c, d)) *
                  (char_projected_rhs_dt_v_minus.get(c, d) +
                   char_speeds[3] * (U3m.get(c, d) - injected_wave.get(c, d) -
                                     normal_dot_three_index_constraint_gamma2));
            } else {
              bc_dt_v_minus->get(a, b) +=
                  (projection_Ab.get(c, a) * projection_Ab.get(d, b) -
                   0.5 * projection_ab.get(a, b) * projection_AB.get(c, d)) *
                  (char_projected_rhs_dt_v_minus.get(c, d) +
                   char_speeds[3] * (U3m.get(c, d) -
                                     normal_dot_three_index_constraint_gamma2 -
                                     mu_phys * U3p.get(c, d)));
            }
          }
        }
      }
    }
  } else {
    for (size_t a = 0; a <= VolumeDim; ++a) {
      for (size_t b = a; b <= VolumeDim; ++b) {
        for (size_t c = 0; c <= VolumeDim; ++c) {
          for (size_t d = 0; d <= VolumeDim; ++d) {
            if constexpr (mu_phys == 0.) {
              bc_dt_v_minus->get(a, b) +=
                  (projection_Ab.get(c, a) * projection_Ab.get(d, b) -
                   0.5 * projection_ab.get(a, b) * projection_AB.get(c, d)) *
                  (char_projected_rhs_dt_v_minus.get(c, d) +
                   char_speeds[3] * (U3m.get(c, d)));
            } else {
              bc_dt_v_minus->get(a, b) +=
                  (projection_Ab.get(c, a) * projection_Ab.get(d, b) -
                   0.5 * projection_ab.get(a, b) * projection_AB.get(c, d)) *
                  (char_projected_rhs_dt_v_minus.get(c, d) +
                   char_speeds[3] * (U3m.get(c, d) - mu_phys * U3p.get(c, d)));
            }
          }
        }
      }
    }
  }
}
}  // namespace detail

template <size_t VolumeDim, typename DataType>
void constraint_preserving_corrections_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const std::array<DataType, 4>& char_speeds) {
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      bc_dt_v_minus->get(a, b) = -char_projected_rhs_dt_v_minus.get(a, b);
    }
  }
  detail::add_constraint_dependent_terms_to_dt_v_minus(
      bc_dt_v_minus, outgoing_null_one_form, incoming_null_vector,
      outgoing_null_vector, projection_ab, projection_Ab, projection_AB,
      constraint_char_zero_plus, constraint_char_zero_minus,
      char_projected_rhs_dt_v_minus, char_speeds);
}

template <size_t VolumeDim, typename DataType>
void constraint_preserving_gauge_corrections_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const std::array<DataType, 4>& char_speeds) {
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      bc_dt_v_minus->get(a, b) = -char_projected_rhs_dt_v_minus.get(a, b);
    }
  }
  detail::add_constraint_dependent_terms_to_dt_v_minus(
      bc_dt_v_minus, outgoing_null_one_form, incoming_null_vector,
      outgoing_null_vector, projection_ab, projection_Ab, projection_AB,
      constraint_char_zero_plus, constraint_char_zero_minus,
      char_projected_rhs_dt_v_minus, char_speeds);
  detail::add_gauge_sommerfeld_terms_to_dt_v_minus(
      bc_dt_v_minus, gamma2, inertial_coords, incoming_null_one_form,
      outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
      projection_Ab, char_projected_rhs_dt_v_psi);
}

template <size_t VolumeDim, typename DataType>
void constraint_preserving_gauge_physical_corrections_dt_v_minus(
    const gsl::not_null<tnsr::aa<DataType, VolumeDim, Frame::Inertial>*>
        bc_dt_v_minus,
    const Scalar<DataType>& gamma2,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>& inertial_coords,
    const double time,
    const tnsr::i<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_one_form,
    const tnsr::I<DataType, VolumeDim, Frame::Inertial>&
        unit_interface_normal_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>&
        spacetime_unit_normal_vector,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& incoming_null_one_form,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>& outgoing_null_one_form,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& incoming_null_vector,
    const tnsr::A<DataType, VolumeDim, Frame::Inertial>& outgoing_null_vector,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& projection_ab,
    const tnsr::Ab<DataType, VolumeDim, Frame::Inertial>& projection_Ab,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>& projection_AB,
    const tnsr::II<DataType, VolumeDim, Frame::Inertial>&
        inverse_spatial_metric,
    const tnsr::ii<DataType, VolumeDim, Frame::Inertial>& extrinsic_curvature,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>& spacetime_metric,
    const tnsr::AA<DataType, VolumeDim, Frame::Inertial>&
        inverse_spacetime_metric,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>&
        three_index_constraint,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_psi,
    const tnsr::aa<DataType, VolumeDim, Frame::Inertial>&
        char_projected_rhs_dt_v_minus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_plus,
    const tnsr::a<DataType, VolumeDim, Frame::Inertial>&
        constraint_char_zero_minus,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& phi,
    const tnsr::ijaa<DataType, VolumeDim, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataType, VolumeDim, Frame::Inertial>& d_pi,
    const std::array<DataType, 4>& char_speeds,
    const MathFunction<1, Frame::Inertial>* const incoming_wave_profile,
    const std::array<double, 6>& incoming_wave_components) {
  for (size_t a = 0; a <= VolumeDim; ++a) {
    for (size_t b = a; b <= VolumeDim; ++b) {
      bc_dt_v_minus->get(a, b) = -char_projected_rhs_dt_v_minus.get(a, b);
    }
  }
  detail::add_constraint_dependent_terms_to_dt_v_minus(
      bc_dt_v_minus, outgoing_null_one_form, incoming_null_vector,
      outgoing_null_vector, projection_ab, projection_Ab, projection_AB,
      constraint_char_zero_plus, constraint_char_zero_minus,
      char_projected_rhs_dt_v_minus, char_speeds);
  detail::add_physical_terms_to_dt_v_minus(
      bc_dt_v_minus, gamma2, unit_interface_normal_one_form,
      unit_interface_normal_vector, spacetime_unit_normal_vector, projection_ab,
      projection_Ab, projection_AB, inverse_spatial_metric, extrinsic_curvature,
      spacetime_metric, inverse_spacetime_metric, three_index_constraint,
      char_projected_rhs_dt_v_minus, phi, d_phi, d_pi, char_speeds, time,
      incoming_wave_profile, incoming_wave_components);
  detail::add_gauge_sommerfeld_terms_to_dt_v_minus(
      bc_dt_v_minus, gamma2, inertial_coords, incoming_null_one_form,
      outgoing_null_one_form, incoming_null_vector, outgoing_null_vector,
      projection_Ab, char_projected_rhs_dt_v_psi);
}

template <size_t Dim>
IntermediateVariables<Dim>::IntermediateVariables(const size_t num_points)
    : spacetime_unit_normal_one_form(num_points, 0.),
      unit_interface_normal_vector(num_points, 0.),
      four_index_constraint(num_points, 0.),
      inverse_spatial_metric(num_points, 0.),
      extrinsic_curvature(num_points, 0.),
      incoming_null_one_form(num_points, 0.),
      outgoing_null_one_form(num_points, 0.),
      incoming_null_vector(num_points, 0.),
      outgoing_null_vector(num_points, 0.),
      projection_ab(num_points, 0.),
      projection_Ab(num_points, 0.),
      projection_AB(num_points, 0.),
      dt_spacetime_metric(num_points, 0.),
      dt_pi(num_points, 0.),
      dt_phi(num_points, 0.),
      char_projected_rhs_dt_v_psi(num_points, 0.),
      char_projected_rhs_dt_v_zero(num_points, 0.),
      char_projected_rhs_dt_v_plus(num_points, 0.),
      char_projected_rhs_dt_v_minus(num_points, 0.),
      constraint_char_zero_plus(num_points, 0.),
      constraint_char_zero_minus(num_points, 0.),
      char_speeds{} {}

template <size_t Dim>
void compute_intermediate_variables(
    const gsl::not_null<IntermediateVariables<Dim>*> vars,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        face_mesh_velocity,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
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
    const tnsr::ijaa<DataVector, Dim, Frame::Inertial>& d_phi) {
  // The unit normal one-form to the spatial slice has only a time component
  get<0>(vars->spacetime_unit_normal_one_form) = -get(lapse);
  for (size_t i = 1; i <= Dim; ++i) {
    vars->spacetime_unit_normal_one_form.get(i) = 0.;
  }

  // Inertial time derivatives: subtract the advective mesh-velocity term from
  // the logical time derivatives
  vars->dt_spacetime_metric = logical_dt_spacetime_metric;
  vars->dt_pi = logical_dt_pi;
  vars->dt_phi = logical_dt_phi;
  if (face_mesh_velocity.has_value()) {
    for (size_t a = 0; a < Dim + 1; ++a) {
      for (size_t b = a; b < Dim + 1; ++b) {
        for (size_t d = 0; d < Dim; ++d) {
          vars->dt_spacetime_metric.get(a, b) -=
              face_mesh_velocity->get(d) * d_spacetime_metric.get(d, a, b);
          vars->dt_pi.get(a, b) -=
              face_mesh_velocity->get(d) * d_pi.get(d, a, b);
        }
        for (size_t i = 0; i < Dim; ++i) {
          for (size_t d = 0; d < Dim; ++d) {
            vars->dt_phi.get(i, a, b) -=
                face_mesh_velocity->get(d) * d_phi.get(d, i, a, b);
          }
        }
      }
    }
  }

  const DataVector one_over_lapse_sqrd = 1.0 / (get(lapse) * get(lapse));
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      vars->inverse_spatial_metric.get(i, j) =
          inverse_spacetime_metric.get(1 + i, 1 + j) +
          (shift.get(i) * shift.get(j) * one_over_lapse_sqrd);
    }
  }

  raise_or_lower_index(make_not_null(&vars->unit_interface_normal_vector),
                       normal_covector, vars->inverse_spatial_metric);

  gh::extrinsic_curvature(make_not_null(&vars->extrinsic_curvature),
                          spacetime_unit_normal_vector, pi, phi);

  if (LIKELY(Dim == 3)) {
    gh::four_index_constraint(make_not_null(&vars->four_index_constraint),
                              d_phi);
  } else if (UNLIKELY(Dim == 2)) {
    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = 0; b <= Dim; ++b) {
        vars->four_index_constraint.get(0, a, b) =
            d_phi.get(0, 1, a, b) - d_phi.get(1, 0, a, b);
        vars->four_index_constraint.get(1, a, b) =
            -vars->four_index_constraint.get(0, a, b);
      }
    }
  } else {
    std::fill(vars->four_index_constraint.begin(),
              vars->four_index_constraint.end(), 0.);
  }

  gr::interface_null_normal(make_not_null(&vars->incoming_null_one_form),
                            vars->spacetime_unit_normal_one_form,
                            normal_covector, shift, -1.);
  gr::interface_null_normal(make_not_null(&vars->outgoing_null_one_form),
                            vars->spacetime_unit_normal_one_form,
                            normal_covector, shift, 1.);
  gr::interface_null_normal(make_not_null(&vars->incoming_null_vector),
                            spacetime_unit_normal_vector,
                            vars->unit_interface_normal_vector, -1.);
  gr::interface_null_normal(make_not_null(&vars->outgoing_null_vector),
                            spacetime_unit_normal_vector,
                            vars->unit_interface_normal_vector, 1.);

  gr::transverse_projection_operator(
      make_not_null(&vars->projection_ab), spacetime_metric,
      vars->spacetime_unit_normal_one_form, normal_covector, shift);
  gr::transverse_projection_operator(
      make_not_null(&vars->projection_Ab), spacetime_unit_normal_vector,
      vars->spacetime_unit_normal_one_form, vars->unit_interface_normal_vector,
      normal_covector, shift);
  gr::transverse_projection_operator(
      make_not_null(&vars->projection_AB), inverse_spacetime_metric,
      spacetime_unit_normal_vector, vars->unit_interface_normal_vector);

  const auto dt_char_fields = gh::characteristic_fields(
      gamma2, vars->inverse_spatial_metric, vars->dt_spacetime_metric,
      vars->dt_pi, vars->dt_phi, normal_covector);
  vars->char_projected_rhs_dt_v_psi =
      get<gh::Tags::VSpacetimeMetric<DataVector, Dim>>(dt_char_fields);
  vars->char_projected_rhs_dt_v_zero =
      get<gh::Tags::VZero<DataVector, Dim>>(dt_char_fields);
  vars->char_projected_rhs_dt_v_plus =
      get<gh::Tags::VPlus<DataVector, Dim>>(dt_char_fields);
  vars->char_projected_rhs_dt_v_minus =
      get<gh::Tags::VMinus<DataVector, Dim>>(dt_char_fields);

  // c^{\hat{0}-}_a = F_a + n^k C_{ka}
  tnsr::ia<DataVector, Dim, Frame::Inertial> two_index_constraint(
      get_size(get(lapse)), 0.);
  gh::two_index_constraint(
      make_not_null(&two_index_constraint), spacetime_deriv_gauge_source,
      vars->spacetime_unit_normal_one_form, spacetime_unit_normal_vector,
      vars->inverse_spatial_metric, inverse_spacetime_metric, pi, phi, d_pi,
      d_phi, gamma2, three_index_constraint);
  gh::f_constraint(make_not_null(&vars->constraint_char_zero_plus),
                   gauge_source, spacetime_deriv_gauge_source,
                   vars->spacetime_unit_normal_one_form,
                   spacetime_unit_normal_vector, vars->inverse_spatial_metric,
                   inverse_spacetime_metric, pi, phi, d_pi, d_phi, gamma2,
                   three_index_constraint);
  for (size_t a = 0; a < Dim + 1; ++a) {
    vars->constraint_char_zero_minus.get(a) =
        vars->constraint_char_zero_plus.get(a);
    for (size_t i = 0; i < Dim; ++i) {
      vars->constraint_char_zero_plus.get(a) -=
          vars->unit_interface_normal_vector.get(i) *
          two_index_constraint.get(i, a);
      vars->constraint_char_zero_minus.get(a) +=
          vars->unit_interface_normal_vector.get(i) *
          two_index_constraint.get(i, a);
    }
  }

  gh::characteristic_speeds(make_not_null(&vars->char_speeds), gamma1, lapse,
                            shift, normal_covector, face_mesh_velocity);
}

double min_characteristic_speed(const std::array<DataVector, 4>& char_speeds) {
  const std::array<double, 4> min_speeds{
      {min(char_speeds[0]), min(char_speeds[1]), min(char_speeds[2]),
       min(char_speeds[3])}};
  return *std::min_element(min_speeds.begin(), min_speeds.end());
}

namespace {
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

template <size_t Dim>
void project_corrections_onto_evolved_variables(
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        dt_spacetime_metric_correction,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        dt_pi_correction,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        dt_phi_correction,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        bc_dt_v_psi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim, Frame::Inertial>*>
        bc_dt_v_zero,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        bc_dt_v_plus,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame::Inertial>*>
        bc_dt_v_minus,
    const std::array<DataVector, 4>& char_speeds,
    const Scalar<DataVector>& gamma2,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector) {
  // Only add corrections at grid points where the char speeds are negative
  set_bc_corr_zero_when_char_speed_is_positive(bc_dt_v_psi, char_speeds[0]);
  set_bc_corr_zero_when_char_speed_is_positive(bc_dt_v_zero, char_speeds[1]);
  set_bc_corr_zero_when_char_speed_is_positive(bc_dt_v_plus, char_speeds[2]);
  set_bc_corr_zero_when_char_speed_is_positive(bc_dt_v_minus, char_speeds[3]);

  // The boundary conditions are imposed as corrections to the characteristic
  // projections of the right-hand sides of the GH evolution equations, Eqs.
  // (63) - (65) of Lindblom et al (2005). Project them back onto corrections
  // to dt<evolved variables>.
  const auto dt_evolved_vars = gh::evolved_fields_from_characteristic_fields(
      gamma2, *bc_dt_v_psi, *bc_dt_v_zero, *bc_dt_v_plus, *bc_dt_v_minus,
      normal_covector);
  *dt_pi_correction = get<gh::Tags::Pi<DataVector, Dim>>(dt_evolved_vars);
  *dt_phi_correction = get<gh::Tags::Phi<DataVector, Dim>>(dt_evolved_vars);
  *dt_spacetime_metric_correction =
      get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(dt_evolved_vars);
}
}  // namespace gh::BoundaryConditions::Bjorhus

// Explicit Instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                                   \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::constraint_preserving_corrections_dt_v_psi( \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>  \
          bc_dt_v_psi,                                                         \
      const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          unit_interface_normal_vector,                                        \
      const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>&                \
          three_index_constraint,                                              \
      const std::array<DTYPE(data), 4>& char_speeds);                          \
  template void gh::BoundaryConditions::Bjorhus::                              \
      constraint_preserving_corrections_dt_v_zero(                             \
          const gsl::not_null<                                                 \
              tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>*>             \
              bc_dt_v_zero,                                                    \
          const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&              \
              unit_interface_normal_vector,                                    \
          const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>&            \
              four_index_constraint,                                           \
          const std::array<DTYPE(data), 4>& char_speeds);                      \
  template void gh::BoundaryConditions::Bjorhus::detail::                      \
      add_gauge_sommerfeld_terms_to_dt_v_minus(                                \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const Scalar<DTYPE(data)>& gamma2,                                   \
          const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&              \
              inertial_coords,                                                 \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_one_form,                                          \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_psi);                                    \
  template void gh::BoundaryConditions::Bjorhus::detail::                      \
      add_constraint_dependent_terms_to_dt_v_minus(                            \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_ab,                                                   \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_AB,                                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_plus,                                       \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_minus,                                      \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_minus,                                   \
          const std::array<DTYPE(data), 4>& char_speeds);                      \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::detail::add_physical_terms_to_dt_v_minus(   \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>  \
          bc_dt_v_minus,                                                       \
      const Scalar<DTYPE(data)>& gamma2,                                       \
      const tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          unit_interface_normal_one_form,                                      \
      const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          unit_interface_normal_vector,                                        \
      const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          spacetime_unit_normal_vector,                                        \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& projection_ab,  \
      const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>& projection_Ab,  \
      const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>& projection_AB,  \
      const tnsr::II<DTYPE(data), DIM(data), Frame::Inertial>&                 \
          inverse_spatial_metric,                                              \
      const tnsr::ii<DTYPE(data), DIM(data), Frame::Inertial>&                 \
          extrinsic_curvature,                                                 \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&                 \
          spacetime_metric,                                                    \
      const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&                 \
          inverse_spacetime_metric,                                            \
      const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>&                \
          three_index_constraint,                                              \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&                 \
          char_projected_rhs_dt_v_minus,                                       \
      const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>& phi,           \
      const tnsr::ijaa<DTYPE(data), DIM(data), Frame::Inertial>& d_phi,        \
      const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>& d_pi,          \
      const std::array<DTYPE(data), 4>& char_speeds, const double time,        \
      const MathFunction<1, Frame::Inertial>* incoming_wave_profile,           \
      const std::array<double, 6>& incoming_wave_components,                   \
      const tnsr::ii<DTYPE(data), DIM(data), Frame::Inertial>*                 \
          incoming_weyl_propagating_minus);                                    \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::detail::add_constraint_sector_projection(   \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>  \
          result,                                                              \
      const DTYPE(data) & scalar_coefficient,                                  \
      const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          outgoing_null_one_form,                                              \
      const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          incoming_null_vector,                                                \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& projection_ab,  \
      const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>& projection_Ab,  \
      const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>& projection_AB,  \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& source);        \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::detail::add_physical_sector_projection(     \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>  \
          result,                                                              \
      const DTYPE(data) & scalar_coefficient,                                  \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& projection_ab,  \
      const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>& projection_Ab,  \
      const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>& projection_AB,  \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& source);        \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::detail::add_gauge_sector_projection(        \
      const gsl::not_null<tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>  \
          result,                                                              \
      const DTYPE(data) & scalar_coefficient,                                  \
      const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          incoming_null_one_form,                                              \
      const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          outgoing_null_one_form,                                              \
      const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          incoming_null_vector,                                                \
      const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&                  \
          outgoing_null_vector,                                                \
      const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>& projection_Ab,  \
      const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& source);        \
  template void gh::BoundaryConditions::Bjorhus::detail::                      \
      add_gauge_sector_terms_to_dt_v_minus(                                    \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const DTYPE(data) & scalar_coefficient,                              \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_one_form,                                          \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>& source);    \
  template void gh::BoundaryConditions::Bjorhus::                              \
      constraint_preserving_corrections_dt_v_minus(                            \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_ab,                                                   \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_AB,                                                   \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_minus,                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_plus,                                       \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_minus,                                      \
          const std::array<DTYPE(data), 4>& char_speeds);                      \
  template void gh::BoundaryConditions::Bjorhus::                              \
      constraint_preserving_gauge_corrections_dt_v_minus(                      \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const Scalar<DTYPE(data)>& gamma2,                                   \
          const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&              \
              inertial_coords,                                                 \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_one_form,                                          \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_ab,                                                   \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_AB,                                                   \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_psi,                                     \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_minus,                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_plus,                                       \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_minus,                                      \
          const std::array<DTYPE(data), 4>& char_speeds);                      \
  template void gh::BoundaryConditions::Bjorhus::                              \
      constraint_preserving_gauge_physical_corrections_dt_v_minus(             \
          const gsl::not_null<                                                 \
              tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>*>              \
              bc_dt_v_minus,                                                   \
          const Scalar<DTYPE(data)>& gamma2,                                   \
          const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&              \
              inertial_coords,                                                 \
          const double time,                                                   \
          const tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>&              \
              unit_interface_normal_one_form,                                  \
          const tnsr::I<DTYPE(data), DIM(data), Frame::Inertial>&              \
              unit_interface_normal_vector,                                    \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              spacetime_unit_normal_vector,                                    \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_one_form,                                          \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_one_form,                                          \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              incoming_null_vector,                                            \
          const tnsr::A<DTYPE(data), DIM(data), Frame::Inertial>&              \
              outgoing_null_vector,                                            \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_ab,                                                   \
          const tnsr::Ab<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_Ab,                                                   \
          const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&             \
              projection_AB,                                                   \
          const tnsr::II<DTYPE(data), DIM(data), Frame::Inertial>&             \
              inverse_spatial_metric,                                          \
          const tnsr::ii<DTYPE(data), DIM(data), Frame::Inertial>&             \
              extrinsic_curvature,                                             \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              spacetime_metric,                                                \
          const tnsr::AA<DTYPE(data), DIM(data), Frame::Inertial>&             \
              inverse_spacetime_metric,                                        \
          const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>&            \
              three_index_constraint,                                          \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_psi,                                     \
          const tnsr::aa<DTYPE(data), DIM(data), Frame::Inertial>&             \
              char_projected_rhs_dt_v_minus,                                   \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_plus,                                       \
          const tnsr::a<DTYPE(data), DIM(data), Frame::Inertial>&              \
              constraint_char_zero_minus,                                      \
          const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>& phi,       \
          const tnsr::ijaa<DTYPE(data), DIM(data), Frame::Inertial>& d_phi,    \
          const tnsr::iaa<DTYPE(data), DIM(data), Frame::Inertial>& d_pi,      \
          const std::array<DTYPE(data), 4>& char_speeds,                       \
          const MathFunction<1, Frame::Inertial>* incoming_wave_profile,       \
          const std::array<double, 6>& incoming_wave_components);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (DataVector))

#undef INSTANTIATE
#undef DTYPE

#define INSTANTIATE_FACE(_, data)                                              \
  template struct gh::BoundaryConditions::Bjorhus::IntermediateVariables<DIM(  \
      data)>;                                                                  \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::compute_intermediate_variables(             \
      const gsl::not_null<                                                     \
          gh::BoundaryConditions::Bjorhus::IntermediateVariables<DIM(data)>*>  \
          vars,                                                                \
      const std::optional<tnsr::I<DataVector, DIM(data), Frame::Inertial>>&    \
          face_mesh_velocity,                                                  \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& normal_covector,  \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&                  \
          spacetime_metric,                                                    \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>& pi,              \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>& phi,            \
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,      \
      const Scalar<DataVector>& lapse,                                         \
      const tnsr::I<DataVector, DIM(data), Frame::Inertial>& shift,            \
      const tnsr::AA<DataVector, DIM(data), Frame::Inertial>&                  \
          inverse_spacetime_metric,                                            \
      const tnsr::A<DataVector, DIM(data), Frame::Inertial>&                   \
          spacetime_unit_normal_vector,                                        \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>&                 \
          three_index_constraint,                                              \
      const tnsr::a<DataVector, DIM(data), Frame::Inertial>& gauge_source,     \
      const tnsr::ab<DataVector, DIM(data), Frame::Inertial>&                  \
          spacetime_deriv_gauge_source,                                        \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&                  \
          logical_dt_spacetime_metric,                                         \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>& logical_dt_pi,   \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>& logical_dt_phi, \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>&                 \
          d_spacetime_metric,                                                  \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>& d_pi,           \
      const tnsr::ijaa<DataVector, DIM(data), Frame::Inertial>& d_phi);        \
  template void                                                                \
  gh::BoundaryConditions::Bjorhus::project_corrections_onto_evolved_variables( \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), Frame::Inertial>*>   \
          dt_spacetime_metric_correction,                                      \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), Frame::Inertial>*>   \
          dt_pi_correction,                                                    \
      const gsl::not_null<tnsr::iaa<DataVector, DIM(data), Frame::Inertial>*>  \
          dt_phi_correction,                                                   \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), Frame::Inertial>*>   \
          bc_dt_v_psi,                                                         \
      const gsl::not_null<tnsr::iaa<DataVector, DIM(data), Frame::Inertial>*>  \
          bc_dt_v_zero,                                                        \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), Frame::Inertial>*>   \
          bc_dt_v_plus,                                                        \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), Frame::Inertial>*>   \
          bc_dt_v_minus,                                                       \
      const std::array<DataVector, 4>& char_speeds,                            \
      const Scalar<DataVector>& gamma2,                                        \
      const tnsr::i<DataVector, DIM(data), Frame::Inertial>& normal_covector);

GENERATE_INSTANTIATIONS(INSTANTIATE_FACE, (1, 2, 3))

#undef INSTANTIATE_FACE
#undef DIM
