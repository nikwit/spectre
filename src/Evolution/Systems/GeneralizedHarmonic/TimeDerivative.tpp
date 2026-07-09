// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/DiscontinuousGalerkin/TimeDerivativeDecisions.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/DuDtTempTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Dispatch.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Gauges.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpacetimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace gh {
namespace TimeDerivative_detail {
// Repoints the components of `view` at [offset, offset + extent) of the
// corresponding components of `full`.
template <typename TensorType>
void set_block_view(const gsl::not_null<TensorType*> view, TensorType& full,
                    const size_t offset, const size_t extent) {
  for (size_t c = 0; c < full.size(); ++c) {
    (*view)[c].set_data_ref(full[c].data() + offset, extent);
  }
}

template <typename TensorType>
void set_const_block_view(const TensorType& view, const TensorType& full,
                          const size_t offset, const size_t extent) {
  for (size_t c = 0; c < full.size(); ++c) {
    make_const_view(make_not_null(&view[c]), full[c], offset, extent);
  }
}
}  // namespace TimeDerivative_detail

template <class AllSolutionsForChristoffelAnalytic, size_t Dim>
evolution::dg::TimeDerivativeDecisions<Dim>
TimeDerivative<AllSolutionsForChristoffelAnalytic, Dim>::apply_impl(
    const gsl::not_null<tnsr::aa<DataVector, Dim>*> dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*> dt_pi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim>*> dt_phi,
    const gsl::not_null<Scalar<DataVector>*> temp_gamma1,
    const gsl::not_null<Scalar<DataVector>*> temp_gamma2,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> gauge_function,
    const gsl::not_null<tnsr::ab<DataVector, Dim>*>
        spacetime_deriv_gauge_function,
    const gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
    const gsl::not_null<Scalar<DataVector>*> half_pi_two_normals,
    const gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
    const gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> pi_one_normal,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> gauge_constraint,
    const gsl::not_null<tnsr::i<DataVector, Dim>*> half_phi_two_normals,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*>
        shift_dot_three_index_constraint,
    const gsl::not_null<tnsr::aa<DataVector, Dim>*>
        mesh_velocity_dot_three_index_constraint,
    const gsl::not_null<tnsr::ia<DataVector, Dim>*> phi_one_normal,
    const gsl::not_null<tnsr::aB<DataVector, Dim>*> pi_2_up,
    const gsl::not_null<tnsr::iaa<DataVector, Dim>*> three_index_constraint,
    const gsl::not_null<tnsr::Iaa<DataVector, Dim>*> phi_1_up,
    const gsl::not_null<tnsr::iaB<DataVector, Dim>*> phi_3_up,
    const gsl::not_null<tnsr::abC<DataVector, Dim>*>
        christoffel_first_kind_3_up,
    const gsl::not_null<Scalar<DataVector>*> lapse,
    const gsl::not_null<tnsr::I<DataVector, Dim>*> shift,
    const gsl::not_null<tnsr::II<DataVector, Dim>*> inverse_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
    const gsl::not_null<tnsr::AA<DataVector, Dim>*> inverse_spacetime_metric,
    const gsl::not_null<tnsr::abb<DataVector, Dim>*> christoffel_first_kind,
    const gsl::not_null<tnsr::a<DataVector, Dim>*> trace_christoffel,
    const gsl::not_null<tnsr::A<DataVector, Dim>*> normal_spacetime_vector,
    const tnsr::iaa<DataVector, Dim>& d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim>& d_pi,
    const tnsr::ijaa<DataVector, Dim>& d_phi,
    const tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::aa<DataVector, Dim>& pi, const tnsr::iaa<DataVector, Dim>& phi,
    const Scalar<DataVector>& gamma0, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2,
    const gauges::GaugeCondition& gauge_condition, const Mesh<Dim>& mesh,
    double time,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        mesh_velocity) {
  const size_t number_of_points = get<0, 0>(*dt_spacetime_metric).size();
  // Need constraint damping on interfaces in DG schemes
  *temp_gamma1 = gamma1;
  *temp_gamma2 = gamma2;

  const tnsr::ii<DataVector, Dim> spatial_metric{};
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = i; j < Dim; ++j) {
      make_const_view(make_not_null(&spatial_metric.get(i, j)),
                      spacetime_metric.get(i + 1, j + 1), 0, number_of_points);
    }
  }
  determinant_and_inverse(det_spatial_metric, inverse_spatial_metric,
                          spatial_metric);
  gr::shift(shift, spacetime_metric, *inverse_spatial_metric);
  gr::lapse(lapse, *shift, spacetime_metric);
  gr::inverse_spacetime_metric(inverse_spacetime_metric, *lapse, *shift,
                               *inverse_spatial_metric);
  // Compute the part of the dt_spacetime_metric equation that doesn't involve
  // constraints so we can use it for da_spacetime_metric to compute Christoffel
  // symbols.
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_spacetime_metric->get(mu, nu) = -get(*lapse) * pi.get(mu, nu);
      for (size_t m = 0; m < Dim; ++m) {
        dt_spacetime_metric->get(mu, nu) += shift->get(m) * phi.get(m, mu, nu);
      }
    }
  }

  const tnsr::abb<DataVector, Dim> da_spacetime_metric{};
  for (size_t a = 0; a < Dim + 1; ++a) {
    for (size_t b = a; b < Dim + 1; ++b) {
      make_const_view(make_not_null(&da_spacetime_metric.get(0, a, b)),
                      dt_spacetime_metric->get(a, b), 0, number_of_points);
      for (size_t i = 0; i < Dim; ++i) {
        make_const_view(make_not_null(&da_spacetime_metric.get(i + 1, a, b)),
                        phi.get(i, a, b), 0, number_of_points);
      }
    }
  }

  gr::christoffel_first_kind(christoffel_first_kind, da_spacetime_metric);
  gr::spacetime_normal_vector(normal_spacetime_vector, *lapse, *shift);

  get(*gamma1gamma2) = get(gamma1) * get(gamma2);
  const DataVector& gamma12 = get(*gamma1gamma2);

  for (size_t m = 0; m < Dim; ++m) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        phi_1_up->get(m, mu, nu) =
            inverse_spatial_metric->get(m, 0) * phi.get(0, mu, nu);
        for (size_t n = 1; n < Dim; ++n) {
          phi_1_up->get(m, mu, nu) +=
              inverse_spatial_metric->get(m, n) * phi.get(n, mu, nu);
        }
      }
    }
  }

  for (size_t m = 0; m < Dim; ++m) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
        phi_3_up->get(m, nu, alpha) =
            inverse_spacetime_metric->get(alpha, 0) * phi.get(m, nu, 0);
        for (size_t beta = 1; beta < Dim + 1; ++beta) {
          phi_3_up->get(m, nu, alpha) +=
              inverse_spacetime_metric->get(alpha, beta) * phi.get(m, nu, beta);
        }
      }
    }
  }

  for (size_t nu = 0; nu < Dim + 1; ++nu) {
    for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
      pi_2_up->get(nu, alpha) =
          inverse_spacetime_metric->get(alpha, 0) * pi.get(nu, 0);
      for (size_t beta = 1; beta < Dim + 1; ++beta) {
        pi_2_up->get(nu, alpha) +=
            inverse_spacetime_metric->get(alpha, beta) * pi.get(nu, beta);
      }
    }
  }

  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
        christoffel_first_kind_3_up->get(mu, nu, alpha) =
            inverse_spacetime_metric->get(alpha, 0) *
            christoffel_first_kind->get(mu, nu, 0);
        for (size_t beta = 1; beta < Dim + 1; ++beta) {
          christoffel_first_kind_3_up->get(mu, nu, alpha) +=
              inverse_spacetime_metric->get(alpha, beta) *
              christoffel_first_kind->get(mu, nu, beta);
        }
      }
      // The trace of the Christoffel symbol of the first kind,
      // \f$g^{bc}\Gamma_{abc} = \Gamma_{ab}{}^{b}\f$, is the diagonal (in the
      // traced-and-raised index) of `christoffel_first_kind_3_up`. Accumulate
      // it here to avoid a separate full contraction over the metric.
      if (nu == 0) {
        gauge_constraint->get(mu) = christoffel_first_kind_3_up->get(mu, 0, 0);
      } else {
        gauge_constraint->get(mu) +=
            christoffel_first_kind_3_up->get(mu, nu, nu);
      }
    }
  }

  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    pi_one_normal->get(mu) = get<0>(*normal_spacetime_vector) * pi.get(0, mu);
    for (size_t nu = 1; nu < Dim + 1; ++nu) {
      pi_one_normal->get(mu) +=
          normal_spacetime_vector->get(nu) * pi.get(nu, mu);
    }
  }

  get(*half_pi_two_normals) =
      get<0>(*normal_spacetime_vector) * get<0>(*pi_one_normal);
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    get(*half_pi_two_normals) +=
        normal_spacetime_vector->get(mu) * pi_one_normal->get(mu);
  }
  get(*half_pi_two_normals) *= 0.5;

  for (size_t n = 0; n < Dim; ++n) {
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      phi_one_normal->get(n, nu) =
          get<0>(*normal_spacetime_vector) * phi.get(n, 0, nu);
      for (size_t mu = 1; mu < Dim + 1; ++mu) {
        phi_one_normal->get(n, nu) +=
            normal_spacetime_vector->get(mu) * phi.get(n, mu, nu);
      }
    }
  }

  for (size_t n = 0; n < Dim; ++n) {
    half_phi_two_normals->get(n) =
        get<0>(*normal_spacetime_vector) * phi_one_normal->get(n, 0);
    for (size_t mu = 1; mu < Dim + 1; ++mu) {
      half_phi_two_normals->get(n) +=
          normal_spacetime_vector->get(mu) * phi_one_normal->get(n, mu);
    }
    half_phi_two_normals->get(n) *= 0.5;
  }

  for (size_t n = 0; n < Dim; ++n) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        three_index_constraint->get(n, mu, nu) =
            d_spacetime_metric.get(n, mu, nu) - phi.get(n, mu, nu);
      }
    }
  }

  get(*gamma1_plus_1) = 1.0 + gamma1.get();
  const DataVector& gamma1p1 = get(*gamma1_plus_1);

  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      shift_dot_three_index_constraint->get(mu, nu) =
          get<0>(*shift) * three_index_constraint->get(0, mu, nu);
      if (mesh_velocity.has_value()) {
        mesh_velocity_dot_three_index_constraint->get(mu, nu) =
            get<0>(*mesh_velocity) * three_index_constraint->get(0, mu, nu);
      }
      for (size_t m = 1; m < Dim; ++m) {
        shift_dot_three_index_constraint->get(mu, nu) +=
            shift->get(m) * three_index_constraint->get(m, mu, nu);
        if (mesh_velocity.has_value()) {
          mesh_velocity_dot_three_index_constraint->get(mu, nu) +=
              mesh_velocity->get(m) * three_index_constraint->get(m, mu, nu);
        }
      }
    }
  }

  const bool using_harmonic_gauge = gauge_condition.is_harmonic();
  if (not using_harmonic_gauge) {
    // Compute gauge condition.
    get(*sqrt_det_spatial_metric) = sqrt(get(*det_spatial_metric));
  }
  gauges::dispatch<AllSolutionsForChristoffelAnalytic, Dim>(
      gauge_function, spacetime_deriv_gauge_function, *lapse, *shift,
      *sqrt_det_spatial_metric, *inverse_spatial_metric, da_spacetime_metric,
      *half_pi_two_normals, *half_phi_two_normals, spacetime_metric, phi, mesh,
      time, inertial_coords, inverse_jacobian, gauge_condition);
  if (not using_harmonic_gauge) {
    // Compute source function last so that we don't need to recompute any of
    // the other temporary tags.
    for (size_t nu = 0; nu < Dim + 1; ++nu) {
      gauge_constraint->get(nu) += gauge_function->get(nu);
    }

    // Reuse `trace_christoffel` as scratch space for 2 H^a.  This avoids
    // constructing the full Christoffel symbol of the second kind just to
    // contract it with H_a below.
    for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
      trace_christoffel->get(alpha) = 2.0 *
                                      inverse_spacetime_metric->get(alpha, 0) *
                                      gauge_function->get(0);
      for (size_t beta = 1; beta < Dim + 1; ++beta) {
        trace_christoffel->get(alpha) +=
            2.0 * inverse_spacetime_metric->get(alpha, beta) *
            gauge_function->get(beta);
      }
    }
  }

  get(*normal_dot_gauge_constraint) =
      get<0>(*normal_spacetime_vector) * get<0>(*gauge_constraint);
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    get(*normal_dot_gauge_constraint) +=
        normal_spacetime_vector->get(mu) * gauge_constraint->get(mu);
  }

  // Here are the actual equations

  // Equation for dt_spacetime_metric
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_spacetime_metric->get(mu, nu) +=
          gamma1p1 * shift_dot_three_index_constraint->get(mu, nu);
      if (mesh_velocity.has_value()) {
        dt_spacetime_metric->get(mu, nu) +=
            get(gamma1) * mesh_velocity_dot_three_index_constraint->get(mu, nu);
      }
    }
  }

  // Equation for dt_pi

  // We first compute the n_a contributions but only for a=0 since n_i=0
  // identically. We also use dt_Pi_{00} as temporary storage to avoid any extra
  // allocations and multiply normal_dot_gauge_constraint=(n^a C_a) by gamma0
  // since it always shows up multiplied by gamma0 in the equations. This
  // reduces the number of multiplications that are needed for the RHS
  // evaluation.
  //
  // WARNING: normal_dot_gauge_constraint is rescaled by gamma0!
  get(*normal_dot_gauge_constraint) *= get(gamma0);

  // Use dt_pi_{00} as temporary storage.
  get<0, 0>(*dt_pi) = -get(gamma0) * get(*lapse);
  for (size_t i = 1; i < Dim + 1; ++i) {
    dt_pi->get(0, i) =
        get<0, 0>(*dt_pi) * gauge_constraint->get(i) -
        get(*normal_dot_gauge_constraint) * spacetime_metric.get(0, i);
  }
  get<0, 0>(*dt_pi) =
      2.0 * get<0, 0>(*dt_pi) * get<0>(*gauge_constraint) -
      get(*normal_dot_gauge_constraint) * get<0, 0>(spacetime_metric);

  // Set space-space components
  for (size_t mu = 1; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_pi->get(mu, nu) =
          -get(*normal_dot_gauge_constraint) * spacetime_metric.get(mu, nu);
    }
  }

  // Add additional pieces to dt_pi that aren't just n_a*(stuff)
  for (size_t mu = 0; mu < Dim + 1; ++mu) {
    for (size_t nu = mu; nu < Dim + 1; ++nu) {
      dt_pi->get(mu, nu) -= get(*half_pi_two_normals) * pi.get(mu, nu);

      if (not using_harmonic_gauge) {
        dt_pi->get(mu, nu) -= spacetime_deriv_gauge_function->get(mu, nu) +
                              spacetime_deriv_gauge_function->get(nu, mu);
      }
      for (size_t delta = 0; delta < Dim + 1; ++delta) {
        dt_pi->get(mu, nu) -= 2 * pi.get(mu, delta) * pi_2_up->get(nu, delta);
        if (not using_harmonic_gauge) {
          dt_pi->get(mu, nu) += christoffel_first_kind->get(delta, mu, nu) *
                                trace_christoffel->get(delta);
        }
        for (size_t n = 0; n < Dim; ++n) {
          dt_pi->get(mu, nu) +=
              2 * phi_1_up->get(n, mu, delta) * phi_3_up->get(n, nu, delta);
        }

        for (size_t alpha = 0; alpha < Dim + 1; ++alpha) {
          dt_pi->get(mu, nu) -=
              2. * christoffel_first_kind_3_up->get(mu, alpha, delta) *
              christoffel_first_kind_3_up->get(nu, delta, alpha);
        }
      }

      for (size_t m = 0; m < Dim; ++m) {
        dt_pi->get(mu, nu) -=
            pi_one_normal->get(m + 1) * phi_1_up->get(m, mu, nu);

        for (size_t n = 0; n < Dim; ++n) {
          dt_pi->get(mu, nu) -=
              inverse_spatial_metric->get(m, n) * d_phi.get(m, n, mu, nu);
        }
      }

      dt_pi->get(mu, nu) *= get(*lapse);

      dt_pi->get(mu, nu) +=
          gamma12 * shift_dot_three_index_constraint->get(mu, nu);
      if (mesh_velocity.has_value()) {
        dt_pi->get(mu, nu) +=
            gamma12 * mesh_velocity_dot_three_index_constraint->get(mu, nu);
      }

      for (size_t m = 0; m < Dim; ++m) {
        // DualFrame term
        dt_pi->get(mu, nu) += shift->get(m) * d_pi.get(m, mu, nu);
      }
    }
  }

  // Equation for dt_phi
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t mu = 0; mu < Dim + 1; ++mu) {
      for (size_t nu = mu; nu < Dim + 1; ++nu) {
        dt_phi->get(i, mu, nu) =
            pi.get(mu, nu) * half_phi_two_normals->get(i) -
            d_pi.get(i, mu, nu) +
            get(gamma2) * three_index_constraint->get(i, mu, nu);
        for (size_t n = 0; n < Dim; ++n) {
          dt_phi->get(i, mu, nu) +=
              phi_one_normal->get(i, n + 1) * phi_1_up->get(n, mu, nu);
        }

        dt_phi->get(i, mu, nu) *= get(*lapse);
        for (size_t m = 0; m < Dim; ++m) {
          dt_phi->get(i, mu, nu) += shift->get(m) * d_phi.get(m, i, mu, nu);
        }
      }
    }
  }
  return {true};
}

template <class AllSolutionsForChristoffelAnalytic, size_t Dim>
evolution::dg::TimeDerivativeDecisions<Dim>
TimeDerivative<AllSolutionsForChristoffelAnalytic, Dim>::apply(
    const gsl::not_null<tnsr::aa<DataVector, Dim> *> dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<DataVector, Dim> *> dt_pi,
    const gsl::not_null<tnsr::iaa<DataVector, Dim> *> dt_phi,
    const gsl::not_null<Scalar<DataVector> *> temp_gamma1,
    const gsl::not_null<Scalar<DataVector> *> temp_gamma2,
    const gsl::not_null<tnsr::a<DataVector, Dim> *> gauge_function,
    const gsl::not_null<tnsr::ab<DataVector, Dim> *>
        spacetime_deriv_gauge_function,
    const gsl::not_null<Scalar<DataVector> *> gamma1gamma2,
    const gsl::not_null<Scalar<DataVector> *> half_pi_two_normals,
    const gsl::not_null<Scalar<DataVector> *> normal_dot_gauge_constraint,
    const gsl::not_null<Scalar<DataVector> *> gamma1_plus_1,
    const gsl::not_null<tnsr::a<DataVector, Dim> *> pi_one_normal,
    const gsl::not_null<tnsr::a<DataVector, Dim> *> gauge_constraint,
    const gsl::not_null<tnsr::i<DataVector, Dim> *> half_phi_two_normals,
    const gsl::not_null<tnsr::aa<DataVector, Dim> *>
        shift_dot_three_index_constraint,
    const gsl::not_null<tnsr::aa<DataVector, Dim> *>
        mesh_velocity_dot_three_index_constraint,
    const gsl::not_null<tnsr::ia<DataVector, Dim> *> phi_one_normal,
    const gsl::not_null<tnsr::aB<DataVector, Dim> *> pi_2_up,
    const gsl::not_null<tnsr::iaa<DataVector, Dim> *> three_index_constraint,
    const gsl::not_null<tnsr::Iaa<DataVector, Dim> *> phi_1_up,
    const gsl::not_null<tnsr::iaB<DataVector, Dim> *> phi_3_up,
    const gsl::not_null<tnsr::abC<DataVector, Dim> *>
        christoffel_first_kind_3_up,
    const gsl::not_null<Scalar<DataVector> *> lapse,
    const gsl::not_null<tnsr::I<DataVector, Dim> *> shift,
    const gsl::not_null<tnsr::II<DataVector, Dim> *> inverse_spatial_metric,
    const gsl::not_null<Scalar<DataVector> *> det_spatial_metric,
    const gsl::not_null<Scalar<DataVector> *> sqrt_det_spatial_metric,
    const gsl::not_null<tnsr::AA<DataVector, Dim> *> inverse_spacetime_metric,
    const gsl::not_null<tnsr::abb<DataVector, Dim> *> christoffel_first_kind,
    const gsl::not_null<tnsr::a<DataVector, Dim> *> trace_christoffel,
    const gsl::not_null<tnsr::A<DataVector, Dim> *> normal_spacetime_vector,
    const tnsr::iaa<DataVector, Dim> &d_spacetime_metric,
    const tnsr::iaa<DataVector, Dim> &d_pi,
    const tnsr::ijaa<DataVector, Dim> &d_phi,
    const tnsr::aa<DataVector, Dim> &spacetime_metric,
    const tnsr::aa<DataVector, Dim> &pi, const tnsr::iaa<DataVector, Dim> &phi,
    const Scalar<DataVector> &gamma0, const Scalar<DataVector> &gamma1,
    const Scalar<DataVector> &gamma2,
    const gauges::GaugeCondition &gauge_condition, const Mesh<Dim> &mesh,
    const double time,
    const tnsr::I<DataVector, Dim, Frame::Inertial> &inertial_coords,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial> &inverse_jacobian,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>
        &mesh_velocity,
    const size_t block_size) {
  const size_t number_of_points = get<0, 0>(*dt_spacetime_metric).size();
  // Evaluating the time derivative on cache-sized blocks of grid points keeps
  // the ~30 temporaries resident in the core-private cache instead of
  // streaming them through shared cache levels between the elementwise
  // statements. This requires every operation to be pointwise, which holds
  // for the harmonic and damped-harmonic gauges but not for gauge conditions
  // that differentiate numerically (AnalyticChristoffel).
  const bool gauge_is_pointwise =
      gauge_condition.is_harmonic() or
      dynamic_cast<const gauges::DampedHarmonic *>(&gauge_condition) != nullptr;
  // Blocking is opt-in: on Apple M4 the fixed per-statement overheads of the
  // blaze expressions at block granularity outweigh the locality gain (see
  // kernel-optimization-findings.md section 13), so the default evaluates the
  // whole grid at once. Pass a block size explicitly to evaluate in blocks.
  const size_t target_block_size = block_size;
  if (block_size == 0 or not gauge_is_pointwise or
      number_of_points < 2 * target_block_size) {
    return apply_impl(
        dt_spacetime_metric, dt_pi, dt_phi, temp_gamma1, temp_gamma2,
        gauge_function, spacetime_deriv_gauge_function, gamma1gamma2,
        half_pi_two_normals, normal_dot_gauge_constraint, gamma1_plus_1,
        pi_one_normal, gauge_constraint, half_phi_two_normals,
        shift_dot_three_index_constraint,
        mesh_velocity_dot_three_index_constraint, phi_one_normal, pi_2_up,
        three_index_constraint, phi_1_up, phi_3_up, christoffel_first_kind_3_up,
        lapse, shift, inverse_spatial_metric, det_spatial_metric,
        sqrt_det_spatial_metric, inverse_spacetime_metric,
        christoffel_first_kind, trace_christoffel, normal_spacetime_vector,
        d_spacetime_metric, d_pi, d_phi, spacetime_metric, pi, phi, gamma0,
        gamma1, gamma2, gauge_condition, mesh, time, inertial_coords,
        inverse_jacobian, mesh_velocity);
  }
  const size_t number_of_blocks =
      (number_of_points + target_block_size - 1) / target_block_size;
  const size_t points_per_block =
      (number_of_points + number_of_blocks - 1) / number_of_blocks;

  using TimeDerivative_detail::set_block_view;
  using TimeDerivative_detail::set_const_block_view;

  // Views into every argument; their components are repointed per block.
  tnsr::aa<DataVector, Dim> dt_spacetime_metric_view{};
  tnsr::aa<DataVector, Dim> dt_pi_view{};
  tnsr::iaa<DataVector, Dim> dt_phi_view{};
  Scalar<DataVector> temp_gamma1_view{};
  Scalar<DataVector> temp_gamma2_view{};
  tnsr::a<DataVector, Dim> gauge_function_view{};
  tnsr::ab<DataVector, Dim> spacetime_deriv_gauge_function_view{};
  Scalar<DataVector> gamma1gamma2_view{};
  Scalar<DataVector> half_pi_two_normals_view{};
  Scalar<DataVector> normal_dot_gauge_constraint_view{};
  Scalar<DataVector> gamma1_plus_1_view{};
  tnsr::a<DataVector, Dim> pi_one_normal_view{};
  tnsr::a<DataVector, Dim> gauge_constraint_view{};
  tnsr::i<DataVector, Dim> half_phi_two_normals_view{};
  tnsr::aa<DataVector, Dim> shift_dot_three_index_constraint_view{};
  tnsr::aa<DataVector, Dim> mesh_velocity_dot_three_index_constraint_view{};
  tnsr::ia<DataVector, Dim> phi_one_normal_view{};
  tnsr::aB<DataVector, Dim> pi_2_up_view{};
  tnsr::iaa<DataVector, Dim> three_index_constraint_view{};
  tnsr::Iaa<DataVector, Dim> phi_1_up_view{};
  tnsr::iaB<DataVector, Dim> phi_3_up_view{};
  tnsr::abC<DataVector, Dim> christoffel_first_kind_3_up_view{};
  Scalar<DataVector> lapse_view{};
  tnsr::I<DataVector, Dim> shift_view{};
  tnsr::II<DataVector, Dim> inverse_spatial_metric_view{};
  Scalar<DataVector> det_spatial_metric_view{};
  Scalar<DataVector> sqrt_det_spatial_metric_view{};
  tnsr::AA<DataVector, Dim> inverse_spacetime_metric_view{};
  tnsr::abb<DataVector, Dim> christoffel_first_kind_view{};
  tnsr::a<DataVector, Dim> trace_christoffel_view{};
  tnsr::A<DataVector, Dim> normal_spacetime_vector_view{};
  const tnsr::iaa<DataVector, Dim> d_spacetime_metric_view{};
  const tnsr::iaa<DataVector, Dim> d_pi_view{};
  const tnsr::ijaa<DataVector, Dim> d_phi_view{};
  const tnsr::aa<DataVector, Dim> spacetime_metric_view{};
  const tnsr::aa<DataVector, Dim> pi_view{};
  const tnsr::iaa<DataVector, Dim> phi_view{};
  const Scalar<DataVector> gamma0_view{};
  const Scalar<DataVector> gamma1_view{};
  const Scalar<DataVector> gamma2_view{};
  const tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords_view{};
  std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> mesh_velocity_view{};
  if (mesh_velocity.has_value()) {
    mesh_velocity_view.emplace();
  }

  for (size_t offset = 0; offset < number_of_points;
       offset += points_per_block) {
    const size_t extent = std::min(points_per_block, number_of_points - offset);
    set_block_view(make_not_null(&dt_spacetime_metric_view),
                   *dt_spacetime_metric, offset, extent);
    set_block_view(make_not_null(&dt_pi_view), *dt_pi, offset, extent);
    set_block_view(make_not_null(&dt_phi_view), *dt_phi, offset, extent);
    set_block_view(make_not_null(&temp_gamma1_view), *temp_gamma1, offset,
                   extent);
    set_block_view(make_not_null(&temp_gamma2_view), *temp_gamma2, offset,
                   extent);
    set_block_view(make_not_null(&gauge_function_view), *gauge_function, offset,
                   extent);
    set_block_view(make_not_null(&spacetime_deriv_gauge_function_view),
                   *spacetime_deriv_gauge_function, offset, extent);
    set_block_view(make_not_null(&gamma1gamma2_view), *gamma1gamma2, offset,
                   extent);
    set_block_view(make_not_null(&half_pi_two_normals_view),
                   *half_pi_two_normals, offset, extent);
    set_block_view(make_not_null(&normal_dot_gauge_constraint_view),
                   *normal_dot_gauge_constraint, offset, extent);
    set_block_view(make_not_null(&gamma1_plus_1_view), *gamma1_plus_1, offset,
                   extent);
    set_block_view(make_not_null(&pi_one_normal_view), *pi_one_normal, offset,
                   extent);
    set_block_view(make_not_null(&gauge_constraint_view), *gauge_constraint,
                   offset, extent);
    set_block_view(make_not_null(&half_phi_two_normals_view),
                   *half_phi_two_normals, offset, extent);
    set_block_view(make_not_null(&shift_dot_three_index_constraint_view),
                   *shift_dot_three_index_constraint, offset, extent);
    set_block_view(
        make_not_null(&mesh_velocity_dot_three_index_constraint_view),
        *mesh_velocity_dot_three_index_constraint, offset, extent);
    set_block_view(make_not_null(&phi_one_normal_view), *phi_one_normal, offset,
                   extent);
    set_block_view(make_not_null(&pi_2_up_view), *pi_2_up, offset, extent);
    set_block_view(make_not_null(&three_index_constraint_view),
                   *three_index_constraint, offset, extent);
    set_block_view(make_not_null(&phi_1_up_view), *phi_1_up, offset, extent);
    set_block_view(make_not_null(&phi_3_up_view), *phi_3_up, offset, extent);
    set_block_view(make_not_null(&christoffel_first_kind_3_up_view),
                   *christoffel_first_kind_3_up, offset, extent);
    set_block_view(make_not_null(&lapse_view), *lapse, offset, extent);
    set_block_view(make_not_null(&shift_view), *shift, offset, extent);
    set_block_view(make_not_null(&inverse_spatial_metric_view),
                   *inverse_spatial_metric, offset, extent);
    set_block_view(make_not_null(&det_spatial_metric_view), *det_spatial_metric,
                   offset, extent);
    set_block_view(make_not_null(&sqrt_det_spatial_metric_view),
                   *sqrt_det_spatial_metric, offset, extent);
    set_block_view(make_not_null(&inverse_spacetime_metric_view),
                   *inverse_spacetime_metric, offset, extent);
    set_block_view(make_not_null(&christoffel_first_kind_view),
                   *christoffel_first_kind, offset, extent);
    set_block_view(make_not_null(&trace_christoffel_view), *trace_christoffel,
                   offset, extent);
    set_block_view(make_not_null(&normal_spacetime_vector_view),
                   *normal_spacetime_vector, offset, extent);
    set_const_block_view(d_spacetime_metric_view, d_spacetime_metric, offset,
                         extent);
    set_const_block_view(d_pi_view, d_pi, offset, extent);
    set_const_block_view(d_phi_view, d_phi, offset, extent);
    set_const_block_view(spacetime_metric_view, spacetime_metric, offset,
                         extent);
    set_const_block_view(pi_view, pi, offset, extent);
    set_const_block_view(phi_view, phi, offset, extent);
    set_const_block_view(gamma0_view, gamma0, offset, extent);
    set_const_block_view(gamma1_view, gamma1, offset, extent);
    set_const_block_view(gamma2_view, gamma2, offset, extent);
    set_const_block_view(inertial_coords_view, inertial_coords, offset, extent);
    if (mesh_velocity.has_value()) {
      set_const_block_view(*mesh_velocity_view, *mesh_velocity, offset, extent);
    }
    // The inverse Jacobian is only used by gauge conditions that
    // differentiate numerically, which are excluded from the blocked path,
    // so it is passed through unblocked.
    apply_impl(
        make_not_null(&dt_spacetime_metric_view), make_not_null(&dt_pi_view),
        make_not_null(&dt_phi_view), make_not_null(&temp_gamma1_view),
        make_not_null(&temp_gamma2_view), make_not_null(&gauge_function_view),
        make_not_null(&spacetime_deriv_gauge_function_view),
        make_not_null(&gamma1gamma2_view),
        make_not_null(&half_pi_two_normals_view),
        make_not_null(&normal_dot_gauge_constraint_view),
        make_not_null(&gamma1_plus_1_view), make_not_null(&pi_one_normal_view),
        make_not_null(&gauge_constraint_view),
        make_not_null(&half_phi_two_normals_view),
        make_not_null(&shift_dot_three_index_constraint_view),
        make_not_null(&mesh_velocity_dot_three_index_constraint_view),
        make_not_null(&phi_one_normal_view), make_not_null(&pi_2_up_view),
        make_not_null(&three_index_constraint_view),
        make_not_null(&phi_1_up_view), make_not_null(&phi_3_up_view),
        make_not_null(&christoffel_first_kind_3_up_view),
        make_not_null(&lapse_view), make_not_null(&shift_view),
        make_not_null(&inverse_spatial_metric_view),
        make_not_null(&det_spatial_metric_view),
        make_not_null(&sqrt_det_spatial_metric_view),
        make_not_null(&inverse_spacetime_metric_view),
        make_not_null(&christoffel_first_kind_view),
        make_not_null(&trace_christoffel_view),
        make_not_null(&normal_spacetime_vector_view), d_spacetime_metric_view,
        d_pi_view, d_phi_view, spacetime_metric_view, pi_view, phi_view,
        gamma0_view, gamma1_view, gamma2_view, gauge_condition, mesh, time,
        inertial_coords_view, inverse_jacobian, mesh_velocity_view);
  }
  return {true};
}
}  // namespace gh
