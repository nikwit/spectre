// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace gh::BoundaryConditions::detail {
/*!
 * The gauge projection of the characteristic RHS that enforces
 * d_t(l^b u^-_ab) = -rate l^b u^-_ab, including the derivatives of the
 * characteristic basis. Here d_t on the left is a full coordinate derivative.
 *
 * This helper assumes a fixed coordinate boundary and time-independent gamma2.
 * The metric RHS must include any boundary correction already selected for
 * v_psi. The returned tensor is purely in the gauge sector, and its contraction
 * with l^b supplies the required frozen-basis characteristic RHS. The caller
 * must subtract the old gauge RHS and mask outgoing characteristics.
 */
template <size_t Dim>
tnsr::aa<DataVector, Dim> algebraic_gauge_rhs(
    const tnsr::aa<DataVector, Dim> &metric,
    const tnsr::aa<DataVector, Dim> &pi, const tnsr::iaa<DataVector, Dim> &phi,
    const Scalar<DataVector> &gamma2,
    const tnsr::II<DataVector, Dim> &inverse_spatial_metric,
    const tnsr::AA<DataVector, Dim> &inverse_spacetime_metric,
    const tnsr::A<DataVector, Dim> &time_normal,
    const tnsr::I<DataVector, Dim> &face_normal,
    const tnsr::a<DataVector, Dim> &incoming_null_one_form,
    const tnsr::A<DataVector, Dim> &outgoing_null_vector,
    const tnsr::aa<DataVector, Dim> &dt_metric, const double rate,
    const tnsr::a<DataVector, Dim> *const initial_difference = nullptr,
    const double outgoing_rate = 0.,
    const tnsr::iaa<DataVector, Dim> *const reference_phi = nullptr,
    const double reference_rate = 0.) {
  const size_t size = get(gamma2).size();
  DataVector dt_metric_nn(size, 0.);
  DataVector dt_metric_tt(size, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      dt_metric_nn +=
          dt_metric.get(i + 1, j + 1) * face_normal.get(i) * face_normal.get(j);
    }
  }
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = 0; b <= Dim; ++b) {
      dt_metric_tt +=
          dt_metric.get(a, b) * time_normal.get(a) * time_normal.get(b);
    }
  }
  tnsr::I<DataVector, Dim> dt_face_normal(size, 0.);
  for (size_t i = 0; i < Dim; ++i) {
    dt_face_normal.get(i) = 0.5 * face_normal.get(i) * dt_metric_nn;
    for (size_t j = 0; j < Dim; ++j) {
      for (size_t k = 0; k < Dim; ++k) {
        dt_face_normal.get(i) -= inverse_spatial_metric.get(i, j) *
                                 dt_metric.get(j + 1, k + 1) *
                                 face_normal.get(k);
      }
    }
  }
  tnsr::A<DataVector, Dim> dt_outgoing_null_vector(size, 0.);
  for (size_t a = 0; a <= Dim; ++a) {
    auto &dt_l = dt_outgoing_null_vector.get(a);
    dt_l = -0.5 * time_normal.get(a) * dt_metric_tt;
    for (size_t b = 0; b <= Dim; ++b) {
      for (size_t c = 0; c <= Dim; ++c) {
        dt_l -= inverse_spacetime_metric.get(a, b) * time_normal.get(c) *
                dt_metric.get(b, c);
      }
    }
    if (a > 0) {
      dt_l += dt_face_normal.get(a - 1);
    }
    dt_l /= std::sqrt(2.);
  }
  // q_target = -rate q - dot(l).u^- + l.dot(n).Phi.
  tnsr::a<DataVector, Dim> q_target(size, 0.);
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = 0; b <= Dim; ++b) {
      DataVector u_minus = pi.get(a, b) - get(gamma2) * metric.get(a, b);
      DataVector normal_term(size, 0.);
      for (size_t i = 0; i < Dim; ++i) {
        u_minus -= face_normal.get(i) * phi.get(i, a, b);
        normal_term += dt_face_normal.get(i) * phi.get(i, a, b);
      }
      q_target.get(a) += outgoing_null_vector.get(b) * normal_term -
                         (rate * outgoing_null_vector.get(b) +
                          dt_outgoing_null_vector.get(b)) *
                             u_minus;
    }
  }
  // q_minus - q_plus = -2 l^b n^i Phi_iab: no gamma2 subtraction.
  // Add the outgoing-driven target to the full-basis derivative terms above.
  if (initial_difference != nullptr) {
    for (size_t a = 0; a <= Dim; ++a) {
      DataVector difference(size, 0.);
      for (size_t b = 0; b <= Dim; ++b) {
        for (size_t i = 0; i < Dim; ++i) {
          difference -= 2. * outgoing_null_vector.get(b) * face_normal.get(i) *
                        phi.get(i, a, b);
        }
      }
      q_target.get(a) -=
          outgoing_rate * (difference - initial_difference->get(a));
    }
  }
  // Stationary-reference subtraction. E_a = q_a + F_a where
  // F_a = gamma2 l^b g_ab + sqrt(2) l^b l^i barPhi_iab.
  // Enforce dt E = -reference_rate E by adding -rate E - dt F to dt q.
  // Both null vectors in F use the live metric; barPhi is fixed.
  if (reference_phi != nullptr) {
    for (size_t a = 0; a <= Dim; ++a) {
      DataVector residual(size, 0.);
      DataVector dt_f(size, 0.);
      for (size_t b = 0; b <= Dim; ++b) {
        DataVector u = pi.get(a, b);
        for (size_t i = 0; i < Dim; ++i) {
          u -= face_normal.get(i) * phi.get(i, a, b);
          u += std::sqrt(2.) * outgoing_null_vector.get(i + 1) *
               reference_phi->get(i, a, b);
          dt_f += std::sqrt(2.) *
                  (dt_outgoing_null_vector.get(b) *
                       outgoing_null_vector.get(i + 1) +
                   outgoing_null_vector.get(b) *
                       dt_outgoing_null_vector.get(i + 1)) *
                  reference_phi->get(i, a, b);
        }
        residual += outgoing_null_vector.get(b) * u;
        dt_f +=
            get(gamma2) * (dt_outgoing_null_vector.get(b) * metric.get(a, b) +
                           outgoing_null_vector.get(b) * dt_metric.get(a, b));
      }
      q_target.get(a) -= reference_rate * residual + dt_f;
    }
  }
  DataVector l_dot_q(size, 0.);
  for (size_t a = 0; a <= Dim; ++a) {
    l_dot_q += outgoing_null_vector.get(a) * q_target.get(a);
  }
  // The unique gauge-sector lift X_ab with l^b X_ab = q_target_a.
  tnsr::aa<DataVector, Dim> result(size, 0.);
  for (size_t a = 0; a <= Dim; ++a) {
    for (size_t b = a; b <= Dim; ++b) {
      result.get(a, b) = -incoming_null_one_form.get(a) * q_target.get(b) -
                         incoming_null_one_form.get(b) * q_target.get(a) -
                         incoming_null_one_form.get(a) *
                             incoming_null_one_form.get(b) * l_dot_q;
    }
  }
  return result;
}
} // namespace gh::BoundaryConditions::detail
