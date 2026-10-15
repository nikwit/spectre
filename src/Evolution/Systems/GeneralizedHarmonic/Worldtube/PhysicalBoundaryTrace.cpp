// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryTrace.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"

namespace gh::worldtube {
PhysicalScreen::PhysicalScreen(const tnsr::aa<DataVector, 3>& metric,
                               const tnsr::AA<DataVector, 3>& inverse_metric,
                               const tnsr::A<DataVector, 3>& time_normal,
                               const tnsr::I<DataVector, 3>& interface_normal)
    : lower(get<0, 0>(metric).size(), 0.),
      mixed(get<0, 0>(metric).size(), 0.),
      upper(get<0, 0>(metric).size(), 0.) {
  const size_t n = get<0, 0>(metric).size();
  tnsr::A<DataVector, 3> s(n, 0.);
  tnsr::a<DataVector, 3> s_lower(n, 0.), n_lower(n, 0.);
  for (size_t i = 0; i < 3; ++i) {
    s.get(i + 1) = interface_normal.get(i);
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      s_lower.get(a) += metric.get(a, b) * s.get(b);
      n_lower.get(a) += metric.get(a, b) * time_normal.get(b);
    }
  }
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      mixed.get(a, b) = (a == b ? 1. : 0.) +
                        time_normal.get(a) * n_lower.get(b) -
                        s.get(a) * s_lower.get(b);
    }
    for (size_t b = a; b < 4; ++b) {
      lower.get(a, b) = metric.get(a, b) + n_lower.get(a) * n_lower.get(b) -
                        s_lower.get(a) * s_lower.get(b);
      upper.get(a, b) = inverse_metric.get(a, b) +
                        time_normal.get(a) * time_normal.get(b) -
                        s.get(a) * s.get(b);
    }
  }
}

tnsr::aa<DataVector, 3> physical_projection(
    const tnsr::aa<DataVector, 3>& tensor, const PhysicalScreen& screen) {
  const size_t n = get<0, 0>(tensor).size();
  tnsr::aa<DataVector, 3> result(n, 0.);
  BoundaryConditions::Bjorhus::detail::add_physical_sector_projection(
      make_not_null(&result), DataVector(n, 1.), screen.lower, screen.mixed,
      screen.upper, tensor);
  return result;
}

tnsr::aa<DataVector, 3> incoming_weyl_to_physical_target(
    const tnsr::ii<DataVector, 3>& incoming_mode,
    const PhysicalScreen& receiver_screen) {
  // Embedding a spatial covariant tensor with zero temporal components and
  // applying Q performs both the spacetime lift and the screen detracing.
  tnsr::aa<DataVector, 3> embedded(get<0, 0>(incoming_mode).size(), 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      embedded.get(i + 1, j + 1) = 2. * incoming_mode.get(i, j);
    }
  }
  return physical_projection(embedded, receiver_screen);
}

void physical_ghost(const gsl::not_null<tnsr::aa<DataVector, 3>*> ghost_pi,
                    const gsl::not_null<tnsr::iaa<DataVector, 3>*> ghost_phi,
                    const tnsr::aa<DataVector, 3>& pi,
                    const tnsr::iaa<DataVector, 3>& phi,
                    const tnsr::aa<DataVector, 3>& target,
                    const tnsr::i<DataVector, 3>& normal_covector,
                    const tnsr::I<DataVector, 3>& normal_vector,
                    const PhysicalScreen& screen) {
  auto difference = target;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      difference.get(a, b) -= pi.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        difference.get(a, b) += normal_vector.get(i) * phi.get(i, a, b);
      }
    }
  }
  const auto jump = physical_projection(difference, screen);
  *ghost_pi = pi;
  *ghost_phi = phi;
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      ghost_pi->get(a, b) += 0.5 * jump.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        ghost_phi->get(i, a, b) -=
            0.5 * normal_covector.get(i) * jump.get(a, b);
      }
    }
  }
}

std::array<DataVector, 2> physical_trace_rhs(
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::aa<DataVector, 3>& rhs_pi,
    const tnsr::iaa<DataVector, 3>& rhs_phi,
    const tnsr::iaa<DataVector, 3>& three_index_constraint,
    const tnsr::aa<DataVector, 3>& weyl,
    const tnsr::aa<DataVector, 3>& target_weyl,
    const Scalar<DataVector>& gamma2, const DataVector& lambda_minus,
    const tnsr::I<DataVector, 3>& normal_vector,
    const tnsr::I<DataVector, 3>& dt_normal_vector,
    const std::array<tnsr::AA<DataVector, 3>, 2>& polarization,
    const std::array<tnsr::AA<DataVector, 3>, 2>& dt_polarization) {
  const size_t n = get<0, 0>(pi).size();
  std::array<DataVector, 2> result{DataVector(n, 0.), DataVector(n, 0.)};
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      DataVector v = pi.get(a, b);
      DataVector rhs = rhs_pi.get(a, b) +
                       lambda_minus * (weyl.get(a, b) - target_weyl.get(a, b));
      for (size_t i = 0; i < 3; ++i) {
        v -= normal_vector.get(i) * phi.get(i, a, b);
        rhs -=
            normal_vector.get(i) * (rhs_phi.get(i, a, b) +
                                    lambda_minus * get(gamma2) *
                                        three_index_constraint.get(i, a, b)) +
            dt_normal_vector.get(i) * phi.get(i, a, b);
      }
      const double multiplicity = a == b ? 1. : 2.;
      for (size_t p = 0; p < 2; ++p) {
        result[p] += multiplicity * (polarization[p].get(a, b) * rhs +
                                     dt_polarization[p].get(a, b) * v);
      }
    }
  }
  return result;
}
}  // namespace gh::worldtube
