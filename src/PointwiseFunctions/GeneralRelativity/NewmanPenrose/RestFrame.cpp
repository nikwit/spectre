// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
namespace {
using Complex4 = std::array<std::complex<double>, 4>;
}  // namespace

void principal_null_pair(const gsl::not_null<AdaptedFourVector*> outgoing,
                         const gsl::not_null<AdaptedFourVector*> incoming,
                         const Scalar<ComplexDataVector>& a_bar,
                         const Scalar<ComplexDataVector>& b) {
  const size_t num_points = get_size(get(a_bar));
  *outgoing = AdaptedFourVector(num_points, 0.);
  *incoming = AdaptedFourVector(num_points, 0.);
  const double one_over_sqrt2 = 1. / sqrt(2.);
  for (size_t p = 0; p < num_points; ++p) {
    // The NR tetrad of eq. surface-tetrad in adapted components
    Complex4 ell{{one_over_sqrt2, one_over_sqrt2, 0., 0.}};
    Complex4 kay{{one_over_sqrt2, -one_over_sqrt2, 0., 0.}};
    Complex4 emm{
        {0., 0., one_over_sqrt2, std::complex<double>{0., one_over_sqrt2}}};
    // Type I with tetrad parameter -a, a = conj(a_bar): l fixed
    {
      const std::complex<double> a = -std::conj(get(a_bar)[p]);
      Complex4 new_kay{};
      Complex4 new_emm{};
      for (size_t mu = 0; mu < 4; ++mu) {
        gsl::at(new_kay, mu) =
            gsl::at(kay, mu) + std::conj(a) * gsl::at(emm, mu) +
            a * std::conj(gsl::at(emm, mu)) + std::norm(a) * gsl::at(ell, mu);
        gsl::at(new_emm, mu) = gsl::at(emm, mu) + a * gsl::at(ell, mu);
      }
      kay = new_kay;
      emm = new_emm;
    }
    // Type II with tetrad parameter -b: k fixed
    {
      const std::complex<double> minus_b = -get(b)[p];
      Complex4 new_ell{};
      for (size_t mu = 0; mu < 4; ++mu) {
        gsl::at(new_ell, mu) = gsl::at(ell, mu) +
                               std::conj(minus_b) * gsl::at(emm, mu) +
                               minus_b * std::conj(gsl::at(emm, mu)) +
                               std::norm(minus_b) * gsl::at(kay, mu);
      }
      ell = new_ell;
    }
    for (size_t mu = 0; mu < 4; ++mu) {
      outgoing->get(mu)[p] = std::real(gsl::at(ell, mu));
      incoming->get(mu)[p] = std::real(gsl::at(kay, mu));
    }
  }
}

TangentBoostMember tangent_boost_member(const Scalar<ComplexDataVector>& a_bar,
                                        const Scalar<ComplexDataVector>& b) {
  const size_t num_points = get_size(get(a_bar));
  AdaptedFourVector ell{};
  AdaptedFourVector kay{};
  principal_null_pair(make_not_null(&ell), make_not_null(&kay), a_bar, b);
  TangentBoostMember member{TriadVector(num_points, 0.),
                            TriadVector(num_points, 0.),
                            Scalar<DataVector>(num_points, 0.)};
  const double one_over_sqrt2 = 1. / sqrt(2.);
  for (size_t p = 0; p < num_points; ++p) {
    if (get<0>(ell)[p] <= 0. or get<0>(kay)[p] <= 0.) {
      ERROR("Principal null legs are not future-directed at point " << p);
    }
    // Rescale the legs so their time components agree: r^0 = 0
    const double boost = sqrt(get<0>(kay)[p] / get<0>(ell)[p]);
    std::array<double, 4> t_leg{};
    std::array<double, 4> r_leg{};
    for (size_t mu = 0; mu < 4; ++mu) {
      const double ell_mu = boost * ell.get(mu)[p];
      const double kay_mu = kay.get(mu)[p] / boost;
      gsl::at(t_leg, mu) = (ell_mu + kay_mu) * one_over_sqrt2;
      gsl::at(r_leg, mu) = (ell_mu - kay_mu) * one_over_sqrt2;
    }
    const double lorentz_factor = t_leg[0];
    const double radial_norm =
        sqrt(square(r_leg[1]) + square(r_leg[2]) + square(r_leg[3]));
    double speed_squared = 0.;
    for (size_t i = 0; i < 3; ++i) {
      member.radial_direction.get(i)[p] = gsl::at(r_leg, i + 1) / radial_norm;
      member.transverse_velocity.get(i)[p] =
          gsl::at(t_leg, i + 1) / lorentz_factor;
      speed_squared += square(member.transverse_velocity.get(i)[p]);
    }
    if (speed_squared >= 1.) {
      ERROR("Decoded transverse velocity is not subluminal at point " << p);
    }
    get(member.lorentz_factor)[p] = lorentz_factor;
  }
  return member;
}

Scalar<DataVector> invariant_tanh_rapidity(
    const TangentBoostMember& member, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& spatial_gradient,
    const Scalar<DataVector>& time_derivative) {
  const size_t num_points = get_size(get(lapse));
  const RealMatrix cholesky_inverse =
      inverse_lower_triangular(cholesky_factor(spatial_metric));
  Scalar<DataVector> tanh_rapidity(num_points, 0.);
  for (size_t p = 0; p < num_points; ++p) {
    double normal_gradient = get(time_derivative)[p];
    for (size_t i = 0; i < 3; ++i) {
      normal_gradient -= shift.get(i)[p] * spatial_gradient.get(i)[p];
    }
    normal_gradient /= get(lapse)[p];
    // Covector components in the adapted triad: R L^{-1} D K
    std::array<double, 3> orthonormal{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        gsl::at(orthonormal, i) +=
            cholesky_inverse.get(i, j)[p] * spatial_gradient.get(j)[p];
      }
    }
    double radial_gradient = 0.;
    double tangential = 0.;
    for (size_t i = 0; i < 3; ++i) {
      double adapted = 0.;
      for (size_t j = 0; j < 3; ++j) {
        adapted += adapted_rotation.get(i, j)[p] * gsl::at(orthonormal, j);
      }
      radial_gradient += member.radial_direction.get(i)[p] * adapted;
      tangential += member.transverse_velocity.get(i)[p] * adapted;
    }
    get(tanh_rapidity)[p] = -get(member.lorentz_factor)[p] *
                            (normal_gradient + tangential) / radial_gradient;
  }
  return tanh_rapidity;
}

Scalar<DataVector> invariant_rapidity(
    const TangentBoostMember& member, const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3, Frame::Inertial>& spatial_gradient,
    const Scalar<DataVector>& time_derivative) {
  Scalar<DataVector> rapidity =
      invariant_tanh_rapidity(member, adapted_rotation, spatial_metric, lapse,
                              shift, spatial_gradient, time_derivative);
  for (size_t p = 0; p < get_size(get(rapidity)); ++p) {
    const double tanh_rapidity = get(rapidity)[p];
    if (std::abs(tanh_rapidity) >= 1.) {
      ERROR("Invariant-boost rapidity is not physical at point "
            << p << ": tanh(eta) = " << tanh_rapidity);
    }
    get(rapidity)[p] = std::atanh(tanh_rapidity);
  }
  return rapidity;
}
}  // namespace gr::np
