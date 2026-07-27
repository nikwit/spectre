// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/KerrSchildGrVars.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace CurvedScalarWave {

void zero_spin_kerr_schild_gr_vars(
    const gsl::not_null<Scalar<DataVector>*> lapse,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> deriv_lapse,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
    const gsl::not_null<tnsr::iJ<DataVector, 3, Frame::Inertial>*> deriv_shift,
    const gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*>
        spatial_metric,
    const gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        trace_spatial_christoffel_second_kind,
    const gsl::not_null<Scalar<DataVector>*> trace_extrinsic_curvature,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const double mass, const std::array<double, 3>& center) {
  const size_t num_points = get<0>(inertial_coords).size();
  set_number_of_grid_points(lapse, num_points);
  set_number_of_grid_points(deriv_lapse, num_points);
  set_number_of_grid_points(shift, num_points);
  set_number_of_grid_points(deriv_shift, num_points);
  set_number_of_grid_points(spatial_metric, num_points);
  set_number_of_grid_points(inverse_spatial_metric, num_points);
  set_number_of_grid_points(trace_spatial_christoffel_second_kind, num_points);
  set_number_of_grid_points(trace_extrinsic_curvature, num_points);
  // The computation is split into passes over the grid points with only a
  // few contiguous data streams each, which the hardware prefetcher handles
  // well. A single fused loop writing all 32 output components at once is
  // an order of magnitude slower as soon as the data does not fit in L1.
  TempBuffer<tmpl::list<
      ::Tags::TempI<0, 3, Frame::Inertial, DataVector>,
      ::Tags::TempScalar<1, DataVector>, ::Tags::TempScalar<2, DataVector>,
      ::Tags::TempScalar<3, DataVector>, ::Tags::TempScalar<4, DataVector>>>
      buffer{num_points};
  auto& n = get<::Tags::TempI<0, 3, Frame::Inertial, DataVector>>(buffer);
  // phi = 2M/r
  DataVector& phi = get(get<::Tags::TempScalar<1, DataVector>>(buffer));
  DataVector& one_over_r = get(get<::Tags::TempScalar<2, DataVector>>(buffer));
  DataVector& one_over_one_plus_phi =
      get(get<::Tags::TempScalar<3, DataVector>>(buffer));
  DataVector& factor = get(get<::Tags::TempScalar<4, DataVector>>(buffer));

  // n_i = (x_i - c_i)/r and the scalar factors
  for (size_t i = 0; i < 3; ++i) {
    n.get(i) = inertial_coords.get(i) - gsl::at(center, i);
  }
  one_over_r =
      1.0 / sqrt(square(get<0>(n)) + square(get<1>(n)) + square(get<2>(n)));
  for (size_t i = 0; i < 3; ++i) {
    n.get(i) *= one_over_r;
  }
  phi = 2.0 * mass * one_over_r;
  one_over_one_plus_phi = 1.0 / (1.0 + phi);

  DataVector& lapse_dv = get(*lapse);
  lapse_dv = sqrt(one_over_one_plus_phi);
  // lapse^3 = lapse / (1 + phi)
  DataVector& trace_k_dv = get(*trace_extrinsic_curvature);
  trace_k_dv = lapse_dv * one_over_one_plus_phi;
  // deriv_lapse_i = phi/(2r) lapse^3 n_i
  factor = 0.5 * phi * one_over_r * trace_k_dv;
  for (size_t i = 0; i < 3; ++i) {
    deriv_lapse->get(i) = factor * n.get(i);
  }
  // K = phi/r (1 + 3 phi/2) lapse^3, reusing lapse^3 stored in trace_k_dv
  trace_k_dv *= phi * one_over_r * (1.0 + 1.5 * phi);
  // trace_christoffel^i = phi(3+4phi)/(2r(1+phi)^2) n^i
  factor = 0.5 * phi * one_over_r * (3.0 + 4.0 * phi) *
           square(one_over_one_plus_phi);
  for (size_t i = 0; i < 3; ++i) {
    trace_spatial_christoffel_second_kind->get(i) = factor * n.get(i);
  }
  // shift^i = phi/(1+phi) n^i
  factor = phi * one_over_one_plus_phi;
  for (size_t i = 0; i < 3; ++i) {
    shift->get(i) = factor * n.get(i);
  }
  // gamma_ij = delta_ij + phi n_i n_j and
  // gamma^ij = delta^ij - phi/(1+phi) n^i n^j, reusing the shift factor
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      if (i == j) {
        spatial_metric->get(i, j) = 1.0 + phi * square(n.get(i));
        inverse_spatial_metric->get(i, j) = 1.0 - factor * square(n.get(i));
      } else {
        spatial_metric->get(i, j) = phi * n.get(i) * n.get(j);
        inverse_spatial_metric->get(i, j) = -factor * n.get(i) * n.get(j);
      }
    }
  }
  // deriv_shift_i^j = phi/(r(1+phi)) delta_i^j
  //                   - phi(2+phi)/(r(1+phi)^2) n_i n^j
  factor = phi * one_over_r * one_over_one_plus_phi;
  // c2 stored in one_over_r, which is no longer needed
  DataVector& c2 = one_over_r;
  c2 = factor * (2.0 + phi) * one_over_one_plus_phi;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      if (i == j) {
        deriv_shift->get(i, j) = factor - c2 * square(n.get(i));
      } else {
        deriv_shift->get(i, j) = -c2 * n.get(i) * n.get(j);
      }
    }
  }
}

}  // namespace CurvedScalarWave
