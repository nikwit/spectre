// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <random>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/DataStructures/RandomUnitNormal.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylScalars.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.GeneralRelativity.WeylScalars.Schwarzschild",
    "[Unit][PointwiseFunctions]") {
  const double mass = 1.98765;
  const size_t num_points = 1000;
  MAKE_GENERATOR(generator);
  const std::uniform_real_distribution<> coord_dist(4.0, 10.);
  const auto coords =
      make_with_random_values<tnsr::I<DataVector, 3, Frame::Inertial>>(
          make_not_null(&generator), coord_dist, num_points);

  const auto radius = magnitude(coords);
  const auto r = get(radius);
  const auto theta =
      atan2(hypot(get<0>(coords), get<1>(coords)), get<2>(coords));

  tnsr::ii<DataVector, 3, Frame::Inertial> spatial_metric(num_points, 0.0);

  get<0, 0>(spatial_metric) = 1. / (1 - 2. * mass / r);
  get<1, 1>(spatial_metric) = r * r;
  get<2, 2>(spatial_metric) = r * r * sin(theta) * sin(theta);

  tnsr::Ij<DataVector, 3, Frame::Inertial> weyl_electric_updown(num_points,
                                                                0.0);
  get<0, 0>(weyl_electric_updown) = -(2.0 * mass) / (r * r * r);
  get<1, 1>(weyl_electric_updown) = mass / (r * r * r);
  get<2, 2>(weyl_electric_updown) = mass / (r * r * r);

  tnsr::ii<DataVector, 3, Frame::Inertial> weyl_electric(num_points, 0.0);
  tenex::evaluate<ti::i, ti::j>(
      make_not_null(&weyl_electric),
      weyl_electric_updown(ti::K, ti::j) * spatial_metric(ti::i, ti::k));

  const auto weyl_magnetic =
      make_with_value<tnsr::ii<DataVector, 3, Frame::Inertial>>(coords, 0.0);

  const auto spatial_normal_vector =
      random_unit_normal(make_not_null(&generator), spatial_metric);

  const auto weyl_scalars = gr::weyl_scalars(
      weyl_electric, weyl_magnetic, spatial_metric, spatial_normal_vector);

  const auto& psi_0 = weyl_scalars.scalars[0];
  const auto& psi_1 = weyl_scalars.scalars[1];
  const auto& psi_2 = weyl_scalars.scalars[2];
  const auto& psi_3 = weyl_scalars.scalars[3];
  const auto& psi_4 = weyl_scalars.scalars[4];

  Scalar<ComplexDataVector> invariant_I{ComplexDataVector(num_points, 0.0)};
  get(invariant_I) = get(psi_0) * get(psi_4) - 4.0 * get(psi_1) * get(psi_3) +
                     3.0 * get(psi_2) * get(psi_2);

  Scalar<ComplexDataVector> invariant_J{ComplexDataVector(num_points, 0.0)};
  get(invariant_J) = get(psi_0) * get(psi_2) * get(psi_4) +
                     2.0 * get(psi_1) * get(psi_2) * get(psi_3) -
                     get(psi_2) * get(psi_2) * get(psi_2) -
                     get(psi_0) * get(psi_3) * get(psi_3) -
                     get(psi_4) * get(psi_1) * get(psi_1);

  const DataVector m_over_r_cubed = mass / cube(r);
  const DataVector expected_I_real = 3.0 * square(m_over_r_cubed);
  const DataVector expected_J_real = -cube(m_over_r_cubed);
  ComplexDataVector expected_I(num_points, 0.0);
  ComplexDataVector expected_J(num_points, 0.0);
  for (size_t i = 0; i < num_points; ++i) {
    expected_I[i] = std::complex<double>(expected_I_real[i], 0.0);
    expected_J[i] = std::complex<double>(expected_J_real[i], 0.0);
  }

  Approx local_approx = Approx::custom().epsilon(1e-14).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(get(invariant_I), expected_I, local_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(get(invariant_J), expected_J, local_approx);
}
