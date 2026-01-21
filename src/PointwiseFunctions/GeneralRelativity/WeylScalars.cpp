// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/WeylScalars.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>

#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/ProjectionOperators.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylPropagating.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace gr {
template <typename Frame>
void weyl_scalars(
    gsl::not_null<WeylScalarsResult<Frame>*> weyl_scalars_result,
                  const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
                  const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
                  const tnsr::ii<DataVector, 3, Frame>& spatial_metric,
                  const tnsr::I<DataVector, 3, Frame>& spatial_normal_vector) {
#ifdef SPECTRE_DEBUG
  const auto magnitude_normal =
      magnitude(spatial_normal_vector, spatial_metric);
  for (size_t i = 0; i < get<0>(spatial_normal_vector).size(); i++) {
    ASSERT(std::abs(magnitude_normal.get()[i] - 1.0) < 1e-12,
           "The spatial normal vector must be normalized. Its magnitude is "
               << magnitude_normal.get()[i] << " at index " << i);
  }
#endif  // SPECTRE_DEBUG
  const size_t num_points = get<0>(spatial_normal_vector).size();
  tnsr::I<DataVector, 3, Frame> vector1(num_points, 0.0);
  tnsr::I<DataVector, 3, Frame> vector2(num_points, 0.0);
  for (size_t i = 0; i < num_points; i++) {
    const auto& x = get<0>(spatial_normal_vector)[i];
    const auto& y = get<1>(spatial_normal_vector)[i];
    const auto& z = get<2>(spatial_normal_vector)[i];
    bool y_greater_x = std::abs(y) > std::abs(x);
    vector1.get(0)[i] = y_greater_x ? 1.0 : 0.0;
    vector1.get(1)[i] = y_greater_x ? 0.0 : 1.0;
    if (std::abs(z) > (y_greater_x ? std::abs(y) : std::abs(x))) {
      const size_t index_not_set = y_greater_x ? 1 : 0;
      vector2.get(index_not_set)[i] = 1.0;
    } else {
      vector2.get(2)[i] = 1.0;
    }
  }

  // Modified Gram-Schmidt to build an orthonormal tetrad.
  const auto& e0 = spatial_normal_vector;
  tnsr::I<DataVector, 3, Frame> e1(num_points, 0.0);
  tnsr::I<DataVector, 3, Frame> e2(num_points, 0.0);
  Scalar<DataVector> dot01{DataVector(num_points, 0.0)};
  Scalar<DataVector> dot02{DataVector(num_points, 0.0)};
  Scalar<DataVector> dot12{DataVector(num_points, 0.0)};
  Scalar<DataVector> mag1{DataVector(num_points, 0.0)};
  Scalar<DataVector> mag2{DataVector(num_points, 0.0)};

  dot_product(make_not_null(&dot01), e0, vector1, spatial_metric);
  tenex::evaluate<ti::I>(make_not_null(&e1),
                         vector1(ti::I) - dot01() * e0(ti::I));
  magnitude(make_not_null(&mag1), e1, spatial_metric);
  tenex::evaluate<ti::I>(make_not_null(&e1), e1(ti::I) / mag1());

  dot_product(make_not_null(&dot02), e0, vector2, spatial_metric);
  dot_product(make_not_null(&dot12), e1, vector2, spatial_metric);
  tenex::evaluate<ti::I>(
      make_not_null(&e2),
      vector2(ti::I) - dot02() * e0(ti::I) - dot12() * e1(ti::I));
  magnitude(make_not_null(&mag2), e2, spatial_metric);
  tenex::evaluate<ti::I>(make_not_null(&e2), e2(ti::I) / mag2());

#ifdef SPECTRE_DEBUG
  const auto mag_e0 = magnitude(e0, spatial_metric);
  const auto mag_e1 = magnitude(e1, spatial_metric);
  const auto mag_e2 = magnitude(e2, spatial_metric);
  const auto dot_e0e1 = dot_product(e0, e1, spatial_metric);
  const auto dot_e0e2 = dot_product(e0, e2, spatial_metric);
  const auto dot_e1e2 = dot_product(e1, e2, spatial_metric);
  const double tol = 1e-8;
  for (size_t i = 0; i < get(mag_e0).size(); i++) {
    ASSERT(std::abs(get(mag_e0)[i] - 1.) < tol,
           "The spatial normal vector must be normalized. Its magnitude is "
               << get(mag_e0)[i] << " at index " << i);
    ASSERT(std::abs(get(mag_e1)[i] - 1.) < tol,
           "The spatial normal vector must be normalized. Its magnitude is "
               << get(mag_e1)[i] << " at index " << i);
    ASSERT(std::abs(get(mag_e2)[i] - 1.) < tol,
           "The spatial normal vector must be normalized. Its magnitude is "
               << get(mag_e2)[i] << " at index " << i);

    ASSERT(std::abs(get(dot_e0e1)[i]) < tol,
           "e0 and e1 must be orthogonal. Their dot product is "
               << get(dot_e0e1)[i] << " at index " << i);
    ASSERT(std::abs(get(dot_e0e2)[i]) < tol,
           "e0 and e2 must be orthogonal. Their dot product is "
               << get(dot_e0e2)[i] << " at index " << i);
    ASSERT(std::abs(get(dot_e1e2)[i]) < tol,
           "e1 and e2 must be orthogonal. Their dot product is "
               << get(dot_e1e2)[i] << " at index " << i);
  }
#endif  // SPECTRE_DEBUG

  const std::complex<double> imag = std::complex<double>(0.0, 1.0);
  const auto electric_minus_magnetic = tenex::evaluate<ti::i, ti::j>(
      weyl_electric(ti::i, ti::j) - imag * weyl_magnetic(ti::i, ti::j));

  const auto m =
      tenex::evaluate<ti::I>(M_SQRT1_2 * (e1(ti::I) + imag * e2(ti::I)));
  const auto mbar =
      tenex::evaluate<ti::I>(M_SQRT1_2 * (e1(ti::I) - imag * e2(ti::I)));
  weyl_scalars_result->m = m;

  auto& psi_0 = weyl_scalars_result->scalars[0];
  psi_0 = tenex::evaluate(-electric_minus_magnetic(ti::i, ti::j) * m(ti::I) *
                          m(ti::J));
  auto& psi_1 = weyl_scalars_result->scalars[1];
  psi_1 = tenex::evaluate(M_SQRT1_2 * electric_minus_magnetic(ti::i, ti::j) *
                          m(ti::I) * e0(ti::J));
  auto& psi_2 = weyl_scalars_result->scalars[2];
  psi_2 = tenex::evaluate(-0.5 * electric_minus_magnetic(ti::i, ti::j) *
                          e0(ti::I) * e0(ti::J));
  auto& psi_3 = weyl_scalars_result->scalars[3];
  psi_3 = tenex::evaluate(-M_SQRT1_2 * electric_minus_magnetic(ti::i, ti::j) *
                          mbar(ti::I) * e0(ti::J));
  auto& psi_4 = weyl_scalars_result->scalars[4];
  psi_4 = tenex::evaluate(-electric_minus_magnetic(ti::i, ti::j) * mbar(ti::I) *
                          mbar(ti::J));
}

template <typename Frame>
WeylScalarsResult<Frame> weyl_scalars(
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3, Frame>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame>& spatial_normal_vector) {
  WeylScalarsResult<Frame> weyl_scalars_result{
      {}, tnsr::I<ComplexDataVector, 3, Frame>{}};
  for (auto& scalar : weyl_scalars_result.scalars) {
    scalar = make_with_value<Scalar<ComplexDataVector>>(
        get<0, 0>(weyl_electric), std::numeric_limits<double>::signaling_NaN());
  }
  weyl_scalars_result.m =
      make_with_value<tnsr::I<ComplexDataVector, 3, Frame>>(
          get<0, 0>(weyl_electric),
          std::numeric_limits<double>::signaling_NaN());
  weyl_scalars(make_not_null(&weyl_scalars_result), weyl_electric,
               weyl_magnetic, spatial_metric, spatial_normal_vector);
  return weyl_scalars_result;
}
}  // namespace gr

#define FRAME(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                             \
  template void gr::weyl_scalars(                                        \
      const gsl::not_null<gr::WeylScalarsResult<FRAME(data)>*>           \
          weyl_scalars_result,                                           \
      const tnsr::ii<DataVector, 3, FRAME(data)>& weyl_electric,         \
      const tnsr::ii<DataVector, 3, FRAME(data)>& weyl_magnetic,         \
      const tnsr::ii<DataVector, 3, FRAME(data)>& spatial_metric,        \
      const tnsr::I<DataVector, 3, FRAME(data)>& spatial_normal_vector); \
  template gr::WeylScalarsResult<FRAME(data)> gr::weyl_scalars(          \
      const tnsr::ii<DataVector, 3, FRAME(data)>& weyl_electric,         \
      const tnsr::ii<DataVector, 3, FRAME(data)>& weyl_magnetic,         \
      const tnsr::ii<DataVector, 3, FRAME(data)>& spatial_metric,        \
      const tnsr::I<DataVector, 3, FRAME(data)>& spatial_normal_vector);

GENERATE_INSTANTIATIONS(INSTANTIATE, (Frame::Grid, Frame::Inertial))

#undef FRAME
#undef INSTANTIATE
