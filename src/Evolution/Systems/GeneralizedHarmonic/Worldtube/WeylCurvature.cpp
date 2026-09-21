// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/SliceTensorToVariables.hpp"
#include "DataStructures/Tensor/EagerMath/Determinant.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/FaceNormal.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/CovariantDerivOfExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::worldtube {
void weyl_electric_magnetic(
    const gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*> electric,
    const gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*> magnetic,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::A<DataVector, 3, Frame::Inertial>& spacetime_unit_normal_vector,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& extrinsic_curvature,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric) {
  const size_t num_points = get_size(get<0, 0>(spatial_metric));
  // d_k gamma_ij = Phi_kij
  tnsr::ijj<DataVector, 3, Frame::Inertial> d_spatial_metric(num_points);
  for (size_t k = 0; k < 3; ++k) {
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = i; j < 3; ++j) {
        d_spatial_metric.get(k, i, j) = phi.get(k, i + 1, j + 1);
      }
    }
  }
  const auto christoffel_second_kind = raise_or_lower_first_index(
      gr::christoffel_first_kind(d_spatial_metric), inverse_spatial_metric);
  const auto cov_deriv_extrinsic_curvature =
      gh::covariant_deriv_of_extrinsic_curvature(
          extrinsic_curvature, spacetime_unit_normal_vector,
          christoffel_second_kind, inverse_spacetime_metric, phi, d_pi, d_phi);
  const auto ricci =
      gh::spatial_ricci_tensor(phi, d_phi, inverse_spatial_metric);
  gr::weyl_electric(electric, ricci, extrinsic_curvature,
                    inverse_spatial_metric);
  const Scalar<DataVector> sqrt_det_spatial_metric{
      sqrt(get(determinant(spatial_metric)))};
  gr::weyl_magnetic(magnetic, cov_deriv_extrinsic_curvature, spatial_metric,
                    sqrt_det_spatial_metric);
}

WeylCurvature weyl_curvature(
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& d_pi,
    const tnsr::ijaa<DataVector, 3, Frame::Inertial>& d_phi) {
  const size_t num_points = get_size(get<0, 0>(spacetime_metric));
  WeylCurvature result{};
  result.spatial_metric =
      tnsr::ii<DataVector, 3, Frame::Inertial>(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      result.spatial_metric.get(i, j) = spacetime_metric.get(i + 1, j + 1);
    }
  }
  result.inverse_spatial_metric =
      determinant_and_inverse(result.spatial_metric).second;
  result.shift = gr::shift(spacetime_metric, result.inverse_spatial_metric);
  result.lapse = gr::lapse(result.shift, spacetime_metric);
  const auto spacetime_unit_normal_vector =
      gr::spacetime_normal_vector(result.lapse, result.shift);
  const auto inverse_spacetime_metric = gr::inverse_spacetime_metric(
      result.lapse, result.shift, result.inverse_spatial_metric);
  const auto extrinsic_curvature =
      gh::extrinsic_curvature(spacetime_unit_normal_vector, pi, phi);
  weyl_electric_magnetic(make_not_null(&result.electric),
                         make_not_null(&result.magnetic), phi, d_phi, d_pi,
                         spacetime_unit_normal_vector, result.spatial_metric,
                         result.inverse_spatial_metric, extrinsic_curvature,
                         inverse_spacetime_metric);
  return result;
}

namespace {
struct ElectricTag : db::SimpleTag {
  using type = tnsr::ii<DataVector, 3, Frame::Inertial>;
};
struct MagneticTag : db::SimpleTag {
  using type = tnsr::ii<DataVector, 3, Frame::Inertial>;
};
}  // namespace

FaceCurvature face_curvature(
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Direction<3>& direction) {
  const size_t sliced_dim = direction.dimension();
  const size_t fixed_index = index_to_slice_at(mesh.extents(), direction);
  const auto d_pi = partial_derivative(pi, mesh, inverse_jacobian);
  const auto d_phi = partial_derivative(phi, mesh, inverse_jacobian);
  const WeylCurvature volume =
      weyl_curvature(spacetime_metric, pi, phi, d_pi, d_phi);
  auto face = data_on_slice<
      gr::Tags::SpatialMetric<DataVector, 3>,
      gr::Tags::InverseSpatialMetric<DataVector, 3>,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>, ElectricTag,
      MagneticTag,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>>(
      mesh.extents(), sliced_dim, fixed_index, volume.spatial_metric,
      volume.inverse_spatial_metric, volume.lapse, volume.shift,
      volume.electric, volume.magnetic, inverse_jacobian);
  FaceCurvature result{};
  result.spatial_metric = get<gr::Tags::SpatialMetric<DataVector, 3>>(face);
  result.inverse_spatial_metric =
      get<gr::Tags::InverseSpatialMetric<DataVector, 3>>(face);
  result.lapse = get<gr::Tags::Lapse<DataVector>>(face);
  result.shift = get<gr::Tags::Shift<DataVector, 3>>(face);
  result.electric = get<ElectricTag>(face);
  result.magnetic = get<MagneticTag>(face);
  result.unit_normal_covector = unnormalized_face_normal(
      mesh.slice_away(sliced_dim),
      get<domain::Tags::InverseJacobian<3, Frame::ElementLogical,
                                        Frame::Inertial>>(face),
      direction);
  const DataVector normal_magnitude = get(
      magnitude(result.unit_normal_covector, result.inverse_spatial_metric));
  for (size_t i = 0; i < 3; ++i) {
    result.unit_normal_covector.get(i) /= normal_magnitude;
  }
  return result;
}

Scalar<DataVector> kretschmann_scalar(
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inverse_spatial_metric) {
  Scalar<DataVector> result(get_size(get<0, 0>(electric)), 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          get(result) += 8. * inverse_spatial_metric.get(i, k) *
                         inverse_spatial_metric.get(j, l) *
                         (electric.get(i, j) * electric.get(k, l) -
                          magnetic.get(i, j) * magnetic.get(k, l));
        }
      }
    }
  }
  return result;
}
}  // namespace gh::worldtube
