// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeRecordingData.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace gh::worldtube {
void compute_recording_data(
    const gsl::not_null<Variables<recording_tags>*> result,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian) {
  result->initialize(mesh.number_of_grid_points());
  get<gr::Tags::SpacetimeMetric<DataVector, 3>>(*result) = metric;
  get<gh::Tags::Pi<DataVector, 3>>(*result) = pi;
  get<gh::Tags::Phi<DataVector, 3>>(*result) = phi;
  get<recording_deriv<gr::Tags::SpacetimeMetric<DataVector, 3>>>(*result) =
      partial_derivative(metric, mesh, inverse_jacobian);
  auto& d_pi = get<recording_deriv<gh::Tags::Pi<DataVector, 3>>>(*result);
  auto& d_phi = get<recording_deriv<gh::Tags::Phi<DataVector, 3>>>(*result);
  d_pi = partial_derivative(pi, mesh, inverse_jacobian);
  d_phi = partial_derivative(phi, mesh, inverse_jacobian);
  const auto curvature = weyl_curvature(metric, pi, phi, d_pi, d_phi);
  auto& kretschmann = get<Tags::Kretschmann>(*result);
  kretschmann = kretschmann_scalar(curvature.electric, curvature.magnetic,
                                   curvature.inverse_spatial_metric);
  get<recording_deriv<Tags::Kretschmann>>(*result) =
      partial_derivative(kretschmann, mesh, inverse_jacobian);
}
}  // namespace gh::worldtube
