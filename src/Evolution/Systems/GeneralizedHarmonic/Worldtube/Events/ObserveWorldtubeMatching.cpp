// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Events/ObserveWorldtubeMatching.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <optional>
#include <pup.h>
#include <string>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/SliceTensorToVariables.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FaceNormal.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gh::worldtube::Events {
namespace {
struct ElectricTag : db::SimpleTag {
  using type = tnsr::ii<DataVector, 3, Frame::Inertial>;
};
struct MagneticTag : db::SimpleTag {
  using type = tnsr::ii<DataVector, 3, Frame::Inertial>;
};

double frobenius_max(const tnsr::ii<DataVector, 3, Frame::Inertial>& tensor,
                     const tnsr::II<DataVector, 3, Frame::Inertial>& inverse) {
  DataVector square_norm(get_size(get<0, 0>(tensor)), 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          square_norm += inverse.get(i, k) * inverse.get(j, l) *
                         tensor.get(i, j) * tensor.get(k, l);
        }
      }
    }
  }
  return sqrt(max(square_norm));
}

double max_abs(const ComplexDataVector& values) {
  double result = 0.;
  for (const auto& value : values) {
    result = std::max(result, std::abs(value));
  }
  return result;
}

// The largest pointwise relative misfit of the type-D solve over the five
// scalars
double type_d_residual(const gr::np::WeylScalars& predicted,
                       const gr::np::WeylScalars& psi) {
  const size_t num_points = get_size(psi.get(0));
  double result = 0.;
  for (size_t p = 0; p < num_points; ++p) {
    double misfit = 0.;
    double norm = 0.;
    for (size_t a = 0; a < 5; ++a) {
      misfit += std::norm(predicted.get(a)[p] - psi.get(a)[p]);
      norm += std::norm(psi.get(a)[p]);
    }
    result = std::max(result, sqrt(misfit / std::max(norm, 1.e-300)));
  }
  return result;
}
}  // namespace

ObserveWorldtubeMatching::ObserveWorldtubeMatching(
    const std::string& subfile_name, const std::optional<double> mass)
    : subfile_path_("/" + subfile_name), mass_(mass) {}

std::vector<std::string> ObserveWorldtubeMatching::legend() {
  return {"Time",
          "NumberOfFaces",
          "MaxAbsPsi0",
          "MaxAbsElectric",
          "MaxAbsMagnetic",
          "MinReCoulomb",
          "MaxReCoulomb",
          "MaxAbsImCoulomb",
          "MaxTypeDResidual",
          "MaxAbsPsi0TypeD",
          "MaxTypeDMismatch",
          "MaxAbsPsi0Quadrupole",
          "MaxQuadrupoleMismatch",
          "Psi4FitRelativeResidual",
          "MaxAbsTanhRapidity",
          "MinRapidity",
          "MaxRapidity",
          "MinMeasuredRadius",
          "MaxMeasuredRadius"};
}

std::optional<MatchingReductionData>
ObserveWorldtubeMatching::compute_reduction_data(
    const double time,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Element<3>& element, const Domain<3>& domain,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& mesh_velocity,
    const KretschmannFaceData<3>& face_data) const {
  const auto direction =
      excision_face_direction(domain.excision_spheres(), element);
  if (not direction.has_value()) {
    return std::nullopt;
  }
  const size_t sliced_dim = direction->dimension();
  const size_t fixed_index = index_to_slice_at(mesh.extents(), *direction);
  const Mesh<2> face_mesh = mesh.slice_away(sliced_dim);

  // Curvature in the volume, sliced to the face
  const auto d_pi = partial_derivative(pi, mesh, inverse_jacobian);
  const auto d_phi = partial_derivative(phi, mesh, inverse_jacobian);
  const WeylCurvature curvature =
      weyl_curvature(spacetime_metric, pi, phi, d_pi, d_phi);
  const auto face = data_on_slice<
      gr::Tags::SpatialMetric<DataVector, 3>,
      gr::Tags::InverseSpatialMetric<DataVector, 3>,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, 3>, ElectricTag,
      MagneticTag,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>>(
      mesh.extents(), sliced_dim, fixed_index, curvature.spatial_metric,
      curvature.inverse_spatial_metric, curvature.lapse, curvature.shift,
      curvature.electric, curvature.magnetic, inverse_jacobian);
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataVector, 3>>(face);
  const auto& inverse_spatial_metric =
      get<gr::Tags::InverseSpatialMetric<DataVector, 3>>(face);
  const auto& electric = get<ElectricTag>(face);
  const auto& magnetic = get<MagneticTag>(face);

  // Outward unit normal covector of the face
  auto normal_covector = unnormalized_face_normal(
      face_mesh,
      get<domain::Tags::InverseJacobian<3, Frame::ElementLogical,
                                        Frame::Inertial>>(face),
      *direction);
  const DataVector normal_magnitude =
      get(magnitude(normal_covector, inverse_spatial_metric));
  for (size_t i = 0; i < 3; ++i) {
    normal_covector.get(i) /= normal_magnitude;
  }

  const MatchingEvaluation type_d = evaluate_matching(
      PhysicalModel::TypeD, std::nullopt, electric, magnetic, spatial_metric,
      normal_covector, get<gr::Tags::Lapse<DataVector>>(face),
      get<gr::Tags::Shift<DataVector, 3>>(face), nullptr);

  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  double max_abs_psi0_quadrupole = nan;
  double max_quadrupole_mismatch = nan;
  double fit_residual = nan;
  double max_abs_tanh_rapidity = nan;
  double min_rapidity = nan;
  double max_rapidity = nan;
  double min_radius = nan;
  double max_radius = nan;
  if (mass_.has_value()) {
    // Face data of the observed state: the backward difference against the
    // stored history, at the observation time
    KretschmannFaceData<3> current_face_data = face_data;
    update_kretschmann_face_data(
        make_not_null(&current_face_data), spacetime_metric, pi, phi, mesh,
        inverse_jacobian, element, domain.excision_spheres(), mesh_velocity,
        time);
    // The order-two model errors on an unphysical boost; a diagnostic must
    // not, so check first
    const gr::np::FrameRegistration registration =
        gr::np::register_frame(type_d.psi, type_d.adapted_rotation, *mass_);
    const auto& lapse = get<gr::Tags::Lapse<DataVector>>(face);
    const auto& shift = get<gr::Tags::Shift<DataVector, 3>>(face);
    max_abs_tanh_rapidity = max(abs(get(gr::np::invariant_tanh_rapidity(
        registration.member, type_d.adapted_rotation, spatial_metric, lapse,
        shift, current_face_data.d_kretschmann,
        current_face_data.dt_kretschmann))));
    min_radius = min(get(registration.measured_radius));
    max_radius = max(get(registration.measured_radius));
    if (max_abs_tanh_rapidity < 1.) {
      const MatchingEvaluation quadrupole = evaluate_matching(
          PhysicalModel::Quadrupole, mass_, electric, magnetic, spatial_metric,
          normal_covector, lapse, shift, &current_face_data);
      max_abs_psi0_quadrupole = max_abs(get(quadrupole.psi0_target));
      max_quadrupole_mismatch =
          max_abs(get(quadrupole.psi0_target) - quadrupole.psi.get(0));
      fit_residual = quadrupole.second_order->fit.relative_residual;
      min_rapidity = min(get(*quadrupole.rapidity));
      max_rapidity = max(get(*quadrupole.rapidity));
    }
  }

  return MatchingReductionData{
      time,
      1_st,
      max_abs(type_d.psi.get(0)),
      frobenius_max(electric, inverse_spatial_metric),
      frobenius_max(magnetic, inverse_spatial_metric),
      min(real(get(type_d.coulomb))),
      max(real(get(type_d.coulomb))),
      max(abs(imag(get(type_d.coulomb)))),
      type_d_residual(type_d.type_d_rotation.predicted_psi, type_d.psi),
      max_abs(get(type_d.psi0_target)),
      max_abs(get(type_d.psi0_target) - type_d.psi.get(0)),
      max_abs_psi0_quadrupole,
      max_quadrupole_mismatch,
      fit_residual,
      max_abs_tanh_rapidity,
      min_rapidity,
      max_rapidity,
      min_radius,
      max_radius};
}

void ObserveWorldtubeMatching::pup(PUP::er& p) {
  Event::pup(p);
  p | subfile_path_;
  p | mass_;
}

PUP::able::PUP_ID ObserveWorldtubeMatching::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::worldtube::Events
