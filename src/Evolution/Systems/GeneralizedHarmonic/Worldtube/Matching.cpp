// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"

#include <cstddef>
#include <optional>
#include <ostream>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Options/ParseError.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/NullRotations.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/RestFrame.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gh::worldtube {
PhysicalModel convert_physical_model_from_yaml(const Options::Option& options) {
  const auto read = options.parse_as<std::string>();
  if (read == "None") {
    return PhysicalModel::None;
  } else if (read == "TypeD") {
    return PhysicalModel::TypeD;
  } else if (read == "Quadrupole") {
    return PhysicalModel::Quadrupole;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to a physical model. Must be "
              "one of None, TypeD or Quadrupole.");
}

std::ostream& operator<<(std::ostream& os, const PhysicalModel model) {
  switch (model) {
    case PhysicalModel::None:
      return os << "None";
    case PhysicalModel::TypeD:
      return os << "TypeD";
    case PhysicalModel::Quadrupole:
      return os << "Quadrupole";
    default:
      ERROR("Unknown PhysicalModel");
  }
}

MatchingEvaluation evaluate_matching(
    const PhysicalModel model, const std::optional<double> mass,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial>& unit_normal_covector,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const KretschmannFaceData<3>* const face_data,
    const std::optional<gr::np::TidalMoments>& imposed_moments) {
  const size_t num_points = get_size(get<0>(unit_normal_covector));
  MatchingEvaluation result{};
  // The adapted triad takes the Euclidean direction of the normal covector;
  // s is then that covector normalized with the metric, i.e. the face normal
  // pointing out of the domain (see gr::np::adapted_triad()).
  tnsr::I<DataVector, 3, Frame::Inertial> directions(num_points, 0.);
  const DataVector euclidean_norm = sqrt(square(get<0>(unit_normal_covector)) +
                                         square(get<1>(unit_normal_covector)) +
                                         square(get<2>(unit_normal_covector)));
  for (size_t i = 0; i < 3; ++i) {
    directions.get(i) = unit_normal_covector.get(i) / euclidean_norm;
  }
  result.adapted_rotation = gr::np::adapted_triad(spatial_metric, directions);
  result.psi = gr::np::weyl_scalars_from_electric_magnetic(
      electric, magnetic, spatial_metric, directions);

  switch (model) {
    case PhysicalModel::TypeD: {
      result.coulomb = gr::np::coulomb_scalar(gr::np::invariant_i(result.psi),
                                              gr::np::invariant_j(result.psi));
      result.type_d_rotation =
          gr::np::solve_type_d_rotation(result.psi, result.coulomb);
      result.psi0_target =
          gr::np::psi0_leading(result.coulomb, result.type_d_rotation.b);
      break;
    }
    case PhysicalModel::Quadrupole: {
      if (not mass.has_value()) {
        ERROR("PhysicalModel: Quadrupole needs the mass of the hole.");
      }
      if (face_data == nullptr or not face_data->direction.has_value()) {
        ERROR(
            "PhysicalModel: Quadrupole needs the Kretschmann face data of the "
            "element, but none is available. The element must abut an "
            "excision sphere and gh::worldtube::UpdateKretschmannFaceData "
            "must run before the boundary condition is applied.");
      }
      if (get(face_data->kretschmann).size() != num_points or
          get(face_data->dt_kretschmann).size() != num_points) {
        ERROR("The Kretschmann face data has "
              << get(face_data->kretschmann).size()
              << " points but the face has " << num_points
              << ". The face data must be updated for the current mesh "
                 "before the boundary condition is applied.");
      }
      result.registration =
          gr::np::register_frame(result.psi, result.adapted_rotation, *mass);
      result.coulomb = result.registration->coulomb;
      result.type_d_rotation = result.registration->rotation;
      result.rapidity = gr::np::invariant_rapidity(
          result.registration->member, result.adapted_rotation, spatial_metric,
          lapse, shift, face_data->d_kretschmann, face_data->dt_kretschmann);
      const std::optional<DataVector> weights =
          face_data->quadrature_weights.size() == num_points
              ? std::optional<DataVector>{face_data->quadrature_weights}
              : std::nullopt;
      result.second_order = gr::np::evaluate_second_order(
          *result.registration, *result.rapidity, result.adapted_rotation,
          *mass, weights, imposed_moments);
      result.psi0_target = result.second_order->psi0_target;
      break;
    }
    case PhysicalModel::None:
      ERROR("PhysicalModel: None supplies no incoming mode to evaluate.");
    default:
      ERROR("Unknown PhysicalModel");
  }

  // U^{8-} = w^- / 2 in covariant coordinate components
  result.incoming_mode = gr::np::orthonormal_to_coordinate_covariant(
      gr::np::incoming_weyl_field(result.psi0_target, result.adapted_rotation),
      gr::np::cholesky_factor(spatial_metric));
  for (auto& component : result.incoming_mode) {
    component *= 0.5;
  }
  return result;
}
}  // namespace gh::worldtube
