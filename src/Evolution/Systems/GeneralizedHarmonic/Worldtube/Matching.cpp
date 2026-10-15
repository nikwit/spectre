// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <ostream>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Options/ParseError.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/CoulombDecode.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/GeometricTide.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/LongitudinalTide.hpp"
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
  } else if (read == "QuadrupoleLongitudinal") {
    return PhysicalModel::QuadrupoleLongitudinal;
  } else if (read == "QuadrupoleGeometric") {
    return PhysicalModel::QuadrupoleGeometric;
  } else if (read == "ThirdOrderGeometric") {
    return PhysicalModel::ThirdOrderGeometric;
  } else if (read == "QuadrupoleCoulomb") {
    return PhysicalModel::QuadrupoleCoulomb;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to a physical model. Must be "
              "one of None, TypeD, Quadrupole, QuadrupoleGeometric, "
              "ThirdOrderGeometric, QuadrupoleLongitudinal or "
              "QuadrupoleCoulomb.");
}

std::ostream& operator<<(std::ostream& os, const PhysicalModel model) {
  switch (model) {
    case PhysicalModel::None:
      return os << "None";
    case PhysicalModel::TypeD:
      return os << "TypeD";
    case PhysicalModel::Quadrupole:
      return os << "Quadrupole";
    case PhysicalModel::QuadrupoleLongitudinal:
      return os << "QuadrupoleLongitudinal";
    case PhysicalModel::QuadrupoleGeometric:
      return os << "QuadrupoleGeometric";
    case PhysicalModel::ThirdOrderGeometric:
      return os << "ThirdOrderGeometric";
    case PhysicalModel::QuadrupoleCoulomb:
      return os << "QuadrupoleCoulomb";
    default:
      ERROR("Unknown PhysicalModel");
  }
}

bool is_order_two(const PhysicalModel model) {
  return model == PhysicalModel::QuadrupoleLongitudinal or
         model == PhysicalModel::QuadrupoleGeometric or
         model == PhysicalModel::Quadrupole or
         model == PhysicalModel::QuadrupoleCoulomb;
}

Scalar<DataVector> normal_derivative_of_coulomb(
    const Scalar<ComplexDataVector>& coulomb,
    const tnsr::i<DataVector, 3, Frame::Inertial>& d_kretschmann,
    const tnsr::I<DataVector, 3, Frame::Inertial>& unit_normal_vector) {
  Scalar<DataVector> result(get_size(get(coulomb)), 0.);
  for (size_t i = 0; i < 3; ++i) {
    get(result) += unit_normal_vector.get(i) * d_kretschmann.get(i);
  }
  get(result) /= 96. * real(get(coulomb));
  return result;
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
    const std::optional<gr::np::TidalMoments>& imposed_moments,
    const bool use_fixed_point) {
  return evaluate_matching_with_fit_input(
      model, mass, electric, magnetic, spatial_metric, unit_normal_covector,
      lapse, shift, face_data, imposed_moments, use_fixed_point, nullptr);
}

MatchingEvaluation evaluate_matching_with_fit_input(
    const PhysicalModel model, const std::optional<double> mass,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& electric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& magnetic,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::i<DataVector, 3, Frame::Inertial>& unit_normal_covector,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const KretschmannFaceData<3>* const face_data,
    const std::optional<gr::np::TidalMoments>& imposed_moments,
    const bool use_fixed_point, const gr::np::WeylScalars* const fit_psi) {
  if (use_fixed_point and
      (model != PhysicalModel::QuadrupoleGeometric or imposed_moments)) {
    ERROR(
        "Geometric Psi0 fixed point requires QuadrupoleGeometric without "
        "imposed moments");
  }
  if (fit_psi != nullptr and (model != PhysicalModel::QuadrupoleGeometric or
                              imposed_moments or use_fixed_point)) {
    ERROR("Independent fit input requires unrelaxed QuadrupoleGeometric.");
  }
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
    case PhysicalModel::QuadrupoleLongitudinal: {
      if (not mass or imposed_moments or face_data == nullptr or
          not face_data->direction or
          face_data->quadrature_weights.size() != num_points or
          get(face_data->dt_kretschmann).size() != num_points) {
        ERROR(
            "QuadrupoleLongitudinal requires mass, current full-sphere face "
            "data and no imposed moments");
      }
      for (auto& component : directions) {
        component *= -1.;
      }
      const size_t grid_lmax = static_cast<size_t>(
          (sqrt(8. * static_cast<double>(num_points) + 1.) - 3.) / 4. + .5);
      const char* setting = std::getenv("NP_GEOMETRIC_LMAX");
      const size_t lm = setting == nullptr
                            ? std::min(size_t{8}, grid_lmax)
                            : static_cast<size_t>(std::stoul(setting));
      if ((grid_lmax + 1) * (2 * grid_lmax + 1) != num_points or lm < 2 or
          lm > grid_lmax) {
        ERROR(
            "Longitudinal matching requires a complete spherical face and "
            "NP_GEOMETRIC_LMAX between 2 and grid l_max");
      }
      auto longitudinal = gr::np::evaluate_longitudinal_tide(
          result.psi, result.adapted_rotation, spatial_metric, directions,
          face_data->quadrature_weights, *mass, lm, lapse, shift,
          face_data->d_kretschmann, face_data->dt_kretschmann);
      result.psi0_target =
          Scalar<ComplexDataVector>{longitudinal.completed.get(0)};
      result.coulomb = longitudinal.registration.coulomb;
      result.type_d_rotation = longitudinal.registration.rotation;
      result.registration = std::move(longitudinal.registration);
      result.rapidity = std::move(longitudinal.rapidity);
      // The longitudinal fit has ten real columns, not five complex columns.
      // Expose its moments and residual, but not the legacy direct columns.
      result.second_order.emplace();
      result.second_order->fit.components = longitudinal.moments;
      result.second_order->fit.relative_residual =
          longitudinal.relative_residual;
      result.second_order->psi0_target = result.psi0_target;
      break;
    }
    case PhysicalModel::ThirdOrderGeometric: {
      if (face_data == nullptr or not face_data->third_order.has_value() or
          face_data->third_order->time != face_data->time or
          get(face_data->third_order->psi0_target).size() != num_points) {
        ERROR(
            "ThirdOrderGeometric needs its current target from "
            "UpdateKretschmannFaceData before the boundary condition");
      }
      result.psi0_target = face_data->third_order->psi0_target;
      break;
    }
    case PhysicalModel::TypeD: {
      result.coulomb = gr::np::coulomb_scalar(gr::np::invariant_i(result.psi),
                                              gr::np::invariant_j(result.psi));
      result.type_d_rotation =
          gr::np::solve_type_d_rotation(result.psi, result.coulomb);
      result.psi0_target =
          gr::np::psi0_leading(result.coulomb, result.type_d_rotation.b);
      break;
    }
    case PhysicalModel::QuadrupoleGeometric:
      if (use_fixed_point) {
        if (face_data == nullptr or not face_data->direction or
            not face_data->psi0_fixed_point or
            face_data->psi0_fixed_point->time != face_data->time or
            get(face_data->psi0_fixed_point->psi0_target).size() !=
                num_points) {
          ERROR(
              "Geometric Psi0 fixed point needs its current converged target "
              "from UpdateKretschmannFaceData");
        }
        result.psi0_target = face_data->psi0_fixed_point->psi0_target;
        break;
      }
      [[fallthrough]];
    case PhysicalModel::Quadrupole:
    case PhysicalModel::QuadrupoleCoulomb: {
      if (not mass.has_value()) {
        ERROR("PhysicalModel: " << model << " needs the mass of the hole.");
      }
      if (face_data == nullptr or not face_data->direction.has_value()) {
        ERROR("PhysicalModel: "
              << model
              << " needs the Kretschmann face data of the element, but none "
                 "is available. The element must abut an excision sphere and "
                 "gh::worldtube::UpdateKretschmannFaceData must run before the "
                 "boundary condition is applied.");
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
      std::optional<gr::np::TidalMoments> moments = imposed_moments;
      if (model == PhysicalModel::QuadrupoleCoulomb and
          not moments.has_value()) {
        // The tide from the Coulomb channel: unit normal vector s^i =
        // gamma^{ij} n_j for the normal derivative of K
        const auto inverse_spatial_metric =
            determinant_and_inverse(spatial_metric).second;
        tnsr::I<DataVector, 3, Frame::Inertial> unit_normal_vector(num_points,
                                                                   0.);
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            unit_normal_vector.get(i) +=
                inverse_spatial_metric.get(i, j) * unit_normal_covector.get(j);
          }
        }
        result.coulomb_decode = gr::np::decode_tidal_moments_from_coulomb(
            *result.registration, *result.rapidity, *mass,
            normal_derivative_of_coulomb(
                result.coulomb, face_data->d_kretschmann, unit_normal_vector),
            weights);
        if (not result.coulomb_decode->valid) {
          ERROR(
              "The Coulomb decode of the tide failed: the radius solve did "
              "not stay on the outer branch of sqrt(1 - 2M/r)/r^4 or did not "
              "converge (largest relative residual "
              << result.coulomb_decode->newton_residual
              << "). The excision must lie outside about 2.4M for "
                 "PhysicalModel: QuadrupoleCoulomb.");
        }
        moments = result.coulomb_decode->components;
      }
      if (model == PhysicalModel::QuadrupoleGeometric) {
        if (not weights.has_value()) {
          ERROR("Geometric matching requires full-sphere quadrature weights");
        }
        // On a round coordinate sphere the covector direction gives the
        // NR angular label. The inner domain normal points towards the hole.
        for (auto& component : directions) {
          component *= -1.;
        }
        const size_t grid_lmax = static_cast<size_t>(
            (sqrt(8. * static_cast<double>(num_points) + 1.) - 3.) / 4. + 0.5);
        if ((grid_lmax + 1) * (2 * grid_lmax + 1) != num_points) {
          ERROR("Geometric matching requires one complete spherical face");
        }
        const char* map_setting = std::getenv("NP_GEOMETRIC_LMAX");
        const size_t map_lmax =
            map_setting == nullptr
                ? std::min(size_t{8}, grid_lmax)
                : static_cast<size_t>(std::stoul(map_setting));
        if (map_lmax > grid_lmax or map_lmax < 2) {
          ERROR("NP_GEOMETRIC_LMAX must lie between 2 and the grid l_max");
        }
        const auto frame = gr::np::geometric_frame(
            *result.registration, *result.rapidity, result.adapted_rotation,
            spatial_metric, directions, *weights, *mass, map_lmax);
        auto geometric =
            fit_psi == nullptr
                ? gr::np::evaluate_geometric_second_order(*result.registration,
                                                          frame, *mass, moments)
                : gr::np::evaluate_geometric_second_order_with_fit_data(
                      *result.registration, frame, *mass,
                      Scalar<ComplexDataVector>{
                          gr::np::pull_back(*fit_psi,
                                            result.registration->rotation.a_bar,
                                            result.registration->rotation.b)
                              .get(4)});
        result.second_order = std::move(geometric.second_order);
      } else {
        result.second_order = gr::np::evaluate_second_order(
            *result.registration, *result.rapidity, result.adapted_rotation,
            *mass, weights, moments);
      }
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
