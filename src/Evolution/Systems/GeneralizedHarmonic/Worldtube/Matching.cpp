// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
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
  } else if (read == "QuadrupoleGeometric") {
    return PhysicalModel::QuadrupoleGeometric;
  } else if (read == "QuadrupoleCoulomb") {
    return PhysicalModel::QuadrupoleCoulomb;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert input option to a physical model. Must be "
              "one of None, TypeD, Quadrupole, QuadrupoleGeometric or "
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
    case PhysicalModel::QuadrupoleGeometric:
      return os << "QuadrupoleGeometric";
    case PhysicalModel::QuadrupoleCoulomb:
      return os << "QuadrupoleCoulomb";
    default:
      ERROR("Unknown PhysicalModel");
  }
}

bool is_order_two(const PhysicalModel model) {
  return model == PhysicalModel::QuadrupoleGeometric or
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
    case PhysicalModel::Quadrupole:
    case PhysicalModel::QuadrupoleGeometric:
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
        // On a centered coordinate sphere the covector direction gives the
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
        auto geometric = gr::np::evaluate_geometric_second_order(
            *result.registration, *result.rapidity, result.adapted_rotation,
            spatial_metric, directions, *weights, *mass, map_lmax, moments);
        // Optional research diagnostics from the actual live prescription.
        // The historical observer still labels its passive legacy fit as
        // Quadrupole; its imposed-moment column must not be used for this
        // model.
        const char* diagnostic_path = std::getenv("NP_GEOMETRIC_DIAGNOSTICS");
        static thread_local double last_output =
            -std::numeric_limits<double>::infinity();
        if (diagnostic_path != nullptr and not moments.has_value() and
            face_data->time >= last_output + 0.099999) {
          std::ofstream out(diagnostic_path, std::ios::app);
          if (not out) {
            ERROR("Cannot open geometric diagnostic file");
          }
          if (out.tellp() == 0) {
            out << "time,lmax,lambda0,lambda1,lambda2,lambda3,lambda4,min_"
                   "jacobian,area_ratio,dyad_error,fit_residual,max_target,max_"
                   "error";
            for (size_t a = 0; a < 5; ++a)
              out << ",H" << a << "re,H" << a << "im";
            out << "\n";
          }
          out << std::setprecision(17) << face_data->time << ',' << map_lmax;
          for (const double value : geometric.map.eigenvalues)
            out << ',' << value;
          out << ',' << geometric.map.minimum_jacobian << ','
              << geometric.map.area_over_label_area << ','
              << geometric.maximum_dyad_error << ','
              << geometric.second_order.fit.relative_residual << ','
              << max(abs(get(geometric.second_order.psi0_target))) << ','
              << max(abs(get(geometric.second_order.psi0_target) -
                         result.psi.get(0)));
          for (const auto value : geometric.second_order.fit.components)
            out << ',' << value.real() << ',' << value.imag();
          out << "\n";
          last_output = face_data->time;
        }
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

  // Research controls for the stationary M=1, R=2.5 single-hole campaign.
  // The prescribed NR spin-2 pattern is fixed analytically; it is never
  // constructed from either model's fitted axes or moments.
  const char* seed_text = std::getenv("NP_CONTROL_SEED");
  const char* diagnostic_file = std::getenv("NP_MATCHING_DIAGNOSTICS");
  if (face_data != nullptr and
      (imposed_moments.has_value() or model == PhysicalModel::TypeD) and
      (seed_text != nullptr or diagnostic_file != nullptr)) {
    const double time = face_data->time;
    const double seed = seed_text == nullptr ? 0. : std::stod(seed_text);
    const char* parity = std::getenv("NP_CONTROL_PARITY");
    const bool magnetic_seed =
        parity != nullptr and std::string(parity) == "magnetic";
    const double envelope =
        time > 0. and time < 2.
            ? seed * std::pow(sin(std::acos(-1.) * time / 2.), 4)
            : 0.;
    std::array<Scalar<ComplexDataVector>, 5> fixed_columns;
    for (auto& column : fixed_columns)
      get(column) = ComplexDataVector(num_points, 0.);
    for (size_t p = 0; p < num_points; ++p) {
      const double nx = -unit_normal_covector.get(0)[p] / euclidean_norm[p];
      const double ny = -unit_normal_covector.get(1)[p] / euclidean_norm[p];
      const double nz = -unit_normal_covector.get(2)[p] / euclidean_norm[p];
      const double st = std::hypot(nx, ny);
      // m=(e_Theta-i e_Phi)/sqrt(2) for the inward NR normal.
      const std::complex<double> mx{nz * nx / (st * sqrt(2.)),
                                    ny / (st * sqrt(2.))};
      const std::complex<double> my{nz * ny / (st * sqrt(2.)),
                                    -nx / (st * sqrt(2.))};
      const std::complex<double> mz{-st / sqrt(2.), 0.};
      get(fixed_columns[0])[p] = 1.8 * (mx * mx - mz * mz);
      get(fixed_columns[1])[p] = 1.8 * (my * my - mz * mz);
      get(fixed_columns[2])[p] = 3.6 * mx * my;
      get(fixed_columns[3])[p] = 3.6 * mx * mz;
      get(fixed_columns[4])[p] = 3.6 * my * mz;
      get(result.psi0_target)[p] +=
          (magnetic_seed ? std::complex<double>{0., envelope}
                         : std::complex<double>{envelope, 0.}) *
          (get(fixed_columns[0])[p] - get(fixed_columns[1])[p]);
    }
    static thread_local std::array<double, 5> last_diagnostic{
        {-1., -1., -1., -1., -1.}};
    auto& last = last_diagnostic[static_cast<size_t>(model)];
    if (diagnostic_file != nullptr and time >= last + 0.049999) {
      const auto fit =
          gr::np::fit_psi4(Scalar<ComplexDataVector>{result.psi.get(0)},
                           fixed_columns, face_data->quadrature_weights);
      double rms0 = 0., rms4 = 0.;
      for (size_t p = 0; p < num_points; ++p) {
        rms0 +=
            face_data->quadrature_weights[p] * std::norm(result.psi.get(0)[p]);
        rms4 +=
            face_data->quadrature_weights[p] * std::norm(result.psi.get(4)[p]);
      }
      std::ofstream out(diagnostic_file, std::ios::app);
      if (not out)
        ERROR("Cannot open matching diagnostic file");
      if (out.tellp() == 0)
        out << "time,model,rms0,rms4,max_target,H0re,H0im,H1re,H1im,H2re,H2im,"
               "H3re,H3im,H4re,H4im\n";
      out << std::setprecision(17) << time << ',' << model << ',' << sqrt(rms0)
          << ',' << sqrt(rms4) << ',' << max(abs(get(result.psi0_target)));
      for (const auto value : fit.components)
        out << ',' << value.real() << ',' << value.imag();
      out << "\n";
      last = time;
    }
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
