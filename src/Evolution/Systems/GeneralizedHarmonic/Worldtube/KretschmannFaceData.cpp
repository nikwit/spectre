// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <numbers>
#include <optional>
#include <pup.h>
#include <string>
#include <unordered_map>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/IndexIterator.hpp"
#include "DataStructures/SliceTensorToVariables.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/WeylCurvature.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

namespace gh::worldtube {
namespace {
struct KretschmannTag : db::SimpleTag {
  using type = Scalar<DataVector>;
};
template <size_t Dim>
struct DerivKretschmannTag : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim, Frame::Inertial>;
};
template <size_t Dim>
struct MeshVelocityTag : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;
};
}  // namespace

void validate_geometric_matching_face(
    const tnsr::I<DataVector, 3, Frame::Inertial>& coordinates,
    const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Direction<3>& direction) {
  if (direction != Direction<3>::lower_xi() or
      mesh.basis(1) != Spectral::Basis::SphericalHarmonic or
      mesh.basis(2) != Spectral::Basis::SphericalHarmonic or
      mesh.quadrature(1) != Spectral::Quadrature::Gauss or
      mesh.quadrature(2) != Spectral::Quadrature::Equiangular or
      mesh.extents(2) != 2 * mesh.extents(1) - 1) {
    ERROR(
        "QuadrupoleGeometric requires a complete spherical-harmonic inner "
        "face. Enable UseSphericalHarmonics on the worldtube object.");
  }
  const auto weights = face_quadrature_weights<3>(mesh.slice_away(0));
  const size_t nr = mesh.extents(0), np = weights.size();
  std::array<double, 3> center{};
  for (size_t p = 0; p < np; ++p) {
    for (size_t i = 0; i < 3; ++i) {
      center[i] += weights[p] * coordinates.get(i)[nr * p];
    }
  }
  DataVector radius(np, 0.);
  double mean_radius = 0.;
  for (size_t p = 0; p < np; ++p) {
    for (size_t i = 0; i < 3; ++i) {
      radius[p] += square(coordinates.get(i)[nr * p] - center[i]);
    }
    radius[p] = sqrt(radius[p]);
    mean_radius += weights[p] * radius[p];
  }
  if (not(mean_radius > 0. and std::isfinite(mean_radius))) {
    ERROR("QuadrupoleGeometric requires a finite nonzero coordinate radius.");
  }
  for (size_t p = 0; p < np; ++p) {
    double normal_norm = 0.;
    for (size_t i = 0; i < 3; ++i) {
      normal_norm += square(inverse_jacobian.get(0, i)[nr * p]);
    }
    normal_norm = sqrt(normal_norm);
    if (not(normal_norm > 0. and std::isfinite(normal_norm)) or
        not(std::abs(radius[p] / mean_radius - 1.) < 1.e-9)) {
      ERROR(
          "QuadrupoleGeometric requires a round inertial-coordinate sphere; "
          "disable non-spherical shape/skew deformations on the worldtube.");
    }
    for (size_t i = 0; i < 3; ++i) {
      if (not(std::abs(inverse_jacobian.get(0, i)[nr * p] / normal_norm -
                       (coordinates.get(i)[nr * p] - center[i]) / radius[p]) <
              1.e-9)) {
        ERROR(
            "QuadrupoleGeometric requires radial coordinate-normal covectors.");
      }
    }
  }
}

void ThirdOrderFaceData::pup(PUP::er& p) {
  p | history;
  p | time;
  p | psi0_target;
  p | moments;
  p | dots;
  p | derivative_order;
  p | clock_rate;
  p | tilt_residual;
  p | fit_residual;
  p | fit_condition;
}
bool operator==(const ThirdOrderFaceData& a, const ThirdOrderFaceData& b) {
  return a.history == b.history and a.time == b.time and
         a.psi0_target == b.psi0_target and a.moments == b.moments and
         a.dots == b.dots and a.derivative_order == b.derivative_order and
         a.clock_rate == b.clock_rate and a.tilt_residual == b.tilt_residual and
         a.fit_residual == b.fit_residual and
         a.fit_condition == b.fit_condition;
}

void update_third_order_matching(
    const gsl::not_null<KretschmannFaceData<3>*> data,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coordinates,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& velocity,
    const double mass, const double time, const TimeStepId& time_id) {
  if (not data->direction.has_value() or time != data->time) {
    ERROR("Update the Kretschmann data before third-order matching");
  }
  if (not time_id.time_runs_forward()) {
    ERROR("ThirdOrderGeometric supports forward evolution only");
  }
  validate_geometric_matching_face(coordinates, mesh, inverse_jacobian,
                                   *data->direction);
  const auto face = face_curvature(spacetime_metric, pi, phi, mesh,
                                   inverse_jacobian, *data->direction);
  const auto& weights = data->quadrature_weights;
  const size_t np = weights.size(), nr = mesh.extents(0);
  gr::np::TriadVector labels(np, 0.), inward(np, 0.), mesh_velocity(np, 0.);
  std::array<double, 3> center{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t p = 0; p < np; ++p) {
      center[i] += weights[p] * coordinates.get(i)[nr * p];
    }
  }
  double coordinate_radius = 0.;
  for (size_t p = 0; p < np; ++p) {
    double radius = 0.;
    for (size_t i = 0; i < 3; ++i) {
      labels.get(i)[p] = coordinates.get(i)[nr * p] - center[i];
      radius += square(labels.get(i)[p]);
      if (velocity)
        mesh_velocity.get(i)[p] = velocity->get(i)[nr * p];
    }
    radius = sqrt(radius);
    coordinate_radius += weights[p] * radius;
    for (size_t i = 0; i < 3; ++i) {
      labels.get(i)[p] /= radius;
      inward.get(i)[p] = -labels.get(i)[p];
    }
  }
  const auto rotation = gr::np::adapted_triad(face.spatial_metric, inward);
  const auto psi = gr::np::weyl_scalars_from_electric_magnetic(
      face.electric, face.magnetic, face.spatial_metric, inward);
  const auto registration = gr::np::register_frame(psi, rotation, mass);
  const auto rapidity = gr::np::invariant_rapidity(
      registration.member, rotation, face.spatial_metric, face.lapse,
      face.shift, data->d_kretschmann, data->dt_kretschmann);
  const size_t grid_lmax = mesh.extents(1) - 1;
  const char* setting = std::getenv("NP_GEOMETRIC_LMAX");
  const size_t lmax =
      setting == nullptr ? std::min(size_t{8}, grid_lmax) : std::stoul(setting);
  if (lmax < 2 or lmax > grid_lmax) {
    ERROR("NP_GEOMETRIC_LMAX must lie between 2 and the grid l_max");
  }
  const auto frame =
      gr::np::geometric_frame(registration, rapidity, rotation,
                              face.spatial_metric, labels, weights, mass, lmax);
  const auto temporal = gr::np::geometric_time_data(
      frame, labels, coordinate_radius, registration.measured_radius, mass,
      rotation, face.spatial_metric, face.lapse, face.shift, mesh_velocity);
  const auto columns = gr::np::third_order_tide_columns(
      frame, registration.measured_radius, mass, temporal.slice_tilt);
  // This preliminary undotted fit is the same history variable as offline.
  // Differentiating its O(epsilon^3) contamination affects fourth order;
  // differentiating the corrected fit would create an implicit feedback loop.
  const auto preliminary = gr::np::fit_third_order_tide(registration, columns,
                                                        frame.map.weights, {});
  gr::np::GeometricTideSample sample{};
  sample.time = time;
  sample.direction = frame.map.direction;
  std::copy_n(preliminary.moments.begin(), 10, sample.quadrupole.begin());
  if (not data->third_order)
    data->third_order.emplace();
  auto& state = *data->third_order;
  const auto derivative = gr::np::causal_tidal_derivative(
      make_not_null(&state.history), sample, frame.map, temporal,
      time_id.step_time().value(), time_id.substep() == 0);
  auto fit = gr::np::fit_third_order_tide(registration, columns,
                                          frame.map.weights, derivative.dots);
  state.time = time;
  state.psi0_target = std::move(fit.psi0_target);
  state.moments = fit.moments;
  state.dots = derivative.dots;
  state.derivative_order = derivative.derivative_order;
  state.clock_rate = temporal.clock_rate;
  state.tilt_residual = temporal.tilt_residual;
  state.fit_residual = fit.relative_residual;
  state.fit_condition = fit.condition_number;
  data->filtered_moments.reset();
}

template <size_t Dim>
void KretschmannFaceData<Dim>::pup(PUP::er& p) {
  p | direction;
  p | time;
  p | kretschmann;
  p | d_kretschmann;
  p | dt_kretschmann;
  p | quadrature_weights;
  p | previous_time;
  p | previous_kretschmann;
  p | filtered_moments;
  p | filtered_moments_time;
  p | third_order;
  p | initial_gauge_difference;
  p | radial_gauge;
  p | replay;
}

template <size_t Dim>
bool operator==(const KretschmannFaceData<Dim>& lhs,
                const KretschmannFaceData<Dim>& rhs) {
  // NaN times compare equal when both are NaN, so that default-constructed
  // objects are equal
  const auto same_time = [](const double a, const double b) {
    return (std::isnan(a) and std::isnan(b)) or a == b;
  };
  return lhs.third_order == rhs.third_order and lhs.replay == rhs.replay and
         lhs.direction == rhs.direction and same_time(lhs.time, rhs.time) and
         lhs.kretschmann == rhs.kretschmann and
         lhs.d_kretschmann == rhs.d_kretschmann and
         lhs.dt_kretschmann == rhs.dt_kretschmann and
         lhs.quadrature_weights == rhs.quadrature_weights and
         same_time(lhs.previous_time, rhs.previous_time) and
         lhs.previous_kretschmann == rhs.previous_kretschmann and
         lhs.radial_gauge == rhs.radial_gauge and
         lhs.initial_gauge_difference == rhs.initial_gauge_difference and
         lhs.filtered_moments == rhs.filtered_moments and
         same_time(lhs.filtered_moments_time, rhs.filtered_moments_time);
}

template <size_t Dim>
bool operator!=(const KretschmannFaceData<Dim>& lhs,
                const KretschmannFaceData<Dim>& rhs) {
  return not(lhs == rhs);
}

template <size_t Dim>
std::optional<Direction<Dim>> excision_face_direction(
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const Element<Dim>& element) {
  for (const auto& [name, excision_sphere] : excision_spheres) {
    (void)name;
    if (const auto direction = excision_sphere.abutting_direction(element.id());
        direction.has_value()) {
      return direction;
    }
  }
  return std::nullopt;
}

template <size_t Dim>
DataVector face_quadrature_weights(const Mesh<Dim - 1>& face_mesh) {
  if constexpr (Dim == 1) {
    (void)face_mesh;
    return DataVector{1.};
  } else {
    DataVector weights(face_mesh.number_of_grid_points(), 1.);
    for (size_t d = 0; d < Dim - 1; ++d) {
      const Mesh<1> mesh_1d = face_mesh.slice_through(d);
      const size_t num_points_1d = mesh_1d.number_of_grid_points();
      DataVector weights_1d(num_points_1d, 0.);
      if (mesh_1d.basis(0) == Spectral::Basis::SphericalHarmonic and
          mesh_1d.quadrature(0) == Spectral::Quadrature::Gauss) {
        // Gauss-Legendre in cos(theta); the weights are symmetric, so the
        // ordering of the theta collocation points does not matter
        weights_1d = Spectral::quadrature_weights<Spectral::Basis::Legendre,
                                                  Spectral::Quadrature::Gauss>(
            num_points_1d);
      } else if (mesh_1d.basis(0) == Spectral::Basis::SphericalHarmonic or
                 mesh_1d.basis(0) == Spectral::Basis::Fourier) {
        weights_1d = 2. * std::numbers::pi / static_cast<double>(num_points_1d);
      } else {
        weights_1d = Spectral::quadrature_weights(mesh_1d);
      }
      for (IndexIterator<Dim - 1> index(face_mesh.extents()); index; ++index) {
        weights[index.collapsed_index()] *= weights_1d[index()[d]];
      }
    }
    double sum = 0.;
    for (const double weight : weights) {
      sum += weight;
    }
    weights /= sum;
    return weights;
  }
}

template <size_t Dim>
void update_kretschmann_face_data(
    const gsl::not_null<KretschmannFaceData<Dim>*> data,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Element<Dim>& element,
    const std::unordered_map<std::string, ExcisionSphere<Dim>>&
        excision_spheres,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        mesh_velocity,
    const double time) {
  const auto direction = excision_face_direction(excision_spheres, element);
  if (not direction.has_value()) {
    *data = KretschmannFaceData<Dim>{};
    return;
  }
  data->direction = direction;
  const size_t sliced_dim = direction->dimension();
  const size_t fixed_index = index_to_slice_at(mesh.extents(), *direction);
  const Mesh<Dim - 1> face_mesh = mesh.slice_away(sliced_dim);
  const size_t num_face_points = face_mesh.number_of_grid_points();
  data->quadrature_weights = face_quadrature_weights<Dim>(face_mesh);

  if constexpr (Dim != 3) {
    (void)spacetime_metric;
    (void)pi;
    (void)phi;
    (void)inverse_jacobian;
    (void)mesh_velocity;
    (void)fixed_index;
    (void)num_face_points;
    data->time = time;
    return;
  } else {
    if (not data->initial_gauge_difference.has_value()) {
      const auto face = face_curvature(spacetime_metric, pi, phi, mesh,
                                       inverse_jacobian, *direction);
      const auto sliced_phi = data_on_slice<gh::Tags::Phi<DataVector, Dim>>(
          mesh.extents(), sliced_dim, fixed_index, phi);
      const auto& face_phi = get<gh::Tags::Phi<DataVector, Dim>>(sliced_phi);
      tnsr::I<DataVector, Dim, Frame::Inertial> normal(num_face_points, 0.);
      for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = 0; j < Dim; ++j) {
          normal.get(i) += face.inverse_spatial_metric.get(i, j) *
                           face.unit_normal_covector.get(j);
        }
      }
      tnsr::A<DataVector, Dim, Frame::Inertial> ell(num_face_points, 0.);
      get<0>(ell) = 1. / (sqrt(2.) * get(face.lapse));
      for (size_t i = 0; i < Dim; ++i) {
        ell.get(i + 1) =
            (-face.shift.get(i) / get(face.lapse) + normal.get(i)) / sqrt(2.);
      }
      data->initial_gauge_difference.emplace(num_face_points, 0.);
      for (size_t a = 0; a <= Dim; ++a) {
        for (size_t b = 0; b <= Dim; ++b) {
          for (size_t i = 0; i < Dim; ++i) {
            data->initial_gauge_difference->get(a) -=
                2. * ell.get(b) * normal.get(i) * face_phi.get(i, a, b);
          }
        }
      }
    }
    // Curvature and Kretschmann scalar in the volume, then its gradient
    const auto d_pi = partial_derivative(pi, mesh, inverse_jacobian);
    const auto d_phi = partial_derivative(phi, mesh, inverse_jacobian);
    const WeylCurvature curvature =
        weyl_curvature(spacetime_metric, pi, phi, d_pi, d_phi);
    const Scalar<DataVector> kretschmann =
        kretschmann_scalar(curvature.electric, curvature.magnetic,
                           curvature.inverse_spatial_metric);
    const auto d_kretschmann =
        partial_derivative(kretschmann, mesh, inverse_jacobian);

    const auto face_vars =
        data_on_slice<KretschmannTag, DerivKretschmannTag<Dim>>(
            mesh.extents(), sliced_dim, fixed_index, kretschmann,
            d_kretschmann);
    const Scalar<DataVector>& face_kretschmann = get<KretschmannTag>(face_vars);

    // Time derivative at fixed inertial coordinates from the backward
    // difference at fixed grid points, corrected for the mesh velocity
    const bool have_history =
        get(data->previous_kretschmann).size() == num_face_points and
        time > data->previous_time;
    const bool same_time = time == data->previous_time;
    const bool reuse_previous_estimate =
        same_time and get(data->dt_kretschmann).size() == num_face_points;
    if (not reuse_previous_estimate) {
      DataVector grid_rate(num_face_points, 0.);
      if (have_history) {
        grid_rate = (get(face_kretschmann) - get(data->previous_kretschmann)) /
                    (time - data->previous_time);
      }
      if (mesh_velocity.has_value()) {
        const auto face_velocity = data_on_slice<MeshVelocityTag<Dim>>(
            mesh.extents(), sliced_dim, fixed_index, *mesh_velocity);
        for (size_t i = 0; i < Dim; ++i) {
          grid_rate -= get<MeshVelocityTag<Dim>>(face_velocity).get(i) *
                       get<DerivKretschmannTag<Dim>>(face_vars).get(i);
        }
      }
      get(data->dt_kretschmann) = std::move(grid_rate);
    }
    if (not same_time) {
      data->previous_kretschmann = face_kretschmann;
      data->previous_time = time;
    }
    data->kretschmann = face_kretschmann;
    data->d_kretschmann = get<DerivKretschmannTag<Dim>>(face_vars);
    data->time = time;
  }
}

void relax_tidal_moments(
    const gsl::not_null<std::optional<gr::np::TidalMoments>*> moments,
    const gsl::not_null<double*> moments_time, const gr::np::TidalMoments& raw,
    const double time, const std::optional<double>& relaxation_time) {
  if (not relaxation_time.has_value() or not moments->has_value() or
      not(time >= *moments_time)) {
    *moments = raw;
    *moments_time = time;
    return;
  }
  if (time == *moments_time) {
    return;
  }
  const double factor = std::min((time - *moments_time) / *relaxation_time, 1.);
  for (size_t a = 0; a < 5; ++a) {
    gsl::at(**moments, a) += factor * (gsl::at(raw, a) - gsl::at(**moments, a));
  }
  *moments_time = time;
}

template <size_t Dim>
void update_filtered_tidal_moments(
    const gsl::not_null<KretschmannFaceData<Dim>*> data,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
    const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
    const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const PhysicalModel model, const double mass,
    const std::optional<double>& relaxation_time, const double time) {
  if constexpr (Dim != 3) {
    (void)data;
    (void)spacetime_metric;
    (void)pi;
    (void)phi;
    (void)mesh;
    (void)inverse_jacobian;
    (void)model;
    (void)mass;
    (void)relaxation_time;
    (void)time;
    ERROR("The order-two worldtube model is only implemented in 3 dimensions.");
  } else {
    if (not data->direction.has_value()) {
      ERROR(
          "The tidal moments can only be updated on an element that abuts an "
          "excision sphere; update the Kretschmann face data first.");
    }
    const FaceCurvature face = face_curvature(
        spacetime_metric, pi, phi, mesh, inverse_jacobian, *data->direction);
    if (not is_order_two(model)) {
      ERROR("Tidal moments are only defined for the order-two models, not "
            << model);
    }
    const MatchingEvaluation raw =
        evaluate_matching(model, mass, face.electric, face.magnetic,
                          face.spatial_metric, face.unit_normal_covector,
                          face.lapse, face.shift, &*data, std::nullopt);
    relax_tidal_moments(make_not_null(&data->filtered_moments),
                        make_not_null(&data->filtered_moments_time),
                        raw.second_order->fit.components, time,
                        relaxation_time);
  }
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                               \
  template struct KretschmannFaceData<DIM(data)>;                            \
  template bool operator==(const KretschmannFaceData<DIM(data)>& lhs,        \
                           const KretschmannFaceData<DIM(data)>& rhs);       \
  template bool operator!=(const KretschmannFaceData<DIM(data)>& lhs,        \
                           const KretschmannFaceData<DIM(data)>& rhs);       \
  template std::optional<Direction<DIM(data)>> excision_face_direction(      \
      const std::unordered_map<std::string, ExcisionSphere<DIM(data)>>&      \
          excision_spheres,                                                  \
      const Element<DIM(data)>& element);                                    \
  template DataVector face_quadrature_weights<DIM(data)>(                    \
      const Mesh<DIM(data) - 1>& face_mesh);                                 \
  template void update_filtered_tidal_moments(                               \
      gsl::not_null<KretschmannFaceData<DIM(data)>*>,                        \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&,               \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&,               \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>&,              \
      const Mesh<DIM(data)>&,                                                \
      const InverseJacobian<DataVector, DIM(data), Frame::ElementLogical,    \
                            Frame::Inertial>&,                               \
      PhysicalModel, double, const std::optional<double>&, double);          \
  template void update_kretschmann_face_data(                                \
      gsl::not_null<KretschmannFaceData<DIM(data)>*>,                        \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&,               \
      const tnsr::aa<DataVector, DIM(data), Frame::Inertial>&,               \
      const tnsr::iaa<DataVector, DIM(data), Frame::Inertial>&,              \
      const Mesh<DIM(data)>&,                                                \
      const InverseJacobian<DataVector, DIM(data), Frame::ElementLogical,    \
                            Frame::Inertial>&,                               \
      const Element<DIM(data)>&,                                             \
      const std::unordered_map<std::string, ExcisionSphere<DIM(data)>>&,     \
      const std::optional<tnsr::I<DataVector, DIM(data), Frame::Inertial>>&, \
      double);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace gh::worldtube
