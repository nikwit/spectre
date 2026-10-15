// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Actions/EvolvePhysicalTrace.hpp"

#include <cmath>

#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeEvolvedPhysicalGhost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/AnalyticChristoffel.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/ConstraintDamping/DampingFunction.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace gh::worldtube {
namespace {
// On the lower radial face the radial index is zero and is the fastest index.
template <typename Tensor>
Tensor lower_face(const Tensor& volume, const Mesh<3>& mesh) {
  const size_t n = mesh.extents(1) * mesh.extents(2);
  Tensor result(n);
  for (size_t c = 0; c < volume.size(); ++c) {
    for (size_t p = 0; p < n; ++p) {
      result[c][p] = volume[c][p * mesh.extents(0)];
    }
  }
  return result;
}

template <typename Tensor>
auto tangential_derivative(
    const Tensor& field, const ylm::Spherepack& spherepack,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& face_inverse_jacobian,
    const PhysicalTraceGeometry& geometry) {
  const size_t n = get(geometry.lapse).size();
  TensorMetafunctions::prepend_spatial_index<Tensor, 3, UpLo::Lo,
                                             Frame::Inertial>
      result(n, 0.);
  for (size_t c = 0; c < field.size(); ++c) {
    // The angular map Jacobian uses the same Pfaffian convention as Spherepack.
    const auto angular = spherepack.gradient(field[c]);
    std::array<DataVector, 3> derivative;
    DataVector normal_derivative(n, 0.);
    for (size_t i = 0; i < 3; ++i) {
      derivative[i] = face_inverse_jacobian.get(1, i) * angular[0] +
                      face_inverse_jacobian.get(2, i) * angular[1];
      normal_derivative += geometry.normal_vector.get(i) * derivative[i];
    }
    for (size_t i = 0; i < 3; ++i) {
      result.get(prepend(field.get_tensor_index(c), i)) =
          derivative[i] - geometry.normal_covector.get(i) * normal_derivative;
    }
  }
  return result;
}
}  // namespace

DataVector physical_trace_on_sphere(
    const bool initialize, const DataVector& amplitudes,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const tnsr::I<DataVector, 3>& coordinates,
    const tnsr::I<DataVector, 3, Frame::Grid>& grid_coordinates,
    const Domain<3>& domain, const Element<3>& element,
    const PhysicalTraceBoundaryConditions& boundary_conditions,
    const double time, const LtsMode lts_mode,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity,
    const gh::gauges::GaugeCondition& gauge,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma0_function,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma1_function,
    const ConstraintDamping::DampingFunction<3, Frame::Grid>& gamma2_function,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  const BoundaryConditions::WorldtubeEvolvedPhysicalGhost* boundary = nullptr;
  for (const auto& direction : element.external_boundaries()) {
    const auto* candidate =
        dynamic_cast<const BoundaryConditions::WorldtubeEvolvedPhysicalGhost*>(
            boundary_conditions.at(element.id().block_id())
                .at(direction)
                .get());
    if (candidate != nullptr) {
      if (direction != Direction<3>::lower_xi() or boundary != nullptr) {
        ERROR("Evolved physical trace requires one lower radial boundary.");
      }
      boundary = candidate;
    }
  }
  if (boundary == nullptr) {
    return {};
  }
  if (lts_mode != LtsMode::Off or
      (mesh_velocity.has_value() and not boundary->uses_recorded_weyl()) or
      mesh.basis(1) != Spectral::Basis::SphericalHarmonic or
      mesh.basis(2) != Spectral::Basis::SphericalHarmonic or
      mesh.quadrature(0) != Spectral::Quadrature::GaussLobatto or
      mesh.extents(2) != 2 * mesh.extents(1) - 1 or
      element.id().segment_id(1).refinement_level() != 0 or
      element.id().segment_id(2).refinement_level() != 0) {
    ERROR(
        "Evolved physical trace requires a full-sphere face and global "
        "time stepping. Moving meshes require RecordedWeyl. Angular "
        "subdomains and LTS are unsupported.");
  }
  const auto* analytic_gauge =
      dynamic_cast<const gh::gauges::AnalyticChristoffel*>(&gauge);
  const auto* damped_harmonic =
      dynamic_cast<const gh::gauges::DampedHarmonic*>(&gauge);
  if (not gauge.is_harmonic() and analytic_gauge == nullptr and
      damped_harmonic == nullptr) {
    ERROR(
        "Evolved physical trace supports Harmonic, DampedHarmonic and "
        "AnalyticChristoffel gauge only.");
  }
  const auto g = lower_face(metric, mesh);
  const auto receiver_pi = lower_face(pi, mesh);
  const auto receiver_phi = lower_face(phi, mesh);
  const auto x = lower_face(coordinates, mesh);
  const size_t n = get<0>(x).size();
  const std::optional<tnsr::I<DataVector, 3>> face_velocity =
      mesh_velocity.has_value()
          ? std::make_optional(lower_face(*mesh_velocity, mesh))
          : std::nullopt;
  const auto geometry_at = [&](const tnsr::aa<DataVector, 3>& dt_g) {
    return boundary->uses_recorded_weyl()
               ? mapped_sphere_trace_geometry(
                     g, dt_g, x,
                     face_velocity.value_or(tnsr::I<DataVector, 3>(n, 0.)),
                     mesh.extents(1) - 1)
               : spherical_trace_geometry(g, dt_g, x);
  };
  auto geometry = geometry_at(tnsr::aa<DataVector, 3>(n, 0.));
  if (initialize) {
    return physical_trace_amplitudes(receiver_pi, receiver_phi, geometry);
  }
  // For gamma1=-1 the receiver metric RHS is local, including mesh
  // advection and its native GH constraint adjustment.
  tnsr::aa<DataVector, 3> dt_g(n);
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = a; b < 4; ++b) {
      dt_g.get(a, b) = -get(geometry.lapse) * receiver_pi.get(a, b);
      for (size_t i = 0; i < 3; ++i) {
        dt_g.get(a, b) += geometry.shift.get(i) * receiver_phi.get(i, a, b);
        if (face_velocity.has_value()) {
          dt_g.get(a, b) += face_velocity->get(i) * receiver_phi.get(i, a, b);
        }
      }
    }
  }
  geometry = geometry_at(dt_g);
  const auto grid_x = lower_face(grid_coordinates, mesh);
  Scalar<DataVector> gamma0(n), gamma1(n), gamma2(n);
  gamma0_function(make_not_null(&gamma0), grid_x, time, functions_of_time);
  gamma1_function(make_not_null(&gamma1), grid_x, time, functions_of_time);
  gamma2_function(make_not_null(&gamma2), grid_x, time, functions_of_time);
  if (min(get(gamma1)) != -1. or max(get(gamma1)) != -1.) {
    ERROR("Evolved physical trace requires gamma1=-1.");
  }
  const PhysicalScreen screen(g, geometry.inverse_metric, geometry.time_normal,
                              geometry.normal_vector);
  tnsr::aa<DataVector, 3> ghost_pi(n);
  tnsr::iaa<DataVector, 3> ghost_phi(n);
  physical_ghost(make_not_null(&ghost_pi), make_not_null(&ghost_phi),
                 receiver_pi, receiver_phi,
                 physical_trace_tensor(amplitudes, geometry),
                 geometry.normal_covector, geometry.normal_vector, screen);
  const auto face_jac = lower_face(inverse_jacobian, mesh);
  const auto& spherepack = ylm::get_spherepack_cache(mesh.extents(1) - 1);
  const auto dg = tangential_derivative(g, spherepack, face_jac, geometry);
  const auto dpi =
      tangential_derivative(ghost_pi, spherepack, face_jac, geometry);
  const auto dphi =
      tangential_derivative(ghost_phi, spherepack, face_jac, geometry);
  tnsr::a<DataVector, 3> gauge_h(n, 0.);
  tnsr::ab<DataVector, 3> d_gauge_h(n, 0.);
  if (analytic_gauge != nullptr) {
    // Differentiate the prescribed gauge on the original VOLUME grid. Setting
    // its radial derivatives to zero would change the gauge forcing.
    tnsr::a<DataVector, 3> volume_h(mesh.number_of_grid_points());
    tnsr::ab<DataVector, 3> volume_dh(mesh.number_of_grid_points());
    analytic_gauge->gauge_and_spacetime_derivative(
        make_not_null(&volume_h), make_not_null(&volume_dh), mesh, time,
        coordinates, inverse_jacobian, gh::Solutions::all_solutions<3>{});
    gauge_h = lower_face(volume_h, mesh);
    d_gauge_h = lower_face(volume_dh, mesh);
  }
  const auto target = boundary->curvature_target(
      time, x, screen, &domain, &element, &functions_of_time, &mesh);
  return physical_boundary_rhs(g, ghost_pi, ghost_phi, dg, dpi, dphi, gamma0,
                               gamma2, gauge_h, d_gauge_h, target, geometry,
                               face_velocity, damped_harmonic, &x);
}
}  // namespace gh::worldtube
