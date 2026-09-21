// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/SphericalShells.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Phi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Pi.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

namespace TestHelpers::gh_worldtube {
using EvolvedVariables = Variables<
    tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
               gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>>;

/// The evolved generalized harmonic variables of a Kerr-Schild solution at
/// `coords` and `time`
inline EvolvedVariables evolved_variables(
    const gr::Solutions::KerrSchild& solution,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords, const double time) {
  using frame = Frame::Inertial;
  using tags = tmpl::list<
      gr::Tags::Lapse<DataVector>, ::Tags::dt<gr::Tags::Lapse<DataVector>>,
      gr::Solutions::KerrSchild::DerivLapse<DataVector, frame>,
      gr::Tags::Shift<DataVector, 3, frame>,
      ::Tags::dt<gr::Tags::Shift<DataVector, 3, frame>>,
      gr::Solutions::KerrSchild::DerivShift<DataVector, frame>,
      gr::Tags::SpatialMetric<DataVector, 3, frame>,
      ::Tags::dt<gr::Tags::SpatialMetric<DataVector, 3, frame>>,
      gr::Solutions::KerrSchild::DerivSpatialMetric<DataVector, frame>>;
  const auto vars = solution.variables(coords, time, tags{});
  const auto& lapse = get<gr::Tags::Lapse<DataVector>>(vars);
  const auto& shift = get<gr::Tags::Shift<DataVector, 3, frame>>(vars);
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataVector, 3, frame>>(vars);
  EvolvedVariables result(get_size(get(lapse)));
  get<gr::Tags::SpacetimeMetric<DataVector, 3>>(result) =
      gr::spacetime_metric(lapse, shift, spatial_metric);
  get<gh::Tags::Phi<DataVector, 3>>(result) = gh::phi(
      lapse,
      get<gr::Solutions::KerrSchild::DerivLapse<DataVector, frame>>(vars),
      shift,
      get<gr::Solutions::KerrSchild::DerivShift<DataVector, frame>>(vars),
      spatial_metric,
      get<gr::Solutions::KerrSchild::DerivSpatialMetric<DataVector, frame>>(
          vars));
  get<gh::Tags::Pi<DataVector, 3>>(result) = gh::pi(
      lapse, get<::Tags::dt<gr::Tags::Lapse<DataVector>>>(vars), shift,
      get<::Tags::dt<gr::Tags::Shift<DataVector, 3, frame>>>(vars),
      spatial_metric,
      get<::Tags::dt<gr::Tags::SpatialMetric<DataVector, 3, frame>>>(vars),
      get<gh::Tags::Phi<DataVector, 3>>(result));
  return result;
}

/// \brief The Kretschmann scalar \f$48 M^2 / r'^6\f$ of a non-spinning
/// Kerr-Schild hole of mass 1 boosted with `velocity` from the origin, and
/// its spatial gradient and time derivative at fixed inertial coordinates.
///
/// Here \f$r'\f$ is the radius in the rest frame of the hole,
/// \f$r'^2 = |x_\perp|^2 + \gamma^2 (x_\parallel - v t)^2\f$; the scalar is
/// stationary in that frame, so \f$\partial_t K = -v^i \partial_i K\f$.
struct BoostedKretschmann {
  Scalar<DataVector> kretschmann;
  tnsr::i<DataVector, 3, Frame::Inertial> d_kretschmann;
  Scalar<DataVector> dt_kretschmann;
};

inline BoostedKretschmann boosted_kretschmann(
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const std::array<double, 3>& velocity, const double time) {
  const size_t num_points = get_size(get<0>(coords));
  const double speed_squared =
      square(velocity[0]) + square(velocity[1]) + square(velocity[2]);
  const double lorentz_factor = 1. / sqrt(1. - speed_squared);
  DataVector parallel(num_points, 0.);
  if (speed_squared > 0.) {
    for (size_t i = 0; i < 3; ++i) {
      parallel += coords.get(i) * gsl::at(velocity, i);
    }
    parallel /= sqrt(speed_squared);
    parallel -= sqrt(speed_squared) * time;
  }
  DataVector rest_radius_squared(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    rest_radius_squared += square(coords.get(i));
  }
  // |x_perp|^2 + gamma^2 (x_par - v t)^2 with |x|^2 = |x_perp|^2 + x_par^2 at
  // t = 0; for t != 0 restore x_par = x.v_hat
  DataVector x_parallel(num_points, 0.);
  if (speed_squared > 0.) {
    for (size_t i = 0; i < 3; ++i) {
      x_parallel += coords.get(i) * gsl::at(velocity, i);
    }
    x_parallel /= sqrt(speed_squared);
  }
  rest_radius_squared +=
      -square(x_parallel) + square(lorentz_factor) * square(parallel);
  BoostedKretschmann result{};
  get(result.kretschmann) = 48. / cube(rest_radius_squared);
  // d_i r'^2 = 2 x_perp_i + 2 gamma^2 (x_par - v t) v_hat_i
  //          = 2 x_i - 2 x_par v_hat_i + 2 gamma^2 (x_par - v t) v_hat_i
  const DataVector d_kretschmann_d_rest_radius_squared =
      -3. * 48. / pow<4>(rest_radius_squared);
  result.d_kretschmann =
      tnsr::i<DataVector, 3, Frame::Inertial>(num_points, 0.);
  get(result.dt_kretschmann) = DataVector(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    DataVector d_rest_radius_squared = 2. * coords.get(i);
    if (speed_squared > 0.) {
      const double v_hat = gsl::at(velocity, i) / sqrt(speed_squared);
      d_rest_radius_squared +=
          (-2. * x_parallel + 2. * square(lorentz_factor) * parallel) * v_hat;
    }
    result.d_kretschmann.get(i) =
        d_kretschmann_d_rest_radius_squared * d_rest_radius_squared;
    get(result.dt_kretschmann) -=
        gsl::at(velocity, i) * result.d_kretschmann.get(i);
  }
  return result;
}

/// \brief A boosted Kerr-Schild hole of mass 1 excised by a domain, and one
/// element of that domain, see `wedge_element()` and `shell_element()`.
struct KerrSchildElement {
  KerrSchildElement(Domain<3> domain_in, const Mesh<3>& mesh_in,
                    const ElementId<3>& element_id,
                    const std::array<double, 3>& velocity_in)
      : velocity(velocity_in),
        solution(1., {{0., 0., 0.}}, {{0., 0., 0.}}, velocity_in),
        domain(std::move(domain_in)),
        element(element_id, {}),
        mesh(mesh_in),
        element_map(element.id(), domain.blocks()[element_id.block_id()]),
        logical_coords(logical_coordinates(mesh)),
        inertial_coords(element_map(logical_coords)),
        inverse_jacobian(element_map.inv_jacobian(logical_coords)) {}

  /// The evolved variables at `time`. SpECTRE's boosted Kerr-Schild solution
  /// is implemented at \f$t = 0\f$ only, so the hole at `time` is the one
  /// whose center has moved to \f$v t\f$, which is exact for a boosted
  /// stationary solution.
  EvolvedVariables evolved_variables(const double time) const {
    const gr::Solutions::KerrSchild solution_at_time{
        1.,
        {{0., 0., 0.}},
        {{velocity[0] * time, velocity[1] * time, velocity[2] * time}},
        velocity};
    return gh_worldtube::evolved_variables(solution_at_time, inertial_coords,
                                           0.);
  }

  std::array<double, 3> velocity;
  gr::Solutions::KerrSchild solution;
  Domain<3> domain;
  Element<3> element;
  Mesh<3> mesh;
  ElementMap<3, Frame::Inertial> element_map;
  tnsr::I<DataVector, 3, Frame::ElementLogical> logical_coords;
  tnsr::I<DataVector, 3, Frame::Inertial> inertial_coords;
  InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian;
};

/// A wedge element of a `Sphere` domain excising the hole, with a
/// Legendre-Gauss-Lobatto mesh
inline KerrSchildElement wedge_element(
    const size_t num_points_per_dim, const std::array<double, 3>& velocity,
    const double inner_radius = 2.5, const double outer_radius = 4.0,
    const size_t block_id = 0,
    const std::array<SegmentId, 3>& segment_ids = {
        {SegmentId{0, 0}, SegmentId{0, 0}, SegmentId{0, 0}}}) {
  return KerrSchildElement{
      domain::creators::Sphere{inner_radius, outer_radius,
                               domain::creators::Sphere::Excision{}, 0_st,
                               num_points_per_dim, true}
          .create_domain(),
      Mesh<3>{num_points_per_dim, Spectral::Basis::Legendre,
              Spectral::Quadrature::GaussLobatto},
      ElementId<3>{block_id, segment_ids}, velocity};
}

/// The single shell element of a `SphericalShells` domain excising the hole,
/// whose spherical-harmonic angular mesh owns the whole excision face
inline KerrSchildElement shell_element(const size_t num_radial_points,
                                       const size_t l_max,
                                       const std::array<double, 3>& velocity,
                                       const double inner_radius = 2.5,
                                       const double outer_radius = 5.0) {
  return KerrSchildElement{
      domain::creators::SphericalShells{inner_radius, outer_radius, 0_st,
                                        num_radial_points, l_max}
          .create_domain(),
      Mesh<3>{{{num_radial_points, l_max + 1, 2 * l_max + 1}},
              {{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
                Spectral::Basis::SphericalHarmonic}},
              {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::Gauss,
                Spectral::Quadrature::Equiangular}}},
      ElementId<3>{0}, velocity};
}
}  // namespace TestHelpers::gh_worldtube
