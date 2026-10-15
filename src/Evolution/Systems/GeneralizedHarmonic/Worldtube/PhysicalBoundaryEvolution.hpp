// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/PhysicalBoundaryTrace.hpp"

namespace gh::gauges {
class DampedHarmonic;
}

namespace gh::worldtube {
namespace Tags {
/// Two independently evolved amplitudes, stored as [B_plus, B_cross].
/// Empty on elements without an evolved physical boundary. This is evolved
/// state, not a cache in the boundary-condition object.
struct PhysicalTrace : db::SimpleTag {
  using type = DataVector;
};
}  // namespace Tags

/// Receiver geometry and its derivative at fixed surface labels. The normal
/// one-form and tangent seeds may move; their derivatives hold surface labels
/// fixed. Omitting these derivatives selects a fixed coordinate surface.
/// The seeds need not be orthonormal, but must be independent and tangential.
/// Gram-Schmidt and its analytic derivative use the receiver spatial metric.
struct PhysicalTraceGeometry {
  PhysicalTraceGeometry(
      const tnsr::aa<DataVector, 3>& metric,
      const tnsr::aa<DataVector, 3>& dt_metric,
      const tnsr::i<DataVector, 3>& normal_one_form,
      const std::array<tnsr::I<DataVector, 3>, 2>& tangent_seeds,
      const std::optional<tnsr::i<DataVector, 3>>& dt_normal_one_form =
          std::nullopt,
      const std::optional<std::array<tnsr::I<DataVector, 3>, 2>>&
          dt_tangent_seeds = std::nullopt);

  Scalar<DataVector> lapse;
  tnsr::I<DataVector, 3> shift;
  tnsr::AA<DataVector, 3> inverse_metric;
  tnsr::II<DataVector, 3> inverse_spatial_metric;
  tnsr::A<DataVector, 3> time_normal;
  tnsr::i<DataVector, 3> normal_covector;
  tnsr::I<DataVector, 3> normal_vector;
  tnsr::I<DataVector, 3> dt_normal_vector;
  std::array<tnsr::AA<DataVector, 3>, 2> polarization;
  std::array<tnsr::AA<DataVector, 3>, 2> dt_polarization;
  std::array<tnsr::aa<DataVector, 3>, 2> lower_polarization;
};

/// Geometry for a fixed sphere centered at the inertial origin, with the
/// outward normal of the retained EXTERIOR domain. There are no polar grid
/// points on the full-sphere Gauss-Legendre angular mesh.
PhysicalTraceGeometry spherical_trace_geometry(
    const tnsr::aa<DataVector, 3>& metric,
    const tnsr::aa<DataVector, 3>& dt_metric,
    const tnsr::I<DataVector, 3>& coordinates);

/// Geometry at fixed logical angles of a mapped full-sphere face. Coordinate
/// tangents and their time derivatives are differentiated with Spherepack.
/// The resulting dyad follows the grid, including translation and rotation.
PhysicalTraceGeometry mapped_sphere_trace_geometry(
    const tnsr::aa<DataVector, 3>& metric,
    const tnsr::aa<DataVector, 3>& dt_metric,
    const tnsr::I<DataVector, 3>& coordinates,
    const tnsr::I<DataVector, 3>& mesh_velocity, size_t l_max);

DataVector physical_trace_amplitudes(const tnsr::aa<DataVector, 3>& pi,
                                     const tnsr::iaa<DataVector, 3>& phi,
                                     const PhysicalTraceGeometry& geometry);

tnsr::aa<DataVector, 3> physical_trace_tensor(
    const DataVector& amplitudes, const PhysicalTraceGeometry& geometry);

/// Native GH and constraint-adjusted incoming Weyl evaluation of the SAME
/// boundary jet, followed by the physical trace equation. Supply the mixed
/// ghost fields and their tangential derivatives, not the receiver jet.
/// Gauge H and dH must correspond to the gauge prescription at those fields.
/// Geometry transport is supplied from the RECEIVER metric time derivative.
/// Requires gamma1=-1. For native DampedHarmonic, pass the gauge and inertial
/// coordinates; explicit H,dH are ignored and evaluated on the ghost state.
DataVector physical_boundary_rhs(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi,
    const tnsr::iaa<DataVector, 3>& d_metric,
    const tnsr::iaa<DataVector, 3>& d_pi,
    const tnsr::ijaa<DataVector, 3>& d_phi, const Scalar<DataVector>& gamma0,
    const Scalar<DataVector>& gamma2, const tnsr::a<DataVector, 3>& gauge_h,
    const tnsr::ab<DataVector, 3>& d_gauge_h,
    const tnsr::aa<DataVector, 3>& target_weyl,
    const PhysicalTraceGeometry& geometry,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity = std::nullopt,
    const gh::gauges::DampedHarmonic* damped_harmonic = nullptr,
    const tnsr::I<DataVector, 3>* coordinates = nullptr);
}  // namespace gh::worldtube
