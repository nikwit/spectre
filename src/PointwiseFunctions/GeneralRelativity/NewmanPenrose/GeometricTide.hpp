// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <cstddef>
#include <optional>

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Psi4Fit.hpp"

namespace gr::np {
/// A sphere map from the first three nonconstant screen-Laplacian modes.
/// Derivatives use the orthonormal NR label basis (dTheta, dPhi/sinTheta).
struct EigenSphereMap {
  TriadVector direction;
  std::array<TriadVector, 2> derivative;
  DataVector weights;
  std::array<double, 5> eigenvalues{};
  double minimum_jacobian{};
  double area_over_label_area{};
};

/// Dense symmetric generalized eigenproblem in real scalar harmonics.
/// h = (h00,h01,h11) in the orthonormal label basis. A constant rescaling
/// of h does not change the map. Eigenvalues are reported after normalizing
/// its area to 4 pi. The returned weights include the map's area Jacobian.
EigenSphereMap laplace_eigenmap(const TriadVector& label_directions,
                                const DataVector& label_weights,
                                const std::array<DataVector, 3>& h,
                                size_t l_max);

struct GeometricTideEvaluation {
  SecondOrderEvaluation second_order;
  EigenSphereMap map;
  double maximum_dyad_error{};
};

/// Geometric quadrupole on a round coordinate sphere. label_directions
/// point away from its center, unlike the domain-outward inner normal.
/// Retains the legacy leading-frame fit, radius, boost, and incoming-slot
/// replacement; only the common angular axes, dyad transport and area
/// weights change. The coordinate radius cancels from both map and polar
/// factor, so unit-radius coordinate tangents suffice.
GeometricTideEvaluation evaluate_geometric_second_order(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const TriadVector& label_directions, const DataVector& label_weights,
    double mass, size_t l_max,
    const std::optional<TidalMoments>& imposed_components = std::nullopt);
}  // namespace gr::np
