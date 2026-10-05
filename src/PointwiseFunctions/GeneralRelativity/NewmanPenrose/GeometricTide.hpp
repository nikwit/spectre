// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <utility>

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

/// Angular identification and temporal data shared by both tidal orders.
/// screen uses unit coordinate-radius tangents in the adapted NR tetrad.
/// observer is the invariant rest observer in that same tetrad; boost is
/// the total type-III rapidity relative to registration.pulled_back.
struct GeometricFrame {
  EigenSphereMap map;
  std::array<ComplexDataVector, 3> dyad;
  std::array<AdaptedFourVector, 2> screen;
  AdaptedFourVector observer;
  DataVector boost;
  double maximum_dyad_error{};
};

GeometricFrame geometric_frame(
    const FrameRegistration& registration, const Scalar<DataVector>& rapidity,
    const RealMatrix& adapted_rotation,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const TriadVector& label_directions, const DataVector& label_weights,
    double mass, size_t l_max);

/// Least-squares potential of a tangential covector on the label sphere.
/// The constant is fixed by the weighted mean; residual measures nonclosure
/// and truncation. The two gradient components use dTheta, dPhi/sinTheta.
std::pair<DataVector, double> sphere_gradient_potential(
    const TriadVector& labels, const DataVector& weights,
    const std::array<DataVector, 2>& gradient, size_t l_max);

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
