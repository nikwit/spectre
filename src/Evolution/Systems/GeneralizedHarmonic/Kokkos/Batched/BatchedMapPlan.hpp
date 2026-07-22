// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Block.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {

/// \brief Classification of a block's grid-to-inertial coordinate map for the
/// batched GH frame-field update.
///
/// The batched executable can evaluate some coordinate maps directly on the
/// device (avoiding the per-stage host map loop and the associated
/// host-to-device transfers of coordinates, mesh velocity, and the composed
/// inverse Jacobian). This enum tags, per block, which device path applies.
///
/// The classification is fixed at initialization and never changes during an
/// evolution, so it can be computed once and reused (and, later, used to group
/// elements by map kind for per-group kernel launches).
enum class BatchedMapKind : int {
  /// No time dependence: the grid-to-inertial map is the identity, so inertial
  /// coordinates equal grid coordinates, the mesh velocity is zero, and the
  /// composed logical-to-inertial inverse Jacobian equals the static
  /// logical-to-grid inverse Jacobian. Evaluated on device.
  Static = 0,
  /// A time-dependent map that does not yet have a device implementation. The
  /// existing host map evaluation path is used for these blocks.
  HostEvaluated = 1,
  /// A rigid rotation about the z-axis through the origin (the map produced by
  /// domain::creators::time_dependence::RotationAboutZAxis in 3D:
  /// ProductOf2Maps<Rotation<2>, Identity<1>>), driven by a single scalar
  /// function of time (the rotation angle). Evaluated on device from the
  /// staged angle and angular velocity; see rigid_rotation_z_frame_fields.
  RigidRotationZ = 2,
  /// A uniform radial scaling about the origin, x^i = a(t) xi^i, driven by a
  /// single scalar function of time a(t). This is
  /// domain::CoordinateMaps::TimeDependent::CubicScale in its linear-scaling
  /// mode (both expansion functions of time equal, e.g. the CubicScale time
  /// dependence with UseLinearScaling). Evaluated on device from the staged
  /// scale factor and its time derivative; see
  /// uniform_radial_scaling_frame_fields.
  UniformRadialScaling = 3,
  /// A general affine map, x^i = A^i_j xi^j + b^i, with affine frame
  /// velocity. The matrix, translation, their time derivatives, and the
  /// inverse Jacobian are inferred and verified from the host map at every
  /// stage, then evaluated for all grid points on the device. This supports
  /// composed rotations (including quaternion-driven rotations), uniform
  /// expansion, and translations without introspecting type-erased maps.
  ComposedAffine = 4,
  /// The nonlinear transition region of a RotScaleTrans map. Rotation,
  /// expansion, and translation parameters are staged exactly from the map
  /// descriptor and evaluated per grid point on the device.
  RotScaleTransTransition = 5,
};

template <size_t Dim>
struct AffineMapData {
  static constexpr size_t matrix_offset = 0;
  static constexpr size_t translation_offset = Dim * Dim;
  static constexpr size_t velocity_matrix_offset = translation_offset + Dim;
  static constexpr size_t velocity_translation_offset =
      velocity_matrix_offset + Dim * Dim;
  static constexpr size_t inverse_jacobian_offset =
      velocity_translation_offset + Dim;
  static constexpr size_t size = inverse_jacobian_offset + Dim * Dim;
};

template <size_t Dim>
struct RotScaleTransMapData {
  static constexpr size_t rotation_offset = 0;
  static constexpr size_t rotation_deriv_offset = rotation_offset + Dim * Dim;
  static constexpr size_t scale_a_offset = rotation_deriv_offset + Dim * Dim;
  static constexpr size_t scale_b_offset = scale_a_offset + 1;
  static constexpr size_t scale_a_deriv_offset = scale_b_offset + 1;
  static constexpr size_t scale_b_deriv_offset = scale_a_deriv_offset + 1;
  static constexpr size_t translation_offset = scale_b_deriv_offset + 1;
  static constexpr size_t translation_deriv_offset = translation_offset + Dim;
  static constexpr size_t inner_radius_offset = translation_deriv_offset + Dim;
  static constexpr size_t outer_radius_offset = inner_radius_offset + 1;
  static constexpr size_t size = outer_radius_offset + 1;
};

KOKKOS_INLINE_FUNCTION void rot_scale_trans_transition_frame_fields(
    const double (&grid_coords)[3], const double* const map_data,
    double (&inertial_coords)[3], double (&mesh_velocity)[3],
    double (&grid_to_inertial_inverse_jacobian)[9]) {
  using offsets = RotScaleTransMapData<3>;
  const double inner_radius = map_data[offsets::inner_radius_offset];
  const double outer_radius = map_data[offsets::outer_radius_offset];
  const double radius =
      sqrt(grid_coords[0] * grid_coords[0] + grid_coords[1] * grid_coords[1] +
           grid_coords[2] * grid_coords[2]);
  const double radius_difference = outer_radius - inner_radius;
  const bool closer_to_outer =
      1.0 - radius / (inner_radius + outer_radius) < 0.5;

  const double scale_a = map_data[offsets::scale_a_offset];
  const double scale_b = map_data[offsets::scale_b_offset];
  const double scale_a_deriv = map_data[offsets::scale_a_deriv_offset];
  const double scale_b_deriv = map_data[offsets::scale_b_deriv_offset];
  const double scale_difference = scale_a - scale_b;
  const double scale_deriv_difference = scale_a_deriv - scale_b_deriv;
  const double radial_scale =
      closer_to_outer ? inner_radius * (outer_radius - radius) *
                            scale_difference / (radius_difference * radius)
                      : outer_radius * (inner_radius - radius) *
                            scale_difference / (radius_difference * radius);
  const double radial_scale_deriv =
      closer_to_outer
          ? inner_radius * (outer_radius - radius) * scale_deriv_difference /
                (radius_difference * radius)
          : outer_radius * (inner_radius - radius) * scale_deriv_difference /
                (radius_difference * radius);
  const double scale = (closer_to_outer ? scale_b : scale_a) + radial_scale;
  const double scale_deriv =
      (closer_to_outer ? scale_b_deriv : scale_a_deriv) + radial_scale_deriv;
  const double translation_weight =
      closer_to_outer ? (outer_radius - radius) / radius_difference
                      : 1.0 + (inner_radius - radius) / radius_difference;

  double rotated_coords[3]{0.0, 0.0, 0.0};
  double rotated_velocity[3]{0.0, 0.0, 0.0};
  for (size_t target_d = 0; target_d < 3; ++target_d) {
    for (size_t source_d = 0; source_d < 3; ++source_d) {
      rotated_coords[target_d] +=
          map_data[offsets::rotation_offset + target_d * 3 + source_d] *
          grid_coords[source_d];
      rotated_velocity[target_d] +=
          map_data[offsets::rotation_deriv_offset + target_d * 3 + source_d] *
          grid_coords[source_d];
    }
    inertial_coords[target_d] =
        scale * rotated_coords[target_d] +
        translation_weight * map_data[offsets::translation_offset + target_d];
    mesh_velocity[target_d] =
        scale * rotated_velocity[target_d] +
        scale_deriv * rotated_coords[target_d] +
        translation_weight *
            map_data[offsets::translation_deriv_offset + target_d];
  }

  const double scale_radial_deriv = -inner_radius * outer_radius *
                                    scale_difference /
                                    (radius_difference * radius * radius);
  const double translation_radial_deriv = -1.0 / radius_difference;
  double jacobian[9]{};
  for (size_t target_d = 0; target_d < 3; ++target_d) {
    for (size_t source_d = 0; source_d < 3; ++source_d) {
      jacobian[target_d * 3 + source_d] =
          scale * map_data[offsets::rotation_offset + target_d * 3 + source_d] +
          (scale_radial_deriv * rotated_coords[target_d] +
           translation_radial_deriv *
               map_data[offsets::translation_offset + target_d]) *
              grid_coords[source_d] / radius;
    }
  }

  const double determinant =
      jacobian[0] * (jacobian[4] * jacobian[8] - jacobian[5] * jacobian[7]) -
      jacobian[1] * (jacobian[3] * jacobian[8] - jacobian[5] * jacobian[6]) +
      jacobian[2] * (jacobian[3] * jacobian[7] - jacobian[4] * jacobian[6]);
  grid_to_inertial_inverse_jacobian[0] =
      (jacobian[4] * jacobian[8] - jacobian[5] * jacobian[7]) / determinant;
  grid_to_inertial_inverse_jacobian[1] =
      (jacobian[2] * jacobian[7] - jacobian[1] * jacobian[8]) / determinant;
  grid_to_inertial_inverse_jacobian[2] =
      (jacobian[1] * jacobian[5] - jacobian[2] * jacobian[4]) / determinant;
  grid_to_inertial_inverse_jacobian[3] =
      (jacobian[5] * jacobian[6] - jacobian[3] * jacobian[8]) / determinant;
  grid_to_inertial_inverse_jacobian[4] =
      (jacobian[0] * jacobian[8] - jacobian[2] * jacobian[6]) / determinant;
  grid_to_inertial_inverse_jacobian[5] =
      (jacobian[2] * jacobian[3] - jacobian[0] * jacobian[5]) / determinant;
  grid_to_inertial_inverse_jacobian[6] =
      (jacobian[3] * jacobian[7] - jacobian[4] * jacobian[6]) / determinant;
  grid_to_inertial_inverse_jacobian[7] =
      (jacobian[1] * jacobian[6] - jacobian[0] * jacobian[7]) / determinant;
  grid_to_inertial_inverse_jacobian[8] =
      (jacobian[0] * jacobian[4] - jacobian[1] * jacobian[3]) / determinant;
}

template <size_t Dim>
KOKKOS_INLINE_FUNCTION void composed_affine_frame_fields(
    const double (&grid_coords)[Dim], const double* const affine_data,
    double (&inertial_coords)[Dim], double (&mesh_velocity)[Dim],
    double (&grid_to_inertial_inverse_jacobian)[Dim * Dim]) {
  for (size_t inertial_d = 0; inertial_d < Dim; ++inertial_d) {
    inertial_coords[inertial_d] =
        affine_data[AffineMapData<Dim>::translation_offset + inertial_d];
    mesh_velocity[inertial_d] =
        affine_data[AffineMapData<Dim>::velocity_translation_offset +
                    inertial_d];
    for (size_t grid_d = 0; grid_d < Dim; ++grid_d) {
      inertial_coords[inertial_d] +=
          affine_data[AffineMapData<Dim>::matrix_offset + inertial_d * Dim +
                      grid_d] *
          grid_coords[grid_d];
      mesh_velocity[inertial_d] +=
          affine_data[AffineMapData<Dim>::velocity_matrix_offset +
                      inertial_d * Dim + grid_d] *
          grid_coords[grid_d];
    }
  }
  for (size_t component = 0; component < Dim * Dim; ++component) {
    grid_to_inertial_inverse_jacobian[component] =
        affine_data[AffineMapData<Dim>::inverse_jacobian_offset + component];
  }
}

/// \brief Closed-form grid-to-inertial frame fields for a rigid rotation about
/// the z-axis through the origin.
///
/// Reproduces domain::CoordinateMaps::TimeDependent::Rotation<2> (composed with
/// an identity in z) for a single grid point, given the rotation angle
/// \p theta and angular velocity \p theta_dot at the current time. Callable on
/// both host (for the initialization-time verification against the host map)
/// and device (in the batched frame-field kernel).
///
/// Outputs (all in the inertial frame, raw arrays to avoid device-side
/// Tensor/DataVector allocation):
/// - \p inertial_coords : the mapped point, inertial = R(theta) . grid.
/// - \p mesh_velocity : d/dt of the mapped point at fixed grid coordinates,
///   theta_dot * R'(theta) . grid.
/// - \p grid_to_inertial_inverse_jacobian : row-major [grid * 3 + inertial],
///   d(grid)/d(inertial) = R(theta)^T (with identity in z), matching the
///   indexing consumed when composing with the static logical-to-grid inverse
///   Jacobian.
KOKKOS_INLINE_FUNCTION void rigid_rotation_z_frame_fields(
    const double grid_x, const double grid_y, const double grid_z,
    const double theta, const double theta_dot, double (&inertial_coords)[3],
    double (&mesh_velocity)[3],
    double (&grid_to_inertial_inverse_jacobian)[9]) {
  const double cos_theta = cos(theta);
  const double sin_theta = sin(theta);

  inertial_coords[0] = cos_theta * grid_x - sin_theta * grid_y;
  inertial_coords[1] = sin_theta * grid_x + cos_theta * grid_y;
  inertial_coords[2] = grid_z;

  mesh_velocity[0] = theta_dot * (-sin_theta * grid_x - cos_theta * grid_y);
  mesh_velocity[1] = theta_dot * (cos_theta * grid_x - sin_theta * grid_y);
  mesh_velocity[2] = 0.0;

  // d(grid_a)/d(inertial_b) = (R^T)_{a,b}, row-major [a * 3 + b].
  grid_to_inertial_inverse_jacobian[0] = cos_theta;   // (x_grid, x_inertial)
  grid_to_inertial_inverse_jacobian[1] = sin_theta;   // (x_grid, y_inertial)
  grid_to_inertial_inverse_jacobian[2] = 0.0;         // (x_grid, z_inertial)
  grid_to_inertial_inverse_jacobian[3] = -sin_theta;  // (y_grid, x_inertial)
  grid_to_inertial_inverse_jacobian[4] = cos_theta;   // (y_grid, y_inertial)
  grid_to_inertial_inverse_jacobian[5] = 0.0;         // (y_grid, z_inertial)
  grid_to_inertial_inverse_jacobian[6] = 0.0;         // (z_grid, x_inertial)
  grid_to_inertial_inverse_jacobian[7] = 0.0;         // (z_grid, y_inertial)
  grid_to_inertial_inverse_jacobian[8] = 1.0;         // (z_grid, z_inertial)
}

/// \brief Closed-form grid-to-inertial frame fields for a uniform radial
/// scaling about the origin, inertial = scale * grid.
///
/// Reproduces domain::CoordinateMaps::TimeDependent::CubicScale in its
/// linear-scaling mode for a single grid point, given the scale factor
/// \p scale and its time derivative \p scale_dot at the current time. Callable
/// on both host (verification) and device (kernel). Outputs use the same
/// conventions and row-major inverse-Jacobian layout as
/// rigid_rotation_z_frame_fields; here d(grid_a)/d(inertial_b) =
/// (1 / scale) delta_{ab}.
KOKKOS_INLINE_FUNCTION void uniform_radial_scaling_frame_fields(
    const double grid_x, const double grid_y, const double grid_z,
    const double scale, const double scale_dot, double (&inertial_coords)[3],
    double (&mesh_velocity)[3],
    double (&grid_to_inertial_inverse_jacobian)[9]) {
  inertial_coords[0] = scale * grid_x;
  inertial_coords[1] = scale * grid_y;
  inertial_coords[2] = scale * grid_z;

  mesh_velocity[0] = scale_dot * grid_x;
  mesh_velocity[1] = scale_dot * grid_y;
  mesh_velocity[2] = scale_dot * grid_z;

  const double inverse_scale = 1.0 / scale;
  for (size_t i = 0; i < 9; ++i) {
    grid_to_inertial_inverse_jacobian[i] = 0.0;
  }
  grid_to_inertial_inverse_jacobian[0] = inverse_scale;  // (x_grid, x_inertial)
  grid_to_inertial_inverse_jacobian[4] = inverse_scale;  // (y_grid, y_inertial)
  grid_to_inertial_inverse_jacobian[8] = inverse_scale;  // (z_grid, z_inertial)
}

/// \brief Classify a block's grid-to-inertial map for the batched device path.
///
/// Only static (non-time-dependent) blocks are currently evaluated on the
/// device. Every time-dependent block falls back to the host map path until
/// its map kind is ported (see the device-map-port milestones).
template <size_t Dim>
BatchedMapKind classify_block_map(const Block<Dim>& block) {
  return block.is_time_dependent() ? BatchedMapKind::HostEvaluated
                                   : BatchedMapKind::Static;
}

}  // namespace gh::Actions
