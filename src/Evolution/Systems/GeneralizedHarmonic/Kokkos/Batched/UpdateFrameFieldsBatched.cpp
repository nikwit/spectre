// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <variant>
#include <vector>

#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/TimeDependent/RotationMatrixHelpers.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/BatchedMapPlan.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions {
namespace {

template <size_t VolumeDim>
void resize_host_array_if_needed(
    const gsl::not_null<std::array<DataVector, VolumeDim>*> host_data,
    const size_t total_points) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    if (gsl::at(*host_data, d).size() != total_points) {
      gsl::at(*host_data, d) = DataVector{total_points, 0.0};
    }
  }
}

template <size_t VolumeDim>
void ensure_frame_field_storage(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t num_elements, const size_t points_per_element,
    const size_t total_points) {
  using packed_geometry_type =
      evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>;

  if (packed_geometry->element_inverse_jacobian_device.extent(0) !=
          num_elements or
      packed_geometry->element_inverse_jacobian_device.extent(1) !=
          points_per_element or
      packed_geometry->element_inverse_jacobian_device.extent(2) !=
          VolumeDim * VolumeDim) {
    packed_geometry->element_inverse_jacobian_device =
        typename packed_geometry_type::device_inverse_jacobian_type{
            "GhBatchedElementInverseJacobian", num_elements, points_per_element,
            VolumeDim * VolumeDim};
  }
  const auto ensure_host_jacobian =
      [num_elements, points_per_element](
          const gsl::not_null<
              typename packed_geometry_type::host_inverse_jacobian_cache_type*>
              jacobian,
          const std::string& label) {
        if (jacobian->extent(0) != num_elements or
            jacobian->extent(1) != points_per_element or
            jacobian->extent(2) != VolumeDim * VolumeDim) {
          *jacobian =
              typename packed_geometry_type::host_inverse_jacobian_cache_type{
                  label, num_elements, points_per_element,
                  VolumeDim * VolumeDim};
        }
      };
  ensure_host_jacobian(
      make_not_null(
          &packed_geometry->logical_to_distorted_inverse_jacobian_host),
      "GhBatchedLogicalToDistortedInverseJacobian");
  ensure_host_jacobian(
      make_not_null(&packed_geometry->grid_to_distorted_inverse_jacobian_host),
      "GhBatchedGridToDistortedInverseJacobian");
  ensure_host_jacobian(
      make_not_null(&packed_geometry->distorted_to_inertial_jacobian_host),
      "GhBatchedDistortedToInertialJacobian");
  ensure_host_jacobian(
      make_not_null(
          &packed_geometry->distorted_to_inertial_inverse_jacobian_host),
      "GhBatchedDistortedToInertialInverseJacobian");

  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->distorted_coordinates_host),
      total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->inertial_coordinates_host), total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->mesh_velocity_host), total_points);
  resize_host_array_if_needed<VolumeDim>(
      make_not_null(&packed_geometry->grid_to_distorted_frame_velocity_host),
      total_points);

  for (size_t d = 0; d < VolumeDim; ++d) {
    if (packed_geometry->inertial_coordinates_device.get(d).extent(0) !=
        total_points) {
      packed_geometry->inertial_coordinates_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedInertialCoordinates", total_points);
    }
    if (packed_geometry->mesh_velocity_device.get(d).extent(0) !=
        total_points) {
      packed_geometry->mesh_velocity_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedMeshVelocity", total_points);
    }
  }
}

template <size_t VolumeDim>
void copy_host_array_to_device(
    const std::array<DataVector, VolumeDim>& host_data,
    const gsl::not_null<typename evolution::Kokkos::PackedGeometry<
        gh::System<VolumeDim>>::device_inertial_coordinates_type*>
        device_data) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    auto host_component = ::Kokkos::create_mirror_view(device_data->get(d));
    for (size_t s = 0; s < host_data[d].size(); ++s) {
      host_component(s) = host_data[d][s];
    }
    ::Kokkos::deep_copy(device_data->get(d), host_component);
  }
}

// Upload the domain-static geometry (grid coordinates and the logical-to-grid
// inverse Jacobian) to the device once. These never change during an
// evolution, so this is a one-time cost guarded by the current device extents.
template <size_t VolumeDim>
void ensure_static_device_geometry(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t num_elements, const size_t points_per_element,
    const size_t total_points) {
  using packed_geometry_type =
      evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>;

  bool coordinates_ready = true;
  for (size_t d = 0; d < VolumeDim; ++d) {
    if (packed_geometry->grid_coordinates_device.get(d).extent(0) !=
        total_points) {
      coordinates_ready = false;
    }
  }
  if (not coordinates_ready) {
    for (size_t d = 0; d < VolumeDim; ++d) {
      packed_geometry->grid_coordinates_device.get(d) =
          ::Kokkos::View<double*>("GhBatchedGridCoordinates", total_points);
    }
    copy_host_array_to_device<VolumeDim>(
        packed_geometry->grid_coordinates_host,
        make_not_null(&packed_geometry->grid_coordinates_device));
  }

  if (packed_geometry->logical_to_grid_inverse_jacobian_device.extent(0) !=
          num_elements or
      packed_geometry->logical_to_grid_inverse_jacobian_device.extent(1) !=
          points_per_element or
      packed_geometry->logical_to_grid_inverse_jacobian_device.extent(2) !=
          VolumeDim * VolumeDim) {
    packed_geometry->logical_to_grid_inverse_jacobian_device =
        typename packed_geometry_type::device_inverse_jacobian_type{
            "GhBatchedLogicalToGridInverseJacobianDevice", num_elements,
            points_per_element, VolumeDim * VolumeDim};
    ::Kokkos::deep_copy(
        packed_geometry->logical_to_grid_inverse_jacobian_device,
        packed_geometry->logical_to_grid_inverse_jacobian_host);
  }
}

// Device frame-field update for a domain whose every block is evaluated on the
// device (a mix of BatchedMapKind::Static and BatchedMapKind::RigidRotationZ).
// Writes the device-resident geometry consumed by the RHS/boundary kernels
// directly on the device, dispatching per element on the map kind. `theta` and
// `theta_dot` are the shared rigid-rotation angle and angular velocity at the
// current time (unused for static-only domains). This is the enum-dispatch form
// of the device map path; the per-element kind is also the natural key for a
// later group-by-kind launch optimization.
template <size_t VolumeDim>
void update_frame_fields_device(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const size_t points_per_element, const size_t total_points,
    const double fot_value, const double fot_value_dot) {
  const auto grid_coordinates = packed_geometry->grid_coordinates_device;
  const auto inertial_coordinates =
      packed_geometry->inertial_coordinates_device;
  const auto mesh_velocity = packed_geometry->mesh_velocity_device;
  const auto logical_to_grid_inverse_jacobian =
      packed_geometry->logical_to_grid_inverse_jacobian_device;
  const auto element_inverse_jacobian =
      packed_geometry->element_inverse_jacobian_device;
  const auto element_map_kind = packed_geometry->element_map_kind_device;
  const auto affine_map_data = packed_geometry->affine_map_data_device;
  const size_t ppe = points_per_element;
  constexpr int static_kind = static_cast<int>(BatchedMapKind::Static);
  constexpr int host_evaluated_kind =
      static_cast<int>(BatchedMapKind::HostEvaluated);
  constexpr int rigid_rotation_z_kind =
      static_cast<int>(BatchedMapKind::RigidRotationZ);
  constexpr int uniform_radial_scaling_kind =
      static_cast<int>(BatchedMapKind::UniformRadialScaling);
  constexpr int composed_affine_kind =
      static_cast<int>(BatchedMapKind::ComposedAffine);
  constexpr int rot_scale_trans_transition_kind =
      static_cast<int>(BatchedMapKind::RotScaleTransTransition);
  ::Kokkos::parallel_for(
      "GhUpdateFrameFields", total_points,
      KOKKOS_LAMBDA(const int linear_index_int) {
        const size_t s = static_cast<size_t>(linear_index_int);
        const size_t e = s / ppe;
        const size_t p = s % ppe;
        const int kind = element_map_kind(e);
        if (kind == host_evaluated_kind) {
          // These entries were copied from the host before a mixed-map launch.
          return;
        } else if (kind == static_kind) {
          for (size_t d = 0; d < VolumeDim; ++d) {
            inertial_coordinates.get(d)(s) = grid_coordinates.get(d)(s);
            mesh_velocity.get(d)(s) = 0.0;
          }
          for (size_t component = 0; component < VolumeDim * VolumeDim;
               ++component) {
            element_inverse_jacobian(e, p, component) =
                logical_to_grid_inverse_jacobian(e, p, component);
          }
        } else {
          double inertial[3];
          double velocity[3];
          double grid_to_inertial_inv_jac[9];
          if (kind == composed_affine_kind) {
            const double grid_coords[3]{grid_coordinates.get(0)(s),
                                        grid_coordinates.get(1)(s),
                                        grid_coordinates.get(2)(s)};
            composed_affine_frame_fields<VolumeDim>(
                grid_coords, &affine_map_data(e, 0), inertial, velocity,
                grid_to_inertial_inv_jac);
          } else if (kind == rot_scale_trans_transition_kind) {
            const double grid_coords[3]{grid_coordinates.get(0)(s),
                                        grid_coordinates.get(1)(s),
                                        grid_coordinates.get(2)(s)};
            rot_scale_trans_transition_frame_fields(
                grid_coords, &affine_map_data(e, 0), inertial, velocity,
                grid_to_inertial_inv_jac);
          } else if (kind == rigid_rotation_z_kind) {
            rigid_rotation_z_frame_fields(
                grid_coordinates.get(0)(s), grid_coordinates.get(1)(s),
                grid_coordinates.get(2)(s), fot_value, fot_value_dot, inertial,
                velocity, grid_to_inertial_inv_jac);
          } else if (kind == uniform_radial_scaling_kind) {
            uniform_radial_scaling_frame_fields(
                grid_coordinates.get(0)(s), grid_coordinates.get(1)(s),
                grid_coordinates.get(2)(s), fot_value, fot_value_dot, inertial,
                velocity, grid_to_inertial_inv_jac);
          } else {
            Kokkos::abort("Unsupported device map kind.");
          }
          for (size_t d = 0; d < VolumeDim; ++d) {
            inertial_coordinates.get(d)(s) = inertial[d];
            mesh_velocity.get(d)(s) = velocity[d];
          }
          // Compose logical->grid (static) with grid->inertial (device map):
          // invjac(logical, inertial) = sum_grid L2G(logical, grid) *
          //                             G2I(grid, inertial).
          for (size_t logical_d = 0; logical_d < VolumeDim; ++logical_d) {
            for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
              double composed = 0.0;
              for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
                composed +=
                    logical_to_grid_inverse_jacobian(
                        e, p, logical_d * VolumeDim + grid_d) *
                    grid_to_inertial_inv_jac[grid_d * VolumeDim + inertial_d];
              }
              element_inverse_jacobian(
                  e, p, logical_d * VolumeDim + inertial_d) = composed;
            }
          }
        }
      });
}

// Copy a device coordinate/velocity array back into its host mirror. Used to
// keep the host geometry arrays current for the host-side observation and
// apparent-horizon consumers when the device path evaluates a moving map.
template <size_t VolumeDim>
void copy_device_array_to_host(
    const typename evolution::Kokkos::PackedGeometry<
        gh::System<VolumeDim>>::device_inertial_coordinates_type& device_data,
    const gsl::not_null<std::array<DataVector, VolumeDim>*> host_data) {
  for (size_t d = 0; d < VolumeDim; ++d) {
    const auto mirror = ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace{}, device_data.get(d));
    const size_t num_points = device_data.get(d).extent(0);
    if (gsl::at(*host_data, d).size() != num_points) {
      gsl::at(*host_data, d) = DataVector{num_points, 0.0};
    }
    for (size_t s = 0; s < num_points; ++s) {
      gsl::at(*host_data, d)[s] = mirror(s);
    }
  }
}

// Verify that a block's grid-to-inertial map is, at the given time, the rigid
// z-rotation closed form used by the device path (rigid_rotation_z_frame_fields
// with the angle and angular velocity from `fot_name`). Compares the host map's
// coordinates, mesh velocity, and grid-to-inertial inverse Jacobian at the
// element's grid points. This is the safety gate for the name-based
// recognizer: only if the host map provably matches the device closed form is
// the block evaluated on the device; otherwise it falls back to the host path.
template <size_t VolumeDim>
bool verify_rigid_rotation_z(
    const Block<VolumeDim>& block,
    const tnsr::I<DataVector, VolumeDim, Frame::Grid>& grid_coords,
    const std::string& fot_name, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  const auto map_data =
      block.moving_mesh_grid_to_inertial_map().coords_frame_velocity_jacobians(
          grid_coords, time, functions_of_time);
  const auto& inertial_coords = std::get<0>(map_data);
  const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
  const auto& mesh_velocity = std::get<3>(map_data);

  const auto angle_and_deriv =
      functions_of_time.at(fot_name)->func_and_deriv(time);
  if (angle_and_deriv[0].size() != 1 or angle_and_deriv[1].size() != 1) {
    return false;
  }
  const double theta = angle_and_deriv[0][0];
  const double theta_dot = angle_and_deriv[1][0];

  const double tolerance = 1.0e-11;
  const auto within_tolerance = [tolerance](const double a, const double b) {
    return std::abs(a - b) <= tolerance * (1.0 + std::abs(b));
  };

  const size_t num_points = grid_coords.get(0).size();
  for (size_t p = 0; p < num_points; ++p) {
    double inertial[3];
    double velocity[3];
    double grid_to_inertial_inv_jac[9];
    rigid_rotation_z_frame_fields(grid_coords.get(0)[p], grid_coords.get(1)[p],
                                  grid_coords.get(2)[p], theta, theta_dot,
                                  inertial, velocity, grid_to_inertial_inv_jac);
    for (size_t d = 0; d < VolumeDim; ++d) {
      if (not within_tolerance(inertial_coords.get(d)[p], inertial[d])) {
        return false;
      }
      if (not within_tolerance(mesh_velocity.get(d)[p], velocity[d])) {
        return false;
      }
    }
    for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
      for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
        if (not within_tolerance(
                grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p],
                grid_to_inertial_inv_jac[grid_d * VolumeDim + inertial_d])) {
          return false;
        }
      }
    }
  }
  return true;
}

// Verify that a block's grid-to-inertial map is, at the given time, the uniform
// radial scaling closed form used by the device path
// (uniform_radial_scaling_frame_fields with the scale factor and its derivative
// from `fot_name`). Same safety-gate role as verify_rigid_rotation_z: a
// non-scaling map (e.g. a rotation, or a non-uniform CubicScale) fails the
// comparison and falls back to the host path.
template <size_t VolumeDim>
bool verify_uniform_radial_scaling(
    const Block<VolumeDim>& block,
    const tnsr::I<DataVector, VolumeDim, Frame::Grid>& grid_coords,
    const std::string& fot_name, const double time,
    const domain::FunctionsOfTimeMap& functions_of_time) {
  const auto map_data =
      block.moving_mesh_grid_to_inertial_map().coords_frame_velocity_jacobians(
          grid_coords, time, functions_of_time);
  const auto& inertial_coords = std::get<0>(map_data);
  const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
  const auto& mesh_velocity = std::get<3>(map_data);

  const auto scale_and_deriv =
      functions_of_time.at(fot_name)->func_and_deriv(time);
  if (scale_and_deriv[0].size() != 1 or scale_and_deriv[1].size() != 1) {
    return false;
  }
  const double scale = scale_and_deriv[0][0];
  const double scale_dot = scale_and_deriv[1][0];
  if (scale <= 0.0) {
    return false;
  }

  const double tolerance = 1.0e-11;
  const auto within_tolerance = [tolerance](const double a, const double b) {
    return std::abs(a - b) <= tolerance * (1.0 + std::abs(b));
  };

  const size_t num_points = grid_coords.get(0).size();
  for (size_t p = 0; p < num_points; ++p) {
    double inertial[3];
    double velocity[3];
    double grid_to_inertial_inv_jac[9];
    uniform_radial_scaling_frame_fields(
        grid_coords.get(0)[p], grid_coords.get(1)[p], grid_coords.get(2)[p],
        scale, scale_dot, inertial, velocity, grid_to_inertial_inv_jac);
    for (size_t d = 0; d < VolumeDim; ++d) {
      if (not within_tolerance(inertial_coords.get(d)[p], inertial[d])) {
        return false;
      }
      if (not within_tolerance(mesh_velocity.get(d)[p], velocity[d])) {
        return false;
      }
    }
    for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
      for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
        if (not within_tolerance(
                grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p],
                grid_to_inertial_inv_jac[grid_d * VolumeDim + inertial_d])) {
          return false;
        }
      }
    }
  }
  return true;
}

// Infer a general affine grid-to-inertial map from four well-conditioned
// points in a representative element, then verify it against every element
// point. Restricting probes to the block avoids singular evaluations of radial
// maps at artificial points such as the origin. This is repeated at every
// stage so a controlled map that ceases to be affine safely returns to the
// host path.
template <size_t VolumeDim>
bool infer_composed_affine_map(
    const gsl::not_null<std::array<double, AffineMapData<VolumeDim>::size>*>
        affine_data,
    const Block<VolumeDim>& block,
    const tnsr::I<DataVector, VolumeDim, Frame::Grid>& representative_coords,
    const double time, const domain::FunctionsOfTimeMap& functions_of_time) {
  static_assert(VolumeDim == 3,
                "The affine map inference currently supports only 3D.");
  const size_t representative_points = representative_coords.get(0).size();
  if (representative_points < VolumeDim + 1) {
    return false;
  }

  std::array<size_t, VolumeDim + 1> fit_points{{0, 0, 0, 0}};
  double longest_displacement_squared = 0.0;
  for (size_t p = 1; p < representative_points; ++p) {
    double displacement_squared = 0.0;
    for (size_t d = 0; d < VolumeDim; ++d) {
      const double displacement =
          representative_coords.get(d)[p] - representative_coords.get(d)[0];
      displacement_squared += displacement * displacement;
    }
    if (displacement_squared > longest_displacement_squared) {
      longest_displacement_squared = displacement_squared;
      fit_points[1] = p;
    }
  }
  if (longest_displacement_squared == 0.0) {
    return false;
  }

  const std::array<double, 3> first_displacement{
      representative_coords.get(0)[fit_points[1]] -
          representative_coords.get(0)[0],
      representative_coords.get(1)[fit_points[1]] -
          representative_coords.get(1)[0],
      representative_coords.get(2)[fit_points[1]] -
          representative_coords.get(2)[0]};
  double largest_cross_product_squared = 0.0;
  std::array<double, 3> plane_normal{};
  for (size_t p = 1; p < representative_points; ++p) {
    const std::array<double, 3> displacement{
        representative_coords.get(0)[p] - representative_coords.get(0)[0],
        representative_coords.get(1)[p] - representative_coords.get(1)[0],
        representative_coords.get(2)[p] - representative_coords.get(2)[0]};
    const std::array<double, 3> cross_product{
        first_displacement[1] * displacement[2] -
            first_displacement[2] * displacement[1],
        first_displacement[2] * displacement[0] -
            first_displacement[0] * displacement[2],
        first_displacement[0] * displacement[1] -
            first_displacement[1] * displacement[0]};
    const double cross_product_squared = cross_product[0] * cross_product[0] +
                                         cross_product[1] * cross_product[1] +
                                         cross_product[2] * cross_product[2];
    if (cross_product_squared > largest_cross_product_squared) {
      largest_cross_product_squared = cross_product_squared;
      fit_points[2] = p;
      plane_normal = cross_product;
    }
  }
  if (largest_cross_product_squared == 0.0) {
    return false;
  }

  double largest_volume = 0.0;
  for (size_t p = 1; p < representative_points; ++p) {
    double signed_volume = 0.0;
    for (size_t d = 0; d < VolumeDim; ++d) {
      signed_volume += plane_normal[d] * (representative_coords.get(d)[p] -
                                          representative_coords.get(d)[0]);
    }
    if (std::abs(signed_volume) > largest_volume) {
      largest_volume = std::abs(signed_volume);
      fit_points[3] = p;
    }
  }
  if (largest_volume == 0.0) {
    return false;
  }

  const auto map_data =
      block.moving_mesh_grid_to_inertial_map().coords_frame_velocity_jacobians(
          representative_coords, time, functions_of_time);
  const auto& inertial_coords = std::get<0>(map_data);
  const auto& inverse_jacobian = std::get<1>(map_data);
  const auto& mesh_velocity = std::get<3>(map_data);

  // Columns of point_displacements span the representative element. Its
  // inverse converts differences at the three fit points into derivatives
  // with respect to the grid coordinates.
  std::array<double, 9> point_displacements{};
  for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
    for (size_t column = 0; column < VolumeDim; ++column) {
      point_displacements[grid_d * VolumeDim + column] =
          representative_coords.get(grid_d)[fit_points[column + 1]] -
          representative_coords.get(grid_d)[fit_points[0]];
    }
  }
  const double a = point_displacements[0];
  const double b = point_displacements[1];
  const double c = point_displacements[2];
  const double d = point_displacements[3];
  const double e = point_displacements[4];
  const double f = point_displacements[5];
  const double g = point_displacements[6];
  const double h = point_displacements[7];
  const double i = point_displacements[8];
  const double determinant =
      a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  const double coordinate_scale = std::sqrt(longest_displacement_squared);
  if (not std::isfinite(determinant) or
      std::abs(determinant) <=
          1.0e-13 * coordinate_scale * coordinate_scale * coordinate_scale) {
    return false;
  }
  const std::array<double, 9> inverse_point_displacements{
      (e * i - f * h) / determinant, (c * h - b * i) / determinant,
      (b * f - c * e) / determinant, (f * g - d * i) / determinant,
      (a * i - c * g) / determinant, (c * d - a * f) / determinant,
      (d * h - e * g) / determinant, (b * g - a * h) / determinant,
      (a * e - b * d) / determinant};

  for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
    for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
      double matrix_component = 0.0;
      double velocity_matrix_component = 0.0;
      for (size_t column = 0; column < VolumeDim; ++column) {
        matrix_component +=
            (inertial_coords.get(inertial_d)[fit_points[column + 1]] -
             inertial_coords.get(inertial_d)[fit_points[0]]) *
            inverse_point_displacements[column * VolumeDim + grid_d];
        velocity_matrix_component +=
            (mesh_velocity.get(inertial_d)[fit_points[column + 1]] -
             mesh_velocity.get(inertial_d)[fit_points[0]]) *
            inverse_point_displacements[column * VolumeDim + grid_d];
      }
      (*affine_data)[AffineMapData<VolumeDim>::matrix_offset +
                     inertial_d * VolumeDim + grid_d] = matrix_component;
      (*affine_data)[AffineMapData<VolumeDim>::velocity_matrix_offset +
                     inertial_d * VolumeDim + grid_d] =
          velocity_matrix_component;
      (*affine_data)[AffineMapData<VolumeDim>::inverse_jacobian_offset +
                     grid_d * VolumeDim + inertial_d] =
          inverse_jacobian.get(grid_d, inertial_d)[fit_points[0]];
    }
    double translation = inertial_coords.get(inertial_d)[fit_points[0]];
    double velocity_translation = mesh_velocity.get(inertial_d)[fit_points[0]];
    for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
      translation -= (*affine_data)[AffineMapData<VolumeDim>::matrix_offset +
                                    inertial_d * VolumeDim + grid_d] *
                     representative_coords.get(grid_d)[fit_points[0]];
      velocity_translation -=
          (*affine_data)[AffineMapData<VolumeDim>::velocity_matrix_offset +
                         inertial_d * VolumeDim + grid_d] *
          representative_coords.get(grid_d)[fit_points[0]];
    }
    (*affine_data)[AffineMapData<VolumeDim>::translation_offset + inertial_d] =
        translation;
    (*affine_data)[AffineMapData<VolumeDim>::velocity_translation_offset +
                   inertial_d] = velocity_translation;
  }

  constexpr double tolerance = 1.0e-11;
  const auto within_tolerance = [tolerance](const double lhs,
                                            const double rhs) {
    return std::isfinite(lhs) and std::isfinite(rhs) and
           std::abs(lhs - rhs) <= tolerance * (1.0 + std::abs(rhs));
  };
  for (size_t p = 0; p < representative_points; ++p) {
    for (size_t inertial_d = 0; inertial_d < VolumeDim; ++inertial_d) {
      double expected_coord =
          (*affine_data)[AffineMapData<VolumeDim>::translation_offset +
                         inertial_d];
      double expected_velocity =
          (*affine_data)[AffineMapData<VolumeDim>::velocity_translation_offset +
                         inertial_d];
      for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
        expected_coord +=
            (*affine_data)[AffineMapData<VolumeDim>::matrix_offset +
                           inertial_d * VolumeDim + grid_d] *
            representative_coords.get(grid_d)[p];
        expected_velocity +=
            (*affine_data)[AffineMapData<VolumeDim>::velocity_matrix_offset +
                           inertial_d * VolumeDim + grid_d] *
            representative_coords.get(grid_d)[p];
      }
      if (not within_tolerance(inertial_coords.get(inertial_d)[p],
                               expected_coord) or
          not within_tolerance(mesh_velocity.get(inertial_d)[p],
                               expected_velocity)) {
        return false;
      }
      for (size_t grid_d = 0; grid_d < VolumeDim; ++grid_d) {
        if (not within_tolerance(
                inverse_jacobian.get(grid_d, inertial_d)[p],
                (*affine_data)
                    [AffineMapData<VolumeDim>::inverse_jacobian_offset +
                     grid_d * VolumeDim + inertial_d])) {
          return false;
        }
      }
    }
  }
  return true;
}

template <size_t VolumeDim>
bool stage_device_map_data(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const evolution::Kokkos::PackedTopology<gh::System<VolumeDim>>&
        packed_topology,
    const Domain<VolumeDim>& domain,
    const domain::FunctionsOfTimeMap& functions_of_time, const double time) {
  const size_t num_elements = packed_topology.local_element_ids.size();
  constexpr size_t affine_data_size = AffineMapData<VolumeDim>::size;
  static_assert(affine_data_size >= RotScaleTransMapData<VolumeDim>::size);
  if (packed_geometry->affine_map_data_device.extent(0) != num_elements or
      packed_geometry->affine_map_data_device.extent(1) != affine_data_size) {
    packed_geometry->affine_map_data_device =
        typename evolution::Kokkos::PackedGeometry<
            gh::System<VolumeDim>>::device_affine_map_data_type{
            "GhBatchedAffineMapData", num_elements, affine_data_size};
  }
  auto host_affine_data =
      ::Kokkos::create_mirror_view(packed_geometry->affine_map_data_device);
  std::map<size_t, std::array<double, affine_data_size>> data_by_block{};
  for (size_t e = 0; e < num_elements; ++e) {
    const int map_kind = packed_geometry->element_map_kind_host[e];
    if (map_kind != static_cast<int>(BatchedMapKind::ComposedAffine) and
        map_kind != static_cast<int>(BatchedMapKind::RotScaleTransTransition)) {
      continue;
    }
    const auto& element_id = packed_topology.local_element_ids[e];
    const size_t block_id = element_id.block_id();
    auto data = data_by_block.find(block_id);
    if (data == data_by_block.end()) {
      const size_t point_begin = packed_topology.element_point_offsets_host[e];
      const size_t points_per_element = packed_topology.points_per_element;
      tnsr::I<DataVector, VolumeDim, Frame::Grid> representative_coords{
          points_per_element};
      for (size_t d = 0; d < VolumeDim; ++d) {
        for (size_t p = 0; p < points_per_element; ++p) {
          representative_coords.get(d)[p] =
              packed_geometry->grid_coordinates_host[d][point_begin + p];
        }
      }
      std::array<double, affine_data_size> inferred_data{};
      bool data_is_ready = false;
      const auto& descriptors =
          packed_geometry->block_map_descriptors_host.at(block_id);
      if (map_kind ==
          static_cast<int>(BatchedMapKind::RotScaleTransTransition)) {
        ASSERT(descriptors.has_value() and descriptors->size() == 1,
               "A RotScaleTrans transition map requires one exact map "
               "descriptor.");
        const auto& rot_scale_trans = std::get<
            domain::CoordinateMaps::TimeDependent::RotScaleTransMapDescriptor>(
            descriptors->front());
        using map_data = RotScaleTransMapData<VolumeDim>;
        inferred_data[map_data::inner_radius_offset] =
            rot_scale_trans.inner_radius;
        inferred_data[map_data::outer_radius_offset] =
            rot_scale_trans.outer_radius;
        for (size_t d = 0; d < VolumeDim; ++d) {
          inferred_data[map_data::rotation_offset + d * VolumeDim + d] = 1.0;
        }
        if (rot_scale_trans.rotation_name.has_value()) {
          const auto& rotation_f_of_t =
              *functions_of_time.at(*rot_scale_trans.rotation_name);
          const Matrix rotation =
              rotation_matrix<VolumeDim>(time, rotation_f_of_t);
          const Matrix rotation_deriv =
              rotation_matrix_deriv<VolumeDim>(time, rotation_f_of_t);
          for (size_t target_d = 0; target_d < VolumeDim; ++target_d) {
            for (size_t source_d = 0; source_d < VolumeDim; ++source_d) {
              inferred_data[map_data::rotation_offset + target_d * VolumeDim +
                            source_d] = rotation(target_d, source_d);
              inferred_data[map_data::rotation_deriv_offset +
                            target_d * VolumeDim + source_d] =
                  rotation_deriv(target_d, source_d);
            }
          }
        }
        inferred_data[map_data::scale_a_offset] = 1.0;
        inferred_data[map_data::scale_b_offset] = 1.0;
        if (rot_scale_trans.scale_names.has_value()) {
          const auto scale_a =
              functions_of_time.at(rot_scale_trans.scale_names->first)
                  ->func_and_deriv(time);
          const auto scale_b =
              functions_of_time.at(rot_scale_trans.scale_names->second)
                  ->func_and_deriv(time);
          inferred_data[map_data::scale_a_offset] = scale_a[0][0];
          inferred_data[map_data::scale_b_offset] = scale_b[0][0];
          inferred_data[map_data::scale_a_deriv_offset] = scale_a[1][0];
          inferred_data[map_data::scale_b_deriv_offset] = scale_b[1][0];
        }
        if (rot_scale_trans.translation_name.has_value()) {
          const auto translation =
              functions_of_time.at(*rot_scale_trans.translation_name)
                  ->func_and_deriv(time);
          ASSERT(translation[0].size() == VolumeDim and
                     translation[1].size() == VolumeDim,
                 "RotScaleTrans translation function must have "
                     << VolumeDim << " components.");
          for (size_t d = 0; d < VolumeDim; ++d) {
            inferred_data[map_data::translation_offset + d] = translation[0][d];
            inferred_data[map_data::translation_deriv_offset + d] =
                translation[1][d];
          }
        }
        data_is_ready = true;
      }
      if (not data_is_ready and descriptors.has_value() and
          descriptors->size() == 1) {
        if (const auto* const rot_scale_trans =
                std::get_if<domain::CoordinateMaps::TimeDependent::
                                RotScaleTransMapDescriptor>(
                    &descriptors->front());
            rot_scale_trans != nullptr and
            rot_scale_trans->region !=
                domain::CoordinateMaps::TimeDependent::
                    RotScaleTransMapDescriptor::Region::Transition) {
          using descriptor_type =
              domain::CoordinateMaps::TimeDependent::RotScaleTransMapDescriptor;
          const bool use_inner_values =
              rot_scale_trans->region == descriptor_type::Region::Inner;
          double scale = 1.0;
          double scale_deriv = 0.0;
          if (rot_scale_trans->scale_names.has_value()) {
            const std::string& scale_name =
                use_inner_values ? rot_scale_trans->scale_names->first
                                 : rot_scale_trans->scale_names->second;
            const auto scale_and_deriv =
                functions_of_time.at(scale_name)->func_and_deriv(time);
            ASSERT(scale_and_deriv[0].size() == 1 and
                       scale_and_deriv[1].size() == 1,
                   "RotScaleTrans expansion function '" << scale_name
                                                        << "' must be scalar.");
            scale = scale_and_deriv[0][0];
            scale_deriv = scale_and_deriv[1][0];
          }

          Matrix rotation{VolumeDim, VolumeDim, 0.0};
          Matrix rotation_deriv{VolumeDim, VolumeDim, 0.0};
          for (size_t d = 0; d < VolumeDim; ++d) {
            rotation(d, d) = 1.0;
          }
          if (rot_scale_trans->rotation_name.has_value()) {
            const auto& rotation_f_of_t =
                *functions_of_time.at(*rot_scale_trans->rotation_name);
            rotation = rotation_matrix<VolumeDim>(time, rotation_f_of_t);
            rotation_deriv =
                rotation_matrix_deriv<VolumeDim>(time, rotation_f_of_t);
          }

          std::array<double, VolumeDim> translation{};
          std::array<double, VolumeDim> translation_deriv{};
          if (use_inner_values and
              rot_scale_trans->translation_name.has_value()) {
            const auto value_and_deriv =
                functions_of_time.at(*rot_scale_trans->translation_name)
                    ->func_and_deriv(time);
            ASSERT(value_and_deriv[0].size() == VolumeDim and
                       value_and_deriv[1].size() == VolumeDim,
                   "RotScaleTrans translation function '"
                       << *rot_scale_trans->translation_name << "' must have "
                       << VolumeDim << " components.");
            for (size_t d = 0; d < VolumeDim; ++d) {
              translation[d] = value_and_deriv[0][d];
              translation_deriv[d] = value_and_deriv[1][d];
            }
          }

          for (size_t target_d = 0; target_d < VolumeDim; ++target_d) {
            inferred_data[AffineMapData<VolumeDim>::translation_offset +
                          target_d] = translation[target_d];
            inferred_data
                [AffineMapData<VolumeDim>::velocity_translation_offset +
                 target_d] = translation_deriv[target_d];
            for (size_t source_d = 0; source_d < VolumeDim; ++source_d) {
              inferred_data[AffineMapData<VolumeDim>::matrix_offset +
                            target_d * VolumeDim + source_d] =
                  scale * rotation(target_d, source_d);
              inferred_data[AffineMapData<VolumeDim>::velocity_matrix_offset +
                            target_d * VolumeDim + source_d] =
                  scale_deriv * rotation(target_d, source_d) +
                  scale * rotation_deriv(target_d, source_d);
              inferred_data[AffineMapData<VolumeDim>::inverse_jacobian_offset +
                            source_d * VolumeDim + target_d] =
                  rotation(target_d, source_d) / scale;
            }
          }
          data_is_ready = true;
        }
      }
      if (not data_is_ready and
          not infer_composed_affine_map<VolumeDim>(
              make_not_null(&inferred_data), domain.blocks()[block_id],
              representative_coords, time, functions_of_time)) {
        return false;
      }
      data = data_by_block.emplace(block_id, inferred_data).first;
    }
    for (size_t component = 0; component < affine_data_size; ++component) {
      host_affine_data(e, component) = data->second[component];
    }
  }
  ::Kokkos::deep_copy(packed_geometry->affine_map_data_device,
                      host_affine_data);
  return true;
}

// Build the per-element device map plan once: classify every element's block,
// verify the specialized and general-affine device representations against the
// host map, record any shared scalar function-of-time name, and upload the
// per-element kind array to the device.
template <size_t VolumeDim>
void build_map_plan(
    const gsl::not_null<
        evolution::Kokkos::PackedGeometry<gh::System<VolumeDim>>*>
        packed_geometry,
    const evolution::Kokkos::PackedTopology<gh::System<VolumeDim>>&
        packed_topology,
    const Domain<VolumeDim>& domain,
    const domain::FunctionsOfTimeMap& functions_of_time, const double time) {
  const size_t num_elements = packed_topology.local_element_ids.size();
  const size_t points_per_element = packed_topology.points_per_element;
  packed_geometry->element_map_kind_host.assign(
      num_elements, static_cast<int>(BatchedMapKind::HostEvaluated));
  packed_geometry->block_map_descriptors_host.resize(domain.blocks().size());

  std::map<size_t, BatchedMapKind> kind_by_block{};
  for (size_t e = 0; e < num_elements; ++e) {
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& block = domain.blocks()[element_id.block_id()];
    const size_t block_id = element_id.block_id();

    BatchedMapKind kind = BatchedMapKind::HostEvaluated;
    const auto cached = kind_by_block.find(block_id);
    if (cached != kind_by_block.end()) {
      kind = cached->second;
    } else {
      kind = classify_block_map(block);
      if (kind != BatchedMapKind::Static and not block.has_distorted_frame()) {
        auto& descriptors =
            packed_geometry->block_map_descriptors_host[block_id];
        descriptors = block.moving_mesh_grid_to_inertial_map()
                          .time_dependent_map_descriptors();
        if (descriptors.has_value() and descriptors->size() == 1) {
          if (const auto* const rot_scale_trans =
                  std::get_if<domain::CoordinateMaps::TimeDependent::
                                  RotScaleTransMapDescriptor>(
                      &descriptors->front());
              rot_scale_trans != nullptr) {
            kind = rot_scale_trans->region ==
                           domain::CoordinateMaps::TimeDependent::
                               RotScaleTransMapDescriptor::Region::Transition
                       ? BatchedMapKind::RotScaleTransTransition
                       : BatchedMapKind::ComposedAffine;
          }
        }
        // Attempt to recognize a device-supported time-dependent map. First
        // check the inexpensive specialized single-function maps, then infer
        // and verify a general affine composition. Distorted-frame maps remain
        // on the host path.
        const auto& names =
            block.moving_mesh_grid_to_inertial_map().function_of_time_names();
        const size_t point_begin =
            packed_topology.element_point_offsets_host[e];
        tnsr::I<DataVector, VolumeDim, Frame::Grid> grid_coords{
            points_per_element};
        for (size_t d = 0; d < VolumeDim; ++d) {
          for (size_t p = 0; p < points_per_element; ++p) {
            grid_coords.get(d)[p] =
                packed_geometry->grid_coordinates_host[d][point_begin + p];
          }
        }
        if (kind == BatchedMapKind::HostEvaluated and names.size() == 1) {
          const std::string& name = *names.begin();
          if (verify_rigid_rotation_z<VolumeDim>(block, grid_coords, name, time,
                                                 functions_of_time)) {
            kind = BatchedMapKind::RigidRotationZ;
          } else if (verify_uniform_radial_scaling<VolumeDim>(
                         block, grid_coords, name, time, functions_of_time)) {
            kind = BatchedMapKind::UniformRadialScaling;
          }
          if (kind != BatchedMapKind::HostEvaluated) {
            if (packed_geometry->device_map_function_of_time_name.empty()) {
              packed_geometry->device_map_function_of_time_name = name;
            } else {
              ASSERT(packed_geometry->device_map_function_of_time_name == name,
                     "Batched device maps currently require a single shared "
                     "scalar function of time, but found both '"
                         << packed_geometry->device_map_function_of_time_name
                         << "' and '" << name << "'.");
            }
          }
        }
        if (kind == BatchedMapKind::HostEvaluated) {
          std::array<double, AffineMapData<VolumeDim>::size> affine_data{};
          if (infer_composed_affine_map<VolumeDim>(make_not_null(&affine_data),
                                                   block, grid_coords, time,
                                                   functions_of_time)) {
            kind = BatchedMapKind::ComposedAffine;
          }
        }
      }
      kind_by_block.emplace(block_id, kind);
    }
    packed_geometry->element_map_kind_host[e] = static_cast<int>(kind);
  }

  packed_geometry->element_map_kind_device =
      ::Kokkos::View<int*>("GhBatchedElementMapKind", num_elements);
  auto host_kind =
      ::Kokkos::create_mirror_view(packed_geometry->element_map_kind_device);
  for (size_t e = 0; e < num_elements; ++e) {
    host_kind(e) = packed_geometry->element_map_kind_host[e];
  }
  ::Kokkos::deep_copy(packed_geometry->element_map_kind_device, host_kind);
  packed_geometry->map_plan_initialized = true;
}

}  // namespace

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache>
void UpdateFrameFieldsBatchedWithTag<FunctionsOfTimeTag,
                                     AddFunctionsOfTimeToConstGlobalCache>::
    apply(const gsl::not_null<typename packed_geometry_tag::type*>
              packed_geometry,
          const typename packed_topology_tag::type& packed_topology,
          const Domain<volume_dim>& domain,
          const domain::FunctionsOfTimeMap& functions_of_time,
          const double time) {
  const size_t total_points = packed_topology.total_points;
  if (total_points == 0) {
    return;
  }

  const size_t num_elements = packed_topology.local_element_ids.size();
  const size_t points_per_element = packed_topology.points_per_element;
  ASSERT(points_per_element > 0,
         "Packed points-per-element must be positive for GH batched updates.");
  ASSERT(num_elements * points_per_element == total_points,
         "Packed topology metadata is inconsistent with total packed points.");

  ASSERT(packed_geometry->logical_to_grid_inverse_jacobian_host.extent(0) ==
                 num_elements and
             packed_geometry->logical_to_grid_inverse_jacobian_host.extent(1) ==
                 points_per_element and
             packed_geometry->logical_to_grid_inverse_jacobian_host.extent(2) ==
                 volume_dim * volume_dim,
         "Packed logical-to-grid inverse Jacobian cache must be initialized by "
         "domain-static initialization before UpdateFrameFieldsBatched.");

  ensure_frame_field_storage<volume_dim>(packed_geometry, num_elements,
                                         points_per_element, total_points);

  for (size_t d = 0; d < volume_dim; ++d) {
    ASSERT(packed_geometry->grid_coordinates_host[d].size() == total_points,
           "Packed grid coordinates are not initialized. Expected "
               << total_points << " points in component " << d << ", got "
               << packed_geometry->grid_coordinates_host[d].size() << ".");
  }

  // Classify each block's map once. Verified static, specialized scalar, and
  // general-affine maps are evaluated on the device; unsupported maps fall
  // back to the host path below.
  if (not packed_geometry->map_plan_initialized) {
    ensure_static_device_geometry<volume_dim>(packed_geometry, num_elements,
                                              points_per_element, total_points);
    build_map_plan<volume_dim>(packed_geometry, packed_topology, domain,
                               functions_of_time, time);
  }

  bool has_staged_map_data = false;
  for (const int kind : packed_geometry->element_map_kind_host) {
    has_staged_map_data =
        has_staged_map_data or
        kind == static_cast<int>(BatchedMapKind::ComposedAffine) or
        kind == static_cast<int>(BatchedMapKind::RotScaleTransTransition);
  }
  const bool staged_maps_ready =
      not has_staged_map_data or
      stage_device_map_data<volume_dim>(packed_geometry, packed_topology,
                                        domain, functions_of_time, time);

  bool has_host_evaluated = not staged_maps_ready;
  bool has_device_evaluated = false;
  bool all_static = true;
  for (size_t e = 0; e < num_elements; ++e) {
    const int kind = packed_geometry->element_map_kind_host[e];
    if (kind == static_cast<int>(BatchedMapKind::HostEvaluated)) {
      has_host_evaluated = true;
    } else {
      has_device_evaluated = true;
    }
    if (kind != static_cast<int>(BatchedMapKind::Static)) {
      all_static = false;
    }
  }
  if (not staged_maps_ready) {
    // A map that ceases to verify as affine falls back to an all-host update
    // for this stage.
    has_device_evaluated = false;
  }

  double fot_value = 0.0;
  double fot_value_dot = 0.0;
  if (has_device_evaluated and
      not packed_geometry->device_map_function_of_time_name.empty()) {
    const auto value_and_deriv =
        functions_of_time.at(packed_geometry->device_map_function_of_time_name)
            ->func_and_deriv(time);
    fot_value = value_and_deriv[0][0];
    fot_value_dot = value_and_deriv[1][0];
  }

  if (not has_host_evaluated) {
    ensure_static_device_geometry<volume_dim>(packed_geometry, num_elements,
                                              points_per_element, total_points);
    update_frame_fields_device<volume_dim>(packed_geometry, points_per_element,
                                           total_points, fot_value,
                                           fot_value_dot);
    // For a static domain the host geometry arrays are constant and were set
    // during initialization, so leave them untouched. When a moving map is
    // evaluated on the device, refresh the host coordinate and mesh-velocity
    // arrays so the host-side observation and apparent-horizon consumers see
    // the current geometry. (Narrowing this to observation/AH times is a later
    // milestone.)
    if (not all_static) {
      copy_device_array_to_host<volume_dim>(
          packed_geometry->inertial_coordinates_device,
          make_not_null(&packed_geometry->inertial_coordinates_host));
      copy_device_array_to_host<volume_dim>(
          packed_geometry->mesh_velocity_device,
          make_not_null(&packed_geometry->mesh_velocity_host));
    }
    return;
  }

  auto host_inverse_jacobian = ::Kokkos::create_mirror_view(
      packed_geometry->element_inverse_jacobian_device);

  for (size_t e = 0; e < num_elements; ++e) {
    if (staged_maps_ready and
        packed_geometry->element_map_kind_host[e] !=
            static_cast<int>(BatchedMapKind::HostEvaluated)) {
      continue;
    }
    const auto& element_id = packed_topology.local_element_ids[e];
    const auto& block = domain.blocks()[element_id.block_id()];

    const size_t point_begin = packed_topology.element_point_offsets_host[e];
    const size_t point_end = packed_topology.element_point_offsets_host[e + 1];
    const size_t element_points = point_end - point_begin;
    ASSERT(element_points == points_per_element,
           "GH batched frame update assumes uniform points per element.");

    tnsr::I<DataVector, volume_dim, Frame::Grid> grid_coords{element_points};
    for (size_t d = 0; d < volume_dim; ++d) {
      for (size_t p = 0; p < element_points; ++p) {
        grid_coords.get(d)[p] =
            packed_geometry->grid_coordinates_host[d][point_begin + p];
      }
    }

    if (block.is_time_dependent()) {
      const auto map_data = block.moving_mesh_grid_to_inertial_map()
                                .coords_frame_velocity_jacobians(
                                    grid_coords, time, functions_of_time);
      const auto& inertial_coords = std::get<0>(map_data);
      const auto& grid_to_inertial_inverse_jacobian = std::get<1>(map_data);
      const auto& grid_to_inertial_jacobian = std::get<2>(map_data);
      const auto& mesh_velocity = std::get<3>(map_data);

      std::optional<std::tuple<
          tnsr::I<DataVector, volume_dim, Frame::Distorted>,
          InverseJacobian<DataVector, volume_dim, Frame::Grid,
                          Frame::Distorted>,
          Jacobian<DataVector, volume_dim, Frame::Grid, Frame::Distorted>,
          tnsr::I<DataVector, volume_dim, Frame::Distorted>>>
          grid_to_distorted_map_data{};
      std::optional<std::tuple<
          tnsr::I<DataVector, volume_dim, Frame::Inertial>,
          InverseJacobian<DataVector, volume_dim, Frame::Distorted,
                          Frame::Inertial>,
          Jacobian<DataVector, volume_dim, Frame::Distorted, Frame::Inertial>,
          tnsr::I<DataVector, volume_dim, Frame::Inertial>>>
          distorted_to_inertial_map_data{};
      if (block.has_distorted_frame()) {
        grid_to_distorted_map_data =
            block.moving_mesh_grid_to_distorted_map()
                .coords_frame_velocity_jacobians(grid_coords, time,
                                                 functions_of_time);
        distorted_to_inertial_map_data =
            block.moving_mesh_distorted_to_inertial_map()
                .coords_frame_velocity_jacobians(
                    std::get<0>(*grid_to_distorted_map_data), time,
                    functions_of_time);
      }

      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              inertial_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] =
              mesh_velocity.get(d)[p];
          packed_geometry->distorted_coordinates_host[d][point_begin + p] =
              block.has_distorted_frame()
                  ? std::get<0>(*grid_to_distorted_map_data).get(d)[p]
                  : grid_coords.get(d)[p];
          packed_geometry
              ->grid_to_distorted_frame_velocity_host[d][point_begin + p] =
              block.has_distorted_frame()
                  ? std::get<3>(*grid_to_distorted_map_data).get(d)[p]
                  : 0.0;
        }

        for (size_t logical_d = 0; logical_d < volume_dim; ++logical_d) {
          for (size_t inertial_d = 0; inertial_d < volume_dim; ++inertial_d) {
            double composed_inverse_jacobian = 0.0;
            for (size_t grid_d = 0; grid_d < volume_dim; ++grid_d) {
              composed_inverse_jacobian +=
                  packed_geometry->logical_to_grid_inverse_jacobian_host(
                      e, p, logical_d * volume_dim + grid_d) *
                  grid_to_inertial_inverse_jacobian.get(grid_d, inertial_d)[p];
            }
            host_inverse_jacobian(e, p, logical_d * volume_dim + inertial_d) =
                composed_inverse_jacobian;
          }
        }
        for (size_t source_d = 0; source_d < volume_dim; ++source_d) {
          for (size_t target_d = 0; target_d < volume_dim; ++target_d) {
            const double grid_to_distorted_inverse_jacobian =
                block.has_distorted_frame()
                    ? std::get<1>(*grid_to_distorted_map_data)
                          .get(source_d, target_d)[p]
                    : (source_d == target_d ? 1.0 : 0.0);
            packed_geometry->grid_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                grid_to_distorted_inverse_jacobian;
            packed_geometry->distorted_to_inertial_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                block.has_distorted_frame()
                    ? std::get<2>(*distorted_to_inertial_map_data)
                          .get(source_d, target_d)[p]
                    : grid_to_inertial_jacobian.get(source_d, target_d)[p];
            packed_geometry->distorted_to_inertial_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                block.has_distorted_frame()
                    ? std::get<1>(*distorted_to_inertial_map_data)
                          .get(source_d, target_d)[p]
                    : grid_to_inertial_inverse_jacobian.get(source_d,
                                                            target_d)[p];

            double composed_logical_to_distorted_inverse_jacobian = 0.0;
            for (size_t grid_d = 0; grid_d < volume_dim; ++grid_d) {
              composed_logical_to_distorted_inverse_jacobian +=
                  packed_geometry->logical_to_grid_inverse_jacobian_host(
                      e, p, source_d * volume_dim + grid_d) *
                  (block.has_distorted_frame()
                       ? std::get<1>(*grid_to_distorted_map_data)
                             .get(grid_d, target_d)[p]
                       : (grid_d == target_d ? 1.0 : 0.0));
            }
            packed_geometry->logical_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                composed_logical_to_distorted_inverse_jacobian;
          }
        }
      }
    } else {
      for (size_t p = 0; p < element_points; ++p) {
        for (size_t d = 0; d < volume_dim; ++d) {
          packed_geometry->distorted_coordinates_host[d][point_begin + p] =
              grid_coords.get(d)[p];
          packed_geometry->inertial_coordinates_host[d][point_begin + p] =
              grid_coords.get(d)[p];
          packed_geometry->mesh_velocity_host[d][point_begin + p] = 0.0;
          packed_geometry
              ->grid_to_distorted_frame_velocity_host[d][point_begin + p] = 0.0;
        }
        for (size_t component = 0; component < volume_dim * volume_dim;
             ++component) {
          host_inverse_jacobian(e, p, component) =
              packed_geometry->logical_to_grid_inverse_jacobian_host(e, p,
                                                                     component);
          packed_geometry->logical_to_distorted_inverse_jacobian_host(
              e, p, component) =
              packed_geometry->logical_to_grid_inverse_jacobian_host(e, p,
                                                                     component);
        }
        for (size_t source_d = 0; source_d < volume_dim; ++source_d) {
          for (size_t target_d = 0; target_d < volume_dim; ++target_d) {
            packed_geometry->grid_to_distorted_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
            packed_geometry->distorted_to_inertial_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
            packed_geometry->distorted_to_inertial_inverse_jacobian_host(
                e, p, source_d * volume_dim + target_d) =
                source_d == target_d ? 1.0 : 0.0;
          }
        }
      }
    }
  }

  copy_host_array_to_device<volume_dim>(
      packed_geometry->inertial_coordinates_host,
      make_not_null(&packed_geometry->inertial_coordinates_device));
  copy_host_array_to_device<volume_dim>(
      packed_geometry->mesh_velocity_host,
      make_not_null(&packed_geometry->mesh_velocity_device));
  ::Kokkos::deep_copy(packed_geometry->element_inverse_jacobian_device,
                      host_inverse_jacobian);

  // In a mixed domain the host-evaluated entries above seed the device
  // arrays, then this kernel overwrites only the supported entries. Refresh
  // the host coordinate and velocity views afterward so host-side consumers
  // see one coherent geometry snapshot.
  if (has_device_evaluated) {
    update_frame_fields_device<volume_dim>(packed_geometry, points_per_element,
                                           total_points, fot_value,
                                           fot_value_dot);
    copy_device_array_to_host<volume_dim>(
        packed_geometry->inertial_coordinates_device,
        make_not_null(&packed_geometry->inertial_coordinates_host));
    copy_device_array_to_host<volume_dim>(
        packed_geometry->mesh_velocity_device,
        make_not_null(&packed_geometry->mesh_velocity_host));
  }
}

template struct UpdateFrameFieldsBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize>;
template struct UpdateFrameFieldsBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize, false>;
template struct UpdateFrameFieldsBatchedWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize>;
template struct UpdateFrameFieldsBatchedWithTag<
    ::control_system::Tags::FunctionsOfTimeInitialize, false>;

}  // namespace gh::Actions
