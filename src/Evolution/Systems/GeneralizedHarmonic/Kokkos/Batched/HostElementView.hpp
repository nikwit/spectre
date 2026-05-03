// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/Kokkos/PackedDataBundles.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Batched {

struct HostElementView {
  using system = gh::System<3>;

  ElementId<3> element_id{};
  Element<3> element{};
  Mesh<3> mesh{};
  typename system::variables_tag::type evolved_vars{};
  tnsr::I<DataVector, 3, Frame::Grid> grid_coordinates{};
  tnsr::I<DataVector, 3, Frame::Distorted> distorted_coordinates{};
  tnsr::I<DataVector, 3, Frame::Inertial> inertial_coordinates{};
  tnsr::I<DataVector, 3, Frame::Inertial> mesh_velocity{};
  tnsr::I<DataVector, 3, Frame::Distorted> grid_to_distorted_frame_velocity{};
  InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Distorted>
      logical_to_distorted_inverse_jacobian{};
  InverseJacobian<DataVector, 3, Frame::Grid, Frame::Distorted>
      grid_to_distorted_inverse_jacobian{};
  Jacobian<DataVector, 3, Frame::Distorted, Frame::Inertial>
      distorted_to_inertial_jacobian{};
  InverseJacobian<DataVector, 3, Frame::Distorted, Frame::Inertial>
      distorted_to_inertial_inverse_jacobian{};
  InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>
      logical_to_inertial_inverse_jacobian{};
};

template <typename Callback>
void for_each_host_element_view(
    const evolution::Kokkos::PackedTopology<gh::System<3>>& topology,
    const evolution::Kokkos::PackedGeometry<gh::System<3>>& geometry,
    const evolution::Kokkos::PackedEvolutionState<gh::System<3>>& state,
    Callback&& callback) {
  using system = gh::System<3>;
  using host_variables_type = typename system::variables_tag::type;

  ASSERT(topology.local_element_ids.size() == topology.local_elements.size(),
         "Packed topology has mismatched element id and element storage.");
  ASSERT(
      topology.local_element_ids.size() == topology.element_extents_host.size(),
      "Packed topology has mismatched element id and extent storage.");
  ASSERT(topology.element_point_offsets_host.size() ==
             topology.local_element_ids.size() + 1,
         "Packed topology must have one more point offset than elements.");
  ASSERT(geometry.element_inverse_jacobian_device.extent(0) ==
             topology.local_element_ids.size(),
         "Packed inverse Jacobian element count does not match topology.");
  ASSERT(geometry.element_inverse_jacobian_device.extent(2) == 9,
         "Packed inverse Jacobian must store 3x3 components.");
  for (size_t d = 0; d < 3; ++d) {
    ASSERT(gsl::at(geometry.grid_coordinates_host, d).size() ==
               topology.total_points,
           "Packed grid coordinates do not match total points.");
    ASSERT(gsl::at(geometry.distorted_coordinates_host, d).size() ==
               topology.total_points,
           "Packed distorted coordinates do not match total points.");
    ASSERT(gsl::at(geometry.inertial_coordinates_host, d).size() ==
               topology.total_points,
           "Packed inertial coordinates do not match total points.");
    ASSERT(
        gsl::at(geometry.mesh_velocity_host, d).size() == topology.total_points,
        "Packed mesh velocity does not match total points.");
    ASSERT(gsl::at(geometry.grid_to_distorted_frame_velocity_host, d).size() ==
               topology.total_points,
           "Packed grid-to-distorted frame velocity does not match total "
           "points.");
  }

  host_variables_type host_vars{topology.total_points};
  copy_to_host(make_not_null(&host_vars), state.device_variables);
  const auto host_inverse_jacobian = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace{}, geometry.element_inverse_jacobian_device);

  for (size_t e = 0; e < topology.local_element_ids.size(); ++e) {
    const size_t point_begin = topology.element_point_offsets_host[e];
    const size_t point_end = topology.element_point_offsets_host[e + 1];
    ASSERT(point_end >= point_begin,
           "Packed element point offsets must be nondecreasing.");
    const size_t element_points = point_end - point_begin;
    ASSERT(element_points == geometry.element_inverse_jacobian_device.extent(1),
           "Element point count does not match inverse-Jacobian storage.");

    HostElementView view{};
    view.element_id = topology.local_element_ids[e];
    view.element = topology.local_elements[e];
    view.mesh =
        Mesh<3>{topology.element_extents_host[e], topology.uniform_basis_host,
                topology.uniform_quadrature_host};
    ASSERT(element_points == view.mesh.number_of_grid_points(),
           "Packed element point count does not match reconstructed mesh.");
    view.evolved_vars = host_variables_type{element_points};
    view.grid_coordinates = tnsr::I<DataVector, 3, Frame::Grid>{element_points};
    view.distorted_coordinates =
        tnsr::I<DataVector, 3, Frame::Distorted>{element_points};
    view.inertial_coordinates =
        tnsr::I<DataVector, 3, Frame::Inertial>{element_points};
    view.mesh_velocity =
        tnsr::I<DataVector, 3, Frame::Inertial>{element_points};
    view.grid_to_distorted_frame_velocity =
        tnsr::I<DataVector, 3, Frame::Distorted>{element_points};
    view.logical_to_distorted_inverse_jacobian =
        InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Distorted>{
            element_points};
    view.grid_to_distorted_inverse_jacobian =
        InverseJacobian<DataVector, 3, Frame::Grid, Frame::Distorted>{
            element_points};
    view.distorted_to_inertial_jacobian =
        Jacobian<DataVector, 3, Frame::Distorted, Frame::Inertial>{
            element_points};
    view.distorted_to_inertial_inverse_jacobian =
        InverseJacobian<DataVector, 3, Frame::Distorted, Frame::Inertial>{
            element_points};
    view.logical_to_inertial_inverse_jacobian =
        InverseJacobian<DataVector, 3, Frame::ElementLogical, Frame::Inertial>{
            element_points};

    for (size_t component = 0;
         component < host_variables_type::number_of_independent_components;
         ++component) {
      for (size_t p = 0; p < element_points; ++p) {
        view.evolved_vars.data()[component * element_points + p] =
            host_vars
                .data()[component * topology.total_points + point_begin + p];
      }
    }
    for (size_t d = 0; d < 3; ++d) {
      for (size_t p = 0; p < element_points; ++p) {
        view.grid_coordinates.get(d)[p] =
            gsl::at(geometry.grid_coordinates_host, d)[point_begin + p];
        view.distorted_coordinates.get(d)[p] =
            gsl::at(geometry.distorted_coordinates_host, d)[point_begin + p];
        view.inertial_coordinates.get(d)[p] =
            gsl::at(geometry.inertial_coordinates_host, d)[point_begin + p];
        view.mesh_velocity.get(d)[p] =
            gsl::at(geometry.mesh_velocity_host, d)[point_begin + p];
        view.grid_to_distorted_frame_velocity.get(d)[p] = gsl::at(
            geometry.grid_to_distorted_frame_velocity_host, d)[point_begin + p];
      }
    }
    for (size_t logical_d = 0; logical_d < 3; ++logical_d) {
      for (size_t inertial_d = 0; inertial_d < 3; ++inertial_d) {
        for (size_t p = 0; p < element_points; ++p) {
          view.logical_to_inertial_inverse_jacobian.get(logical_d,
                                                        inertial_d)[p] =
              host_inverse_jacobian(e, p, logical_d * 3 + inertial_d);
        }
      }
    }
    for (size_t source_d = 0; source_d < 3; ++source_d) {
      for (size_t target_d = 0; target_d < 3; ++target_d) {
        for (size_t p = 0; p < element_points; ++p) {
          view.logical_to_distorted_inverse_jacobian.get(source_d,
                                                         target_d)[p] =
              geometry.logical_to_distorted_inverse_jacobian_host(
                  e, p, source_d * 3 + target_d);
          view.grid_to_distorted_inverse_jacobian.get(source_d, target_d)[p] =
              geometry.grid_to_distorted_inverse_jacobian_host(
                  e, p, source_d * 3 + target_d);
          view.distorted_to_inertial_jacobian.get(source_d, target_d)[p] =
              geometry.distorted_to_inertial_jacobian_host(
                  e, p, source_d * 3 + target_d);
          view.distorted_to_inertial_inverse_jacobian.get(source_d,
                                                          target_d)[p] =
              geometry.distorted_to_inertial_inverse_jacobian_host(
                  e, p, source_d * 3 + target_d);
        }
      }
    }

    callback(view);
  }
}

Variables<ah::source_vars<3>> make_ah_source_vars(const HostElementView& view);

}  // namespace gh::Batched
