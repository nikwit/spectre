// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <exception>
#include <optional>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/MortarDataPrimitives.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ScalarWave::Actions {

struct ComputeTimeDerivativeKokkos {
 private:
  static constexpr size_t volume_dim = 3;
  using System = ScalarWave::System<volume_dim>;
  static constexpr bool use_nodegroup_dg_elements = false;

  using device_variables_tag = ScalarWave::KokkosTags::DeviceVariables<System>;
  using device_dt_variables_tag =
      ScalarWave::KokkosTags::DeviceDtVariables<System>;
  using device_inverse_jacobian_tag =
      ScalarWave::KokkosTags::DeviceInverseJacobian<volume_dim>;
  using device_constraint_gamma2_tag =
      ScalarWave::KokkosTags::DeviceConstraintGamma2;
  using device_face_to_volume_index_map_tag =
      ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<volume_dim>;
  using device_face_unit_normal_covector_tag =
      ScalarWave::KokkosTags::DeviceFaceUnitNormalCovector<volume_dim>;
  using device_mortar_data_tag =
      ScalarWave::KokkosTags::DeviceMortarData<volume_dim>;
  using outgoing_boundary_data_tag =
      ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<volume_dim>;
  using external_boundary_data_tag =
      ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<volume_dim>;
  using mortar_mesh_tag = evolution::dg::Tags::MortarMesh<volume_dim>;
  using device_variables_type = device_variables_tag::type;
  using device_dt_type = device_dt_variables_tag::type;
  using device_inverse_jacobian_type = device_inverse_jacobian_tag::type;
  using device_constraint_gamma2_type = device_constraint_gamma2_tag::type;
  using device_face_to_volume_index_map_type =
      device_face_to_volume_index_map_tag::type;
  using device_face_unit_normal_covector_type =
      device_face_unit_normal_covector_tag::type;
  using device_mortar_data_type = device_mortar_data_tag::type;
  using outgoing_boundary_data_type = outgoing_boundary_data_tag::type;
  using external_boundary_data_type = external_boundary_data_tag::type;
  using volume_time_derivative_terms =
      typename System::compute_volume_time_derivative_terms;
  using package_field_tags =
      typename ScalarWave::BoundaryCorrections::UpwindPenalty<
          volume_dim>::dg_package_field_tags;
  template <size_t I>
  using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;
  using device_package_field_tags =
      db::wrap_tags_in<::Tags::MirrorView, package_field_tags>;

 public:
  static void compute_volume_terms_and_package_boundary_data(
      gsl::not_null<device_dt_type*> device_dt,
      gsl::not_null<outgoing_boundary_data_type*> outgoing_boundary_data,
      gsl::not_null<external_boundary_data_type*> external_boundary_data,
      const device_variables_type& device_vars,
      const device_inverse_jacobian_type& device_inverse_jacobian,
      const device_constraint_gamma2_type& device_constraint_gamma2,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_unit_normal_covector_type&
          device_face_unit_normal_covector,
      const device_mortar_data_type& device_mortar_data,
      const typename mortar_mesh_tag::type& mortar_meshes,
      const ScalarWave::Tags::ConstraintGamma2::type& host_constraint_gamma2,
      const typename domain::Tags::ExternalBoundaryConditions<volume_dim>::type&
          external_boundary_conditions_by_block,
      const tnsr::I<DataVector, volume_dim, Frame::Inertial>&
          inertial_coordinates,
      const double time, const Mesh<volume_dim>& mesh,
      const Element<volume_dim>& element, const TimeStepId& time_step_id);
  static void orient_boundary_data_for_send(
      gsl::not_null<Variables<device_package_field_tags>*>
          oriented_boundary_data,
      const Variables<device_package_field_tags>& boundary_data,
      const ::Kokkos::View<size_t*>& oriented_mortar_grid_point_source_index);

 private:
  template <typename ParallelComponent, typename Metavariables>
  static void send_packaged_boundary_data(
      const outgoing_boundary_data_type& outgoing_boundary_data,
      const device_mortar_data_type& device_mortar_data,
      const Element<volume_dim>& element, const TimeStepId& time_step_id,
      Parallel::GlobalCache<Metavariables>& cache) {
    ::Kokkos::fence("ComputeTimeDerivativeKokkosBeforeBoundarySend");
    using inbox_tag = ScalarWave::KokkosTags::BoundaryCorrectionInbox<
        volume_dim, use_nodegroup_dg_elements>;
    auto& receiver_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);

    for (const auto& [direction, neighbors] : element.neighbors()) {
      for (const auto& neighbor : neighbors) {
        const auto& orientation = neighbors.orientation(neighbor);
        const auto direction_from_neighbor = orientation(direction.opposite());
        const DirectionalId<volume_dim> mortar_id{direction, neighbor};

        ASSERT(outgoing_boundary_data.count(mortar_id) == 1,
               "Missing outgoing Kokkos boundary data for mortar " << mortar_id
                                                                   << ".");

        auto data_for_neighbor = outgoing_boundary_data.at(mortar_id);
        if (not orientation.is_aligned()) {
          Variables<device_package_field_tags> oriented_boundary_data{
              data_for_neighbor.boundary_correction_data
                  .number_of_grid_points()};
          const auto& mortar_data = device_mortar_data.at(mortar_id);
          orient_boundary_data_for_send(
              make_not_null(&oriented_boundary_data),
              data_for_neighbor.boundary_correction_data,
              mortar_data.oriented_mortar_grid_point_source_index);
          data_for_neighbor.boundary_correction_data =
              std::move(oriented_boundary_data);
        }

        auto* const local_receiver = Parallel::local(receiver_proxy[neighbor]);
        ASSERT(local_receiver != nullptr,
               "ComputeTimeDerivativeKokkos boundary data exchange currently "
               "requires sender/receiver elements to be local to the same "
               "Charm++ PE. Could not find local receiver for neighbor "
                   << neighbor << ".");

        local_receiver->template receive_data<inbox_tag>(
            time_step_id,
            std::make_pair(DirectionalId{direction_from_neighbor, element.id()},
                           std::move(data_for_neighbor)));
      }
    }
  }

 public:
  using return_tags =
      tmpl::list<device_dt_variables_tag, outgoing_boundary_data_tag,
                 external_boundary_data_tag>;
  using argument_tags =
      tmpl::list<device_variables_tag, device_inverse_jacobian_tag,
                 device_constraint_gamma2_tag,
                 device_face_to_volume_index_map_tag,
                 device_face_unit_normal_covector_tag, device_mortar_data_tag,
                 mortar_mesh_tag, ScalarWave::Tags::ConstraintGamma2,
                 domain::Tags::Coordinates<volume_dim, Frame::Inertial>,
                 ::Tags::Time, domain::Tags::Mesh<volume_dim>,
                 domain::Tags::Element<volume_dim>, ::Tags::TimeStepId>;
  using inbox_tags = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::ExternalBoundaryConditions<volume_dim>>;

  static void apply(
      gsl::not_null<device_dt_type*> device_dt,
      gsl::not_null<outgoing_boundary_data_type*> outgoing_boundary_data,
      gsl::not_null<external_boundary_data_type*> external_boundary_data,
      const device_variables_type& device_vars,
      const device_inverse_jacobian_type& device_inverse_jacobian,
      const device_constraint_gamma2_type& device_constraint_gamma2,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_unit_normal_covector_type&
          device_face_unit_normal_covector,
      const device_mortar_data_type& device_mortar_data,
      const typename mortar_mesh_tag::type& mortar_meshes,
      const ScalarWave::Tags::ConstraintGamma2::type& host_constraint_gamma2,
      const typename domain::Tags::ExternalBoundaryConditions<volume_dim>::type&
          external_boundary_conditions_by_block,
      const tnsr::I<DataVector, volume_dim, Frame::Inertial>&
          inertial_coordinates,
      const double time, const Mesh<volume_dim>& mesh,
      const Element<volume_dim>& element, const TimeStepId& time_step_id) {
    compute_volume_terms_and_package_boundary_data(
        device_dt, outgoing_boundary_data, external_boundary_data, device_vars,
        device_inverse_jacobian, device_constraint_gamma2,
        device_face_to_volume_index_map, device_face_unit_normal_covector,
        device_mortar_data, mortar_meshes, host_constraint_gamma2,
        external_boundary_conditions_by_block, inertial_coordinates, time, mesh,
        element, time_step_id);
  }

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent,
            typename Metavariables>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    ASSERT(Parallel::number_of_nodes<size_t>(cache) == 1,
           "ComputeTimeDerivativeKokkos currently supports only a single "
           "node (no inter-node communication).");

    const auto& element = db::get<domain::Tags::Element<volume_dim>>(box);
    const auto& external_boundary_conditions_by_block =
        Parallel::get<domain::Tags::ExternalBoundaryConditions<volume_dim>>(
            cache);
    const auto& time_step_id = db::get<::Tags::TimeStepId>(box);
    try {
      db::mutate<device_dt_variables_tag, outgoing_boundary_data_tag,
                 external_boundary_data_tag>(
          [&device_vars = db::get<device_variables_tag>(box),
           &device_inverse_jacobian = db::get<device_inverse_jacobian_tag>(box),
           &device_constraint_gamma2 =
               db::get<device_constraint_gamma2_tag>(box),
           &device_face_to_volume_index_map =
               db::get<device_face_to_volume_index_map_tag>(box),
           &device_face_unit_normal_covector =
               db::get<device_face_unit_normal_covector_tag>(box),
           &device_mortar_data = db::get<device_mortar_data_tag>(box),
           &mortar_meshes = db::get<mortar_mesh_tag>(box),
           &host_constraint_gamma2 =
               db::get<ScalarWave::Tags::ConstraintGamma2>(box),
           &external_boundary_conditions_by_block,
           &inertial_coordinates =
               db::get<domain::Tags::Coordinates<volume_dim, Frame::Inertial>>(
                   box),
           &time = db::get<::Tags::Time>(box),
           &mesh = db::get<domain::Tags::Mesh<volume_dim>>(box), &time_step_id,
           &element](gsl::not_null<device_dt_type*> device_dt,
                     gsl::not_null<outgoing_boundary_data_type*>
                         outgoing_boundary_data,
                     gsl::not_null<external_boundary_data_type*>
                         external_boundary_data) {
            compute_volume_terms_and_package_boundary_data(
                device_dt, outgoing_boundary_data, external_boundary_data,
                device_vars, device_inverse_jacobian, device_constraint_gamma2,
                device_face_to_volume_index_map,
                device_face_unit_normal_covector, device_mortar_data,
                mortar_meshes, host_constraint_gamma2,
                external_boundary_conditions_by_block, inertial_coordinates,
                time, mesh, element, time_step_id);
          },
          make_not_null(&box));

      send_packaged_boundary_data<ParallelComponent>(
          db::get<outgoing_boundary_data_tag>(box),
          db::get<device_mortar_data_tag>(box), element, time_step_id, cache);
    } catch (const std::exception& e) {
      ERROR_NO_TRACE("ComputeTimeDerivativeKokkos action failed on element "
                     << element.id() << " at " << time_step_id << ": "
                     << e.what());
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace ScalarWave::Actions
