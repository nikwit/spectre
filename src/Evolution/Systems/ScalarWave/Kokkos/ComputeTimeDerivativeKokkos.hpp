// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/InterfaceHelpers.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
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
  using outgoing_boundary_data_tag =
      ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<volume_dim>;
  using device_variables_type = device_variables_tag::type;
  using device_dt_type = device_dt_variables_tag::type;
  using device_inverse_jacobian_type = device_inverse_jacobian_tag::type;
  using device_constraint_gamma2_type = device_constraint_gamma2_tag::type;
  using device_face_to_volume_index_map_type =
      device_face_to_volume_index_map_tag::type;
  using device_face_unit_normal_covector_type =
      device_face_unit_normal_covector_tag::type;
  using outgoing_boundary_data_type = outgoing_boundary_data_tag::type;
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
      const device_variables_type& device_vars,
      const device_inverse_jacobian_type& device_inverse_jacobian,
      const device_constraint_gamma2_type& device_constraint_gamma2,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_unit_normal_covector_type&
          device_face_unit_normal_covector,
      const Mesh<volume_dim>& mesh, const Element<volume_dim>& element,
      const TimeStepId& time_step_id);

 private:
  template <typename ParallelComponent, typename Metavariables>
  static void send_packaged_boundary_data(
      const outgoing_boundary_data_type& outgoing_boundary_data,
      const Element<volume_dim>& element, const TimeStepId& time_step_id,
      Parallel::GlobalCache<Metavariables>& cache) {
    Kokkos::fence("ComputeTimeDerivativeKokkosBeforeBoundarySend");
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

        // TODO(#kokkos-dg): Reorient payload data for non-aligned neighbors.
        // For now only metadata orientation is adjusted.
        data_for_neighbor.volume_mesh = orientation(data_for_neighbor.volume_mesh);
        data_for_neighbor.boundary_correction_mesh =
            orient_mesh_on_slice(data_for_neighbor.boundary_correction_mesh,
                                 direction.dimension(), orientation);

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
      tmpl::list<device_dt_variables_tag, outgoing_boundary_data_tag>;
  using argument_tags =
      tmpl::list<device_variables_tag, device_inverse_jacobian_tag,
                 device_constraint_gamma2_tag,
                 device_face_to_volume_index_map_tag,
                 device_face_unit_normal_covector_tag,
                 domain::Tags::Mesh<volume_dim>,
                 domain::Tags::Element<volume_dim>,
                 ::Tags::TimeStepId>;
  using inbox_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;

  static void apply(
      gsl::not_null<device_dt_type*> device_dt,
      gsl::not_null<outgoing_boundary_data_type*> outgoing_boundary_data,
      const device_variables_type& device_vars,
      const device_inverse_jacobian_type& device_inverse_jacobian,
      const device_constraint_gamma2_type& device_constraint_gamma2,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_unit_normal_covector_type&
          device_face_unit_normal_covector,
      const Mesh<volume_dim>& mesh, const Element<volume_dim>& element,
      const TimeStepId& time_step_id) {
    compute_volume_terms_and_package_boundary_data(
        device_dt, outgoing_boundary_data, device_vars, device_inverse_jacobian,
        device_constraint_gamma2, device_face_to_volume_index_map,
        device_face_unit_normal_covector, mesh, element, time_step_id);
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
    db::mutate<device_dt_variables_tag, outgoing_boundary_data_tag>(
        [&device_vars = db::get<device_variables_tag>(box),
         &device_inverse_jacobian = db::get<device_inverse_jacobian_tag>(box),
         &device_constraint_gamma2 = db::get<device_constraint_gamma2_tag>(box),
         &device_face_to_volume_index_map =
             db::get<device_face_to_volume_index_map_tag>(box),
         &device_face_unit_normal_covector =
             db::get<device_face_unit_normal_covector_tag>(box),
         &mesh = db::get<domain::Tags::Mesh<volume_dim>>(box),
         &time_step_id = db::get<::Tags::TimeStepId>(box),
         &element](
            gsl::not_null<device_dt_type*> device_dt,
            gsl::not_null<outgoing_boundary_data_type*> outgoing_boundary_data) {
          compute_volume_terms_and_package_boundary_data(
              device_dt, outgoing_boundary_data, device_vars,
              device_inverse_jacobian, device_constraint_gamma2,
              device_face_to_volume_index_map, device_face_unit_normal_covector,
              mesh, element, time_step_id);
        },
        make_not_null(&box));

    send_packaged_boundary_data<ParallelComponent>(
        db::get<outgoing_boundary_data_tag>(box), element,
        db::get<::Tags::TimeStepId>(box), cache);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace ScalarWave::Actions
