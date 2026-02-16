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
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ScalarWave::Actions {

namespace detail {

template <typename DbTagsList, typename... InboxTags>
bool receive_boundary_data_global_time_stepping_kokkos(
    gsl::not_null<db::DataBox<DbTagsList>*> box,
    gsl::not_null<tuples::TaggedTuple<InboxTags...>*> inboxes) {
  const auto number_of_neighbors =
      db::get<domain::Tags::Element<3>>(*box).number_of_neighbors();

  if (number_of_neighbors == 0) {
    db::mutate<ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<3>>(
        [](gsl::not_null<typename ScalarWave::KokkosTags::
                             IncomingBoundaryCorrectionData<3>::type*>
               incoming_boundary_data) { incoming_boundary_data->clear(); },
        box);
    return true;
  }

  const TimeStepId& temporal_id = db::get<::Tags::TimeStepId>(*box);

  auto& inbox =
      tuples::get<ScalarWave::KokkosTags::BoundaryCorrectionInbox<3, false>>(
          *inboxes);
  const auto received_record = inbox.find(temporal_id);
  if (received_record == inbox.end()) {
    return false;
  }

  auto& received_neighbor_data = received_record->second;
  if (received_neighbor_data.size() != number_of_neighbors) {
    ASSERT(received_neighbor_data.size() < number_of_neighbors,
           "Received too many messages: " << received_neighbor_data);
    return false;
  }

  db::mutate<ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<3>>(
      [&received_neighbor_data](
          gsl::not_null<typename ScalarWave::KokkosTags::
                            IncomingBoundaryCorrectionData<3>::type*>
              incoming_boundary_data) {
        *incoming_boundary_data = std::move(received_neighbor_data);
      },
      box);
  inbox.erase(received_record);
  return true;
}

}  // namespace detail

struct ApplyBoundaryCorrectionsToTimeDerivativeKokkos {
 private:
  static constexpr size_t volume_dim = 3;
  using System = ScalarWave::System<3>;
  using device_dt_variables_tag =
      ScalarWave::KokkosTags::DeviceDtVariables<System>;
  using boundary_correction_data_type =
      ScalarWave::KokkosTags::BoundaryCorrectionData<3>;
  using outgoing_boundary_data_type =
      DirectionalIdMap<3, boundary_correction_data_type>;
  using incoming_boundary_data_type =
      DirectionalIdMap<3, boundary_correction_data_type>;
  using external_boundary_data_type =
      DirectionalIdMap<3, boundary_correction_data_type>;
  using device_face_to_volume_index_map_tag =
      ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<3>;
  using device_face_normal_magnitude_tag =
      ScalarWave::KokkosTags::DeviceFaceNormalMagnitude<3>;
  using device_mortar_data_tag =
      ScalarWave::KokkosTags::DeviceMortarData<volume_dim>;
  using mortar_mesh_tag = evolution::dg::Tags::MortarMesh<volume_dim>;
  using device_dt_type = device_dt_variables_tag::type;
  using device_face_to_volume_index_map_type =
      device_face_to_volume_index_map_tag::type;
  using device_face_normal_magnitude_type =
      device_face_normal_magnitude_tag::type;
  using device_mortar_data_type = device_mortar_data_tag::type;
  using package_field_tags =
      typename ScalarWave::BoundaryCorrections::UpwindPenalty<
          3>::dg_package_field_tags;
  template <size_t I>
  using package_field_tag = tmpl::at<package_field_tags, tmpl::size_t<I>>;

 public:
  static void apply_boundary_corrections_on_device(
      gsl::not_null<device_dt_type*> device_dt,
      const outgoing_boundary_data_type& outgoing_boundary_data,
      const incoming_boundary_data_type& incoming_boundary_data,
      const external_boundary_data_type& external_boundary_data,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_normal_magnitude_type& device_face_normal_magnitude,
      const device_mortar_data_type& device_mortar_data,
      const typename mortar_mesh_tag::type& mortar_meshes, const Mesh<3>& mesh,
      const Element<3>& element);

  using inbox_tags =
      tmpl::list<ScalarWave::KokkosTags::BoundaryCorrectionInbox<3, false>>;
  using const_global_cache_tags = tmpl::list<>;

  static void apply(
      gsl::not_null<device_dt_type*> device_dt,
      const outgoing_boundary_data_type& outgoing_boundary_data,
      const incoming_boundary_data_type& incoming_boundary_data,
      const external_boundary_data_type& external_boundary_data,
      const device_face_to_volume_index_map_type&
          device_face_to_volume_index_map,
      const device_face_normal_magnitude_type& device_face_normal_magnitude,
      const device_mortar_data_type& device_mortar_data,
      const typename mortar_mesh_tag::type& mortar_meshes, const Mesh<3>& mesh,
      const Element<3>& element) {
    apply_boundary_corrections_on_device(
        device_dt, outgoing_boundary_data, incoming_boundary_data,
        external_boundary_data, device_face_to_volume_index_map,
        device_face_normal_magnitude, device_mortar_data, mortar_meshes, mesh,
        element);
  }

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    ASSERT(Parallel::number_of_nodes<size_t>(cache) == 1,
           "ApplyBoundaryCorrectionsToTimeDerivativeKokkos currently supports "
           "only a single node (no inter-node communication).");
    if (not detail::receive_boundary_data_global_time_stepping_kokkos(
            make_not_null(&box), make_not_null(&inboxes))) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    const auto& element = db::get<domain::Tags::Element<3>>(box);
    const auto& time_step_id = db::get<::Tags::TimeStepId>(box);
    try {
      db::mutate<device_dt_variables_tag>(
          [&outgoing_boundary_data = db::get<
               ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<3>>(box),
           &incoming_boundary_data = db::get<
               ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<3>>(box),
           &external_boundary_data = db::get<
               ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<3>>(box),
           &device_face_to_volume_index_map =
               db::get<device_face_to_volume_index_map_tag>(box),
           &device_face_normal_magnitude =
               db::get<device_face_normal_magnitude_tag>(box),
           &device_mortar_data = db::get<device_mortar_data_tag>(box),
           &mortar_meshes = db::get<mortar_mesh_tag>(box),
           &mesh = db::get<domain::Tags::Mesh<3>>(box),
           &element](gsl::not_null<device_dt_type*> device_dt) {
            apply_boundary_corrections_on_device(
                device_dt, outgoing_boundary_data, incoming_boundary_data,
                external_boundary_data, device_face_to_volume_index_map,
                device_face_normal_magnitude, device_mortar_data, mortar_meshes,
                mesh, element);
          },
          make_not_null(&box));
    } catch (const std::exception& e) {
      ERROR_NO_TRACE(
          "ApplyBoundaryCorrectionsToTimeDerivativeKokkos action failed on "
          "element "
          << element.id() << " at " << time_step_id << ": " << e.what());
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace ScalarWave::Actions
