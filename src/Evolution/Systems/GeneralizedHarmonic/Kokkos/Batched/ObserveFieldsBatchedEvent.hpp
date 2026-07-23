// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateFrameFieldsBatched.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Events {

/*!
 * \brief Write the packed batched GH evolved variables and inertial
 * coordinates as volume data, one grid per packed element.
 *
 * The output layout matches the production `ObserveFields` event: grids are
 * named by their `ElementId`, tensor components use the standard component
 * suffixes, and the observation id is keyed on `<subfile>.vol`. This makes
 * batched volume dumps directly comparable, point by point, with volume
 * dumps from the host DG executables on the same domain.
 */
class ObserveFieldsBatched : public Event {
 public:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using host_variables_type = typename system::variables_tag::type;
  using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;

  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "Name of volume data subfile (without leading slash)."};
  };

  explicit ObserveFieldsBatched(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ObserveFieldsBatched);  // NOLINT

  using options = tmpl::list<SubfileName>;
  static constexpr Options::String help =
      "Write batched GH evolved variables and inertial coordinates as "
      "volume data, one grid per packed element.";

  ObserveFieldsBatched() = default;
  explicit ObserveFieldsBatched(const std::string& subfile_name)
      : subfile_path_("/" + subfile_name) {}

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<packed_geometry_tag>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<system>,
                 evolution::Kokkos::Tags::PackedEvolutionState<system>,
                 ::Tags::Time>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const gsl::not_null<typename packed_geometry_tag::type*> packed_geometry,
      const evolution::Kokkos::PackedTopology<system>& packed_topology,
      const evolution::Kokkos::PackedEvolutionState<system>&
          packed_evolution_state,
      const double time, Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& /*observation_value*/) const {
    const size_t num_elements = packed_topology.local_element_ids.size();
    if (num_elements == 0) {
      return;
    }
    gh::Actions::synchronize_frame_fields_to_host<volume_dim>(packed_geometry);
    const size_t points_per_element = packed_topology.points_per_element;
    ASSERT(packed_topology.total_points == num_elements * points_per_element,
           "ObserveFieldsBatched requires a uniform mesh: total points "
               << packed_topology.total_points << ", elements " << num_elements
               << ", points per element " << points_per_element);
    for (size_t d = 0; d < volume_dim; ++d) {
      ASSERT(packed_geometry->inertial_coordinates_host[d].size() ==
                 packed_topology.total_points,
             "Packed inertial-coordinate sizes do not match total packed "
             "points.");
    }

    host_variables_type host_vars{packed_topology.total_points};
    copy_to_host(make_not_null(&host_vars),
                 packed_evolution_state.device_variables);
    const auto& spacetime_metric =
        get<gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>>(
            host_vars);
    const auto& pi =
        get<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>(host_vars);
    const auto& phi =
        get<gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>(host_vars);

    const Mesh<volume_dim> mesh{packed_topology.uniform_extents_host,
                                packed_topology.uniform_basis_host,
                                packed_topology.uniform_quadrature_host};

    std::vector<ElementVolumeData> volume_data{};
    volume_data.reserve(num_elements);
    const auto slice_component = [&points_per_element](
                                     const DataVector& component,
                                     const size_t point_offset) {
      DataVector result{points_per_element};
      for (size_t p = 0; p < points_per_element; ++p) {
        result[p] = component[point_offset + p];
      }
      return result;
    };
    for (size_t e = 0; e < num_elements; ++e) {
      const size_t point_offset = packed_topology.element_point_offsets_host[e];
      std::vector<TensorComponent> components{};
      components.reserve(spacetime_metric.size() + pi.size() + phi.size() +
                         volume_dim);
      for (size_t d = 0; d < volume_dim; ++d) {
        components.emplace_back(
            "InertialCoordinates" +
                tnsr::I<DataVector, volume_dim,
                        Frame::Inertial>::component_suffix(d),
            slice_component(packed_geometry->inertial_coordinates_host[d],
                            point_offset));
      }
      for (size_t i = 0; i < spacetime_metric.size(); ++i) {
        components.emplace_back(
            "SpacetimeMetric" + spacetime_metric.component_suffix(i),
            slice_component(spacetime_metric[i], point_offset));
      }
      for (size_t i = 0; i < pi.size(); ++i) {
        components.emplace_back("Pi" + pi.component_suffix(i),
                                slice_component(pi[i], point_offset));
      }
      for (size_t i = 0; i < phi.size(); ++i) {
        components.emplace_back("Phi" + phi.component_suffix(i),
                                slice_component(phi[i], point_offset));
      }
      volume_data.emplace_back(packed_topology.local_element_ids[e],
                               std::move(components), mesh);
    }

    auto& volume_writer = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    Parallel::threaded_action<observers::ThreadedActions::WriteVolumeData>(
        volume_writer[0], Parallel::get<observers::Tags::VolumeFileName>(cache),
        subfile_path_, observers::ObservationId{time, subfile_path_ + ".vol"},
        std::move(volume_data));
  }

  using is_ready_argument_tags = tmpl::list<>;
  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return false; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Event::pup(p);
    p | subfile_path_;
  }

 private:
  std::string subfile_path_ = "/BatchedVolumeData";
};

}  // namespace gh::Events
