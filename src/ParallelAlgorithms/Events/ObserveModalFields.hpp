// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <functional>
#include <initializer_list>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataBox/ValidateSelection.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/ModalVector.hpp"
#include "Domain/Structure/BlockGroups.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/GetSectionObservationKey.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Auto.hpp"
#include "Options/String.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/OptionalHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <size_t Dim>
class Mesh;
namespace Frame {
struct Inertial;
}  // namespace Frame
/// \endcond

namespace dg {
namespace Events {
/// \cond
template <size_t VolumeDim, typename Tensors,
          typename NonTensorComputeTagsList = tmpl::list<>,
          typename ArraySectionIdTag = void>
class ObserveModalFields;
/// \endcond

/*!
 * \ingroup DiscontinuousGalerkinGroup
 * \brief Observe volume tensor fields in a modal basis.
 *
 * Writes modal expansion coefficients for the selected tensor components. If
 * an `interpolation_mesh` is provided the coefficients are truncated to the
 * specified extents, requiring that the mesh share the element's basis and
 * quadrature and that no extent increases.
 */
template <size_t VolumeDim, typename... Tensors,
          typename... NonTensorComputeTags, typename ArraySectionIdTag>
class ObserveModalFields<VolumeDim, tmpl::list<Tensors...>,
                         tmpl::list<NonTensorComputeTags...>, ArraySectionIdTag>
    : public Event {
 public:
  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "The name of the subfile inside the HDF5 file without an extension and "
        "without a preceding '/'."};
  };

  /// \cond
  explicit ObserveModalFields(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ObserveModalFields);  // NOLINT
  /// \endcond

  struct VariablesToObserve {
    static constexpr Options::String help = "Subset of variables to observe";
    using type = std::vector<std::string>;
    static size_t lower_bound_on_size() { return 1; }
  };

  struct InterpolateToMesh {
    using type = Options::Auto<Mesh<VolumeDim>, Options::AutoLabel::None>;
    static constexpr Options::String help =
        "An optional mesh defining the modal truncation extents. The mesh must "
        "share the element's basis and quadrature and cannot increase the "
        "number of points in any dimension.";
  };

  struct BlocksToObserve {
    using type =
        Options::Auto<std::vector<std::string>, Options::AutoLabel::All>;
    static constexpr Options::String help = {
        "A list of block and group names on which to observe."};
  };

  using options = tmpl::list<SubfileName, VariablesToObserve, BlocksToObserve,
                             InterpolateToMesh>;

  static constexpr Options::String help =
      "Observe volume tensor fields in a modal basis.\n"
      "\n"
      "Writes modal coefficients for:\n"
      " * InertialCoordinates\n"
      " * Tensors listed in the 'VariablesToObserve' option\n"
      "\n"
      "The coefficients are written as double-precision real values.";

  ObserveModalFields() = default;

  ObserveModalFields(
      const std::string& subfile_name,
      const std::vector<std::string>& variables_to_observe,
      std::optional<std::vector<std::string>> active_block_or_block_groups = {},
      std::optional<Mesh<VolumeDim>> truncation_mesh = {},
      std::optional<std::string> dependency = {},
      const Options::Context& context = {});

  ObserveModalFields(
      const std::string& subfile_name,
      const std::vector<std::string>& variables_to_observe,
      std::optional<std::vector<std::string>> active_block_or_block_groups = {},
      std::optional<Mesh<VolumeDim>> truncation_mesh = {},
      const Options::Context& context = {})
      : ObserveModalFields(subfile_name, variables_to_observe,
                           std::move(active_block_or_block_groups),
                           std::move(truncation_mesh), std::nullopt, context) {}

  using compute_tags_for_observation_box =
      tmpl::list<Tensors..., NonTensorComputeTags...>;

  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<::Tags::ObservationBox,
                                   ::Events::Tags::ObserverMesh<VolumeDim>>;

  template <typename DataBoxType, typename ComputeTagsList,
            typename Metavariables, typename ParallelComponent>
  void operator()(const ObservationBox<DataBoxType, ComputeTagsList>& box,
                  const Mesh<VolumeDim>& mesh,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ElementId<VolumeDim>& array_index,
                  const ParallelComponent* const component,
                  const ObservationValue& observation_value) const {
    if (not active_block(get<domain::Tags::Domain<VolumeDim>>(box),
                         array_index)) {
      return;
    }
    const std::optional<std::string> section_observation_key =
        observers::get_section_observation_key<ArraySectionIdTag>(box);
    if (not section_observation_key.has_value()) {
      return;
    }
    call_operator_impl(subfile_path_ + *section_observation_key,
                       variables_to_observe_, truncation_mesh_, mesh, box,
                       cache, array_index, component, observation_value,
                       dependency_);
  }

  template <typename DataBoxType, typename ComputeTagsList,
            typename Metavariables, typename ParallelComponent>
  static void call_operator_impl(
      const std::string& subfile_path,
      const std::unordered_set<std::string>& variables_to_observe,
      const std::optional<Mesh<VolumeDim>>& truncation_mesh,
      const Mesh<VolumeDim>& mesh,
      const ObservationBox<DataBoxType, ComputeTagsList>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<VolumeDim>& element_id,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value,
      const std::optional<std::string>& dependency) {
    const Mesh<VolumeDim>& mesh_for_output =
        truncation_mesh.has_value() ? truncation_mesh.value() : mesh;
    if (truncation_mesh.has_value()) {
      for (size_t d = 0; d < VolumeDim; ++d) {
        ASSERT(mesh_for_output.basis(d) == mesh.basis(d),
               "Modal truncation requires identical bases. Element basis "
                   << mesh.basis(d) << ", requested basis "
                   << mesh_for_output.basis(d));
        ASSERT(mesh_for_output.quadrature(d) == mesh.quadrature(d),
               "Modal truncation requires identical quadratures. Element "
               "quadrature "
                   << mesh.quadrature(d) << ", requested quadrature "
                   << mesh_for_output.quadrature(d));
        ASSERT(mesh_for_output.extents(d) <= mesh.extents(d),
               "Cannot increase resolution when truncating modal data. "
               "Requested extent "
                   << mesh_for_output.extents(d) << " exceeds element extent "
                   << mesh.extents(d) << " in dimension " << d);
      }
    }

    const auto modal_to_datavector = [](const ModalVector& modal) {
      DataVector result(modal.size());
      for (size_t i = 0; i < modal.size(); ++i) {
        result[i] = modal[i];
      }
      return result;
    };

    const auto truncate_modal = [&mesh,
                                 &mesh_for_output](const auto& modal_in) {
      using ModalType = std::decay_t<decltype(modal_in)>;
      if (mesh_for_output.number_of_grid_points() ==
          mesh.number_of_grid_points()) {
        return ModalType{modal_in};
      }
      ModalType modal_out(mesh_for_output.number_of_grid_points());
      const Index<VolumeDim> source_extents(mesh.extents());
      const Index<VolumeDim> target_extents(mesh_for_output.extents());
      for (size_t target_linear = 0;
           target_linear < mesh_for_output.number_of_grid_points();
           ++target_linear) {
        const Index<VolumeDim> target_multi =
            expanded_index(target_linear, target_extents);
        Index<VolumeDim> source_multi{0};
        for (size_t d = 0; d < VolumeDim; ++d) {
          source_multi[d] = target_multi[d];
        }
        const size_t source_linear =
            collapsed_index(source_multi, source_extents);
        modal_out[target_linear] = modal_in[source_linear];
      }
      return modal_out;
    };

    std::vector<TensorComponent> components;
    components.reserve(alg::accumulate(
        std::initializer_list<size_t>{
            std::decay_t<decltype(value(typename Tensors::type{}))>::size()...},
        0_st));

    const auto record_tensor_component_impl =
        [&components](DataVector&& tensor_component,
                      const std::string& component_name) {
          components.emplace_back(component_name, std::move(tensor_component));
        };

    const auto record_tensor_components_impl =
        [&record_tensor_component_impl, &mesh, &truncate_modal,
         &modal_to_datavector](const auto& tensor,
                               const std::string& tag_name) {
          using TensorType = std::decay_t<decltype(tensor)>;
          using VectorType = typename TensorType::type;
          for (size_t i = 0; i < tensor.size(); ++i) {
            const std::string component_name =
                tag_name + tensor.component_suffix(i);
            static_assert(std::is_same_v<VectorType, DataVector>,
                          "ObserveModalFields assumes real DataVector data.");
            const ModalVector element_modal =
                to_modal_coefficients(tensor[i], mesh);
            const ModalVector truncated_modal = truncate_modal(element_modal);
            record_tensor_component_impl(modal_to_datavector(truncated_modal),
                                         component_name);
          }
        };
    const auto record_tensor_components = [&box, &record_tensor_components_impl,
                                           &variables_to_observe](
                                              const auto tensor_tag_v) {
      using tensor_tag = tmpl::type_from<decltype(tensor_tag_v)>;
      const std::string tag_name = db::tag_name<tensor_tag>();
      if (const auto var_to_observe = variables_to_observe.find(tag_name);
          var_to_observe != variables_to_observe.end()) {
        const auto& tensor = get<tensor_tag>(box);
        if (not has_value(tensor)) {
          [[maybe_unused]] static bool t =
              ObserveModalFields::print_warning_about_optional<tensor_tag>();
          return;
        }
        record_tensor_components_impl(value(tensor), tag_name);
      }
    };
    EXPAND_PACK_LEFT_TO_RIGHT(record_tensor_components(tmpl::type_<Tensors>{}));

    const Parallel::ArrayComponentId array_component_id{
        std::add_pointer_t<ParallelComponent>{nullptr},
        Parallel::ArrayIndex<ElementId<VolumeDim>>{element_id}};
    ElementVolumeData element_volume_data{element_id, std::move(components),
                                          mesh_for_output};
    observers::ObservationId observation_id{observation_value.value,
                                            subfile_path + ".vol"};

    auto& local_observer = *Parallel::local_branch(
        Parallel::get_parallel_component<
            tmpl::conditional_t<Parallel::is_nodegroup_v<ParallelComponent>,
                                observers::ObserverWriter<Metavariables>,
                                observers::Observer<Metavariables>>>(cache));

    if constexpr (Parallel::is_nodegroup_v<ParallelComponent>) {
      std::unordered_map<Parallel::ArrayComponentId,
                         std::vector<ElementVolumeData>>
          data_to_send{};
      data_to_send[array_component_id] =
          std::vector{std::move(element_volume_data)};
      Parallel::threaded_action<
          observers::ThreadedActions::ContributeVolumeDataToWriter>(
          local_observer, std::move(observation_id), array_component_id,
          subfile_path, std::move(data_to_send), dependency);
    } else {
      Parallel::simple_action<observers::Actions::ContributeVolumeData>(
          local_observer, std::move(observation_id), subfile_path,
          array_component_id, std::move(element_volume_data), dependency);
    }
  }

  using observation_registration_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTagsList>
  std::optional<
      std::pair<observers::TypeOfObservation, observers::ObservationKey>>
  get_observation_type_and_key_for_registration(
      const db::DataBox<DbTagsList>& box) const {
    if (not active_block(db::get<domain::Tags::Domain<VolumeDim>>(box),
                         db::get<domain::Tags::Element<VolumeDim>>(box).id())) {
      return std::nullopt;
    }
    const std::optional<std::string> section_observation_key =
        observers::get_section_observation_key<ArraySectionIdTag>(box);
    if (not section_observation_key.has_value()) {
      return std::nullopt;
    }
    return {{observers::TypeOfObservation::Volume,
             observers::ObservationKey{
                 subfile_path_ + section_observation_key.value() + ".vol"}}};
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Event::pup(p);
    p | subfile_path_;
    p | variables_to_observe_;
    p | active_block_or_block_groups_;
    p | truncation_mesh_;
    p | dependency_;
  }

 private:
  template <typename Tag>
  static bool print_warning_about_optional() {
    Parallel::printf(
        "Warning: ObserveModalFields is trying to dump the tag %s "
        "but it is stored as a std::optional and has not been "
        "evaluated. This most commonly occurs when you are "
        "trying to either observe an analytic solution or errors when "
        "no analytic solution is available.\n",
        db::tag_name<Tag>());
    return false;
  }

  bool active_block(const Domain<VolumeDim>& domain,
                    const ElementId<VolumeDim>& element_id) const {
    if (not active_block_or_block_groups_.has_value()) {
      return true;
    }
    const std::unordered_set<std::string> block_names =
        domain::expand_block_groups_to_block_names(
            active_block_or_block_groups_.value(), domain.block_names(),
            domain.block_groups());
    return alg::found(block_names,
                      domain.blocks().at(element_id.block_id()).name());
  }

  std::string subfile_path_;
  std::unordered_set<std::string> variables_to_observe_{};
  std::optional<std::vector<std::string>> active_block_or_block_groups_{};
  std::optional<Mesh<VolumeDim>> truncation_mesh_{};
  std::optional<std::string> dependency_;
};

template <size_t VolumeDim, typename... Tensors,
          typename... NonTensorComputeTags, typename ArraySectionIdTag>
ObserveModalFields<VolumeDim, tmpl::list<Tensors...>,
                   tmpl::list<NonTensorComputeTags...>, ArraySectionIdTag>::
    ObserveModalFields(
        const std::string& subfile_name,
        const std::vector<std::string>& variables_to_observe,
        std::optional<std::vector<std::string>> active_block_or_block_groups,
        std::optional<Mesh<VolumeDim>> truncation_mesh,
        std::optional<std::string> dependency, const Options::Context& context)
    : subfile_path_("/" + subfile_name),
      variables_to_observe_([&context, &variables_to_observe]() {
        if (variables_to_observe.empty()) {
          PARSE_ERROR(context,
                      "ObserveModalFields requires at least one variable.");
        }
        std::unordered_set<std::string> result{};
        result.reserve(variables_to_observe.size() + 1);
        for (const auto& name : variables_to_observe) {
          result.insert(name);
        }
        return result;
      }()),
      active_block_or_block_groups_(std::move(active_block_or_block_groups)),
      truncation_mesh_(std::move(truncation_mesh)),
      dependency_(std::move(dependency)) {
  ASSERT(
      (... or (db::tag_name<Tensors>() == "InertialCoordinates")),
      "There is no tag with name 'InertialCoordinates' specified "
      "for the observer. Please make sure you specify a tag in the 'Tensors' "
      "list that has the 'db::tag_name()' 'InertialCoordinates'.");
  db::validate_selection<tmpl::list<Tensors...>>(variables_to_observe, context);
  variables_to_observe_.insert("InertialCoordinates");
}

/// \cond
template <size_t VolumeDim, typename... Tensors,
          typename... NonTensorComputeTags, typename ArraySectionIdTag>
PUP::able::PUP_ID ObserveModalFields<VolumeDim, tmpl::list<Tensors...>,
                                     tmpl::list<NonTensorComputeTags...>,
                                     ArraySectionIdTag>::my_PUP_ID =
    0;  // NOLINT
/// \endcond
}  // namespace Events
}  // namespace dg
