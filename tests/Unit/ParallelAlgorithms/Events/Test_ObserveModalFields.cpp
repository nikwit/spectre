// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/MetavariablesTag.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Creators/Rectilinear.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/ParallelAlgorithms/Events/ObserveFields.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Events/ObserveModalFields.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {

namespace observe_fields = TestHelpers::dg::Events::ObserveFields;

template <typename System, typename = void>
struct prim_vars_impl {
  using type = tmpl::list<>;
  static constexpr bool has_prims = false;
};

template <typename System>
struct prim_vars_impl<System,
                      std::void_t<typename System::primitive_variables_tag>> {
  using type = typename System::primitive_variables_tag::tags_list;
  static constexpr bool has_prims = true;
};

template <typename System>
using prim_vars_list = typename prim_vars_impl<System>::type;

template <typename System, typename ObserveEvent>
void test_modal_observe(
    const ObserveEvent& observe,
    const std::optional<Mesh<System::volume_dim>>& truncation_mesh) {
  using metavariables = observe_fields::Metavariables<System, false>;
  constexpr size_t volume_dim = System::volume_dim;
  using DataType = typename System::data_type;
  using element_component = observe_fields::ElementComponent<metavariables>;
  using observer_component =
      observe_fields::MockObserverComponent<metavariables>;
  using coordinates_tag =
      domain::Tags::Coordinates<volume_dim, Frame::Inertial>;

  const ElementId<volume_dim> element_id(0);
  const Element<volume_dim> element(element_id, {});
  const domain::creators::Rectilinear<volume_dim> rectilinear{
      make_array<volume_dim>(-2.0), make_array<volume_dim>(2.0),
      make_array<volume_dim>(0_st), make_array<volume_dim>(5_st),
      make_array<volume_dim>(true)};
  const Mesh<volume_dim> mesh(5, Spectral::Basis::Legendre,
                              Spectral::Quadrature::GaussLobatto);

  const double observation_time = 2.0;
  Variables<typename System::variables_tag::tags_list> vars(
      mesh.number_of_grid_points());
  Variables<tmpl::list<coordinates_tag>> coordinate_vars(
      mesh.number_of_grid_points());
  Variables<prim_vars_list<System>> prim_vars(mesh.number_of_grid_points());
  std::iota(vars.data(), vars.data() + vars.size(), 1.0);
  std::iota(coordinate_vars.data(),
            coordinate_vars.data() + coordinate_vars.size(),
            static_cast<double>(vars.size()));
  if constexpr (not std::is_same_v<tmpl::list<>, prim_vars_list<System>>) {
    std::iota(prim_vars.data(), prim_vars.data() + prim_vars.size(),
              static_cast<double>(vars.size() + coordinate_vars.size()));
  }

  const typename System::solution_for_test analytic_solution{};
  using solution_variables = typename System::solution_for_test::vars_for_test;
  using analytic_solution_variables =
      db::wrap_tags_in<::Tags::detail::AnalyticImpl, solution_variables>;

  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  MockRuntimeSystem runner(
      tuples::TaggedTuple<
          ::Tags::AnalyticSolution<typename System::solution_for_test>>{
          std::move(analytic_solution)});
  ActionTesting::emplace_component<element_component>(make_not_null(&runner),
                                                      element_id);
  ActionTesting::emplace_group_component<observer_component>(&runner);

  auto box = db::create<db::AddSimpleTags<
      Parallel::Tags::MetavariablesImpl<metavariables>,
      domain::Tags::Domain<volume_dim>, domain::Tags::Element<volume_dim>,
      domain::Tags::Mesh<volume_dim>,
      ::Tags::Variables<typename decltype(vars)::tags_list>,
      ::Tags::Variables<typename decltype(prim_vars)::tags_list>,
      coordinates_tag, ::Tags::AnalyticSolutions<solution_variables>,
      observers::Tags::ObservationKey<void>>>(
      metavariables{}, rectilinear.create_domain(), element, mesh, vars,
      prim_vars, get<coordinates_tag>(coordinate_vars),
      std::optional<Variables<analytic_solution_variables>>{},
      std::optional<std::string>{});

  observe_fields::MockContributeVolumeData::results = {};

  auto obs_box = make_observation_box<tmpl::push_back<
      tmpl::filter<typename ObserveEvent::compute_tags_for_observation_box,
                   db::is_compute_tag<tmpl::_1>>,
      ::Events::Tags::ObserverMeshCompute<volume_dim>>>(make_not_null(&box));

  const Event::ObservationValue observation_value{"TestObservation",
                                                  observation_time};

  observe(obs_box, mesh,
          ActionTesting::cache<element_component>(runner, element_id),
          element_id, static_cast<const element_component*>(nullptr),
          observation_value);
  runner.template invoke_queued_simple_action<observer_component>(0);

  const auto& results = observe_fields::MockContributeVolumeData::results;
  const Mesh<volume_dim>& mesh_for_output =
      truncation_mesh.has_value() ? truncation_mesh.value() : mesh;
  const auto mesh_extents = mesh_for_output.extents();
  const std::vector<size_t> expected_extents{mesh_extents.begin(),
                                             mesh_extents.end()};
  CHECK(results.received_volume_data.extents == expected_extents);

  const auto truncate_modal = [&mesh, &mesh_for_output](const ModalVector& in) {
    if (mesh.number_of_grid_points() ==
        mesh_for_output.number_of_grid_points()) {
      return ModalVector{in};
    }
    ModalVector out(mesh_for_output.number_of_grid_points());
    const Index<volume_dim> source_extents(mesh.extents());
    const Index<volume_dim> target_extents(mesh_for_output.extents());
    for (size_t target_linear = 0;
         target_linear < mesh_for_output.number_of_grid_points();
         ++target_linear) {
      const Index<volume_dim> target_multi =
          expanded_index(target_linear, target_extents);
      Index<volume_dim> source_multi{0};
      for (size_t d = 0; d < volume_dim; ++d) {
        source_multi[d] = target_multi[d];
      }
      const size_t source_linear =
          collapsed_index(source_multi, source_extents);
      out[target_linear] = in[source_linear];
    }
    return out;
  };

  const auto nodal_to_modal_data = [&mesh,
                                    &truncate_modal](const DataVector& nodal) {
    ModalVector modal = to_modal_coefficients(nodal, mesh);
    ModalVector truncated = truncate_modal(modal);
    DataVector result(truncated.size());
    for (size_t i = 0; i < truncated.size(); ++i) {
      result[i] = truncated[i];
    }
    return result;
  };

  std::unordered_map<std::string, DataVector> expected_components{};
  const auto record_component = [&expected_components, &nodal_to_modal_data](
                                    const std::string& name,
                                    const DataVector& nodal) {
    expected_components[name] = nodal_to_modal_data(nodal);
  };
  for (size_t i = 0; i < volume_dim; ++i) {
    record_component(
        "InertialCoordinates_" + std::string(1, gsl::at({'x', 'y', 'z'}, i)),
        get<coordinates_tag>(coordinate_vars).get(i));
  }

  try {
    System::check_data([&record_component, &prim_vars, &vars](
                           const std::string& name, auto tag_v,
                           const auto... indices) {
      using tag = decltype(tag_v);
      if constexpr (std::is_same_v<tag, observe_fields::Tags::ScalarVarTimesTwo<
                                            DataType>>) {
        record_component(
            name, DataType{2.0 * get<typename System::ScalarVar>(vars).get()});
      } else if constexpr (std::is_same_v<
                               tag, observe_fields::Tags::ScalarVarTimesThree<
                                        DataType>>) {
        record_component(
            name, DataType{3.0 * get<typename System::ScalarVar>(vars).get()});
      } else {
        if constexpr (tmpl::list_contains_v<
                          typename std::decay_t<decltype(prim_vars)>::tags_list,
                          tag>) {
          record_component(name, get<tag>(prim_vars).get(indices...));
        } else {
          record_component(name, get<tag>(vars).get(indices...));
        }
      }
    });
  } catch (const std::exception& e) {
    INFO("check_data threw exception: " << e.what());
    throw;
  }

  REQUIRE(results.received_volume_data.tensor_components.size() ==
          expected_components.size());
  for (const auto& component : results.received_volume_data.tensor_components) {
    const auto expected_it = expected_components.find(component.name);
    REQUIRE(expected_it != expected_components.end());
    DataVector observed;
    if (std::holds_alternative<DataVector>(component.data)) {
      observed = std::get<DataVector>(component.data);
    } else {
      const auto& as_float = std::get<std::vector<float>>(component.data);
      observed = DataVector(as_float.size());
      for (size_t i = 0; i < as_float.size(); ++i) {
        observed[i] = as_float[i];
      }
    }
    CHECK_ITERABLE_APPROX(observed, expected_it->second);
  }
}

template <typename System>
void test_modal_factory_creation() {
  const auto obs = TestHelpers::test_creation<
      std::unique_ptr<Event>, observe_fields::Metavariables<System, false>>(
      "ObserveModalFields:\n"
      "  SubfileName: element_data\n"
      "  VariablesToObserve: [Scalar, ScalarVarTimesTwo, ScalarVarTimesThree]\n"
      "  BlocksToObserve: All\n"
      "  InterpolateToMesh: None\n");
  namespace helper_tags = observe_fields::Tags;
  CHECK(obs != nullptr);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Events.ObserveModalFields",
                  "[Unit][ParallelAlgorithms]") {
  using System = observe_fields::ScalarSystem<dg::Events::ObserveModalFields>;

  const typename System::ObserveEvent observe{
      "element_data",
      std::vector<std::string>{"Scalar", "ScalarVarTimesTwo",
                               "ScalarVarTimesThree"},
      std::nullopt, std::nullopt, std::nullopt};

  test_modal_observe<System>(observe, std::nullopt);

  const Mesh<System::volume_dim> truncation_mesh(
      3, Spectral::Basis::Legendre, Spectral::Quadrature::GaussLobatto);
  const typename System::ObserveEvent observe_truncated{
      "element_data",
      std::vector<std::string>{"Scalar", "ScalarVarTimesTwo",
                               "ScalarVarTimesThree"},
      std::nullopt, truncation_mesh, std::nullopt};
  test_modal_observe<System>(observe_truncated, truncation_mesh);

  test_modal_factory_creation<System>();
}
