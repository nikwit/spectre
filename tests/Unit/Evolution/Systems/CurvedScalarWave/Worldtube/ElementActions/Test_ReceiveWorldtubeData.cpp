// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <random>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Block.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/ReceiveWorldtubeData.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/InitializeElementFacesGridCoordinates.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/SendToElements.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonChare.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {

template <typename Metavariables>
struct MockElementArray {
  using metavariables = Metavariables;
  static constexpr size_t Dim = metavariables::volume_dim;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<Dim>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              db::AddSimpleTags<
                  domain::Tags::Element<Dim>, domain::Tags::Mesh<Dim>,
                  Tags::GeodesicPunctureField<Dim>,
                  Tags::IteratedPunctureField<Dim>,
                  gr::Tags::Shift<DataVector, Dim>, gr::Tags::Lapse<DataVector>,
                  ::Tags::TimeStepId, Tags::WorldtubeSolution<Dim>,
                  Tags::FaceCoordinates<Dim, Frame::Inertial, true>>,
              db::AddComputeTags<>>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             tmpl::list<Actions::ReceiveWorldtubeData>>>;
};
template <typename Metavariables>
struct MockWorldtubeSingleton {
  using metavariables = Metavariables;
  static constexpr size_t Dim = metavariables::volume_dim;
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              db::AddSimpleTags<
                  Tags::ElementFacesGridCoordinates<Dim>, ::Tags::TimeStepId,
                  Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim,
                                       Frame::Inertial>,
                  Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim,
                                       Frame::Inertial>,
                  Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim,
                                       Frame::Inertial>,
                  Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                                       Frame::Inertial>,
                  Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim,
                                       Frame::Inertial>,
                  Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 2, Dim,
                                       Frame::Inertial>,
                  Tags::Psi0, Tags::dtPsi0,
                  Tags::ParticlePositionVelocity<Dim>>,
              db::AddComputeTags<>>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<Actions::SendToElements<Metavariables>>>>;
  using component_being_mocked = WorldtubeSingleton<Metavariables>;
};

template <size_t Dim>
struct MockMetavariables {
  static constexpr size_t volume_dim = Dim;
  using dg_element_array = MockElementArray<MockMetavariables>;
  using component_list = tmpl::list<MockWorldtubeSingleton<MockMetavariables>,
                                    MockElementArray<MockMetavariables>>;
  using const_global_cache_tags =
      tmpl::list<Tags::ExcisionSphere<Dim>, Tags::ExpansionOrder,
                 Tags::MaxIterations, Tags::WorldtubeRadius>;
};

// This test checks that `SendToElements` sends the coefficients of the Taylor
// series of the regular field and its time derivative to the elements
// abutting the worldtube and that `ReceiveWorldtubeData` evaluates the series
// at the centered inertial coordinates of the element faces and assembles the
// worldtube solution from it: the puncture field (set to zero here) is added
// and the time derivative is converted to Pi with lapse and shift (set to
// Minkowski here). At second order, the constant coefficient is the
// separately evolved `Tags::Psi0`, the constant coefficient of the
// time-derivative field is corrected by the motion of the expansion center
// (checked here with a nonzero particle velocity) and the trace of the
// second-order coefficient is reconstructed from the boundary monopole and
// added to the sent quadrupole coefficients.
SPECTRE_TEST_CASE("Unit.CurvedScalarWave.Worldtube.ReceiveWorldtubeData",
                  "[Unit]") {
  static constexpr size_t Dim = 3;
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> dist(-10., 10.);
  using metavars = MockMetavariables<Dim>;
  domain::creators::register_derived_with_charm();
  using element_chare = MockElementArray<metavars>;
  using worldtube_chare = MockWorldtubeSingleton<metavars>;
  const size_t initial_extent = 5;
  const size_t face_size = initial_extent * initial_extent;
  const auto basis = Spectral::Basis::Legendre;
  const auto quadrature = Spectral::Quadrature::GaussLobatto;
  // we create several differently refined shells so a different number of
  // elements sends data
  for (const auto& [expansion_order, initial_refinement, worldtube_radius] :
       cartesian_product(std::array<size_t, 3>{0, 1, 2},
                         std::array<size_t, 3>{0, 1, 2},
                         make_array(0.07, 1., 2.8))) {
    CAPTURE(expansion_order);
    CAPTURE(worldtube_radius);
    CAPTURE(initial_refinement);
    const domain::creators::Sphere shell{worldtube_radius,
                                         3.,
                                         domain::creators::Sphere::Excision{},
                                         initial_refinement,
                                         initial_extent,
                                         true};
    const auto shell_domain = shell.create_domain();
    const auto excision_sphere =
        shell_domain.excision_spheres().at("ExcisionSphere");

    const auto& initial_refinements = shell.initial_refinement_levels();
    const auto& initial_extents = shell.initial_extents();
    ActionTesting::MockRuntimeSystem<metavars> runner{
        {excision_sphere, expansion_order, static_cast<size_t>(0),
         excision_sphere.radius()}};
    const auto element_ids = initial_element_ids(initial_refinements);
    const auto& blocks = shell_domain.blocks();

    using puncture_field_type =
        Variables<tmpl::list<CurvedScalarWave::Tags::Psi,
                             ::Tags::dt<CurvedScalarWave::Tags::Psi>,
                             ::Tags::deriv<CurvedScalarWave::Tags::Psi,
                                           tmpl::size_t<3>, Frame::Inertial>>>;

    // The puncture field will get added to the regular field. Here, we set
    // the puncture field to 0, so the evaluated Taylor series is passed on
    // directly and we can check the analytical result.
    const puncture_field_type puncture_field{face_size, 0.};
    const double psi_monopole = dist(generator);
    const double dt_psi_monopole = dist(generator);
    const auto psi_dipole =
        make_with_random_values<tnsr::i<double, Dim, Frame::Inertial>>(
            make_not_null(&generator), dist, 0.);
    const auto dt_psi_dipole =
        make_with_random_values<tnsr::i<double, Dim, Frame::Inertial>>(
            make_not_null(&generator), dist, 0.);
    const auto psi_quadrupole =
        make_with_random_values<tnsr::ii<double, Dim, Frame::Inertial>>(
            make_not_null(&generator), dist, 0.);
    const auto dt_psi_quadrupole =
        make_with_random_values<tnsr::ii<double, Dim, Frame::Inertial>>(
            make_not_null(&generator), dist, 0.);
    // at second order the constant coefficient of the regular field and its
    // time derivative are evolved separately by the worldtube
    const double psi_0 = dist(generator);
    const double dt_psi_0 = dist(generator);
    // the particle sits at the center of the excision sphere so the centered
    // inertial face coordinates are equal to the grid coordinates of the
    // faces. The velocity is nonzero to check the correction to the constant
    // coefficient of the time-derivative field from the motion of the
    // expansion center.
    const tnsr::I<double, Dim> particle_position(0.);
    const auto particle_velocity =
        make_with_random_values<tnsr::I<double, Dim>>(make_not_null(&generator),
                                                      dist, 0.);
    const std::array<tnsr::I<double, Dim>, 2> particle_pos_vel{
        {particle_position, particle_velocity}};

    // The Taylor coefficients that `SendToElements` should assemble and send,
    // ordered as [monopole, dipole x/y/z, quadrupole xx, xy, xz, yy, yz, zz],
    // where the quadrupole coefficients hold the independent components of
    // the full (trace-included) symmetric second-order coefficient.
    const size_t num_coefs =
        expansion_order == 0 ? 1 : (expansion_order == 1 ? 4 : 10);
    DataVector expected_psi_coefs(num_coefs);
    DataVector expected_dt_psi_coefs(num_coefs);
    expected_psi_coefs[0] = psi_monopole;
    expected_dt_psi_coefs[0] = dt_psi_monopole;
    if (expansion_order > 0) {
      for (size_t i = 0; i < Dim; ++i) {
        expected_psi_coefs[i + 1] = psi_dipole.get(i);
        expected_dt_psi_coefs[i + 1] = dt_psi_dipole.get(i);
      }
    }
    if (expansion_order > 1) {
      // the constant coefficient is evolved separately at second order
      expected_psi_coefs[0] = psi_0;
      // the constant coefficient of the time-derivative field gets a
      // correction from the motion of the expansion center:
      // (dt Psi)_0 = dtPsi0 - Psi_i * v^i
      expected_dt_psi_coefs[0] = dt_psi_0;
      for (size_t i = 0; i < Dim; ++i) {
        expected_dt_psi_coefs[0] -=
            psi_dipole.get(i) * particle_velocity.get(i);
      }
      // the traces of the second-order coefficients are reconstructed from
      // the boundary monopoles and added back to the trace-free quadrupoles
      const double trace_psi_over_3 =
          (psi_monopole - expected_psi_coefs[0]) / square(worldtube_radius);
      const double trace_dt_psi_over_3 =
          (dt_psi_monopole - expected_dt_psi_coefs[0]) /
          square(worldtube_radius);
      size_t index = 4;
      for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = i; j < Dim; ++j, ++index) {
          expected_psi_coefs[index] =
              psi_quadrupole.get(i, j) + (i == j ? trace_psi_over_3 : 0.);
          expected_dt_psi_coefs[index] =
              dt_psi_quadrupole.get(i, j) + (i == j ? trace_dt_psi_over_3 : 0.);
        }
      }
    }
    // the full (trace-included) second-order coefficients used to evaluate
    // the expected Taylor series
    tnsr::ii<double, Dim, Frame::Inertial> psi_coef_2_full{};
    tnsr::ii<double, Dim, Frame::Inertial> dt_psi_coef_2_full{};
    if (expansion_order > 1) {
      size_t index = 4;
      for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = i; j < Dim; ++j, ++index) {
          psi_coef_2_full.get(i, j) = expected_psi_coefs[index];
          dt_psi_coef_2_full.get(i, j) = expected_dt_psi_coefs[index];
        }
      }
    }
    const Time dummy_time{{1., 2.}, {1, 2}};
    const TimeStepId dummy_time_step_id{true, 123, dummy_time};
    std::unordered_map<ElementId<Dim>, tnsr::I<DataVector, Dim, Frame::Grid>>
        element_faces_grid_coords{};
    Initialization::InitializeElementFacesGridCoordinates<Dim>::apply(
        make_not_null(&element_faces_grid_coords), initial_extents,
        initial_refinements, quadrature, shell_domain, excision_sphere);
    for (const auto& element_id : element_ids) {
      auto element = domain::create_initial_element(element_id, blocks,
                                                    initial_refinements);
      auto mesh = domain::create_initial_mesh(initial_extents, element, basis,
                                              quadrature);
      const size_t grid_size = mesh.number_of_grid_points();

      // we set lapse and shift to Minkowski so dt Psi = - Pi
      Scalar<DataVector> lapse(grid_size, 1.);
      tnsr::I<DataVector, Dim, Frame::Inertial> shift(grid_size, 0.);
      // not initialized, the action does that for us
      typename CurvedScalarWave::System<Dim>::variables_tag::type
          worldtube_solution{};
      std::optional<puncture_field_type> optional_puncture_field =
          excision_sphere.abutting_direction(element_id).has_value()
              ? std::make_optional<puncture_field_type>(puncture_field)
              : std::nullopt;
      std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>
          inertial_face_coords{};

      if (element_faces_grid_coords.count(element_id) > 0) {
        const auto& grid_face_coords = element_faces_grid_coords.at(element_id);
        inertial_face_coords.emplace(get<0>(grid_face_coords).size());
        for (size_t i = 0; i < Dim; ++i) {
          inertial_face_coords.value().get(i) = grid_face_coords.get(i);
        }
      }
      ActionTesting::emplace_array_component_and_initialize<element_chare>(
          &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
          element_id,
          {std::move(element), std::move(mesh),
           std::move(optional_puncture_field), std::nullopt, std::move(shift),
           std::move(lapse), dummy_time_step_id, worldtube_solution,
           inertial_face_coords});
    }

    ActionTesting::emplace_singleton_component_and_initialize<worldtube_chare>(
        &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
        {element_faces_grid_coords, dummy_time_step_id,
         Scalar<double>(psi_monopole), Scalar<double>(dt_psi_monopole),
         psi_dipole, dt_psi_dipole, psi_quadrupole, dt_psi_quadrupole,
         Scalar<DataVector>(DataVector(1, psi_0)),
         Scalar<DataVector>(DataVector(1, dt_psi_0)), particle_pos_vel});

    ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

    // If it is an abutting element, we don't have data yet and should not be
    // ready. Else, there is nothing to do so it's ready.
    for (const auto& element_id : element_ids) {
      if (element_faces_grid_coords.count(element_id)) {
        CHECK(not ActionTesting::next_action_if_ready<element_chare>(
            make_not_null(&runner), element_id));
      } else {
        CHECK(ActionTesting::next_action_if_ready<element_chare>(
            make_not_null(&runner), element_id));
      }
    }
    // SendToElements
    ActionTesting::next_action<worldtube_chare>(make_not_null(&runner), 0);

    for (const auto& element_id : element_ids) {
      using inbox_tag = Tags::RegularFieldInbox<Dim>;
      const auto& element_inbox =
          ActionTesting::get_inbox_tag<element_chare, inbox_tag>(runner,
                                                                 element_id);
      if (element_faces_grid_coords.count(element_id)) {
        // the Taylor series and its coefficients are assembled in a
        // different operation order than the expected values, so roundoff
        // accumulates slightly differently
        Approx local_approx = Approx::custom().epsilon(1.e-10).scale(1.);
        CHECK(element_inbox.count(dummy_time_step_id));
        const auto& inbox_data = element_inbox.at(dummy_time_step_id);
        CHECK_ITERABLE_CUSTOM_APPROX(
            get(get<CurvedScalarWave::Tags::Psi>(inbox_data)),
            expected_psi_coefs, local_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(inbox_data)),
            expected_dt_psi_coefs, local_approx);

        const auto& inertial_coords =
            ActionTesting::get_databox_tag<
                element_chare,
                Tags::FaceCoordinates<Dim, Frame::Inertial, true>>(runner,
                                                                   element_id)
                .value();
        DataVector expected_solution_psi(face_size, expected_psi_coefs[0]);
        DataVector expected_solution_dt_psi(face_size,
                                            expected_dt_psi_coefs[0]);
        if (expansion_order > 0) {
          for (size_t i = 0; i < Dim; ++i) {
            expected_solution_psi +=
                expected_psi_coefs[i + 1] * inertial_coords.get(i);
            expected_solution_dt_psi +=
                expected_dt_psi_coefs[i + 1] * inertial_coords.get(i);
          }
        }
        if (expansion_order > 1) {
          for (size_t i = 0; i < Dim; ++i) {
            for (size_t j = 0; j < Dim; ++j) {
              expected_solution_psi += psi_coef_2_full.get(i, j) *
                                       inertial_coords.get(i) *
                                       inertial_coords.get(j);
              expected_solution_dt_psi += dt_psi_coef_2_full.get(i, j) *
                                          inertial_coords.get(i) *
                                          inertial_coords.get(j);
            }
          }
        }

        // ReceiveWorldtubeData
        CHECK(ActionTesting::next_action_if_ready<element_chare>(
            make_not_null(&runner), element_id));
        CHECK(element_inbox.empty());
        const auto& worldtube_solution =
            ActionTesting::get_databox_tag<element_chare,
                                           Tags::WorldtubeSolution<Dim>>(
                runner, element_id);
        CHECK_ITERABLE_CUSTOM_APPROX(
            get(get<CurvedScalarWave::Tags::Psi>(worldtube_solution)),
            expected_solution_psi, local_approx);
        CHECK_ITERABLE_CUSTOM_APPROX(
            get(get<CurvedScalarWave::Tags::Pi>(worldtube_solution)),
            -expected_solution_dt_psi, local_approx);
        for (size_t i = 0; i < Dim; ++i) {
          DataVector expected_di_psi(
              face_size, expansion_order > 0 ? expected_psi_coefs[i + 1] : 0.);
          if (expansion_order > 1) {
            for (size_t j = 0; j < Dim; ++j) {
              expected_di_psi +=
                  2. * psi_coef_2_full.get(i, j) * inertial_coords.get(j);
            }
          }
          CHECK_ITERABLE_CUSTOM_APPROX(
              get<CurvedScalarWave::Tags::Phi<Dim>>(worldtube_solution).get(i),
              expected_di_psi, local_approx);
        }
      } else {
        CHECK(element_inbox.empty());
      }
    }
  }
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
