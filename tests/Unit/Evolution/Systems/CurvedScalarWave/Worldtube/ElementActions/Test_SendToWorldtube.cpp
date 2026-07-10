// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <random>
#include <unordered_map>
#include <vector>

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
#include "Domain/ElementMap.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/IteratePunctureField.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/ReceiveWorldtubeData.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/SendToWorldtube.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/InitializeElementFacesGridCoordinates.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/IterateAccelerationTerms.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/ReceiveElementData.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/UpdateAcceleration.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonChare.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/RealSphericalHarmonics.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Spherepack.hpp"
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
                  domain::Tags::Coordinates<Dim, Frame::Grid>,
                  Tags::GeodesicPunctureField<Dim>,
                  gr::Tags::Shift<DataVector, Dim>, gr::Tags::Lapse<DataVector>,
                  domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                                Frame::Inertial>,
                  domain::Tags::InverseJacobian<Dim, Frame::Grid,
                                                Frame::Inertial>,
                  typename CurvedScalarWave::System<Dim>::variables_tag,
                  domain::Tags::MeshVelocity<Dim>, ::Tags::TimeStepId,
                  Tags::ParticlePositionVelocity<Dim>, Tags::CurrentIteration,
                  Tags::FaceCoordinates<Dim, Frame::Inertial, true>>,
              db::AddComputeTags<
                  Tags::FaceCoordinatesCompute<Dim, Frame::Grid, true>,
                  Tags::FaceQuantitiesCompute>>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<Actions::SendToWorldtube, Actions::IteratePunctureField,
                     CurvedScalarWave::Worldtube::Actions::
                         ReceiveWorldtubeData>>>;
};
// Mock element for the spherical-harmonic shell case where a single element
// covers the entire worldtube boundary. The face quantities and face
// coordinates are simple tags here so the test can fill them directly with
// analytic data on the spherical-harmonic collocation grid.
template <typename Metavariables>
struct MockElementArraySphericalShell {
  using metavariables = Metavariables;
  static constexpr size_t Dim = metavariables::volume_dim;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<Dim>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              db::AddSimpleTags<domain::Tags::Element<Dim>,
                                domain::Tags::Mesh<Dim>, Tags::FaceQuantities,
                                Tags::FaceCoordinates<Dim, Frame::Grid, true>,
                                Tags::GeodesicPunctureField<Dim>,
                                Tags::IteratedPunctureField<Dim>,
                                CurvedScalarWave::Tags::Phi<Dim>,
                                domain::Tags::MeshVelocity<Dim>,
                                domain::Tags::InverseJacobian<Dim, Frame::Grid,
                                                              Frame::Inertial>,
                                gr::Tags::Shift<DataVector, Dim>,
                                gr::Tags::Lapse<DataVector>, ::Tags::TimeStepId,
                                Tags::CurrentIteration>,
              db::AddComputeTags<>>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             tmpl::list<Actions::SendToWorldtube,
                                        CurvedScalarWave::Worldtube::Actions::
                                            ReceiveWorldtubeData>>>;
};

template <typename Metavariables>
struct MockWorldtubeSingleton {
  using metavariables = Metavariables;
  static constexpr size_t Dim = metavariables::volume_dim;
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::EvolvedPosition<Dim>, Tags::EvolvedVelocity<Dim>,
                 Tags::Psi0, Tags::dtPsi0>>;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              db::AddSimpleTags<
                  Tags::ElementFacesGridCoordinates<Dim>, ::Tags::TimeStepId,
                  Tags::CurrentIteration, Tags::GeodesicAcceleration<Dim>,
                  CurvedScalarWave::Worldtube::Tags::ParticlePositionVelocity<
                      Dim>,
                  Tags::BackgroundQuantities<Dim>, variables_tag,
                  dt_variables_tag,
                  gr::Tags::InverseSpacetimeMetric<double, Dim, Frame::Grid>,
                  gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim,
                                                                Frame::Grid>>,
              db::AddComputeTags<>>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<Actions::ReceiveElementData,
                     ::Actions::MutateApply<IterateAccelerationTerms>,
                     Actions::SendAccelerationTerms<Metavariables>,
                     ::Actions::MutateApply<UpdateAcceleration>>>>;
  using component_being_mocked = WorldtubeSingleton<Metavariables>;
};

template <size_t Dim>
struct MockMetavariables {
  static constexpr size_t volume_dim = Dim;

  using component_list =
      tmpl::list<MockWorldtubeSingleton<MockMetavariables>,
                 MockElementArray<MockMetavariables>,
                 MockElementArraySphericalShell<MockMetavariables>>;
  using dg_element_array = MockElementArray<MockMetavariables>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Dim>, Tags::ExcisionSphere<Dim>,
                 Tags::WorldtubeRadius, Tags::PunctureFieldConfig,
                 Tags::ExpansionOrder, Tags::MaxIterations, Tags::Charge,
                 Tags::Mass, ::Tags::Time, Tags::SelfForceTurnOnTime,
                 Tags::SelfForceTurnOnInterval>;
};

// Checks the projection onto spherical harmonics for a face mesh with
// `Spectral::Basis::SphericalHarmonic`, where a single spherical-shell element
// covers the entire worldtube boundary and the projection integrals are
// evaluated with the Gauss quadrature of the spherical-harmonic collocation
// grid. The regular field and its time derivative are set to linear
// combinations of real spherical harmonics with random coefficients for
// l <= 2, plus a high-l mode below the exactness limit of the quadrature
// which must not alias into the projected coefficients. Because
// `ylm::real_spherical_harmonic` is orthonormal over the unit sphere, the
// sent coefficients must equal the input coefficients times the squared
// worldtube radius, which is the same convention as the wedge path tested
// below (the worldtube singleton divides by the squared radius). The
// Euclidean area element of `Tags::FaceQuantities` is filled with garbage to
// check that it is unused on this code path, and the mesh velocity is zero
// so the advective term vanishes.
void test_s2_projection() {
  static constexpr size_t Dim = 3;
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> dist(-10., 10.);
  using metavars = MockMetavariables<Dim>;
  using element_chare = MockElementArraySphericalShell<metavars>;
  using worldtube_chare = MockWorldtubeSingleton<metavars>;
  using puncture_field_type =
      Variables<tmpl::list<CurvedScalarWave::Tags::Psi,
                           ::Tags::dt<CurvedScalarWave::Tags::Psi>,
                           ::Tags::deriv<CurvedScalarWave::Tags::Psi,
                                         tmpl::size_t<3>, Frame::Inertial>>>;
  const size_t l_max = 7;
  const size_t n_theta = l_max + 1;
  const size_t n_phi = 2 * l_max + 1;
  const size_t n_r = 4;
  const Mesh<Dim> mesh{
      {{n_r, n_theta, n_phi}},
      {{Spectral::Basis::Legendre, Spectral::Basis::SphericalHarmonic,
        Spectral::Basis::SphericalHarmonic}},
      {{Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::Gauss,
        Spectral::Quadrature::Equiangular}}};
  const size_t face_size = n_theta * n_phi;
  const size_t grid_size = mesh.number_of_grid_points();

  // The collocation angles of the spherical-harmonic grid: Gauss-Legendre
  // points in theta from SPHEREPACK and equiangular points in phi. Theta
  // varies fastest in the face DataVector.
  std::vector<double> gauss_points(n_theta + 1);
  std::vector<double> gauss_weights(n_theta + 1);
  std::vector<double> gaqd_work(n_theta);
  int gaqd_err = 0;
  gaqd_(static_cast<int>(n_theta), gauss_points.data(), gauss_weights.data(),
        gaqd_work.data(), static_cast<int>(gauss_weights.size()), &gaqd_err);
  REQUIRE(gaqd_err == 0);
  DataVector theta(face_size);
  DataVector phi(face_size);
  for (size_t j = 0; j < n_phi; ++j) {
    for (size_t i = 0; i < n_theta; ++i) {
      theta[i + j * n_theta] = gauss_points[i];
      phi[i + j * n_theta] =
          2. * M_PI * static_cast<double>(j) / static_cast<double>(n_phi);
    }
  }

  // unused but the tag is needed to compile
  const double time = std::numeric_limits<double>::signaling_NaN();
  const Time dummy_time{{1., 2.}, {1, 2}};
  const TimeStepId dummy_time_step_id{true, 123, dummy_time};
  const ElementId<Dim> element_id{0};
  for (const auto& [expansion_order, worldtube_radius] : cartesian_product(
           std::array<size_t, 3>{0, 1, 2}, make_array(0.07, 1.6))) {
    CAPTURE(expansion_order);
    CAPTURE(worldtube_radius);
    const ExcisionSphere<Dim> excision_sphere{
        worldtube_radius,
        tnsr::I<double, Dim, Frame::Grid>{{0., 0., 0.}},
        {{0, Direction<Dim>::lower_xi()}}};

    std::array<double, 9> psi_coefs{};
    std::array<double, 9> dt_psi_coefs{};
    DataVector psi_face(face_size, 0.);
    DataVector dt_psi_face(face_size, 0.);
    {
      size_t index = 0;
      for (size_t l = 0; l <= 2; ++l) {
        for (int m = -static_cast<int>(l); m <= static_cast<int>(l);
             ++m, ++index) {
          gsl::at(psi_coefs, index) = dist(generator);
          gsl::at(dt_psi_coefs, index) = dist(generator);
          const DataVector harmonic =
              ylm::real_spherical_harmonic(theta, phi, l, m);
          psi_face += gsl::at(psi_coefs, index) * harmonic;
          dt_psi_face += gsl::at(dt_psi_coefs, index) * harmonic;
        }
      }
    }
    // Contaminate the fields with high-l modes below the exactness limit of
    // the quadrature. They are orthogonal to the projected harmonics, so
    // they must not alias into the sent coefficients.
    psi_face += dist(generator) *
                ylm::real_spherical_harmonic(theta, phi, l_max - 1, 3);
    dt_psi_face += dist(generator) *
                   ylm::real_spherical_harmonic(theta, phi, l_max - 1, -4);

    // The puncture field is set to zero so psi and dt_psi are projected
    // directly and we can check the analytical result.
    typename Tags::FaceQuantities::type face_quantities{};
    face_quantities.emplace(face_size);
    get(get<CurvedScalarWave::Tags::Psi>(*face_quantities)) = psi_face;
    get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*face_quantities)) =
        dt_psi_face;
    // the Euclidean area element must be unused when projecting on the
    // spherical-harmonic grid, so we fill it with garbage
    get(get<gr::surfaces::Tags::AreaElement<DataVector>>(*face_quantities)) =
        1.2345e6;

    // the particle is at the center of the excision sphere, so the centered
    // face coordinates are just the coordinates of the worldtube sphere
    tnsr::I<DataVector, Dim, Frame::Grid> face_coords(face_size);
    get<0>(face_coords) = worldtube_radius * sin(theta) * cos(phi);
    get<1>(face_coords) = worldtube_radius * sin(theta) * sin(phi);
    get<2>(face_coords) = worldtube_radius * cos(theta);

    const auto puncture_field_options =
        CurvedScalarWave::Worldtube::PunctureField{
            CurvedScalarWave::Worldtube::PunctureField::Schwarzschild{
                expansion_order, 1.}};
    tuples::TaggedTuple<domain::Tags::Domain<Dim>, Tags::ExcisionSphere<Dim>,
                        Tags::WorldtubeRadius, Tags::PunctureFieldConfig,
                        Tags::ExpansionOrder, Tags::MaxIterations, Tags::Charge,
                        Tags::Mass, ::Tags::Time, Tags::SelfForceTurnOnTime,
                        Tags::SelfForceTurnOnInterval>
        tuple_of_opts{Domain<Dim>{},
                      excision_sphere,
                      worldtube_radius,
                      puncture_field_options,
                      expansion_order,
                      static_cast<size_t>(1),
                      0.1,
                      std::nullopt,
                      time,
                      std::nullopt,
                      std::nullopt};
    ActionTesting::MockRuntimeSystem<metavars> runner{std::move(tuple_of_opts)};

    // we set the mesh velocity to zero so the advective term vanishes and
    // we can recover the exact time derivative
    ActionTesting::emplace_array_component_and_initialize<element_chare>(
        &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
        element_id,
        {Element<Dim>{element_id, {}}, mesh, std::move(face_quantities),
         std::make_optional(face_coords),
         std::make_optional(puncture_field_type{face_size, 0.}),
         std::optional<puncture_field_type>{},
         tnsr::i<DataVector, Dim>{grid_size, 0.},
         std::make_optional(
             tnsr::I<DataVector, Dim, Frame::Inertial>{grid_size, 0.}),
         InverseJacobian<DataVector, Dim, Frame::Grid, Frame::Inertial>{
             grid_size, 0.},
         tnsr::I<DataVector, Dim, Frame::Inertial>{grid_size, 0.},
         Scalar<DataVector>{grid_size, 1.}, dummy_time_step_id,
         static_cast<size_t>(0)});

    // these are all unused
    tnsr::I<double, Dim> particle_position(0.);
    auto particle_velocity = particle_position;
    const std::array<tnsr::I<double, Dim>, 2> particle_pos_vel{
        {std::move(particle_position), std::move(particle_velocity)}};
    tuples::TaggedTuple<
        gr::Tags::SpacetimeMetric<double, Dim>,
        gr::Tags::InverseSpacetimeMetric<double, Dim>,
        gr::Tags::SpacetimeChristoffelSecondKind<double, Dim>,
        gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim>,
        Tags::TimeDilationFactor>
        background_quantities{};
    typename MockWorldtubeSingleton<metavars>::variables_tag::type
        evolved_vars_wt{};
    typename MockWorldtubeSingleton<metavars>::dt_variables_tag::type
        dt_variables{};
    tnsr::I<double, Dim> geodesic_acc{};
    tnsr::AA<double, Dim, Frame::Grid> inverse_spacetime_metric{};
    tnsr::A<double, Dim, Frame::Grid> trace_spacetime_christoffel{};
    ActionTesting::emplace_singleton_component_and_initialize<worldtube_chare>(
        &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
        {std::unordered_map<ElementId<Dim>,
                            tnsr::I<DataVector, Dim, Frame::Grid>>{},
         dummy_time_step_id, static_cast<size_t>(0), geodesic_acc,
         particle_pos_vel, background_quantities, evolved_vars_wt, dt_variables,
         inverse_spacetime_metric, trace_spacetime_christoffel});

    ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
    ActionTesting::next_action<element_chare>(make_not_null(&runner),
                                              element_id);

    // the mesh velocity is zero, so the advective term should vanish
    const auto& regular_field_advective_term =
        ActionTesting::get_databox_tag<element_chare,
                                       Tags::RegularFieldAdvectiveTerm<Dim>>(
            runner, element_id);
    CHECK_ITERABLE_APPROX(regular_field_advective_term,
                          Scalar<DataVector>(face_size, 0.));

    const auto& worldtube_inbox =
        ActionTesting::get_inbox_tag<worldtube_chare,
                                     Tags::SphericalHarmonicsInbox<Dim>>(runner,
                                                                         0);
    // the single spherical-shell element covers the entire worldtube
    // boundary, so exactly one message must have been sent
    REQUIRE(worldtube_inbox.size() == 1);
    REQUIRE(worldtube_inbox.count(dummy_time_step_id) == 1);
    const auto& time_step_data = worldtube_inbox.at(dummy_time_step_id);
    REQUIRE(time_step_data.size() == 1);
    REQUIRE(time_step_data.count(element_id) == 1);
    const auto& sent_coefs = time_step_data.at(element_id);
    const auto& sent_psi_coefs =
        get(get<CurvedScalarWave::Tags::Psi>(sent_coefs));
    const auto& sent_dt_psi_coefs =
        get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(sent_coefs));
    const size_t num_modes = square(expansion_order + 1);
    REQUIRE(sent_psi_coefs.size() == num_modes);
    REQUIRE(sent_dt_psi_coefs.size() == num_modes);
    // the quadrature is exact for these band-limited integrands, so we can
    // check to tight tolerance
    Approx s2_approx = Approx::custom().epsilon(1e-12).scale(1.0);
    const double r_squared = square(worldtube_radius);
    for (size_t index = 0; index < num_modes; ++index) {
      CAPTURE(index);
      CHECK(sent_psi_coefs[index] ==
            s2_approx(gsl::at(psi_coefs, index) * r_squared));
      CHECK(sent_dt_psi_coefs[index] ==
            s2_approx(gsl::at(dt_psi_coefs, index) * r_squared));
    }
  }
}

// This test checks that `SendToWorldtube` integrates the regular field on the
// element surfaces abutting the worldtube and sends this data to the worldtube
// which reduces it in `ReceiveElementData`. The projection is done in the
// co-moving grid frame, so the time derivative is transformed with an
// advective term which is checked to vanish here for zero mesh velocity.
// There are several other actions in the action lists which are needed for
// the test to compile but are never used. The iterative scheme is tested in
// Test_Iterations.cpp.
SPECTRE_TEST_CASE("Unit.CurvedScalarWave.Worldtube.SendToWorldtube", "[Unit]") {
  static constexpr size_t Dim = 3;
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> dist(-10., 10.);
  using metavars = MockMetavariables<Dim>;
  domain::creators::register_derived_with_charm();
  using element_chare = MockElementArray<metavars>;
  using worldtube_chare = MockWorldtubeSingleton<metavars>;
  const size_t initial_extent = 10;
  const size_t face_size = initial_extent * initial_extent;
  const auto basis = Spectral::Basis::Legendre;
  const auto quadrature = Spectral::Quadrature::GaussLobatto;
  // unused but the tag is needed to compile
  const double time = std::numeric_limits<double>::signaling_NaN();
  // we create several differently refined shells so a different number of
  // elements sends data
  for (const auto& [expansion_order, initial_refinement, worldtube_radius] :
       cartesian_product(std::array<size_t, 3>{0, 1, 2},
                         std::array<size_t, 2>{0, 1},
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
    // self force and therefore iterative scheme is turned off
    const auto puncture_field_options =
        CurvedScalarWave::Worldtube::PunctureField{
            CurvedScalarWave::Worldtube::PunctureField::Schwarzschild{
                expansion_order, 1.}};
    tuples::TaggedTuple<domain::Tags::Domain<Dim>, Tags::ExcisionSphere<Dim>,
                        Tags::WorldtubeRadius, Tags::PunctureFieldConfig,
                        Tags::ExpansionOrder, Tags::MaxIterations, Tags::Charge,
                        Tags::Mass, ::Tags::Time, Tags::SelfForceTurnOnTime,
                        Tags::SelfForceTurnOnInterval>
        tuple_of_opts{shell.create_domain(),
                      excision_sphere,
                      excision_sphere.radius(),
                      puncture_field_options,
                      expansion_order,
                      static_cast<size_t>(0),
                      0.1,
                      std::nullopt,
                      time,
                      std::nullopt,
                      std::nullopt};
    ActionTesting::MockRuntimeSystem<metavars> runner{std::move(tuple_of_opts)};
    const auto element_ids = initial_element_ids(initial_refinements);
    const auto& blocks = shell_domain.blocks();

    using puncture_field_type =
        Variables<tmpl::list<CurvedScalarWave::Tags::Psi,
                             ::Tags::dt<CurvedScalarWave::Tags::Psi>,
                             ::Tags::deriv<CurvedScalarWave::Tags::Psi,
                                           tmpl::size_t<3>, Frame::Inertial>>>;

    // The puncture field will get subtracted from the DG field. Here, we set
    // the puncture field to 0, so psi and dt_psi are integrated directly
    // and we can check the analytical result.
    const puncture_field_type puncture_field{face_size, 0.};
    const double psi_coefs_0 = dist(generator);
    const double pi_coefs_0 = dist(generator);
    const auto psi_coefs_1 =
        make_with_random_values<tnsr::i<double, Dim, Frame::Grid>>(
            make_not_null(&generator), dist, 0.);
    const auto pi_coefs_1 =
        make_with_random_values<tnsr::i<double, Dim, Frame::Grid>>(
            make_not_null(&generator), dist, 0.);
    const auto psi_coefs_2 =
        make_with_random_values<tnsr::ii<double, Dim, Frame::Grid>>(
            make_not_null(&generator), dist, 0.);
    const auto pi_coefs_2 =
        make_with_random_values<tnsr::ii<double, Dim, Frame::Grid>>(
            make_not_null(&generator), dist, 0.);
    double psi_coefs_2_trace = 0.;
    double pi_coefs_2_trace = 0.;
    for (size_t i = 0; i < Dim; ++i) {
      psi_coefs_2_trace += psi_coefs_2.get(i, i);
      pi_coefs_2_trace += pi_coefs_2.get(i, i);
    }
    const Time dummy_time{{1., 2.}, {1, 2}};
    const TimeStepId dummy_time_step_id{true, 123, dummy_time};
    // the particle is at the center of the excision sphere which is the
    // origin in the Shell domain
    tnsr::I<double, Dim> particle_position(0.);
    auto particle_velocity = particle_position;
    const std::array<tnsr::I<double, Dim>, 2> particle_pos_vel{
        {std::move(particle_position), std::move(particle_velocity)}};
    for (const auto& element_id : element_ids) {
      auto element = domain::create_initial_element(element_id, blocks,
                                                    initial_refinements);
      auto mesh = domain::create_initial_mesh(initial_extents, element, basis,
                                              quadrature);
      const auto& my_block = blocks.at(element_id.block_id());
      const ElementMap inertial_element_map(
          element_id, my_block.stationary_map().get_clone());
      const ElementMap<Dim, Frame::Grid> grid_element_map(
          element_id, my_block.stationary_map().get_to_grid_frame());
      const auto logical_coords = logical_coordinates(mesh);
      const auto grid_coords = grid_element_map(logical_coords);

      auto inertial_inv_jacobian =
          inertial_element_map.inv_jacobian(logical_coords);
      const size_t grid_size = mesh.number_of_grid_points();
      // the domain is stationary, so the grid to inertial map is the
      // identity. This tag is only needed to compile `ReceiveWorldtubeData`
      // which is never called in this test.
      InverseJacobian<DataVector, Dim, Frame::Grid, Frame::Inertial>
          grid_to_inertial_inv_jacobian(grid_size, 0.);
      for (size_t i = 0; i < Dim; ++i) {
        grid_to_inertial_inv_jacobian.get(i, i) = 1.;
      }
      // we set lapse and shift to Minkowski so dt Psi = - Pi, and the value we
      // pass in for Pi will get integrated directly
      Scalar<DataVector> lapse(grid_size, 1.);
      tnsr::I<DataVector, Dim, Frame::Inertial> shift(grid_size, 0.);
      typename CurvedScalarWave::System<Dim>::variables_tag::type evolved_vars(
          grid_size, 0.);
      const bool is_abutting =
          excision_sphere.abutting_direction(element_id).has_value();
      get(get<CurvedScalarWave::Tags::Psi>(evolved_vars)) = psi_coefs_0;
      get(get<CurvedScalarWave::Tags::Pi>(evolved_vars)) = pi_coefs_0;
      if (expansion_order > 0) {
        for (size_t i = 0; i < Dim; ++i) {
          get(get<CurvedScalarWave::Tags::Psi>(evolved_vars)) +=
              psi_coefs_1.get(i) * grid_coords.get(i);
          get(get<CurvedScalarWave::Tags::Pi>(evolved_vars)) +=
              pi_coefs_1.get(i) * grid_coords.get(i);
        }
      }
      if (expansion_order > 1) {
        for (size_t i = 0; i < Dim; ++i) {
          for (size_t j = 0; j < Dim; ++j) {
            get(get<CurvedScalarWave::Tags::Psi>(evolved_vars)) +=
                psi_coefs_2.get(i, j) * grid_coords.get(i) * grid_coords.get(j);
            get(get<CurvedScalarWave::Tags::Pi>(evolved_vars)) +=
                pi_coefs_2.get(i, j) * grid_coords.get(i) * grid_coords.get(j);
          }
        }
      }
      std::optional<puncture_field_type> optional_puncture_field =
          is_abutting ? std::make_optional<puncture_field_type>(puncture_field)
                      : std::nullopt;
      // we set the mesh velocity to zero so the advective term vanishes and
      // we can recover the exact time derivative
      std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> mesh_velocity =
          std::make_optional<tnsr::I<DataVector, Dim, Frame::Inertial>>(
              grid_size, 0.);

      ActionTesting::emplace_array_component_and_initialize<element_chare>(
          &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
          element_id,
          {std::move(element), std::move(mesh), grid_coords,
           std::move(optional_puncture_field), std::move(shift),
           std::move(lapse), std::move(inertial_inv_jacobian),
           std::move(grid_to_inertial_inv_jacobian), evolved_vars,
           std::move(mesh_velocity), dummy_time_step_id, particle_pos_vel,
           static_cast<size_t>(0), std::nullopt});
    }

    std::unordered_map<ElementId<Dim>, tnsr::I<DataVector, Dim, Frame::Grid>>
        element_faces_grid_coords{};
    Initialization::InitializeElementFacesGridCoordinates<Dim>::apply(
        make_not_null(&element_faces_grid_coords), initial_extents,
        initial_refinements, quadrature, shell_domain, excision_sphere);

    // these are all unused
    tuples::TaggedTuple<
        gr::Tags::SpacetimeMetric<double, Dim>,
        gr::Tags::InverseSpacetimeMetric<double, Dim>,
        gr::Tags::SpacetimeChristoffelSecondKind<double, Dim>,
        gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim>,
        Tags::TimeDilationFactor>
        background_quantities{};
    typename MockWorldtubeSingleton<MockMetavariables<Dim>>::variables_tag::type
        evolved_vars_wt{};
    typename MockWorldtubeSingleton<
        MockMetavariables<Dim>>::dt_variables_tag::type dt_variables{};
    tnsr::I<double, Dim> geodesic_acc{};
    tnsr::AA<double, Dim, Frame::Grid> inverse_spacetime_metric{};
    tnsr::A<double, Dim, Frame::Grid> trace_spacetime_christoffel{};
    ActionTesting::emplace_singleton_component_and_initialize<worldtube_chare>(
        &runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
        {element_faces_grid_coords, dummy_time_step_id, static_cast<size_t>(0),
         geodesic_acc, particle_pos_vel, background_quantities, evolved_vars_wt,
         dt_variables, inverse_spacetime_metric, trace_spacetime_christoffel});

    ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

    // ReceiveElementData should not be ready yet as the worldtube has not
    // received any data
    CHECK(not ActionTesting::next_action_if_ready<worldtube_chare>(
        make_not_null(&runner), 0));
    // SendToWorldtube called on all elements
    for (const auto& element_id : element_ids) {
      ActionTesting::next_action<element_chare>(make_not_null(&runner),
                                                element_id);
      if (excision_sphere.abutting_direction(element_id).has_value()) {
        // the mesh velocity is zero, so the advective term should vanish
        const auto& regular_field_advective_term =
            ActionTesting::get_databox_tag<
                element_chare, Tags::RegularFieldAdvectiveTerm<Dim>>(
                runner, element_id);
        CHECK_ITERABLE_APPROX(regular_field_advective_term,
                              Scalar<DataVector>(face_size, 0.));
      }
    }

    using inbox_tag = Tags::SphericalHarmonicsInbox<Dim>;
    const auto& worldtube_inbox =
        ActionTesting::get_inbox_tag<worldtube_chare, inbox_tag>(runner, 0);
    CHECK(worldtube_inbox.count(dummy_time_step_id));
    auto time_step_data = worldtube_inbox.at(dummy_time_step_id);
    // these are all the element ids of elements abutting the worldtube, we
    // check that these are the ones that were sent.
    for (const auto& [element_id, _] : element_faces_grid_coords) {
      CHECK(time_step_data.count(element_id));
      time_step_data.erase(element_id);
    }
    // Check that have received only data from elements abutting the worldtube
    CHECK(time_step_data.empty());
    // ReceiveElementData called
    CHECK(ActionTesting::next_action_if_ready<worldtube_chare>(
        make_not_null(&runner), 0));
    CHECK(worldtube_inbox.empty());
    const auto& psi_monopole_worldtube = ActionTesting::get_databox_tag<
        worldtube_chare,
        Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>>(runner,
                                                                       0);
    const auto& dt_psi_monopole_worldtube = ActionTesting::get_databox_tag<
        worldtube_chare, Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0,
                                              Dim, Frame::Grid>>(runner, 0);
    // the integral is over a low resolution DG grid which introduces a large
    // error
    Approx apprx = Approx::custom().epsilon(1e-8).scale(1.0);
    if (expansion_order < 2) {
      CHECK(get(psi_monopole_worldtube) == apprx(psi_coefs_0));
      CHECK(get(dt_psi_monopole_worldtube) == -apprx(pi_coefs_0));
    }
    if (expansion_order > 0) {
      const auto& psi_dipole_worldtube = ActionTesting::get_databox_tag<
          worldtube_chare,
          Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>>(runner,
                                                                         0);
      const auto& dt_psi_dipole_worldtube = ActionTesting::get_databox_tag<
          worldtube_chare, Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>,
                                                1, Dim, Frame::Grid>>(runner,
                                                                      0);
      CHECK_ITERABLE_CUSTOM_APPROX(psi_dipole_worldtube, psi_coefs_1, apprx);
      for (size_t i = 0; i < Dim; ++i) {
        CHECK(dt_psi_dipole_worldtube.get(i) == apprx(-pi_coefs_1.get(i)));
      }
      if (expansion_order > 1) {
        // the trace of the second order coefficients gets absorbed by the
        // monopole
        const double expected_psi_monopole =
            psi_coefs_0 + square(worldtube_radius) * psi_coefs_2_trace / 3.;
        const double expected_pi_monopole =
            pi_coefs_0 + square(worldtube_radius) * pi_coefs_2_trace / 3.;
        CHECK(get(psi_monopole_worldtube) == apprx(expected_psi_monopole));
        CHECK(get(dt_psi_monopole_worldtube) == -apprx(expected_pi_monopole));
        const auto& psi_quadrupole_worldtube = ActionTesting::get_databox_tag<
            worldtube_chare,
            Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Grid>>(
            runner, 0);
        const auto& dt_psi_quadrupole_worldtube =
            ActionTesting::get_databox_tag<
                worldtube_chare,
                Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 2, Dim,
                                     Frame::Grid>>(runner, 0);
        tnsr::ii<double, Dim, Frame::Grid> expected_psi_quadrupole{};
        tnsr::ii<double, Dim, Frame::Grid> expected_dt_psi_quadrupole{};
        for (size_t i = 0; i < Dim; ++i) {
          for (size_t j = 0; j < Dim; ++j) {
            expected_psi_quadrupole.get(i, j) = psi_coefs_2.get(i, j);
            expected_dt_psi_quadrupole.get(i, j) = -pi_coefs_2.get(i, j);
          }
          expected_psi_quadrupole.get(i, i) -= psi_coefs_2_trace / 3.;
          expected_dt_psi_quadrupole.get(i, i) += pi_coefs_2_trace / 3.;
        }
        CHECK_ITERABLE_CUSTOM_APPROX(psi_quadrupole_worldtube,
                                     expected_psi_quadrupole, apprx);
        CHECK_ITERABLE_CUSTOM_APPROX(dt_psi_quadrupole_worldtube,
                                     expected_dt_psi_quadrupole, apprx);
      }
    }
  }
  test_s2_projection();
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
