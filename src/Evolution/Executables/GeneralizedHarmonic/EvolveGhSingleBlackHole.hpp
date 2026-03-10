// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

#include "ControlSystem/Actions/InitializeMeasurements.hpp"
#include "ControlSystem/CleanFunctionsOfTime.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Size/Factory.hpp"
#include "ControlSystem/ControlErrors/Size/State.hpp"
#include "ControlSystem/Measurements/NonFactoryCreatable.hpp"
#include "ControlSystem/Measurements/SingleHorizon.hpp"
#include "ControlSystem/Metafunctions.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "ControlSystem/Systems/Shape.hpp"
#include "ControlSystem/Systems/Size.hpp"
#include "ControlSystem/Systems/Translation.hpp"
#include "ControlSystem/Trigger.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Evolution/Actions/RunEventsAndTriggers.hpp"
#include "Evolution/Executables/GeneralizedHarmonic/Deadlock.hpp"
#include "Evolution/Executables/GeneralizedHarmonic/GeneralizedHarmonicBase.hpp"
#include "Evolution/Systems/Cce/Callbacks/DumpBondiSachsOnWorldtube.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/SetInitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BjorhusImpl.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Options/FactoryHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/ArrayCollection/DgElementCollection.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/Protocols/RegistrationMetavariables.hpp"
#include "ParallelAlgorithms/Actions/FunctionsOfTimeAreReady.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/Amr/Events/ObserveAmrCriteria.hpp"
#include "ParallelAlgorithms/Amr/Events/ObserveAmrStats.hpp"
#include "ParallelAlgorithms/Amr/Events/RefineMesh.hpp"
#include "ParallelAlgorithms/Amr/Projectors/CopyFromCreatorOrLeaveAsIs.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FailedHorizonFind.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/IgnoreFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ObserveFieldsOnHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Callbacks/ObserveTimeSeriesOnHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Component.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeExcisionBoundaryVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeHorizonVolumeQuantities.tpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Criteria/Criterion.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Criteria/Factory.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Events/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/InterpolationTarget.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Protocols/HorizonMetavars.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsOnFailure.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/CleanUpInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementInitInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetReceiveVars.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceivePoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveSurfaceData.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Events/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/ComputeVarsToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/CubicCurvatureScalars.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/CovariantDerivOfExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/QuadraticCurvatureScalars.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/AdvanceTime.hpp"
#include "Time/ChangeSlabSize/Action.hpp"
#include "Time/ChangeSlabSize/Tags.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/Tags/StepperErrors.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeAndPrevious.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/Numeric.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube_diagnostics {
using diagnostic_frame = Frame::Inertial;
using volume_source_tags =
    tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3, diagnostic_frame>,
               gh::Tags::Pi<DataVector, 3, diagnostic_frame>,
               gh::Tags::Phi<DataVector, 3, diagnostic_frame>,
               ::Tags::deriv<gh::Tags::Pi<DataVector, 3, diagnostic_frame>,
                             tmpl::size_t<3>, diagnostic_frame>,
               ::Tags::deriv<gh::Tags::Phi<DataVector, 3, diagnostic_frame>,
                             tmpl::size_t<3>, diagnostic_frame>>;

using volume_target_tags =
    tmpl::list<gr::Tags::SpatialMetric<DataVector, 3, diagnostic_frame>,
               gr::Tags::InverseSpatialMetric<DataVector, 3, diagnostic_frame>,
               gr::Tags::WeylElectric<DataVector, 3, diagnostic_frame>,
               gr::Tags::WeylMagnetic<DataVector, 3, diagnostic_frame>>;
using worldtube_volume_target_tags =
    tmpl::push_back<volume_target_tags,
                    gr::Tags::SpacetimeMetric<DataVector, 3, diagnostic_frame>>;

struct ComputeVolumeQuantities
    : tt::ConformsTo<intrp::protocols::ComputeVarsToInterpolate> {
  using allowed_src_tags = volume_source_tags;
  using required_src_tags = volume_source_tags;

  template <typename TargetFrame>
  using allowed_dest_tags =
      tmpl::conditional_t<std::is_same_v<TargetFrame, diagnostic_frame>,
                          worldtube_volume_target_tags, tmpl::list<>>;
  template <typename TargetFrame>
  using required_dest_tags =
      tmpl::conditional_t<std::is_same_v<TargetFrame, diagnostic_frame>,
                          volume_target_tags, tmpl::list<>>;

  template <typename SrcTagList, typename DestTagList>
  static void apply(const gsl::not_null<Variables<DestTagList>*> target_vars,
                    const Variables<SrcTagList>& src_vars, const Mesh<3>&) {
    static_assert(
        std::is_same_v<tmpl::list_difference<SrcTagList, allowed_src_tags>,
                       tmpl::list<>>,
        "Found a source tag that is not allowed");
    static_assert(
        std::is_same_v<tmpl::list_difference<required_src_tags, SrcTagList>,
                       tmpl::list<>>,
        "A required source tag is missing");
    static_assert(
        std::is_same_v<tmpl::list_difference<
                           DestTagList, allowed_dest_tags<diagnostic_frame>>,
                       tmpl::list<>>,
        "Found a destination tag that is not allowed");

    if (target_vars->number_of_grid_points() !=
        src_vars.number_of_grid_points()) {
      target_vars->initialize(src_vars.number_of_grid_points());
    }

    auto& spatial_metric =
        get<gr::Tags::SpatialMetric<DataVector, 3, diagnostic_frame>>(
            *target_vars);
    auto& inverse_spatial_metric =
        get<gr::Tags::InverseSpatialMetric<DataVector, 3, diagnostic_frame>>(
            *target_vars);
    auto& weyl_electric =
        get<gr::Tags::WeylElectric<DataVector, 3, diagnostic_frame>>(
            *target_vars);
    auto& weyl_magnetic =
        get<gr::Tags::WeylMagnetic<DataVector, 3, diagnostic_frame>>(
            *target_vars);

    const auto& spacetime_metric =
        get<gr::Tags::SpacetimeMetric<DataVector, 3, diagnostic_frame>>(
            src_vars);
    const auto& pi =
        get<gh::Tags::Pi<DataVector, 3, diagnostic_frame>>(src_vars);
    const auto& phi =
        get<gh::Tags::Phi<DataVector, 3, diagnostic_frame>>(src_vars);
    const auto& d_pi =
        get<::Tags::deriv<gh::Tags::Pi<DataVector, 3, diagnostic_frame>,
                          tmpl::size_t<3>, diagnostic_frame>>(src_vars);
    const auto& d_phi =
        get<::Tags::deriv<gh::Tags::Phi<DataVector, 3, diagnostic_frame>,
                          tmpl::size_t<3>, diagnostic_frame>>(src_vars);

    if constexpr (tmpl::list_contains_v<DestTagList,
                                        gr::Tags::SpacetimeMetric<
                                            DataVector, 3, diagnostic_frame>>) {
      get<gr::Tags::SpacetimeMetric<DataVector, 3, diagnostic_frame>>(
          *target_vars) = spacetime_metric;
    }

    Scalar<DataVector> det_spatial_metric{src_vars.number_of_grid_points()};
    Scalar<DataVector> lapse{src_vars.number_of_grid_points()};
    Scalar<DataVector> sqrt_det_spatial_metric{
        src_vars.number_of_grid_points()};
    tnsr::I<DataVector, 3, diagnostic_frame> shift{
        src_vars.number_of_grid_points()};
    tnsr::AA<DataVector, 3, diagnostic_frame> inverse_spacetime_metric{
        src_vars.number_of_grid_points()};
    tnsr::A<DataVector, 3, diagnostic_frame> spacetime_normal_vector{
        src_vars.number_of_grid_points()};
    tnsr::ii<DataVector, 3, diagnostic_frame> extrinsic_curvature{
        src_vars.number_of_grid_points()};
    tnsr::Ijj<DataVector, 3, diagnostic_frame> spatial_christoffel{
        src_vars.number_of_grid_points()};
    tnsr::ijj<DataVector, 3, diagnostic_frame>
        covariant_deriv_extrinsic_curvature{src_vars.number_of_grid_points()};
    tnsr::ii<DataVector, 3, diagnostic_frame> spatial_ricci{
        src_vars.number_of_grid_points()};

    gr::spatial_metric(make_not_null(&spatial_metric), spacetime_metric);
    determinant_and_inverse(make_not_null(&det_spatial_metric),
                            make_not_null(&inverse_spatial_metric),
                            spatial_metric);
    gr::shift(make_not_null(&shift), spacetime_metric, inverse_spatial_metric);
    gr::lapse(make_not_null(&lapse), shift, spacetime_metric);
    gr::inverse_spacetime_metric(make_not_null(&inverse_spacetime_metric),
                                 lapse, shift, inverse_spatial_metric);
    gr::spacetime_normal_vector(make_not_null(&spacetime_normal_vector), lapse,
                                shift);
    gh::extrinsic_curvature(make_not_null(&extrinsic_curvature),
                            spacetime_normal_vector, pi, phi);
    gh::christoffel_second_kind(make_not_null(&spatial_christoffel), phi,
                                inverse_spatial_metric);
    gh::covariant_deriv_of_extrinsic_curvature(
        make_not_null(&covariant_deriv_extrinsic_curvature),
        extrinsic_curvature, spacetime_normal_vector, spatial_christoffel,
        inverse_spacetime_metric, phi, d_pi, d_phi);
    gh::spatial_ricci_tensor(make_not_null(&spatial_ricci), phi, d_phi,
                             inverse_spatial_metric);
    get(sqrt_det_spatial_metric) = sqrt(get(det_spatial_metric));
    gr::weyl_electric(make_not_null(&weyl_electric), spatial_ricci,
                      extrinsic_curvature, inverse_spatial_metric);
    gr::weyl_magnetic(make_not_null(&weyl_magnetic),
                      covariant_deriv_extrinsic_curvature, spatial_metric,
                      sqrt_det_spatial_metric);
  }
};

template <typename SurfaceFrame>
inline tnsr::I<DataVector, 3, diagnostic_frame> inward_unit_interface_normal(
    const tnsr::I<DataVector, 3, SurfaceFrame>& coords,
    const ylm::Strahlkorper<SurfaceFrame>& strahlkorper,
    const tnsr::II<DataVector, 3, diagnostic_frame>& inverse_spatial_metric) {
  const size_t num_points = get<0>(coords).size();
  tnsr::i<DataVector, 3, diagnostic_frame> inward_normal_one_form{num_points};
  Scalar<DataVector> euclidean_radius{num_points, 0.0};
  Scalar<DataVector> normal_magnitude{num_points, 0.0};
  tnsr::I<DataVector, 3, diagnostic_frame> unit_interface_normal_vector{
      num_points};

  const auto center = strahlkorper.expansion_center();
  for (size_t i = 0; i < 3; ++i) {
    inward_normal_one_form.get(i) = center[i] - coords.get(i);
    get(euclidean_radius) += square(inward_normal_one_form.get(i));
  }
  get(euclidean_radius) = sqrt(get(euclidean_radius));
  for (size_t i = 0; i < 3; ++i) {
    inward_normal_one_form.get(i) /= get(euclidean_radius);
  }

  for (size_t i = 0; i < 3; ++i) {
    unit_interface_normal_vector.get(i) = 0.0;
    for (size_t j = 0; j < 3; ++j) {
      unit_interface_normal_vector.get(i) +=
          inverse_spatial_metric.get(i, j) * inward_normal_one_form.get(j);
      get(normal_magnitude) += inverse_spatial_metric.get(i, j) *
                               inward_normal_one_form.get(i) *
                               inward_normal_one_form.get(j);
    }
  }
  get(normal_magnitude) = sqrt(get(normal_magnitude));
  for (size_t i = 0; i < 3; ++i) {
    unit_interface_normal_vector.get(i) /= get(normal_magnitude);
  }
  return unit_interface_normal_vector;
}

struct WorldtubeWeylScalarDiagnostics : db::SimpleTag {
  using type =
      gh::BoundaryConditions::Bjorhus::detail::WorldtubeWeylScalarDiagnostics<
          DataVector, diagnostic_frame>;
  static std::string name() { return "WorldtubeWeylScalarDiagnostics"; }
};

struct WorldtubeWeylScalarDiagnosticsCompute : WorldtubeWeylScalarDiagnostics,
                                               db::ComputeTag {
  using base = WorldtubeWeylScalarDiagnostics;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      gr::Tags::SpatialMetric<DataVector, 3, diagnostic_frame>,
      gr::Tags::InverseSpatialMetric<DataVector, 3, diagnostic_frame>,
      gr::Tags::WeylElectric<DataVector, 3, diagnostic_frame>,
      gr::Tags::WeylMagnetic<DataVector, 3, diagnostic_frame>,
      ylm::Tags::CartesianCoords<Frame::Grid>,
      ylm::Tags::Strahlkorper<Frame::Grid>>;

  static void function(
      const gsl::not_null<return_type*> diagnostics,
      const tnsr::ii<DataVector, 3, diagnostic_frame>& spatial_metric,
      const tnsr::II<DataVector, 3, diagnostic_frame>& inverse_spatial_metric,
      const tnsr::ii<DataVector, 3, diagnostic_frame>& weyl_electric,
      const tnsr::ii<DataVector, 3, diagnostic_frame>& weyl_magnetic,
      const tnsr::I<DataVector, 3, Frame::Grid>& coords,
      const ylm::Strahlkorper<Frame::Grid>& strahlkorper) {
    const auto unit_interface_normal_vector = inward_unit_interface_normal(
        coords, strahlkorper, inverse_spatial_metric);
    gh::BoundaryConditions::Bjorhus::detail::worldtube_weyl_scalar_diagnostics(
        diagnostics, weyl_electric, weyl_magnetic, spatial_metric,
        inverse_spatial_metric, unit_interface_normal_vector);
  }
};

struct WorldtubeGaussBonnetScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "WorldtubeGaussBonnetScalar"; }
};

struct WorldtubeGaussBonnetScalarCompute : WorldtubeGaussBonnetScalar,
                                           db::ComputeTag {
  using base = WorldtubeGaussBonnetScalar;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<WorldtubeWeylScalarDiagnostics>;
  static void function(
      const gsl::not_null<return_type*> result,
      const WorldtubeWeylScalarDiagnostics::type& diagnostics) {
    *result = diagnostics.gauss_bonnet_scalar;
  }
};

template <size_t PsiIndex, bool ImaginaryPart>
struct WorldtubePsiComponent : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() {
    return MakeString{} << "WorldtubePsi" << PsiIndex
                        << (ImaginaryPart ? "Imag" : "Real");
  }
};

template <size_t PsiIndex, bool ImaginaryPart>
struct WorldtubePsiComponentCompute
    : WorldtubePsiComponent<PsiIndex, ImaginaryPart>,
      db::ComputeTag {
  using base = WorldtubePsiComponent<PsiIndex, ImaginaryPart>;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<WorldtubeWeylScalarDiagnostics>;
  static void function(
      const gsl::not_null<return_type*> result,
      const WorldtubeWeylScalarDiagnostics::type& diagnostics) {
    if constexpr (ImaginaryPart) {
      get(*result) = imag(get(diagnostics.weyl_scalars[PsiIndex]));
    } else {
      get(*result) = real(get(diagnostics.weyl_scalars[PsiIndex]));
    }
  }
};

template <bool ImaginaryPart>
struct WorldtubePsi0BcComponent : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() {
    return ImaginaryPart ? "WorldtubePsi0BcImag" : "WorldtubePsi0BcReal";
  }
};

template <bool ImaginaryPart>
struct WorldtubePsi0BcComponentCompute
    : WorldtubePsi0BcComponent<ImaginaryPart>,
      db::ComputeTag {
  using base = WorldtubePsi0BcComponent<ImaginaryPart>;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<WorldtubeWeylScalarDiagnostics>;
  static void function(
      const gsl::not_null<return_type*> result,
      const WorldtubeWeylScalarDiagnostics::type& diagnostics) {
    if constexpr (ImaginaryPart) {
      get(*result) = imag(get(diagnostics.inferred_psi0));
    } else {
      get(*result) = real(get(diagnostics.inferred_psi0));
    }
  }
};

struct WorldtubePsi2Kinnersley : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "WorldtubePsi2Kinnersley"; }
};

struct WorldtubePsi2KinnersleyCompute : WorldtubePsi2Kinnersley,
                                        db::ComputeTag {
  using base = WorldtubePsi2Kinnersley;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<WorldtubeWeylScalarDiagnostics>;
  static void function(
      const gsl::not_null<return_type*> result,
      const WorldtubeWeylScalarDiagnostics::type& diagnostics) {
    *result = diagnostics.psi2_kinnersley;
  }
};

using surface_observe_tags = tmpl::list<
    WorldtubeGaussBonnetScalar, gr::Tags::PontryaginScalar<DataVector>,
    gr::Tags::CubicInvariantReal<DataVector>,
    gr::Tags::CubicInvariantImag<DataVector>, WorldtubePsiComponent<0, false>,
    WorldtubePsiComponent<0, true>, WorldtubePsiComponent<1, false>,
    WorldtubePsiComponent<1, true>, WorldtubePsiComponent<2, false>,
    WorldtubePsiComponent<2, true>, WorldtubePsiComponent<3, false>,
    WorldtubePsiComponent<3, true>, WorldtubePsiComponent<4, false>,
    WorldtubePsiComponent<4, true>, WorldtubePsi0BcComponent<false>,
    WorldtubePsi0BcComponent<true>, WorldtubePsi2Kinnersley>;

using surface_compute_items = tmpl::list<
    WorldtubeWeylScalarDiagnosticsCompute, WorldtubeGaussBonnetScalarCompute,
    gr::Tags::PontryaginScalarCompute<DataVector, 3, diagnostic_frame>,
    gr::Tags::CubicInvariantRealCompute<DataVector, 3, diagnostic_frame>,
    gr::Tags::CubicInvariantImagCompute<DataVector, 3, diagnostic_frame>,
    WorldtubePsiComponentCompute<0, false>,
    WorldtubePsiComponentCompute<0, true>,
    WorldtubePsiComponentCompute<1, false>,
    WorldtubePsiComponentCompute<1, true>,
    WorldtubePsiComponentCompute<2, false>,
    WorldtubePsiComponentCompute<2, true>,
    WorldtubePsiComponentCompute<3, false>,
    WorldtubePsiComponentCompute<3, true>,
    WorldtubePsiComponentCompute<4, false>,
    WorldtubePsiComponentCompute<4, true>,
    WorldtubePsi0BcComponentCompute<false>,
    WorldtubePsi0BcComponentCompute<true>, WorldtubePsi2KinnersleyCompute>;
}  // namespace gh::worldtube_diagnostics

namespace gh::gb_translation {
struct Measurement : tt::ConformsTo<control_system::protocols::Measurement> {
  struct Submeasurement
      : tt::ConformsTo<control_system::protocols::Submeasurement> {
    static std::string name() { return Measurement::name(); }

   private:
    template <typename ControlSystems>
    struct InterpolationTarget
        : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
      static std::string name() { return "ControlSystemGaussBonnetDipole"; }
      using temporal_id = ::Tags::TimeAndPrevious<0>;
      using vars_to_interpolate_to_target =
          gh::worldtube_diagnostics::volume_target_tags;
      using compute_vars_to_interpolate =
          gh::worldtube_diagnostics::ComputeVolumeQuantities;
      using compute_items_on_source =
          tmpl::list<::Tags::TimeAndPreviousCompute<0>>;
      using compute_items_on_target = tmpl::list<
          gr::Tags::WeylElectricScalarCompute<DataVector, 3, Frame::Inertial>,
          gr::Tags::WeylMagneticScalarCompute<DataVector, 3, Frame::Inertial>,
          gr::Tags::GaussBonnetScalarCompute<DataVector>>;
      using compute_target_points =
          intrp::TargetPoints::Sphere<InterpolationTarget, ::Frame::Grid>;
      using post_interpolation_callbacks = tmpl::list<
          control_system::RunCallbacks<Submeasurement, ControlSystems>>;

      template <typename Metavariables>
      using interpolating_component =
          typename Metavariables::gh_dg_element_array;
    };

   public:
    template <typename ControlSystems>
    using interpolation_target_tag = InterpolationTarget<ControlSystems>;
    template <typename ControlSystems>
    using horizon_metavars = void;

    template <typename ControlSystems>
    using event = NonFactoryCreatableWrapper<
        intrp::Events::InterpolateWithoutInterpComponent<
            3, InterpolationTarget<ControlSystems>,
            gh::worldtube_diagnostics::volume_source_tags>>;
  };

  static std::string name() { return "GaussBonnetDipole"; }
  using submeasurements = tmpl::list<Submeasurement>;
};

struct ControlError : tt::ConformsTo<control_system::protocols::ControlError> {
  using object_centers = domain::object_list<>;
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Computes the translation control error from the Gauss-Bonnet dipole. "
      "This should not take any options."};

  // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
  std::optional<double> get_suggested_timescale() const { return std::nullopt; }

  void reset() {}
  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const ::TimescaleTuner<true>& /*unused*/,
                        const Parallel::GlobalCache<Metavariables>& /*cache*/,
                        const double /*time*/,
                        const std::string& /*function_of_time_name*/,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    return get<control_system::QueueTags::Center<::domain::ObjectLabel::None>>(
        measurements);
  }
};

template <size_t DerivOrder>
struct Translation : tt::ConformsTo<control_system::protocols::ControlSystem> {
  static constexpr size_t deriv_order = DerivOrder;

  static std::string name() { return "Translation"; }

  static std::optional<std::string> component_name(
      const size_t component, const size_t num_components) {
    ASSERT(num_components == 3,
           "Translation control expects 3 components but there are "
               << num_components << " instead.");
    return component == 0 ? "x" : component == 1 ? "y" : "z";
  }

  using measurement = Measurement;
  using control_error = ControlError;

  struct MeasurementQueue : db::SimpleTag {
    using type =
        LinkedMessageQueue<double, tmpl::list<control_system::QueueTags::Center<
                                       ::domain::ObjectLabel::None>>>;
  };

  using simple_tags = tmpl::list<MeasurementQueue>;

  struct process_measurement {
    template <typename Submeasurement>
    using argument_tags = tmpl::list<ylm::Tags::Strahlkorper<Frame::Grid>,
                                     gr::Tags::GaussBonnetScalar<DataVector>>;

    template <typename Metavariables>
    static void apply(Measurement::Submeasurement /*submeasurement*/,
                      const ylm::Strahlkorper<Frame::Grid>& strahlkorper,
                      const Scalar<DataVector>& gauss_bonnet_scalar,
                      Parallel::GlobalCache<Metavariables>& cache,
                      const LinkedMessageId<double>& measurement_id) {
      auto& control_sys_proxy = Parallel::get_parallel_component<
          ControlComponent<Metavariables, Translation>>(cache);
      Parallel::simple_action<::Actions::UpdateMessageQueue<
          MeasurementQueue, control_system::UpdateControlSystem<Translation>,
          control_system::QueueTags::Center<::domain::ObjectLabel::None>>>(
          control_sys_proxy, measurement_id,
          gauss_bonnet_dipole(strahlkorper, gauss_bonnet_scalar));
    }

   private:
    static DataVector gauss_bonnet_dipole(
        const ylm::Strahlkorper<Frame::Grid>& strahlkorper,
        const Scalar<DataVector>& gauss_bonnet_scalar) {
      ASSERT(strahlkorper.l_max() > 0 and strahlkorper.m_max() > 0,
             "Need l_max >= 1 and m_max >= 1 to compute a dipole.");
      const auto& ylm = strahlkorper.ylm_spherepack();
      const size_t points_per_sphere = ylm.physical_size();
      const DataVector& all_values = get(gauss_bonnet_scalar);
      ASSERT(all_values.size() >= points_per_sphere and
                 all_values.size() % points_per_sphere == 0,
             "Unexpected number of Gauss-Bonnet points: "
                 << all_values.size() << ", expected a multiple of "
                 << points_per_sphere << ".");

      // If multiple radii are observed, use the outermost sphere.
      const size_t offset = all_values.size() - points_per_sphere;
      DataVector values_on_control_sphere{points_per_sphere};
      for (size_t i = 0; i < points_per_sphere; ++i) {
        values_on_control_sphere[i] = sqrt(all_values[offset + i] / 16. / 3.);
      }

      const DataVector gb_coefs = ylm.phys_to_spec(values_on_control_sphere);
      ylm::SpherepackIterator iterator(strahlkorper.l_max(),
                                       strahlkorper.m_max());
      ModalVector l0_coefs{1, 0.0};
      l0_coefs[0] = gb_coefs[iterator.set(0, 0)()] * sqrt(M_PI / 2.);
      ModalVector l1_coefs{3, 0.0};
      l1_coefs[0] = gb_coefs[iterator.set(1, 1)()] * sqrt(M_PI);
      l1_coefs[1] = -gb_coefs[iterator.set(1, -1)()] * sqrt(M_PI);
      l1_coefs[2] = gb_coefs[iterator.set(1, 0)()] * sqrt(M_PI / 2.);

      const double monopole = l0_coefs[0];
      const double dipole_magnitude =
          sqrt(square(l1_coefs[0]) + square(l1_coefs[1]) + square(l1_coefs[2]));

      if (dipole_magnitude == 0.0) {
        return DataVector{3, 0.0};
      }

      ASSERT(monopole != 0.0,
             "Gauss-Bonnet monopole coefficient is zero, so the control "
             "error normalization is singular.");

      const double radius = 1.9;
      const double delta = dipole_magnitude * radius / (sqrt(3.0) * monopole);
      DataVector result{3};
      result[0] = delta * l1_coefs[0] / dipole_magnitude;
      result[1] = delta * l1_coefs[1] / dipole_magnitude;
      result[2] = delta * l1_coefs[2] / dipole_magnitude;
      Parallel::printf(
          "Gauss-Bonnet monopole: %e, dipole magnitude: %e, control error: "
          "(%e, %e, %e)\n",
          monopole, dipole_magnitude, result[0], result[1], result[2]);
      return result;
    }
  };
};
}  // namespace gh::gb_translation



template <bool UseLts>
struct EvolutionMetavars : public GeneralizedHarmonicTemplateBase<3, UseLts> {
  static constexpr bool local_time_stepping = UseLts;
  static constexpr size_t volume_dim = 3;
  using gh_base = GeneralizedHarmonicTemplateBase<volume_dim, UseLts>;
  using typename gh_base::initialize_initial_data_dependent_quantities_actions;
  using typename gh_base::system;

  static constexpr Options::String help{
      "Evolve the Einstein field equations using the Generalized Harmonic "
      "formulation,\n"
      "on a domain with a single horizon and corresponding excised region"};

  struct ApparentHorizon : tt::ConformsTo<ah::protocols::HorizonMetavars> {
    using time_tag = ah::Tags::ObservationTime<0>;

    using frame = ::Frame::Inertial;

    using horizon_find_callbacks = tmpl::list<
        ah::callbacks::ObserveTimeSeriesOnHorizon<
            ::ah::tags_for_observing<Frame::Inertial>, ApparentHorizon>,
        ah::callbacks::ObserveFieldsOnHorizon<::ah::surface_tags_for_observing,
                                              ApparentHorizon>>;
    using horizon_find_failure_callbacks =
        tmpl::list<ah::callbacks::FailedHorizonFind<ApparentHorizon, false>>;

    using compute_tags_on_element =
        tmpl::list<ah::Tags::ObservationTimeCompute<0>>;

    static constexpr ah::Destination destination = ah::Destination::Observation;

    static std::string name() { return "ApparentHorizon"; }
  };

  struct ExcisionBoundary
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe =
        tmpl::list<gr::Tags::Lapse<DataVector>,
                   gr::Tags::Shift<DataVector, 3, Frame::Grid>>;
    using compute_vars_to_interpolate =
        ah::ComputeExcisionBoundaryVolumeQuantities;
    using vars_to_interpolate_to_target = tags_to_observe;
    using compute_items_on_source = tmpl::list<>;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        intrp::TargetPoints::Sphere<ExcisionBoundary, ::Frame::Grid>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::ObserveSurfaceData<
            tags_to_observe, ExcisionBoundary, ::Frame::Grid>>;
    // run_callbacks
    template <typename metavariables>
    using interpolating_component = typename metavariables::gh_dg_element_array;
  };

  using control_systems =
      tmpl::list<control_system::Systems::Shape<
                     ::domain::ObjectLabel::None, 2,
                     control_system::measurements::SingleHorizon<
                         ::domain::ObjectLabel::None>>,
                 gh::gb_translation::Translation<2>,
                 control_system::Systems::Size<::domain::ObjectLabel::None, 2>>;

  static constexpr bool use_control_systems =
      tmpl::size<control_systems>::value > 0;

  struct BondiSachs;
  struct GaussBonnetSpheres;
  struct WorldtubeDiagnostics;

  using interpolation_target_tags = tmpl::push_back<
      control_system::metafunctions::interpolation_target_tags<control_systems>,
      ExcisionBoundary, BondiSachs, GaussBonnetSpheres, WorldtubeDiagnostics>;
  using source_vars_no_deriv =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, volume_dim>,
                 gh::Tags::Pi<DataVector, volume_dim>,
                 gh::Tags::Phi<DataVector, volume_dim>>;
  using curvature_surface_source_vars =
      gh::worldtube_diagnostics::volume_source_tags;
  using curvature_surface_observe_tags =
      tmpl::list<gr::Tags::GaussBonnetScalar<DataVector>,
                 gr::Tags::PontryaginScalar<DataVector>,
                 gr::Tags::CubicInvariantReal<DataVector>,
                 gr::Tags::CubicInvariantImag<DataVector>,
                 gr::Tags::SpacetimeMetric<DataVector, 3, Frame::Inertial>>;
  using curvature_surface_compute_items = tmpl::list<
      gr::Tags::WeylElectricScalarCompute<DataVector, 3, Frame::Inertial>,
      gr::Tags::WeylMagneticScalarCompute<DataVector, 3, Frame::Inertial>,
      gr::Tags::GaussBonnetScalarCompute<DataVector>,
      gr::Tags::PontryaginScalarCompute<DataVector, 3, Frame::Inertial>,
      gr::Tags::CubicInvariantRealCompute<DataVector, 3, Frame::Inertial>,
      gr::Tags::CubicInvariantImagCompute<DataVector, 3, Frame::Inertial>>;
  using worldtube_diagnostic_source_vars =
      gh::worldtube_diagnostics::volume_source_tags;

  struct BondiSachs : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    static std::string name() { return "BondiSachsInterpolation"; }
    using temporal_id = ::Tags::Time;
    using vars_to_interpolate_to_target = source_vars_no_deriv;
    using compute_target_points =
        intrp::TargetPoints::Sphere<BondiSachs, ::Frame::Inertial>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::DumpBondiSachsOnWorldtube<BondiSachs>>;
    using compute_items_on_target = tmpl::list<>;
    template <typename Metavariables>
    using interpolating_component = typename Metavariables::gh_dg_element_array;
  };

  struct GaussBonnetSpheres
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    static std::string name() { return "GaussBonnetSpheres"; }
    using temporal_id = ::Tags::Time;
    using compute_vars_to_interpolate =
        gh::worldtube_diagnostics::ComputeVolumeQuantities;
    using vars_to_interpolate_to_target =
        gh::worldtube_diagnostics::worldtube_volume_target_tags;
    using compute_items_on_target = curvature_surface_compute_items;
    using compute_target_points =
        intrp::TargetPoints::Sphere<GaussBonnetSpheres, ::Frame::Grid>;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::ObserveSurfaceData<
            curvature_surface_observe_tags, GaussBonnetSpheres, ::Frame::Grid>>;
    template <typename Metavariables>
    using interpolating_component = typename Metavariables::gh_dg_element_array;
  };

  struct WorldtubeDiagnostics
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    static std::string name() { return "WorldtubeDiagnostics"; }
    using temporal_id = ::Tags::Time;
    using compute_vars_to_interpolate =
        gh::worldtube_diagnostics::ComputeVolumeQuantities;
    using vars_to_interpolate_to_target =
        gh::worldtube_diagnostics::volume_target_tags;
    using compute_target_points =
        intrp::TargetPoints::Sphere<WorldtubeDiagnostics, ::Frame::Grid>;
    using compute_items_on_target =
        gh::worldtube_diagnostics::surface_compute_items;
    using post_interpolation_callbacks =
        tmpl::list<intrp::callbacks::ObserveSurfaceData<
            gh::worldtube_diagnostics::surface_observe_tags,
            WorldtubeDiagnostics, ::Frame::Grid>>;
    template <typename Metavariables>
    using interpolating_component = typename Metavariables::gh_dg_element_array;
  };

  // The interpolator_source_vars need to be the same in both the Interpolate
  // event and the InterpolateWithoutInterpComponent event.  The Interpolate
  // event interpolates to the horizon, and the
  // InterpolateWithoutInterpComponent event interpolates to the excision
  // boundary. Every Target gets the same interpolator_source_vars, so they need
  // to be made the same. Otherwise a static assert is triggered.
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = Options::add_factory_classes<
        // Restrict to monotonic time steppers in LTS to avoid control
        // systems deadlocking.
        tmpl::insert<
            tmpl::erase<typename gh_base::factory_creation::factory_classes,
                        LtsTimeStepper>,
            tmpl::pair<LtsTimeStepper,
                       TimeSteppers::monotonic_lts_time_steppers>>,
        tmpl::pair<ah::Criterion, ah::Criteria::standard_criteria>,
        tmpl::pair<
            Event,
            tmpl::flatten<tmpl::list<
                ah::Events::FindApparentHorizon<ApparentHorizon>,
                control_system::metafunctions::control_system_events<
                    control_systems>,
                control_system::CleanFunctionsOfTime,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, BondiSachs, source_vars_no_deriv>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, GaussBonnetSpheres, curvature_surface_source_vars>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, ExcisionBoundary, ::ah::source_vars<volume_dim>>,
                intrp::Events::InterpolateWithoutInterpComponent<
                    3, WorldtubeDiagnostics, worldtube_diagnostic_source_vars>,
                amr::Events::RefineMesh,
                amr::Events::ObserveAmrStats<volume_dim>>>>,
        tmpl::pair<DenseTrigger,
                   control_system::control_system_triggers<control_systems>>,
        tmpl::pair<control_system::size::State,
                   control_system::size::States::factory_creatable_states>>;
  };

  using typename gh_base::const_global_cache_tags;

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::append<
          tmpl::at<typename factory_creation::factory_classes, Event>,
          typename ExcisionBoundary::post_interpolation_callbacks,
          typename GaussBonnetSpheres::post_interpolation_callbacks,
          typename WorldtubeDiagnostics::post_interpolation_callbacks>>;

  using dg_registration_list = typename gh_base::dg_registration_list;

  using step_actions =
      typename gh_base::template step_actions<EvolutionMetavars,
                                              control_systems>;

  using initialization_actions = tmpl::push_back<
      tmpl::pop_back<typename gh_base::template initialization_actions<
          EvolutionMetavars, use_control_systems>>,
      control_system::Actions::InitializeMeasurements<control_systems>,
      intrp::Actions::ElementInitInterpPoints<volume_dim,
                                              interpolation_target_tags>,
      tmpl::back<typename gh_base::template initialization_actions<
          EvolutionMetavars, use_control_systems>>>;

  using gh_dg_element_array = DgElementArray<
      EvolutionMetavars,
      tmpl::flatten<tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::RegisterWithElementDataReader,
              tmpl::list<importers::Actions::RegisterWithElementDataReader,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::ImportInitialData,
              tmpl::list<gh::Actions::SetInitialData,
                         gh::Actions::ReceiveNumericInitialData,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions, system>>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<Parallel::Phase::Restart,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::WriteCheckpoint,
              tmpl::list<evolution::Actions::RunEventsAndTriggers<
                             Triggers::WhenToCheck::AtCheckpoints>,
                         Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<Parallel::Phase::CheckDomain,
                                 tmpl::list<::amr::Actions::SendAmrDiagnostics,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::flatten<tmpl::list<
                  ::domain::Actions::CheckFunctionsOfTimeAreReady<volume_dim>,
                  std::conditional_t<local_time_stepping,
                                     evolution::Actions::RunEventsAndTriggers<
                                         Triggers::WhenToCheck::AtSteps>,
                                     tmpl::list<>>,
                  evolution::Actions::RunEventsAndTriggers<
                      Triggers::WhenToCheck::AtSlabs>,
                  Actions::ChangeSlabSize, step_actions,
                  Actions::MutateApply<AdvanceTime<>>,
                  PhaseControl::Actions::ExecutePhaseChange>>>,
          Parallel::PhaseActions<
              Parallel::Phase::PostFailureCleanup,
              tmpl::list<Actions::RunEventsOnFailure<::Tags::Time>,
                         Parallel::Actions::TerminatePhase>>>>>;

  struct amr : tt::ConformsTo<::amr::protocols::AmrMetavariables> {
    using element_array = gh_dg_element_array;
    using projectors = tmpl::list<
        Initialization::ProjectTimeStepping<volume_dim>,
        evolution::dg::Initialization::ProjectDomain<volume_dim>,
        ::amr::projectors::ProjectVariables<volume_dim,
                                            typename system::variables_tag>,
        evolution::dg::Initialization::ProjectMortars<volume_dim,
                                                      local_time_stepping>,
        dg::Actions::InitializeFilters<
            typename gh_base::FilterEvolvedVariables>,
        Initialization::ProjectTimeStepperHistory<EvolutionMetavars>,
        evolution::Actions::ProjectRunEventsAndDenseTriggers,
        ::amr::projectors::DefaultInitialize<
            Initialization::Tags::InitialTimeDelta,
            Initialization::Tags::InitialSlabSize<gh_base::local_time_stepping>,
            ::domain::Tags::InitialExtents<volume_dim>,
            ::domain::Tags::InitialRefinementLevels<volume_dim>,
            evolution::dg::Tags::Quadrature,
            Tags::StepperErrors<typename system::variables_tag>,
            SelfStart::Tags::InitialValue<typename system::variables_tag>,
            SelfStart::Tags::InitialValue<Tags::TimeStep>>,
        ::amr::projectors::CopyFromCreatorOrLeaveAsIs<tmpl::push_back<
            tmpl::append<
                typename control_system::Actions::InitializeMeasurements<
                    control_systems>::simple_tags,
                tmpl::transform<
                    intrp::InterpolationTarget_detail::
                        get_non_sequential_target_tags<
                            interpolation_target_tags>,
                    tmpl::bind<intrp::Tags::PointInfo, tmpl::_1,
                               tmpl::pin<tmpl::size_t<volume_dim>>>>>,
            Tags::ChangeSlabSize::NumberOfExpectedMessages,
            Tags::ChangeSlabSize::NewSlabSize,
            ::Filters::Tags::Filter<Filters::Exponential<volume_dim, 0>>,
            ::Filters::Tags::Filter<ylm::TensorYlm::TensorYlmFilter>>>>;
    static constexpr bool keep_coarse_grids = false;
    static constexpr bool p_refine_only_in_event = true;
  };

  struct registration
      : tt::ConformsTo<Parallel::protocols::RegistrationMetavariables> {
    using element_registrars =
        tmpl::map<tmpl::pair<gh_dg_element_array, dg_registration_list>>;
  };

  using control_system_horizon_metavars =
      control_system::metafunctions::horizon_metavars<control_systems>;
  using control_components =
      control_system::control_components<EvolutionMetavars, control_systems>;

  static void run_deadlock_analysis_simple_actions(
      Parallel::GlobalCache<EvolutionMetavars>& cache,
      const std::vector<std::string>& deadlocked_components) {
    gh::deadlock::run_deadlock_analysis_simple_actions<
        gh_dg_element_array, control_components, interpolation_target_tags,
        tmpl::list<ApparentHorizon>, false>(cache, deadlocked_components);
  }

  using component_list = tmpl::flatten<tmpl::list<
      ::amr::Component<EvolutionMetavars>,
      observers::Observer<EvolutionMetavars>,
      observers::ObserverWriter<EvolutionMetavars>,
      mem_monitor::MemoryMonitor<EvolutionMetavars>,
      importers::ElementDataReader<EvolutionMetavars>, gh_dg_element_array,
      ah::Component<EvolutionMetavars, ApparentHorizon>, control_components,
      tmpl::transform<
          control_system_horizon_metavars,
          tmpl::bind<ah::Component, tmpl::pin<EvolutionMetavars>, tmpl::_1>>,
      tmpl::transform<interpolation_target_tags,
                      tmpl::bind<intrp::InterpolationTarget,
                                 tmpl::pin<EvolutionMetavars>, tmpl::_1>>>>;
};
