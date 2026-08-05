// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <type_traits>

#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Translation.hpp"
#include "ControlSystem/Measurements/GaussBonnetCenter.hpp"
#include "ControlSystem/Measurements/NonFactoryCreatable.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "ControlSystem/Tags/QueueTags.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Actions/UpdateMessageQueue.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/ComputeVarsToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
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
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Time/Tags/TimeAndPrevious.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::single_bh_gauss_bonnet {

using diagnostic_frame = Frame::Inertial;

// Compute the electric and magnetic Weyl tensors at volume points before
// interpolating them to the control sphere.  This follows the same route used
// by the binary-black-hole Gauss--Bonnet tracker, and avoids differentiating
// data after interpolation.
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

struct ComputeVolumeQuantities
    : tt::ConformsTo<intrp::protocols::ComputeVarsToInterpolate> {
  using allowed_src_tags = volume_source_tags;
  using required_src_tags = volume_source_tags;

  template <typename TargetFrame>
  using allowed_dest_tags =
      tmpl::conditional_t<std::is_same_v<TargetFrame, diagnostic_frame>,
                          volume_target_tags, tmpl::list<>>;
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

    Scalar<DataVector> det_spatial_metric{src_vars.number_of_grid_points()};
    Scalar<DataVector> lapse{src_vars.number_of_grid_points()};
    tnsr::I<DataVector, 3, diagnostic_frame> shift{
        src_vars.number_of_grid_points()};
    tnsr::AA<DataVector, 3, diagnostic_frame> inverse_spacetime_metric{
        src_vars.number_of_grid_points()};
    tnsr::A<DataVector, 3, diagnostic_frame> spacetime_normal_vector{
        src_vars.number_of_grid_points()};
    tnsr::ii<DataVector, 3, diagnostic_frame> extrinsic_curvature{
        src_vars.number_of_grid_points()};
    Scalar<DataVector> sqrt_det_spatial_metric{
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

struct Measurement : tt::ConformsTo<control_system::protocols::Measurement> {
  struct GaussBonnet
      : tt::ConformsTo<control_system::protocols::Submeasurement> {
    static std::string name() { return Measurement::name() + "::GaussBonnet"; }

   private:
    template <typename ControlSystems>
    struct InterpolationTarget
        : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
      static std::string name() { return "ControlSystemGaussBonnetCenter"; }

      using temporal_id = ::Tags::TimeAndPrevious<0>;
      using vars_to_interpolate_to_target = volume_target_tags;
      using compute_vars_to_interpolate = ComputeVolumeQuantities;
      using compute_items_on_source =
          tmpl::list<::Tags::TimeAndPreviousCompute<0>>;
      using compute_items_on_target = tmpl::list<
          gr::Tags::WeylElectricScalarCompute<DataVector, 3, Frame::Inertial>,
          gr::Tags::WeylMagneticScalarCompute<DataVector, 3, Frame::Inertial>,
          gr::Tags::GaussBonnetScalarCompute<DataVector>>;
      using compute_target_points =
          intrp::TargetPoints::Sphere<InterpolationTarget, ::Frame::Grid>;
      using post_interpolation_callbacks =
          tmpl::list<control_system::RunCallbacks<GaussBonnet, ControlSystems>>;

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
            3, InterpolationTarget<ControlSystems>, volume_source_tags>>;
  };

  static std::string name() { return "GaussBonnetCenter"; }
  using submeasurements = tmpl::list<GaussBonnet>;
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
  using control_error = control_system::ControlErrors::Translation<1>;

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
    static void apply(typename Measurement::GaussBonnet /*submeasurement*/,
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
          control_system::measurements::gauss_bonnet_center(
              strahlkorper, gauss_bonnet_scalar));
    }
  };
};

}  // namespace gh::single_bh_gauss_bonnet
