// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/ComputeRecordingData.hpp"
#include "ParallelAlgorithms/Interpolation/ComputeExcisionBoundaryVolumeQuantities.hpp"
#include "ParallelAlgorithms/Interpolation/ComputeExcisionBoundaryVolumeQuantities.tpp"

namespace gh::worldtube {
/// Inertial GH tensors and volume derivatives, followed by the existing lapse
/// and target-frame shift diagnostics. The frame of the interpolation points
/// does not change the frame of the recorded GH tensor components.
template <typename TargetFrame>
using excision_recording_tags =
    tmpl::append<recording_tags,
                 tmpl::list<gr::Tags::Lapse<DataVector>,
                            gr::Tags::Shift<DataVector, 3, TargetFrame>>>;

struct ComputeExcisionRecordingData
    : tt::ConformsTo<intrp::protocols::ComputeVarsToInterpolate> {
  using allowed_src_tags = recording_source_tags;
  using required_src_tags = recording_source_tags;
  template <typename TargetFrame>
  using allowed_dest_tags = excision_recording_tags<TargetFrame>;
  template <typename TargetFrame>
  using required_dest_tags = excision_recording_tags<TargetFrame>;

  // The interpolation framework supplies the map Jacobians and mesh velocities
  // when a target-frame shift must be computed. The recording derivatives use
  // the observer's logical-to-inertial inverse Jacobian in either case.
  template <typename SrcTags, typename DestTags, typename... Geometry>
  static void apply(const gsl::not_null<Variables<DestTags>*> result,
                    const Variables<SrcTags>& source, const Mesh<3>& mesh,
                    const Geometry&... geometry) {
    static_assert(std::is_same_v<SrcTags, recording_source_tags>);
    using extra_tags = tmpl::list_difference<DestTags, recording_tags>;
    static_assert(
        std::is_same_v<DestTags, tmpl::append<recording_tags, extra_tags>>);
    const size_t n = source.number_of_grid_points();
    result->initialize(n);
    const size_t recording_size =
        n * Variables<recording_tags>::number_of_independent_components;
    Variables<recording_tags> recording(result->data(), recording_size);
    ComputeRecordingData::apply(make_not_null(&recording), source, mesh);

    Variables<recording_evolved_tags> evolved(n);
    tmpl::for_each<recording_evolved_tags>([&](auto tag_v) {
      using Tag = tmpl::type_from<decltype(tag_v)>;
      get<Tag>(evolved) = get<Tag>(source);
    });
    Variables<extra_tags> extra(result->data() + recording_size,
                                result->size() - recording_size);
    intrp::ComputeExcisionBoundaryVolumeQuantities::apply(
        make_not_null(&extra), evolved, mesh, geometry...);
  }
};
}  // namespace gh::worldtube
