// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/ComputeVarsToInterpolate.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube {
namespace Tags {
/// Vacuum curvature invariant C_abcd C^abcd, computed with the same convention
/// as the worldtube matching models (8(E_ij E^ij - B_ij B^ij)).
struct Kretschmann : db::SimpleTag {
  using type = Scalar<DataVector>;
};
}  // namespace Tags

using recording_evolved_tags =
    tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
               gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>;
template <typename Tag>
using recording_deriv = ::Tags::deriv<Tag, tmpl::size_t<3>, Frame::Inertial>;
using recording_tags = tmpl::append<
    recording_evolved_tags,
    tmpl::list<recording_deriv<gr::Tags::SpacetimeMetric<DataVector, 3>>,
               recording_deriv<gh::Tags::Pi<DataVector, 3>>,
               recording_deriv<gh::Tags::Phi<DataVector, 3>>, Tags::Kretschmann,
               recording_deriv<Tags::Kretschmann>>>;
using recording_inverse_jacobian =
    ::Events::Tags::ObserverInverseJacobian<3, Frame::ElementLogical,
                                            Frame::Inertial>;
using recording_source_tags =
    tmpl::push_back<recording_evolved_tags, recording_inverse_jacobian>;

/// Compute all recording data on the donor element's volume mesh. Inertial
/// spatial derivatives use the full logical-to-inertial inverse Jacobian at
/// the observation time, including the moving maps. K is computed from the
/// volume derivatives of Pi and Phi, then differentiated on that same mesh.
/// Only the resulting fields are interpolated to the recording sphere.
void compute_recording_data(
    gsl::not_null<Variables<recording_tags>*> result,
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian);

struct ComputeRecordingData
    : tt::ConformsTo<intrp::protocols::ComputeVarsToInterpolate> {
  using allowed_src_tags = recording_source_tags;
  using required_src_tags = recording_source_tags;
  template <typename Frame>
  using allowed_dest_tags = recording_tags;
  template <typename Frame>
  using required_dest_tags = recording_tags;

  template <typename SrcTags, typename DestTags>
  static void apply(const gsl::not_null<Variables<DestTags>*> result,
                    const Variables<SrcTags>& source, const Mesh<3>& mesh) {
    static_assert(std::is_same_v<SrcTags, recording_source_tags>);
    static_assert(std::is_same_v<DestTags, recording_tags>);
    compute_recording_data(
        result, get<gr::Tags::SpacetimeMetric<DataVector, 3>>(source),
        get<gh::Tags::Pi<DataVector, 3>>(source),
        get<gh::Tags::Phi<DataVector, 3>>(source), mesh,
        get<recording_inverse_jacobian>(source));
  }
};
}  // namespace gh::worldtube
