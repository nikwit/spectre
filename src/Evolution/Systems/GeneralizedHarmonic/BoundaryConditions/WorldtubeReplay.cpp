// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeReplay.hpp"

#include <utility>
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/Serialization/PupStlCpp11.hpp"

namespace gh::BoundaryConditions {
WorldtubeReplay::WorldtubeReplay(std::string filename, std::string subfile,
                                 const double radius,
                                 const bool read_physical_data)
    : filename_(std::move(filename)),
      subfile_(std::move(subfile)),
      radius_(radius),
      read_physical_data_(read_physical_data),
      tape_(worldtube::cached_sphere_replay(filename_, subfile_, radius_,
                                            read_physical_data_)) {}
WorldtubeReplay::WorldtubeReplay(CkMigrateMessage* msg)
    : BoundaryCondition<3>(msg) {}
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
WorldtubeReplay::get_clone() const {
  return std::make_unique<WorldtubeReplay>(*this);
}
void WorldtubeReplay::pup(PUP::er& p) {
  BoundaryCondition<3>::pup(p);
  p | filename_;
  p | subfile_;
  p | radius_;
  p | read_physical_data_;
  if (p.isUnpacking())
    tape_ = worldtube::cached_sphere_replay(filename_, subfile_, radius_,
                                            read_physical_data_);
}
std::optional<std::string> WorldtubeReplay::dg_ghost(
    gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> spacetime_metric,
    gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
    gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> phi,
    gsl::not_null<Scalar<DataVector>*> gamma1,
    gsl::not_null<Scalar<DataVector>*> gamma2,
    gsl::not_null<Scalar<DataVector>*> lapse,
    gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
    gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*> inv_spatial_metric,
    const std::optional<
        tnsr::I<DataVector, 3, Frame::Inertial>>& /*face_mesh_velocity*/,
    const tnsr::i<DataVector, 3, Frame::Inertial>& /*normal_covector*/,
    const tnsr::I<DataVector, 3, Frame::Inertial>& /*normal_vector*/,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const Scalar<DataVector>& interior_gamma1,
    const Scalar<DataVector>& interior_gamma2, double time,
    const Domain<3>& domain, const Element<3>& element,
    const domain::FunctionsOfTimeMap& functions_of_time) const {
  const auto center_at_time = [&](const double sample_time) {
    return detail::excision_sphere_center(domain.excision_spheres(),
                                          element.id(), sample_time,
                                          functions_of_time);
  };
  const auto values = tape_->evaluate(time, coords, center_at_time);
  *spacetime_metric = get<gr::Tags::SpacetimeMetric<DataVector, 3>>(values);
  *pi = get<gh::Tags::Pi<DataVector, 3>>(values);
  *phi = get<gh::Tags::Phi<DataVector, 3>>(values);
  *gamma1 = interior_gamma1;
  *gamma2 = interior_gamma2;
  const auto spatial_metric = gr::spatial_metric(*spacetime_metric);
  determinant_and_inverse(lapse, inv_spatial_metric, spatial_metric);
  gr::shift(shift, *spacetime_metric, *inv_spatial_metric);
  gr::lapse(lapse, *shift, *spacetime_metric);
  return {};
}
PUP::able::PUP_ID WorldtubeReplay::my_PUP_ID = 0;  // NOLINT
}  // namespace gh::BoundaryConditions
