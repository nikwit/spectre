// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/SphereReplay.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "Time/Tags/Time.hpp"

namespace gh::BoundaryConditions {
/// All-reference moving-sphere boundary. Native DG packaging and lifting handle
/// mesh velocity and all characteristic families; no model sectors are
/// replaced.
class WorldtubeReplay final : public BoundaryCondition<3> {
 public:
  struct H5Filename {
    using type = std::string;
    static constexpr Options::String help = "GhWorldtube HDF5 tape.";
  };
  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help =
        "Sphere volume subfile, usually GhWorldtube.";
  };
  struct Radius {
    using type = double;
    static constexpr Options::String help =
        "Selected constant inertial coordinate radius.";
    static double lower_bound() { return 0.; }
  };
  using options = tmpl::list<H5Filename, SubfileName, Radius>;
  static constexpr Options::String help =
      "Replay one recorded inertial sphere through the native GH numerical "
      "flux.";
  WorldtubeReplay() = default;
  WorldtubeReplay(std::string filename, std::string subfile, double radius,
                  bool read_physical_data = false);
  explicit WorldtubeReplay(CkMigrateMessage* msg);
  const worldtube::SphereReplay& tape() const { return *tape_; }
  WRAPPED_PUPable_decl_template(WorldtubeReplay);
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override;
  void pup(PUP::er& p) override;
  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;
  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags =
      tmpl::list<domain::Tags::Coordinates<3, Frame::Inertial>,
                 gh::Tags::ConstraintGamma1, gh::Tags::ConstraintGamma2>;
  using dg_gridless_tags =
      tmpl::list<::Tags::Time, domain::Tags::Domain<3>,
                 domain::Tags::Element<3>, domain::Tags::FunctionsOfTime>;
  std::optional<std::string> dg_ghost(
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, 3, Frame::Inertial>*> pi,
      gsl::not_null<tnsr::iaa<DataVector, 3, Frame::Inertial>*> phi,
      gsl::not_null<Scalar<DataVector>*> gamma1,
      gsl::not_null<Scalar<DataVector>*> gamma2,
      gsl::not_null<Scalar<DataVector>*> lapse,
      gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> shift,
      gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
          inv_spatial_metric,
      const std::optional<
          tnsr::I<DataVector, 3, Frame::Inertial>>& /*face_mesh_velocity*/,
      const tnsr::i<DataVector, 3, Frame::Inertial>& /*normal_covector*/,
      const tnsr::I<DataVector, 3, Frame::Inertial>& /*normal_vector*/,
      const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
      const Scalar<DataVector>& interior_gamma1,
      const Scalar<DataVector>& interior_gamma2, double time,
      const Domain<3>& domain, const Element<3>& element,
      const domain::FunctionsOfTimeMap& functions_of_time) const;

 private:
  std::string filename_{};
  std::string subfile_{};
  double radius_{0.};
  bool read_physical_data_{false};
  std::shared_ptr<const worldtube::SphereReplay> tape_{};
};
}  // namespace gh::BoundaryConditions
