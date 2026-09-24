// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/WorldtubeTypeD.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matching.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/ConstraintDamping/GaussianPlusConstant.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube {
namespace Initialization {
/// \brief Add an empty `Tags::KretschmannFaceData` to the DataBox, to be
/// used with `Initialization::Actions::InitializeItems`.
template <size_t Dim>
struct KretschmannFaceData {
  using const_global_cache_tags = tmpl::list<>;
  using mutable_global_cache_tags = tmpl::list<>;
  using simple_tags_from_options = tmpl::list<>;
  using simple_tags = tmpl::list<Tags::KretschmannFaceData<Dim>>;
  using compute_tags = tmpl::list<>;

  using return_tags = simple_tags;
  using argument_tags = tmpl::list<>;

  static void apply(
      const gsl::not_null<worldtube::KretschmannFaceData<Dim>*> data) {
    *data = worldtube::KretschmannFaceData<Dim>{};
  }
};
}  // namespace Initialization

/*!
 * \brief Mutator updating the `Tags::KretschmannFaceData` of the element
 * from the current evolved variables, see `update_kretschmann_face_data()`.
 *
 * \details Apply with `Actions::MutateApply` immediately before
 * `evolution::dg::Actions::ComputeTimeDerivative`, so that the worldtube
 * boundary condition with `PhysicalModel: Quadrupole` finds the face data of
 * the state whose time derivative it corrects. Elements that abut no excision
 * sphere do no work.
 */
template <size_t Dim>
struct UpdateKretschmannFaceData {
  using return_tags = tmpl::list<Tags::KretschmannFaceData<Dim>>;
  using argument_tags =
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, Dim>,
                 gh::Tags::Pi<DataVector, Dim>, gh::Tags::Phi<DataVector, Dim>,
                 domain::Tags::Mesh<Dim>,
                 domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                               Frame::Inertial>,
                 domain::Tags::Element<Dim>, domain::Tags::Domain<Dim>,
                 domain::Tags::MeshVelocity<Dim>, ::Tags::Time,
                 domain::Tags::ExternalBoundaryConditions<Dim>,
                 domain::Tags::Coordinates<Dim, Frame::Inertial>, ::Tags::TimeStepId,
                 gh::Tags::DampingFunctionGamma2<Dim, Frame::Grid>,
                 domain::Tags::FunctionsOfTime>;

  static void apply(
      const gsl::not_null<worldtube::KretschmannFaceData<Dim>*> data,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, Dim, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, Dim, Frame::Inertial>& phi,
      const Mesh<Dim>& mesh,
      const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const Element<Dim>& element, const Domain<Dim>& domain,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          mesh_velocity,
      const double time,
      const std::vector<DirectionMap<
          Dim, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>&
          external_boundary_conditions,
      const tnsr::I<DataVector,Dim,Frame::Inertial>& coordinates,
      const TimeStepId& time_id,
      const ConstraintDamping::DampingFunction<Dim, Frame::Grid>& damping_gamma2,
      const domain::FunctionsOfTimeMap& functions_of_time) {
    update_kretschmann_face_data(
        data, spacetime_metric, pi, phi, mesh, inverse_jacobian, element,
        domain.excision_spheres(), mesh_velocity, time);
    if constexpr (Dim == 3) {
      // The relaxed tidal moments are only needed, and only defined, when the
      // excision face carries the order-two worldtube condition
      if (data->direction.has_value()) {
        const auto& conditions =
            external_boundary_conditions[element.id().block_id()];
        const auto* const worldtube =
            dynamic_cast<const BoundaryConditions::WorldtubeTypeD<Dim>*>(
                conditions.at(*data->direction).get());
        if (worldtube != nullptr and worldtube->radial_response().has_value()) {
          if (mesh_velocity.has_value() or domain.blocks()[element.id().block_id()].is_time_dependent() or
              element.id().refinement_levels() != std::array<size_t,Dim>{}) {
            ERROR("RadialResponse requires a static unrefined spherical-shell element.");
          }
          if (dynamic_cast<const ConstraintDamping::Constant<Dim,Frame::Grid>*>(&damping_gamma2) == nullptr and
              dynamic_cast<const ConstraintDamping::GaussianPlusConstant<Dim,Frame::Grid>*>(&damping_gamma2) == nullptr) {
            ERROR("RadialResponse requires time-independent gamma2.");
          }
          const auto center = BoundaryConditions::detail::excision_sphere_center(
              domain.excision_spheres(), element.id(), time, functions_of_time);
          tnsr::I<DataVector,Dim,Frame::Grid> grid_coordinates(mesh.number_of_grid_points(),0.);
          for(size_t i=0;i<Dim;++i) { grid_coordinates.get(i)=coordinates.get(i); }
          Scalar<DataVector> gamma2(mesh.number_of_grid_points(),0.);
          damping_gamma2(make_not_null(&gamma2),grid_coordinates,time,functions_of_time);
          update_radial_gauge_face_data(make_not_null(&data->radial_gauge),
              spacetime_metric,pi,phi,gamma2,coordinates,center,mesh,inverse_jacobian,*data->direction,time);
        }
        if (worldtube != nullptr and worldtube->face_replay().has_value()) {
          if (mesh_velocity.has_value() or element.id().refinement_levels() != std::array<size_t,Dim>{}) {
            ERROR("Face replay requires a static unrefined spherical-shell element.");
          }
          update_face_replay(make_not_null(&data->replay), *worldtube->face_replay(),
              spacetime_metric, pi, phi, coordinates, mesh, inverse_jacobian, time_id, time);
        } else {
          data->replay.reset();
        }
        if (worldtube != nullptr and
            is_order_two(worldtube->physical_model())) {
          update_filtered_tidal_moments(
              data, spacetime_metric, pi, phi, mesh, inverse_jacobian,
              worldtube->physical_model(), *worldtube->mass(),
              worldtube->moment_relaxation_time(), time);
          return;
        }
      }
      data->filtered_moments.reset();
    } else {
      (void)external_boundary_conditions;
    }
  }
};
}  // namespace gh::worldtube
