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
#include "Time/Tags/Time.hpp"
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
                 domain::Tags::ExternalBoundaryConditions<Dim>>;

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
          external_boundary_conditions) {
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
