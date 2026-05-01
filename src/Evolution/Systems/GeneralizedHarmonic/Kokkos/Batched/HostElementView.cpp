// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/HostElementView.hpp"

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Batched {

Variables<ah::source_vars<3>> make_ah_source_vars(const HostElementView& view) {
  using spacetime_metric_tag =
      gr::Tags::SpacetimeMetric<DataVector, 3, Frame::Inertial>;
  using pi_tag = gh::Tags::Pi<DataVector, 3, Frame::Inertial>;
  using phi_tag = gh::Tags::Phi<DataVector, 3, Frame::Inertial>;
  using deriv_phi_tag =
      ::Tags::deriv<phi_tag, tmpl::size_t<3>, Frame::Inertial>;

  Variables<ah::source_vars<3>> source_vars{view.mesh.number_of_grid_points()};
  get<spacetime_metric_tag>(source_vars) =
      get<spacetime_metric_tag>(view.evolved_vars);
  get<pi_tag>(source_vars) = get<pi_tag>(view.evolved_vars);
  get<phi_tag>(source_vars) = get<phi_tag>(view.evolved_vars);
  partial_derivative(make_not_null(&get<deriv_phi_tag>(source_vars)),
                     get<phi_tag>(view.evolved_vars), view.mesh,
                     view.logical_to_inertial_inverse_jacobian);
  return source_vars;
}

}  // namespace gh::Batched
