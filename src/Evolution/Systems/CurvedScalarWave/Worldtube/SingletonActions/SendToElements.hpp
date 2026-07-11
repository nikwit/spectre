// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct TimeStepId;
}  // namespace Tags
/// \endcond

namespace CurvedScalarWave::Worldtube::Actions {
/*!
 * \brief Sends the coefficients of the Taylor expansion of the regular field
 * to each element abutting the worldtube.
 *
 * \details The expansion is defined in the inertial frame, centered on the
 * current particle position. The coefficients are ordered as the monopole,
 * the three dipole components and, at second order, the six independent
 * quadrupole components (xx, xy, xz, yy, yz, zz), first for \f$\Psi^R\f$ and
 * then for \f$\partial_t \Psi^R\f$.
 *
 * At second order the constant coefficient \f$\Psi^R_0\f$ is evolved
 * separately (see `UpdateAcceleration`) because the monopole of the regular
 * field on the worldtube boundary also contains the trace of the second-order
 * coefficient, which is reconstructed here, see Eq. (10a) of the matching
 * scheme. Because the expansion center moves with the particle, the constant
 * coefficient of the time-derivative field differs from the time derivative
 * of the constant coefficient by the motion of the center,
 *
 * \f{equation}{
 * (\partial_t \Psi^R)_0 = \dot\Psi^R_0 - \Psi^R_i \dot x_p^i,
 * \f}
 *
 * which is used to reconstruct the trace of the second-order coefficient of
 * the time-derivative field.
 */
template <typename Metavariables>
struct SendToElements {
  static constexpr size_t Dim = Metavariables::volume_dim;
  using psi_tag = CurvedScalarWave::Tags::Psi;
  using dt_psi_tag = ::Tags::dt<CurvedScalarWave::Tags::Psi>;
  using tags_to_send = tmpl::list<psi_tag, dt_psi_tag>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const size_t order = db::get<Tags::ExpansionOrder>(box);
    auto& element_proxies = Parallel::get_parallel_component<
        typename Metavariables::dg_element_array>(cache);
    const auto& faces_grid_coords =
        get<Tags::ElementFacesGridCoordinates<Dim>>(box);
    const auto& psi_l0 =
        get<Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Inertial>>(
            box);
    const auto& dt_psi_l0 =
        get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim,
                                 Frame::Inertial>>(box);
    const auto& psi_l1 =
        get<Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Inertial>>(
            box);
    const auto& dt_psi_l1 =
        get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                                 Frame::Inertial>>(box);
    const size_t num_coefs = order == 0 ? 1 : (order == 1 ? 4 : 10);
    Variables<tags_to_send> vars_to_send(num_coefs);
    DataVector& psi_coefs = get(get<psi_tag>(vars_to_send));
    DataVector& dt_psi_coefs = get(get<dt_psi_tag>(vars_to_send));
    psi_coefs[0] = get(psi_l0);
    dt_psi_coefs[0] = get(dt_psi_l0);
    if (order > 0) {
      for (size_t i = 0; i < Dim; ++i) {
        psi_coefs[i + 1] = psi_l1.get(i);
        dt_psi_coefs[i + 1] = dt_psi_l1.get(i);
      }
    }
    if (order > 1) {
      const auto& psi_l2 = get<
          Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Inertial>>(
          box);
      const auto& dt_psi_l2 =
          get<Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 2, Dim,
                                   Frame::Inertial>>(box);
      const auto& psi_0 = get<Tags::Psi0>(box);
      const auto& dt_psi_0 = get<Tags::dtPsi0>(box);
      const auto& particle_velocity =
          get<Tags::ParticlePositionVelocity<Dim>>(box)[1];
      const double wt_radius = db::get<Tags::WorldtubeRadius>(box);

      // the constant coefficient of the regular field is evolved separately
      // at second order, see the class documentation
      psi_coefs[0] = get(psi_0)[0];
      // the constant coefficient of the time-derivative field gets a
      // correction from the motion of the expansion center
      double dt_psi_coef_0 = get(dt_psi_0)[0];
      for (size_t i = 0; i < Dim; ++i) {
        dt_psi_coef_0 -= psi_l1.get(i) * particle_velocity.get(i);
      }
      dt_psi_coefs[0] = dt_psi_coef_0;

      // reconstruct the traces of the second-order coefficients from the
      // boundary monopoles, see Eq. (10a) of the matching scheme
      const double trace_psi_2_over_3 =
          (get(psi_l0) - psi_coefs[0]) / square(wt_radius);
      const double trace_dt_psi_2_over_3 =
          (get(dt_psi_l0) - dt_psi_coefs[0]) / square(wt_radius);
      size_t index = 4;
      for (size_t i = 0; i < Dim; ++i) {
        for (size_t j = i; j < Dim; ++j, ++index) {
          psi_coefs[index] =
              psi_l2.get(i, j) + (i == j ? trace_psi_2_over_3 : 0.);
          dt_psi_coefs[index] =
              dt_psi_l2.get(i, j) + (i == j ? trace_dt_psi_2_over_3 : 0.);
        }
      }
      ASSERT(index == num_coefs, "Internal indexing error");
    }
    for (const auto& [element_id, _] : faces_grid_coords) {
      auto vars_to_send_copy = vars_to_send;
      Parallel::receive_data<Tags::RegularFieldInbox<Dim>>(
          element_proxies[element_id], db::get<::Tags::TimeStepId>(box),
          std::move(vars_to_send_copy));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
