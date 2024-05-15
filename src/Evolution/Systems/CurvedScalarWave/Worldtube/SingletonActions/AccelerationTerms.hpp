// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {

/// @{
/*!
 * \brief Calculates the time derivative of `Psi0`, the constant coefficient of
 * the expansion of `Psi`.
 *
 * \details The derivation comes from expanding the scalar wave equation to
 * second order and reads
 *
 * \f{equation}
 * g^{00}_0 \ddot{\Psi}^R_0(t_s) + 2 g_0^{0i}
 * \dot{\Psi}^N_i(t_s) + 2 g_0^{ij} \Psi^N_{\langle ij \rangle}(t_s)
 * + \frac{2 \delta_{ij}  g_0^{ij}}{R^2}
 * \left(\Psi^N_0(t_s) - \Psi^R_0(t_s) \right) -
 * \Gamma_0^0\dot{\Psi}R_0(t_s) -  \Gamma_0^i \Psi_i^N(t_s) = 0.
 *\f}
 *
 * Here, \f$ \Gamma^\mu_0 \f$ and \f$ g^{\mu \nu}_0 \f$ are the trace of the
 * spacetime Christoffel symbol and the inverse spacetime metric, respectively,
 * evaluated at the position of the particle; \f$\Psi^N_0\f$, \f$\Psi^N_i\f$,
 * \f$\Psi^N_\langle i j \rangle\f$ are the monopole, dipole and quadrupole of
 * the regular field on the worldtube boundary transformed to symmetric
 * trace-free tensors and \f$ R\f$ is the worldtube radius.
 */
struct AccelerationTermsMutator {
  static constexpr size_t Dim = 3;

  using return_tags = tmpl::list<Tags::SelfForce>;
  using argument_tags = tmpl::list<
      Tags::Position, Tags::Velocity,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim, Frame::Grid>,
      ::Tags::Time, Tags::ParticleMass, Tags::ParticleCharge, Tags::TurnOnTime,
      Tags::TurnOnInterval, Tags::CurrentIteration,
      CurvedScalarWave::Tags::BackgroundSpacetime<gr::Solutions::KerrSchild>>;

  static void apply(const gsl::not_null<Scalar<DataVector>*> self_force,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& position,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& velocity,
                    const Scalar<double>& psi_monopole,
                    const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
                    const Scalar<double>& dt_psi_monopole,
                    const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole,
                    const double time, const double mass, const double charge,
                    const double turn_on_time, const double turn_on_interval,
                    const size_t iteration,
                    const gr::Solutions::KerrSchild& kerr_schild);
};

namespace Actions {
struct ComputeAccelerationTerms {
  static constexpr size_t Dim = 3;
  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::Psi0, Tags::dtPsi0, Tags::Position, Tags::Velocity>>;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
  using simple_tags = tmpl::list<
      dt_variables_tag, variables_tag,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim, Frame::Grid>,
      Tags::ParticlePositionVelocity<3>, Tags::SelfForce>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* /*meta*/) {
    db::mutate_apply<AccelerationTermsMutator>(make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
/// @}

}  // namespace Actions
}  // namespace CurvedScalarWave::Worldtube
