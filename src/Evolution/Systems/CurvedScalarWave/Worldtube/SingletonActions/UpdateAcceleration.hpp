// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/TaggedTuple.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {

/*!
 * \brief Computes the time derivative of the evolved variables of the
 * worldtube singleton at this time step.
 *
 * \details If `max_iterations` is 0, the acceleration of the particle will
 * simply be geodesic, see `gr::geodesic_acceleration`. Otherwise, the
 * acceleration due to the scalar self-force is additionally applied to it,
 * see `self_force_acceleration`.
 *
 * At expansion order 2, the constant coefficient of the regular field
 * `Tags::Psi0` is additionally evolved with an ODE derived from the
 * \f$O(\rho^0)\f$ part of the Klein-Gordon equation expanded in inertial
 * coordinates around the moving particle position,
 *
 * \f{align*}{
 * g^{tt}_0 \ddot\Psi^R_0 ={}& \Gamma^t_0 \dot\Psi^R_0
 * + \left(2 g^{tt}_0 \dot x_p^i - 2 g^{it}_0\right) \dot\Psi^R_i
 * + \left(g^{tt}_0 \ddot x_p^i - \dot x_p^i \Gamma^t_0 + \Gamma^i_0\right)
 *   \Psi^R_i
 * - \left(2 g^{tt}_0 \dot x_p^i \dot x_p^j - 4 g^{it}_0 \dot x_p^j
 *   + 2 g^{ij}_0\right) \Psi^R_{ij},
 * \f}
 *
 * where the time derivative of the dipole coefficient is obtained from the
 * projection of the time-derivative field using
 * \f$\dot\Psi^R_i = (\partial_t\Psi^R)_i + 2 \Psi^R_{ij} \dot x_p^j\f$, and
 * the full second-order coefficient is
 * \f$\Psi^R_{ij} = \Psi^R_{\langle ij\rangle} + \delta_{ij}
 * (\Psi^{N,R}_{\langle 0 \rangle} - \Psi^R_0)/\rho^2\f$. At lower expansion
 * orders the time derivatives of `Tags::Psi0` and `Tags::dtPsi0` are set to
 * zero and the values are unused.
 */
struct UpdateAcceleration {
  static constexpr size_t Dim = 3;
  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::EvolvedPosition<Dim>, Tags::EvolvedVelocity<Dim>,
                 Tags::Psi0, Tags::dtPsi0>>;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;
  using return_tags = tmpl::list<dt_variables_tag>;
  using argument_tags = tmpl::list<
      variables_tag, Tags::ParticlePositionVelocity<Dim>,
      Tags::BackgroundQuantities<Dim>, Tags::GeodesicAcceleration<Dim>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Inertial>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim,
                           Frame::Inertial>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Inertial>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                           Frame::Inertial>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Inertial>,
      Tags::Charge, Tags::Mass, Tags::MaxIterations, ::Tags::Time,
      Tags::SelfForceTurnOnTime, Tags::SelfForceTurnOnInterval,
      Tags::ExpansionOrder, Tags::WorldtubeRadius>;
  static void apply(
      gsl::not_null<Variables<
          tmpl::list<::Tags::dt<Tags::EvolvedPosition<Dim>>,
                     ::Tags::dt<Tags::EvolvedVelocity<Dim>>,
                     ::Tags::dt<Tags::Psi0>, ::Tags::dt<Tags::dtPsi0>>>*>
          dt_evolved_vars,
      const Variables<
          tmpl::list<Tags::EvolvedPosition<Dim>, Tags::EvolvedVelocity<Dim>,
                     Tags::Psi0, Tags::dtPsi0>>& evolved_vars,
      const std::array<tnsr::I<double, Dim>, 2>& pos_vel,
      const tuples::TaggedTuple<
          gr::Tags::SpacetimeMetric<double, Dim>,
          gr::Tags::InverseSpacetimeMetric<double, Dim>,
          gr::Tags::SpacetimeChristoffelSecondKind<double, Dim>,
          gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim>,
          Tags::TimeDilationFactor>& background,
      const tnsr::I<double, Dim, Frame::Inertial>& geodesic_acc,
      const Scalar<double>& psi_monopole, const Scalar<double>& dt_psi_monopole,
      const tnsr::i<double, Dim, Frame::Inertial>& psi_dipole,
      const tnsr::i<double, Dim, Frame::Inertial>& dt_psi_dipole,
      const tnsr::ii<double, Dim, Frame::Inertial>& psi_quadrupole,
      double charge, std::optional<double> mass, size_t max_iterations,
      double time, std::optional<double> turn_on_time,
      std::optional<double> turn_on_interval, size_t expansion_order,
      double worldtube_radius);
};

}  // namespace CurvedScalarWave::Worldtube
