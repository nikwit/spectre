// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
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
 * \details The particle moves on a geodesic, see
 * `gr::geodesic_acceleration`. The scalar self-force is not supported by the
 * grid-frame worldtube scheme and requesting iterations of the acceleration
 * is an error.
 *
 * At expansion order 2, the constant coefficient of the regular field
 * `Tags::Psi0` is additionally evolved with an ODE derived from an expansion
 * of the Klein-Gordon equation in the co-rotating grid frame, see Eq. (39)
 * of \cite Wittek:2023nyi. At lower expansion orders its time derivatives
 * are set to zero.
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
      Tags::GeodesicAcceleration<Dim>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<Tags::PsiWorldtube, 2, Dim, Frame::Grid>,
      Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim, Frame::Grid>,
      gr::Tags::InverseSpacetimeMetric<double, Dim, Frame::Grid>,
      gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim, Frame::Grid>,
      Tags::ExpansionOrder, Tags::WorldtubeRadius, Tags::MaxIterations>;
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
      const tnsr::I<double, Dim, Frame::Inertial>& geodesic_acc,
      const Scalar<double>& psi_monopole,
      const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
      const tnsr::ii<double, Dim, Frame::Grid>& psi_quadrupole,
      const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole,
      const tnsr::AA<double, Dim, Frame::Grid>& inverse_spacetime_metric,
      const tnsr::A<double, Dim, Frame::Grid>& trace_spacetime_christoffel,
      size_t expansion_order, double worldtube_radius, size_t max_iterations);
};

}  // namespace CurvedScalarWave::Worldtube
