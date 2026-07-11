// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/UpdateAcceleration.hpp"

#include <cstddef>
#include <optional>

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SelfForce.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {
void UpdateAcceleration::apply(
    const gsl::not_null<Variables<
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
    const double charge, const std::optional<double> mass,
    const size_t max_iterations, const double time,
    const std::optional<double> turn_on_time,
    const std::optional<double> turn_on_interval, const size_t expansion_order,
    const double worldtube_radius) {
  tnsr::I<double, Dim> self_force_acc(0.);
  const auto& particle_velocity = pos_vel.at(1);
  double roll_on = 0.;
  if (max_iterations > 0 and time > turn_on_time.value()) {
    roll_on =
        turn_on_function(time - turn_on_time.value(), turn_on_interval.value());
    const auto& inverse_metric =
        get<gr::Tags::InverseSpacetimeMetric<double, Dim>>(background);
    const auto& dilation_factor = get<Tags::TimeDilationFactor>(background);
    const double evolved_mass = mass.value() - charge * get(psi_monopole);
    self_force_acceleration(make_not_null(&self_force_acc), dt_psi_monopole,
                            psi_dipole, particle_velocity, charge, evolved_mass,
                            inverse_metric, dilation_factor);
  }
  for (size_t i = 0; i < Dim; ++i) {
    get<::Tags::dt<Tags::EvolvedPosition<Dim>>>(*dt_evolved_vars).get(i)[0] =
        particle_velocity.get(i);
    get<::Tags::dt<Tags::EvolvedVelocity<Dim>>>(*dt_evolved_vars).get(i)[0] =
        geodesic_acc.get(i) + roll_on * self_force_acc.get(i);
  }

  auto& dt_psi0 = get(get<::Tags::dt<Tags::Psi0>>(*dt_evolved_vars));
  auto& dt2_psi0 = get(get<::Tags::dt<Tags::dtPsi0>>(*dt_evolved_vars));
  if (expansion_order < 2) {
    // Psi0 is only evolved at expansion order 2, at lower orders the
    // monopole is fully determined by the boundary data
    dt_psi0 = 0.;
    dt2_psi0 = 0.;
    return;
  }
  // ODE for the constant coefficient of the regular field derived from the
  // O(rho^0) part of the Klein-Gordon equation expanded in inertial
  // coordinates around the moving particle position, see the class
  // documentation.
  const auto& imetric =
      get<gr::Tags::InverseSpacetimeMetric<double, Dim>>(background);
  const auto& trace_christoffel =
      get<gr::Tags::TraceSpacetimeChristoffelSecondKind<double, Dim>>(
          background);
  const auto& psi0 = get(get<Tags::Psi0>(evolved_vars));
  const auto& evolved_dt_psi0 = get(get<Tags::dtPsi0>(evolved_vars));
  dt_psi0 = evolved_dt_psi0;

  // the full coordinate acceleration of the expansion center, including the
  // self-force if it is enabled
  tnsr::I<double, Dim> particle_acceleration = geodesic_acc;
  for (size_t i = 0; i < Dim; ++i) {
    particle_acceleration.get(i) += roll_on * self_force_acc.get(i);
  }

  // full second-order coefficient: STF part plus the trace reconstructed
  // from the boundary monopole and the evolved Psi0
  tnsr::ii<double, Dim> psi_ij = psi_quadrupole;
  const double trace_psi_2_over_3 =
      (get(psi_monopole) - psi0[0]) / square(worldtube_radius);
  for (size_t i = 0; i < Dim; ++i) {
    psi_ij.get(i, i) += trace_psi_2_over_3;
  }

  // time derivative of the dipole coefficient: the projection of the
  // time-derivative field is corrected for the motion of the expansion
  // center
  tnsr::i<double, Dim> psi_dot_i = dt_psi_dipole;
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      psi_dot_i.get(i) += 2. * psi_ij.get(i, j) * particle_velocity.get(j);
    }
  }

  dt2_psi0 = get<0>(trace_christoffel) * evolved_dt_psi0;
  for (size_t i = 0; i < Dim; ++i) {
    dt2_psi0 += (2. * imetric.get(0, 0) * particle_velocity.get(i) -
                 2. * imetric.get(0, i + 1)) *
                psi_dot_i.get(i);
    dt2_psi0 += (imetric.get(0, 0) * particle_acceleration.get(i) -
                 particle_velocity.get(i) * get<0>(trace_christoffel) +
                 trace_christoffel.get(i + 1)) *
                psi_dipole.get(i);
    for (size_t j = 0; j < Dim; ++j) {
      dt2_psi0 -= (2. * imetric.get(0, 0) * particle_velocity.get(i) *
                       particle_velocity.get(j) -
                   4. * imetric.get(0, i + 1) * particle_velocity.get(j) +
                   2. * imetric.get(i + 1, j + 1)) *
                  psi_ij.get(i, j);
    }
  }
  dt2_psi0 /= imetric.get(0, 0);
}
}  // namespace CurvedScalarWave::Worldtube
