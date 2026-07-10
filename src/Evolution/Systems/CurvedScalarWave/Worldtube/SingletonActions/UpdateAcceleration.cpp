// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/UpdateAcceleration.hpp"

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
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
    const tnsr::I<double, Dim, Frame::Inertial>& geodesic_acc,
    const Scalar<double>& psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
    const tnsr::ii<double, Dim, Frame::Grid>& psi_quadrupole,
    const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole,
    const tnsr::AA<double, Dim, Frame::Grid>& inverse_spacetime_metric,
    const tnsr::A<double, Dim, Frame::Grid>& trace_spacetime_christoffel,
    const size_t expansion_order, const double worldtube_radius,
    const size_t max_iterations) {
  if (UNLIKELY(max_iterations > 0)) {
    ERROR(
        "The scalar self-force is not supported by the grid-frame worldtube "
        "scheme. Select `SelfForce: None`.");
  }
  const auto& particle_velocity = pos_vel.at(1);
  for (size_t i = 0; i < Dim; ++i) {
    get<::Tags::dt<Tags::EvolvedPosition<Dim>>>(*dt_evolved_vars).get(i)[0] =
        particle_velocity.get(i);
    get<::Tags::dt<Tags::EvolvedVelocity<Dim>>>(*dt_evolved_vars).get(i)[0] =
        geodesic_acc.get(i);
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
  // O(rho^0) part of the Klein-Gordon equation expanded in the co-rotating
  // grid frame, see Eq. (39) of https://arxiv.org/abs/2304.05329
  const auto& psi0 = get(get<Tags::Psi0>(evolved_vars));
  const auto& evolved_dt_psi0 = get(get<Tags::dtPsi0>(evolved_vars));
  dt_psi0 = evolved_dt_psi0;
  double trace_inverse_spatial_metric = 0.;
  dt2_psi0 = get<0>(trace_spacetime_christoffel) * evolved_dt_psi0;
  for (size_t i = 0; i < Dim; ++i) {
    dt2_psi0 -=
        2. * inverse_spacetime_metric.get(0, i + 1) * dt_psi_dipole.get(i);
    dt2_psi0 += trace_spacetime_christoffel.get(i + 1) * psi_dipole.get(i);
    trace_inverse_spatial_metric += inverse_spacetime_metric.get(i + 1, i + 1);
    for (size_t j = 0; j < Dim; ++j) {
      dt2_psi0 -= 2. * inverse_spacetime_metric.get(i + 1, j + 1) *
                  psi_quadrupole.get(i, j);
    }
  }
  dt2_psi0 -= 2. * trace_inverse_spatial_metric * (get(psi_monopole) - psi0) /
              (worldtube_radius * worldtube_radius);
  dt2_psi0 /= inverse_spacetime_metric.get(0, 0);
}
}  // namespace CurvedScalarWave::Worldtube
