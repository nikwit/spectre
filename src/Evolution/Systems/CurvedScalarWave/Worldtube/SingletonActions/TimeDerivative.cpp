// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/TimeDerivative.hpp"

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {

void TimeDerivativeMutator::apply(
    const gsl::not_null<Variables<
        tmpl::list<::Tags::dt<Tags::Psi0>, ::Tags::dt<Tags::dtPsi0>,
                   ::Tags::dt<Tags::Position>, ::Tags::dt<Tags::Velocity>>>*>
        dt_evolved_vars,
    const Variables<tmpl::list<Tags::Psi0, Tags::dtPsi0, Tags::Position,
                               Tags::Velocity>>& evolved_vars,
    const Scalar<double>& psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
    const tnsr::ii<double, Dim, Frame::Grid>& psi_quadrupole,
    const Scalar<double>& dt_psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole,
    const tnsr::AA<double, Dim, Frame::Grid>& inverse_spacetime_metric,
    const tnsr::A<double, Dim, Frame::Grid>& trace_spacetime_christoffel,
    const ExcisionSphere<Dim>& excision_sphere, const double time,
    const std::array<double, 2>& worldtube_radius_and_velocity,
    const double mass, const double charge,
    const gr::Solutions::KerrSchild& kerr_schild) {
  const double wt_radius = worldtube_radius_and_velocity.at(0);
  const auto& psi0 = get(get<Tags::Psi0>(evolved_vars));
  const auto& dt_psi0 = get(get<Tags::dtPsi0>(evolved_vars));
  get(get<::Tags::dt<Tags::Psi0>>(*dt_evolved_vars)) = dt_psi0;
  double trace_inverse_spatial_metric = 0.;
  auto& dt2_psi0 = get(get<::Tags::dt<Tags::dtPsi0>>(*dt_evolved_vars));
  dt2_psi0 = get<0>(trace_spacetime_christoffel) * dt_psi0;
  for (size_t i = 0; i < Dim; ++i) {
    dt2_psi0 -=
        2. * inverse_spacetime_metric.get(0, i + 1) * dt_psi_dipole.get(i);
    dt2_psi0 += trace_spacetime_christoffel.get(i + 1) * psi_dipole.get(i);
    trace_inverse_spatial_metric += inverse_spacetime_metric.get(i + 1, i + 1);
    for (size_t j = 0; j < 3; ++j) {
      dt2_psi0 -= 2. * inverse_spacetime_metric.get(i + 1, j + 1) *
                  psi_quadrupole.get(i, j);
    }
  }
  dt2_psi0 -= 2. * trace_inverse_spatial_metric * (get(psi_monopole) - psi0) /
              (wt_radius * wt_radius);
  dt2_psi0 /= inverse_spacetime_metric.get(0, 0);
  tnsr::I<double, 3> inertial_particle_position{};
  tnsr::I<double, 3> particle_velocity{};
  for (size_t i = 0; i < 3; ++i) {
    inertial_particle_position.get(i) =
        get<Tags::Position>(evolved_vars).get(i)[0];
    particle_velocity.get(i) = get<Tags::Velocity>(evolved_vars).get(i)[0];
  }

  const auto spacetime_vars = kerr_schild.variables(
      inertial_particle_position, time,
      tmpl::list<
          gr::Tags::Lapse<double>,
          gr::Tags::Shift<double, Dim, Frame::Inertial>,
          gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
          gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>,
          ::Tags::dt<gr::Tags::Lapse<double>>,
          ::Tags::deriv<gr::Tags::Lapse<double>, tmpl::size_t<Dim>,
                        Frame::Inertial>,
          ::Tags::dt<gr::Tags::Shift<double, Dim, Frame::Inertial>>,
          ::Tags::deriv<gr::Tags::Shift<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>,
          ::Tags::dt<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>,
          ::Tags::deriv<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>{});
  const auto inverse_spacetime_metric_inertial = gr::inverse_spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars));
  const auto spacetime_metric_inertial = gr::spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars));
  const auto d_spacetime_metric = gr::derivatives_of_spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<::Tags::dt<gr::Tags::Lapse<double>>>(spacetime_vars),
      get<::Tags::deriv<gr::Tags::Lapse<double>, tmpl::size_t<Dim>,
                        Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<::Tags::dt<gr::Tags::Shift<double, Dim, Frame::Inertial>>>(
          spacetime_vars),
      get<::Tags::deriv<gr::Tags::Shift<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars),
      get<::Tags::dt<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>>(
          spacetime_vars),
      get<::Tags::deriv<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>(spacetime_vars));

  const auto christoffel = gr::christoffel_second_kind(
      d_spacetime_metric, inverse_spacetime_metric_inertial);

  tnsr::I<double, Dim> particle_acceleration{};
  double u0_squared = spacetime_metric_inertial.get(0, 0);
  for (size_t i = 0; i < Dim; ++i) {
    particle_acceleration.get(i) =
        particle_velocity.get(i) * christoffel.get(0, 0, 0) -
        christoffel.get(i + 1, 0, 0);
    u0_squared +=
        2. * spacetime_metric_inertial.get(i + 1, 0) * particle_velocity.get(i);

    for (size_t j = 0; j < Dim; ++j) {
      u0_squared += spacetime_metric_inertial.get(i + 1, j + 1) *
                    particle_velocity.get(i) * particle_velocity.get(j);
      particle_acceleration.get(i) +=
          2. * particle_velocity.get(j) *
          (particle_velocity.get(i) * christoffel.get(0, j + 1, 0) -
           christoffel.get(i + 1, j + 1, 0));
      for (size_t k = 0; k < Dim; ++k) {
        particle_acceleration.get(i) +=
            particle_velocity.get(j) * particle_velocity.get(k) *
            (particle_velocity.get(i) * christoffel.get(0, j + 1, k + 1) -
             christoffel.get(i + 1, j + 1, k + 1));
      }
    }
  }

  const double turn_up_time = 1200.;
  if (time > turn_up_time) {
    ::InverseJacobian<double, Dim, Frame::Grid, Frame::Inertial> inv_jacobian{};
    const double angle = atan2(inertial_particle_position.get(1),
                               inertial_particle_position.get(0));
    inv_jacobian.get(0, 0) = cos(angle);
    inv_jacobian.get(0, 1) = sin(angle);
    inv_jacobian.get(1, 0) = -sin(angle);
    inv_jacobian.get(1, 1) = cos(angle);
    inv_jacobian.get(2, 2) = 1.;

    tnsr::I<double, Dim> di_psi_inertial{};
    for (size_t i = 0; i < Dim; ++i) {
      di_psi_inertial.get(i) = get<0>(psi_dipole) * inv_jacobian.get(0, i) +
                               get<1>(psi_dipole) * inv_jacobian.get(1, i) +
                               get<2>(psi_dipole) * inv_jacobian.get(2, i);
    }
    double v_dot_di_psi = 0.;
    for (size_t i = 0; i < Dim; ++i) {
      v_dot_di_psi += particle_velocity.get(i) * di_psi_inertial.get(i);
    }
    u0_squared = -1. / u0_squared;
    const double t_minus_turnup = time - turn_up_time;
    const double roll_on = t_minus_turnup < 800. ? t_minus_turnup / 800. : 1.;
    for (size_t i = 0; i < Dim; ++i) {
      particle_acceleration.get(i) +=
          (inverse_spacetime_metric_inertial.get(i + 1, 0) -
           particle_velocity.get(i) *
               inverse_spacetime_metric_inertial.get(0, 0)) *
          (get(dt_psi_monopole) - v_dot_di_psi) * roll_on * charge / mass /
          u0_squared;
      for (size_t j = 0; j < Dim; ++j) {
        particle_acceleration.get(i) +=
            (inverse_spacetime_metric_inertial.get(i + 1, j + 1) -
             particle_velocity.get(i) *
                 inverse_spacetime_metric_inertial.get(0, j + 1)) *
            di_psi_inertial.get(j) * roll_on * charge / mass / u0_squared;
      }
    }
  }
  for (size_t i = 0; i < Dim; ++i) {
    get<::Tags::dt<Tags::Position>>(*dt_evolved_vars).get(i)[0] =
        get<Tags::Velocity>(evolved_vars).get(i)[0];
    get<::Tags::dt<Tags::Velocity>>(*dt_evolved_vars).get(i)[0] =
        particle_acceleration.get(i);
  }
}

}  // namespace CurvedScalarWave::Worldtube
