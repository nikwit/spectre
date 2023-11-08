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
    const gsl::not_null<Scalar<DataVector>*> self_force,
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
    const double mass, const double charge, const double turn_on_time,
    const double turn_on_interval, const size_t expansion_order,
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
  get(*self_force) = DataVector(9, 0.);
  if (time > turn_on_time) {
    u0_squared = -1. / u0_squared;
    const double t_minus_turnup = time - turn_on_time;
    double roll_on =
        t_minus_turnup < turn_on_interval
            ? t_minus_turnup /
                  turn_on_interval  // square(sin(M_PI_2 * t_minus_turnup /
                                    // turn_on_interval))
            : 1.;

    for (size_t i = 0; i < Dim; ++i) {
      particle_acceleration.get(i) +=
          (inverse_spacetime_metric_inertial.get(i + 1, 0) -
           particle_velocity.get(i) *
               inverse_spacetime_metric_inertial.get(0, 0)) *
          get(dt_psi_monopole) * roll_on * charge / mass / u0_squared;
      if (expansion_order > 0) {
        for (size_t j = 0; j < Dim; ++j) {
          particle_acceleration.get(i) +=
              (inverse_spacetime_metric_inertial.get(i + 1, j + 1) -
               particle_velocity.get(i) *
                   inverse_spacetime_metric_inertial.get(0, j + 1)) *
              psi_dipole.get(j) * roll_on * charge / mass / u0_squared;
        }
      }
    }
    const double u0 = sqrt(u0_squared);
    const tnsr::A<double, Dim> four_velocity{
        {u0, get<0>(particle_velocity) * u0, get<1>(particle_velocity) * u0,
         get<2>(particle_velocity) * u0}};
    const tnsr::a<double, Dim> d_psiR{{get(dt_psi_monopole), get<0>(psi_dipole),
                                       get<1>(psi_dipole), get<2>(psi_dipole)}};
    const auto contracted_christoffel_inertial =
        trace_last_indices(christoffel, inverse_spacetime_metric_inertial);
    double dt2_psiR = -get<0>(contracted_christoffel_inertial) * get<0>(d_psiR);
    for (size_t i = 0; i < Dim; ++i) {
      dt2_psiR +=
          2. * inverse_spacetime_metric_inertial.get(0, i + 1) *
              dt_psi_dipole.get(i) -
          contracted_christoffel_inertial.get(i + 1) * d_psiR.get(i + 1);
      for (size_t j = 0; j < Dim; ++j) {
      }
    }

    dt2_psiR /= get<0, 0>(inverse_spacetime_metric_inertial);
    const auto& pos = inertial_particle_position;
    const auto& vel = particle_velocity;
    const auto& acc = particle_acceleration;

    const double r = magnitude(pos).get();
    tnsr::a<double, Dim> dt_d_psiR{{dt2_psiR, get<0>(dt_psi_dipole),
                                    get<1>(dt_psi_dipole),
                                    get<2>(dt_psi_dipole)}};
    const auto& metric = spacetime_metric_inertial;
    tnsr::A<double, Dim> dt_four_velocity{};
    tnsr::iaa<double, Dim> d_metric{};
    tnsr::iAA<double, Dim> d_imetric{};
    tnsr::ijAA<double, Dim> d2_metric{};
    tnsr::ii<double, Dim> delta_ll{0.};
    tnsr::Ij<double, Dim> delta_ul{0.};
    tnsr::i<double, Dim> pos_lower{};

    for (size_t i = 0; i < Dim; ++i) {
      delta_ll.get(i, i) = 1.;
      delta_ul.get(i, i) = 1.;
      pos_lower.get(i) = pos.get(i);
    }

    const auto d_imetric_ij = tenex::evaluate<ti::i, ti::J, ti::K>(
        6. * pos(ti::J) * pos(ti::K) * pos_lower(ti::i) /
            (square(r) * cube(r)) -
        2. * delta_ul(ti::J, ti::i) * pos(ti::K) / cube(r) -
        2. * delta_ul(ti::K, ti::i) * pos(ti::J) / cube(r));
    const auto d_imetric_i0 = tenex::evaluate<ti::i, ti::J>(
        -4. * pos_lower(ti::i) * pos(ti::J) / square(square(r)) +
        2. * delta_ul(ti::J, ti::i) / square(r));
    const auto d_imetric_00 =
        tenex::evaluate<ti::i>(2. * pos_lower(ti::i) / cube(r));

    /*const auto d2_metric_ij = tenex::evaluate<ti::i, ti::j, ti::K, ti::L>(
        -30. * pos_lower(ti::i) * pos_lower(ti::j) * pos(ti::K) * pos(ti::L) /
            (cube(r) * square(square(r))) +
        6. *
            (delta_ll(ti::i, ti::j) * pos(ti::K) * pos(ti::L) +
             delta_ul(ti::K, ti::i) * pos_lower(ti::j) * pos(ti::L) +
             delta_ul(ti::K, ti::j) * pos_lower(ti::i) * pos(ti::L) +
             delta_ul(ti::L, ti::i) * pos_lower(ti::j) * pos(ti::K) +
             delta_ul(ti::L, ti::j) * pos_lower(ti::i) * pos(ti::K)) /
            (square(r) * cube(r)) -
        2. *
            (delta_ul(ti::L, ti::i) * delta_ul(ti::K, ti::j) +
             delta_ul(ti::K, ti::i) * delta_ul(ti::L, ti::j)) /
            cube(r));

    const auto d2_metric_i0 = tenex::evaluate<ti::j, ti::k, ti::I>(
        16. * pos(ti::I) * pos_lower(ti::j) * pos_lower(ti::k) /
            (square(cube(r))) -
        4. *
            (delta_ll(ti::k, ti::j) * pos(ti::i) +
             delta_ul(ti::I, ti::k) * pos(ti::j) +
             delta_ul(ti::I, ti::j) * pos(ti::k)) /
            square(square(r)));
    const auto d2_metric_00 = tenex::evaluate<ti::i, ti::k>(
        2. * delta_ll(ij) / cube(r) -
        6. * pos_lower(ti::i) * pos_lower(ti::j) / (cube(r) * square(r)));*/
    for (size_t i = 0; i < Dim; ++i) {
      d_imetric.get(i, 0, 0) = d_imetric_00.get(i);
      d_metric.get(i, 0, 0) = -d_imetric_00.get(i);
      get<0>(dt_d_psiR) += particle_velocity.get(i) * dt_psi_dipole.get(i);
      for (size_t j = 0; j < Dim; ++j) {
        d_imetric.get(i, j + 1, 0) = d_imetric_i0.get(i, j);
        d_metric.get(i, j + 1, 0) = d_imetric_i0.get(i, j);
        // d2_metric.get(i, j, 0, 0) = d2_metric_00.get(i, j);
        for (size_t k = 0; k < Dim; ++k) {
          d_imetric.get(i, j + 1, k + 1) = d_imetric_ij.get(i, j, k);
          d_metric.get(i, j + 1, k + 1) = -d_imetric_ij.get(i, j, k);
          // d2_metric.get(i, j, k + 1, 0) = d2_metric_i0.get(i, j, k);
          for (size_t l = 0; l < Dim; ++l) {
            // d2_metric.get(i, j, k + 1, l + 1) = d2_metric_ij.get(i, j, k, l);
          }
        }
      }
    }

    const auto dt_four_velocity = tenex::evaluate<ti::A>(
        roll_on * charge / mass / u0 *
            inverse_spacetime_metric_inertial(ti::A, ti::B) * d_psiR(ti::b) -
        christoffel(ti::A, ti::b, ti::c) * four_velocity(ti::B) *
            four_velocity(ti::C) / u0);

    const auto dt_metric = tenex::evaluate<ti::A, ti::B>(
        particle_velocity(ti::I) * d_imetric(ti::i, ti::A, ti::B));
    /*const auto dt2_metric = tenex::evaluate<ti::A, ti::B>(
        particle_velocity(ti::I) * particle_velocity(ti::J) *
            d2_metric(ti::i, ti::j, ti::A, ti::B) +
        particle_acceleration(ti::I) * d_imetric(ti::i, ti::A, ti::B));*/
    const auto f = tenex::evaluate<ti::B>(
        roll_on * charge / mass * d_psiR(ti::a) *
        (inverse_spacetime_metric_inertial(ti::A, ti::B) +
         four_velocity(ti::A) * four_velocity(ti::B)));
    const auto dt_f = tenex::evaluate<ti::A>(
        roll_on * charge / mass *
        ((dt_metric(ti::A, ti::B) +
          four_velocity(ti::A) * dt_four_velocity(ti::B) +
          dt_four_velocity(ti::A) * four_velocity(ti::B)) *
             d_psiR(ti::b) +
         (inverse_spacetime_metric_inertial(ti::A, ti::B) +
          four_velocity(ti::A) * four_velocity(ti::B)) *
             dt_d_psiR(ti::b)));

    for (size_t i = 0; i < Dim; ++i) {
      get(*self_force)[i + 3] = f.get(i);
      get(*self_force)[i + 6] = dt_f.get(i);
    }
  }
  for (size_t i = 0; i < Dim; ++i) {
    get(*self_force)[i] = particle_acceleration.get(i);
  }

  for (size_t i = 0; i < Dim; ++i) {
    get<::Tags::dt<Tags::Position>>(*dt_evolved_vars).get(i)[0] =
        get<Tags::Velocity>(evolved_vars).get(i)[0];
    get<::Tags::dt<Tags::Velocity>>(*dt_evolved_vars).get(i)[0] =
        particle_acceleration.get(i);
  }
}

}  // namespace CurvedScalarWave::Worldtube
