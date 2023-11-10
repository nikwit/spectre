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
  const auto imetric = gr::inverse_spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars));
  const auto metric = gr::spacetime_metric(
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

  const auto christoffel =
      gr::christoffel_second_kind(d_spacetime_metric, imetric);

  tnsr::I<double, Dim> acc{};
  const auto& pos = inertial_particle_position;
  const auto& vel = particle_velocity;

  double u0_squared = metric.get(0, 0);
  for (size_t i = 0; i < Dim; ++i) {
    acc.get(i) =
        vel.get(i) * christoffel.get(0, 0, 0) - christoffel.get(i + 1, 0, 0);
    u0_squared += 2. * metric.get(i + 1, 0) * vel.get(i);

    for (size_t j = 0; j < Dim; ++j) {
      u0_squared += metric.get(i + 1, j + 1) * vel.get(i) * vel.get(j);
      acc.get(i) += 2. * vel.get(j) *
                    (vel.get(i) * christoffel.get(0, j + 1, 0) -
                     christoffel.get(i + 1, j + 1, 0));
      for (size_t k = 0; k < Dim; ++k) {
        acc.get(i) += vel.get(j) * vel.get(k) *
                      (vel.get(i) * christoffel.get(0, j + 1, k + 1) -
                       christoffel.get(i + 1, j + 1, k + 1));
      }
    }
  }
  get(*self_force) = DataVector(15, 0.);
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
      acc.get(i) += (imetric.get(i + 1, 0) - vel.get(i) * imetric.get(0, 0)) *
                    get(dt_psi_monopole) * roll_on * charge / mass / u0_squared;
      if (expansion_order > 0) {
        for (size_t j = 0; j < Dim; ++j) {
          acc.get(i) +=
              (imetric.get(i + 1, j + 1) - vel.get(i) * imetric.get(0, j + 1)) *
              psi_dipole.get(j) * roll_on * charge / mass / u0_squared;
        }
      }
    }
    const double u0 = sqrt(u0_squared);
    const tnsr::A<double, Dim> u{
        {u0, get<0>(vel) * u0, get<1>(vel) * u0, get<2>(vel) * u0}};
    const tnsr::a<double, Dim> d_psiR{{get(dt_psi_monopole), get<0>(psi_dipole),
                                       get<1>(psi_dipole), get<2>(psi_dipole)}};
    const auto contracted_christoffel_inertial =
        trace_last_indices(christoffel, imetric);
    double dt2_psiR = get<0>(contracted_christoffel_inertial) * get<0>(d_psiR);
    for (size_t i = 0; i < Dim; ++i) {
      dt2_psiR +=
          -2. * imetric.get(0, i + 1) * dt_psi_dipole.get(i) +
          contracted_christoffel_inertial.get(i + 1) * d_psiR.get(i + 1);
    }

    dt2_psiR /= get<0, 0>(imetric);

    const double r = magnitude(pos).get();
    tnsr::a<double, Dim> dt_d_psiR{{dt2_psiR, get<0>(dt_psi_dipole),
                                    get<1>(dt_psi_dipole),
                                    get<2>(dt_psi_dipole)}};

    tnsr::iAA<double, Dim> di_imetric{};
    tnsr::ijAA<double, Dim> dij_imetric{};
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

    const auto d2_imetric_ij = tenex::evaluate<ti::i, ti::j, ti::K, ti::L>(
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

    const auto d2_imetric_i0 = tenex::evaluate<ti::j, ti::k, ti::I>(
        16. * pos(ti::I) * pos_lower(ti::j) * pos_lower(ti::k) /
            (square(cube(r))) -
        4. *
            (delta_ll(ti::k, ti::j) * pos(ti::I) +
             delta_ul(ti::I, ti::k) * pos_lower(ti::j) +
             delta_ul(ti::I, ti::j) * pos_lower(ti::k)) /
            square(square(r)));
    const auto d2_imetric_00 = tenex::evaluate<ti::i, ti::j>(
        2. * delta_ll(ti::i, ti::j) / cube(r) -
        6. * pos_lower(ti::i) * pos_lower(ti::j) / (cube(r) * square(r)));
    for (size_t i = 0; i < Dim; ++i) {
      di_imetric.get(i, 0, 0) = d_imetric_00.get(i);
      get<0>(dt_d_psiR) += vel.get(i) * dt_psi_dipole.get(i);
      for (size_t j = 0; j < Dim; ++j) {
        di_imetric.get(i, j + 1, 0) = d_imetric_i0.get(i, j);
        dij_imetric.get(i, j, 0, 0) = d2_imetric_00.get(i, j);
        for (size_t k = 0; k < Dim; ++k) {
          di_imetric.get(i, j + 1, k + 1) = d_imetric_ij.get(i, j, k);
          dij_imetric.get(i, j, k + 1, 0) = d2_imetric_i0.get(i, j, k);
          for (size_t l = 0; l < Dim; ++l) {
            dij_imetric.get(i, j, k + 1, l + 1) = d2_imetric_ij.get(i, j, k, l);
          }
        }
      }
    }

    tnsr::iaa<double, Dim> di_metric{};
    tnsr::ijaa<double, Dim> dij_metric{};
    tnsr::iAbb<double, Dim> di_christoffel{};
    tnsr::abb<double, Dim> d_metric{};
    tnsr::iabb<double, Dim> di_d_metric{};

    tenex::evaluate<ti::i, ti::a, ti::b>(make_not_null(&di_metric),
                                         -metric(ti::a, ti::c) *
                                             metric(ti::b, ti::d) *
                                             di_imetric(ti::i, ti::C, ti::D));
    tenex::evaluate<ti::i, ti::j, ti::a, ti::b>(
        make_not_null(&dij_metric),
        -metric(ti::a, ti::c) * metric(ti::b, ti::d) *
                dij_imetric(ti::i, ti::j, ti::C, ti::D) -
            2. * metric(ti::a, ti::c) * d_metric(ti::i, ti::b, ti::d) *
                di_imetric(ti::j, ti::C, ti::D));

    for (size_t a = 0; a <= Dim; ++a) {
      for (size_t b = 0; b <= Dim; ++b) {
        d_metric.get(0, a, b) = 0.;
        for (size_t i = 0; i < Dim; ++i) {
          d_metric.get(i + 1, a, b) = di_metric.get(i, a, b);
          di_d_metric.get(i, 0, a, b) = 0.;
          for (size_t j = 0; j < Dim; ++j) {
            di_d_metric.get(i, j + 1, a, b) = dij_metric.get(i, j, a, b);
          }
        }
      }
    }

    tenex::evaluate<ti::i, ti::A, ti::b, ti::c>(
        make_not_null(&di_christoffel),
        0.5 * di_imetric(ti::i, ti::A, ti::D) *
                (d_metric(ti::b, ti::c, ti::d) + d_metric(ti::c, ti::b, ti::d) -
                 d_metric(ti::d, ti::b, ti::c)) +
            0.5 * imetric(ti::A, ti::D) *
                (di_d_metric(ti::i, ti::b, ti::c, ti::d) +
                 di_d_metric(ti::i, ti::c, ti::b, ti::d) -
                 di_d_metric(ti::i, ti::d, ti::b, ti::c)));

    const auto dt_christoffel = tenex::evaluate<ti::A, ti::b, ti::c>(
        vel(ti::I) * di_christoffel(ti::i, ti::A, ti::b, ti::c));

    const auto dt_u = tenex::evaluate<ti::A>(
        roll_on * charge / mass / u0 * imetric(ti::A, ti::B) * d_psiR(ti::b) -
        christoffel(ti::A, ti::b, ti::c) * u(ti::B) * u(ti::C) / u0);

    const auto dt_imetric = tenex::evaluate<ti::A, ti::B>(
        vel(ti::I) * di_imetric(ti::i, ti::A, ti::B));
    const auto dt_metric = tenex::evaluate<ti::a, ti::b>(
        vel(ti::I) * di_metric(ti::i, ti::a, ti::b));
    const auto dt2_imetric = tenex::evaluate<ti::A, ti::B>(
        vel(ti::I) * vel(ti::J) * dij_imetric(ti::i, ti::j, ti::A, ti::B) +
        acc(ti::I) * di_imetric(ti::i, ti::A, ti::B));

    const auto dt2_u = tenex::evaluate<ti::A>(
        roll_on * charge / mass / u0 *
            (dt_imetric(ti::A, ti::B) * d_psiR(ti::b) +
             imetric(ti::A, ti::B) * dt_d_psiR(ti::b)) -
        (dt_christoffel(ti::A, ti::b, ti::c) * u(ti::B) * u(ti::C) +
         2. * christoffel(ti::A, ti::b, ti::c) * dt_u(ti::B) * u(ti::C) +
         get<0>(dt_u) * dt_u(ti::A)) /
            u0);
    tnsr::iA<double, Dim> d_contracted_christoffel;
    for (size_t i = 0; i < Dim; ++i) {
      d_contracted_christoffel.get(i, 0) = 4. * pos.get(i) / square(square(r));
      for (size_t j = 0; j < Dim; ++j) {
        d_contracted_christoffel.get(i, j + 1) =
            -6. * pos.get(i) * pos.get(j) / (square(r) * cube(r));
      }
      d_contracted_christoffel.get(i, i + 1) += 2 / cube(r);
    }

    tnsr::i<double, Dim> d_dt2_psiR{0.};
    for (size_t i = 0; i < Dim; ++i) {
      d_dt2_psiR.get(i) +=
          -di_imetric.get(i, 0, 0) * dt2_psiR +
          d_contracted_christoffel.get(i, 0) * get(dt_psi_monopole) +
          contracted_christoffel_inertial.get(0) * dt_d_psiR.get(i + 1);
      for (size_t j = 0; j < Dim; ++j) {
        d_dt2_psiR.get(i) +=
            -2. * di_imetric.get(i, 0, j + 1) * dt_d_psiR.get(j + 1) +
            d_contracted_christoffel.get(i, j + 1) * psi_dipole.get(j);
      }
      d_dt2_psiR.get(i) /= get<0, 0>(imetric);
    }
    double dt3_psiR = contracted_christoffel_inertial.get(0) * dt2_psiR;
    for (size_t i = 0; i < Dim; ++i) {
      dt3_psiR +=
          -2. * imetric.get(0, i + 1) * d_dt2_psiR.get(i) +
          contracted_christoffel_inertial.get(i + 1) * dt_d_psiR.get(i + 1);
    }

    tnsr::a<double, Dim> dt2_d_psiR;
    dt2_d_psiR.get(0) = dt3_psiR;
    for (size_t i = 0; i < Dim; ++i) {
      dt2_d_psiR.get(0) +=
          2. * d_dt2_psiR.get(i) * vel.get(i) + dt_d_psiR.get(i) * acc.get(i);
      dt2_d_psiR.get(i) = d_dt2_psiR.get(i);
    }

    const auto f =
        tenex::evaluate<ti::B>(roll_on * charge / mass * d_psiR(ti::a) *
                               (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B)));
    const auto dt_f = tenex::evaluate<ti::A>(
        roll_on * charge / mass *
        ((dt_imetric(ti::A, ti::B) + u(ti::A) * dt_u(ti::B) +
          dt_u(ti::A) * u(ti::B)) *
             d_psiR(ti::b) +
         (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B)) * dt_d_psiR(ti::b)));
    const auto dt2_f = tenex::evaluate<ti::A>(
        (dt2_imetric(ti::A, ti::B) + dt2_u(ti::A) * u(ti::B) +
         dt2_u(ti::B) * u(ti::A) + 2. * dt_u(ti::A) * dt_u(ti::B)) *
            d_psiR(ti::b) +
        2. * dt_d_psiR(ti::b) *
            (dt_imetric(ti::A, ti::B) + dt_u(ti::A) * u(ti::B) +
             dt_u(ti::B) * u(ti::A)) +
        dt2_d_psiR(ti::b) * (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B)));

    const auto cov_f = tenex::evaluate<ti::A>(dt_f(ti::A) * u0 +
                                              christoffel(ti::A, ti::b, ti::c) *
                                                  u(ti::B) * f(ti::C));

    const auto dt_cov_f = tenex::evaluate<ti::A>(
        dt2_f(ti::A) * u0 + dt_f(ti::A) * get<0>(dt_u) +
        dt_christoffel(ti::A, ti::b, ti::c) * u(ti::B) * f(ti::C) +
        christoffel(ti::A, ti::b, ti::c) * dt_u(ti::B) * f(ti::C) +
        christoffel(ti::A, ti::b, ti::c) * u(ti::B) * dt_f(ti::C));

    for (size_t i = 0; i < Dim; ++i) {
      get(*self_force)[i + 3] = f.get(i);
      get(*self_force)[i + 6] = dt_f.get(i);
      get(*self_force)[i + 9] = cov_f.get(i);
      get(*self_force)[i + 12] = dt_cov_f.get(i);
    }
  }
  for (size_t i = 0; i < Dim; ++i) {
    get(*self_force)[i] = acc.get(i);
  }

  for (size_t i = 0; i < Dim; ++i) {
    get<::Tags::dt<Tags::Position>>(*dt_evolved_vars).get(i)[0] =
        get<Tags::Velocity>(evolved_vars).get(i)[0];
    get<::Tags::dt<Tags::Velocity>>(*dt_evolved_vars).get(i)[0] = acc.get(i);
  }
}

}  // namespace CurvedScalarWave::Worldtube
