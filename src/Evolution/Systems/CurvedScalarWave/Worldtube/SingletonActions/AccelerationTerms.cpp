// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/AccelerationTerms.hpp"

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Worldtube.hpp"
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

void AccelerationTermsMutator::apply(
    const gsl::not_null<Scalar<DataVector>*> self_force,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& position,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& velocity,
    const Scalar<double>& psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& psi_dipole,
    const Scalar<double>& dt_psi_monopole,
    const tnsr::i<double, Dim, Frame::Grid>& dt_psi_dipole, const double time,
    const double mass, const double charge, const double turn_on_time,
    const double turn_on_interval, const size_t iteration,
    const gr::Solutions::KerrSchild& kerr_schild) {
  tnsr::I<double, 3> pos{};
  tnsr::I<double, 3> vel{};
  for (size_t i = 0; i < 3; ++i) {
    pos.get(i) = position.get(i)[0];
    vel.get(i) = velocity.get(i)[0];
  }

  const auto [metric, imetric, di_metric, di_imetric, dij_metric, dij_imetric,
              christoffel, di_christoffel, contracted_christoffel,
              di_contracted_christoffel] =
      kerr_schild_quantities(kerr_schild, pos);

  tnsr::I<double, Dim> acc{};

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

    const double t_over_tsigma = t_minus_turnup / turn_on_interval;
    const double t_over_tsigma_pow4 = square(square(t_over_tsigma));
    const double roll_on = 1. - exp(-t_over_tsigma_pow4);
    const double dt_roll_on = 4. * cube(t_minus_turnup) *
                              exp(-t_over_tsigma_pow4) /
                              square(square(turn_on_interval));
    const double dt2_roll_on = 4. * square(t_minus_turnup) *
                               exp(-t_over_tsigma_pow4) *
                               (3. * square(square(turn_on_interval)) -
                                4. * square(square(t_minus_turnup))) /
                               square(square(square(turn_on_interval)));
    const double evolved_mass =
        iteration == 0 ? mass : mass - charge * get(psi_monopole);

    for (size_t i = 0; i < Dim; ++i) {
      acc.get(i) += (imetric.get(i + 1, 0) - vel.get(i) * imetric.get(0, 0)) *
                    get(dt_psi_monopole) * roll_on * charge / evolved_mass /
                    u0_squared;
      for (size_t j = 0; j < Dim; ++j) {
        acc.get(i) +=
            (imetric.get(i + 1, j + 1) - vel.get(i) * imetric.get(0, j + 1)) *
            psi_dipole.get(j) * roll_on * charge / evolved_mass / u0_squared;
      }
    }
    const double u0 = sqrt(u0_squared);
    const tnsr::A<double, Dim> u{
        {u0, get<0>(vel) * u0, get<1>(vel) * u0, get<2>(vel) * u0}};
    const tnsr::a<double, Dim> d_psiR{{get(dt_psi_monopole), get<0>(psi_dipole),
                                       get<1>(psi_dipole), get<2>(psi_dipole)}};
    double dt2_psiR = get<0>(contracted_christoffel) * get<0>(d_psiR);
    for (size_t i = 0; i < Dim; ++i) {
      dt2_psiR += -2. * imetric.get(0, i + 1) * dt_psi_dipole.get(i) +
                  contracted_christoffel.get(i + 1) * psi_dipole.get(i);
    }

    dt2_psiR /= get<0, 0>(imetric);

    tnsr::a<double, Dim> dt_d_psiR{{dt2_psiR, get<0>(dt_psi_dipole),
                                    get<1>(dt_psi_dipole),
                                    get<2>(dt_psi_dipole)}};

    for (size_t i = 0; i < Dim; ++i) {
      get<0>(dt_d_psiR) += vel.get(i) * dt_psi_dipole.get(i);
    }

    double dt_evolved_mass = get(dt_psi_monopole);
    double dt2_evolved_mass = dt2_psiR;

    for (size_t i = 0; i < Dim; ++i) {
      dt_evolved_mass += psi_dipole.get(i) * vel.get(i);
      dt2_evolved_mass += 2. * dt_psi_dipole.get(i) * vel.get(i) +
                          psi_dipole.get(i) * acc.get(i);
    }
    dt_evolved_mass *= iteration == 0 ? 0. : -charge;
    dt2_evolved_mass *= iteration == 0 ? 0. : -charge;

    const double dt_mass_factor =
        -dt_evolved_mass * charge / (evolved_mass * evolved_mass);
    const double dt2_mass_factor = charge *
                                   (2. * dt_evolved_mass * dt_evolved_mass -
                                    evolved_mass * dt2_evolved_mass) /
                                   cube(evolved_mass);

    const auto dt_christoffel = tenex::evaluate<ti::A, ti::b, ti::c>(
        vel(ti::I) * di_christoffel(ti::i, ti::A, ti::b, ti::c));

    const auto dt_imetric = tenex::evaluate<ti::A, ti::B>(
        vel(ti::I) * di_imetric(ti::i, ti::A, ti::B));
    const auto dt_metric = tenex::evaluate<ti::a, ti::b>(
        vel(ti::I) * di_metric(ti::i, ti::a, ti::b));
    const auto dt2_imetric = tenex::evaluate<ti::A, ti::B>(
        vel(ti::I) * vel(ti::J) * dij_imetric(ti::i, ti::j, ti::A, ti::B) +
        acc(ti::I) * di_imetric(ti::i, ti::A, ti::B));

    const auto dt_u = tenex::evaluate<ti::A>(
        charge / evolved_mass / u0 * imetric(ti::A, ti::B) * d_psiR(ti::b) -
        christoffel(ti::A, ti::b, ti::c) * u(ti::B) * u(ti::C) / u0);
    const auto dt2_u = tenex::evaluate<ti::A>(
        charge / evolved_mass / u0 *
            (dt_imetric(ti::A, ti::B) * d_psiR(ti::b) +
             imetric(ti::A, ti::B) * dt_d_psiR(ti::b)) +
        dt_mass_factor / u0 * imetric(ti::A, ti::B) * d_psiR(ti::b) -
        (dt_christoffel(ti::A, ti::b, ti::c) * u(ti::B) * u(ti::C) +
         2. * christoffel(ti::A, ti::b, ti::c) * dt_u(ti::B) * u(ti::C) +
         get<0>(dt_u) * dt_u(ti::A)) /
            u0);
    const auto dt_u_rollon = tenex::evaluate<ti::A>(roll_on * dt_u(ti::A));
    const auto dt2_u_rollon = tenex::evaluate<ti::A>(dt_roll_on * dt_u(ti::A) +
                                                     roll_on * dt2_u(ti::A));
    tnsr::i<double, Dim> d_dt2_psiR{0.};
    for (size_t i = 0; i < Dim; ++i) {
      d_dt2_psiR.get(i) +=
          -di_imetric.get(i, 0, 0) * dt2_psiR +
          di_contracted_christoffel.get(i, 0) * get(dt_psi_monopole) +
          contracted_christoffel.get(0) * dt_psi_dipole.get(i);
      for (size_t j = 0; j < Dim; ++j) {
        d_dt2_psiR.get(i) +=
            -2. * di_imetric.get(i, 0, j + 1) * dt_psi_dipole.get(j) +
            di_contracted_christoffel.get(i, j + 1) * psi_dipole.get(j);
      }
      d_dt2_psiR.get(i) /= get<0, 0>(imetric);
    }
    double dt3_psiR = contracted_christoffel.get(0) * dt2_psiR;
    for (size_t i = 0; i < Dim; ++i) {
      dt3_psiR += -2. * imetric.get(0, i + 1) * d_dt2_psiR.get(i) +
                  contracted_christoffel.get(i + 1) * dt_psi_dipole.get(i);
    }

    tnsr::a<double, Dim> dt2_d_psiR;
    dt2_d_psiR.get(0) = dt3_psiR;
    for (size_t i = 0; i < Dim; ++i) {
      dt2_d_psiR.get(0) += 2. * d_dt2_psiR.get(i) * vel.get(i) +
                           dt_psi_dipole.get(i) * acc.get(i);
      dt2_d_psiR.get(i) = d_dt2_psiR.get(i);
    }

    const auto f =
        tenex::evaluate<ti::B>(charge / evolved_mass * d_psiR(ti::a) *
                               (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B)));
    auto dt_f = tenex::evaluate<ti::A>(
        charge / evolved_mass *
        ((dt_imetric(ti::A, ti::B) + u(ti::A) * dt_u_rollon(ti::B) +
          dt_u_rollon(ti::A) * u(ti::B)) *
             d_psiR(ti::b) +
         (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B)) * dt_d_psiR(ti::b)));
    auto dt2_f = tenex::evaluate<ti::A>(
        charge / evolved_mass *
        ((dt2_imetric(ti::A, ti::B) + dt2_u_rollon(ti::A) * u(ti::B) +
          dt2_u_rollon(ti::B) * u(ti::A) +
          2. * dt_u_rollon(ti::A) * dt_u_rollon(ti::B)) *
             d_psiR(ti::b) +
         2. * dt_d_psiR(ti::b) *
             (dt_imetric(ti::A, ti::B) + dt_u_rollon(ti::A) * u(ti::B) +
              dt_u_rollon(ti::B) * u(ti::A)) +
         dt2_d_psiR(ti::b) * (imetric(ti::A, ti::B) + u(ti::A) * u(ti::B))));
    for (size_t i = 0; i < 4; ++i) {
      dt_f.get(i) += dt_mass_factor * f.get(i);
      dt2_f.get(i) +=
          2. * dt_mass_factor * dt_f.get(i) + dt2_mass_factor * f.get(i);
    }
    const auto f_roll_on = tenex::evaluate<ti::A>(roll_on * f(ti::A));
    const auto dt_f_roll_on =
        tenex::evaluate<ti::A>(dt_roll_on * f(ti::A) + roll_on * dt_f(ti::A));
    const auto dt2_f_roll_on = tenex::evaluate<ti::A>(
        dt2_roll_on * f(ti::A) + 2. * dt_roll_on * dt_f(ti::A) +
        roll_on * dt2_f(ti::A));

    const auto cov_f = tenex::evaluate<ti::A>(dt_f_roll_on(ti::A) * u0 +
                                              christoffel(ti::A, ti::b, ti::c) *
                                                  u(ti::B) * f_roll_on(ti::C));

    const auto dt_cov_f = tenex::evaluate<ti::A>(
        dt2_f_roll_on(ti::A) * u0 + dt_f_roll_on(ti::A) * get<0>(dt_u_rollon) +
        dt_christoffel(ti::A, ti::b, ti::c) * u(ti::B) * f_roll_on(ti::C) +
        christoffel(ti::A, ti::b, ti::c) * dt_u_rollon(ti::B) *
            f_roll_on(ti::C) +
        christoffel(ti::A, ti::b, ti::c) * u(ti::B) * dt_f_roll_on(ti::C));

    for (size_t i = 0; i < Dim; ++i) {
      get(*self_force)[i + 3] = f_roll_on.get(i);
      get(*self_force)[i + 6] = dt_f_roll_on.get(i);
      get(*self_force)[i + 9] = cov_f.get(i);
      get(*self_force)[i + 12] = dt_cov_f.get(i);
    }
  }
  for (size_t i = 0; i < Dim; ++i) {
    get(*self_force)[i] = acc.get(i);
  }
}

}  // namespace CurvedScalarWave::Worldtube
