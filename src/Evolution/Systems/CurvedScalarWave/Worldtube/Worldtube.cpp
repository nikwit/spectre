// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/Worldtube.hpp"

#include <tuple>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"

namespace CurvedScalarWave::Worldtube {

double worldtube_shrink_factor(const double orbit_radius,
                               const double original_orbit_radius) {
  const double orbit_radius_fraction = orbit_radius / original_orbit_radius;
  return orbit_radius_fraction * sqrt(orbit_radius_fraction);
}

double worldtube_shrink_factor_derivative(const double orbit_radius,
                                          const double orbit_velocity,
                                          const double original_orbit_radius) {
  const double orbit_radius_fraction = orbit_radius / original_orbit_radius;
  return 1.5 * sqrt(orbit_radius_fraction) * orbit_velocity /
         original_orbit_radius;
}

double broken_power_shrink(const double orbit_radius, const double amp,
                           const double rb) {
  const double delta = 0.05;
  const double r_by_rb = orbit_radius / rb;
  return amp * r_by_rb * sqrt(r_by_rb) *
         pow(0.5 * (1. + pow(r_by_rb, 1. / delta)), -1.5 * delta);
}
double broken_power_shrink_derivative(const double orbit_radius,
                                      const double amp, const double rb,
                                      const double orbit_radius_derivative) {
  const double delta = 0.05;
  const double r_by_rb = orbit_radius / rb;
  double sol = amp * 1.5 * sqrt(r_by_rb) * orbit_radius_derivative / rb *
               pow(0.5 * (1. + pow(r_by_rb, 1. / delta)), -1.5 * delta);
  sol += amp * r_by_rb * sqrt(r_by_rb) * -1.5 * delta *
         pow(0.5 * (1. + pow(r_by_rb, 1. / delta)), -1.5 * delta - 1.) * 0.5 /
         delta * pow(r_by_rb, 1. / delta - 1.) * orbit_radius_derivative / rb;
  return sol;
}

template <size_t Dim>
std::tuple<tnsr::aa<double, Dim>, tnsr::AA<double, Dim>, tnsr::iaa<double, Dim>,
           tnsr::iAA<double, Dim>, tnsr::ijaa<double, Dim>,
           tnsr::ijAA<double, Dim>, tnsr::Abb<double, Dim>,
           tnsr::iAbb<double, Dim>, tnsr::A<double, Dim>, tnsr::iA<double, Dim>>
kerr_schild_quantities(const gr::Solutions::KerrSchild& kerr_schild,
                       const tnsr::I<double, Dim>& pos) {
  const auto spacetime_vars = kerr_schild.variables(
      pos, 0.,
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

  const double r = magnitude(pos).get();

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
      6. * pos(ti::J) * pos(ti::K) * pos_lower(ti::i) / (square(r) * cube(r)) -
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

  tenex::evaluate<ti::i, ti::a, ti::b>(
      make_not_null(&di_metric), -metric(ti::a, ti::c) * metric(ti::b, ti::d) *
                                     di_imetric(ti::i, ti::C, ti::D));
  tenex::evaluate<ti::j, ti::i, ti::a, ti::b>(
      make_not_null(&dij_metric),
      -metric(ti::a, ti::c) * metric(ti::b, ti::d) *
              dij_imetric(ti::j, ti::i, ti::C, ti::D) -
          2. * metric(ti::a, ti::c) * di_metric(ti::j, ti::b, ti::d) *
              di_imetric(ti::i, ti::C, ti::D));

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
  const auto christoffel =
      gr::christoffel_second_kind(d_spacetime_metric, imetric);
  tenex::evaluate<ti::i, ti::A, ti::b, ti::c>(
      make_not_null(&di_christoffel),
      0.5 * di_imetric(ti::i, ti::A, ti::D) *
              (d_metric(ti::b, ti::c, ti::d) + d_metric(ti::c, ti::b, ti::d) -
               d_metric(ti::d, ti::b, ti::c)) +
          0.5 * imetric(ti::A, ti::D) *
              (di_d_metric(ti::i, ti::b, ti::c, ti::d) +
               di_d_metric(ti::i, ti::c, ti::b, ti::d) -
               di_d_metric(ti::i, ti::d, ti::b, ti::c)));

  const auto contracted_christoffel = trace_last_indices(christoffel, imetric);
  tnsr::iA<double, Dim> di_contracted_christoffel;
  for (size_t i = 0; i < Dim; ++i) {
    di_contracted_christoffel.get(i, 0) = 4. * pos.get(i) / square(square(r));
    for (size_t j = 0; j < Dim; ++j) {
      di_contracted_christoffel.get(i, j + 1) =
          -6. * pos.get(i) * pos.get(j) / (square(r) * cube(r));
    }
    di_contracted_christoffel.get(i, i + 1) += 2. / cube(r);
  }

  return std::make_tuple(metric, imetric, di_metric, di_imetric, dij_metric,
                         dij_imetric, christoffel, di_christoffel,
                         contracted_christoffel, di_contracted_christoffel);
}

template std::tuple<tnsr::aa<double, 3>, tnsr::AA<double, 3>,
                    tnsr::iaa<double, 3>, tnsr::iAA<double, 3>,
                    tnsr::ijaa<double, 3>, tnsr::ijAA<double, 3>,
                    tnsr::Abb<double, 3>, tnsr::iAbb<double, 3>,
                    tnsr::A<double, 3>, tnsr::iA<double, 3>>
kerr_schild_quantities<3>(const gr::Solutions::KerrSchild&,
                          const tnsr::I<double, 3>&);

}  // namespace CurvedScalarWave::Worldtube
