// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/Context.hpp"
#include "Options/Options.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/PupStlCpp11.hpp"
#include "Utilities/TMPL.hpp"
#include <array>
#include <cmath>
#include <pup.h>
#include <utility>

namespace gh::BoundaryConditions::detail {
/// Stationary Schwarzschild metric in a frozen quadratic radial/height map.
/// r=r0+p x+A x^2/2, tau=a T+b x+B x^2/2, x=R-R0.
/// This is a local extension; it is not required to be DH away from the face.
struct SchwarzschildReferenceParameters {
  struct Mass {
    using type = double;
    static constexpr Options::String help{"Schwarzschild mass."};
  };
  struct CoordinateRadius {
    using type = double;
    static constexpr Options::String help{"Calibration coordinate radius."};
  };
  struct ArealRadius {
    using type = double;
    static constexpr Options::String help{"Areal radius at calibration."};
  };
  struct RadialSlope {
    using type = double;
    static constexpr Options::String help{"First derivative of areal radius."};
  };
  struct RadialCurvature {
    using type = double;
    static constexpr Options::String help{"Second derivative of areal radius."};
  };
  struct ClockRate {
    using type = double;
    static constexpr Options::String help{
        "Kerr-Schild time derivative with respect to inertial time."};
  };
  struct HeightSlope {
    using type = double;
    static constexpr Options::String help{
        "Radial derivative of Kerr-Schild time."};
  };
  struct HeightCurvature {
    using type = double;
    static constexpr Options::String help{
        "Second radial derivative of Kerr-Schild time."};
  };
  struct RelaxationTime {
    using type = double;
    static constexpr Options::String help{
        "Positive time scale for decay of the subtracted radiation residual."};
  };
  using options = tmpl::list<Mass, CoordinateRadius, ArealRadius, RadialSlope,
                             RadialCurvature, ClockRate, HeightSlope,
                             HeightCurvature, RelaxationTime>;
  static constexpr Options::String help{
      "Frozen local Schwarzschild reference. Static mesh only."};
  SchwarzschildReferenceParameters() = default;
  SchwarzschildReferenceParameters(double mass, double coordinate_radius,
                                   double areal_radius, double radial_slope,
                                   double radial_curvature, double clock_rate,
                                   double height_slope, double height_curvature,
                                   double relaxation_time,
                                   const Options::Context &context = {})
      : values{{mass, coordinate_radius, areal_radius, radial_slope,
                radial_curvature, clock_rate, height_slope, height_curvature,
                relaxation_time}} {
    for (const auto value : values) {
      if (not std::isfinite(value)) {
        PARSE_ERROR(context,
                    "SchwarzschildReference parameters must be finite.");
      }
    }
    if (mass <= 0. or coordinate_radius <= 0. or areal_radius <= 0. or
        radial_slope <= 0. or clock_rate <= 0. or relaxation_time <= 0.) {
      PARSE_ERROR(context, "SchwarzschildReference mass, radii, radial slope, "
                           "clock rate and relaxation time must be positive.");
    }
    const double z = 2. * mass / areal_radius;
    if ((1. + z) * radial_slope * radial_slope +
            2. * z * height_slope * radial_slope -
            (1. - z) * height_slope * height_slope <=
        0.) {
      PARSE_ERROR(
          context,
          "SchwarzschildReference slice must be spacelike at calibration.");
    }
  }
  std::array<double, 9> values{{1., 2.5, 2.5, 1., 0., 1., 0., 0., 20.}};
  void pup(PUP::er &p) { p | values; }
};
inline bool operator==(const SchwarzschildReferenceParameters &a,
                       const SchwarzschildReferenceParameters &b) {
  return a.values == b.values;
}
inline bool operator!=(const SchwarzschildReferenceParameters &a,
                       const SchwarzschildReferenceParameters &b) {
  return not(a == b);
}
struct SchwarzschildReferenceGauge {
  struct SchwarzschildReference {
    using type = SchwarzschildReferenceParameters;
    static constexpr Options::String help{
        "Prescribe stationary-reference-subtracted radiation data, relaxing "
        "its full live-frame residual to zero."};
  };
  using options = tmpl::list<SchwarzschildReference>;
  static constexpr Options::String help{
      "Fixed analytic Schwarzschild gauge reference."};
  SchwarzschildReferenceGauge() = default;
  explicit SchwarzschildReferenceGauge(SchwarzschildReferenceParameters p)
      : parameters(std::move(p)) {}
  SchwarzschildReferenceParameters parameters{};
};

/// Analytic Cartesian metric and spatial derivatives, including angular terms.
/// The center and map coefficients are held fixed in inertial coordinates.
inline auto schwarzschild_reference_fields(
    const SchwarzschildReferenceParameters &parameters,
    const tnsr::I<DataVector, 3> &coords, const tnsr::I<double, 3> &center) {
  const auto &v = parameters.values;
  const size_t size = get<0>(coords).size();
  DataVector R(size, 0.);
  for (size_t i = 0; i < 3; ++i) {
    R += square(coords.get(i) - center.get(i));
  }
  R = sqrt(R);
  const DataVector x = R - v[1];
  const DataVector r = v[2] + v[3] * x + 0.5 * v[4] * square(x);
  const DataVector p = v[3] + v[4] * x;
  const DataVector b = v[6] + v[7] * x;
  const DataVector z = 2. * v[0] / r;
  const DataVector dz = -z * p / r;
  const DataVector tt = -(1. - z) * square(v[5]);
  const DataVector tr = v[5] * (z * p - (1. - z) * b);
  const DataVector rr =
      (1. + z) * square(p) + 2. * z * b * p - (1. - z) * square(b);
  if (min(R) <= 0. or min(r) <= 0. or min(p) <= 0. or min(rr) <= 0.) {
    ERROR("SchwarzschildReference map is singular or its slice is not "
          "spacelike here.");
  }
  const DataVector d = square(r / R);
  const DataVector dtt = dz * square(v[5]);
  const DataVector dtr = v[5] * (dz * (p + b) + z * v[4] - (1. - z) * v[7]);
  const DataVector drr = dz * square(p + b) + 2. * (1. + z) * p * v[4] +
                         2. * z * (v[7] * p + b * v[4]) -
                         2. * (1. - z) * b * v[7];
  const DataVector dd = 2. * d * (p / r - 1. / R);
  tnsr::I<DataVector, 3> n(size);
  for (size_t i = 0; i < 3; ++i) {
    n.get(i) = (coords.get(i) - center.get(i)) / R;
  }
  tnsr::aa<DataVector, 3> g(size, 0.);
  tnsr::iaa<DataVector, 3> phi(size, 0.);
  get<0, 0>(g) = tt;
  for (size_t j = 0; j < 3; ++j) {
    g.get(0, j + 1) = tr * n.get(j);
    for (size_t k = j; k < 3; ++k) {
      g.get(j + 1, k + 1) =
          d * (j == k ? 1. : 0.) + (rr - d) * n.get(j) * n.get(k);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    phi.get(i, 0, 0) = dtt * n.get(i);
    for (size_t j = 0; j < 3; ++j) {
      const DataVector nij = (i == j ? 1. : 0.) - n.get(i) * n.get(j);
      phi.get(i, 0, j + 1) = dtr * n.get(i) * n.get(j) + tr / R * nij;
      for (size_t k = j; k < 3; ++k) {
        const DataVector nik = (i == k ? 1. : 0.) - n.get(i) * n.get(k);
        phi.get(i, j + 1, k + 1) =
            dd * n.get(i) * (j == k ? 1. : 0.) +
            (drr - dd) * n.get(i) * n.get(j) * n.get(k) +
            (rr - d) / R * (nij * n.get(k) + nik * n.get(j));
      }
    }
  }
  return std::make_pair(std::move(g), std::move(phi));
}
} // namespace gh::BoundaryConditions::detail
