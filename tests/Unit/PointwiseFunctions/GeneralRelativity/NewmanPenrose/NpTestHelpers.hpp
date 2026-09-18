// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <random>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/Pypp.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TypeD.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Types.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace TestHelpers::gr::np {
using frame = Frame::Inertial;

/// A random face: perturbed flat metric, random unit sphere directions, the
/// adapted triad, and lapse and shift near the static values.
struct RandomGeometry {
  tnsr::ii<DataVector, 3, frame> spatial_metric;
  tnsr::I<DataVector, 3, frame> directions;
  ::gr::np::RealMatrix rotation;
  Scalar<DataVector> lapse;
  tnsr::I<DataVector, 3, frame> shift;
};

template <typename Generator>
RandomGeometry random_geometry(const gsl::not_null<Generator*> generator,
                               const size_t num_points) {
  std::uniform_real_distribution<> small(-0.1, 0.1);
  std::uniform_real_distribution<> unit(-1.0, 1.0);
  const DataVector used_for_size(num_points);
  RandomGeometry geometry{};
  geometry.spatial_metric =
      make_with_random_values<tnsr::ii<DataVector, 3, frame>>(
          generator, make_not_null(&small), used_for_size);
  for (size_t i = 0; i < 3; ++i) {
    geometry.spatial_metric.get(i, i) += 1.;
  }
  geometry.directions = make_with_random_values<tnsr::I<DataVector, 3, frame>>(
      generator, make_not_null(&unit), used_for_size);
  const DataVector norm = sqrt(square(get<0>(geometry.directions)) +
                               square(get<1>(geometry.directions)) +
                               square(get<2>(geometry.directions)));
  for (size_t i = 0; i < 3; ++i) {
    geometry.directions.get(i) /= norm;
  }
  geometry.rotation =
      ::gr::np::adapted_triad(geometry.spatial_metric, geometry.directions);
  geometry.lapse = make_with_value<Scalar<DataVector>>(used_for_size, 1.);
  get(geometry.lapse) += make_with_random_values<DataVector>(
      generator, make_not_null(&small), used_for_size);
  geometry.shift = make_with_random_values<tnsr::I<DataVector, 3, frame>>(
      generator, make_not_null(&small), used_for_size);
  return geometry;
}

template <typename Generator>
ComplexDataVector random_complex(const gsl::not_null<Generator*> generator,
                                 const size_t num_points, const double scale) {
  std::uniform_real_distribution<> dist(-scale, scale);
  const DataVector used_for_size(num_points);
  const auto re = make_with_random_values<DataVector>(
      generator, make_not_null(&dist), used_for_size);
  const auto im = make_with_random_values<DataVector>(
      generator, make_not_null(&dist), used_for_size);
  ComplexDataVector result(num_points);
  for (size_t p = 0; p < num_points; ++p) {
    result[p] = std::complex<double>{re[p], im[p]};
  }
  return result;
}

/// Arbitrary complex Weyl scalars
template <typename Generator>
::gr::np::WeylScalars random_scalars(const gsl::not_null<Generator*> generator,
                                     const size_t num_points,
                                     const double scale) {
  ::gr::np::WeylScalars psi(num_points, std::complex<double>{0., 0.});
  for (size_t a = 0; a < 5; ++a) {
    psi.get(a) = random_complex(generator, num_points, scale);
  }
  return psi;
}

/// Type-D scalars in a slightly misaligned tetrad plus a small tide-like
/// perturbation: the regime the type-D solve is made for.
struct NearTypeD {
  Scalar<ComplexDataVector> coulomb;
  Scalar<ComplexDataVector> a_bar;
  Scalar<ComplexDataVector> b;
  ::gr::np::WeylScalars psi;
};

template <typename Generator>
NearTypeD near_type_d(const gsl::not_null<Generator*> generator,
                      const size_t num_points, const double misalignment,
                      const double perturbation) {
  NearTypeD result{};
  std::uniform_real_distribution<> unit(0.5, 1.5);
  const DataVector used_for_size(num_points);
  get(result.coulomb) = ComplexDataVector(num_points);
  const auto magnitude = make_with_random_values<DataVector>(
      generator, make_not_null(&unit), used_for_size);
  const ComplexDataVector phase =
      random_complex(generator, num_points, 0.05 * misalignment);
  for (size_t p = 0; p < num_points; ++p) {
    get(result.coulomb)[p] = -magnitude[p] * (1. + phase[p]);
  }
  get(result.a_bar) = random_complex(generator, num_points, misalignment);
  get(result.b) = random_complex(generator, num_points, misalignment);
  result.psi = ::gr::np::type_d_scalars(result.coulomb, result.a_bar, result.b);
  const auto tide = random_scalars(generator, num_points, perturbation);
  for (size_t a = 0; a < 5; ++a) {
    result.psi.get(a) += tide.get(a);
  }
  return result;
}

/// The face data as the lists the sphere-wide Python adapters expect
inline std::vector<DataVector> pack_reals(const RandomGeometry& geometry) {
  std::vector<DataVector> reals{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      reals.push_back(geometry.spatial_metric.get(i, j));
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      reals.push_back(geometry.rotation.get(i, j));
    }
  }
  reals.push_back(get(geometry.lapse));
  for (size_t i = 0; i < 3; ++i) {
    reals.push_back(geometry.shift.get(i));
  }
  for (size_t i = 0; i < 3; ++i) {
    reals.push_back(geometry.directions.get(i));
  }
  return reals;
}

inline std::vector<ComplexDataVector> pack_scalars(
    const ::gr::np::WeylScalars& psi) {
  std::vector<ComplexDataVector> result{};
  for (size_t a = 0; a < 5; ++a) {
    result.push_back(psi.get(a));
  }
  return result;
}

/// The manufactured slice of the minimal study: exact type-D background with
/// a generic quadrupole tide, boosted at speed 0.22, on the flat unit-lapse
/// slice. `truth` holds the ten real moment components
/// (E11, E22, E12, E13, E23, B11, ...).
struct ManufacturedSlice {
  RandomGeometry geometry;
  ::gr::np::WeylScalars psi;
  std::array<double, 10> truth;
  double mass;
};

inline ManufacturedSlice manufactured_slice(const int seed,
                                            const size_t num_points) {
  ManufacturedSlice slice{};
  const auto psi_list = pypp::call<std::vector<ComplexDataVector>>(
      "NpMatching", "manufactured_psi", seed, num_points);
  const auto reals = pypp::call<std::vector<DataVector>>(
      "NpMatching", "manufactured_reals", seed, num_points);
  slice.truth = pypp::call<std::array<double, 10>>(
      "NpMatching", "manufactured_truth", seed, num_points);
  slice.mass = pypp::call<double>("NpMatching", "manufactured_mass");
  slice.psi = ::gr::np::WeylScalars(num_points, std::complex<double>{0., 0.});
  for (size_t a = 0; a < 5; ++a) {
    slice.psi.get(a) = psi_list[a];
  }
  auto& geometry = slice.geometry;
  geometry.spatial_metric = tnsr::ii<DataVector, 3, frame>(num_points, 0.);
  geometry.rotation = ::gr::np::RealMatrix(num_points, 0.);
  geometry.shift = tnsr::I<DataVector, 3, frame>(num_points, 0.);
  geometry.directions = tnsr::I<DataVector, 3, frame>(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      if (j >= i) {
        geometry.spatial_metric.get(i, j) = reals[3 * i + j];
      }
      geometry.rotation.get(i, j) = reals[9 + 3 * i + j];
    }
  }
  get(geometry.lapse) = reals[18];
  for (size_t i = 0; i < 3; ++i) {
    geometry.shift.get(i) = reals[19 + i];
    geometry.directions.get(i) = reals[22 + i];
  }
  return slice;
}
}  // namespace TestHelpers::gr::np
