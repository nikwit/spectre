// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/DuDtTempTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Harmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/TimeDerivative.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "Utilities/Gsl.hpp"

extern "C" void CkRegisterMainModule(void) {}

namespace {
constexpr size_t Dim = 3;

DataVector field(const DataVector& x, const DataVector& y, const DataVector& z,
                 const double offset, const double ax, const double ay,
                 const double az) {
  return offset + ax * x + ay * square(y) + az * sin(z + 0.25 * x);
}

void set_identity_inverse_jacobian(
    const gsl::not_null<InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                                        Frame::Inertial>*>
        inverse_jacobian,
    const size_t number_of_points) {
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      inverse_jacobian->get(i, j) =
          DataVector(number_of_points, i == j ? 1.0 : 0.0);
    }
  }
}

template <typename TensorType>
void fill_tensor(const gsl::not_null<TensorType*> tensor,
                 const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
                 const double scale) {
  for (size_t storage_index = 0; storage_index < tensor->size();
       ++storage_index) {
    (*tensor)[storage_index] =
        field(get<0>(coords), get<1>(coords), get<2>(coords),
              scale * (1.0 + 0.03 * static_cast<double>(storage_index)),
              0.01 * static_cast<double>(storage_index + 1),
              -0.008 * static_cast<double>(storage_index + 2),
              0.006 * static_cast<double>(storage_index + 3));
  }
}

struct BenchmarkState {
  explicit BenchmarkState(const size_t points_per_dimension)
      : mesh(points_per_dimension, Spectral::Basis::Legendre,
             Spectral::Quadrature::GaussLobatto),
        number_of_points(mesh.number_of_grid_points()),
        inertial_coords(number_of_points),
        spacetime_metric(number_of_points),
        pi(number_of_points),
        phi(number_of_points),
        d_spacetime_metric(number_of_points),
        d_pi(number_of_points),
        d_phi(number_of_points),
        gamma0(number_of_points),
        gamma1(number_of_points),
        gamma2(number_of_points),
        dt_spacetime_metric(number_of_points),
        dt_pi(number_of_points),
        dt_phi(number_of_points),
        buffer(number_of_points) {
    const auto logical_coords = logical_coordinates(mesh);
    for (size_t i = 0; i < Dim; ++i) {
      inertial_coords.get(i) = logical_coords.get(i);
    }
    set_identity_inverse_jacobian(make_not_null(&inverse_jacobian),
                                  number_of_points);

    tnsr::ii<DataVector, Dim> spatial_metric(number_of_points);
    get<0, 0>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), 1.30, 0.04, 0.02, 0.01);
    get<1, 1>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), 1.25, -0.03, 0.015, -0.02);
    get<2, 2>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), 1.20, 0.02, -0.01, 0.015);
    get<0, 1>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), 0.03, 0.01, 0.0, 0.004);
    get<0, 2>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), -0.02, 0.0, 0.006, -0.003);
    get<1, 2>(spatial_metric) =
        field(get<0>(inertial_coords), get<1>(inertial_coords),
              get<2>(inertial_coords), 0.015, -0.004, 0.003, 0.002);

    Scalar<DataVector> lapse(number_of_points);
    get(lapse) = field(get<0>(inertial_coords), get<1>(inertial_coords),
                       get<2>(inertial_coords), 1.1, 0.03, -0.02, 0.01);
    tnsr::I<DataVector, Dim> shift(number_of_points);
    get<0>(shift) = field(get<0>(inertial_coords), get<1>(inertial_coords),
                          get<2>(inertial_coords), 0.05, 0.015, 0.0, 0.01);
    get<1>(shift) = field(get<0>(inertial_coords), get<1>(inertial_coords),
                          get<2>(inertial_coords), -0.04, 0.0, -0.012, 0.005);
    get<2>(shift) = field(get<0>(inertial_coords), get<1>(inertial_coords),
                          get<2>(inertial_coords), 0.03, -0.01, 0.006, 0.0);
    gr::spacetime_metric(make_not_null(&spacetime_metric), lapse, shift,
                         spatial_metric);

    fill_tensor(make_not_null(&pi), inertial_coords, 0.08);
    fill_tensor(make_not_null(&phi), inertial_coords, 0.04);
    fill_tensor(make_not_null(&d_spacetime_metric), inertial_coords, 0.02);
    fill_tensor(make_not_null(&d_pi), inertial_coords, 0.03);
    fill_tensor(make_not_null(&d_phi), inertial_coords, 0.025);

    get(gamma0) = DataVector(number_of_points, 1.0);
    get(gamma1) = DataVector(number_of_points, -1.0);
    get(gamma2) = DataVector(number_of_points, 1.0);
  }

  size_t block_size = 0;

  void apply() {
    gh::TimeDerivative<gh::Solutions::all_solutions<Dim>, Dim>::apply(
        make_not_null(&dt_spacetime_metric), make_not_null(&dt_pi),
        make_not_null(&dt_phi),
        make_not_null(&get<gh::Tags::ConstraintGamma1>(buffer)),
        make_not_null(&get<gh::Tags::ConstraintGamma2>(buffer)),
        make_not_null(&get<gh::Tags::GaugeH<DataVector, Dim>>(buffer)),
        make_not_null(
            &get<gh::Tags::SpacetimeDerivGaugeH<DataVector, Dim>>(buffer)),
        make_not_null(&get<gh::Tags::Gamma1Gamma2>(buffer)),
        make_not_null(&get<gh::Tags::HalfPiTwoNormals>(buffer)),
        make_not_null(&get<gh::Tags::NormalDotOneIndexConstraint>(buffer)),
        make_not_null(&get<gh::Tags::Gamma1Plus1>(buffer)),
        make_not_null(&get<gh::Tags::PiOneNormal<Dim>>(buffer)),
        make_not_null(&get<gh::Tags::GaugeConstraint<DataVector, Dim>>(buffer)),
        make_not_null(&get<gh::Tags::HalfPhiTwoNormals<Dim>>(buffer)),
        make_not_null(
            &get<gh::Tags::ShiftDotThreeIndexConstraint<Dim>>(buffer)),
        make_not_null(
            &get<gh::Tags::MeshVelocityDotThreeIndexConstraint<Dim>>(buffer)),
        make_not_null(&get<gh::Tags::PhiOneNormal<Dim>>(buffer)),
        make_not_null(&get<gh::Tags::PiSecondIndexUp<Dim>>(buffer)),
        make_not_null(
            &get<gh::Tags::ThreeIndexConstraint<DataVector, Dim>>(buffer)),
        make_not_null(&get<gh::Tags::PhiFirstIndexUp<Dim>>(buffer)),
        make_not_null(&get<gh::Tags::PhiThirdIndexUp<Dim>>(buffer)),
        make_not_null(
            &get<gh::Tags::SpacetimeChristoffelFirstKindThirdIndexUp<Dim>>(
                buffer)),
        make_not_null(&get<gr::Tags::Lapse<DataVector>>(buffer)),
        make_not_null(&get<gr::Tags::Shift<DataVector, Dim>>(buffer)),
        make_not_null(
            &get<gr::Tags::InverseSpatialMetric<DataVector, Dim>>(buffer)),
        make_not_null(&get<gr::Tags::DetSpatialMetric<DataVector>>(buffer)),
        make_not_null(&get<gr::Tags::SqrtDetSpatialMetric<DataVector>>(buffer)),
        make_not_null(
            &get<gr::Tags::InverseSpacetimeMetric<DataVector, Dim>>(buffer)),
        make_not_null(
            &get<gr::Tags::SpacetimeChristoffelFirstKind<DataVector, Dim>>(
                buffer)),
        make_not_null(
            &get<gr::Tags::SpacetimeChristoffelSecondKind<DataVector, Dim>>(
                buffer)),
        make_not_null(
            &get<gr::Tags::TraceSpacetimeChristoffelFirstKind<DataVector, Dim>>(
                buffer)),
        make_not_null(
            &get<gr::Tags::SpacetimeNormalVector<DataVector, Dim>>(buffer)),
        d_spacetime_metric, d_pi, d_phi, spacetime_metric, pi, phi, gamma0,
        gamma1, gamma2,
        use_harmonic
            ? static_cast<const gh::gauges::GaugeCondition&>(harmonic_gauge)
            : gauge_condition,
        mesh, time, inertial_coords, inverse_jacobian, {});
  }

  double checksum() const {
    return get<0, 0>(dt_spacetime_metric)[0] + get<1, 2>(dt_pi)[3] +
           get<2, 0, 3>(dt_phi)[7];
  }

  Mesh<Dim> mesh;
  size_t number_of_points;
  tnsr::I<DataVector, Dim, Frame::Inertial> inertial_coords;
  InverseJacobian<DataVector, Dim, Frame::ElementLogical, Frame::Inertial>
      inverse_jacobian{};
  tnsr::aa<DataVector, Dim, Frame::Inertial> spacetime_metric;
  tnsr::aa<DataVector, Dim, Frame::Inertial> pi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> phi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> d_spacetime_metric;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> d_pi;
  tnsr::ijaa<DataVector, Dim, Frame::Inertial> d_phi;
  Scalar<DataVector> gamma0;
  Scalar<DataVector> gamma1;
  Scalar<DataVector> gamma2;
  gh::gauges::DampedHarmonic gauge_condition{100.0, std::array{1.2, 1.5, 1.7},
                                             std::array{2, 4, 6}};
  gh::gauges::Harmonic harmonic_gauge{};
  bool use_harmonic = std::getenv("GH_BENCH_HARMONIC") != nullptr;
  double time = 1.3;
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_spacetime_metric;
  tnsr::aa<DataVector, Dim, Frame::Inertial> dt_pi;
  tnsr::iaa<DataVector, Dim, Frame::Inertial> dt_phi;
  Variables<tmpl::list<
      gh::Tags::ConstraintGamma1, gh::Tags::ConstraintGamma2,
      gh::Tags::GaugeH<DataVector, Dim>,
      gh::Tags::SpacetimeDerivGaugeH<DataVector, Dim>, gh::Tags::Gamma1Gamma2,
      gh::Tags::HalfPiTwoNormals, gh::Tags::NormalDotOneIndexConstraint,
      gh::Tags::Gamma1Plus1, gh::Tags::PiOneNormal<Dim>,
      gh::Tags::GaugeConstraint<DataVector, Dim>,
      gh::Tags::HalfPhiTwoNormals<Dim>,
      gh::Tags::ShiftDotThreeIndexConstraint<Dim>,
      gh::Tags::MeshVelocityDotThreeIndexConstraint<Dim>,
      gh::Tags::PhiOneNormal<Dim>, gh::Tags::PiSecondIndexUp<Dim>,
      gh::Tags::ThreeIndexConstraint<DataVector, Dim>,
      gh::Tags::PhiFirstIndexUp<Dim>, gh::Tags::PhiThirdIndexUp<Dim>,
      gh::Tags::SpacetimeChristoffelFirstKindThirdIndexUp<Dim>,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
      gr::Tags::InverseSpatialMetric<DataVector, Dim>,
      gr::Tags::DetSpatialMetric<DataVector>,
      gr::Tags::SqrtDetSpatialMetric<DataVector>,
      gr::Tags::InverseSpacetimeMetric<DataVector, Dim>,
      gr::Tags::SpacetimeChristoffelFirstKind<DataVector, Dim>,
      gr::Tags::SpacetimeChristoffelSecondKind<DataVector, Dim>,
      gr::Tags::TraceSpacetimeChristoffelFirstKind<DataVector, Dim>,
      gr::Tags::SpacetimeNormalVector<DataVector, Dim>>>
      buffer;
};
}  // namespace

int main(const int argc, const char* const* const argv) {
  const size_t points_per_dimension =
      argc > 1 ? static_cast<size_t>(std::strtoull(argv[1], nullptr, 10)) : 10;
  const size_t iterations =
      argc > 2 ? static_cast<size_t>(std::strtoull(argv[2], nullptr, 10))
               : 10000;
  const size_t warmup_iterations =
      argc > 3 ? static_cast<size_t>(std::strtoull(argv[3], nullptr, 10)) : 200;
  // 0 = automatic; pass a value >= number of points to disable blocking
  const size_t block_size =
      argc > 4 ? static_cast<size_t>(std::strtoull(argv[4], nullptr, 10)) : 0;

  BenchmarkState state{points_per_dimension};
  state.block_size = block_size;
  for (size_t i = 0; i < warmup_iterations; ++i) {
    state.apply();
  }

  const auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < iterations; ++i) {
    state.apply();
  }
  const auto stop = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = stop - start;

  std::cout << "points_per_dimension: " << points_per_dimension << "\n"
            << "number_of_points: " << state.number_of_points << "\n"
            << "iterations: " << iterations << "\n"
            << "seconds: " << elapsed.count() << "\n"
            << "microseconds_per_call: "
            << 1.0e6 * elapsed.count() / static_cast<double>(iterations) << "\n"
            << "checksum: " << state.checksum() << "\n";
  return 0;
}
