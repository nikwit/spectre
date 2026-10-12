// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/ModalVector.hpp"
#include "NumericalAlgorithms/LinearOperators/Chalmers.hpp"
#include "NumericalAlgorithms/LinearOperators/CoefficientTransforms.hpp"
#include "NumericalAlgorithms/Spectral/CollocationPoints.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/QuadratureWeights.hpp"
#include "Utilities/MakeArray.hpp"

SPECTRE_TEST_CASE("Unit.LinearOperators.ChalmersParity",
                  "[Unit][NumericalAlgorithms]") {
  // Independent direct sums check the stride-two tail, both terminal
  // parities, and the small-sigma analytic summation branch.
  for (const size_t n : {6_st, 7_st, 12_st, 13_st}) {
    for (const double sigma : {0.005, 0.1, 0.10001, 0.7, 2.0}) {
      for (const size_t parity : {0_st, 1_st, 2_st}) {
        DataVector modes(n, 0.0);
        for (size_t k = 0; k < n; ++k) {
          if (parity == 2 or k % 2 == parity) {
            modes[k] = 2.3 * std::exp(-sigma * k);
          }
        }
        const auto fit = Chalmers::estimate_sequence_parity(modes);
        CHECK(fit.sigma == approx(sigma));
        CHECK(fit.max_log_fit_rms == approx(0.0));
        double squared_error = 0.0;
        for (size_t k = n - 2; k < 10000; ++k) {
          if (parity == 2 or k % 2 == parity) {
            squared_error +=
                2.0 * 2.3 * 2.3 * std::exp(-2.0 * sigma * k) / (2.0 * k + 1.0);
          }
        }
        CHECK(fit.error == approx(std::sqrt(squared_error)));
        for (const double scale : {1.e-100, 1.e100}) {
          const auto scaled = Chalmers::estimate_sequence_parity(scale * modes);
          CHECK(scaled.sigma == approx(sigma));
          CHECK(scaled.error / scale == approx(fit.error));
        }
      }
    }
  }
  // Roundoff in a symmetry-zero sequence cannot set the active decay rate.
  DataVector modes(11, 0.0);
  for (size_t k = 1; k < modes.size(); k += 2) {
    modes[k] = std::exp(-0.8 * k);
  }
  const auto clean = Chalmers::estimate_sequence_parity(modes, 3.e-5);
  for (const double noise : {1.e-15, 1.e-14, 1.e-13, 1.e-12}) {
    auto noisy = modes;
    for (size_t k = 0; k < noisy.size(); k += 2) {
      noisy[k] = noise;
    }
    const auto fit = Chalmers::estimate_sequence_parity(noisy, 3.e-5);
    CHECK(fit.sigma == approx(0.8));
    CHECK(fit.ignored_sequences == 1);
    CHECK(fit.error >= clean.error);
  }
  // A weaker but significant slow sequence must still control smoothness.
  modes = DataVector(16);
  for (size_t k = 0; k < modes.size(); ++k) {
    modes[k] = k % 2 == 0 ? std::exp(-1.2 * k) : 1.e-3 * std::exp(-0.15 * k);
  }
  CHECK(Chalmers::estimate_sequence_parity(modes, 1.e-5).sigma == approx(0.15));
  modes = DataVector(12, 0.0);
  modes[9] = 0.01;
  CHECK(std::isinf(Chalmers::estimate_sequence_parity(modes, 1.e-5).error));
  CHECK(
      Chalmers::estimate_sequence_parity(modes, 1.e-5).insufficient_sequences ==
      1);
  // Constant backgrounds do not keep a roundoff-only even variation active.
  modes = DataVector(10, 1.e-14);
  modes[0] = 1.0;
  CHECK(Chalmers::estimate_sequence_parity(modes, 1.e-5).error < 1.e-12);
  CHECK_THROWS_WITH(Chalmers::estimate_sequence_parity(modes, -1.0),
                    Catch::Matchers::ContainsSubstring("finite nonnegative"));
}

SPECTRE_TEST_CASE("Unit.LinearOperators.Chalmers",
                  "[Unit][NumericalAlgorithms]") {
  for (const double sigma : {0.005, 0.1, 0.10001, 0.5, 1.0, 2.0}) {
    DataVector coefficients(10);
    for (size_t k = 0; k < coefficients.size(); ++k) {
      coefficients[k] = 2.3 * std::exp(-sigma * static_cast<double>(k));
    }
    const auto fit = Chalmers::estimate_sequence(coefficients);
    CHECK(fit.sigma == approx(sigma));
    CHECK(fit.log_fit_rms == approx(0.0));
    CHECK(fit.significant_modes == 10);
    CHECK(fit.status == Chalmers::FitStatus::Fit);
    double squared_error = 0.0;
    for (size_t k = 9; k < 10000; ++k) {
      squared_error += 2.0 * 2.3 * 2.3 *
                       std::exp(-2.0 * sigma * static_cast<double>(k)) /
                       static_cast<double>(2 * k + 1);
    }
    CHECK(fit.error == approx(std::sqrt(squared_error)));
    for (const double scale : {1.e-100, 1.e100}) {
      const auto scaled = Chalmers::estimate_sequence(scale * coefficients);
      CHECK(scaled.sigma == approx(fit.sigma));
      CHECK(scaled.error / scale == approx(fit.error));
    }
  }
  CHECK(Chalmers::estimate_sequence(DataVector(8, 0.0)).error == 0.0);
  CHECK(Chalmers::estimate_sequence(DataVector{1.0, 2.0, 3.0, 0.0, 1.e-16})
            .status == Chalmers::FitStatus::TerminatedTail);
  const auto insufficient =
      Chalmers::estimate_sequence(DataVector{1.0, 0.0, 2.0});
  CHECK(insufficient.status == Chalmers::FitStatus::InsufficientModes);
  CHECK(std::isnan(insufficient.sigma));
  CHECK(std::isinf(insufficient.error));
  CHECK(Chalmers::estimate_sequence(DataVector(8, 1.0)).status ==
        Chalmers::FitStatus::Nondecaying);
  for (size_t parity = 0; parity < 2; ++parity) {
    DataVector coefficients(14, 0.0);
    for (size_t k = parity; k < coefficients.size(); k += 2) {
      coefficients[k] = std::exp(-0.4 * static_cast<double>(k));
    }
    const auto fit = Chalmers::estimate_sequence(coefficients);
    CHECK(fit.sigma == approx(0.4));
    CHECK(fit.significant_modes == 7);
  }
  DataVector plateau(16);
  for (size_t k = 0; k < plateau.size(); ++k) {
    plateau[k] = std::exp(-1.5 * static_cast<double>(std::min(k, size_t{8})));
  }
  const auto noisy_fit = Chalmers::estimate_sequence(plateau);
  CHECK(noisy_fit.sigma < 1.0);
  CHECK(noisy_fit.log_fit_rms > 1.0);

  const Mesh<3> mesh{{{6, 7, 8}},
                     Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  ModalVector coefficients(mesh.number_of_grid_points());
  for (size_t k = 0; k < 8; ++k) {
    for (size_t j = 0; j < 7; ++j) {
      for (size_t i = 0; i < 6; ++i) {
        coefficients[i + 6 * (j + 7 * k)] =
            std::exp(-0.3 * i - 1.2 * j - 0.7 * k);
      }
    }
  }
  const auto directional =
      Chalmers::estimate(to_nodal_coefficients(coefficients, mesh), mesh);
  const auto expected = std::array{0.3, 1.2, 0.7};
  for (size_t d = 0; d < 3; ++d) {
    // Small transverse stripes amplify nodal-to-modal roundoff in their logs.
    CHECK(directional[d].sigma == approx(expected[d]).epsilon(1.e-10));
    CHECK(directional[d].insufficient_sequences == 0);
  }
  const auto energy = Chalmers::estimate_directional_energy(
      to_nodal_coefficients(coefficients, mesh), mesh);
  for (size_t d = 0; d < 3; ++d) {
    double transverse_norm_squared = 1.0;
    for (size_t j = 0; j < 3; ++j) {
      if (j == d) {
        continue;
      }
      double norm_squared = 0.0;
      for (size_t k = 0; k < mesh.extents(j); ++k) {
        norm_squared +=
            2.0 * std::exp(-2.0 * expected[j] * k) / (2.0 * k + 1.0);
      }
      transverse_norm_squared *= norm_squared;
    }
    DataVector sequence(mesh.extents(d));
    for (size_t k = 0; k < sequence.size(); ++k) {
      sequence[k] =
          std::sqrt(transverse_norm_squared) * std::exp(-expected[d] * k);
    }
    const auto exact = Chalmers::estimate_sequence(sequence);
    CHECK(energy[d].sigma == approx(exact.sigma).epsilon(1.e-10));
    CHECK(energy[d].error == approx(exact.error).epsilon(1.e-10));
  }
  // Coarsening must account for the discarded mode even if the retained
  // spectrum is exactly constant and therefore has no fitted tail.
  const Mesh<1> coarsen_mesh{8, Spectral::Basis::Legendre,
                             Spectral::Quadrature::GaussLobatto};
  ModalVector coarsen_modes(8, 0.0);
  coarsen_modes[0] = 1.0;
  coarsen_modes[7] = 0.2;
  CHECK(Chalmers::p_coarsening_error(
            to_nodal_coefficients(coarsen_modes, coarsen_mesh),
            coarsen_mesh)[0] == approx(0.2 * std::sqrt(2.0 / 15.0)));
  // For an exponential, the candidate is the retained-spectrum indicator
  // plus the exact norm of the removed modal slab.
  for (size_t k = 0; k < 8; ++k) {
    coarsen_modes[k] = std::exp(-0.7 * k);
  }
  DataVector retained(7);
  for (size_t k = 0; k < 7; ++k) {
    retained[k] = coarsen_modes[k];
  }
  CHECK(Chalmers::p_coarsening_error(
            to_nodal_coefficients(coarsen_modes, coarsen_mesh),
            coarsen_mesh)[0] ==
        approx(Chalmers::estimate_sequence(retained).error +
               std::sqrt(2.0 / 15.0) * coarsen_modes[7]));
  // Full transverse Legendre weights are included, including for anisotropic
  // meshes. A tiny budget does not conceal resolved content.
  ModalVector slab(mesh.number_of_grid_points(), 0.0);
  slab[0] = 1.0;
  slab[2 + 6 * (3 + 7 * 7)] = 0.2;
  const auto slab_values = to_nodal_coefficients(slab, mesh);
  const auto slab_errors = Chalmers::p_coarsening_error(slab_values, mesh);
  CHECK(slab_errors[2] == approx(0.2 * std::sqrt(8.0 / (5.0 * 7.0 * 15.0))));
  CHECK(Chalmers::p_coarsening_error(slab_values, mesh, 1.e-8)[2] ==
        approx(slab_errors[2]));
  const Mesh<1> minimum_mesh{2, Spectral::Basis::Legendre,
                             Spectral::Quadrature::GaussLobatto};
  CHECK(std::isinf(
      Chalmers::p_coarsening_error(DataVector(2, 0.0), minimum_mesh)[0]));
  const Mesh<1> fit_floor_mesh{6, Spectral::Basis::Legendre,
                               Spectral::Quadrature::GaussLobatto};
  CHECK(std::isinf(
      Chalmers::p_coarsening_error(DataVector(6, 1.0), fit_floor_mesh)[0]));
  const Mesh<1> above_floor_mesh{7, Spectral::Basis::Legendre,
                                 Spectral::Quadrature::GaussLobatto};
  CHECK(Chalmers::p_coarsening_error(DataVector(7, 1.0), above_floor_mesh,
                                     1.e-12)[0] < 1.e-12);
  CHECK_THROWS_WITH(Chalmers::estimate_sequence(DataVector{1.0}),
                    Catch::Matchers::ContainsSubstring("at least two modes"));
  // Tiny transverse modal noise can give nondecaying or insufficient fits
  // even though its full norm is far below the requested accuracy.
  const Mesh<3> noise_mesh{8, Spectral::Basis::Legendre,
                           Spectral::Quadrature::GaussLobatto};
  ModalVector noisy_modes(noise_mesh.number_of_grid_points(), 0.0);
  noisy_modes[0] = 1.0;
  noisy_modes[1] = 0.5;
  noisy_modes[2] = 0.1;
  noisy_modes[7 + 8 * 7] = 1.e-10;
  const auto noisy_values = to_nodal_coefficients(noisy_modes, noise_mesh);
  const auto unguarded = Chalmers::estimate(noisy_values, noise_mesh);
  CHECK(std::isinf(unguarded[0].error));
  CHECK(std::isinf(unguarded[1].error));
  const auto guarded = Chalmers::estimate(noisy_values, noise_mesh, 1.e-7);
  for (const auto& direction : guarded) {
    CHECK(direction.error == 0.0);
    CHECK(direction.insufficient_sequences == 0);
    CHECK(direction.ignored_sequences > 0);
  }
  // A budget smaller than the perturbation must not hide it.
  CHECK(std::isinf(
      Chalmers::estimate(noisy_values, noise_mesh, 1.e-12)[0].error));
  // A large constant does not make a tiny transverse variation important.
  // Three points reproduce the formerly undetermined transverse fit.
  const Mesh<3> constant_mesh{{{8, 3, 3}},
                              Spectral::Basis::Legendre,
                              Spectral::Quadrature::GaussLobatto};
  ModalVector constant_modes(constant_mesh.number_of_grid_points(), 0.0);
  constant_modes[0] = 1.0;
  constant_modes[8 * 2] = 1.e-10;
  const auto constant_values =
      to_nodal_coefficients(constant_modes, constant_mesh);
  CHECK(
      std::isinf(Chalmers::estimate(constant_values, constant_mesh)[1].error));
  CHECK(Chalmers::estimate(constant_values, constant_mesh, 1.e-7)[1].error ==
        0.0);
  CHECK(std::isinf(
      Chalmers::estimate(constant_values, constant_mesh, 1.e-12)[1].error));
  CHECK(Chalmers::estimate_directional_energy(constant_values, constant_mesh,
                                              1.e-7)[1]
            .error == 0.0);
  CHECK(std::isinf(Chalmers::estimate_directional_energy(
                       constant_values, constant_mesh, 1.e-12)[1]
                       .error));
  CHECK_THROWS_WITH(
      Chalmers::estimate_directional_energy(constant_values, constant_mesh,
                                            -1.0),
      Catch::Matchers::ContainsSubstring("finite and nonnegative"));
  CHECK_THROWS_WITH(Chalmers::estimate_directional_energy(
                        DataVector(constant_mesh.number_of_grid_points(),
                                   std::numeric_limits<double>::quiet_NaN()),
                        constant_mesh, 1.e-7),
                    Catch::Matchers::ContainsSubstring("non-finite"));
  CHECK_THROWS_WITH(
      Chalmers::estimate(noisy_values, noise_mesh, -1.0),
      Catch::Matchers::ContainsSubstring("finite and nonnegative"));
}

SPECTRE_TEST_CASE("Unit.LinearOperators.ChalmersJoin",
                  "[Unit][NumericalAlgorithms]") {
  using Spectral::SegmentSize;
  const Mesh<1> low_mesh{6, Spectral::Basis::Legendre,
                         Spectral::Quadrature::GaussLobatto};
  const Mesh<1> high_mesh{10, Spectral::Basis::Legendre,
                          Spectral::Quadrature::GaussLobatto};
  const DataVector x_low = (Spectral::collocation_points(low_mesh) - 1.0) / 2.0;
  const DataVector x_high =
      (Spectral::collocation_points(high_mesh) + 1.0) / 2.0;
  const auto polynomial = Chalmers::h_coarsening_error(
      DataVector{1.0 + x_low * x_low * x_low}, low_mesh,
      DataVector{1.0 + x_high * x_high * x_high}, high_mesh, 0, 1.e-10);
  CHECK(polynomial.projection_error < 1.e-12);
  CHECK(polynomial.error < 1.e-12);
  CHECK(polynomial.parent_amplitude == approx(2.0));

  const Mesh<1> mesh{8, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  // Two individually constant children can contain a large unresolved jump.
  const auto jump = Chalmers::h_coarsening_error(
      DataVector(8, -1.0), mesh, DataVector(8, 1.0), mesh, 0, 1.e-10);
  const std::array<double, 4> exact_modes{{1.5, -0.875, 0.6875, -0.5859375}};
  double lost_norm_squared = 2.0;
  for (size_t j = 0; j < exact_modes.size(); ++j) {
    lost_norm_squared -=
        2.0 / (4.0 * j + 3.0) * exact_modes[j] * exact_modes[j];
  }
  CHECK(jump.projection_error == approx(std::sqrt(lost_norm_squared)));
  CHECK(jump.error >= jump.projection_error);
  const auto parity_jump = Chalmers::h_coarsening_error(
      DataVector(8, -1.0), mesh, DataVector(8, 1.0), mesh, 0, 1.e-10, 0, true);
  CHECK(parity_jump.projection_error == approx(jump.projection_error));
  CHECK(parity_jump.error >= parity_jump.projection_error);
  const Mesh<3> mesh3{8, Spectral::Basis::Legendre,
                      Spectral::Quadrature::GaussLobatto};
  for (size_t d = 0; d < 3; ++d) {
    const auto jump3 = Chalmers::h_coarsening_error(
        DataVector(512, -1.0), mesh3, DataVector(512, 1.0), mesh3, d, 1.e-10);
    CHECK(jump3.projection_error == approx(2.0 * jump.projection_error));
  }

  // Mixed transverse modes and unequal anisotropic child extents reproduce
  // a resolved parent along every join axis.
  const Mesh<3> low3{{{6, 8, 7}},
                     Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  const Mesh<3> high3{{{9, 7, 10}},
                      Spectral::Basis::Legendre,
                      Spectral::Quadrature::GaussLobatto};
  // Generate the data from a smaller exact polynomial mesh: parent-to-child
  // interpolation requires the child to have at least as many points.
  const Mesh<3> parent3{4, Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  ModalVector mixed_modes(parent3.number_of_grid_points(), 0.0);
  mixed_modes[0] = 1.0;
  mixed_modes[2 + 4 * (3 + 4 * 2)] = 0.2;
  mixed_modes[3 + 4 * (1 + 4 * 3)] = -0.1;
  const auto mixed = to_nodal_coefficients(mixed_modes, parent3);
  const auto full = make_array<3>(SegmentSize::Full);
  for (size_t d = 0; d < 3; ++d) {
    auto lo = full;
    auto hi = full;
    lo[d] = SegmentSize::LowerHalf;
    hi[d] = SegmentSize::UpperHalf;
    const auto reconstructed = Chalmers::h_coarsening_error(
        Spectral::project(mixed, parent3, low3, full, lo), low3,
        Spectral::project(mixed, parent3, high3, full, hi), high3, d, 1.e-10);
    CHECK(reconstructed.projection_error < 1.e-12);
    CHECK(reconstructed.error < 1.e-12);
    CHECK(reconstructed.parent_amplitude == approx(max(abs(mixed))));
  }

  // A lossless projection alone does not establish parent resolution.
  ModalVector modes(8, 0.0);
  modes[7] = 1.0;
  const auto parent = to_nodal_coefficients(modes, mesh);
  const auto lower =
      Spectral::project(parent, mesh, mesh, std::array{SegmentSize::Full},
                        std::array{SegmentSize::LowerHalf});
  const auto upper =
      Spectral::project(parent, mesh, mesh, std::array{SegmentSize::Full},
                        std::array{SegmentSize::UpperHalf});
  const auto underresolved =
      Chalmers::h_coarsening_error(lower, mesh, upper, mesh, 0, 1.e-10);
  CHECK(underresolved.projection_error < 1.e-12);
  CHECK(underresolved.error > 0.1);

  // Independently integrate the error with Gauss quadrature when the child
  // extents differ. Restricting to the smaller child's mesh would lose terms.
  ModalVector high_modes(10, 0.0);
  high_modes[9] = 1.0;
  const DataVector zero(6, 0.0);
  const auto high = to_nodal_coefficients(high_modes, high_mesh);
  const auto unequal =
      Chalmers::h_coarsening_error(zero, low_mesh, high, high_mesh, 0);
  const auto projected = Spectral::project(high, high_mesh, high_mesh,
                                           std::array{SegmentSize::UpperHalf},
                                           std::array{SegmentSize::Full});
  const Mesh<1> gauss{20, Spectral::Basis::Legendre,
                      Spectral::Quadrature::Gauss};
  double quadrature_error = 0.0;
  for (const auto half : {SegmentSize::LowerHalf, SegmentSize::UpperHalf}) {
    DataVector delta =
        Spectral::project(projected, high_mesh, gauss,
                          std::array{SegmentSize::Full}, std::array{half});
    if (half == SegmentSize::UpperHalf) {
      delta -= Spectral::project(high, high_mesh, gauss,
                                 std::array{SegmentSize::Full},
                                 std::array{SegmentSize::Full});
    }
    quadrature_error +=
        0.5 * sum(Spectral::quadrature_weights(gauss) * delta * delta);
  }
  CHECK(unequal.projection_error == approx(std::sqrt(quadrature_error)));
}

SPECTRE_TEST_CASE("Unit.LinearOperators.ChalmersFiltered",
                  "[Unit][NumericalAlgorithms]") {
  const Mesh<1> mesh{10, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  ModalVector modes(10, 0.0);
  DataVector usable(9);
  for (size_t k = 0; k < 9; ++k) {
    modes[k] = std::exp(-0.7 * k);
    usable[k] = modes[k];
  }
  const auto values = to_nodal_coefficients(modes, mesh);
  const auto expected = Chalmers::estimate_sequence(usable);
  CHECK(Chalmers::estimate(values, mesh, 0.0, 1)[0].sigma == approx(0.7));
  CHECK(Chalmers::estimate(values, mesh, 0.0, 1)[0].error ==
        approx(expected.error));
  CHECK(Chalmers::estimate_directional_energy(values, mesh, 0.0, 1)[0].error ==
        approx(expected.error));
  // Coarsening removes the zero highest mode, but filtering on the smaller
  // grid also removes the previously nonzero penultimate mode.
  DataVector retained(8);
  std::copy_n(usable.begin(), 8, retained.begin());
  CHECK(Chalmers::p_coarsening_error(values, mesh, 0.0, 1)[0] ==
        approx(std::sqrt(2.0 / 17.0) * usable[8] +
               Chalmers::estimate_sequence(retained).error));
  modes = 0.0;
  modes[0] = 1.0;
  modes[8] = 0.2;
  CHECK(Chalmers::p_coarsening_error(to_nodal_coefficients(modes, mesh), mesh,
                                     0.0, 1)[0] ==
        approx(0.2 * std::sqrt(2.0 / 17.0)));
  const Mesh<1> minimum{7, Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  CHECK(std::isinf(
      Chalmers::p_coarsening_error(DataVector(7, 1.0), minimum, 0.0, 1)[0]));
  // A parent can acquire a highest mode from restriction of its children.
  // The future parent filter loss must be counted even for a terminated fit.
  ModalVector parent_modes(10, 0.0);
  parent_modes[0] = 1.0;
  parent_modes[9] = 0.2;
  const auto parent = to_nodal_coefficients(parent_modes, mesh);
  const auto lower = Spectral::project(
      parent, mesh, mesh, std::array{Spectral::SegmentSize::Full},
      std::array{Spectral::SegmentSize::LowerHalf});
  const auto upper = Spectral::project(
      parent, mesh, mesh, std::array{Spectral::SegmentSize::Full},
      std::array{Spectral::SegmentSize::UpperHalf});
  const auto joined =
      Chalmers::h_coarsening_error(lower, mesh, upper, mesh, 0, 1.e-12, 1);
  CHECK(joined.projection_error < 1.e-12);
  CHECK(joined.error == approx(0.2 * std::sqrt(2.0 / 19.0)));
  // Orthogonal products have known norms. A mixed highest-mode product
  // belongs to two filtered slabs but must contribute to the loss only once.
  const auto check_intersecting_slabs = [](const auto dim_v) {
    constexpr size_t dim = decltype(dim_v)::value;
    const Mesh<dim> tensor_mesh{10, Spectral::Basis::Legendre,
                                Spectral::Quadrature::GaussLobatto};
    ModalVector tensor_modes(tensor_mesh.number_of_grid_points(), 0.0);
    tensor_modes[0] = 1.0;
    tensor_modes[9] = 0.2;
    tensor_modes[90] = 0.3;
    tensor_modes[99] = 0.4;
    const auto tensor_parent = to_nodal_coefficients(tensor_modes, tensor_mesh);
    const auto full = make_array<dim>(Spectral::SegmentSize::Full);
    auto lower_half = full;
    auto upper_half = full;
    lower_half[0] = Spectral::SegmentSize::LowerHalf;
    upper_half[0] = Spectral::SegmentSize::UpperHalf;
    const auto result = Chalmers::h_coarsening_error(
        Spectral::project(tensor_parent, tensor_mesh, tensor_mesh, full,
                          lower_half),
        tensor_mesh,
        Spectral::project(tensor_parent, tensor_mesh, tensor_mesh, full,
                          upper_half),
        tensor_mesh, 0, 1.e-12, 1);
    CHECK(result.projection_error < 1.e-12);
    CHECK(result.error ==
          approx(std::sqrt(static_cast<double>(size_t{1} << dim) *
                           ((0.04 + 0.09) / 19.0 + 0.16 / (19.0 * 19.0)))));
  };
  check_intersecting_slabs(std::integral_constant<size_t, 2>{});
  check_intersecting_slabs(std::integral_constant<size_t, 3>{});
  // Parity zeros plus a filtered highest mode must not look like a
  // terminated tail: the penultimate usable mode still carries a tail.
  modes = 0.0;
  for (size_t k = 1; k < 9; k += 2) {
    modes[k] = std::exp(-0.7 * k);
  }
  CHECK(Chalmers::estimate(to_nodal_coefficients(modes, mesh), mesh)[0].error ==
        0.0);
  CHECK(Chalmers::estimate(to_nodal_coefficients(modes, mesh), mesh, 0.0, 1)[0]
            .error > 0.0);
}
