// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"

/*!
 * \brief The set of utilities for performing CurvedScalarWave evolution with a
 * worldtube excision scheme.
 *
 * \details The worldtube excision scheme is a method that aims to enable NR
 * evolutions of intermediate mass ratio binary black hole simulations. In
 * standard BBH simulations two excision spheres are cut out from the domain
 * within the apparent horizons of the respective black holes. For larger mass
 * ratios, this introduces a scale disparity in the evolution system because the
 * small grid spacing in the elements near the smaller black hole constrain the
 * time step to be orders of magnitude smaller than near the larger black hole
 * due to the CFL condition. The worldtube excision scheme avoids this by
 * excising a much larger region (the worldtube) around the smaller black hole.
 * Since the excision boundary no longer lies within the apparent horizon,
 * boundary conditions are required. These are derived by approximating the
 * solution inside the worldtube using a perturbative solution - a black hole
 * perturbed by another black hole. The solution is calibrated by the evolved
 * metric on the worldtube boundary and in turn provides boundary conditions to
 * the evolution system.
 *
 * Here, we test this scheme using a toy problem of a scalar charge in
 * circular orbit around a Schwarzschild black hole.
 */
namespace CurvedScalarWave::Worldtube {

double worldtube_shrink_factor(const double orbit_radius,
                               const double original_orbit_radius);

double worldtube_shrink_factor_derivative(const double orbit_radius,
                                          const double orbit_velocity,
                                          const double original_orbit_radius);

template <size_t Dim>
std::tuple<tnsr::aa<double, Dim>, tnsr::AA<double, Dim>, tnsr::iaa<double, Dim>,
           tnsr::iAA<double, Dim>, tnsr::ijaa<double, Dim>,
           tnsr::ijAA<double, Dim>, tnsr::Abb<double, Dim>,
           tnsr::iAbb<double, Dim>, tnsr::A<double, Dim>, tnsr::iA<double, Dim>>
kerr_schild_quantities(const gr::Solutions::KerrSchild& kerr_schild,
                       const tnsr::I<double, Dim>& pos);

double broken_power_shrink(const double orbit_radius, const double amp,
                           const double rb);
double broken_power_shrink_derivative(const double orbit_radius,
                                      const double amp, const double rb,
                                      const double orbit_radius_derivative);

}  // namespace CurvedScalarWave::Worldtube
