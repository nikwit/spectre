// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/Factory.hpp"
#include "NumericalAlgorithms/LinearOperators/Filters/SphericalShell.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave {
/*!
 * \ingroup DiscontinuousGalerkinGroup
 * \brief A `tmpl::list` of all concrete filter types for the curved scalar
 * wave system.
 *
 * Extends `Filters::all_filters` with `Filters::SphericalShell` for `Dim == 3`.
 */
template <size_t Dim>
using all_filters = tmpl::append<
    Filters::all_filters<Dim, typename System<Dim>::variables_tag::tags_list>,
    tmpl::conditional_t<Dim == 3,
                        tmpl::list<Filters::SphericalShell<
                            typename System<Dim>::variables_tag::tags_list>>,
                        tmpl::list<>>>;
}  // namespace CurvedScalarWave
