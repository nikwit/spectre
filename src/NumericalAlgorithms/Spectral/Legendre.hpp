// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/ModalVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace Spectral {

/// \brief Evaluate a tensor-product Legendre expansion at a logical point.
template <size_t Dim>
double evaluate_legendre_series(
    const ModalVector& coefficients, const Mesh<Dim>& mesh,
    const tnsr::I<double, Dim, Frame::ElementLogical>& logical_coords);

}  // namespace Spectral

