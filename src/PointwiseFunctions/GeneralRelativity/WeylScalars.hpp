// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace gr {

template <typename Frame>
struct WeylScalarsResult {
  std::array<Scalar<ComplexDataVector>, 5> scalars;
  tnsr::I<ComplexDataVector, 3, Frame> m;
};

/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Computes Newman Penrose quantity \f$\Psi_4\f$ using the characteristic
 * field U\f$^{8+}\f$ and complex vector \f$\bar{m}^i\f$.
 *
 * \details Computes \f$\Psi_4\f$ as: \f$\Psi_4 =
 * U^{8+}_{ij}\bar{m}^i\bar{m}^j\f$ with the characteristic field
 * \f$U^{8+} = (P^{(a}_i P^{b)}_j - \frac{1}{2}P_{ij}P^{ab})
 * (E_{ab} - \epsilon_a^{cd}n_dB_{cb}\f$)
 * and \f$\bar{m}^i\f$ = \f$\frac{(x^i + iy^i)}{\sqrt{2}}\f$. \f$x^i\f$ and
 * \f$y^i\f$ are normalized unit vectors in the frame Frame.
 *
 */
template <typename Frame>
void weyl_scalars(
    gsl::not_null<WeylScalarsResult<Frame>*> weyl_scalars_result,
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3, Frame>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame>& spatial_normal_vector);

template <typename Frame>
WeylScalarsResult<Frame> weyl_scalars(
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3, Frame>& spatial_metric,
    const tnsr::I<DataVector, 3, Frame>& spatial_normal_vector);

/// @}

}  // namespace gr
