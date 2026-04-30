// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/DiscontinuousGalerkin/Initialization/DisableLtsOnNonconformingBoundaries.hpp"

#include <bit>
#include <cstddef>
#include <optional>

#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/FaceType.hpp"
#include "Time/Time.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace evolution::dg::Initialization {
template <size_t Dim>
void DisableLtsOnNonconformingBoundaries<Dim>::apply(
    const gsl::not_null<std::optional<size_t>*> fixed_lts_ratio,
    const Element<Dim>& element, const TimeDelta& time_step) {
  const bool touches_nonconforming_boundary =
      alg::any_of(element.face_types(), [](const auto& direction_and_type) {
        return direction_and_type.second ==
                   domain::FaceType::SingleNonconforming or
               direction_and_type.second ==
                   domain::FaceType::MultipleNonconforming;
      });

  if (not touches_nonconforming_boundary) {
    fixed_lts_ratio->reset();
    return;
  }

  const auto fraction = time_step.fraction();
  const auto numerator = fraction.numerator();
  const auto denominator = fraction.denominator();
  const auto positive_numerator = numerator < 0 ? -numerator : numerator;
  if (positive_numerator == 0) {
    ERROR("Cannot set a fixed LTS ratio from a zero time step.");
  }
  if (denominator % positive_numerator != 0) {
    ERROR("The time step " << time_step
                           << " is not an integer fraction of its slab.");
  }

  const auto ratio = denominator / positive_numerator;
  const auto ratio_as_size_t = static_cast<size_t>(ratio);
  if (ratio <= 0 or
      static_cast<decltype(ratio)>(ratio_as_size_t) != ratio) {
    ERROR("Invalid fixed LTS ratio " << ratio << " from time step "
                                     << time_step);
  }

  if (std::popcount(ratio_as_size_t) != 1) {
    ERROR("Fixed LTS ratio must be a power of 2, not " << ratio_as_size_t);
  }
  fixed_lts_ratio->emplace(ratio_as_size_t);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) \
  template struct DisableLtsOnNonconformingBoundaries<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE
#undef DIM
}  // namespace evolution::dg::Initialization
