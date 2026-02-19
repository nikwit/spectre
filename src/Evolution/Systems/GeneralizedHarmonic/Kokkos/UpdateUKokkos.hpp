// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Kokkos/UpdateU.hpp"

namespace gh::Actions {

template <typename System>
using UpdateUKokkos =
    ::Actions::Kokkos::UpdateU<System, KokkosTags::DeviceVariables,
                               KokkosTags::DeviceStepStart,
                               KokkosTags::DeviceDerivativeHistory>;

}  // namespace gh::Actions
