// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Kokkos/CleanHistory.hpp"

namespace gh::Actions {

template <typename System>
using CleanHistoryKokkos = ::Actions::Kokkos::CleanHistory<
    System, KokkosTags::DeviceVariables, KokkosTags::DeviceStepStart>;

}  // namespace gh::Actions
