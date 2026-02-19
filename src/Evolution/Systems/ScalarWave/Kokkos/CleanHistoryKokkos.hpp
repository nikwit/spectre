// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Kokkos/CleanHistory.hpp"

namespace ScalarWave::Actions {

template <typename System>
using CleanHistoryKokkos = ::Actions::Kokkos::CleanHistory<
    System, KokkosTags::DeviceVariables, KokkosTags::DeviceStepStart>;

}  // namespace ScalarWave::Actions
