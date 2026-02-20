// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"

namespace gh::Events {

template <typename System>
using SyncKokkosToHost = ::Events::Kokkos::SyncKokkosToHost<
    System, gh::KokkosTags::DeviceVariables<System>>;

}  // namespace gh::Events
