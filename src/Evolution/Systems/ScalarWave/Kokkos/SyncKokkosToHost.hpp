// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Kokkos/SyncKokkosToHost.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"

namespace ScalarWave::Events {

template <typename System>
using SyncKokkosToHost = ::Events::Kokkos::SyncKokkosToHost<
    System, ScalarWave::KokkosTags::DeviceVariables<System>>;

}  // namespace ScalarWave::Events
