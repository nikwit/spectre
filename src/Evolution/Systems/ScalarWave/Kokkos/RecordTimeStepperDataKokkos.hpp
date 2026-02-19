// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Kokkos/RecordTimeStepperData.hpp"

namespace ScalarWave::Actions {

template <typename System>
using RecordTimeStepperDataKokkos = ::Actions::Kokkos::RecordTimeStepperData<
    System, KokkosTags::DeviceDtVariables, KokkosTags::DeviceDerivativeHistory>;

}  // namespace ScalarWave::Actions
