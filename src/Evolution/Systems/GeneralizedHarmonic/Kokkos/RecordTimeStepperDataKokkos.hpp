// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Time/Kokkos/RecordTimeStepperData.hpp"

namespace gh::Actions {

template <typename System>
using RecordTimeStepperDataKokkos = ::Actions::Kokkos::RecordTimeStepperData<
    System, KokkosTags::DeviceDtVariables, KokkosTags::DeviceDerivativeHistory>;

}  // namespace gh::Actions
