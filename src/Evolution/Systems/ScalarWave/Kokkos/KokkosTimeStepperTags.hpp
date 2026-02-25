// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tags/MirrorView.hpp"
#include "Evolution/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"

namespace ScalarWave::KokkosTags {

template <typename System>
using DeviceVariables = evolution::Kokkos::Tags::DeviceVariables<System>;

template <typename System>
using DeviceDtVariables = evolution::Kokkos::Tags::DeviceDtVariables<System>;

template <typename System>
using DeviceStepStart = evolution::Kokkos::Tags::DeviceStepStart<System>;

template <typename System>
using DeviceDerivativeHistory =
    evolution::Kokkos::Tags::DeviceDerivativeHistory<System>;

template <size_t Dim>
using DeviceInverseJacobian =
    evolution::Kokkos::Tags::DeviceInverseJacobian<Dim>;

template <size_t Dim>
using DeviceFaceToVolumeIndexMap =
    evolution::Kokkos::Tags::DeviceFaceToVolumeIndexMap<Dim>;

template <size_t Dim>
using DeviceFaceUnitNormalCovector =
    evolution::Kokkos::Tags::DeviceFaceUnitNormalCovector<Dim>;

template <size_t Dim>
using DeviceFaceNormalMagnitude =
    evolution::Kokkos::Tags::DeviceFaceNormalMagnitude<Dim>;

struct DeviceConstraintGamma2 : db::SimpleTag {
  using type =
      typename ::Tags::MirrorView<ScalarWave::Tags::ConstraintGamma2>::type;
};

}  // namespace ScalarWave::KokkosTags
