// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Executables/ScalarWave/Batched/DeviceData.hpp"

namespace ScalarWave::Batched::Tags {

struct DeviceData : db::SimpleTag {
  using type = ScalarWave::Batched::DeviceData;
};

}  // namespace ScalarWave::Batched::Tags
