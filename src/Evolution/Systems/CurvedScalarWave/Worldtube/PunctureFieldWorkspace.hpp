// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/DynamicBuffer.hpp"

namespace CurvedScalarWave::Worldtube {

/*!
 * \brief Reusable, element-local scratch storage for puncture evaluation.
 *
 * Contents are overwritten on every evaluation. Only a change in the number
 * of buffers or grid points allocates storage. No field values or orbit data
 * are cached, so reuse across time-stepper stages needs no invalidation.
 * Scratch contents are not checkpointed or migrated; storage is allocated
 * lazily on the receiving element's next evaluation.
 */
class PunctureFieldWorkspace {
 public:
  DynamicBuffer<DataVector>& buffer(const size_t number_of_vectors,
                                    const size_t number_of_grid_points) {
    if (not buffer_.has_value() or buffer_->size() != number_of_vectors or
        (number_of_vectors > 0 and
         buffer_->at(0).size() != number_of_grid_points)) {
      buffer_.emplace(number_of_vectors, number_of_grid_points);
    }
    return buffer_.value();
  }

  bool empty() const { return not buffer_.has_value(); }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    if (p.isUnpacking()) {
      buffer_.reset();
    }
  }

 private:
  std::optional<DynamicBuffer<DataVector>> buffer_{};
};

}  // namespace CurvedScalarWave::Worldtube
