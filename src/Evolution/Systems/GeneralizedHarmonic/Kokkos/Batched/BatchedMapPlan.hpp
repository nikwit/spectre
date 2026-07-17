// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Block.hpp"

namespace gh::Actions {

/// \brief Classification of a block's grid-to-inertial coordinate map for the
/// batched GH frame-field update.
///
/// The batched executable can evaluate some coordinate maps directly on the
/// device (avoiding the per-stage host map loop and the associated
/// host-to-device transfers of coordinates, mesh velocity, and the composed
/// inverse Jacobian). This enum tags, per block, which device path applies.
///
/// The classification is fixed at initialization and never changes during an
/// evolution, so it can be computed once and reused (and, later, used to group
/// elements by map kind for per-group kernel launches).
enum class BatchedMapKind : int {
  /// No time dependence: the grid-to-inertial map is the identity, so inertial
  /// coordinates equal grid coordinates, the mesh velocity is zero, and the
  /// composed logical-to-inertial inverse Jacobian equals the static
  /// logical-to-grid inverse Jacobian. Evaluated on device.
  Static = 0,
  /// A time-dependent map that does not yet have a device implementation. The
  /// existing host map evaluation path is used for these blocks.
  HostEvaluated = 1,
};

/// \brief Classify a block's grid-to-inertial map for the batched device path.
///
/// Only static (non-time-dependent) blocks are currently evaluated on the
/// device. Every time-dependent block falls back to the host map path until
/// its map kind is ported (see the device-map-port milestones).
template <size_t Dim>
BatchedMapKind classify_block_map(const Block<Dim>& block) {
  return block.is_time_dependent() ? BatchedMapKind::HostEvaluated
                                   : BatchedMapKind::Static;
}

}  // namespace gh::Actions
