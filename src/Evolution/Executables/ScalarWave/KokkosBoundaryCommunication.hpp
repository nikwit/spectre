// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/Structure/DirectionalId.hpp"
#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Evolution/Systems/ScalarWave/BoundaryCorrections/UpwindPenalty.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace ScalarWave::KokkosTags {

template <size_t Dim>
using BoundaryCorrectionDataTags =
    db::wrap_tags_in<::Tags::MirrorView,
                     typename ScalarWave::BoundaryCorrections::UpwindPenalty<
                         Dim>::dg_package_field_tags>;

template <size_t Dim>
using BoundaryCorrectionDataBuffer = Variables<BoundaryCorrectionDataTags<Dim>>;

// Device-resident boundary-correction payload exchanged between elements.
// This type is intentionally same-node only for now.
template <size_t Dim>
struct BoundaryCorrectionData {
  Mesh<Dim> volume_mesh{};
  Mesh<Dim - 1> boundary_correction_mesh{};
  BoundaryCorrectionDataBuffer<Dim> boundary_correction_data{};
  TimeStepId validity_range{};
  size_t integration_order{std::numeric_limits<size_t>::max()};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | volume_mesh;
    p | boundary_correction_mesh;
    p | validity_range;
    p | integration_order;
    size_t boundary_data_size = boundary_correction_data.number_of_grid_points();
    p | boundary_data_size;
    if (boundary_data_size != 0) {
      ERROR("PUP for non-empty Kokkos boundary-correction payloads is "
            "currently unsupported. This path is intentionally same-node only "
            "for now.");
    }
  }
};

template <size_t Dim>
struct InboxBoundaryCorrectionData {
  using mapped_type =
      DirectionalIdMap<Dim, ScalarWave::KokkosTags::BoundaryCorrectionData<Dim>>;

  std::map<TimeStepId, mapped_type> messages;

  int missing_messages = 0;

  bool empty() const { return messages.empty(); }

  void collect_messages() { missing_messages = 0; }

  bool set_missing_messages(const size_t count) {
    missing_messages += static_cast<int>(count);
    return missing_messages <= 0;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | messages;
    p | missing_messages;
  }
};

template <size_t Dim>
struct OutgoingBoundaryCorrectionData : db::SimpleTag {
  using type =
      DirectionalIdMap<Dim, ScalarWave::KokkosTags::BoundaryCorrectionData<Dim>>;
};

template <size_t Dim>
struct IncomingBoundaryCorrectionData : db::SimpleTag {
  using type =
      DirectionalIdMap<Dim, ScalarWave::KokkosTags::BoundaryCorrectionData<Dim>>;
};

// Inbox tag used by the Kokkos-only scalar-wave DG communication path.
template <size_t Dim, bool UseNodegroupDgElements>
struct BoundaryCorrectionInbox {
  using stored_type = ScalarWave::KokkosTags::BoundaryCorrectionData<Dim>;

 public:
  using temporal_id = TimeStepId;
  using type = ScalarWave::KokkosTags::InboxBoundaryCorrectionData<Dim>;
  using value_type = type;

  static bool insert_into_inbox(
      const gsl::not_null<type*> inbox, const temporal_id& time_step_id,
      std::pair<DirectionalId<Dim>, stored_type> data) {
    auto& current_inbox = inbox->messages[time_step_id];
    if (not current_inbox.insert(std::move(data)).second) {
      ERROR("Failed to insert Kokkos boundary data into inbox at temporal id "
            << time_step_id << ".");
    }
    --inbox->missing_messages;
    return inbox->missing_messages == 0;
  }

  static std::string output_inbox(const type& inbox, const size_t padding_size) {
    std::stringstream ss{};
    const std::string pad(padding_size, ' ');
    ss << std::scientific << std::setprecision(16);
    ss << pad << "KokkosBoundaryCorrectionInbox:\n";

    for (const auto& [current_time_step_id, hash_map] : inbox.messages) {
      ss << pad << " Current time: " << current_time_step_id << "\n";
      for (const auto& [key, boundary_data] : hash_map) {
        ss << pad << "  Key: " << key
           << ", next time: " << boundary_data.validity_range << "\n";
      }
    }

    return ss.str();
  }
};

}  // namespace ScalarWave::KokkosTags
