// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosBoundaryCommunication.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Actions {

template <size_t Dim>
struct InitializeKokkosBoundaryCommunication {
  using outgoing_boundary_data_tag =
      ScalarWave::KokkosTags::OutgoingBoundaryCorrectionData<Dim>;
  using incoming_boundary_data_tag =
      ScalarWave::KokkosTags::IncomingBoundaryCorrectionData<Dim>;
  using external_boundary_data_tag =
      ScalarWave::KokkosTags::ExternalBoundaryCorrectionData<Dim>;

  using simple_tags =
      tmpl::list<outgoing_boundary_data_tag, incoming_boundary_data_tag,
                 external_boundary_data_tag>;
  using return_tags = simple_tags;
  using argument_tags = tmpl::list<>;

  static void apply(
      const gsl::not_null<typename outgoing_boundary_data_tag::type*>
          outgoing_boundary_data,
      const gsl::not_null<typename incoming_boundary_data_tag::type*>
          incoming_boundary_data,
      const gsl::not_null<typename external_boundary_data_tag::type*>
          external_boundary_data) {
    outgoing_boundary_data->clear();
    incoming_boundary_data->clear();
    external_boundary_data->clear();
  }
};

}  // namespace ScalarWave::Actions
