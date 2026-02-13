// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tags/MirrorView.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::KokkosTags {

template <typename System>
struct DeviceVariables : db::SimpleTag {
 private:
  using mirrored_variables_tag =
      db::add_tag_prefix<::Tags::MirrorView, typename System::variables_tag>;

 public:
  using type = typename mirrored_variables_tag::type;
};

template <typename System>
struct DeviceDtVariables : db::SimpleTag {
 private:
  using dt_variables_tag =
      db::add_tag_prefix<::Tags::dt, typename System::variables_tag>;
  using mirrored_dt_variables_tag =
      db::add_tag_prefix<::Tags::MirrorView, dt_variables_tag>;

 public:
  using type = typename mirrored_dt_variables_tag::type;
};

template <typename System>
struct DeviceStepStart : db::SimpleTag {
  using type = typename DeviceVariables<System>::type;
};

template <typename System>
struct DeviceDerivativeHistory : db::SimpleTag {
 private:
  using dt_vars_type = typename DeviceDtVariables<System>::type;

 public:
  using type = Kokkos::View<typename dt_vars_type::value_type***>;
};

template <size_t Dim>
struct DeviceInverseJacobian : db::SimpleTag {
 private:
  using host_tag = domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                                 Frame::Inertial>;

 public:
  using type = typename ::Tags::MirrorView<host_tag>::type;
};

struct DeviceConstraintGamma2 : db::SimpleTag {
  using type =
      typename ::Tags::MirrorView<ScalarWave::Tags::ConstraintGamma2>::type;
};

template <size_t Dim>
struct DeviceFaceToVolumeIndexMap : db::SimpleTag {
  using type =
      std::array<std::pair<Kokkos::View<size_t*>, Kokkos::View<size_t*>>, Dim>;
};

}  // namespace ScalarWave::KokkosTags
