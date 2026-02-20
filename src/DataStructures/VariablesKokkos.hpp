// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>

#include "DataStructures/Tags/MirrorView.hpp"
#include "DataStructures/Tensor/AtIndex.hpp"
#include "DataStructures/Variables.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \brief Variables specialization that is backed by a Kokkos::View, meaning
 * it can be allocated on device (GPU) memory.
 *
 * Given a list of tags that hold `Tensor<Kokkos::View<double*>>`, this class
 * allocates a larger `Kokkos::View<double*[num_components]>` to hold all the
 * tensor components in contiguous memory. The memory is either allocated on the
 * host or on a device, depending on the `Kokkos::View` (see Kokkos
 * documentation for details). Individual tensors can be accessed with reference
 * semantics, meaning the `Kokkos::View` returned by `get` are subviews that
 * point into the large memory allocation.
 *
 * This class doesn't support all the arithmetic operations that the `Variables`
 * class backed by a `DataVector` supports. This is because arithmetic
 * operations should be done on device using Kokkos parallel algorithms with
 * precise control over launching kernels on device.
 */
template <typename... Tags, typename DataType, typename... Properties>
class Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>> {
 public:
  using tags_list = tmpl::list<Tags...>;
  static constexpr auto number_of_variables = sizeof...(Tags);
  static constexpr size_t number_of_independent_components =
      (... + Tags::type::size());

  using vector_type = Kokkos::View<DataType, Properties...>;
  static_assert(vector_type::rank() == 1,
                "Variables can only be used with a 1D Kokkos::View. "
                "A second dimension over tensor components is added "
                "automatically.");
  // Store the data on device in a 2D Kokkos::View over grid points and tensor
  // components. The two common layouts are:
  // - SoA: tensor components x grid points
  //   (typically better for derivatives because neighboring grid points
  //   are contiguous for each tensor component)
  // - AoS: grid points x tensor components
  //   (typically better for pointwise operations because all tensor components
  //   for a grid point are contiguous)
  // We choose AoS for now, but this can be changed if needed and should be
  // measured for performance.
  using storage_type =
      Kokkos::View<DataType[number_of_independent_components], Properties...>;
  using value_type = typename storage_type::value_type;

  Variables() = default;
  explicit Variables(size_t number_of_grid_points);
  explicit Variables(storage_type&& rhs);
  ~Variables() = default;

  void initialize(size_t number_of_grid_points) {
    Kokkos::resize(storage_, number_of_grid_points);
    set_reference_variable_data();
  }

  constexpr SPECTRE_ALWAYS_INLINE size_t size() const {
    return storage_.size();
  }
  constexpr SPECTRE_ALWAYS_INLINE size_t number_of_grid_points() const {
    return storage_.extent(0);
  }

  storage_type& view() { return storage_; }
  const storage_type& view() const { return storage_; }

  template <typename Tag, typename... FriendTags, typename FriendDataType,
            typename... FriendProperties>
  friend KOKKOS_INLINE_FUNCTION constexpr typename Tag::type& get(  // NOLINT
      Variables<tmpl::list<FriendTags...>,
                Kokkos::View<FriendDataType, FriendProperties...>>& v);
  template <typename Tag, typename... FriendTags, typename FriendDataType,
            typename... FriendProperties>
  friend KOKKOS_INLINE_FUNCTION constexpr const typename Tag::type&
  get(  // NOLINT
      const Variables<tmpl::list<FriendTags...>,
                      Kokkos::View<FriendDataType, FriendProperties...>>& v);

  template <typename Space>
  auto create_mirror_view(const Space& space) const {
    return Variables<tmpl::list<::Tags::MirrorView<Tags, Space>...>>(
        Kokkos::create_mirror_view(space, storage_));
  }

  template <typename Space>
  auto create_mirror_view_and_copy(const Space& space) const {
    return Variables<tmpl::list<::Tags::MirrorView<Tags, Space>...>>(
        Kokkos::create_mirror_view_and_copy(space, storage_));
  }

  void pup(PUP::er& p);  // NOLINT

 private:
  storage_type storage_{};
  tuples::TaggedTuple<Tags...> reference_variable_data_;

  void set_reference_variable_data();
};

template <typename... Tags, typename DataType, typename... Properties>
Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>>::
    Variables(const size_t number_of_grid_points)
    : storage_("Variables", number_of_grid_points) {
  set_reference_variable_data();
}

template <typename... Tags, typename DataType, typename... Properties>
Variables<tmpl::list<Tags...>,
          Kokkos::View<DataType, Properties...>>::Variables(storage_type&& rhs)
    : storage_(std::move(rhs)) {
  set_reference_variable_data();
}

template <typename... Tags, typename DataType, typename... Properties>
void Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>>::
    set_reference_variable_data() {
  size_t variable_offset = 0;
  tmpl::for_each<tags_list>([this, &variable_offset](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;
    auto& var = tuples::get<Tag>(reference_variable_data_);
    for (size_t i = 0; i < Tag::type::size(); ++i) {
      var[i] = Kokkos::subview(storage_, Kokkos::ALL(), variable_offset);
      ++variable_offset;
    }
  });
}

template <typename Tag, typename... Tags, typename DataType,
          typename... Properties>
KOKKOS_INLINE_FUNCTION constexpr typename Tag::type& get(
    Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>>& v) {
  static_assert(tmpl::list_contains_v<tmpl::list<Tags...>, Tag>,
                "Could not retrieve Tag from Variables. See the first "
                "template parameter of the instantiation for what Tag is "
                "being retrieved and the second template parameter for "
                "what Tags are available.");
  return tuples::get<Tag>(v.reference_variable_data_);
}

template <typename Tag, typename... Tags, typename DataType,
          typename... Properties>
KOKKOS_INLINE_FUNCTION constexpr const typename Tag::type& get(
    const Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>>&
        v) {
  static_assert(tmpl::list_contains_v<tmpl::list<Tags...>, Tag>,
                "Could not retrieve Tag from Variables. See the first "
                "template parameter of the instantiation for what Tag is "
                "being retrieved and the second template parameter for "
                "what Tags are available.");
  return tuples::get<Tag>(v.reference_variable_data_);
}

template <typename... Tags, typename... Is>
KOKKOS_FUNCTION auto make_at_index(const Variables<tmpl::list<Tags...>>& vars,
                                   const Is&... i) {
  tuples::TaggedTuple<::Tags::AtIndex<Tags>...> result{};
  (...,
   (get<::Tags::AtIndex<Tags>>(result) = make_at_index(get<Tags>(vars), i...)));
  return result;
}

template <typename... Tags>
auto copy_to_device(const Variables<tmpl::list<Tags...>>& vars) {
  Variables<tmpl::list<::Tags::MirrorView<Tags>...>> vars_on_device{
      vars.number_of_grid_points()};
  if constexpr (sizeof...(Tags) > 0) {
    const auto num_points = vars.number_of_grid_points();
    const auto num_components = vars_on_device.view().extent(1);
    ASSERT(vars.size() == num_points * num_components,
           "Host Variables size mismatch in copy_to_device: size="
               << vars.size() << " points=" << num_points
               << " components=" << num_components);

    auto vars_on_host = vars_on_device.create_mirror_view(Kokkos::HostSpace{});
    auto host_view = vars_on_host.view();
    const auto* const host_data = vars.data();
    for (size_t component = 0; component < num_components; ++component) {
      for (size_t point = 0; point < num_points; ++point) {
        host_view(point, component) = host_data[component * num_points + point];
      }
    }
    // Then copy to device
    Kokkos::deep_copy(vars_on_device.view(), vars_on_host.view());
  }
  return vars_on_device;
}

template <typename... Properties>
auto copy_to_device(const Tensor<DataVector, Properties...>& tensor) {
  Tensor<Kokkos::View<double*>, Properties...> tensor_on_device{
      "Tensor", tensor.begin()->size()};
  for (size_t i = 0; i < tensor.size(); ++i) {
    // First copy to a Kokkos::View on the host
    auto tensor_on_host = Kokkos::create_mirror_view(tensor_on_device[i]);
    std::copy_n(tensor[i].data(), tensor[i].size(), tensor_on_host.data());
    // Then copy to device
    Kokkos::deep_copy(tensor_on_device[i], tensor_on_host);
  }
  return tensor_on_device;
}

template <typename Space, typename... Properties>
auto copy_to_device(const Tensor<DataVector, Properties...>& tensor,
                    tmpl::type_<Space> /*meta*/) {
  Tensor<Kokkos::View<double*, Space>, Properties...> tensor_on_device{
      "Tensor", tensor.begin()->size()};
  for (size_t i = 0; i < tensor.size(); ++i) {
    // First copy to a Kokkos::View on the host
    auto tensor_on_host = Kokkos::create_mirror_view(tensor_on_device[i]);
    std::copy_n(tensor[i].data(), tensor[i].size(), tensor_on_host.data());
    // Then copy to device
    Kokkos::deep_copy(tensor_on_device[i], tensor_on_host);
  }
  return tensor_on_device;
}

template <typename HostTags, typename DeviceTags>
void copy_to_host(const gsl::not_null<Variables<HostTags>*> vars,
                  const Variables<DeviceTags>& vars_on_device) {
  if constexpr (tmpl::size<DeviceTags>::value > 0) {
    const auto num_points = vars_on_device.number_of_grid_points();
    const auto num_components = vars_on_device.view().extent(1);
    ASSERT(vars->number_of_grid_points() == num_points,
           "Host/device point count mismatch in copy_to_host: host="
               << vars->number_of_grid_points() << " device=" << num_points);
    ASSERT(vars->size() == num_points * num_components,
           "Host Variables size mismatch in copy_to_host: host size="
               << vars->size() << " points=" << num_points
               << " components=" << num_components);

    auto vars_on_host =
        vars_on_device.create_mirror_view_and_copy(Kokkos::HostSpace{});
    const auto host_view = vars_on_host.view();
    auto* const host_data = vars->data();
    for (size_t component = 0; component < num_components; ++component) {
      for (size_t point = 0; point < num_points; ++point) {
        host_data[component * num_points + point] = host_view(point, component);
      }
    }
  }
}

template <typename... Tags, typename DataType, typename... Properties>
void Variables<tmpl::list<Tags...>, Kokkos::View<DataType, Properties...>>::pup(
    PUP::er& p) {  // NOLINT
  (void)p;
  ERROR(
      "PUP for Variables<..., Kokkos::View<...>> is currently unsupported. "
      "This executable currently assumes no migration/checkpointing with "
      "device-resident DataBox items.");
}

namespace PUP {
template <typename DataType, typename... Properties>
void operator|(PUP::er& p,
               Kokkos::View<DataType, Properties...>& view) {  // NOLINT
  (void)p;
  (void)view;
  ERROR("PUP for Kokkos::View is currently unsupported in this build.");
}
}  // namespace PUP
