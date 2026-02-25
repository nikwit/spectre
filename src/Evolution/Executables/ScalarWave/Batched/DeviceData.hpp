// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

#include <pup.h>

#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Evolution/Systems/ScalarWave/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/MortarHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/SegmentSize.hpp"
#include "Time/History.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace ScalarWave::Batched {

class DeviceData {
 public:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using device_variables_type =
      typename ScalarWave::KokkosTags::DeviceVariables<system>::type;
  using device_dt_variables_type =
      typename ScalarWave::KokkosTags::DeviceDtVariables<system>::type;
  using device_step_start_type =
      typename ScalarWave::KokkosTags::DeviceStepStart<system>::type;
  using device_derivative_history_type =
      typename ScalarWave::KokkosTags::DeviceDerivativeHistory<system>::type;
  using device_constraint_gamma2_type =
      typename ScalarWave::KokkosTags::DeviceConstraintGamma2::type;
  using device_face_to_volume_index_map_type =
      typename ScalarWave::KokkosTags::DeviceFaceToVolumeIndexMap<
          volume_dim>::type;
  using oriented_remote_face_index_map_type =
      DirectionalIdMap<volume_dim, ::Kokkos::View<size_t*>>;
  struct MortarMetadata {
    Mesh<volume_dim - 1> mortar_mesh{};
    std::array<Spectral::SegmentSize, volume_dim - 1> mortar_size{};
    bool needs_projection{false};
  };
  using mortar_metadata_map_type = DirectionalIdMap<volume_dim, MortarMetadata>;
  using oriented_remote_face_index_maps_type =
      std::vector<oriented_remote_face_index_map_type>;
  using mortar_metadata_maps_type = std::vector<mortar_metadata_map_type>;

  void initialize_from_domain(
      const Domain<volume_dim>& domain,
      const std::vector<std::array<size_t, volume_dim>>&
          initial_refinement_levels,
      const std::vector<std::array<size_t, volume_dim>>& initial_extents) {
    ASSERT(initial_refinement_levels.size() == domain.blocks().size(),
           "Initial refinement levels must contain one entry per block. Have "
               << initial_refinement_levels.size() << " levels for "
               << domain.blocks().size() << " blocks.");
    ASSERT(initial_extents.size() == domain.blocks().size(),
           "Initial extents must contain one entry per block. Have "
               << initial_extents.size() << " extents for "
               << domain.blocks().size() << " blocks.");

    local_element_ids_.clear();
    local_elements_.clear();
    element_index_by_id_.clear();
    element_extents_host_.clear();
    element_point_offsets_host_.clear();
    element_point_offsets_host_.push_back(0);
    points_per_element_ = 0;
    std::array<size_t, volume_dim> first_extents{};
    bool have_first_extents = false;

    for (const auto& block : domain.blocks()) {
      const size_t block_id = block.id();
      const auto block_element_ids = initial_element_ids(
          block_id, gsl::at(initial_refinement_levels, block_id));
      const auto& block_extents = gsl::at(initial_extents, block_id);
      if (not have_first_extents) {
        first_extents = block_extents;
        have_first_extents = true;
      } else {
        ASSERT(
            block_extents == first_extents,
            "EvolveScalarWaveKokkosBatched currently requires uniform extents "
            "across all blocks/elements. Found extents mismatch: first="
                << first_extents[0] << "," << first_extents[1] << ","
                << first_extents[2] << " current=" << block_extents[0] << ","
                << block_extents[1] << "," << block_extents[2] << ".");
      }
      const size_t block_points_per_element =
          block_extents[0] * block_extents[1] * block_extents[2];
      if (points_per_element_ == 0) {
        points_per_element_ = block_points_per_element;
      } else {
        ASSERT(points_per_element_ == block_points_per_element,
               "EvolveScalarWaveKokkosBatched currently requires uniform "
               "points per element.");
      }

      for (const auto& element_id : block_element_ids) {
        const auto element = ::domain::Initialization::create_initial_element(
            element_id, domain.blocks(), initial_refinement_levels);
        element_index_by_id_.insert_or_assign(element_id,
                                              local_element_ids_.size());
        local_element_ids_.push_back(element_id);
        local_elements_.push_back(element);
        element_extents_host_.push_back(block_extents);
        element_point_offsets_host_.push_back(
            element_point_offsets_host_.back() + points_per_element_);
      }
    }
    if (have_first_extents) {
      uniform_extents_host_ = first_extents;
    }

    total_points_ = element_point_offsets_host_.empty()
                        ? 0
                        : element_point_offsets_host_.back();
    allocate_metadata_on_device();
    allocate_fields_on_device();
  }

  size_t num_elements() const { return local_element_ids_.size(); }
  size_t total_points() const { return total_points_; }
  size_t points_per_element() const { return points_per_element_; }
  const std::array<size_t, volume_dim>& uniform_extents_host() const {
    return uniform_extents_host_;
  }
  const std::array<Spectral::Basis, volume_dim>& uniform_basis_host() const {
    return uniform_basis_host_;
  }
  const std::array<Spectral::Quadrature, volume_dim>&
  uniform_quadrature_host() const {
    return uniform_quadrature_host_;
  }
  void set_uniform_mesh_metadata(
      const std::array<size_t, volume_dim>& extents,
      const std::array<Spectral::Basis, volume_dim>& basis,
      const std::array<Spectral::Quadrature, volume_dim>& quadrature) {
    ASSERT(extents == uniform_extents_host_,
           "Uniform mesh extents mismatch in DeviceData metadata. Stored="
               << uniform_extents_host_[0] << "," << uniform_extents_host_[1]
               << "," << uniform_extents_host_[2] << " requested=" << extents[0]
               << "," << extents[1] << "," << extents[2] << ".");
    uniform_basis_host_ = basis;
    uniform_quadrature_host_ = quadrature;
  }
  void set_inertial_coordinates_host(
      std::array<DataVector, volume_dim> inertial_coordinates_host) {
    for (size_t d = 0; d < volume_dim; ++d) {
      ASSERT(inertial_coordinates_host[d].size() == total_points_,
             "Packed inertial-coordinate size mismatch for component " << d
                                                                        << ".");
    }
    inertial_coordinates_host_ = std::move(inertial_coordinates_host);
  }
  const std::array<DataVector, volume_dim>& inertial_coordinates_host() const {
    return inertial_coordinates_host_;
  }

  std::pair<size_t, size_t> element_point_range(
      const size_t element_index) const {
    ASSERT(element_index + 1 < element_point_offsets_host_.size(),
           "Element index " << element_index
                            << " is out of range for " << num_elements()
                            << " packed elements.");
    return {element_point_offsets_host_[element_index],
            element_point_offsets_host_[element_index + 1]};
  }

  const std::vector<ElementId<volume_dim>>& local_element_ids() const {
    return local_element_ids_;
  }
  const std::vector<Element<volume_dim>>& local_elements() const {
    return local_elements_;
  }
  size_t element_index(const ElementId<volume_dim>& element_id) const {
    const auto iter = element_index_by_id_.find(element_id);
    ASSERT(iter != element_index_by_id_.end(),
           "Packed-element index lookup failed for element id " << element_id
                                                                 << ".");
    return iter->second;
  }

  const std::vector<std::array<size_t, volume_dim>>& element_extents_host()
      const {
    return element_extents_host_;
  }

  const std::vector<size_t>& element_point_offsets_host() const {
    return element_point_offsets_host_;
  }

  const ::Kokkos::View<size_t*>& element_point_offsets_device() const {
    return element_point_offsets_device_;
  }

  const ::Kokkos::View<size_t*[volume_dim]>& element_extents_device() const {
    return element_extents_device_;
  }
  const device_face_to_volume_index_map_type& device_face_to_volume_index_map()
      const {
    return device_face_to_volume_index_map_;
  }
  const oriented_remote_face_index_map_type&
  oriented_remote_face_index_for_local(const size_t element_index) const {
    ASSERT(element_index < oriented_remote_face_index_for_local_.size(),
           "Element index " << element_index
                            << " out of range for oriented-face index maps.");
    return oriented_remote_face_index_for_local_[element_index];
  }
  const mortar_metadata_map_type& mortar_metadata(
      const size_t element_index) const {
    ASSERT(element_index < mortar_metadata_.size(),
           "Element index " << element_index
                            << " out of range for mortar metadata maps.");
    return mortar_metadata_[element_index];
  }
  ::Kokkos::View<double***>& element_inverse_jacobian_device() {
    return element_inverse_jacobian_device_;
  }
  const ::Kokkos::View<double***>& element_inverse_jacobian_device() const {
    return element_inverse_jacobian_device_;
  }

  device_variables_type& device_variables() { return device_variables_; }
  const device_variables_type& device_variables() const {
    return device_variables_;
  }
  auto device_variables_for_element(const size_t element_index) {
    const auto [begin, end] = element_point_range(element_index);
    return ::Kokkos::subview(device_variables_.view(),
                             ::Kokkos::pair<size_t, size_t>{begin, end},
                             ::Kokkos::ALL());
  }
  auto device_variables_for_element(const size_t element_index) const {
    const auto [begin, end] = element_point_range(element_index);
    return ::Kokkos::subview(device_variables_.view(),
                             ::Kokkos::pair<size_t, size_t>{begin, end},
                             ::Kokkos::ALL());
  }

  device_dt_variables_type& device_dt_variables() { return device_dt_variables_; }
  const device_dt_variables_type& device_dt_variables() const {
    return device_dt_variables_;
  }

  device_step_start_type& device_step_start() { return device_step_start_; }
  const device_step_start_type& device_step_start() const {
    return device_step_start_;
  }

  device_derivative_history_type& device_derivative_history() {
    return device_derivative_history_;
  }
  const device_derivative_history_type& device_derivative_history() const {
    return device_derivative_history_;
  }
  device_constraint_gamma2_type& device_constraint_gamma2() {
    return device_constraint_gamma2_;
  }
  const device_constraint_gamma2_type& device_constraint_gamma2() const {
    return device_constraint_gamma2_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | local_element_ids_;
    p | local_elements_;
    p | element_extents_host_;
    p | element_point_offsets_host_;
    p | total_points_;
    p | points_per_element_;
    p | uniform_extents_host_;
    p | uniform_basis_host_;
    p | uniform_quadrature_host_;
    p | inertial_coordinates_host_;
    p | element_index_by_id_;

    size_t device_points = device_variables_.number_of_grid_points();
    p | device_points;
    if (device_points != 0) {
      ERROR("PUP for non-empty ScalarWave::Batched::DeviceData is currently "
            "unsupported.");
    }
  }

 private:
  void allocate_metadata_on_device() {
    element_point_offsets_device_ =
        ::Kokkos::View<size_t*>("BatchedElementPointOffsets",
                                element_point_offsets_host_.size());
    auto host_point_offsets =
        ::Kokkos::create_mirror_view(element_point_offsets_device_);
    for (size_t i = 0; i < element_point_offsets_host_.size(); ++i) {
      host_point_offsets(i) = element_point_offsets_host_[i];
    }
    ::Kokkos::deep_copy(element_point_offsets_device_, host_point_offsets);

    element_extents_device_ =
        ::Kokkos::View<size_t*[volume_dim]>("BatchedElementExtents",
                                            element_extents_host_.size());
    auto host_extents = ::Kokkos::create_mirror_view(element_extents_device_);
    for (size_t e = 0; e < element_extents_host_.size(); ++e) {
      for (size_t d = 0; d < volume_dim; ++d) {
        host_extents(e, d) = element_extents_host_[e][d];
      }
    }
    ::Kokkos::deep_copy(element_extents_device_, host_extents);

    ASSERT(uniform_extents_host_[0] > 0 and uniform_extents_host_[1] > 0 and
               uniform_extents_host_[2] > 0,
           "Cannot initialize batched face-to-volume index maps before uniform "
           "extents are known.");
    const Index<volume_dim> extents_index{
        uniform_extents_host_[0], uniform_extents_host_[1],
        uniform_extents_host_[2]};
    const auto [indices_buffer, volume_and_slice_index_map] =
        volume_and_slice_indices(extents_index);
    (void)indices_buffer;
    for (size_t d = 0; d < volume_dim; ++d) {
      const size_t num_face_points = extents_index.slice_away(d).product();
      auto& lower_face_to_volume_index =
          gsl::at(device_face_to_volume_index_map_, d).first;
      auto& upper_face_to_volume_index =
          gsl::at(device_face_to_volume_index_map_, d).second;
      lower_face_to_volume_index = ::Kokkos::View<size_t*>(
          "BatchedLowerFaceToVolumeIndexMap", num_face_points);
      upper_face_to_volume_index = ::Kokkos::View<size_t*>(
          "BatchedUpperFaceToVolumeIndexMap", num_face_points);
      auto host_lower_face_to_volume_index =
          ::Kokkos::create_mirror_view(lower_face_to_volume_index);
      auto host_upper_face_to_volume_index =
          ::Kokkos::create_mirror_view(upper_face_to_volume_index);
      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).first) {
        host_lower_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }
      for (const auto& volume_and_slice_index :
           gsl::at(volume_and_slice_index_map, d).second) {
        host_upper_face_to_volume_index(volume_and_slice_index.second) =
            volume_and_slice_index.first;
      }
      ::Kokkos::deep_copy(lower_face_to_volume_index,
                          host_lower_face_to_volume_index);
      ::Kokkos::deep_copy(upper_face_to_volume_index,
                          host_upper_face_to_volume_index);
    }

    oriented_remote_face_index_for_local_.clear();
    oriented_remote_face_index_for_local_.resize(local_elements_.size());
    mortar_metadata_.clear();
    mortar_metadata_.resize(local_elements_.size());
    const Mesh<volume_dim> mesh{uniform_extents_host_, uniform_basis_host_,
                                uniform_quadrature_host_};
    for (size_t e = 0; e < local_elements_.size(); ++e) {
      const auto& element = local_elements_[e];
      for (const auto& [direction, neighbors] : element.neighbors()) {
        const size_t sliced_dim = direction.dimension();
        const auto slice_extents = extents_index.slice_away(sliced_dim);
        const Mesh<volume_dim - 1> local_face_mesh = mesh.slice_away(sliced_dim);
        for (const auto& neighbor : neighbors) {
          const auto& orientation = neighbors.orientation(neighbor);
          std::vector<double> local_grid_point_indices(slice_extents.product());
          std::iota(local_grid_point_indices.begin(),
                    local_grid_point_indices.end(), 0.0);
          const std::vector<double> local_index_in_remote_order =
              orient_variables_on_slice(local_grid_point_indices, slice_extents,
                                        sliced_dim, orientation);

          std::vector<size_t> remote_index_for_local(slice_extents.product(), 0);
          std::vector<bool> seen_local_index(slice_extents.product(), false);
          for (size_t remote_index = 0;
               remote_index < local_index_in_remote_order.size();
               ++remote_index) {
            const size_t local_index =
                static_cast<size_t>(local_index_in_remote_order[remote_index]);
            ASSERT(local_index < remote_index_for_local.size(),
                   "Oriented face-index map produced out-of-range local index.");
            ASSERT(not seen_local_index[local_index],
                   "Oriented face-index map produced duplicate local index.");
            remote_index_for_local[local_index] = remote_index;
            seen_local_index[local_index] = true;
          }
          for (size_t local_index = 0; local_index < seen_local_index.size();
               ++local_index) {
            ASSERT(seen_local_index[local_index],
                   "Oriented face-index map is missing local index "
                       << local_index << ".");
          }

          ::Kokkos::View<size_t*> device_map(
              "BatchedRemoteFaceIndexForLocalFaceIndex",
              remote_index_for_local.size());
          auto host_map = ::Kokkos::create_mirror_view(device_map);
          for (size_t i = 0; i < remote_index_for_local.size(); ++i) {
            host_map(i) = remote_index_for_local[i];
          }
          ::Kokkos::deep_copy(device_map, host_map);
          const DirectionalId<volume_dim> mortar_id{direction, neighbor};
          oriented_remote_face_index_for_local_[e].insert_or_assign(
              mortar_id, std::move(device_map));

          const Direction<volume_dim> remote_direction =
              orientation(direction.opposite());
          const Mesh<volume_dim - 1> remote_face_mesh =
              mesh.slice_away(remote_direction.dimension());
          const Mesh<volume_dim - 1> oriented_remote_face_mesh =
              orient_mesh_on_slice(remote_face_mesh, remote_direction.dimension(),
                                   orientation);
          MortarMetadata mortar_metadata{};
          mortar_metadata.mortar_mesh =
              ::dg::mortar_mesh(local_face_mesh, oriented_remote_face_mesh);
          mortar_metadata.mortar_size = ::dg::mortar_size(
              element.id(), neighbor, sliced_dim, orientation);
          mortar_metadata.needs_projection = Spectral::needs_projection(
              local_face_mesh, mortar_metadata.mortar_mesh,
              mortar_metadata.mortar_size);
          mortar_metadata_[e].insert_or_assign(mortar_id,
                                               std::move(mortar_metadata));
        }
      }
    }
  }

  void allocate_fields_on_device() {
    device_variables_ = device_variables_type(total_points_);
    device_dt_variables_ = device_dt_variables_type(total_points_);
    device_step_start_ = device_step_start_type(total_points_);
    device_constraint_gamma2_ =
        device_constraint_gamma2_type("BatchedConstraintGamma2", total_points_);
    device_derivative_history_ = device_derivative_history_type(
        "BatchedDerivativeHistory", TimeSteppers::history_max_substeps,
        total_points_,
        device_dt_variables_type::number_of_independent_components);
    element_inverse_jacobian_device_ = ::Kokkos::View<double***>(
        "BatchedElementInverseJacobian", num_elements(), points_per_element_, 9);

    if (total_points_ > 0) {
      ::Kokkos::deep_copy(device_variables_.view(), 0.0);
      ::Kokkos::deep_copy(device_dt_variables_.view(), 0.0);
      ::Kokkos::deep_copy(device_step_start_.view(), 0.0);
      ::Kokkos::deep_copy(get(device_constraint_gamma2_), 0.0);
      ::Kokkos::deep_copy(device_derivative_history_, 0.0);
      ::Kokkos::deep_copy(element_inverse_jacobian_device_, 0.0);
    }
  }

  std::vector<ElementId<volume_dim>> local_element_ids_{};
  std::vector<Element<volume_dim>> local_elements_{};
  std::map<ElementId<volume_dim>, size_t> element_index_by_id_{};
  std::vector<std::array<size_t, volume_dim>> element_extents_host_{};
  std::vector<size_t> element_point_offsets_host_{};
  size_t total_points_{0};
  size_t points_per_element_{0};
  std::array<size_t, volume_dim> uniform_extents_host_{{0, 0, 0}};
  std::array<Spectral::Basis, volume_dim> uniform_basis_host_{
      {Spectral::Basis::Legendre, Spectral::Basis::Legendre,
       Spectral::Basis::Legendre}};
  std::array<Spectral::Quadrature, volume_dim> uniform_quadrature_host_{
      {Spectral::Quadrature::GaussLobatto, Spectral::Quadrature::GaussLobatto,
       Spectral::Quadrature::GaussLobatto}};
  std::array<DataVector, volume_dim> inertial_coordinates_host_{};

  ::Kokkos::View<size_t*> element_point_offsets_device_{};
  ::Kokkos::View<size_t*[volume_dim]> element_extents_device_{};
  device_face_to_volume_index_map_type device_face_to_volume_index_map_{};
  oriented_remote_face_index_maps_type oriented_remote_face_index_for_local_{};
  mortar_metadata_maps_type mortar_metadata_{};
  ::Kokkos::View<double***> element_inverse_jacobian_device_{};
  device_variables_type device_variables_{};
  device_dt_variables_type device_dt_variables_{};
  device_step_start_type device_step_start_{};
  device_derivative_history_type device_derivative_history_{};
  device_constraint_gamma2_type device_constraint_gamma2_{};
};

}  // namespace ScalarWave::Batched
