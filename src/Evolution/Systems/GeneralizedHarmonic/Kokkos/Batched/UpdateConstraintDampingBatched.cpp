// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/Batched/UpdateConstraintDampingBatched.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesKokkos.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache>
void UpdateConstraintDampingBatchedWithTag<
    FunctionsOfTimeTag, AddFunctionsOfTimeToConstGlobalCache>::
    apply(const gsl::not_null<typename device_constraint_gamma0_tag::type*>
              device_constraint_gamma0,
          const gsl::not_null<typename device_constraint_gamma1_tag::type*>
              device_constraint_gamma1,
          const gsl::not_null<typename device_constraint_gamma2_tag::type*>
              device_constraint_gamma2,
          const typename packed_geometry_tag::type& packed_geometry,
          const domain::FunctionsOfTimeMap& functions_of_time,
          const double time,
          const typename damping_function_gamma0_tag::DampingFunctionType&
              damping_function_gamma0,
          const typename damping_function_gamma1_tag::DampingFunctionType&
              damping_function_gamma1,
          const typename damping_function_gamma2_tag::DampingFunctionType&
              damping_function_gamma2) {
  const size_t total_points = packed_geometry.grid_coordinates_host[0].size();
  for (size_t d = 1; d < volume_dim; ++d) {
    ASSERT(packed_geometry.grid_coordinates_host[d].size() == total_points,
           "Packed grid-coordinate component sizes must match when updating "
           "constraint damping.");
  }

  Scalar<DataVector> host_gamma0{total_points};
  Scalar<DataVector> host_gamma1{total_points};
  Scalar<DataVector> host_gamma2{total_points};
  tnsr::I<DataVector, volume_dim, Frame::Grid> packed_grid_coordinates{
      total_points};
  for (size_t d = 0; d < volume_dim; ++d) {
    packed_grid_coordinates.get(d) = packed_geometry.grid_coordinates_host[d];
  }

  damping_function_gamma0(make_not_null(&host_gamma0), packed_grid_coordinates,
                          time, functions_of_time);
  damping_function_gamma1(make_not_null(&host_gamma1), packed_grid_coordinates,
                          time, functions_of_time);
  damping_function_gamma2(make_not_null(&host_gamma2), packed_grid_coordinates,
                          time, functions_of_time);

  using gamma0_space =
      typename device_constraint_gamma0_tag::type::value_type::memory_space;
  using gamma1_space =
      typename device_constraint_gamma1_tag::type::value_type::memory_space;
  using gamma2_space =
      typename device_constraint_gamma2_tag::type::value_type::memory_space;
  *device_constraint_gamma0 =
      copy_to_device(host_gamma0, tmpl::type_<gamma0_space>{});
  *device_constraint_gamma1 =
      copy_to_device(host_gamma1, tmpl::type_<gamma1_space>{});
  *device_constraint_gamma2 =
      copy_to_device(host_gamma2, tmpl::type_<gamma2_space>{});
}

template struct UpdateConstraintDampingBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize>;
template struct UpdateConstraintDampingBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize, false>;

}  // namespace gh::Actions
