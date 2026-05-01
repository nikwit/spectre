// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Kokkos/KokkosTimeStepperTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Actions {

template <typename FunctionsOfTimeTag,
          bool AddFunctionsOfTimeToConstGlobalCache = true>
struct UpdateConstraintDampingBatchedWithTag {
 private:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
  using packed_geometry_tag = evolution::Kokkos::Tags::PackedGeometry<system>;
  using device_constraint_gamma0_tag = gh::KokkosTags::DeviceConstraintGamma0;
  using device_constraint_gamma1_tag = gh::KokkosTags::DeviceConstraintGamma1;
  using device_constraint_gamma2_tag = gh::KokkosTags::DeviceConstraintGamma2;
  using damping_function_gamma0_tag =
      gh::Tags::DampingFunctionGamma0<volume_dim, Frame::Grid>;
  using damping_function_gamma1_tag =
      gh::Tags::DampingFunctionGamma1<volume_dim, Frame::Grid>;
  using damping_function_gamma2_tag =
      gh::Tags::DampingFunctionGamma2<volume_dim, Frame::Grid>;

 public:
  using return_tags =
      tmpl::list<device_constraint_gamma0_tag, device_constraint_gamma1_tag,
                 device_constraint_gamma2_tag>;
  using argument_tags =
      tmpl::list<packed_geometry_tag, FunctionsOfTimeTag, ::Tags::Time,
                 damping_function_gamma0_tag, damping_function_gamma1_tag,
                 damping_function_gamma2_tag>;
  using const_global_cache_tags = tmpl::conditional_t<
      AddFunctionsOfTimeToConstGlobalCache,
      tmpl::list<FunctionsOfTimeTag, damping_function_gamma0_tag,
                 damping_function_gamma1_tag, damping_function_gamma2_tag>,
      tmpl::list<damping_function_gamma0_tag, damping_function_gamma1_tag,
                 damping_function_gamma2_tag>>;

  static void apply(
      gsl::not_null<typename device_constraint_gamma0_tag::type*>
          device_constraint_gamma0,
      gsl::not_null<typename device_constraint_gamma1_tag::type*>
          device_constraint_gamma1,
      gsl::not_null<typename device_constraint_gamma2_tag::type*>
          device_constraint_gamma2,
      const typename packed_geometry_tag::type& packed_geometry,
      const domain::FunctionsOfTimeMap& functions_of_time, double time,
      const typename damping_function_gamma0_tag::DampingFunctionType&
          damping_function_gamma0,
      const typename damping_function_gamma1_tag::DampingFunctionType&
          damping_function_gamma1,
      const typename damping_function_gamma2_tag::DampingFunctionType&
          damping_function_gamma2);
};

using UpdateConstraintDampingBatched = UpdateConstraintDampingBatchedWithTag<
    ::domain::Tags::FunctionsOfTimeInitialize>;

}  // namespace gh::Actions
