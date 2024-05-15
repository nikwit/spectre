// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/YlmSpherepack.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/PostInterpolationCallback.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp {
namespace callbacks {

template <typename TagToObserve, typename InterpolationTargetTag,
          typename Frame>
struct ObserveYlms
    : tt::ConformsTo<intrp::protocols::PostInterpolationCallback> {
  static constexpr double fill_invalid_points_with =
      std::numeric_limits<double>::quiet_NaN();

  using const_global_cache_tags = tmpl::list<observers::Tags::SurfaceFileName>;

  template <typename DbTags, typename Metavariables, typename TemporalId>
  static void apply(const db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const TemporalId& temporal_id) {
    const auto& strahlkorper = get<StrahlkorperTags::Strahlkorper<Frame>>(box);
    const auto& ylm = strahlkorper.ylm_spherepack();
    const DataVector& collocation_values = get(get<TagToObserve>(box));
    const std::string& surface_name =
        pretty_type::name<InterpolationTargetTag>();
    const double time =
        InterpolationTarget_detail::get_temporal_id_value(temporal_id);
    auto& proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    std::vector<std::string> ylm_legend{};
    std::vector<double> ylm_data{};
    const DataVector spectral_data = ylm.phys_to_spec(collocation_values);
    const size_t l_max_output = 10;
    const size_t num_coefficients =
        YlmSpherepack::spectral_size(l_max_output, l_max_output) / 2;
    const size_t num_columns = num_coefficients + 1;
    ylm_legend.reserve(num_columns);
    ylm_data.reserve(num_columns);
    ylm_legend.emplace_back("Time");
    ylm_data.emplace_back(time);
    SpherepackIterator iter(strahlkorper.l_max(), strahlkorper.l_max());
    for (size_t l = 0; l <= l_max_output; l++) {
      for (int m = -static_cast<int>(l); m <= static_cast<int>(l); m++) {
        ylm_legend.emplace_back(MakeString{} << db::tag_name<TagToObserve>()
                                             << "(" << l << "," << m << ")");
        iter.set(l, m);
        ylm_data.emplace_back(spectral_data[iter()]);
      }
    }
    const std::string ylm_subfile_name{std::string{"/"} + surface_name +
                                       "_Ylm"};

    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        proxy[0], ylm_subfile_name, std::move(ylm_legend),
        std::make_tuple(std::move(ylm_data)));
  }
};
}  // namespace callbacks
}  // namespace intrp
