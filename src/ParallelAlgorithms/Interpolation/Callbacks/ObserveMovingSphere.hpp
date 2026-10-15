// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/PostInterpolationCallback.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"

namespace intrp::callbacks {
namespace observe_moving_sphere_detail {
/// One spherical grid per radius, with explicit inertial geometry. The angular
/// ordering follows the Sphere options; no angular filtering or
/// precision reduction is applied to the tensor components.
template <typename TagsList>
std::vector<ElementVolumeData> volume_data(
    const OptionHolders::Sphere& options, const std::array<double, 3>& center,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coordinates,
    const Variables<TagsList>& vars) {
  const size_t l_max = options.l_max;
  const size_t n = (l_max + 1) * (2 * l_max + 1);
  if (vars.number_of_grid_points() != n * options.radii.size()) {
    ERROR("MovingSphere output has an inconsistent number of points.");
  }
  const bool theta_fastest =
      options.angular_ordering == ylm::AngularOrdering::Strahlkorper;
  const std::vector<size_t> extents =
      theta_fastest ? std::vector<size_t>{l_max + 1, 2 * l_max + 1}
                    : std::vector<size_t>{2 * l_max + 1, l_max + 1};
  const std::vector<Spectral::Quadrature> quadratures =
      theta_fastest
          ? std::vector<Spectral::Quadrature>{Spectral::Quadrature::Gauss,
                                              Spectral::Quadrature::Equiangular}
          : std::vector<Spectral::Quadrature>{Spectral::Quadrature::Equiangular,
                                              Spectral::Quadrature::Gauss};
  std::vector<ElementVolumeData> result;
  size_t shell = 0;
  for (const double radius : options.radii) {
    std::vector<TensorComponent> components;
    for (size_t i = 0; i < 3; ++i) {
      DataVector coords(n);
      for (size_t p = 0; p < n; ++p) {
        coords[p] = coordinates.get(i)[shell * n + p];
      }
      const std::string suffix = std::string{"_"} + "xyz"[i];
      components.emplace_back("InertialCoordinates" + suffix,
                              std::move(coords));
      components.emplace_back("InertialCenter" + suffix,
                              DataVector(n, gsl::at(center, i)));
    }
    components.emplace_back("InertialRadius", DataVector(n, radius));
    tmpl::for_each<TagsList>([&](auto tag_v) {
      using Tag = tmpl::type_from<decltype(tag_v)>;
      const auto& tensor = get<Tag>(vars);
      for (size_t c = 0; c < tensor.size(); ++c) {
        DataVector values(n);
        for (size_t p = 0; p < n; ++p) {
          values[p] = tensor[c][shell * n + p];
          if (not std::isfinite(values[p])) {
            ERROR("MovingSphere received non-finite " << db::tag_name<Tag>());
          }
        }
        components.emplace_back(
            db::tag_name<Tag>() + tensor.component_suffix(c),
            std::move(values));
      }
    });
    result.emplace_back(
        "Sphere" + std::to_string(shell), std::move(components), extents,
        std::vector<Spectral::Basis>(2, Spectral::Basis::SphericalHarmonic),
        quadratures);
    ++shell;
  }
  return result;
}
}  // namespace observe_moving_sphere_detail

/// Write raw interpolated tensors and geometry to SurfaceFileName.h5, under
/// /<InterpolationTargetTag>.vol, using the standard volume format (version 0).
template <typename InterpolationTargetTag>
struct ObserveMovingSphere
    : tt::ConformsTo<intrp::protocols::PostInterpolationCallback> {
  // The interpolation framework fills missing points before invoking us. The
  // finite-data check below rejects incomplete surfaces instead of writing
  // them.
  static constexpr double fill_invalid_points_with =
      std::numeric_limits<double>::quiet_NaN();
  using const_global_cache_tags = tmpl::list<observers::Tags::SurfaceFileName>;

  template <typename DbTags, typename Metavariables, typename TemporalId>
  static void apply(const db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const TemporalId& temporal_id) {
    const double time =
        InterpolationTarget_detail::get_temporal_id_value(temporal_id);
    using target = typename InterpolationTargetTag::compute_target_points;
    const auto center = target::center(cache, time);
    const auto& options =
        Parallel::get<Tags::Sphere<InterpolationTargetTag>>(cache);
    const auto& cached_points = db::get<Tags::AllCoords<Frame::Inertial>>(box);
    const auto moved_points =
        target::points_at_time(cached_points, cache, time);
    auto data = observe_moving_sphere_detail::volume_data(
        options, center, moved_points ? *moved_points : cached_points,
        db::get<::Tags::Variables<
            typename InterpolationTargetTag::vars_to_interpolate_to_target>>(
            box));
    const std::string subfile =
        "/" + pretty_type::name<InterpolationTargetTag>();
    auto& writer = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    Parallel::threaded_action<observers::ThreadedActions::WriteVolumeData>(
        writer[0], Parallel::get<observers::Tags::SurfaceFileName>(cache),
        subfile, observers::ObservationId(time, subfile + ".vol"),
        std::move(data));
  }
};
}  // namespace intrp::callbacks
