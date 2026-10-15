// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/SphereReplay.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <tuple>
#include <variant>
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/H5/VolumeData.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gh::worldtube {
namespace {
std::vector<std::string> component_names(const bool physical) {
  std::vector<std::string> result;
  const auto append = [&result](const auto& tensor, const std::string& name) {
    for (size_t c = 0; c < tensor.size(); ++c) {
      result.push_back(name + tensor.component_suffix(c));
    }
  };
  append(tnsr::aa<DataVector, 3>{}, "SpacetimeMetric");
  append(tnsr::aa<DataVector, 3>{}, "Pi");
  append(tnsr::iaa<DataVector, 3>{}, "Phi");
  if (physical) {
    append(tnsr::iaa<DataVector, 3>{}, "deriv(Pi)");
    append(tnsr::ijaa<DataVector, 3>{}, "deriv(Phi)");
    append(Scalar<DataVector>{}, "Kretschmann");
    append(tnsr::i<DataVector, 3>{}, "deriv(Kretschmann)");
  }
  return result;
}
DataVector read_component(const h5::VolumeData& volume, const size_t id,
                          const std::string& name) {
  auto component = volume.get_tensor_component(id, name);
  if (not std::holds_alternative<DataVector>(component.data)) {
    ERROR("Sphere replay requires double precision: " << name);
  }
  return std::get<DataVector>(std::move(component.data));
}
}  // namespace

SphereReplay::SphereReplay(const std::string& filename,
                           const std::string& subfile, const double radius,
                           const bool read_physical_data)
    : radius_(radius), read_physical_data_(read_physical_data) {
  if (not std::isfinite(radius_) or radius_ <= 0.) {
    ERROR("Sphere replay requires a positive finite radius.");
  }
  const h5::H5File<h5::AccessType::ReadOnly> file(filename);
  const auto& volume = file.get<h5::VolumeData>(subfile);
  const auto ids = volume.list_observation_ids();
  if (ids.size() < 6) {
    ERROR("Sphere replay needs at least six time samples.");
  }
  const auto names = component_names(read_physical_data_);
  std::string selected_grid;
  for (const auto id : ids) {
    const double time = volume.get_observation_value(id);
    if (not std::isfinite(time) or
        (not times_.empty() and time <= times_.back())) {
      ERROR("Sphere replay times must be finite and strictly increasing.");
    }
    const auto grids = volume.get_grid_names(id);
    const auto extents = volume.get_extents(id);
    const auto bases = volume.get_bases(id);
    const auto quadratures = volume.get_quadratures(id);
    const auto radii = read_component(volume, id, "InertialRadius");
    size_t offset = 0, selected_offset = 0, selected = grids.size(), total = 0;
    for (size_t g = 0; g < grids.size(); ++g) {
      const size_t size = std::accumulate(extents[g].begin(), extents[g].end(),
                                          size_t{1}, std::multiplies<>{});
      if (offset + size > radii.size()) {
        ERROR("Sphere replay has inconsistent grid extents.");
      }
      if (std::abs(radii[offset] - radius_) <= 1.e-12 * std::max(1., radius_)) {
        if (selected != grids.size())
          ERROR("Sphere replay radius is ambiguous.");
        selected = g;
        selected_offset = offset;
      }
      offset += size;
    }
    total = offset;
    if (selected == grids.size() or total != radii.size()) {
      ERROR("Sphere replay radius " << radius_ << " is absent or malformed.");
    }
    if (selected_grid.empty())
      selected_grid = grids[selected];
    if (grids[selected] != selected_grid or extents[selected].size() != 2 or
        bases[selected] != std::vector{Spectral::Basis::SphericalHarmonic,
                                       Spectral::Basis::SphericalHarmonic}) {
      ERROR("Sphere replay grid changed or is not spherical harmonic.");
    }
    const bool theta_fast =
        quadratures[selected] == std::vector{Spectral::Quadrature::Gauss,
                                             Spectral::Quadrature::Equiangular};
    if (not theta_fast and
        quadratures[selected] != std::vector{Spectral::Quadrature::Equiangular,
                                             Spectral::Quadrature::Gauss}) {
      ERROR("Sphere replay angular quadratures are unsupported.");
    }
    const size_t nt = extents[selected][theta_fast ? 0 : 1];
    const size_t np = extents[selected][theta_fast ? 1 : 0];
    if (nt < 3 or np != 2 * nt - 1)
      ERROR("Invalid sphere replay angular extents.");
    if (spherepack_ == nullptr)
      spherepack_ = std::make_unique<ylm::Spherepack>(nt - 1, nt - 1);
    if (spherepack_->l_max() != nt - 1)
      ERROR("Sphere replay resolution changed.");
    const size_t points = nt * np;
    const auto extract = [&](const std::string& name) {
      const auto raw = read_component(volume, id, name);
      if (raw.size() != total)
        ERROR("Sphere replay component has wrong size: " << name);
      DataVector result(points);
      for (size_t j = 0; j < np; ++j) {
        for (size_t i = 0; i < nt; ++i) {
          const double value =
              raw[selected_offset + (theta_fast ? i + nt * j : j + np * i)];
          if (not std::isfinite(value))
            ERROR("Non-finite sphere replay component: " << name);
          result[i + nt * j] = value;
        }
      }
      return result;
    };
    const auto recorded_radius = extract("InertialRadius");
    const auto& angles = spherepack_->theta_phi_points();
    std::array<double, 3> center{};
    for (size_t d = 0; d < 3; ++d) {
      const std::string suffix(1, "xyz"[d]);
      const auto centers = extract("InertialCenter_" + suffix);
      const auto coordinates = extract("InertialCoordinates_" + suffix);
      center[d] = centers[0];
      for (size_t q = 0; q < points; ++q) {
        const double unit =
            d == 2 ? cos(angles[0][q])
                   : sin(angles[0][q]) *
                         (d == 0 ? cos(angles[1][q]) : sin(angles[1][q]));
        if (std::abs(centers[q] - center[d]) > 1.e-12 or
            std::abs(recorded_radius[q] - radius_) > 1.e-12 or
            std::abs(coordinates[q] - center[d] - radius_ * unit) >
                1.e-10 * std::max(1., radius_)) {
          ERROR(
              "Sphere replay coordinates do not match the selected radius and "
              "angular ordering.");
        }
      }
    }
    times_.push_back(time);
    centers_.push_back(center);
    coefficients_.emplace_back(names.size());
    for (size_t c = 0; c < names.size(); ++c) {
      coefficients_.back()[c] = spherepack_->phys_to_spec(extract(names[c]));
    }
  }
  // GhWorldtube records on a uniform dense trigger. Detect missing records
  // rather than silently interpolating across a hole in the tape.
  const double interval = times_[1] - times_[0];
  const double tolerance =
      128. * std::numeric_limits<double>::epsilon() *
      std::max({1., std::abs(times_.front()), std::abs(times_.back())});
  for (size_t i = 2; i < times_.size(); ++i) {
    if (std::abs(times_[i] - times_[i - 1] - interval) > tolerance) {
      ERROR("Sphere replay has a missing sample or nonuniform cadence.");
    }
  }
}
std::array<double, 2> SphereReplay::time_bounds() const {
  return {times_.front(), times_.back()};
}
SphereReplay::TimeStencil SphereReplay::time_stencil(
    double time, const CenterTrajectory& center_at_time) const {
  const double tolerance =
      64. * std::numeric_limits<double>::epsilon() *
      std::max({1., std::abs(times_.front()), std::abs(times_.back())});
  if (not std::isfinite(time) or time < times_.front() - tolerance or
      time > times_.back() + tolerance) {
    ERROR("Sphere replay would extrapolate at time "
          << time << "; tape covers [" << times_.front() << ", "
          << times_.back() << "].");
  }
  time = std::clamp(time, times_.front(), times_.back());
  const auto upper = static_cast<size_t>(
      std::lower_bound(times_.begin(), times_.end(), time) - times_.begin());
  const size_t start = std::min(upper > 3 ? upper - 3 : 0, times_.size() - 6);
  std::array<double, 6> weights{};
  for (size_t i = 0; i < 6; ++i) {
    weights[i] = 1.;
    for (size_t j = 0; j < 6; ++j) {
      if (i != j)
        weights[i] *= (time - times_[start + j]) /
                      (times_[start + i] - times_[start + j]);
    }
  }
  std::array<double, 6> derivative_weights{};
  // Differentiate the Lagrange products directly. This is regular even when
  // time coincides with a recorded sample (including either tape endpoint).
  for (size_t i = 0; i < 6; ++i) {
    for (size_t k = 0; k < 6; ++k) {
      if (k == i)
        continue;
      double term = 1. / (times_[start + i] - times_[start + k]);
      for (size_t j = 0; j < 6; ++j) {
        if (j != i and j != k) {
          term *= (time - times_[start + j]) /
                  (times_[start + i] - times_[start + j]);
        }
      }
      derivative_weights[i] += term;
    }
  }
  // A coarse time interpolant need not reproduce the mapped trajectory,
  // especially across control-system updates. Compare at the recorded times
  // instead, and use the exact current center for the angular coordinates.
  for (size_t i = 0; i < 6; ++i) {
    const double sample_time = times_[start + i];
    const auto sample_center = center_at_time(sample_time);
    for (size_t d = 0; d < 3; ++d) {
      if (not std::isfinite(sample_center.get(d)) or
          std::abs(centers_[start + i][d] - sample_center.get(d)) >
              1.e-7 * std::max(1., radius_)) {
        ERROR("Receiver center does not follow the recorded sphere at time "
              << sample_time << ", component " << d << ": mapped "
              << sample_center.get(d) << ", recorded " << centers_[start + i][d]
              << ".");
      }
    }
  }
  return {start, weights, derivative_weights, center_at_time(time)};
}
ylm::Spherepack::InterpolationInfo<DataVector> SphereReplay::angular_info(
    const tnsr::I<DataVector, 3>& coordinates,
    const tnsr::I<double, 3>& center) const {
  const size_t points = get<0>(coordinates).size();
  std::array<DataVector, 2> angles{DataVector(points), DataVector(points)};
  for (size_t p = 0; p < points; ++p) {
    const double x = get<0>(coordinates)[p] - get<0>(center);
    const double y = get<1>(coordinates)[p] - get<1>(center);
    const double z = get<2>(coordinates)[p] - get<2>(center);
    const double r = std::hypot(x, y, z);
    if (not std::isfinite(r) or
        std::abs(r - radius_) > 1.e-9 * std::max(1., radius_)) {
      ERROR("Receiver face is not on the recorded inertial sphere: radius "
            << r << ", expected " << radius_);
    }
    angles[0][p] = std::acos(std::clamp(z / r, -1., 1.));
    angles[1][p] = std::atan2(y, x);
  }
  return spherepack_->set_up_interpolation_info(angles);
}
DataVector SphereReplay::interpolate_component(
    const size_t component, const TimeStencil& stencil,
    const ylm::Spherepack::InterpolationInfo<DataVector>& info,
    const size_t points, const bool time_derivative) const {
  DataVector coefs(spherepack_->spectral_size(), 0.);
  const auto& weights =
      time_derivative ? stencil.derivative_weights : stencil.weights;
  for (size_t i = 0; i < 6; ++i) {
    if (time_derivative) {
      // The derivative annihilates constants. Subtract a nearby sample
      // before summing to avoid cancellation of the large Coulomb background.
      coefs += weights[i] * (coefficients_[stencil.start + i][component] -
                             coefficients_[stencil.start + 2][component]);
    } else {
      coefs += weights[i] * coefficients_[stencil.start + i][component];
    }
  }
  DataVector result(points);
  spherepack_->interpolate_from_coefs(make_not_null(&result), coefs, info);
  return result;
}
Variables<ReplayTags> SphereReplay::evaluate(
    const double time, const tnsr::I<DataVector, 3>& coordinates,
    const CenterTrajectory& center_at_time) const {
  const auto stencil = time_stencil(time, center_at_time);
  const auto& center = stencil.center;
  {
    const std::lock_guard lock(evaluation_mutex_);
    if (last_evaluation_.has_value() and last_evaluation_->time == time and
        last_evaluation_->center == center and
        get<0>(last_evaluation_->coordinates).size() ==
            get<0>(coordinates).size() and
        last_evaluation_->coordinates == coordinates) {
      return last_evaluation_->fields;
    }
  }
  const size_t points = get<0>(coordinates).size();
  const auto info = angular_info(coordinates, center);
  Variables<ReplayTags> result(points);
  for (size_t c = 0; c < 50; ++c) {
    DataVector component(result.data() + c * points, points);
    component = interpolate_component(c, stencil, info, points);
  }
  {
    const std::lock_guard lock(evaluation_mutex_);
    last_evaluation_.emplace(Evaluation{time, coordinates, center, result});
  }
  return result;
}
RecordedPhysicalData SphereReplay::evaluate_physical(
    const double time, const tnsr::I<DataVector, 3>& coordinates,
    const CenterTrajectory& center_at_time,
    const tnsr::I<double, 3>& center_velocity) const {
  if (not read_physical_data_) {
    ERROR("Physical sphere replay data was not loaded.");
  }
  RecordedPhysicalData result{};
  result.fields = evaluate(time, coordinates, center_at_time);
  const auto stencil = time_stencil(time, center_at_time);
  const auto info = angular_info(coordinates, stencil.center);
  const size_t points = get<0>(coordinates).size();
  size_t component = 50;
  const auto fill = [&](auto& tensor) {
    for (auto& c : tensor) {
      c = interpolate_component(component++, stencil, info, points);
    }
  };
  fill(result.d_pi);
  fill(result.d_phi);
  const size_t k_component = component;
  fill(result.kretschmann);
  fill(result.d_kretschmann);
  get(result.dt_kretschmann) =
      interpolate_component(k_component, stencil, info, points, true);
  for (size_t d = 0; d < 3; ++d) {
    if (not std::isfinite(center_velocity.get(d))) {
      ERROR("Non-finite recording center velocity.");
    }
    get(result.dt_kretschmann) -=
        center_velocity.get(d) * result.d_kretschmann.get(d);
  }
  return result;
}
std::shared_ptr<const SphereReplay> cached_sphere_replay(
    const std::string& filename, const std::string& subfile,
    const double radius, const bool read_physical_data) {
  static std::mutex mutex;
  static std::map<std::tuple<std::string, std::string, double, bool>,
                  std::shared_ptr<const SphereReplay>>
      tapes;
  const std::lock_guard lock(mutex);
  const auto key =
      std::make_tuple(filename, subfile, radius, read_physical_data);
  const auto found = tapes.find(key);
  if (found != tapes.end())
    return found->second;
  auto tape = std::make_shared<SphereReplay>(filename, subfile, radius,
                                             read_physical_data);
  tapes.emplace(key, tape);
  return tape;
}
}  // namespace gh::worldtube
