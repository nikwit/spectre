// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/MomentFeedback.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <map>
#include <mutex>
#include <numbers>
#include <sstream>
#include <vector>
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Options/ParseError.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace gh::BoundaryConditions::detail {
// NPFBG001: 8-byte magic, uint64 sample count and point count, N weights,
// then (time, N point-major coordinate triples, N row-major 3x3 Cholesky
// factors) for every geometry sample. Native x86-64 doubles, generated and
// hashed by the research control_geometry extractor. Times are NR times.
// NPFBG002 appends uint64 map_bytes, fots_bytes, serialized coordinate map
// and functions of time, then N point-major grid-coordinate triples. The
// exact map avoids interpolation error in the coordinate identity check.
struct MomentFeedbackReference {
  std::vector<double> times, geometry_times;
  std::vector<gr::np::TidalMoments> moments;
  std::vector<std::vector<double>> geometry;
  DataVector weights;
  std::unique_ptr<domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>>
      map;
  domain::FunctionsOfTimeMap functions_of_time;
  tnsr::I<DataVector, 3, Frame::Grid> grid_coordinates;
};
namespace {
using Matrix = std::array<std::array<double, 3>, 3>;
Matrix polar(const Matrix& a) {
  double u[9], v[9], singular[3], work[3];
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      u[3 * i + j] = a[i][j];
  auto um = gsl_matrix_view_array(u, 3, 3), vm = gsl_matrix_view_array(v, 3, 3);
  auto sv = gsl_vector_view_array(singular, 3),
       wv = gsl_vector_view_array(work, 3);
  if (gsl_linalg_SV_decomp(&um.matrix, &vm.matrix, &sv.vector, &wv.vector) !=
          0 or
      not(*std::min_element(singular, singular + 3) >
          1.e-12 * (*std::max_element(singular, singular + 3))))
    ERROR("Singular feedback axis map");
  Matrix result{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        result[i][j] += u[3 * i + k] * v[3 * j + k];
  const double det =
      result[0][0] *
          (result[1][1] * result[2][2] - result[1][2] * result[2][1]) -
      result[0][1] *
          (result[1][0] * result[2][2] - result[1][2] * result[2][0]) +
      result[0][2] *
          (result[1][0] * result[2][1] - result[1][1] * result[2][0]);
  if (not(det > 0.))
    ERROR("Improper feedback axis map");
  return result;
}
gr::np::TidalMoments rotate(const gr::np::TidalMoments& q, const Matrix& r) {
  using ComplexMatrix = std::array<std::array<std::complex<double>, 3>, 3>;
  ComplexMatrix h{
      {{q[0], q[2], q[3]}, {q[2], q[1], q[4]}, {q[3], q[4], -q[0] - q[1]}}},
      rh{};
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          rh[i][j] += r[i][k] * h[k][l] * r[j][l];
  return {rh[0][0], rh[1][1], rh[0][1], rh[0][2], rh[1][2]};
}
size_t interval(const std::vector<double>& t, const double time) {
  if (not std::isfinite(time) or time < t.front() or time > t.back())
    ERROR("Feedback reference does not cover time "
          << time << " in [" << t.front() << "," << t.back() << "]");
  return std::min(
      static_cast<size_t>(std::upper_bound(t.begin(), t.end(), time) -
                          t.begin() - 1),
      t.size() - 2);
}
}  // namespace
MomentFeedback::MomentFeedback(Strength::type strength, const double start,
                               const double duration, std::string moments_file,
                               std::string geometry_file,
                               const Options::Context& context)
    : start_(start),
      duration_(duration),
      moments_file_(std::move(moments_file)),
      geometry_file_(std::move(geometry_file)) {
  strengths_ = std::holds_alternative<double>(strength)
                   ? std::array<double, 2>{{std::get<double>(strength),
                                            std::get<double>(strength)}}
                   : std::get<std::array<double, 2>>(strength);
  if (not std::isfinite(strengths_[0]) or strengths_[0] < 0. or
      strengths_[0] > 1. or not std::isfinite(strengths_[1]) or
      strengths_[1] < 0. or strengths_[1] > 1. or not std::isfinite(start_) or
      not std::isfinite(duration_) or duration_ <= 0. or
      not std::isfinite(start_ + duration_))
    PARSE_ERROR(context,
                "MomentFeedback requires finite Strength in [0,1], StartTime "
                "and positive RampDuration.");
  load();
}
void MomentFeedback::load() {
  static std::mutex mutex;
  static std::map<std::pair<std::string, std::string>,
                  std::weak_ptr<const MomentFeedbackReference>>
      cache;
  const std::lock_guard<std::mutex> lock(mutex);
  const auto key = std::make_pair(moments_file_, geometry_file_);
  if (const auto p = cache[key].lock()) {
    reference_ = p;
    return;
  }
  auto ref = std::make_shared<MomentFeedbackReference>();
  std::ifstream csv(moments_file_);
  if (not csv)
    ERROR("Cannot open feedback moments " << moments_file_);
  std::string line;
  std::getline(csv, line);
  if (line != "time,E11,E22,E12,E13,E23,B11,B22,B12,B13,B23")
    ERROR("Invalid feedback CSV schema");
  while (std::getline(csv, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    std::array<double, 11> values{};
    for (auto& v : values)
      if (not(row >> v) or not std::isfinite(v))
        ERROR("Invalid feedback CSV row");
    std::string extra;
    if (row >> extra)
      ERROR("Extra feedback CSV columns");
    if (not ref->times.empty() and values[0] <= ref->times.back())
      ERROR("Unordered feedback moments");
    ref->times.push_back(values[0]);
    gr::np::TidalMoments q;
    for (size_t a = 0; a < 5; ++a) {
      q[a] = {values[a + 1], values[a + 6]};
    }
    ref->moments.push_back(q);
  }
  if (ref->times.size() < 2)
    ERROR("Too few feedback moment samples");
  std::ifstream bin(geometry_file_, std::ios::binary);
  char magic[8];
  uint64_t ns = 0, np = 0;
  bin.read(magic, 8);
  bin.read(reinterpret_cast<char*>(&ns), 8);
  bin.read(reinterpret_cast<char*>(&np), 8);
  if (not bin or
      (std::string(magic, 8) != "NPFBG001" and
       std::string(magic, 8) != "NPFBG002") or
      ns < 6 or ns > 100000 or np < 6 or np > 100000)
    ERROR("Invalid feedback geometry header");
  const uint64_t expected_bytes = 24 + 8 * np + ns * (8 + 96 * np);
  const bool exact_coordinates = std::string(magic, 8) == "NPFBG002";
  const auto file_bytes = std::filesystem::file_size(geometry_file_);
  if ((not exact_coordinates and file_bytes != expected_bytes) or
      (exact_coordinates and file_bytes < expected_bytes + 16 + 24 * np)) {
    ERROR("Truncated or extra feedback geometry data");
  }
  ref->weights = DataVector(static_cast<size_t>(np));
  bin.read(reinterpret_cast<char*>(ref->weights.data()),
           static_cast<std::streamsize>(8 * np));
  if (not bin)
    ERROR("Truncated feedback weights");
  double total = 0.;
  for (double w : ref->weights) {
    if (not(w > 0. and std::isfinite(w)))
      ERROR("Invalid feedback weights");
    total += w;
  }
  if (abs(total - 1.) > 1.e-10)
    ERROR("Feedback weights must sum to one");
  for (uint64_t k = 0; k < ns; ++k) {
    double time = 0.;
    bin.read(reinterpret_cast<char*>(&time), 8);
    if (not bin)
      ERROR("Truncated feedback geometry time");
    if (not std::isfinite(time) or (not ref->geometry_times.empty() and
                                    time <= ref->geometry_times.back()))
      ERROR("Invalid feedback geometry times");
    ref->geometry_times.push_back(time);
    ref->geometry.emplace_back(static_cast<size_t>(12 * np));
    bin.read(reinterpret_cast<char*>(ref->geometry.back().data()),
             static_cast<std::streamsize>(96 * np));
    if (not bin)
      ERROR("Truncated feedback geometry");
    for (double x : ref->geometry.back())
      if (not std::isfinite(x))
        ERROR("Nonfinite feedback geometry");
  }
  if (exact_coordinates) {
    uint64_t map_size = 0, fots_size = 0;
    bin.read(reinterpret_cast<char*>(&map_size), 8);
    bin.read(reinterpret_cast<char*>(&fots_size), 8);
    if (not bin or map_size == 0 or fots_size == 0 or map_size > 134217728 or
        fots_size > 134217728 or
        file_bytes != expected_bytes + 16 + map_size + fots_size + 24 * np)
      ERROR("Invalid exact feedback coordinate metadata");
    std::vector<char> map_data(map_size), fots_data(fots_size);
    bin.read(map_data.data(), static_cast<std::streamsize>(map_size));
    bin.read(fots_data.data(), static_cast<std::streamsize>(fots_size));
    ref->map = deserialize<decltype(ref->map)>(map_data.data());
    ref->functions_of_time =
        deserialize<domain::FunctionsOfTimeMap>(fots_data.data());
    if (not ref->map)
      ERROR("Missing exact feedback coordinate map");
    ref->grid_coordinates = tnsr::I<DataVector, 3, Frame::Grid>(np);
    for (size_t p = 0; p < np; ++p)
      for (size_t d = 0; d < 3; ++d) {
        double x = 0.;
        bin.read(reinterpret_cast<char*>(&x), 8);
        if (not bin or not std::isfinite(x))
          ERROR("Invalid feedback grid coordinate");
        ref->grid_coordinates.get(d)[p] = x;
      }
    // Validate the appended map against every original coordinate sample.
    for (size_t k = 0; k < ref->geometry_times.size(); ++k) {
      const auto x = (*ref->map)(ref->grid_coordinates, ref->geometry_times[k],
                                 ref->functions_of_time);
      for (size_t p = 0; p < np; ++p)
        for (size_t d = 0; d < 3; ++d)
          if (abs(x.get(d)[p] - ref->geometry[k][3 * p + d]) > 2.e-8)
            ERROR("Exact feedback map disagrees with tabulated coordinates");
    }
  }
  if (bin.peek() != std::char_traits<char>::eof())
    ERROR("Extra feedback geometry data");
  reference_ = ref;
  cache[key] = ref;
}
double MomentFeedback::strength(const double time) const {
  const auto gains = strengths(time);
  return std::max(gains[0], gains[1]);
}
std::array<double, 2> MomentFeedback::strengths(const double time) const {
  if (time <= start_)
    return {{0., 0.}};
  if (time >= start_ + duration_)
    return strengths_;
  const double s = sin(0.5 * std::numbers::pi * (time - start_) / duration_);
  return {{strengths_[0] * s * s, strengths_[1] * s * s}};
}
gr::np::TidalMoments MomentFeedback::imposed(
    const double time, const gr::np::TidalMoments& live,
    const gr::np::TidalMoments& recorded, const gr::np::TidalMoments& pulse,
    const tnsr::ii<DataVector, 3>& metric,
    const tnsr::I<DataVector, 3>& coordinates,
    const DataVector& weights) const {
  if (not reference_)
    ERROR("Uninitialized feedback reference");
  const auto& ref = *reference_;
  const size_t np = ref.weights.size();
  if (weights.size() != np or get<0>(coordinates).size() != np or
      get<0, 0>(metric).size() != np)
    ERROR("Feedback control face resolution does not match");
  if (max(abs(weights - ref.weights)) > 1.e-12)
    ERROR("Feedback face quadrature does not match");
  const size_t qi = interval(ref.times, time),
               gi = interval(ref.geometry_times, time);
  const double fraction =
      (time - ref.times[qi]) / (ref.times[qi + 1] - ref.times[qi]);
  gr::np::TidalMoments control;
  for (size_t a = 0; a < 5; ++a)
    control[a] = (1. - fraction) * ref.moments[qi][a] +
                 fraction * ref.moments[qi + 1][a];
  const size_t start =
      std::min(gi > 2 ? gi - 2 : 0, ref.geometry_times.size() - 6);
  std::array<double, 6> w{};
  for (size_t i = 0; i < 6; ++i) {
    w[i] = 1.;
    for (size_t j = 0; j < 6; ++j)
      if (i != j)
        w[i] *= (time - ref.geometry_times[start + j]) /
                (ref.geometry_times[start + i] - ref.geometry_times[start + j]);
  }
  std::vector<double> geometry(12 * np, 0.);
  for (size_t j = 0; j < 6; ++j)
    for (size_t k = 0; k < geometry.size(); ++k)
      geometry[k] += w[j] * ref.geometry[start + j][k];
  tnsr::I<DataVector, 3> reference_coordinates(np);
  if (ref.map) {
    reference_coordinates =
        (*ref.map)(ref.grid_coordinates, time, ref.functions_of_time);
  } else {
    for (size_t p = 0; p < np; ++p)
      for (size_t d = 0; d < 3; ++d)
        reference_coordinates.get(d)[p] = geometry[3 * p + d];
  }
  const auto inverse =
      gr::np::inverse_lower_triangular(gr::np::cholesky_factor(metric));
  Matrix mean{};
  for (size_t p = 0; p < np; ++p) {
    for (size_t d = 0; d < 3; ++d)
      if (abs(coordinates.get(d)[p] - reference_coordinates.get(d)[p]) > 2.e-8)
        ERROR(
            "Feedback control coordinates or point ordering do not match the "
            "live face at point "
            << p << ", component " << d << ": live=" << coordinates.get(d)[p]
            << ", reference=" << reference_coordinates.get(d)[p]
            << ", time=" << time);
    Matrix a{};
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        for (size_t k = 0; k < 3; ++k)
          a[i][j] +=
              inverse.get(i, k)[p] * geometry[3 * np + 9 * p + 3 * k + j];
    const auto r = polar(a);
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        mean[i][j] += weights[p] * r[i][j];
  }
  const auto r = polar(mean);
  control = rotate(control, r);
  const auto probe = rotate(pulse, r);
  auto result = recorded;
  const auto gains = strengths(time);
  for (size_t a = 0; a < 5; ++a) {
    const auto delta = live[a] - control[a];
    // Retain the original arithmetic exactly for equal gains.
    result[a] +=
        (gains[0] == gains[1] ? gains[0] * delta
                              : std::complex<double>{gains[0] * delta.real(),
                                                     gains[1] * delta.imag()}) +
        probe[a];
  }
  return result;
}
void MomentFeedback::pup(PUP::er& p) {
  p | strengths_;
  p | start_;
  p | duration_;
  p | moments_file_;
  p | geometry_file_;
  if (p.isUnpacking())
    load();
}
}  // namespace gh::BoundaryConditions::detail
