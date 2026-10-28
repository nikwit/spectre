// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/SliceIterator.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/Tags/GaugeCondition.hpp"
#include "NumericalAlgorithms/Interpolation/LinearRegression.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/LinearOperators/PowerMonitors.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/NodalToModalMatrix.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "ParallelAlgorithms/Amr/Criteria/Constraints.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Factory.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/Math.hpp"

namespace gh::spec_amr {
RefinementVariables::RefinementVariables(const Mesh<3>& mesh) {
  for (size_t d = 0; d < 3; ++d) {
    if (mesh.basis(d) == Spectral::Basis::Legendre or
        mesh.basis(d) == Spectral::Basis::Chebyshev) {
      if (mesh.quadrature(d) != Spectral::Quadrature::Gauss and
          mesh.quadrature(d) != Spectral::Quadrature::GaussLobatto) {
        ERROR("Unsupported polynomial quadrature for SpecSpectral: " << mesh);
      }
      continue;
    }
    if (d == 1 and mesh.basis(1) == Spectral::Basis::SphericalHarmonic and
        mesh.basis(2) == Spectral::Basis::SphericalHarmonic and
        mesh.quadrature(1) == Spectral::Quadrature::Gauss and
        mesh.quadrature(2) == Spectral::Quadrature::Equiangular and
        mesh.extents(2) == 2 * mesh.extents(1) - 1) {
      number_of_variables_ = 2;
      kinds_[1] = Kind::SphericalHarmonic;
      variable_for_dimension_[2] = 1;
      return;
    }
    ERROR(
        "Unsupported basis combination for SpecSpectral refinement "
        "variables (S2 requires m_max == l_max): "
        << mesh);
  }
}

RefinementVariables::Kind RefinementVariables::kind(
    const size_t variable) const {
  ASSERT(variable < size(), "Invalid independent refinement variable");
  return kinds_[variable];
}

size_t RefinementVariables::variable_for_dimension(
    const size_t dimension) const {
  ASSERT(dimension < 3, "Invalid logical dimension");
  return variable_for_dimension_[dimension];
}

std::array<size_t, 3> RefinementVariables::mode_counts(
    const Mesh<3>& mesh) const {
  ASSERT(*this == RefinementVariables{mesh},
         "Mesh and independent refinement variables disagree");
  return size() == 2
             ? std::array<size_t, 3>{{mesh.extents(0), mesh.extents(1), 0}}
             : mesh.extents().indices();
}

std::array<amr::Flag, 3> RefinementVariables::logical_flags(
    const std::array<amr::Flag, 3>& independent_flags) const {
  if (size() == 2 and (independent_flags[1] == amr::Flag::Split or
                       independent_flags[1] == amr::Flag::Join)) {
    ERROR("Cannot h refine the S2 refinement variable");
  }
  std::array<amr::Flag, 3> result{};
  for (size_t d = 0; d < 3; ++d) {
    result[d] = independent_flags[variable_for_dimension(d)];
  }
  return result;
}

std::array<amr::Flag, 3> RefinementVariables::independent_flags(
    const std::array<amr::Flag, 3>& logical_flags) const {
  if (size() == 2 and logical_flags[1] != logical_flags[2]) {
    ERROR("S2 must have a single angular refinement decision");
  }
  auto result = logical_flags;
  if (size() == 2) {
    result[2] = amr::Flag::DoNothing;
    // Validate that the angular request does not change topology.
    (void)this->logical_flags(result);
  }
  return result;
}

void RefinementVariables::pup(PUP::er& p) {
  p | number_of_variables_;
  p | kinds_;
  p | variable_for_dimension_;
}

namespace {
constexpr double roundoff = 100.0 * std::numeric_limits<double>::epsilon();
double tail(const DataVector& power, const size_t count, const double floor) {
  if (count < 2) {
    return 1.0;
  }
  if (power[count - 1] <= floor and power[count - 2] <= floor) {
    return 0.0;
  }
  // A missing leading pair cannot demonstrate a small relative error.
  const double leading = std::max(power[0], power[1]);
  if (leading <= floor) {
    return 1.0;
  }
  double weight_sum = 0.0;
  double log_sum = 0.0;
  for (size_t k = 0; k < count; ++k) {
    if (power[k] <= floor) {
      continue;
    }
    const double weight =
        exp(-square(static_cast<double>(count - 1 - k) - 0.5));
    weight_sum += weight;
    log_sum += weight * log10(power[k] / leading);
  }
  return weight_sum > 0.0 ? pow(10.0, log_sum / weight_sum) : 1.0;
}
double predict(const double current, const double old) {
  return old <= 1.e-4 * current ? 2.0 * current : current * (current / old);
}
}  // namespace

Spectrum diagnose_power(const DataVector& power, const size_t filtered_modes) {
  ASSERT(power.size() >= filtered_modes + 3, "Insufficient retained modes");
  Spectrum result{};
  const size_t retained = power.size() - filtered_modes;
  double scale = 0.0;
  for (size_t k = 0; k < retained; ++k) {
    ASSERT(std::isfinite(power[k]) and power[k] >= 0.0, "Invalid spectrum");
    scale = std::max(scale, power[k]);
  }
  if (scale == 0.0) {
    return result;
  }
  const double floor = roundoff * scale;
  result.error = tail(power, retained, floor);
  result.next_error = std::max(result.error, tail(power, retained - 1, floor));
  std::vector<double> x{}, y{};
  for (size_t k = 0; k < retained; ++k) {
    if (power[k] > floor) {
      x.push_back(static_cast<double>(k));
      y.push_back(log10(power[k] / scale));
    }
  }
  if (x.size() < 2) {
    result.degenerate = result.next_error == 0.0;
    return result;
  }
  result.degenerate = false;
  if (x.size() < 4) {
    return result;  // insufficient evidence for a pile-up override
  }
  const auto fit = [&x, &y](const size_t begin, const size_t end) {
    if (end - begin == 2) {
      const double slope =
          (y[begin + 1] - y[begin]) / (x[begin + 1] - x[begin]);
      return intrp::LinearRegressionResult{y[begin] - slope * x[begin], slope,
                                           0.0, 0.0};
    }
    const auto first = static_cast<std::ptrdiff_t>(begin);
    const auto last = static_cast<std::ptrdiff_t>(end);
    return intrp::linear_regression(
        std::vector<double>(x.begin() + first, x.begin() + last),
        std::vector<double>(y.begin() + first, y.begin() + last));
  };
  std::vector<intrp::LinearRegressionResult> fits{
      fit(0, std::min(x.size(), std::max(x.size() / 3, size_t{4})))};
  for (size_t start = 0; start < 3; ++start) {
    for (size_t end = start + 5; end <= x.size(); ++end) {
      fits.push_back(fit(start, end));
    }
  }
  double max_error = 0.0;
  for (const auto& f : fits) {
    max_error = std::max(max_error, std::abs(f.delta_slope));
  }
  const double regularizer = std::max(1.e-15, 1.e-3 * max_error);
  double weight_sum = 0.0;
  for (const auto& f : fits) {
    const double weight = 1.0 / (regularizer + std::abs(f.delta_slope));
    result.convergence -= weight * f.slope;
    weight_sum += weight;
  }
  result.convergence /= weight_sum;
  // A flat or growing entire spectrum is unresolved, not excess resolution.
  // Do not count roundoff tails of a resolved polynomial as physical pile-up.
  if (result.convergence <= 1.e-10 or result.error == 0.0) {
    return result;
  }
  for (size_t start = 1; start + 1 < x.size(); ++start) {
    const double local =
        std::max(0.0, -fit(start, std::min(x.size(), start + 5)).slope);
    const double ratio = square(local / result.convergence);
    result.pile_up += ratio < 16.0 ? exp(-32.0 * ratio) : 0.0;
  }
  return result;
}

Spectrum diagnose_spec_power(
    const DataVector& power, const double filtered,
    const std::optional<double> after_p_down_filtered) {
  ASSERT(std::isfinite(filtered) and filtered >= 0.0 and
             (not after_p_down_filtered.has_value() or
              (std::isfinite(*after_p_down_filtered) and
               *after_p_down_filtered >= 0.0)),
         "Invalid filter count");
  Spectrum result{};
  std::vector<double> x{}, y{};
  for (size_t k = 0; k < power.size(); ++k) {
    ASSERT(std::isfinite(power[k]) and power[k] >= 0.0, "Invalid spectrum");
    if (power[k] != 0.0) {
      x.push_back(static_cast<double>(k + 1));
      y.push_back(log10(power[k]));
    }
  }
  // SpEC uses the count of nonzero coefficients here, even in a spectrum
  // containing exact parity gaps, rather than the original vector length.
  const double cutoff = static_cast<double>(y.size()) - round(filtered);
  const auto retained = static_cast<size_t>(
      std::count_if(x.begin(), x.end(),
                    [cutoff](const double index) { return index <= cutoff; }));
  if (retained == 0) {
    return result;
  }
  const auto [lo, hi] = std::minmax_element(y.begin(), y.end());
  if (*lo == *hi or *hi - *lo < 1.e-10 * std::abs(*lo + *hi)) {
    return result;
  }
  result.degenerate = false;
  const auto fit = [&x, &y](const size_t begin, const size_t end) {
    if (end < begin + 2) {
      return intrp::LinearRegressionResult{};
    }
    if (end == begin + 2) {
      const double slope =
          (y[begin + 1] - y[begin]) / (x[begin + 1] - x[begin]);
      return intrp::LinearRegressionResult{y[begin] - slope * x[begin], slope,
                                           0.0, 0.0};
    }
    const auto first = static_cast<std::ptrdiff_t>(begin);
    const auto last = static_cast<std::ptrdiff_t>(end);
    return intrp::linear_regression(
        std::vector<double>(x.begin() + first, x.begin() + last),
        std::vector<double>(y.begin() + first, y.begin() + last));
  };
  result.convergence = -fit(0, retained).slope;
  if (retained > 2) {
    std::vector<intrp::LinearRegressionResult> fits{
        fit(0, std::min(y.size(), std::max(retained / 3, size_t{4})))};
    for (size_t start = 0; start < 3; ++start) {
      for (size_t end = start + 5; end <= retained; ++end) {
        fits.push_back(fit(start, end));
      }
    }
    double max_error = 0.0;
    for (const auto& f : fits) {
      max_error = std::max(max_error, std::abs(f.delta_slope));
    }
    double weight_sum = 0.0;
    result.convergence = 0.0;
    for (const auto& f : fits) {
      // The limiting equal-weight case avoids SpEC's 0/0 on an exactly
      // representable straight line, without imposing an amplitude floor.
      const double weight =
          max_error == 0.0
              ? 1.0
              : 1.0 / (1.e-3 * max_error + std::abs(f.delta_slope));
      result.convergence -= weight * f.slope;
      weight_sum += weight;
    }
    result.convergence /= weight_sum;
    for (size_t start = 1; start + 1 < retained; ++start) {
      const double local =
          std::max(0.0, -fit(start, std::min(retained, start + 5)).slope);
      const double ratio = result.convergence == 0.0
                               ? std::numeric_limits<double>::infinity()
                               : square(local / result.convergence);
      result.pile_up += ratio < 16.0 ? exp(-32.0 * ratio) : 0.0;
    }
  }
  const auto relative_tail = [&y](const size_t count) {
    if (count == 0) {
      return 1.0;
    }
    if (count == 1) {
      return pow(10.0, y[0]);
    }
    if (count == 2) {
      return pow(10.0, y[1] - y[0]);
    }
    double log_sum = 0.0, weight_sum = 0.0;
    for (size_t k = 0; k < count; ++k) {
      const double weight = exp(
          -square(static_cast<double>(k) - (static_cast<double>(count) - 1.5)));
      log_sum += weight * y[k];
      weight_sum += weight;
    }
    return pow(10.0, log_sum / weight_sum - std::max(y[0], y[1]));
  };
  result.error = relative_tail(retained);
  // Generalize SpEC's one-less-unfiltered-mode forecast if this filter's
  // effective count changes with resolution. With a fixed count this is
  // exactly SpEC's retained-1, including in spectra with exact zeros.
  const auto next_retained = static_cast<size_t>(
      std::clamp(static_cast<double>(retained) - 1.0 + round(filtered) -
                     round(after_p_down_filtered.value_or(filtered)),
                 0.0, static_cast<double>(y.size())));
  result.next_error = std::max(result.error, relative_tail(next_retained));
  return result;
}

FilteredModes filtered_modes(
    const Mesh<3>& mesh,
    const Filters::Filter<
        3, ylm::TensorYlm::filter_detail::gh_spacetime_vars_list>& filter) {
  const RefinementVariables variables{mesh};
  const auto counts =
      [&](const Mesh<3>& local_mesh) -> std::optional<std::array<double, 3>> {
    const auto response = filter.volume_filter_response(local_mesh);
    if (not response.has_value()) {
      return std::nullopt;
    }
    const auto modes = variables.mode_counts(local_mesh);
    std::array<double, 3> result{};
    for (size_t d = 0; d < 3; ++d) {
      const size_t v = variables.variable_for_dimension(d);
      if ((*response)[d].size() != modes[v]) {
        ERROR("Inconsistent local filter response size for SpecSpectral");
      }
      double count = 0.0;
      for (const double r : (*response)[d]) {
        if (not std::isfinite(r) or r < 0.0 or r > 1.0) {
          ERROR("Invalid local modal filter response for SpecSpectral");
        }
        count += 1.0 - r;
      }
      if (d == 2 and variables.size() == 2 and
          (*response)[1] != (*response)[2]) {
        ERROR("S2 must have a single angular filter response");
      }
      result[v] = count;
    }
    return result;
  };
  FilteredModes result{};
  const auto current = counts(mesh);
  if (not current.has_value()) {
    ERROR("SpecSpectral requires a modal volume response from local filter "
          << filter.name() << " on mesh " << mesh);
  }
  result.current = *current;
  for (size_t v = 0; v < variables.size(); ++v) {
    auto extents = mesh.extents();
    ASSERT(extents[v] > 1, "Insufficient modes for p-down forecast");
    --extents[v];
    if (variables.kind(v) == RefinementVariables::Kind::SphericalHarmonic) {
      extents[2] -= 2;
    }
    const auto lower =
        counts(Mesh<3>{extents.indices(), mesh.basis(), mesh.quadrature()});
    result.can_decrease[v] = lower.has_value();
    result.after_p_down[v] =
        lower.has_value() ? (*lower)[v] : result.current[v];
  }
  return result;
}

std::array<DataVector, 3> polynomial_power(const DataVector& component,
                                           const Mesh<3>& mesh) {
  ASSERT(RefinementVariables{mesh}.size() == 3,
         "Polynomial power requires three polynomial directions");
  std::array<DataVector, 3> result{};
  const Matrix identity{};
  for (size_t d = 0; d < 3; ++d) {
    auto matrices = make_array<3>(std::cref(identity));
    matrices[d] =
        std::cref(Spectral::nodal_to_modal_matrix(mesh.slice_through(d)));
    const auto modal = apply_matrices(matrices, component, mesh.extents());
    result[d] = DataVector(mesh.extents(d), 0.0);
    for (size_t k = 0; k < mesh.extents(d); ++k) {
      for (SliceIterator si(mesh.extents(), d, k); si; ++si) {
        result[d][k] += square(modal[si.volume_offset()]);
      }
    }
    result[d] = sqrt(result[d] / static_cast<double>(
                                     mesh.extents().slice_away(d).product()));
  }
  return result;
}

std::array<double, 3> aggregate_constraints(
    const std::array<double, 3>& directional_rms,
    const RefinementVariables& variables) {
  std::array<double, 3> result{};
  std::array<size_t, 3> counts{};
  for (size_t d = 0; d < 3; ++d) {
    const size_t v = variables.variable_for_dimension(d);
    result[v] += directional_rms[d];
    ++counts[v];
  }
  for (size_t v = 0; v < variables.size(); ++v) {
    result[v] = sqrt(result[v] / static_cast<double>(counts[v]));
  }
  return result;
}

namespace {
std::pair<tnsr::ia<DataVector, 3>, tnsr::iaa<DataVector, 3>> gh_constraints(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const tnsr::I<DataVector, 3>& coordinates, const double time,
    const gauges::GaugeCondition& gauge, const Scalar<DataVector>& gamma2) {
  const auto spatial = gr::spatial_metric(metric);
  const auto [det_spatial, inverse_spatial] = determinant_and_inverse(spatial);
  const auto shift = gr::shift(metric, inverse_spatial);
  const auto lapse = gr::lapse(shift, metric);
  const auto normal = gr::spacetime_normal_vector(lapse, shift);
  const auto normal_form =
      gr::spacetime_normal_one_form<DataVector, 3, Frame::Inertial>(lapse);
  const auto inverse_metric =
      gr::inverse_spacetime_metric(lapse, shift, inverse_spatial);
  using GaugeCompute =
      gauges::Tags::GaugeAndDerivativeCompute<3,
                                              gh::Solutions::all_solutions<3>>;
  GaugeCompute::return_type gauge_data{};
  GaugeCompute::function(make_not_null(&gauge_data), lapse, shift, normal,
                         Scalar<DataVector>{sqrt(get(det_spatial))},
                         inverse_spatial, metric, pi, phi, mesh, time,
                         coordinates, inverse_jacobian, gauge);
  auto c3 = gh::three_index_constraint(
      partial_derivative(metric, mesh, inverse_jacobian), phi);
  auto c2 = gh::two_index_constraint(
      get<gh::Tags::SpacetimeDerivGaugeH<DataVector, 3>>(gauge_data),
      normal_form, normal, inverse_spatial, inverse_metric, pi, phi,
      partial_derivative(pi, mesh, inverse_jacobian),
      partial_derivative(phi, mesh, inverse_jacobian), gamma2, c3);
  return {std::move(c2), std::move(c3)};
}
}  // namespace

Diagnostics diagnose(const tnsr::aa<DataVector, 3>& metric,
                     const tnsr::aa<DataVector, 3>& pi,
                     const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
                     const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                                           Frame::Inertial>& inverse_jacobian,
                     const tnsr::I<DataVector, 3>& coordinates,
                     const double time, const gauges::GaugeCondition& gauge,
                     const Scalar<DataVector>& gamma2) {
  Diagnostics result{};
  result.refinement_variables = RefinementVariables{mesh};
  if (result.refinement_variables.size() != 3) {
    ERROR(
        "Legacy GH diagnostics support polynomial meshes; use diagnose_spec "
        "for spherical shells");
  }
  std::array<std::array<DataVector, 3>, 2> power{};
  for (auto& group : power) {
    for (size_t d = 0; d < 3; ++d) {
      group[d] = DataVector(mesh.extents(d), 0.0);
    }
  }
  const auto accumulate = [&mesh, &power](const auto& tensor,
                                          const size_t group) {
    for (const auto& component : tensor) {
      const auto modes = PowerMonitors::power_monitors(component, mesh);
      for (size_t d = 0; d < 3; ++d) {
        power[group][d] += square(modes[d]);
      }
    }
  };
  accumulate(metric, 0);
  accumulate(pi, 1);
  accumulate(phi, 1);
  for (size_t group = 0; group < 2; ++group) {
    for (size_t d = 0; d < 3; ++d) {
      result.spectra[group][d] = diagnose_power(
          DataVector(sqrt(power[group][d] / (group == 0 ? 10.0 : 40.0))));
    }
  }

  const auto [c2, c3] = gh_constraints(metric, pi, phi, mesh, inverse_jacobian,
                                       coordinates, time, gauge, gamma2);
  const auto jacobian = determinant_and_inverse(inverse_jacobian).second;
  tnsr::i<DataVector, 3, Frame::ElementLogical> normalization{}, projected{};
  DataVector buffer{};
  amr::Criteria::Constraints_detail::normalization_factor_square(
      make_not_null(&normalization), jacobian);
  const auto project = [&](const auto& constraint, const size_t group) {
    amr::Criteria::Constraints_detail::logical_constraints(
        make_not_null(&projected), make_not_null(&buffer), constraint, jacobian,
        normalization);
    for (size_t d = 0; d < 3; ++d) {
      result.constraints[group][d] =
          sqrt(sum(square(projected.get(d))) /
               static_cast<double>(mesh.number_of_grid_points()));
    }
  };
  project(c2, 0);
  project(c3, 1);
  return result;
}

Diagnostics diagnose_spec(
    const tnsr::aa<DataVector, 3>& metric, const tnsr::aa<DataVector, 3>& pi,
    const tnsr::iaa<DataVector, 3>& phi, const Mesh<3>& mesh,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const Jacobian<DataVector, 3, Frame::Grid, Frame::Inertial>& grid_jacobian,
    const tnsr::I<DataVector, 3>& coordinates, const double time,
    const gauges::GaugeCondition& gauge, const Scalar<DataVector>& gamma2,
    const FilteredModes& filter,
    const gsl::not_null<power_monitor::CartToSphereMatrices*> transform_cache) {
  Diagnostics result{};
  result.refinement_variables = RefinementVariables{mesh};
  result.filter = filter;
  std::array<std::array<DataVector, 3>, 2> power{};
  if (result.refinement_variables.size() == 2) {
    Variables<ylm::TensorYlm::filter_detail::gh_spacetime_vars_list> fields(
        mesh.number_of_grid_points());
    get<gr::Tags::SpacetimeMetric<DataVector, 3>>(fields) = metric;
    get<gh::Tags::Pi<DataVector, 3>>(fields) = pi;
    get<gh::Tags::Phi<DataVector, 3>>(fields) = phi;
    auto shell = power_monitor::gh_shell_power_monitors(
        transform_cache, fields, mesh, grid_jacobian,
        power_monitor::AngularNormalization::SpEC);
    for (size_t v = 0; v < 2; ++v) {
      power[0][v] = std::move(shell.spacetime_metric[v]);
      power[1][v] = std::move(shell.pi_phi[v]);
    }
  } else {
    for (auto& group : power) {
      for (size_t d = 0; d < 3; ++d) {
        group[d] = DataVector(mesh.extents(d), 0.0);
      }
    }
    const auto accumulate = [&mesh, &power](const auto& tensor,
                                            const size_t group) {
      for (const auto& component : tensor) {
        const auto modes = polynomial_power(component, mesh);
        for (size_t d = 0; d < 3; ++d) {
          power[group][d] += square(modes[d]);
        }
      }
    };
    accumulate(metric, 0);
    accumulate(pi, 1);
    accumulate(phi, 1);
    for (size_t group = 0; group < 2; ++group) {
      for (size_t d = 0; d < 3; ++d) {
        power[group][d] = sqrt(power[group][d] / (group == 0 ? 10.0 : 40.0));
      }
    }
  }
  for (size_t group = 0; group < 2; ++group) {
    for (size_t v = 0; v < result.refinement_variables.size(); ++v) {
      result.spectra[group][v] = diagnose_spec_power(
          power[group][v], filter.current[v], filter.after_p_down[v]);
    }
  }

  const auto [c2, c3] = gh_constraints(metric, pi, phi, mesh, inverse_jacobian,
                                       coordinates, time, gauge, gamma2);
  const auto jacobian = determinant_and_inverse(inverse_jacobian).second;
  const auto inverse_grid = determinant_and_inverse(grid_jacobian).second;
  const size_t n = mesh.number_of_grid_points();
  Jacobian<DataVector, 3, Frame::ElementLogical, Frame::Grid> element_grid(n,
                                                                           0.0);
  for (size_t j = 0; j < 3; ++j) {
    for (size_t d = 0; d < 3; ++d) {
      for (size_t k = 0; k < 3; ++k) {
        element_grid.get(j, d) += inverse_grid.get(j, k) * jacobian.get(k, d);
      }
    }
  }
  tnsr::i<DataVector, 3, Frame::ElementLogical> normalization(n, 0.0),
      projected{};
  for (size_t d = 0; d < 3; ++d) {
    for (size_t j = 0; j < 3; ++j) {
      for (size_t k = 0; k < 3; ++k) {
        // Keep the product factors separate, as in SpEC, rather than
        // squaring the composed Jacobian (which introduces cross terms).
        normalization.get(d) +=
            square(grid_jacobian.get(k, j) * element_grid.get(j, d));
      }
    }
  }
  DataVector buffer{};
  const auto project = [&](const auto& constraint, const size_t group) {
    amr::Criteria::Constraints_detail::logical_constraints(
        make_not_null(&projected), make_not_null(&buffer), constraint, jacobian,
        normalization);
    std::array<double, 3> rms{};
    for (size_t d = 0; d < 3; ++d) {
      rms[d] = sqrt(sum(square(projected.get(d))) / static_cast<double>(n));
    }
    result.constraints[group] =
        aggregate_constraints(rms, result.refinement_variables);
  };
  project(c2, 0);
  project(c3, 1);
  return result;
}

std::array<double, 3> updated_target_shifts(const History& history,
                                            const Diagnostics& diagnostics,
                                            const Settings& settings) {
  ASSERT(std::isfinite(settings.adjustment_limit) and
             settings.adjustment_limit <= 0.0,
         "The target adjustment limit must be finite and nonpositive");
  std::array<double, 3> result{};
  if (settings.adjustment_limit == 0.0) {
    return result;
  }
  for (size_t d = 0; d < diagnostics.refinement_variables.size(); ++d) {
    const double shift = history.log10_target_shift[d];
    result[d] = shift;
    double convergence = std::numeric_limits<double>::infinity();
    double error = 0.0;
    for (const auto& group : diagnostics.spectra) {
      if (not group[d].degenerate) {
        convergence = std::min(convergence, group[d].convergence);
        error = std::max(error, group[d].error);
      }
    }
    // A negative slope cannot define a positive number of spectral modes to
    // tighten by. Retain the previous target until a usable fit is available.
    if (not std::isfinite(convergence) or convergence <= 1.e-10) {
      continue;
    }
    const double difference =
        error > 0.0 ? log10(error) - log10(settings.spectral_target)
                    : -std::numeric_limits<double>::infinity();
    double most_restrictive = std::numeric_limits<double>::infinity();
    for (size_t c = 0; c < diagnostics.constraints.size(); ++c) {
      const double current = diagnostics.constraints[c][d];
      const double previous = history.constraints[c][d];
      const bool have_previous = history.calls > 0;
      // Retain the previous constraint even when a p change was realized.
      const double predicted =
          have_previous ? predict(current, previous) : current;
      double next_shift = shift;
      if (predicted > settings.constraint_targets[c]) {
        if (shift == 0.0 and difference > 0.0) {
          next_shift = std::min(0.0, difference - 0.8 * convergence);
        } else if (have_previous and current >= previous and
                   shift > settings.adjustment_limit) {
          next_shift = std::min(shift - 0.5 * convergence,
                                difference - 0.8 * convergence);
        } else if (shift > settings.adjustment_limit) {
          next_shift = shift - 0.5 * convergence;
        }
      } else if (have_previous and current < previous and
                 std::max(current, predicted) <
                     0.5 * settings.constraint_targets[c] and
                 shift < 0.0) {
        next_shift = shift + 0.1 * convergence;
      }
      most_restrictive = std::min(most_restrictive, next_shift);
    }
    result[d] = std::clamp(most_restrictive, settings.adjustment_limit, 0.0);
  }
  return result;
}

namespace {
std::array<amr::Flag, 3> decide_impl(
    const gsl::not_null<History*> history, const Diagnostics& diagnostics,
    const std::array<size_t, 3>& extents, const double time,
    const std::array<amr::Flag, 3>& logical_flags, const Settings& settings,
    const bool standalone) {
  const auto& variables = diagnostics.refinement_variables;
  if (history->calls == 0 and history->initial_meshes.empty()) {
    history->refinement_variables = variables;
  }
  ASSERT(history->refinement_variables == variables,
         "Cannot change independent refinement variables without resetting "
         "AMR history");
  for (size_t d = variables.size(); d < 3; ++d) {
    ASSERT(extents[d] == 0, "Unused refinement extents must be zero");
  }
  auto flags = variables.independent_flags(logical_flags);
  // Repeated observations/evaluations at a fixed time must be idempotent.
  const bool initializing = history->initializing;
  if (initializing) {
    if (history->initial_cycle_detected) {
      return make_array<3>(amr::Flag::DoNothing);
    }
    if (not history->initial_meshes.empty() and
        history->initial_meshes.back() == extents) {
      return history->decision;
    }
    if (std::find(history->initial_meshes.begin(),
                  history->initial_meshes.end(),
                  extents) != history->initial_meshes.end()) {
      history->initial_cycle_detected = true;
      history->decision = make_array<3>(amr::Flag::DoNothing);
      return history->decision;
    }
    history->initial_meshes.push_back(extents);
  }
  if (not initializing and history->calls > 0 and time == history->time and
      extents == history->extents) {
    return history->decision;
  }
  ASSERT(time >= history->time, "AMR history ran backwards");
  const size_t call = history->calls;
  const bool joining =
      std::any_of(flags.begin(), flags.end(),
                  [](const auto f) { return f == amr::Flag::Join; });
  bool join_allowed = true;
  for (size_t group = 0; group < 2; ++group) {
    for (size_t d = 0; d < variables.size(); ++d) {
      const double c = diagnostics.constraints[group][d];
      if (not std::isfinite(c)) {
        ERROR("Nonfinite projected GH constraint");
      }
      join_allowed &= 1.25 * c <= settings.constraint_targets[group];
    }
  }
  history->log10_target_shift =
      updated_target_shifts(*history, diagnostics, settings);
  for (size_t d = 0; d < variables.size(); ++d) {
    const double target =
        settings.spectral_target * pow(10.0, history->log10_target_shift[d]);
    if (not initializing and call > 0 and extents[d] != history->extents[d]) {
      // Change took place following the previous evaluation, which is when
      // SpEC starts its cooldown. Use realized extents, not requested flags.
      history->last_change[d] = call - 1;
      history->last_change_time[d] = history->time;
    }
    if (flags[d] == amr::Flag::Split) {
      continue;
    }
    if (flags[d] == amr::Flag::Join) {
      if (not join_allowed) {
        flags[d] = amr::Flag::DoNothing;
      }
      continue;
    }
    if (joining) {
      // The sibling criterion checked a parent at the current child extents.
      // Simultaneous p-down would invalidate that prospective-parent check.
      continue;
    }
    if (not initializing and
        ((call > 0 and call < history->last_change[d] + settings.cooldown) or
         time < history->last_change_time[d] +
                    settings.minimum_time_between_changes)) {
      flags[d] = amr::Flag::DoNothing;
      continue;
    }
    double error = 0.0, next = 0.0;
    double pile_up = std::numeric_limits<double>::infinity();
    bool has_spectrum = false;
    bool constraints_ok = true;
    bool constraint_excess = false;
    // The derivative spectrum predicts the two-index constraint and the
    // metric spectrum predicts the three-index constraint, as in SpEC.
    for (size_t group = 0; group < 2; ++group) {
      const auto& spectrum = diagnostics.spectra[group][d];
      error = std::max(error, spectrum.error);
      next = std::max(next, spectrum.next_error);
      if (not spectrum.degenerate) {
        has_spectrum = true;
        pile_up = std::min(pile_up, spectrum.pile_up);
      }
      const size_t constraint = 1 - group;
      const double current = diagnostics.constraints[constraint][d];
      const double predicted =
          call > 0 ? predict(current, history->constraints[constraint][d])
                   : current;
      constraint_excess |= predicted >= settings.constraint_targets[constraint];
      const double ratio =
          spectrum.error > 0.0
              ? std::max(1.0, spectrum.next_error / spectrum.error)
              : (spectrum.next_error > 0.0
                     ? std::numeric_limits<double>::infinity()
                     : 1.0);
      const double predicted_next = predicted == 0.0 ? 0.0 : predicted * ratio;
      constraints_ok &=
          1.25 * predicted_next <= settings.constraint_targets[constraint];
    }
    if (not has_spectrum) {
      pile_up = 0.0;
    }
    const bool ordinary =
        pow(10.0, 0.1) * next <= settings.spectral_target and constraints_ok;
    const bool forced = settings.force_pile_up and pile_up > 3.0 and
                        (pile_up - 1.1 > 3.0 or next < target);
    const bool excessive = pile_up > 3.0;
    if (standalone or flags[d] == amr::Flag::IncreaseResolution) {
      // Forced pile-up has priority over p-addition. Otherwise apply SpEC's
      // pile-up hysteresis before adding points.
      if (not forced) {
        // Retain SpEC's base-target/constraint permission gate, and also
        // require LargeTrunc against the shifted target. Only the legacy
        // unadjusted hybrid relies on Chalmers to supply that accuracy check.
        if (pile_up <= (settings.force_pile_up ? 1.5 : 3.0) and
            (error >= settings.spectral_target or constraint_excess) and
            ((not standalone and settings.adjustment_limit == 0.0) or
             error > target)) {
          flags[d] = amr::Flag::IncreaseResolution;
          continue;
        }
        flags[d] = amr::Flag::DoNothing;
      }
    }
    const bool small = error < target and next < target;
    const bool pile_drop =
        excessive and ordinary and (pile_up - 1.1 > 3.0 or next < target);
    flags[d] = (initializing or call > 0) and
                       extents[d] > settings.minimum_modes and
                       (forced or pile_drop or (ordinary and small))
                   ? amr::Flag::DecreaseResolution
                   : amr::Flag::DoNothing;
  }
  if (not initializing) {
    ++history->calls;
    history->time = time;
    history->extents = extents;
    history->constraints = diagnostics.constraints;
  }
  history->decision = variables.logical_flags(flags);
  return history->decision;
}
}  // namespace

std::array<amr::Flag, 3> decide(const gsl::not_null<History*> history,
                                const Diagnostics& diagnostics,
                                const std::array<size_t, 3>& extents,
                                const double time,
                                const std::array<amr::Flag, 3> field_decision,
                                const Settings& settings) {
  return decide_impl(history, diagnostics, extents, time, field_decision,
                     settings, false);
}

std::array<amr::Flag, 3> decide_p(const gsl::not_null<History*> history,
                                  const Diagnostics& diagnostics,
                                  const std::array<size_t, 3>& extents,
                                  const double time, const Settings& settings) {
  return decide_impl(
      history, diagnostics, extents, time,
      {{amr::Flag::DoNothing, amr::Flag::DoNothing, amr::Flag::DoNothing}},
      settings, true);
}
}  // namespace gh::spec_amr
