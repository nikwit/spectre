// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/WorldtubeExpansion.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <ostream>
#include <pup.h>
#include <pup_stl.h>
#include <string>
#include <unordered_map>
#include <utility>

#include "DataStructures/Blaze/StepFunction.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TypeTraits/RemoveReferenceWrapper.hpp"

namespace domain::CoordinateMaps::TimeDependent {

WorldtubeExpansion::WorldtubeExpansion(double inner_boundary,
                                       double outer_boundary,
                                       std::string function_of_time_name)
    : f_of_t_name_(std::move(function_of_time_name)),
      r_in_(inner_boundary),
      r_out_(outer_boundary) {
  const double denom = cube(r_in_ - r_out_);
  a_ = 2. / denom;
  b_ = -3. * (r_in_ + r_out_) / denom;
  c_ = 6. * r_in_ * r_out_ / denom;
  d_ = r_in_ * r_in_ * (r_in_ - 3. * r_out_) / denom;
}

template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, 3> WorldtubeExpansion::operator()(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  const double f_of_t = functions_of_time.at(f_of_t_name_)->func(time)[0][0];

  tt::remove_cvref_wrap_t<T> rho_sq =
      square(dereference_wrapper(source_coords[0]));
  for (size_t i = 1; i < Dim; ++i) {
    rho_sq += square(dereference_wrapper(gsl::at(source_coords, i)));
  }
  const auto rho = sqrt(rho_sq);
  const auto cubic = a_ * rho_sq * rho + b_ * rho_sq + c_ * rho + d_;

  std::array<tt::remove_cvref_wrap_t<T>, Dim> result{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(result, i) = gsl::at(source_coords, i) +
                         (1. + step_function(r_out_ - rho) * (cubic - 1.)) *
                             step_function(rho - r_in_) * f_of_t *
                             gsl::at(source_coords, i) / rho;
  }
  return result;
}

std::optional<std::array<double, 3>> WorldtubeExpansion::inverse(
    const std::array<double, 3>& target_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));
  const double f = functions_of_time.at(f_of_t_name_)->func(time)[0][0];

  std::optional<std::array<double, Dim>> result{std::array<double, Dim>{}};
  const double r = magnitude(target_coords);
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(*result, i) = gsl::at(target_coords, i);
  }
  if (r - f > r_out_) {
    for (size_t i = 0; i < Dim; ++i) {
      gsl::at(*result, i) -= f * gsl::at(target_coords, i) / r;
    }
  } else if (r > r_in_) {
    const double inner_minus_outer_cubed = cube(r_in_ - r_out_);
    const double cbrt_3 = std::cbrt(3.);
    const double big_cbrt = std::cbrt(
        -9 * cube(r_in_ - r_out_) * square(f) * (r_in_ + r_out_ + f - 2 * r) +
        sqrt(3.) *
            sqrt(
                square(inner_minus_outer_cubed) * cube(f) *
                (8 * cube(r_in_) - 3 * square(r_in_) * (8 * r_out_ + 3 * f) -
                 square(r_out_) * (8 * r_out_ + 9 * f) +
                 6 * r_in_ *
                     (4 * square(r_out_) + 21 * r_out_ * f + 18 * f * (f - r)) -
                 108 * f * (r_out_ + f) * r + 108 * f * square(r))));
    const double rho = (18 * (r_in_ + r_out_) * f -
                        (6 * square(cbrt_3) * square(r_in_ - r_out_) *
                         (2 * r_in_ - 2 * r_out_ - 3 * f) * f) /
                            big_cbrt +
                        6 * cbrt_3 * big_cbrt) /
                       (36. * f);
    for (size_t i = 0; i < Dim; ++i) {
      gsl::at(*result, i) = rho * gsl::at(target_coords, i) / r;
    }
  }
  return result;
}

template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, 3> WorldtubeExpansion::frame_velocity(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));
  const double dt_a_of_t =
      functions_of_time.at(f_of_t_name_)->func_and_deriv(time)[1][0];

  tt::remove_cvref_wrap_t<T> rho_sq =
      square(dereference_wrapper(source_coords[0]));
  for (size_t i = 1; i < Dim; ++i) {
    rho_sq += square(dereference_wrapper(gsl::at(source_coords, i)));
  }
  const auto rho = sqrt(rho_sq);
  const auto cubic = a_ * rho_sq * rho + b_ * rho_sq + c_ * rho + d_;

  std::array<tt::remove_cvref_wrap_t<T>, Dim> result{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(result, i) = (1. + step_function(r_out_ - rho) * (cubic - 1.)) *
                         step_function(rho - r_in_) * dt_a_of_t *
                         gsl::at(source_coords, i) / rho;
  }
  return result;
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
WorldtubeExpansion::jacobian(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  const double f_of_t = functions_of_time.at(f_of_t_name_)->func(time)[0][0];
  auto jac{make_with_value<
      tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0)};
  tt::remove_cvref_wrap_t<T> rho_sq =
      square(dereference_wrapper(source_coords[0]));
  for (size_t i = 1; i < 3; ++i) {
    rho_sq += square(dereference_wrapper(gsl::at(source_coords, i)));
  }
  const auto rho = sqrt(rho_sq);
  const auto in_outer_region = step_function(rho - r_out_);
  const auto in_trans_region =
      step_function(r_out_ - rho) * step_function(rho - r_in_);

  for (size_t i = 0; i < 3; ++i) {
    jac.get(i, i) = 1.;
    jac.get(i, i) += in_outer_region * (1 + f_of_t / rho);
    jac.get(i, i) +=
        in_trans_region * f_of_t * (a_ * rho_sq + b_ * rho + c_ + d_ / rho);
    for (size_t j = 0; j < 3; ++j) {
      jac.get(i, j) -= in_outer_region * source_coords.at(i) *
                       source_coords.at(j) * f_of_t / (rho_sq * rho);
      jac.get(i, j) -= in_trans_region * source_coords.at(i) *
                       source_coords.at(j) * f_of_t *
                       (2. * a_ + b_ / rho - d_ / (rho_sq * rho));
    }
  }

  return jac;
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
WorldtubeExpansion::inv_jacobian(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  return determinant_and_inverse(
             jacobian(source_coords, time, functions_of_time))
      .second;
}

void WorldtubeExpansion::pup(PUP::er& p) {
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | f_of_t_name_;
    p | r_in_;
    p | r_out_;
    p | a_;
    p | b_;
    p | c_;
    p | d_;
  }
}

bool operator==(const WorldtubeExpansion& lhs, const WorldtubeExpansion& rhs) {
  return lhs.r_in_ == rhs.r_in_ and lhs.r_out_ == rhs.r_out_;
}

// Explicit instantiations
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                 \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, 3>               \
  WorldtubeExpansion::operator()(                                            \
      const std::array<DTYPE(data), 3>& source_coords, const double time,    \
      const std::unordered_map<                                              \
          std::string,                                                       \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&         \
          functions_of_time) const;                                          \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, 3>               \
  WorldtubeExpansion::frame_velocity(                                        \
      const std::array<DTYPE(data), 3>& source_coords, const double time,    \
      const std::unordered_map<                                              \
          std::string,                                                       \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&         \
          functions_of_time) const;                                          \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame> \
  WorldtubeExpansion::jacobian(                                              \
      const std::array<DTYPE(data), 3>& source_coords, const double time,    \
      const std::unordered_map<                                              \
          std::string,                                                       \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&         \
          functions_of_time) const;                                          \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, 3, Frame::NoFrame> \
  WorldtubeExpansion::inv_jacobian(                                          \
      const std::array<DTYPE(data), 3>& source_coords, const double time,    \
      const std::unordered_map<                                              \
          std::string,                                                       \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&         \
          functions_of_time) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector,
                                      std::reference_wrapper<const double>,
                                      std::reference_wrapper<const DataVector>))
#undef DIM
#undef DTYPE
#undef INSTANTIATE
}  // namespace domain::CoordinateMaps::TimeDependent
