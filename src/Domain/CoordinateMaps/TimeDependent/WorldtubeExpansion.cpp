// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"

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

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "NumericalAlgorithms/RootFinding/TOMS748.hpp"
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

template <size_t Dim>
WorldtubeExpansion<Dim>::WorldtubeExpansion(double inner_boundary,
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
std::array<tt::remove_cvref_wrap_t<T>, Dim> WorldtubeExpansion<Dim>::operator()(
    const std::array<T, Dim>& source_coords, const double time,
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

template <size_t Dim>
std::optional<std::array<double, Dim>> CubicScale<Dim>::inverse(
    const std::array<double, Dim>& target_coords, const double time,
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
  if (r - f > outer_boundary_) {
    for (size_t i = 0; i < Dim; ++i) {
      gsl::at(*result, i) -= f * gsl::at(target_coords, i) / r;
    }
  } else if (r > inner_boundary_) {
    const double inner_minus_outer_cubed = cube(r_in_ - r_out_);
    const double cbrt_3 = std::cbrt(3.);
    const double big_cbrt = std::cbrt(
        -9 * cube(r_in_ - r_out_) * square(f) * (r_in_ + r_out_ + f - 2 * r) +
        sqrt(3.) *
            sqrt(square(inner_minus_outer_cubed) * cube(f) *
                 (8 * cube(r_in_) - 3 * square(r_in_) * (8 * r_out_ + 3 * f) -
                  square(r_out_) * (8 * r_out_ + 9 * f) +
                  6 * r_in_ *
                      (4 * square(b) + 21 * r_out_ * f + 18 * f * (f - r)) -
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

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim> CubicScale<Dim>::frame_velocity(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));
  const double dt_a_of_t =
      functions_of_time.at(f_of_t_a_)->func_and_deriv(time)[1][0];

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

template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>
CubicScale<Dim>::jacobian(
    const std::array<T, Dim>& source_coords, const double time,
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
  tt::remove_cvref_wrap_t<T> rho =
      square(dereference_wrapper(source_coords[0]));
  for (size_t i = 1; i < Dim; ++i) {
    rho += square(dereference_wrapper(gsl::at(source_coords, i)));
  }
  rho = sqrt(rho);
  for (const size_t i = 0; i < Dim; ++i) {
    jac.get(i, i) = 1.;
    const double cube_jac = 3. * a_ * rho * source_coords.get(i) +
                            2. * b_ * source_coords.get(i) +
                            c_ * source_coords.get(i) / rho;
    for (const size_t j = 0; j < Dim; ++j) {
      jac.get(j, i) += step_function(r_out_ - rho) *
                       step_function(rho - r_in_) * f * cube_jac;
    }
  }

  return jac;
}

template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>
CubicScale<Dim>::inv_jacobian(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  ASSERT(functions_of_time.find(f_of_t_a_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_a_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));
  ASSERT(functions_of_time.find(f_of_t_b_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_b_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  const double a_of_t = functions_of_time.at(f_of_t_a_)->func(time)[0][0];

  if (functions_of_time_equal_) {
    // optimization for linear radial scaling
    auto inv_jac{make_with_value<
        tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>>(
        dereference_wrapper(source_coords[0]), 0.0)};
    const double one_over_a = 1.0 / a_of_t;
    for (size_t i = 0; i < Dim; ++i) {
      inv_jac.get(i, i) = one_over_a;
    }
    return inv_jac;
  }

  const double b_of_t = functions_of_time.at(f_of_t_b_)->func(time)[0][0];

  tt::remove_cvref_wrap_t<T> rho_squared =
      square(dereference_wrapper(source_coords[0]));
  for (size_t i = 1; i < Dim; ++i) {
    rho_squared += square(dereference_wrapper(gsl::at(source_coords, i)));
  }
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame> inv_jac{};
  get<0, 0>(inv_jac) =
      1.0 / (a_of_t + (b_of_t - a_of_t) * square(one_over_outer_boundary_) *
                          rho_squared);
  for (size_t i = 1; i < Dim; ++i) {
    inv_jac.get(i, i) = get<0, 0>(inv_jac);
  }

  // Factor out `double` computations to ensure minimal DataVector operations
  const double denom_constant_a = a_of_t / square(one_over_outer_boundary_);
  const double denom_constant_b = 3.0 * (b_of_t - a_of_t);
  const double numerator_constant = -2.0 * (b_of_t - a_of_t);
  if (Dim == 1) {
    get<0, 0>(inv_jac) *=
        (1.0 + numerator_constant /
                   (denom_constant_a + denom_constant_b * rho_squared) *
                   square(source_coords[0]));
  } else {
    // Reuse rho^2 allocation
    rho_squared = numerator_constant /
                  (denom_constant_a + denom_constant_b * rho_squared) *
                  get<0, 0>(inv_jac);

    for (size_t i = 0; i < Dim; ++i) {
      for (size_t j = 0; j < Dim; ++j) {
        if (i == j) {
          inv_jac.get(i, j) += rho_squared * gsl::at(source_coords, i) *
                               gsl::at(source_coords, j);
        } else {
          inv_jac.get(i, j) = rho_squared * gsl::at(source_coords, i) *
                              gsl::at(source_coords, j);
        }
      }
    }
  }

  return inv_jac;
}

template <size_t Dim>
void CubicScale<Dim>::pup(PUP::er& p) {
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | f_of_t_a_;
    p | f_of_t_b_;
    p | one_over_outer_boundary_;
    p | functions_of_time_equal_;
  }
}

template <size_t Dim>
bool operator==(const CubicScale<Dim>& lhs, const CubicScale<Dim>& rhs) {
  return lhs.f_of_t_a_ == rhs.f_of_t_a_ and lhs.f_of_t_b_ == rhs.f_of_t_b_ and
         lhs.one_over_outer_boundary_ == rhs.one_over_outer_boundary_ and
         lhs.functions_of_time_equal_ == rhs.functions_of_time_equal_;
}

// Explicit instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                             \
  template class CubicScale<DIM(data)>;                  \
  template bool operator==(const CubicScale<DIM(data)>&, \
                           const CubicScale<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                           \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)> \
  CubicScale<DIM(data)>::operator()(                                   \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)> \
  CubicScale<DIM(data)>::frame_velocity(                               \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  CubicScale<DIM(data)>::jacobian(                                     \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;                                    \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),   \
                    Frame::NoFrame>                                    \
  CubicScale<DIM(data)>::inv_jacobian(                                 \
      const std::array<DTYPE(data), DIM(data)>& source_coords,         \
      const double time,                                               \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&   \
          functions_of_time) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3),
                        (double, DataVector,
                         std::reference_wrapper<const double>,
                         std::reference_wrapper<const DataVector>))
#undef DIM
#undef DTYPE
#undef INSTANTIATE
}  // namespace domain::CoordinateMaps::TimeDependent
