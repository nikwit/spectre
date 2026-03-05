// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <utility>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Kokkos/KokkosCore.hpp"

namespace gh::Actions::detail {

template <std::size_t N, class F, std::size_t... Is>
KOKKOS_INLINE_FUNCTION constexpr void static_for_impl(
    F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t N, class F>
KOKKOS_INLINE_FUNCTION constexpr void static_for(F&& f) {
  static_for_impl<N>(static_cast<F&&>(f), std::make_index_sequence<N>{});
}

template <size_t VolumeDim>
KOKKOS_INLINE_FUNCTION void inverse_spatial_metric_and_det(
    const gsl::not_null<tnsr::II<double, VolumeDim, Frame::Inertial>*>
        inverse_spatial_metric,
    const gsl::not_null<double*> det_spatial_metric,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& spacetime_metric) {
  static_assert(VolumeDim == 3,
                "GH Kokkos batched boundary helpers currently support 3D.");

  const double g00 = spacetime_metric.get(1, 1);
  const double g01 = spacetime_metric.get(1, 2);
  const double g02 = spacetime_metric.get(1, 3);
  const double g11 = spacetime_metric.get(2, 2);
  const double g12 = spacetime_metric.get(2, 3);
  const double g22 = spacetime_metric.get(3, 3);

  *det_spatial_metric = g00 * (g11 * g22 - g12 * g12) -
                        g01 * (g01 * g22 - g12 * g02) +
                        g02 * (g01 * g12 - g11 * g02);
  const double inv_det = 1.0 / *det_spatial_metric;

  inverse_spatial_metric->get(0, 0) = (g11 * g22 - g12 * g12) * inv_det;
  inverse_spatial_metric->get(0, 1) = (g02 * g12 - g01 * g22) * inv_det;
  inverse_spatial_metric->get(0, 2) = (g01 * g12 - g02 * g11) * inv_det;
  inverse_spatial_metric->get(1, 1) = (g00 * g22 - g02 * g02) * inv_det;
  inverse_spatial_metric->get(1, 2) = (g02 * g01 - g00 * g12) * inv_det;
  inverse_spatial_metric->get(2, 2) = (g00 * g11 - g01 * g01) * inv_det;
}

template <size_t VolumeDim>
KOKKOS_INLINE_FUNCTION void compute_packaged_boundary_data_at_point(
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*>
        char_speed_v_spacetime_metric,
    const gsl::not_null<tnsr::iaa<double, VolumeDim, Frame::Inertial>*>
        char_speed_v_zero,
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*>
        char_speed_v_plus,
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*>
        char_speed_v_minus,
    const gsl::not_null<tnsr::iaa<double, VolumeDim, Frame::Inertial>*>
        char_speed_n_times_v_plus,
    const gsl::not_null<tnsr::iaa<double, VolumeDim, Frame::Inertial>*>
        char_speed_n_times_v_minus,
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*>
        char_speed_gamma2_v_spacetime_metric,
    const gsl::not_null<tnsr::a<double, VolumeDim, Frame::Inertial>*>
        char_speeds,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& spacetime_metric,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& pi,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>& phi,
    const double gamma1, const double gamma2,
    const tnsr::i<double, VolumeDim, Frame::Inertial>&
        unnormalized_normal_covector) {
  tnsr::II<double, VolumeDim, Frame::Inertial> inverse_spatial_metric{};
  double det_spatial_metric = 0.0;
  inverse_spatial_metric_and_det<VolumeDim>(
      make_not_null(&inverse_spatial_metric),
      make_not_null(&det_spatial_metric), spacetime_metric);
  (void)det_spatial_metric;

  tnsr::I<double, VolumeDim, Frame::Inertial> shift{};
  for (size_t i = 0; i < VolumeDim; ++i) {
    shift.get(i) = 0.0;
    for (size_t j = 0; j < VolumeDim; ++j) {
      shift.get(i) +=
          inverse_spatial_metric.get(i, j) * spacetime_metric.get(0, j + 1);
    }
  }
  double lapse_squared = -spacetime_metric.get(0, 0);
  for (size_t i = 0; i < VolumeDim; ++i) {
    lapse_squared += shift.get(i) * spacetime_metric.get(0, i + 1);
  }
  const double lapse = sqrt(lapse_squared);

  tnsr::I<double, VolumeDim, Frame::Inertial> normal_vector{};
  double normal_magnitude_squared = 0.0;
  for (size_t i = 0; i < VolumeDim; ++i) {
    normal_vector.get(i) = 0.0;
    for (size_t j = 0; j < VolumeDim; ++j) {
      normal_vector.get(i) += inverse_spatial_metric.get(i, j) *
                              unnormalized_normal_covector.get(j);
    }
    normal_magnitude_squared +=
        normal_vector.get(i) * unnormalized_normal_covector.get(i);
  }
  const double one_over_normal_magnitude = 1.0 / sqrt(normal_magnitude_squared);
  tnsr::i<double, VolumeDim, Frame::Inertial> normal_covector{};
  for (size_t i = 0; i < VolumeDim; ++i) {
    normal_vector.get(i) *= one_over_normal_magnitude;
    normal_covector.get(i) =
        unnormalized_normal_covector.get(i) * one_over_normal_magnitude;
  }

  double shift_dot_normal = 0.0;
  for (size_t i = 0; i < VolumeDim; ++i) {
    shift_dot_normal += shift.get(i) * normal_covector.get(i);
  }
  shift_dot_normal *= -1.0;

  char_speeds->get(0) = (1.0 + gamma1) * shift_dot_normal;
  char_speeds->get(1) = shift_dot_normal;
  char_speeds->get(2) = lapse + shift_dot_normal;
  char_speeds->get(3) = -lapse + shift_dot_normal;

  for (size_t a = 0; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      char_speed_gamma2_v_spacetime_metric->get(a, b) =
          gamma2 * spacetime_metric.get(a, b);
    }
  }

  tnsr::aa<double, VolumeDim, Frame::Inertial> normal_dot_phi{};
  for (size_t a = 0; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      normal_dot_phi.get(a, b) = normal_vector.get(0) * phi.get(0, a, b);
      for (size_t i = 1; i < VolumeDim; ++i) {
        normal_dot_phi.get(a, b) += normal_vector.get(i) * phi.get(i, a, b);
      }
    }
  }

  for (size_t a = 0; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      char_speed_v_plus->get(a, b) =
          char_speeds->get(2) *
          (pi.get(a, b) + normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));
      char_speed_v_minus->get(a, b) =
          char_speeds->get(3) *
          (pi.get(a, b) - normal_dot_phi.get(a, b) -
           char_speed_gamma2_v_spacetime_metric->get(a, b));

      for (size_t i = 0; i < VolumeDim; ++i) {
        char_speed_v_zero->get(i, a, b) =
            char_speeds->get(1) *
            (phi.get(i, a, b) -
             normal_covector.get(i) * normal_dot_phi.get(a, b));
      }
    }
  }

  for (size_t a = 0; a < VolumeDim + 1; ++a) {
    for (size_t b = a; b < VolumeDim + 1; ++b) {
      for (size_t i = 0; i < VolumeDim; ++i) {
        char_speed_n_times_v_plus->get(i, a, b) =
            char_speed_v_plus->get(a, b) * normal_covector.get(i);
        char_speed_n_times_v_minus->get(i, a, b) =
            char_speed_v_minus->get(a, b) * normal_covector.get(i);
      }
      char_speed_v_spacetime_metric->get(a, b) =
          char_speeds->get(0) * spacetime_metric.get(a, b);
      char_speed_gamma2_v_spacetime_metric->get(a, b) *= char_speeds->get(0);
    }
  }
}

KOKKOS_INLINE_FUNCTION double step_function_double(const double value) {
  return value < 0.0 ? 0.0 : 1.0;
}

template <size_t VolumeDim>
KOKKOS_INLINE_FUNCTION void compute_boundary_terms_at_point(
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*>
        dt_spacetime_metric,
    const gsl::not_null<tnsr::aa<double, VolumeDim, Frame::Inertial>*> dt_pi,
    const gsl::not_null<tnsr::iaa<double, VolumeDim, Frame::Inertial>*> dt_phi,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>&
        local_v_spacetime_metric,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>& local_v_zero,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& local_v_plus,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& local_v_minus,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>&
        local_normal_times_v_plus,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>&
        local_normal_times_v_minus,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>&
        local_gamma2_v_spacetime_metric,
    const tnsr::a<double, VolumeDim, Frame::Inertial>& local_char_speeds,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>&
        remote_v_spacetime_metric,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>& remote_v_zero,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& remote_v_plus,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>& remote_v_minus,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>&
        remote_normal_times_v_plus,
    const tnsr::iaa<double, VolumeDim, Frame::Inertial>&
        remote_normal_times_v_minus,
    const tnsr::aa<double, VolumeDim, Frame::Inertial>&
        remote_gamma2_v_spacetime_metric,
    const tnsr::a<double, VolumeDim, Frame::Inertial>& remote_char_speeds) {
  const double weighted_lambda_spacetime_metric_int =
      step_function_double(-get<0>(local_char_speeds));
  const double weighted_lambda_spacetime_metric_ext =
      -step_function_double(get<0>(remote_char_speeds));
  const double weighted_lambda_zero_int =
      step_function_double(-get<1>(local_char_speeds));
  const double weighted_lambda_zero_ext =
      -step_function_double(get<1>(remote_char_speeds));
  const double weighted_lambda_plus_int =
      step_function_double(-get<2>(local_char_speeds));
  const double weighted_lambda_plus_ext =
      -step_function_double(get<2>(remote_char_speeds));
  const double weighted_lambda_minus_int =
      step_function_double(-get<3>(local_char_speeds));
  const double weighted_lambda_minus_ext =
      -step_function_double(get<3>(remote_char_speeds));

  const double half_weighted_lambda_plus_ext = 0.5 * weighted_lambda_plus_ext;
  const double half_weighted_lambda_minus_ext = 0.5 * weighted_lambda_minus_ext;
  const double half_weighted_lambda_plus_int = 0.5 * weighted_lambda_plus_int;
  const double half_weighted_lambda_minus_int = 0.5 * weighted_lambda_minus_int;

  static_for<VolumeDim + 1>([&](auto a_c) {
    constexpr size_t a = static_cast<size_t>(a_c);
    static_for<VolumeDim + 1 - a>([&](auto off_c) {
      constexpr size_t b = a + static_cast<size_t>(off_c);

      get<a, b>(*dt_spacetime_metric) =
          fma(weighted_lambda_spacetime_metric_ext,
              get<a, b>(remote_v_spacetime_metric),
              -weighted_lambda_spacetime_metric_int *
                  get<a, b>(local_v_spacetime_metric));

      const double remote_pi_term =
          fma(half_weighted_lambda_plus_ext, get<a, b>(remote_v_plus),
              half_weighted_lambda_minus_ext * get<a, b>(remote_v_minus));
      const double local_pi_term =
          fma(half_weighted_lambda_plus_int, get<a, b>(local_v_plus),
              half_weighted_lambda_minus_int * get<a, b>(local_v_minus));
      const double remote_gamma_term =
          weighted_lambda_spacetime_metric_ext *
          get<a, b>(remote_gamma2_v_spacetime_metric);
      const double local_gamma_term =
          weighted_lambda_spacetime_metric_int *
          get<a, b>(local_gamma2_v_spacetime_metric);
      get<a, b>(*dt_pi) = (remote_pi_term + remote_gamma_term) -
                          (local_pi_term + local_gamma_term);

      static_for<VolumeDim>([&](auto d_c) {
        constexpr size_t d = static_cast<size_t>(d_c);
        const double remote_flux =
            fma(half_weighted_lambda_plus_ext,
                get<d, a, b>(remote_normal_times_v_plus),
                -half_weighted_lambda_minus_ext *
                    get<d, a, b>(remote_normal_times_v_minus));
        const double local_flux =
            fma(half_weighted_lambda_plus_int,
                get<d, a, b>(local_normal_times_v_plus),
                -half_weighted_lambda_minus_int *
                    get<d, a, b>(local_normal_times_v_minus));
        const double remote_zero =
            weighted_lambda_zero_ext * get<d, a, b>(remote_v_zero);
        const double local_zero =
            weighted_lambda_zero_int * get<d, a, b>(local_v_zero);
        get<d, a, b>(*dt_phi) =
            (remote_flux + remote_zero) - (local_flux + local_zero);
      });
    });
  });
}

}  // namespace gh::Actions::detail
