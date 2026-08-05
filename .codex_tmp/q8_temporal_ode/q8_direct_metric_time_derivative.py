#!/usr/bin/env python3
"""Recover first-order coefficient rates from the evolved metric derivative.

This script uses Pi, Phi, lapse, and shift to construct the simulation-time
derivative of the inverse spacetime metric at each downloaded q=8 slice.  It
adds the velocity of the curvature-centered virtual worldtube, fits the same
V1+C3 first-order response blocks to the derivative, and compares the result
with:

1. finite differences of the independently fitted first-order coefficients;
2. the rates freely fitted from the second-order metric residual.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/q8-direct-time-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/q8-direct-time-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

HERE = Path(__file__).resolve().parent
ANALYSIS_PATH = HERE / "q8_temporal_ode_consistency.py"
spec = importlib.util.spec_from_file_location(
    "q8_temporal_ode", ANALYSIS_PATH)
ode = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(ode)

PREVIOUS_RESULTS = HERE / "q8_temporal_ode_consistency.npz"
OUTPUT_FILE = HERE / "q8_direct_metric_time_derivative.npz"

PI_COMPONENT_NAMES = tuple(
    f"Pi_{suffix}" for suffix in ode.COVARIANT_COMPONENT_SUFFIXES)
SHIFT_COMPONENT_NAMES = ("Shift_x", "Shift_y", "Shift_z")
ALL_GH_COMPONENT_NAMES = (
    ode.INVERSE_METRIC_COMPONENT_NAMES
    + ode.PHI_COMPONENT_NAMES
    + PI_COMPONENT_NAMES
    + ("Lapse",)
    + SHIFT_COMPONENT_NAMES
)

VECTOR_CLOSURE_BLOCKS = (("GTI", 2), ("GIJ", 1))
CLOCK_CLOSURE_BLOCKS = (("GTI", 1),)


def covariant_symmetric_matrix(component_values):
    """Unpack ten symmetric covariant components into 4x4 matrices."""
    component_values = np.asarray(component_values)
    result = np.zeros(component_values.shape[1:] + (4, 4))
    for values, (a, b) in zip(
            component_values, ode.COVARIANT_DATA_COMPONENTS):
        result[..., a, b] = values
        result[..., b, a] = values
    return result


def local_polynomial_derivative(
        sample_times, sample_values, target_times,
        number_of_points=41, degree=3):
    """Differentiate nonuniform high-cadence data with a local polynomial."""
    sample_times = np.asarray(sample_times)
    sample_values = np.asarray(sample_values)
    result = np.empty((len(target_times),) + sample_values.shape[1:])
    for target_index, target in enumerate(target_times):
        nearest = np.argsort(
            np.abs(sample_times - target))[:number_of_points]
        nearest = nearest[np.argsort(sample_times[nearest])]
        offsets = sample_times[nearest] - target
        design = np.column_stack([
            offsets**power for power in range(degree + 1)
        ])
        coefficients, *_ = np.linalg.lstsq(
            design, sample_values[nearest], rcond=None)
        result[target_index] = coefficients[1]
    return result


def control_system_center_data():
    with ode.spectre_h5.H5File(
            ode.REDUCTIONS_FILE, "r") as h5file:
        dataset = h5file.get_dat(
            "/ApparentHorizons/ControlSystemAhB_Centers")
        data = np.asarray(dataset.get_data())
    return data[:, 0], data[:, 4:7]


def center_velocities(times, ah_centers, curvature_centers):
    control_times, control_centers = control_system_center_data()
    ah_velocity = local_polynomial_derivative(
        control_times, control_centers, times)
    offset = curvature_centers - ah_centers
    offset_velocity = savgol_filter(
        offset,
        window_length=7,
        polyorder=3,
        deriv=1,
        delta=float(np.median(np.diff(times))),
        axis=0,
        mode="interp",
    )
    curvature_velocity = ah_velocity + offset_velocity
    direct_low_cadence_velocity = savgol_filter(
        curvature_centers,
        window_length=7,
        polyorder=3,
        deriv=1,
        delta=float(np.median(np.diff(times))),
        axis=0,
        mode="interp",
    )
    return curvature_velocity, ah_velocity, direct_low_cadence_velocity


def unpack_gh_fields(values):
    """Unpack one batched interpolation of all required GH fields."""
    number_of_radii = len(ode.RADII)
    number_of_points = len(ode.DIRECTIONS)
    start = 0

    inverse_components = values[
        start:start + len(ode.INVERSE_METRIC_COMPONENT_NAMES)]
    start += len(ode.INVERSE_METRIC_COMPONENT_NAMES)
    inverse_components = inverse_components.reshape(
        len(ode.INVERSE_METRIC_COMPONENT_NAMES),
        number_of_radii,
        number_of_points,
    )
    inverse_metric = np.asarray([
        ode.inverse_metric_matrix(
            inverse_components[:, radius_index])
        for radius_index in range(number_of_radii)
    ])

    phi_components = values[
        start:start + len(ode.PHI_COMPONENT_NAMES)]
    start += len(ode.PHI_COMPONENT_NAMES)
    phi_components = phi_components.reshape(
        3,
        len(ode.COVARIANT_DATA_COMPONENTS),
        number_of_radii,
        number_of_points,
    )
    phi = np.zeros(
        (number_of_radii, number_of_points, 3, 4, 4))
    for derivative in range(3):
        phi[:, :, derivative] = covariant_symmetric_matrix(
            phi_components[derivative])

    pi_components = values[
        start:start + len(PI_COMPONENT_NAMES)]
    start += len(PI_COMPONENT_NAMES)
    pi_components = pi_components.reshape(
        len(PI_COMPONENT_NAMES), number_of_radii, number_of_points)
    pi = covariant_symmetric_matrix(pi_components)

    lapse = values[start].reshape(number_of_radii, number_of_points)
    start += 1
    shift = np.moveaxis(
        values[start:start + 3].reshape(
            3, number_of_radii, number_of_points),
        0,
        -1,
    )
    return inverse_metric, phi, pi, lapse, shift


def interpolate_gh_fields(observation_id, center):
    points = np.concatenate([
        center[None, :] + radius * ode.DIRECTIONS
        for radius in ode.RADII
    ])
    values = ode.interpolate_components(
        observation_id, points, ALL_GH_COMPONENT_NAMES)
    return unpack_gh_fields(values)


def inverse_metric_derivatives(
        inverse_metric, phi, pi, lapse, shift, center_velocity):
    """Construct fixed-inertial and co-moving derivatives of g^{ab}."""
    # Pi_ab = -(partial_T g_ab - beta^i Phi_iab) / alpha.
    covariant_time_derivative = (
        np.einsum("rpi,rpiab->rpab", shift, phi)
        - lapse[..., None, None] * pi
    )
    inverse_fixed_time_derivative = -np.einsum(
        "rpac,rpbd,rpcd->rpab",
        inverse_metric,
        inverse_metric,
        covariant_time_derivative,
    )
    inverse_spatial_derivative = -np.einsum(
        "rpac,rpbd,rpicd->rpiab",
        inverse_metric,
        inverse_metric,
        phi,
    )
    inverse_comoving_time_derivative = (
        inverse_fixed_time_derivative
        + np.einsum(
            "i,rpiab->rpab",
            center_velocity,
            inverse_spatial_derivative,
        )
    )
    return (
        inverse_fixed_time_derivative,
        inverse_comoving_time_derivative,
        inverse_spatial_derivative,
    )


def fit_rate_from_metric_derivative(
        metric_derivative, radius_index):
    """Fit dp_(1)/dT from the co-moving inverse-metric derivative."""
    result = np.zeros(13)
    for blocks, indices in (
        (ode.VECTOR_BLOCKS, ode.VECTOR_PARAMETER_INDICES),
        (ode.CLOCK_STRAIN_BLOCKS,
         ode.CLOCK_STRAIN_PARAMETER_INDICES),
    ):
        matrix = ode.projected_design(
            ode.FIRST_ORDER_COLUMNS[radius_index],
            blocks,
            indices,
        )
        data = ode.projected_data(metric_derivative, blocks)
        estimate, rank, _ = ode.scaled_least_squares(matrix, data)
        if rank != len(indices):
            raise RuntimeError("Metric-time-derivative block lost rank")
        result[indices] = estimate
    return result


def closure_fraction(
        metric_derivative, estimate, radius_index, blocks):
    prediction = ode.metric_linear_combination(
        ode.FIRST_ORDER_COLUMNS[radius_index], estimate)
    data = ode.projected_data(metric_derivative, blocks)
    predicted = ode.projected_data(prediction, blocks)
    return np.linalg.norm(data - predicted) / max(
        np.linalg.norm(data), 1.0e-300)


def symmetric_matrix_series(values):
    return np.asarray([
        ode.symmetric_matrix_from_six(row) for row in values])


def symmetric_rate_from_second_order(second_order):
    return np.einsum(
        "...a,aij->...ij",
        second_order[..., :6],
        ode.SYMMETRIC_BASIS,
    )


def group_relative_mismatch(left, right, axis):
    return np.sqrt(np.sum((left - right)**2, axis=axis)) / np.maximum(
        np.maximum(
            np.sqrt(np.sum(left**2, axis=axis)),
            np.sqrt(np.sum(right**2, axis=axis)),
        ),
        1.0e-300,
    )


def group_alignment(left, right, axis):
    numerator = np.sum(left * right, axis=axis)
    denominator = (
        np.sqrt(np.sum(left**2, axis=axis))
        * np.sqrt(np.sum(right**2, axis=axis))
    )
    return numerator / np.maximum(denominator, 1.0e-300)


def make_plots(
        times,
        first_order,
        second_order,
        direct_harmonic_rates,
        finite_difference_harmonic_rates,
        derivative_validation,
        vector_closure,
        clock_closure,
        stage=3,
):
    central = slice(2, -2)
    radius_scale = ode.RADII / ode.SMALL_BLACK_HOLE_MASS

    direct_clock = direct_harmonic_rates[
        central, :, stage, 0]
    finite_clock = finite_difference_harmonic_rates[
        central, :, stage, 0]
    free_clock = second_order[central, :, stage, -1]
    direct_strain = (
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(
            direct_harmonic_rates[
                central, :, stage, 7:13].reshape(-1, 6)
        ).reshape(
            len(times[central]), len(ode.RADII), 3, 3)
    )
    finite_strain = (
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(
            finite_difference_harmonic_rates[
                central, :, stage, 7:13].reshape(-1, 6)
        ).reshape(
            len(times[central]), len(ode.RADII), 3, 3)
    )
    free_strain = symmetric_rate_from_second_order(
        second_order[central, :, stage])

    fig, axes = plt.subplots(
        2, 2, figsize=(14.0, 9.0), constrained_layout=True)
    axes[0, 0].loglog(
        radius_scale,
        np.sqrt(np.mean(direct_clock**2, axis=0)),
        "o-",
        label=r"direct from $\Pi,\Phi$",
    )
    axes[0, 0].loglog(
        radius_scale,
        np.sqrt(np.mean(finite_clock**2, axis=0)),
        "s--",
        label="finite difference of first-order fit",
    )
    axes[0, 0].loglog(
        radius_scale,
        np.sqrt(np.mean(free_clock**2, axis=0)),
        "^:",
        label="free second-order residual fit",
    )
    axes[0, 0].set_title(r"RMS time acceleration $|\ddot q^0|$")
    axes[0, 0].legend(fontsize=9)

    axes[0, 1].loglog(
        radius_scale,
        np.sqrt(np.mean(
            np.sum(direct_strain**2, axis=(-2, -1)), axis=0)),
        "o-",
        label=r"direct $M_B\dot\Lambda$",
    )
    axes[0, 1].loglog(
        radius_scale,
        np.sqrt(np.mean(
            np.sum(finite_strain**2, axis=(-2, -1)), axis=0)),
        "s--",
        label=r"finite-difference $M_B\dot\Lambda$",
    )
    axes[0, 1].loglog(
        radius_scale,
        np.sqrt(np.mean(
            np.sum(free_strain**2, axis=(-2, -1)), axis=0)),
        "^:",
        label=r"free $\dot L$ residual fit",
    )
    axes[0, 1].set_title("RMS symmetric strain-rate norm")
    axes[0, 1].legend(fontsize=9)

    group_specs = (
        ("clock", np.asarray([0])),
        ("time gradient", np.arange(1, 4)),
        ("spatial acceleration", np.arange(4, 7)),
        ("strain", np.arange(7, 13)),
    )
    for name, indices in group_specs:
        mismatch = np.asarray([
            np.linalg.norm(
                direct_harmonic_rates[
                    central, radius_index, stage][:, indices]
                - finite_difference_harmonic_rates[
                    central, radius_index, stage][:, indices]
            )
            / max(
                np.linalg.norm(
                    direct_harmonic_rates[
                        central, radius_index, stage][:, indices]),
                np.linalg.norm(
                    finite_difference_harmonic_rates[
                        central, radius_index, stage][:, indices]),
                1.0e-300,
            )
            for radius_index in range(len(ode.RADII))
        ])
        axes[1, 0].semilogy(
            radius_scale, mismatch, "o-", label=name)
    axes[1, 0].set_title(
        "Direct versus finite-difference coefficient rates")
    axes[1, 0].set_ylabel("symmetric relative mismatch")
    axes[1, 0].legend(fontsize=9)

    axes[1, 1].semilogy(
        radius_scale,
        np.mean(derivative_validation[central], axis=0),
        "o-",
        label=r"$\Pi,\Phi$ derivative vs metric time series",
    )
    axes[1, 1].semilogy(
        radius_scale,
        np.mean(vector_closure[central], axis=0),
        "s--",
        label="held-out vector closure",
    )
    axes[1, 1].semilogy(
        radius_scale,
        np.mean(clock_closure[central], axis=0),
        "^:",
        label="held-out clock/strain closure",
    )
    axes[1, 1].set_title("Derivative construction and block closure")
    axes[1, 1].set_ylabel("relative residual")
    axes[1, 1].legend(fontsize=9)

    for axis in axes.flat:
        axis.set_xlabel(r"$R/M_B$")
    fig.suptitle(
        f"Per-slice inverse-metric time-derivative matching, "
        f"alternation {stage}")
    fig.savefig(HERE / "q8_direct_metric_time_derivative.png", dpi=180)
    fig.savefig(HERE / "q8_direct_metric_time_derivative.pdf")
    plt.close(fig)

    # Representative time series at the radius minimizing the sum of direct
    # versus finite clock and strain mismatch.
    combined = np.empty(len(ode.RADII))
    for radius_index in range(len(ode.RADII)):
        combined[radius_index] = (
            group_relative_mismatch(
                direct_clock[:, radius_index],
                finite_clock[:, radius_index],
                axis=0,
            )
            + group_relative_mismatch(
                direct_strain[:, radius_index],
                finite_strain[:, radius_index],
                axis=(0, 1, 2),
            )
        )
    best_radius = int(np.argmin(combined))

    fig, axes = plt.subplots(
        2, 2, figsize=(15.0, 9.0), constrained_layout=True)
    axes[0, 0].semilogy(
        times,
        np.abs(direct_harmonic_rates[:, best_radius, stage, 0]),
        "o-",
        label=r"direct from $\Pi,\Phi$",
    )
    axes[0, 0].semilogy(
        times,
        np.abs(finite_difference_harmonic_rates[
            :, best_radius, stage, 0]),
        "s--",
        label="finite difference",
    )
    axes[0, 0].semilogy(
        times,
        np.abs(second_order[:, best_radius, stage, -1]),
        "^:",
        label="free second-order fit",
    )
    axes[0, 0].set_title(r"Time-acceleration magnitude $|\ddot q^0|$")
    axes[0, 0].legend(fontsize=9)

    direct_norm = np.linalg.norm(
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(
            direct_harmonic_rates[
                :, best_radius, stage, 7:13]),
        axis=(1, 2),
    )
    finite_norm = np.linalg.norm(
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(
            finite_difference_harmonic_rates[
                :, best_radius, stage, 7:13]),
        axis=(1, 2),
    )
    free_norm = np.linalg.norm(
        symmetric_rate_from_second_order(
            second_order[:, best_radius, stage]),
        axis=(1, 2),
    )
    axes[0, 1].semilogy(
        times, direct_norm, "o-", label=r"direct $M_B\dot\Lambda$")
    axes[0, 1].semilogy(
        times, finite_norm, "s--", label="finite difference")
    axes[0, 1].semilogy(
        times, free_norm, "^:", label="free second-order fit")
    axes[0, 1].set_title("Symmetric strain-rate norm")
    axes[0, 1].legend(fontsize=9)

    for name, indices in group_specs:
        axes[1, 0].semilogy(
            times,
            np.linalg.norm(
                direct_harmonic_rates[
                    :, best_radius, stage][:, indices],
                axis=1,
            ),
            "o-",
            label=name,
        )
    axes[1, 0].set_title("Norms of all 13 direct harmonic-time rates")
    axes[1, 0].legend(fontsize=9)

    axes[1, 1].semilogy(
        times,
        derivative_validation[:, best_radius],
        "o-",
        label=r"$\Pi,\Phi$ vs time-series metric",
    )
    axes[1, 1].semilogy(
        times,
        vector_closure[:, best_radius],
        "s--",
        label="vector closure",
    )
    axes[1, 1].semilogy(
        times,
        clock_closure[:, best_radius],
        "^:",
        label="clock/strain closure",
    )
    axes[1, 1].set_title("Per-slice validation")
    axes[1, 1].legend(fontsize=9)

    for axis in axes.flat:
        axis.set_xlabel(r"$T/M_{\rm tot}$")
    fig.suptitle(
        f"Representative direct time-derivative fit: "
        f"R={ode.RADII[best_radius]:.2f}M_tot "
        f"({radius_scale[best_radius]:.2f}M_B)")
    fig.savefig(
        HERE / "q8_direct_metric_time_derivative_timeseries.png",
        dpi=180,
    )
    fig.savefig(
        HERE / "q8_direct_metric_time_derivative_timeseries.pdf")
    plt.close(fig)
    return best_radius


def main():
    previous = np.load(PREVIOUS_RESULTS)
    times = previous["times"]
    observations = [
        item for item in ode.volume_observations()
        if times[0] - 1.0 <= item[0] <= times[-1] + 1.0
    ]
    observation_ids = [item[1] for item in observations]
    ah_centers = previous["ah_centers"]
    curvature_centers = previous["curvature_centers"]
    first_order = previous["first_order_history"]
    second_order = previous["second_order_history"]

    (
        curvature_velocity,
        ah_velocity,
        low_cadence_curvature_velocity,
    ) = center_velocities(times, ah_centers, curvature_centers)
    print(
        "center-velocity high/low-cadence maximum difference:",
        np.max(np.linalg.norm(
            curvature_velocity - low_cadence_curvature_velocity,
            axis=1,
        )),
    )

    inverse_metrics = []
    fixed_derivatives = []
    comoving_derivatives = []
    base_rate_estimates = []
    vector_closure = []
    clock_closure = []

    for time_index, (time, observation_id) in enumerate(observations):
        (
            inverse_metric,
            phi,
            pi,
            lapse,
            shift,
        ) = interpolate_gh_fields(
            observation_id, curvature_centers[time_index])
        (
            fixed_derivative,
            comoving_derivative,
            _,
        ) = inverse_metric_derivatives(
            inverse_metric,
            phi,
            pi,
            lapse,
            shift,
            curvature_velocity[time_index],
        )
        rates_at_radii = []
        vector_closure_at_radii = []
        clock_closure_at_radii = []
        for radius_index in range(len(ode.RADII)):
            derivative_metric = ode.metric_tuple_from_matrix(
                comoving_derivative[radius_index])
            estimate = fit_rate_from_metric_derivative(
                derivative_metric, radius_index)
            rates_at_radii.append(estimate)
            vector_closure_at_radii.append(closure_fraction(
                derivative_metric,
                estimate,
                radius_index,
                VECTOR_CLOSURE_BLOCKS,
            ))
            clock_closure_at_radii.append(closure_fraction(
                derivative_metric,
                estimate,
                radius_index,
                CLOCK_CLOSURE_BLOCKS,
            ))

        inverse_metrics.append(inverse_metric)
        fixed_derivatives.append(fixed_derivative)
        comoving_derivatives.append(comoving_derivative)
        base_rate_estimates.append(rates_at_radii)
        vector_closure.append(vector_closure_at_radii)
        clock_closure.append(clock_closure_at_radii)
        print(
            f"t={time:.1f}: |v_center|="
            f"{np.linalg.norm(curvature_velocity[time_index]):.6f}, "
            f"median closures vector="
            f"{np.median(vector_closure_at_radii):.3e}, "
            f"clock={np.median(clock_closure_at_radii):.3e}"
        )

    inverse_metrics = np.asarray(inverse_metrics)
    fixed_derivatives = np.asarray(fixed_derivatives)
    comoving_derivatives = np.asarray(comoving_derivatives)
    base_rate_estimates = np.asarray(base_rate_estimates)
    vector_closure = np.asarray(vector_closure)
    clock_closure = np.asarray(clock_closure)

    # Direct validation of the GH time derivative at fixed local coordinates.
    time_series_metric_derivative = savgol_filter(
        inverse_metrics,
        window_length=7,
        polyorder=3,
        deriv=1,
        delta=float(np.median(np.diff(times))),
        axis=0,
        mode="interp",
    )
    derivative_validation = np.sqrt(np.sum(
        (comoving_derivatives - time_series_metric_derivative)**2,
        axis=(2, 3, 4),
    )) / np.maximum(
        np.maximum(
            np.sqrt(np.sum(
                comoving_derivatives**2, axis=(2, 3, 4))),
            np.sqrt(np.sum(
                time_series_metric_derivative**2, axis=(2, 3, 4))),
        ),
        1.0e-300,
    )

    # Convert the direct simulation-time rates to centered harmonic time using
    # dt/dT=1/(1+qdot0).  Each alternation has its own clock-rate estimate.
    number_of_stages = second_order.shape[2]
    direct_harmonic_rates = np.empty(
        (len(times), len(ode.RADII), number_of_stages, 13))
    finite_difference_harmonic_rates = np.empty_like(
        direct_harmonic_rates)
    for stage in range(number_of_stages):
        qdot = first_order[:, :, stage, 0]
        direct_harmonic_rates[:, :, stage] = (
            (1.0 + qdot)[..., None] * base_rate_estimates
        )
        coefficient_derivative = savgol_filter(
            first_order[:, :, stage],
            window_length=7,
            polyorder=3,
            deriv=1,
            delta=float(np.median(np.diff(times))),
            axis=0,
            mode="interp",
        )
        finite_difference_harmonic_rates[:, :, stage] = (
            (1.0 + qdot)[..., None] * coefficient_derivative
        )

    # Save the expensive interpolated fields before plotting/post-processing.
    # This checkpoint is intentionally overwritten by the complete archive
    # below once all diagnostics succeed.
    np.savez(
        HERE / "q8_direct_metric_time_derivative_checkpoint.npz",
        times=times,
        radii=ode.RADII,
        inverse_metrics=inverse_metrics,
        inverse_metric_comoving_time_derivatives=comoving_derivatives,
        inverse_metric_time_series_derivatives=time_series_metric_derivative,
        derivative_validation=derivative_validation,
        base_simulation_time_rate_estimates=base_rate_estimates,
        direct_harmonic_time_rate_estimates=direct_harmonic_rates,
        finite_difference_harmonic_time_rate_estimates=(
            finite_difference_harmonic_rates),
        vector_closure=vector_closure,
        clock_strain_closure=clock_closure,
    )

    best_radius = make_plots(
        times,
        first_order,
        second_order,
        direct_harmonic_rates,
        finite_difference_harmonic_rates,
        derivative_validation,
        vector_closure,
        clock_closure,
    )

    np.savez(
        OUTPUT_FILE,
        times=times,
        radii=ode.RADII,
        ah_centers=ah_centers,
        curvature_centers=curvature_centers,
        ah_center_velocities=ah_velocity,
        curvature_center_velocities=curvature_velocity,
        low_cadence_curvature_center_velocities=(
            low_cadence_curvature_velocity),
        inverse_metrics=inverse_metrics,
        inverse_metric_fixed_time_derivatives=fixed_derivatives,
        inverse_metric_comoving_time_derivatives=comoving_derivatives,
        inverse_metric_time_series_derivatives=(
            time_series_metric_derivative),
        derivative_validation=derivative_validation,
        base_simulation_time_rate_estimates=base_rate_estimates,
        direct_harmonic_time_rate_estimates=direct_harmonic_rates,
        finite_difference_harmonic_time_rate_estimates=(
            finite_difference_harmonic_rates),
        vector_closure=vector_closure,
        clock_strain_closure=clock_closure,
        first_order_history=first_order,
        second_order_history=second_order,
    )

    central = slice(2, -2)
    stage = 3
    direct = direct_harmonic_rates[central, best_radius, stage]
    finite = finite_difference_harmonic_rates[
        central, best_radius, stage]
    free_second = second_order[central, best_radius, stage]
    direct_strain_rate = (
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(direct[:, 7:13]))
    finite_strain_rate = (
        ode.SMALL_BLACK_HOLE_MASS
        * symmetric_matrix_series(finite[:, 7:13]))
    free_strain_rate = symmetric_rate_from_second_order(free_second)

    print("\nRepresentative comparison:")
    print(
        f"  radius={ode.RADII[best_radius]:.2f} "
        f"({ode.RADII[best_radius]/ode.SMALL_BLACK_HOLE_MASS:.2f} M_B)"
    )
    print(
        "  metric derivative validation=",
        np.mean(derivative_validation[central, best_radius]),
    )
    print(
        "  direct/finite clock mismatch=",
        group_relative_mismatch(direct[:, 0], finite[:, 0], axis=0),
    )
    print(
        "  direct/free clock mismatch=",
        group_relative_mismatch(
            direct[:, 0], free_second[:, -1], axis=0),
    )
    print(
        "  direct/finite strain mismatch=",
        group_relative_mismatch(
            direct_strain_rate, finite_strain_rate,
            axis=(0, 1, 2)),
    )
    print(
        "  direct/free strain mismatch=",
        group_relative_mismatch(
            direct_strain_rate, free_strain_rate,
            axis=(0, 1, 2)),
    )
    print(
        "  RMS direct clock=",
        np.sqrt(np.mean(direct[:, 0]**2)),
        "finite=",
        np.sqrt(np.mean(finite[:, 0]**2)),
        "free=",
        np.sqrt(np.mean(free_second[:, -1]**2)),
    )
    print(
        "  RMS direct strain=",
        np.sqrt(np.mean(direct_strain_rate**2)),
        "finite=",
        np.sqrt(np.mean(finite_strain_rate**2)),
        "free=",
        np.sqrt(np.mean(free_strain_rate**2)),
    )


if __name__ == "__main__":
    main()
