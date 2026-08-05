#!/usr/bin/env python3
"""Compare minimal block estimators for the direct metric time derivative."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "q8_temporal_ode", HERE / "q8_temporal_ode_consistency.py")
ode = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(ode)

DATA = np.load(HERE / "q8_direct_metric_time_derivative.npz")
TIMES = DATA["times"]
RADII = DATA["radii"]
RAW_DERIVATIVE = DATA["inverse_metric_comoving_time_derivatives"]
FIRST_ORDER = DATA["first_order_history"]
FINITE_RATES = DATA["finite_difference_harmonic_time_rate_estimates"]
STAGE = 3
CENTRAL = slice(2, -2)

VECTOR_CANDIDATES = (
    ("GTT", 1), ("GTI", 0), ("GTI", 2), ("GIJ", 1))
VECTOR_ESTIMATORS = (
    (("GTT", 1), ("GTI", 0)),
    (("GTT", 1), ("GIJ", 1)),
    (("GTI", 0), ("GTI", 2)),
    (("GTI", 0), ("GIJ", 1)),
    (("GTI", 2), ("GIJ", 1)),
)
CLOCK_CANDIDATES = (
    ("GTT", 0), ("GTT", 2), ("GTI", 1), ("GTI", 3),
    ("GIJ", 0), ("GIJ", 2), ("GIJ", 4))
CLOCK_ESTIMATORS = (
    (("GTT", 0), ("GTI", 1)),
    (("GTT", 0), ("GIJ", 0)),
    (("GTT", 0), ("GIJ", 2)),
    (("GTI", 1), ("GIJ", 0)),
    (("GTI", 1), ("GIJ", 2)),
)


def fit_estimator_family(estimators, candidates, parameter_indices):
    estimates = np.empty((
        len(estimators), len(TIMES), len(RADII),
        len(parameter_indices)))
    closures = np.empty((
        len(estimators), len(TIMES), len(RADII)))
    for estimator_index, blocks in enumerate(estimators):
        omitted = tuple(
            block for block in candidates if block not in blocks)
        for time_index in range(len(TIMES)):
            for radius_index in range(len(RADII)):
                derivative_metric = ode.metric_tuple_from_matrix(
                    RAW_DERIVATIVE[time_index, radius_index])
                matrix = ode.projected_design(
                    ode.FIRST_ORDER_COLUMNS[radius_index],
                    blocks,
                    parameter_indices,
                )
                data = ode.projected_data(derivative_metric, blocks)
                estimate, rank, _ = ode.scaled_least_squares(
                    matrix, data)
                if rank != len(parameter_indices):
                    raise RuntimeError("Estimator lost rank")
                estimates[
                    estimator_index, time_index, radius_index] = estimate

                prediction = ode.metric_linear_combination(
                    tuple(
                        ode.FIRST_ORDER_COLUMNS[radius_index][index]
                        for index in parameter_indices
                    ),
                    estimate,
                )
                held_out_data = ode.projected_data(
                    derivative_metric, omitted)
                held_out_prediction = ode.projected_data(
                    prediction, omitted)
                closures[
                    estimator_index, time_index, radius_index] = (
                    np.linalg.norm(
                        held_out_data - held_out_prediction)
                    / max(np.linalg.norm(held_out_data), 1.0e-300)
                )

    qdot = FIRST_ORDER[:, :, STAGE, 0]
    harmonic_estimates = (
        estimates
        * (1.0 + qdot)[None, ..., None]
    )
    mismatch = np.empty((len(estimators), len(RADII)))
    for estimator_index in range(len(estimators)):
        for radius_index in range(len(RADII)):
            direct = harmonic_estimates[
                estimator_index, CENTRAL, radius_index]
            finite = FINITE_RATES[
                CENTRAL, radius_index, STAGE][:, parameter_indices]
            mismatch[estimator_index, radius_index] = (
                np.linalg.norm(direct - finite)
                / max(
                    np.linalg.norm(direct),
                    np.linalg.norm(finite),
                    1.0e-300,
                )
            )
    return (
        harmonic_estimates,
        mismatch,
        np.mean(closures[:, CENTRAL], axis=1),
    )


(
    vector_estimates,
    vector_mismatch,
    vector_closure,
) = fit_estimator_family(
    VECTOR_ESTIMATORS,
    VECTOR_CANDIDATES,
    ode.VECTOR_PARAMETER_INDICES,
)
(
    clock_estimates,
    clock_mismatch,
    clock_closure,
) = fit_estimator_family(
    CLOCK_ESTIMATORS,
    CLOCK_CANDIDATES,
    ode.CLOCK_STRAIN_PARAMETER_INDICES,
)


def plot_heatmap(axis, values, labels, title):
    image = axis.imshow(
        values,
        origin="lower",
        aspect="auto",
        norm=colors.LogNorm(
            vmin=max(np.nanmin(values), 1.0e-2),
            vmax=max(np.nanmax(values), 1.0),
        ),
        extent=(
            RADII[0] / ode.SMALL_BLACK_HOLE_MASS,
            RADII[-1] / ode.SMALL_BLACK_HOLE_MASS,
            0.5,
            len(labels) + 0.5,
        ),
    )
    axis.set_yticks(np.arange(1, len(labels) + 1), labels)
    axis.set_xlabel(r"$R/M_B$")
    axis.set_title(title)
    plt.colorbar(image, ax=axis)


fig, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.0), constrained_layout=True)
plot_heatmap(
    axes[0, 0],
    vector_mismatch,
    [f"V{index}" for index in range(1, 6)],
    "Vector: direct vs finite-difference rates",
)
plot_heatmap(
    axes[0, 1],
    clock_mismatch,
    [f"C{index}" for index in range(1, 6)],
    "Clock/strain: direct vs finite-difference rates",
)
plot_heatmap(
    axes[1, 0],
    vector_closure,
    [f"V{index}" for index in range(1, 6)],
    "Vector held-out derivative closure",
)
plot_heatmap(
    axes[1, 1],
    clock_closure,
    [f"C{index}" for index in range(1, 6)],
    "Clock/strain held-out derivative closure",
)
fig.suptitle("Minimal estimators for the inverse-metric time derivative")
fig.savefig(HERE / "q8_direct_derivative_estimators.png", dpi=180)
fig.savefig(HERE / "q8_direct_derivative_estimators.pdf")
plt.close(fig)

np.savez(
    HERE / "q8_direct_derivative_estimators.npz",
    times=TIMES,
    radii=RADII,
    vector_estimates=vector_estimates,
    vector_mismatch=vector_mismatch,
    vector_closure=vector_closure,
    clock_strain_estimates=clock_estimates,
    clock_strain_mismatch=clock_mismatch,
    clock_strain_closure=clock_closure,
)

print("Best vector direct/finite agreement:")
index = np.unravel_index(
    np.argmin(vector_mismatch), vector_mismatch.shape)
print(
    f"  V{index[0] + 1}, R={RADII[index[1]]:.2f}: "
    f"{vector_mismatch[index]:.6g}, "
    f"closure={vector_closure[index]:.6g}")
print("Best clock/strain direct/finite agreement:")
index = np.unravel_index(
    np.argmin(clock_mismatch), clock_mismatch.shape)
print(
    f"  C{index[0] + 1}, R={RADII[index[1]]:.2f}: "
    f"{clock_mismatch[index]:.6g}, "
    f"closure={clock_closure[index]:.6g}")
print("Best vector held-out closure:")
index = np.unravel_index(
    np.argmin(vector_closure), vector_closure.shape)
print(
    f"  V{index[0] + 1}, R={RADII[index[1]]:.2f}: "
    f"{vector_closure[index]:.6g}, "
    f"rate mismatch={vector_mismatch[index]:.6g}")
print("Best clock/strain held-out closure:")
index = np.unravel_index(
    np.argmin(clock_closure), clock_closure.shape)
print(
    f"  C{index[0] + 1}, R={RADII[index[1]]:.2f}: "
    f"{clock_closure[index]:.6g}, "
    f"rate mismatch={clock_mismatch[index]:.6g}")
