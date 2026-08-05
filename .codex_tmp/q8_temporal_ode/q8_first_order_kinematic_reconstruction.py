#!/usr/bin/env python3
"""Cross-check qdot^0 and qdot^i using the direct metric time derivative.

The algebraic first-order metric fit supplies qdot^mu on every slice.  The
Pi/Phi construction instead supplies d(qdot^mu)/dT, where T is simulation
time.  Integrating these rates therefore predicts the *change* in qdot^mu,
with one integration constant per component.

For the spatial translation there is a second, absolute measurement.  If

    T = t + q^0(t),

then dT/dt = 1 + qdot^0 and the curvature-center velocity measured per unit
simulation time obeys

    qdot^i = (1 + qdot^0) dz_center^i/dT.

This script compares all three pieces of information:

1. the algebraic metric fit on each slice;
2. the integral of the Pi/Phi-derived rates;
3. the kinematic curvature-center velocity.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/q8-qdot-reconstruction-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/q8-qdot-reconstruction-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid


HERE = Path(__file__).resolve().parent
INPUT_FILE = HERE / "q8_direct_metric_time_derivative.npz"
OUTPUT_FILE = HERE / "q8_first_order_kinematic_reconstruction.npz"
STAGE = 3
REFERENCE_RADIUS_INDEX = 0


def relative_history_mismatch(left, right):
    """Relative Euclidean mismatch over all times and components."""
    numerator = np.linalg.norm(left - right)
    denominator = max(
        np.linalg.norm(left), np.linalg.norm(right), 1.0e-300)
    return numerator / denominator


def relative_change_mismatch(left, right):
    """Mismatch of histories after removing their initial values."""
    left_change = left - left[0]
    right_change = right - right[0]
    return relative_history_mismatch(left_change, right_change)


def main():
    data = np.load(INPUT_FILE)
    times = data["times"]
    radii = data["radii"]
    first_order = data["first_order_history"][:, :, STAGE]
    simulation_time_rates = data[
        "base_simulation_time_rate_estimates"]
    center_velocity = data["curvature_center_velocities"]

    integrated_from_metric_anchor = np.empty_like(first_order)
    integrated_from_kinematic_anchor = np.empty(
        first_order[:, :, 4:7].shape)
    kinematic_spatial_velocity = np.empty_like(
        integrated_from_kinematic_anchor)

    clock_change_mismatch = np.empty(len(radii))
    spatial_change_mismatch = np.empty(len(radii))
    raw_center_mismatch = np.empty(len(radii))
    kinematic_center_mismatch = np.empty(len(radii))
    kinematically_anchored_integral_mismatch = np.empty(len(radii))

    for radius_index in range(len(radii)):
        integrated = (
            first_order[0, radius_index]
            + cumulative_trapezoid(
                simulation_time_rates[:, radius_index],
                times,
                axis=0,
                initial=0.0,
            )
        )
        integrated_from_metric_anchor[:, radius_index] = integrated

        # Use the directly reconstructed clock history in the T-to-t
        # conversion.  Only its value on the first slice comes from the
        # algebraic metric fit.
        reconstructed_clock = integrated[:, 0]
        kinematic_velocity = (
            (1.0 + reconstructed_clock[:, None]) * center_velocity)
        kinematic_spatial_velocity[:, radius_index] = kinematic_velocity

        # Anchor the integrated spatial acceleration to the independently
        # measured center velocity instead of to the fitted qdot^i.
        integrated_from_kinematic_anchor[:, radius_index] = (
            kinematic_velocity[0]
            + cumulative_trapezoid(
                simulation_time_rates[:, radius_index, 4:7],
                times,
                axis=0,
                initial=0.0,
            )
        )

        clock_change_mismatch[radius_index] = relative_change_mismatch(
            first_order[:, radius_index, 0], integrated[:, 0])
        spatial_change_mismatch[radius_index] = relative_change_mismatch(
            first_order[:, radius_index, 4:7], integrated[:, 4:7])
        raw_center_mismatch[radius_index] = relative_history_mismatch(
            first_order[:, radius_index, 4:7], center_velocity)
        kinematic_center_mismatch[radius_index] = relative_history_mismatch(
            first_order[:, radius_index, 4:7], kinematic_velocity)
        kinematically_anchored_integral_mismatch[
            radius_index
        ] = relative_history_mismatch(
            first_order[:, radius_index, 4:7],
            integrated_from_kinematic_anchor[:, radius_index],
        )

    radius_index = REFERENCE_RADIUS_INDEX
    radius = radii[radius_index]
    fitted = first_order[:, radius_index]
    integrated = integrated_from_metric_anchor[:, radius_index]
    kinematic = kinematic_spatial_velocity[:, radius_index]
    kinematic_integral = integrated_from_kinematic_anchor[:, radius_index]

    fig, axes = plt.subplots(
        2, 2, figsize=(14.0, 9.5), constrained_layout=True)

    axes[0, 0].plot(
        times, fitted[:, 0], "o", label="algebraic metric fit")
    axes[0, 0].plot(
        times,
        integrated[:, 0],
        "-",
        label=r"integrated direct $d\dot q^0/dT$",
    )
    axes[0, 0].set_title(
        rf"Clock rate at $R={radius:.2f}M$: one initial calibration")
    axes[0, 0].set_xlabel(r"$T/M$")
    axes[0, 0].set_ylabel(r"$\dot q^0$")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    colors = ("tab:blue", "tab:orange")
    component_labels = ("x", "y")
    for component, (color, label) in enumerate(
            zip(colors, component_labels)):
        parameter_index = 4 + component
        axes[0, 1].plot(
            times,
            fitted[:, parameter_index],
            "o",
            color=color,
            label=rf"metric fit $\dot q^{label}$",
        )
        axes[0, 1].plot(
            times,
            kinematic_integral[:, component],
            "-",
            color=color,
            label=rf"integrated $\ddot q^{label}$",
        )
        axes[0, 1].plot(
            times,
            kinematic[:, component],
            "--",
            color=color,
            alpha=0.75,
            label=rf"$(1+\dot q^0)v_{{\rm center}}^{label}$",
        )
    axes[0, 1].set_title(
        rf"Spatial velocity at $R={radius:.2f}M$")
    axes[0, 1].set_xlabel(r"$T/M$")
    axes[0, 1].set_ylabel(r"$\dot q^i$")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend(ncol=2, fontsize=9)

    axes[1, 0].loglog(
        radii,
        clock_change_mismatch,
        "o-",
        label=r"$\dot q^0$: integrated-rate change",
    )
    axes[1, 0].loglog(
        radii,
        spatial_change_mismatch,
        "s-",
        label=r"$\dot q^i$: integrated-rate change",
    )
    axes[1, 0].set_title("Does the direct derivative reproduce the evolution?")
    axes[1, 0].set_xlabel(r"worldtube radius $R/M$")
    axes[1, 0].set_ylabel("relative mismatch of changes")
    axes[1, 0].grid(which="both", alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].loglog(
        radii,
        raw_center_mismatch,
        "o-",
        label=r"raw $v_{\rm center}^i$",
    )
    axes[1, 1].loglog(
        radii,
        kinematic_center_mismatch,
        "s-",
        label=r"$(1+\dot q^0)v_{\rm center}^i$",
    )
    axes[1, 1].loglog(
        radii,
        kinematically_anchored_integral_mismatch,
        "^-",
        label=r"integrated $\ddot q^i$, center-velocity anchor",
    )
    axes[1, 1].set_title("Absolute spatial-velocity cross-check")
    axes[1, 1].set_xlabel(r"worldtube radius $R/M$")
    axes[1, 1].set_ylabel("relative history mismatch")
    axes[1, 1].grid(which="both", alpha=0.25)
    axes[1, 1].legend()

    fig.suptitle(
        r"Independent first-order velocity checks from $\Pi,\Phi$ "
        "and center motion",
        fontsize=15,
    )
    figure_png = HERE / "q8_first_order_kinematic_reconstruction.png"
    figure_pdf = HERE / "q8_first_order_kinematic_reconstruction.pdf"
    fig.savefig(figure_png, dpi=180)
    fig.savefig(figure_pdf)
    plt.close(fig)

    np.savez(
        OUTPUT_FILE,
        times=times,
        radii=radii,
        stage=STAGE,
        fitted_first_order=first_order,
        simulation_time_rate_estimates=simulation_time_rates,
        integrated_from_metric_anchor=integrated_from_metric_anchor,
        curvature_center_velocity=center_velocity,
        kinematic_spatial_velocity=kinematic_spatial_velocity,
        integrated_from_kinematic_anchor=(
            integrated_from_kinematic_anchor),
        clock_change_mismatch=clock_change_mismatch,
        spatial_change_mismatch=spatial_change_mismatch,
        raw_center_mismatch=raw_center_mismatch,
        kinematic_center_mismatch=kinematic_center_mismatch,
        kinematically_anchored_integral_mismatch=(
            kinematically_anchored_integral_mismatch),
    )

    best_spatial = int(np.argmin(spatial_change_mismatch))
    best_kinematic = int(np.argmin(kinematic_center_mismatch))
    print("First-order velocity cross-check")
    print(
        f"  R={radius:.2f}: clock evolution mismatch = "
        f"{clock_change_mismatch[radius_index]:.3%}")
    print(
        f"  R={radius:.2f}: spatial evolution mismatch = "
        f"{spatial_change_mismatch[radius_index]:.3%}")
    print(
        f"  R={radius:.2f}: harmonic-time center mismatch = "
        f"{kinematic_center_mismatch[radius_index]:.3%}")
    print(
        "  best spatial evolution mismatch: "
        f"{spatial_change_mismatch[best_spatial]:.3%} at "
        f"R={radii[best_spatial]:.2f}")
    print(
        "  best harmonic-time center mismatch: "
        f"{kinematic_center_mismatch[best_kinematic]:.3%} at "
        f"R={radii[best_kinematic]:.2f}")
    print(f"  wrote {OUTPUT_FILE}")
    print(f"  wrote {figure_png}")


if __name__ == "__main__":
    main()
