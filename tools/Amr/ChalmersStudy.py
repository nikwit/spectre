# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Compare Chalmers' prediction with independently resampled h/p candidates.

Run with the build's python-spectre, as described in README.md. No evolution is
performed. Added grid points are a storage/work proxy, NOT measured DG runtime.
"""

import argparse
import csv
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legvander
from scipy.special import roots_legendre

from spectre.tools.Amr.Chalmers import (
    estimate_sequence,
    estimate_tensor,
    refinement_choice,
    sample_modal,
)


@lru_cache(maxsize=None)
def reference_grid(order):
    """Independent Gauss integration nodes, not the LGL interpolation grid."""
    return roots_legendre(order)


def approximation_error(function, extents, bounds, order):
    """Physical L2 error, integrated independently on each candidate element."""
    coefficients = sample_modal(function, extents, bounds)
    nodes, weights = reference_grid(order)
    approximation = coefficients
    for axis, extent in enumerate(extents):
        approximation = np.moveaxis(
            np.tensordot(
                legvander(nodes, extent - 1), approximation, axes=(1, axis)
            ),
            0,
            axis,
        )
    coordinates = [
        lower + (nodes + 1) * (upper - lower) / 2 for lower, upper in bounds
    ]
    exact = function(*np.meshgrid(*coordinates, indexing="ij"))
    error_squared = (approximation - exact) ** 2
    for lower, upper in reversed(bounds):
        error_squared = np.tensordot(
            error_squared, weights * (upper - lower) / 2, axes=(-1, 0)
        )
    return float(np.sqrt(error_squared))


def compare(function, extents, *, order, tolerance, metadata):
    """One p+1 and one binary split candidate for each logical direction."""
    bounds = [(-1, 1)] * len(extents)
    coefficients = sample_modal(function, extents, bounds)
    estimates = estimate_tensor(coefficients)
    baseline = approximation_error(function, extents, bounds, order)
    baseline_check = approximation_error(function, extents, bounds, 2 * order)
    points = int(np.prod(extents))
    rows = []
    for axis, estimate in enumerate(estimates):
        enriched = list(extents)
        enriched[axis] += 1
        p_error = approximation_error(function, enriched, bounds, order)
        p_check = approximation_error(function, enriched, bounds, 2 * order)
        children = []
        for child_interval in [(-1, 0), (0, 1)]:
            child = list(bounds)
            child[axis] = child_interval
            children.append(child)
        h_error = float(
            np.hypot(
                *(
                    approximation_error(function, extents, b, order)
                    for b in children
                )
            )
        )
        h_check = float(
            np.hypot(
                *(
                    approximation_error(function, extents, b, 2 * order)
                    for b in children
                )
            )
        )
        p_added = points // extents[axis]
        h_added = points
        # Clip at a numerical accuracy floor before taking logarithms. These
        # gains are descriptive and do not replace the paper's sigma rule.
        p_gain = np.log(max(baseline, 1e-13) / max(p_error, 1e-13)) / p_added
        h_gain = np.log(max(baseline, 1e-13) / max(h_error, 1e-13)) / h_added
        fitted = [e for e in estimate.sequences if np.isfinite(e.log_fit_rms)]
        rows.append(
            {
                **metadata,
                "dim": len(extents),
                "extents": "x".join(map(str, extents)),
                "axis": axis,
                "sigma": estimate.sigma,
                "estimated_error": estimate.error,
                "fit_rms_max": max(
                    (e.log_fit_rms for e in fitted), default=np.nan
                ),
                "insufficient_sequences": sum(
                    e.status == "insufficient_modes" for e in estimate.sequences
                ),
                "choice": refinement_choice(estimate, tolerance),
                "error_base": baseline,
                "error_p": p_error,
                "error_h": h_error,
                "points_base": points,
                "points_p": points + p_added,
                "points_h": points + h_added,
                "log_gain_per_added_point_p": float(p_gain),
                "log_gain_per_added_point_h": float(h_gain),
                "better_gain_per_point": "p" if p_gain >= h_gain else "h",
                "quadrature_abs_change": max(
                    abs(baseline - baseline_check),
                    abs(p_error - p_check),
                    abs(h_error - h_check),
                ),
                "quadrature_scaled_change": max(
                    abs(a - b) / max(b, tolerance)
                    for a, b in [
                        (baseline, baseline_check),
                        (p_error, p_check),
                        (h_error, h_check),
                    ]
                ),
            }
        )
    return rows


def study_1d(tolerance):
    rows = []
    for extent in [3, 4, 5, 6, 8, 10, 12, 16]:
        for width in [0.03, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]:
            for center in [-0.8, -0.37, 0.0, 0.23, 0.75]:
                rows.extend(
                    compare(
                        lambda x: np.exp(-(((x - center) / width) ** 2)),
                        (extent,),
                        order=256,
                        tolerance=tolerance,
                        metadata={
                            "family": "gaussian",
                            "width": width,
                            "center": center,
                        },
                    )
                )
        for frequency in [0.5, 1, 2, 4, 8]:
            for phase in [0, 0.37]:
                rows.extend(
                    compare(
                        lambda x: np.sin(np.pi * frequency * x + phase),
                        (extent,),
                        order=256,
                        tolerance=tolerance,
                        metadata={
                            "family": "sinusoid",
                            "frequency": frequency,
                            "phase": phase,
                        },
                    )
                )
        for name, function in [
            ("constant", lambda x: np.ones_like(x)),
            ("quadratic", lambda x: 1 + x + x**2),
        ]:
            rows.extend(
                compare(
                    function,
                    (extent,),
                    order=256,
                    tolerance=tolerance,
                    metadata={"family": name},
                )
            )
    return rows


def study_3d(tolerance):
    def anisotropic(x, y, z):
        return np.exp(
            -(((x - 0.13) / 0.12) ** 2)
            - ((y + 0.2) / 0.5) ** 2
            - (z / 1.3) ** 2
        )

    def rotated(x, y, z):
        return anisotropic((x + y) / np.sqrt(2), (y - x) / np.sqrt(2), z)

    functions = [
        ("anisotropic_gaussian", anisotropic),
        ("rotated_gaussian", rotated),
        (
            "oblique_wave",
            lambda x, y, z: np.sin(np.pi * (x + 2 * y + 0.5 * z) + 0.37),
        ),
    ]
    rows = []
    for extents in [(4, 6, 8), (8, 8, 8), (12, 10, 8)]:
        for name, function in functions:
            rows.extend(
                compare(
                    function,
                    extents,
                    order=64,
                    tolerance=tolerance,
                    metadata={"family": name},
                )
            )
    return rows


def summarize(rows, tolerance):
    result = {}
    for dim in [1, 3]:
        subset = [r for r in rows if r["dim"] == dim]
        above_target = [r for r in subset if r["error_base"] > tolerance]
        decided = [r for r in above_target if r["choice"] in ["h", "p"]]
        result[f"{dim}d"] = {
            "directional_comparisons": len(subset),
            "choices": dict(Counter(r["choice"] for r in subset)),
            "above_target": len(above_target),
            # In 3D the global error can be dominated by another direction;
            # a hold in one direction is not necessarily a false negative.
            "holds_with_global_error_above_target": sum(
                r["choice"] == "hold" for r in above_target
            ),
            "decided_above_target": len(decided),
            "selected_reduces_error": sum(
                r[f"error_{r['choice']}"] < r["error_base"] for r in decided
            ),
            "selected_reduces_error_by_two": sum(
                r[f"error_{r['choice']}"] < 0.5 * r["error_base"]
                for r in decided
            ),
            "agrees_with_gain_per_point": sum(
                r["choice"] == r["better_gain_per_point"] for r in decided
            ),
            "max_quadrature_abs_change": max(
                r["quadrature_abs_change"] for r in subset
            ),
            "max_quadrature_scaled_change": max(
                r["quadrature_scaled_change"] for r in subset
            ),
        }
    # Report actual h/p predictions at different thresholds without altering
    # the reference run, restricted to measured above-target cases with fits.
    fitted = [
        r
        for r in rows
        if r["dim"] == 1
        and r["error_base"] > tolerance
        and np.isfinite(r["sigma"])
    ]
    result["threshold_sensitivity_1d"] = {
        str(threshold): {
            "cases": len(fitted),
            "agrees_with_gain_per_point": sum(
                ("p" if r["sigma"] > threshold else "h")
                == r["better_gain_per_point"]
                for r in fitted
            ),
        }
        for threshold in [0.5, 1.0, 1.5]
    }
    return result


def plot_results(rows, output, tolerance):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for choice, color in [("h", "#cc6633"), ("p", "#3366aa")]:
        subset = [
            r
            for r in rows
            if r["dim"] == 1
            and r["choice"] == choice
            and r["error_base"] > tolerance
        ]
        axes[0].scatter(
            [r["sigma"] for r in subset],
            [
                r["log_gain_per_added_point_p"]
                - r["log_gain_per_added_point_h"]
                for r in subset
            ],
            s=15,
            alpha=0.6,
            color=color,
            label=f"Rule chooses {choice}",
        )
    axes[0].axvline(1, color="0.4", linestyle="--")
    axes[0].axhline(0, color="0.4", linestyle=":")
    axes[0].set(
        xlabel=r"Fitted $\sigma$ (natural logarithm)",
        ylabel="p gain per added point − h gain per added point",
        title="1D accuracy gain per added grid point",
    )
    axes[0].legend()
    subset = [
        r for r in rows if r["dim"] == 1 and 0 < r["estimated_error"] < np.inf
    ]
    axes[1].loglog(
        [r["error_base"] for r in subset],
        [r["estimated_error"] for r in subset],
        ".",
        alpha=0.5,
    )
    axes[1].plot([1e-16, 10], [1e-16, 10], "--", color="0.4")
    axes[1].set(
        xlabel="Measured L2 approximation error",
        ylabel="Estimated error",
        title="Estimator calibration; not an error bound",
    )
    fig.savefig(output / "refinement_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), layout="constrained")
    for axis, width in zip(axes, [0.05, 0.2, 0.8]):
        subset = [
            r
            for r in rows
            if r["family"] == "gaussian"
            and r["width"] == width
            and r["center"] == 0.23
        ]
        for key, label, style in [
            ("base", "Original", "o-"),
            ("p", "p+1", "s-"),
            ("h", "Two children", "^-"),
        ]:
            axis.semilogy(
                [int(r["extents"]) for r in subset],
                [r[f"error_{key}"] for r in subset],
                style,
                label=label,
            )
        axis.set(
            xlabel="Original number of points",
            ylabel="L2 approximation error",
            title=f"Gaussian width {width}, center 0.23",
        )
    axes[-1].legend()
    fig.savefig(output / "gaussian_errors.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout="constrained")
    modes = np.arange(16)
    spectra = {
        "Clean exponential": np.exp(-1.5 * modes),
        "Same exponential + plateau": np.maximum(
            np.exp(-1.5 * modes), np.exp(-12)
        ),
    }
    for name, coefficients in spectra.items():
        estimate = estimate_sequence(coefficients)
        axes[0].semilogy(
            modes, coefficients, "o-", label=f"{name}: σ={estimate.sigma:.2f}"
        )
    axes[0].set(
        xlabel="Legendre degree",
        ylabel="Coefficient amplitude",
        title="A plateau can turn a p decision into h",
    )
    axes[0].legend(fontsize=8)
    selected = [
        r
        for r in rows
        if r["family"] == "anisotropic_gaussian" and r["extents"] == "8x8x8"
    ]
    axes[1].bar(
        ["x", "y", "z"],
        [r["sigma"] for r in selected],
        color=["#cc6633", "#cc6633", "#3366aa"],
    )
    axes[1].axhline(1, color="0.4", linestyle="--")
    axes[1].set(
        ylabel=r"Directional $\sigma$",
        title="Anisotropic Gaussian on an 8³ cube",
    )
    fig.savefig(output / "spectra_and_directions.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()
    if not np.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error("--tolerance must be finite and positive")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = study_1d(args.tolerance)
    print(f"Completed {len(rows)} one-dimensional comparisons", flush=True)
    rows.extend(study_3d(args.tolerance))
    with (args.output / "comparisons.csv").open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=sorted(set().union(*(r.keys() for r in rows)))
        )
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "tolerance": args.tolerance,
        "threshold": 1.0,
        "fit": "natural-log least squares, all significant modes",
        "cost_proxy": "log(error_before/error_after) / added grid points",
        "quadrature_verified": all(
            r["quadrature_scaled_change"] <= 1e-7 for r in rows
        ),
        **summarize(rows, args.tolerance),
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    plot_results(rows, args.output, args.tolerance)
    print(json.dumps(summary, indent=2))
    # A quadrature failure invalidates the comparison instead of silently
    # presenting underintegrated reference errors as evidence.
    if not summary["quadrature_verified"]:
        raise RuntimeError(
            "Reference quadrature did not converge; increase its order"
        )


if __name__ == "__main__":
    main()
