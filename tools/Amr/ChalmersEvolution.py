# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Reproducible native scalar-wave experiment for the Chalmers hp criterion.

Run with python-spectre. Each run has its own input, log and HDF5 files. Outputs
are never overwritten. Subcommands separate execution from postprocessing.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import time
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np
import yaml
from numpy.polynomial.legendre import legder, leggauss, legval

from spectre.DataStructures import DataVector
from spectre.NumericalAlgorithms.LinearOperators import chalmers_estimate
from spectre.Spectral import Basis, Mesh, Quadrature
from spectre.tools.Amr.Chalmers import to_modal


def input_file(
    strategy,
    tolerance,
    dt,
    end_time,
    volume,
    underresolved=False,
    coarsening_factor=0.0,
    h_coarsening=False,
):
    def interval(period):
        count = round(period / dt)
        if not np.isclose(count * dt, period, rtol=0, atol=1e-12):
            raise ValueError(
                "Time step must divide all observation/AMR periods"
            )
        return {"Slabs": {"EvenlySpaced": {"Interval": count, "Offset": 0}}}

    final_slab = round(end_time / dt)
    if not np.isclose(final_slab * dt, end_time, rtol=0, atol=1e-12):
        raise ValueError("Time step must divide end time")
    data = {
        "PlaneWave": {
            "WaveVector": [1.0],
            "Center": [-2.0],
            "Profile": {
                "Gaussian": {"Amplitude": 1.0, "Width": 0.3, "Center": 0.0}
            },
        }
    }
    criterion = {
        "Chalmers": {
            "VariablesToMonitor": ["Psi", "Pi", "Phi"],
            "AbsoluteTarget": tolerance,
            "RelativeTarget": 0.0,
            "SmoothnessThreshold": 0.5 if strategy == "hp05" else 1.0,
            "PRefineOnly": strategy == "p",
            "CoarseningFactor": coarsening_factor,
            "HCoarsening": h_coarsening,
        }
    }
    events = [
        {
            "ObserveAmrStats": {
                "PrintToTerminal": False,
                "ObservePerCore": False,
            }
        },
        {
            "ObserveTimeStep": {
                "SubfileName": "TimeSteps",
                "PrintTimeToTerminal": False,
                "ObservePerCore": False,
            }
        },
        {
            "ObserveNorms": {
                "SubfileName": "Errors",
                "TensorsToObserve": [
                    {
                        "Name": f"Error({field})",
                        "NormType": "L2Norm",
                        "Components": "Sum",
                    }
                    for field in ["Psi", "Pi", "Phi"]
                ],
            }
        },
    ]
    if volume:
        events.append(
            {
                "ObserveFields": {
                    "SubfileName": "Fields",
                    "VariablesToObserve": ["Psi", "Pi", "Phi"],
                    "InterpolateToMesh": "None",
                    "ProjectToMesh": "None",
                    "CoordinatesFloatingPointType": "Double",
                    "FloatingPointTypes": ["Double"],
                    "BlocksToObserve": "All",
                }
            }
        )
    return {
        "Parallelization": {"ElementDistribution": "NumGridPoints"},
        "ResourceInfo": {"AvoidGlobalProc0": False, "Singletons": "Auto"},
        "InitialData": data,
        "Amr": {
            "Criteria": [] if strategy == "fixed" else [criterion],
            "EnableInBlocks": "All",
            "Policies": {
                "EnforceTwoToOneBalanceInNormalDirection": True,
                "Isotropy": "Anisotropic",
                "Limits": {
                    "RefinementLevel": [0, 7],
                    "NumPolynomialModes": [5, 20],
                    "FourierM": "Auto",
                    "SphericalHarmonicL": "Auto",
                    "ErrorBeyondLimits": True,
                },
                "AllowCoarsening": coarsening_factor > 0,
            },
            "Verbosity": "Quiet",
        },
        "PhaseChangeAndTriggers": (
            []
            if strategy == "fixed"
            else [
                {
                    "Trigger": interval(0.02),
                    "PhaseChanges": [
                        "VisitAndReturn(EvaluateAmrCriteria)",
                        "VisitAndReturn(AdjustDomain)",
                        "VisitAndReturn(CheckDomain)",
                    ],
                }
            ]
        ),
        "Evolution": {
            "InitialTime": 0.0,
            "InitialTimeStep": dt,
            "InitialSlabSize": dt,
            "MinimumTimeStep": 1e-8,
            "LocalTimeStepping": "Off",
            "TimeStepper": {"DormandPrince5": None},
            "LtsStepChoosers": None,
            "VariableOrderAlgorithm": {"GoalOrder": 5},
        },
        "DomainCreator": {
            "RotatedIntervals": {
                "LowerBound": [-4.0],
                "Midpoint": [0.0],
                "UpperBound": [8.0],
                "InitialRefinement": [4 if strategy == "fixed" else 3],
                "InitialGridPoints": [
                    (
                        [20, 20]
                        if strategy == "fixed"
                        else [6 if underresolved else 16, 5]
                    )
                ],
                "TimeDependence": "None",
                "BoundaryConditions": {
                    side: {"DirichletAnalytic": {"AnalyticPrescription": data}}
                    for side in ["LowerBoundary", "UpperBoundary"]
                },
            }
        },
        "SpatialDiscretization": {
            "BoundaryCorrection": {"UpwindPenalty": None},
            "DiscontinuousGalerkin": {
                "Formulation": "StrongInertial",
                "Quadrature": "GaussLobatto",
                "Filtering": [{"None": {"BlocksToFilter": "All"}}],
            },
        },
        "EventsAndTriggersAtSlabs": [
            {"Trigger": interval(0.2), "Events": events},
            {
                "Trigger": {"Slabs": {"Specified": {"Values": [final_slab]}}},
                "Events": ["Completion"],
            },
        ],
        "EventsAndTriggersAtSteps": None,
        "EventsAndTriggersAtCheckpoints": None,
        "EventsAndDenseTriggers": None,
        "Observers": {
            "VolumeFileName": "Volume",
            "ReductionFileName": "Reductions",
        },
    }


def run(args):
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    configuration = input_file(
        args.strategy,
        args.tolerance,
        args.dt,
        args.end_time,
        not args.no_volume,
        args.underresolved,
        coarsening_factor=getattr(args, "coarsening_factor", 0.0),
        h_coarsening=getattr(args, "h_coarsening", False),
    )
    if getattr(args, "initial_points", None) is not None:
        configuration["DomainCreator"]["RotatedIntervals"][
            "InitialGridPoints"
        ] = [[args.initial_points, args.initial_points]]
    (directory / "Input.yaml").write_text(
        "Executable: EvolveScalarWave1D\n---\n"
        + yaml.safe_dump(configuration, sort_keys=False)
    )
    command = [
        str(args.executable.resolve()),
        f"+p{getattr(args, 'processes', 1)}",
        "--input-file",
        "Input.yaml",
    ]
    started = time.perf_counter()
    with (directory / "Run.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=directory,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "OPENBLAS_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
            },
        )
    metadata = {
        "strategy": args.strategy,
        "tolerance": args.tolerance,
        "dt": args.dt,
        "end_time": args.end_time,
        "volume": not args.no_volume,
        "underresolved": args.underresolved,
        "coarsening_factor": getattr(args, "coarsening_factor", 0.0),
        "h_coarsening": getattr(args, "h_coarsening", False),
        "initial_points": getattr(args, "initial_points", None),
        "processes": getattr(args, "processes", 1),
        "executable_sha256": (
            hashlib.sha256(args.executable.read_bytes()).hexdigest()
        ),
        "returncode": result.returncode,
        "process_seconds": time.perf_counter() - started,
        "command": command,
    }
    (directory / "Run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"run": str(directory), **metadata}), flush=True)
    if result.returncode:
        print((directory / "Run.log").read_text()[-7000:])
        raise SystemExit(result.returncode)


@lru_cache(None)
def gauss(order):
    return leggauss(order)


def exact(x, time):
    u = x + 2.0 - time
    psi = np.exp(-((u / 0.3) ** 2))
    derivative = -2.0 * u / 0.3**2 * psi
    return {"Psi": psi, "Pi": derivative, "Phi_x": derivative}


def elements(group):
    names = bytes(group["grid_names"][:]).decode().split(":")
    offset = 0
    for name, extent in zip(names, group["total_extents"][:]):
        extent = int(extent)
        selection = slice(offset, offset + extent)
        x = group["InertialCoordinates_x"][selection]
        values = {
            field: group[field][selection] for field in ["Psi", "Pi", "Phi_x"]
        }
        yield name, x, values
        offset += extent


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze_run(directory):
    metadata = json.loads((directory / "Run.json").read_text())
    result = {"run": directory.name, **metadata}
    if metadata["returncode"] != 0:
        log = (directory / "Run.log").read_text()
        reason = re.search(
            r"Tried refining beyond[^\n]+|Chalmers fit is undetermined[^\n]+",
            log,
        )
        result["failure_reason"] = reason[0] if reason else "See Run.log"
        reductions_path = directory / "Reductions.h5"
        if reductions_path.exists():
            with h5py.File(reductions_path) as reductions:
                if "Amr.dat" in reductions and len(reductions["Amr.dat"]):
                    row = reductions["Amr.dat"][-1]
                    result.update(
                        last_observed_time=float(row[0]),
                        last_elements=int(row[1]),
                        last_points=int(row[2]),
                        last_max_degree=int(row[-1] - 1),
                    )
        (directory / "Analysis.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result
    with h5py.File(directory / "Reductions.h5") as reductions:
        stats = reductions["Amr.dat"][:]
        timing = reductions["TimeSteps.dat"][:]
        result.update(
            evolution_seconds=float(timing[-1, -1] - timing[0, -1]),
            final_elements=int(stats[-1, 1]),
            final_points=int(stats[-1, 2]),
            min_degree=int(stats[-1, -2] - 1),
            max_degree=int(stats[-1, -1] - 1),
            net_element_change=int(stats[-1, 1] - stats[0, 1]),
        )
    volume = directory / "Volume0.h5"
    if not volume.exists():
        (directory / "Analysis.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result
    history, diagnostics, spectra, meshes = [], [], [], []
    worst_quadrature_change = 0.0
    with h5py.File(volume) as f:
        observations = sorted(
            (g for g in f["Fields.vol"].values() if isinstance(g, h5py.Group)),
            key=lambda g: g.attrs["observation_value"],
        )
        for group in observations:
            time = float(group.attrs["observation_value"])
            errors = {order: np.zeros(3) for order in [64, 128]}
            constraint_squared = 0.0
            num_elements, num_points = 0, 0
            for name, x, values in elements(group):
                n = len(x)
                num_elements += 1
                num_points += n
                center, half_width = (x[0] + x[-1]) / 2, (x[-1] - x[0]) / 2
                jacobian = abs(half_width)
                coefficients = {
                    field: to_modal(v) for field, v in values.items()
                }
                meshes.append(
                    {
                        "time": time,
                        "element": name,
                        "left": min(x),
                        "right": max(x),
                        "degree": n - 1,
                    }
                )
                for order in errors:
                    xi, weights = gauss(order)
                    analytic = exact(center + half_width * xi, time)
                    for i, (field, modes) in enumerate(coefficients.items()):
                        difference = legval(xi, modes) - analytic[field]
                        errors[order][i] += jacobian * np.dot(
                            weights, difference**2
                        )
                    if order == 128:
                        constraint = legval(
                            xi, legder(coefficients["Psi"])
                        ) / half_width - legval(xi, coefficients["Phi_x"])
                        constraint_squared += jacobian * np.dot(
                            weights, constraint**2
                        )
                mesh = Mesh[1]([n], Basis.Legendre, Quadrature.GaussLobatto)
                for field, v in values.items():
                    fit = chalmers_estimate(DataVector(v), mesh)[0]
                    diagnostics.append(
                        {
                            "time": time,
                            "element": name,
                            "field": field,
                            "degree": n - 1,
                            "amplitude": float(np.max(np.abs(v))),
                            **fit,
                        }
                    )
                # Keep the full Pi spectrum, including the analytic nodal
                # interpolant and numerical contamination on the same mesh.
                analytic_modes = to_modal(exact(x, time)["Pi"])
                for k, (mode, analytic_mode) in enumerate(
                    zip(coefficients["Pi"], analytic_modes)
                ):
                    spectra.append(
                        {
                            "time": time,
                            "element": name,
                            "center": center,
                            "mode": k,
                            "amplitude": abs(mode),
                            "analytic_amplitude": abs(analytic_mode),
                            "error_amplitude": abs(mode - analytic_mode),
                        }
                    )
            error = np.sqrt(errors[128])
            change = float(np.max(np.abs(error - np.sqrt(errors[64]))))
            worst_quadrature_change = max(worst_quadrature_change, change)
            if change > max(2e-13, 1e-7 * max(error)):
                raise RuntimeError(
                    f"Quadrature did not converge in {directory} at {time}:"
                    f" {change}"
                )
            history.append(
                {
                    "time": time,
                    "elements": num_elements,
                    "points": num_points,
                    "error_psi": error[0],
                    "error_pi": error[1],
                    "error_phi": error[2],
                    "error_state": np.linalg.norm(error),
                    "constraint": np.sqrt(constraint_squared),
                }
            )
    if not np.isclose(history[-1]["time"], metadata["end_time"], atol=1e-10):
        raise RuntimeError(f"Missing final snapshot in {directory}")
    result.update(
        {
            f"final_{key}": float(value)
            for key, value in history[-1].items()
            if key.startswith("error") or key == "constraint"
        }
    )
    result["initial_error_state"] = float(history[0]["error_state"])
    result["max_error_state"] = float(
        max(row["error_state"] for row in history)
    )
    result["quadrature_change"] = worst_quadrature_change
    write_csv(directory / "History.csv", history)
    write_csv(directory / "Fits.csv", diagnostics)
    write_csv(directory / "Spectra.csv", spectra)
    write_csv(directory / "Mesh.csv", meshes)
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def analyze(args):
    results = []
    for directory in args.runs:
        result = analyze_run(directory)
        results.append(result)
        print(
            json.dumps({k: v for k, v in result.items() if k != "command"}),
            flush=True,
        )
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "Summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    plot_results(args.runs, results, args.output)


def final_polynomials(directory):
    with h5py.File(directory / "Volume0.h5") as f:
        group = max(
            (g for g in f["Fields.vol"].values() if isinstance(g, h5py.Group)),
            key=lambda g: g.attrs["observation_value"],
        )
        return float(group.attrs["observation_value"]), [
            (x[0], x[-1], np.array([to_modal(v) for v in values.values()]))
            for _, x, values in elements(group)
        ]


def difference(args):
    """Integrate the difference on a common mesh, without matching nodal grids."""
    time_a, elements_a = final_polynomials(args.first)
    time_b, elements_b = final_polynomials(args.second)
    if not np.isclose(time_a, time_b, atol=1e-10, rtol=0):
        raise ValueError("Snapshots are at different times")
    boundaries = sorted(set(x for e in elements_a + elements_b for x in e[:2]))
    xi, weights = gauss(64)
    squared = np.zeros(3)
    error_difference_squared = np.zeros(3)

    def evaluate(elements, center, points):
        for left, right, coefficients in elements:
            if min(left, right) <= center <= max(left, right):
                logical = (2 * points - left - right) / (right - left)
                return np.array([legval(logical, c) for c in coefficients])
        raise ValueError("No containing element")

    for left, right in zip(boundaries[:-1], boundaries[1:]):
        center = (left + right) / 2
        half_width = (right - left) / 2
        x = center + half_width * xi
        delta = evaluate(elements_a, center, x) - evaluate(
            elements_b, center, x
        )
        squared += half_width * np.sum(weights * delta**2, axis=1)
        # Slab accumulation gives slightly different floating-point final
        # times. Remove only the known analytic phase shift, not any numerical
        # error, and report both comparisons.
        analytic_shift = np.array(list(exact(x, time_a).values())) - np.array(
            list(exact(x, time_b).values())
        )
        error_difference_squared += half_width * np.sum(
            weights * (delta - analytic_shift) ** 2, axis=1
        )
    result = {
        "first": str(args.first),
        "second": str(args.second),
        "time_difference": time_a - time_b,
        "l2_difference": np.sqrt(squared).tolist(),
        "state_difference": float(np.sqrt(sum(squared))),
        "phase_aligned_state_difference": float(
            np.sqrt(sum(error_difference_squared))
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


def plot_results(directories, results, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"p": "#2271b2", "hp1": "#d55e00", "hp05": "#009e73"}
    labels = {"p": "p only", "hp1": "hp, sigma > 1", "hp05": "hp, sigma > 0.5"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for strategy in colors:
        rows = sorted(
            (
                r
                for r in results
                if r["strategy"] == strategy and "final_error_state" in r
            ),
            key=lambda r: r["tolerance"],
        )
        if not rows:
            continue
        axes[0, 0].loglog(
            [r["evolution_seconds"] for r in rows],
            [r["final_error_state"] for r in rows],
            "o-",
            color=colors[strategy],
            label=labels[strategy],
        )
        axes[0, 1].loglog(
            [r["tolerance"] for r in rows],
            [r["final_error_state"] for r in rows],
            "o-",
            color=colors[strategy],
            label=labels[strategy],
        )
    for directory, result in zip(directories, results):
        strategy = result["strategy"]
        if (
            strategy not in colors
            or result["tolerance"] != 1e-8
            or "final_error_state" not in result
        ):
            continue
        history = np.genfromtxt(
            directory / "History.csv", delimiter=",", names=True
        )
        axes[1, 0].semilogy(
            history["time"],
            history["error_state"],
            color=colors[strategy],
            label=labels[strategy],
        )
        axes[1, 1].plot(
            history["time"],
            history["points"],
            color=colors[strategy],
            label=labels[strategy],
        )
    for ax in axes.flat:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
    axes[0, 0].set(
        xlabel="Evolution wall time with diagnostics [s]",
        ylabel="Final physical L2 error, combined fields",
        title="Local build: accuracy and measured cost",
    )
    axes[0, 1].set(
        xlabel="Absolute indicator target",
        ylabel="Final physical L2 error, combined fields",
        title="Tolerance sweep (failed runs omitted)",
    )
    axes[1, 0].set(
        xlabel="Time",
        ylabel="Physical L2 error, combined fields",
        title="Resolved initial pulse, target 1e-8",
    )
    axes[1, 1].set(
        xlabel="Time",
        ylabel="Total grid points",
        title="Refinement without coarsening, target 1e-8",
    )
    fig.savefig(output / "EvolutionComparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    for directory, result in zip(directories, results):
        strategy = result["strategy"]
        if (
            strategy not in colors
            or result["tolerance"] != 1e-8
            or "final_error_state" not in result
        ):
            continue
        _, polynomials = final_polynomials(directory)
        for left, right, modes in polynomials:
            axes[0].plot(
                [left, right],
                [modes.shape[1] - 1] * 2,
                color=colors[strategy],
                alpha=0.7,
            )
            axes[0].plot(
                [left, right],
                [modes.shape[1] - 1] * 2,
                "|",
                color=colors[strategy],
                ms=8,
            )
        with (directory / "Spectra.csv").open() as stream:
            spectra = list(csv.DictReader(stream))
        final_time = max(float(row["time"]) for row in spectra)
        final = [row for row in spectra if float(row["time"]) == final_time]
        selected = min(final, key=lambda row: abs(float(row["center"]) - 3.8))[
            "element"
        ]
        subset = [row for row in final if row["element"] == selected]
        k = [int(row["mode"]) for row in subset]
        axes[1].semilogy(
            k,
            [max(float(row["amplitude"]), 1e-17) for row in subset],
            "o-",
            color=colors[strategy],
            label=labels[strategy] + ": Pi",
        )
        axes[1].semilogy(
            k,
            [max(float(row["error_amplitude"]), 1e-17) for row in subset],
            "--",
            color=colors[strategy],
            label=labels[strategy] + ": modal error",
        )
    axes[0].set(
        xlabel="Physical coordinate",
        ylabel="Polynomial degree",
        title="Final mesh, target 1e-8; ticks mark element boundaries",
    )
    axes[1].set(
        xlabel="Legendre mode",
        ylabel="Coefficient amplitude",
        title="Element nearest x=3.8 at final time; target 1e-8",
    )
    for ax in axes:
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8, ncol=2)
    fig.savefig(output / "MeshAndSpectra.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    runner = commands.add_parser("run")
    runner.add_argument("--executable", type=Path, required=True)
    runner.add_argument("--output", type=Path, required=True)
    runner.add_argument(
        "--strategy", choices=["fixed", "p", "hp1", "hp05"], required=True
    )
    runner.add_argument("--tolerance", type=float, default=1e-6)
    runner.add_argument("--dt", type=float, default=0.0005)
    runner.add_argument("--end-time", type=float, default=6.0)
    runner.add_argument("--no-volume", action="store_true")
    runner.add_argument("--underresolved", action="store_true")
    runner.add_argument("--coarsening-factor", type=float, default=0.0)
    runner.add_argument("--h-coarsening", action="store_true")
    runner.add_argument("--initial-points", type=int)
    runner.add_argument("--processes", type=int, default=1)
    runner.set_defaults(func=run)
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("runs", nargs="+", type=Path)
    analyzer.add_argument("--output", type=Path, required=True)
    analyzer.set_defaults(func=analyze)
    differ = commands.add_parser("difference")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.set_defaults(func=difference)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
