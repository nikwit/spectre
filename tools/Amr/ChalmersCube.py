# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Native 3D Gaussian plane-wave AMR experiments on a Cartesian cube.

The exact solution is exp(-((n.x-t+1.8)/0.35)^2). Tensor-product quadrature
checks the evolved polynomial against the analytic solution independently of
the evolution nodes. All paths are local and runs never overwrite data.
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
from itertools import product
from pathlib import Path

import h5py
import numpy as np
import yaml
from numpy.polynomial.legendre import legder, leggauss, legvander

from spectre.DataStructures import DataVector
from spectre.NumericalAlgorithms.LinearOperators import (
    chalmers_energy_estimate,
    chalmers_estimate,
)
from spectre.Spectral import Basis, Mesh, Quadrature
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import input_file as input_1d
from spectre.tools.Amr.ChalmersEvolution import write_csv

FIELDS = ["Psi", "Pi", "Phi_x", "Phi_y", "Phi_z"]
WIDTH = 0.35
CENTER = -1.8
AMR_INTERVAL = 0.04
OUTPUT_INTERVAL = 0.4


def normal(orientation):
    return (
        np.array([1.0, 0.0, 0.0])
        if orientation == "axis"
        else np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
    )


def exact(xyz, t, orientation):
    n = normal(orientation)
    u = sum(k * x for k, x in zip(n, xyz)) - t - CENTER
    psi = np.exp(-((u / WIDTH) ** 2))
    derivative = -2 * u / WIDTH**2 * psi
    return np.array([psi, derivative, *(k * derivative for k in n)])


def configuration(args):
    amr_interval = getattr(args, "amr_interval", AMR_INTERVAL)
    config = input_1d(
        args.strategy,
        args.tolerance,
        args.dt,
        args.end_time,
        True,
        coarsening_factor=getattr(args, "coarsening_factor", 0.0),
        h_coarsening=getattr(args, "h_coarsening", False),
    )
    for criterion in config["Amr"]["Criteria"]:
        if "Chalmers" in criterion:
            criterion["Chalmers"]["UseDirectionalEnergy"] = (
                getattr(args, "aggregation", "energy") == "energy"
            )
    n = normal(args.orientation)
    data = {
        "PlaneWave": {
            "WaveVector": n.tolist(),
            "Center": (CENTER * n).tolist(),
            "Profile": {
                "Gaussian": {"Amplitude": 1.0, "Width": WIDTH, "Center": 0.0}
            },
        }
    }
    config["InitialData"] = data
    points = [args.fixed_p + 1] * 3 if args.strategy == "fixed" else [14] * 3
    if args.orientation == "axis":
        points[1:] = [3, 3]
    regions = (
        []
        if args.strategy == "fixed"
        else [
            {
                "LowerCornerIndex": (
                    [1, 0, 0] if args.orientation == "axis" else [1, 1, 1]
                ),
                "UpperCornerIndex": [2, 2, 2],
                "Refinement": (
                    [5, 3, 3] if args.orientation == "axis" else [5, 5, 5]
                ),
            }
        ]
    )
    config["DomainCreator"] = {
        "AlignedLattice": {
            "BlockBounds": [[-1.0, 0.0, 1.0]] * 3,
            "Distributions": [[], [], []],
            "SingularityPositions": [[], [], []],
            "InitialLevels": [0, 0, 0],
            "InitialGridPoints": points,
            "RefinedLevels": [],
            "RefinedGridPoints": regions,
            "BlocksToExclude": [],
            "BoundaryConditions": (
                [{"DirichletAnalytic": {"AnalyticPrescription": data}}] * 3
            ),
        }
    }
    if args.orientation == "axis" and getattr(args, "axis_periodic", False):
        config["DomainCreator"]["AlignedLattice"]["BoundaryConditions"][1:] = [
            "Periodic",
            "Periodic",
        ]
    config["Amr"]["Policies"]["Limits"]["RefinementLevel"] = [0, args.max_level]
    config["Amr"]["Policies"]["Limits"]["NumPolynomialModes"] = [2, 20]
    # Land exactly on the physical schedule in both fixed-step and CFL runs.
    times = [
        round(i * amr_interval, 12)
        for i in range(round(args.end_time / amr_interval) + 1)
    ]
    observe = [
        round(i * args.observation_interval, 12)
        for i in range(round(args.end_time / args.observation_interval) + 1)
    ]
    if observe[-1] != args.end_time:
        raise ValueError("End time must be divisible by output interval")
    config["PhaseChangeAndTriggers"] = (
        []
        if args.strategy == "fixed"
        else [
            {
                "Trigger": {"Times": {"Specified": {"Values": times[1:-1]}}},
                "PhaseChanges": [
                    "VisitAndReturn(EvaluateAmrCriteria)",
                    "VisitAndReturn(AdjustDomain)",
                    "VisitAndReturn(CheckDomain)",
                ],
            }
        ]
    )
    events = config["EventsAndTriggersAtSlabs"][0]["Events"]
    events = [e for e in events if "ObserveTimeStep" not in e]
    # Keep fields only at endpoints in timing runs, preserving independent
    # final error measurements even when intermediate I/O is disabled.
    fields_event = events.pop()
    config["EventsAndTriggersAtSlabs"] = [
        {
            "Trigger": {"Times": {"Specified": {"Values": observe}}},
            "Events": events,
        },
        {
            "Trigger": {
                "Times": {
                    "Specified": {
                        "Values": (
                            [0.0, args.end_time]
                            if args.endpoints_only
                            else observe
                        )
                    }
                }
            },
            "Events": [fields_event],
        },
        {
            "Trigger": {"Times": {"Specified": {"Values": [args.end_time]}}},
            "Events": ["Completion"],
        },
        {
            "Trigger": "Always",
            "Events": [
                {
                    "ObserveTimeStep": {
                        "SubfileName": "TimeSteps",
                        "PrintTimeToTerminal": False,
                        "ObservePerCore": False,
                    }
                },
                {
                    "ChangeSlabSize": {
                        "DelayChange": 0,
                        "StepChoosers": [
                            (
                                {"Cfl": {"SafetyFactor": args.cfl}}
                                if args.stepping == "cfl"
                                else {"Constant": args.dt}
                            ),
                            {
                                "StepToTimes": {
                                    "Times": {
                                        "Specified": {"Values": times[1:]}
                                    }
                                }
                            },
                        ],
                    }
                },
            ],
        },
    ]
    return config


def run(
    args, *, config=None, analyzer=None, executable_name="EvolveScalarWave3D"
):
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    config = configuration(args) if config is None else config
    (directory / "Input.yaml").write_text(
        f"Executable: {executable_name}\n---\n"
        + yaml.safe_dump(config, sort_keys=False)
    )
    command = [
        str(args.executable.resolve()),
        f"+p{getattr(args, 'processes', 1)}",
        "--input-file",
        "Input.yaml",
    ]
    executable_hash = hashlib.sha256(args.executable.read_bytes()).hexdigest()
    start = time.perf_counter()
    timed_out = False
    with (directory / "Run.log").open("w") as log:
        try:
            result = subprocess.run(
                command,
                cwd=directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                env={
                    **os.environ,
                    "OPENBLAS_NUM_THREADS": "1",
                    "OMP_NUM_THREADS": "1",
                },
            )
            returncode = result.returncode
        except subprocess.TimeoutExpired:
            returncode = -1
            timed_out = True
    metadata = {
        key: value
        for key, value in vars(args).items()
        if key not in ["func", "output", "command", "executable"]
    }
    metadata.update(
        returncode=returncode,
        timed_out=timed_out,
        process_seconds=time.perf_counter() - start,
        executable=str(args.executable.resolve()),
        estimator_variant=(
            "one-dimensional-stripe"
            if executable_name == "EvolveScalarWave1D"
            else (
                "directional-energy"
                if getattr(args, "aggregation", "energy") == "energy"
                else "variation-budget"
            )
        ),
        executable_sha256=executable_hash,
        processes=getattr(args, "processes", 1),
    )
    if hasattr(args, "estimator_variant"):
        metadata["estimator_variant"] = args.estimator_variant
    (directory / "Run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"run": str(directory), **metadata}), flush=True)
    if returncode:
        print((directory / "Run.log").read_text()[-5000:], flush=True)
    else:
        print(json.dumps((analyzer or analyze_run)(directory)), flush=True)
    return returncode


@lru_cache(None)
def quadrature(order):
    xi, weights = leggauss(order)
    return (
        xi,
        weights,
        weights[:, None, None]
        * weights[None, :, None]
        * weights[None, None, :],
    )


@lru_cache(None)
def evaluation_matrix(extent, order):
    return legvander(quadrature(order)[0], extent - 1)


def evaluate(modes, order):
    result = modes
    for axis, extent in enumerate(modes.shape):
        result = np.moveaxis(
            np.tensordot(
                evaluation_matrix(extent, order), result, axes=(1, axis)
            ),
            0,
            axis,
        )
    return result


def elements(group):
    extents = np.asarray(group["total_extents"]).reshape(-1, 3)
    names = bytes(group["grid_names"][:]).decode().split(":")
    offset = 0
    for name, ext in zip(names, extents):
        size = int(np.prod(ext))
        selection = slice(offset, offset + size)
        xyz = [
            np.asarray(group["InertialCoordinates_" + d][selection]).reshape(
                ext, order="F"
            )
            for d in "xyz"
        ]
        bounds = [(float(np.min(x)), float(np.max(x))) for x in xyz]
        fields = np.array(
            [
                np.asarray(group[field][selection]).reshape(ext, order="F")
                for field in FIELDS
            ]
        )
        yield name, ext, bounds, xyz, fields
        offset += size


def element_errors(modes, bounds, t, orientation, order):
    xi, _, weights = quadrature(order)
    xyz = np.meshgrid(
        *[(a + b) / 2 + (b - a) / 2 * xi for a, b in bounds], indexing="ij"
    )
    numerical = np.array([evaluate(c, order) for c in modes])
    jacobian = np.prod([(b - a) / 2 for a, b in bounds])
    squared = jacobian * np.sum(
        weights * (numerical - exact(xyz, t, orientation)) ** 2, axis=(1, 2, 3)
    )
    constraint = np.array(
        [
            evaluate(legder(modes[0], axis=d), order)
            * 2
            / (bounds[d][1] - bounds[d][0])
            - numerical[d + 2]
            for d in range(3)
        ]
    )
    constraint_sq = jacobian * np.sum(weights * constraint**2)
    return squared, constraint_sq


def analyze_run(directory, *, diagnostics=True):
    metadata = json.loads((directory / "Run.json").read_text())
    result = {"run": directory.name, **metadata}
    if metadata["returncode"]:
        reason = re.search(
            r"Tried refining beyond[^\n]+|Chalmers fit is undetermined[^\n]+",
            (directory / "Run.log").read_text(),
        )
        result["failure_reason"] = reason[0] if reason else "See Run.log"
        if (directory / "Reductions.h5").exists():
            with h5py.File(directory / "Reductions.h5") as reductions:
                if "TimeSteps.dat" in reductions and len(
                    reductions["TimeSteps.dat"]
                ):
                    result["last_time"] = float(
                        reductions["TimeSteps.dat"][-1, 0]
                    )
                if "Amr.dat" in reductions and len(reductions["Amr.dat"]):
                    row = reductions["Amr.dat"][-1]
                    result.update(
                        last_mesh_time=float(row[0]),
                        last_elements=int(row[1]),
                        last_points=int(row[2]),
                    )
        (directory / "Analysis.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result
    if diagnostics and metadata["estimator_variant"] not in [
        "directional-energy",
        "variation-budget",
    ]:
        raise ValueError(
            "This run uses the earlier whole-stripe guard. Keep its saved "
            "diagnostics or use the original analysis version; current fits "
            "must not silently overwrite the earlier experiment."
        )
    with h5py.File(directory / "Reductions.h5") as f:
        stats, times = f["Amr.dat"][:], f["TimeSteps.dat"][:]
        result.update(
            evolution_seconds=float(times[-1, -1] - times[0, -1]),
            steps=len(times) - 1,
            minimum_dt=float(np.min(times[:-1, 2])),
            maximum_dt=float(np.max(times[:-1, 2])),
            final_elements=int(stats[-1, 1]),
            final_points=int(stats[-1, 2]),
        )
    rows, fits, spectra, meshes = [], [], [], []
    min_spacing = float("inf")
    with h5py.File(directory / "Volume0.h5") as f:
        groups = sorted(
            (g for g in f["Fields.vol"].values() if isinstance(g, h5py.Group)),
            key=lambda g: g.attrs["observation_value"],
        )
        for group in groups:
            t = float(group.attrs["observation_value"])
            sums = {order: np.zeros(5) for order in [24, 32]}
            constraint = 0.0
            npoints, nelements = 0, 0
            for name, ext, bounds, xyz, fields in elements(group):
                npoints += int(np.prod(ext))
                nelements += 1
                modes = [to_modal(v) for v in fields]
                for order in sums:
                    error, cons = element_errors(
                        modes, bounds, t, metadata["orientation"], order
                    )
                    sums[order] += error
                    if order == 32:
                        constraint += cons
                for d, (a, b) in enumerate(bounds):
                    min_spacing = min(
                        min_spacing,
                        (b - a)
                        / 2
                        * np.min(np.diff(legendre_grid(int(ext[d]))[0])),
                    )
                meshes.append(
                    {
                        "time": t,
                        "element": name,
                        **{f"degree_{d}": int(ext[d] - 1) for d in range(3)},
                        **{f"left_{d}": bounds[d][0] for d in range(3)},
                        **{f"right_{d}": bounds[d][1] for d in range(3)},
                    }
                )
                if diagnostics:
                    mesh = Mesh[3](
                        ext.tolist(), Basis.Legendre, Quadrature.GaussLobatto
                    )
                    for field, v in zip(FIELDS, fields):
                        amplitude = float(np.max(abs(v)))
                        for d, fit in enumerate(
                            (
                                chalmers_energy_estimate
                                if metadata["estimator_variant"]
                                == "directional-energy"
                                else chalmers_estimate
                            )(
                                DataVector(v.ravel(order="F")),
                                mesh,
                                metadata["tolerance"],
                            )
                        ):
                            fits.append(
                                {
                                    "time": t,
                                    "element": name,
                                    "field": field,
                                    "direction": d,
                                    "amplitude": amplitude,
                                    "component_skipped": (
                                        np.sqrt(8) * amplitude
                                        <= metadata["tolerance"]
                                    ),
                                    **fit,
                                }
                            )
                    analytic_modes = to_modal(
                        exact(xyz, t, metadata["orientation"])[1]
                    )
                    for d in range(3):
                        other = tuple(k for k in range(3) if k != d)
                        # Directional maximum over transverse stripes, plus
                        # numerical modal error on the identical mesh.
                        for k, (a, b, c) in enumerate(
                            zip(
                                np.max(abs(modes[1]), axis=other),
                                np.max(abs(analytic_modes), axis=other),
                                np.max(
                                    abs(modes[1] - analytic_modes), axis=other
                                ),
                            )
                        ):
                            spectra.append(
                                {
                                    "time": t,
                                    "element": name,
                                    "direction": d,
                                    "mode": k,
                                    "amplitude": a,
                                    "analytic_amplitude": b,
                                    "error_amplitude": c,
                                }
                            )
            errors = np.sqrt(sums[32])
            delta = float(np.max(abs(errors - np.sqrt(sums[24]))))
            if delta > max(1e-12, 1e-5 * np.max(errors)):
                raise RuntimeError(
                    f"Quadrature not converged at {t} in {directory}: {delta}"
                )
            rows.append(
                {
                    "time": t,
                    "elements": nelements,
                    "points": npoints,
                    **{
                        f"error_{field}": float(e)
                        for field, e in zip(FIELDS, errors)
                    },
                    "error_state": float(np.linalg.norm(errors)),
                    "constraint": float(np.sqrt(constraint)),
                    "quadrature_change": delta,
                }
            )
    if not np.isclose(
        rows[-1]["time"], metadata["end_time"], rtol=0, atol=1e-12
    ):
        raise RuntimeError("Missing final observation")
    result.update(
        initial_error=rows[0]["error_state"],
        final_error=rows[-1]["error_state"],
        final_constraint=rows[-1]["constraint"],
        max_error=max(r["error_state"] for r in rows),
        quadrature_change=max(r["quadrature_change"] for r in rows),
        min_nodal_spacing=min_spacing,
    )
    final_mesh = [r for r in meshes if r["time"] == rows[-1]["time"]]
    result["final_max_degree"] = [
        max(r[f"degree_{d}"] for r in final_mesh) for d in range(3)
    ]
    result["final_min_width"] = [
        min(r[f"right_{d}"] - r[f"left_{d}"] for r in final_mesh)
        for d in range(3)
    ]
    write_csv(directory / "History.csv", rows)
    write_csv(directory / "Mesh.csv", meshes)
    write_csv(directory / "Fits.csv", fits)
    write_csv(directory / "Spectra.csv", spectra)
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def analyze(args):
    results = [analyze_run(p) for p in args.runs]
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "Summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    for result in results:
        print(json.dumps(result), flush=True)


def polynomial_difference(a, b):
    """Integrate the difference of two piecewise polynomials on affine boxes.

    The common partition consists of all nonempty box intersections. Gauss
    quadrature with p+1 points integrates the squared degree-p difference
    exactly, including when the two adaptive h meshes differ.
    """
    squared = np.zeros(len(FIELDS))
    covered_volume = 0.0
    for bounds_a, modes_a in a.values():
        for bounds_b, modes_b in b.values():
            overlap = [
                (max(aa, ba), min(ab, bb))
                for (aa, ab), (ba, bb) in zip(bounds_a, bounds_b)
            ]
            if any(hi <= lo for lo, hi in overlap):
                continue
            ext = tuple(
                max(na, nb)
                for na, nb in zip(modes_a[0].shape, modes_b[0].shape)
            )
            grids = [leggauss(n) for n in ext]
            weights = np.einsum("i,j,k->ijk", *(g[1] for g in grids))
            jacobian = np.prod([(hi - lo) / 2 for lo, hi in overlap])
            covered_volume += 8 * jacobian

            def restricted(modes, bounds):
                result = modes
                for d, ((lo, hi), (a, b), (xi, _)) in enumerate(
                    zip(overlap, bounds, grids)
                ):
                    physical = (lo + hi) / 2 + (hi - lo) / 2 * xi
                    logical = (2 * physical - a - b) / (b - a)
                    matrix = legvander(logical, modes.shape[d] - 1)
                    result = np.moveaxis(
                        np.tensordot(matrix, result, axes=(1, d)), 0, d
                    )
                return result

            for field, (ma, mb) in enumerate(zip(modes_a, modes_b)):
                delta = restricted(ma, bounds_a) - restricted(mb, bounds_b)
                squared[field] += jacobian * np.sum(weights * delta**2)
    volumes = [
        sum(
            np.prod([hi - lo for lo, hi in bounds])
            for bounds, _ in mesh.values()
        )
        for mesh in [a, b]
    ]
    if not np.allclose(volumes, covered_volume, rtol=1e-12, atol=1e-12):
        raise ValueError("Meshes do not cover the same volume")
    return np.sqrt(squared)


def compare(args):
    """Exact polynomial L2 difference, including different adaptive meshes."""

    def final(directory):
        with h5py.File(directory / "Volume0.h5") as f:
            group = max(
                (
                    g
                    for g in f["Fields.vol"].values()
                    if isinstance(g, h5py.Group)
                ),
                key=lambda g: g.attrs["observation_value"],
            )
            return float(group.attrs["observation_value"]), {
                name: (bounds, [to_modal(v) for v in fields])
                for name, _, bounds, _, fields in elements(group)
            }

    ta, a = final(args.first)
    tb, b = final(args.second)
    if abs(ta - tb) > 1e-12:
        raise ValueError("Final times differ")
    differences = polynomial_difference(a, b)
    result = {
        "first": str(args.first),
        "second": str(args.second),
        "same_h_mesh": (
            a.keys() == b.keys() and all(a[name][0] == b[name][0] for name in a)
        ),
        "field_l2_differences": differences.tolist(),
        "state_l2_difference": float(np.linalg.norm(differences)),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


def summarize(args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    results = []
    for directory in args.runs:
        path = directory / "Analysis.json"
        result = (
            json.loads(path.read_text())
            if path.exists()
            else analyze_run(directory)
        )
        result["directory"] = str(directory.resolve())
        results.append(result)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "Summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    colors = {
        "fixed": "#555555",
        "p": "#2271b2",
        "hp1": "#d55e00",
        "hp05": "#009e73",
    }
    labels = {
        "fixed": "Fixed grids",
        "p": "p only",
        "hp1": "hp, threshold 1",
        "hp05": "hp, threshold 0.5",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for column, orientation in enumerate(["axis", "oblique"]):
        for strategy in colors:
            rows = [
                r
                for r in results
                if r["orientation"] == orientation
                and r["strategy"] == strategy
                and "final_error" in r
                and (orientation != "axis" or r.get("axis_periodic", False))
            ]
            for row_index, stepping in enumerate(["fixed", "cfl"]):
                selected = sorted(
                    (
                        r
                        for r in rows
                        if r["stepping"] == stepping
                        and (stepping != "cfl" or r["cfl"] == 0.5)
                    ),
                    key=lambda r: r["max_error"],
                )
                if selected:
                    axes[row_index, column].semilogy(
                        [r["evolution_seconds"] for r in selected],
                        [r["max_error"] for r in selected],
                        "o-",
                        color=colors[strategy],
                        label=labels[strategy],
                    )
        for row_index, stepping in enumerate(
            ["Fixed step schedule", "CFL factor 0.5"]
        ):
            ax = axes[row_index, column]
            ax.set(
                title=f"{orientation.capitalize()} wave: {stepping}",
                xlabel="Evolution wall time [s]",
                ylabel="Maximum observed combined L2 error",
            )
            ax.grid(alpha=0.25)
            ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
            ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:g}"))
            ax.xaxis.set_minor_locator(ticker.NullLocator())
            ax.legend(fontsize=8)
            failures = [
                r
                for r in results
                if r["orientation"] == orientation
                and r["stepping"] == ("fixed" if row_index == 0 else "cfl")
                and (row_index == 0 or r["cfl"] == 0.5)
                and (orientation != "axis" or r.get("axis_periodic", False))
                and r.get("failure_reason", "").startswith(
                    ("Tried refining", "Chalmers fit")
                )
            ]
            if failures:
                ax.text(
                    0.03,
                    0.03,
                    f"{len(failures)} runs stopped: refinement limit or"
                    " undetermined fit",
                    transform=ax.transAxes,
                    fontsize=8,
                    color="#a04000",
                )
    fig.savefig(args.output / "AccuracyCost.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for column, orientation in enumerate(["axis", "oblique"]):
        choices = [
            r
            for r in results
            if r["orientation"] == orientation
            and r["strategy"] == "hp05"
            and r["tolerance"] == 1e-4
            and r["stepping"] == "fixed"
            and (Path(r["directory"]) / "Volume0.h5").exists()
            and (orientation != "axis" or r.get("axis_periodic", False))
        ]
        if not choices:
            axes[0, column].text(
                0.5,
                0.5,
                "Run did not complete",
                ha="center",
                transform=axes[0, column].transAxes,
            )
            continue
        choice = choices[0]
        directory = Path(choice["directory"])
        mesh_path = directory / "Mesh.csv"
        if "final_error" not in choice:
            # Preserve partial mesh histories without giving a failed run a
            # final error or including it in the successful cost comparisons.
            partial = []
            with h5py.File(directory / "Volume0.h5") as f:
                for group in f["Fields.vol"].values():
                    if not isinstance(group, h5py.Group):
                        continue
                    for name, ext, bounds, _, _ in elements(group):
                        partial.append(
                            {
                                "time": float(group.attrs["observation_value"]),
                                "element": name,
                                **{
                                    f"degree_{d}": int(ext[d] - 1)
                                    for d in range(3)
                                },
                                **{f"left_{d}": bounds[d][0] for d in range(3)},
                                **{
                                    f"right_{d}": bounds[d][1] for d in range(3)
                                },
                            }
                        )
            mesh_path = directory / "PartialMesh.csv"
            write_csv(mesh_path, partial)
        with mesh_path.open() as stream:
            mesh = list(csv.DictReader(stream))
        times = sorted(set(float(r["time"]) for r in mesh))
        widths, degrees = [], []
        for t in times:
            cells = [
                r
                for r in mesh
                if float(r["time"]) == t
                and (
                    float(r["left_0"]) >= 0
                    if orientation == "axis"
                    else all(float(r[f"left_{d}"]) >= 0 for d in range(3))
                )
            ]
            volumes = np.array(
                [
                    np.prod(
                        [
                            float(r[f"right_{d}"]) - float(r[f"left_{d}"])
                            for d in range(3)
                        ]
                    )
                    for r in cells
                ]
            )
            widths.append(
                [
                    min(
                        float(r[f"right_{d}"]) - float(r[f"left_{d}"])
                        for r in cells
                    )
                    for d in range(3)
                ]
            )
            degrees.append(
                [
                    np.average(
                        [float(r[f"degree_{d}"]) for r in cells],
                        weights=volumes,
                    )
                    for d in range(3)
                ]
            )
        for d, style in enumerate(["-", "--", ":"]):
            axes[0, column].plot(
                times, np.array(widths)[:, d], style, marker="o", label="xyz"[d]
            )
            axes[1, column].plot(
                times,
                np.array(degrees)[:, d],
                style,
                marker="o",
                label="xyz"[d],
            )
        axes[0, column].set(
            title=f"{orientation.capitalize()}: hp threshold 0.5, target 1e-4",
            ylabel="Minimum element width",
        )
        axes[1, column].set(
            xlabel="Time", ylabel="Volume-weighted mean polynomial degree"
        )
        axes[0, column].set_ylim(0, 1.05)
        axes[1, column].set_ylim(0, max(6, np.max(degrees) + 1))
        for ax in axes[:, column]:
            ax.grid(alpha=0.25)
            ax.legend()
            ax.set_xlim(0, choice["end_time"])
            if "final_error" not in choice:
                ax.axvline(choice["last_time"], color="#d55e00", alpha=0.7)
                ax.text(
                    0.98,
                    0.92,
                    f"Stopped at t≈{choice['last_time']:.2f}: refinement cap",
                    transform=ax.transAxes,
                    ha="right",
                    fontsize=8,
                    color="#a04000",
                )
    fig.savefig(args.output / "DirectionalRefinement.png", dpi=180)
    plt.close(fig)

    # Compare only actually attained accuracies, without interpolating across
    # tolerance settings or counting incomplete runs as successful.
    matched = []
    for orientation, stepping, target, strategy in product(
        ["axis", "oblique"], ["fixed", "cfl"], [1e-3, 1e-4, 1e-5], colors
    ):
        candidates = [
            r
            for r in results
            if "max_error" in r
            and r["orientation"] == orientation
            and r["stepping"] == stepping
            and r["strategy"] == strategy
            and r["max_error"] <= target
            and (stepping != "cfl" or r["cfl"] == 0.5)
            and (orientation != "axis" or r.get("axis_periodic", False))
        ]
        best = (
            min(candidates, key=lambda r: r["evolution_seconds"])
            if candidates
            else None
        )
        matched.append(
            {
                "orientation": orientation,
                "stepping": stepping,
                "required_max_observed_error": target,
                "strategy": strategy,
                "seconds": (
                    best["evolution_seconds"] if best else "not attained"
                ),
                "case": best["run"] if best else "",
            }
        )
    write_csv(args.output / "MatchedAccuracy.csv", matched)
    print(
        f"Wrote {len(results)} run records and comparison figures to"
        f" {args.output}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    runner = commands.add_parser("run")
    runner.add_argument("--executable", type=Path, required=True)
    runner.add_argument("--output", type=Path, required=True)
    runner.add_argument(
        "--orientation", choices=["axis", "oblique"], required=True
    )
    runner.add_argument(
        "--strategy", choices=["fixed", "p", "hp1", "hp05"], required=True
    )
    runner.add_argument(
        "--aggregation", choices=["energy", "stripes"], default="energy"
    )
    runner.add_argument("--tolerance", type=float, default=1e-5)
    runner.add_argument("--dt", type=float, default=0.0005)
    runner.add_argument("--end-time", type=float, default=2.0)
    runner.add_argument("--fixed-p", type=int, default=19)
    runner.add_argument("--max-level", type=int, default=3)
    runner.add_argument("--stepping", choices=["fixed", "cfl"], default="fixed")
    runner.add_argument("--cfl", type=float, default=0.5)
    runner.add_argument("--endpoints-only", action="store_true")
    runner.add_argument(
        "--axis-periodic",
        action="store_true",
        help="Periodic transverse boundaries for the axis symmetry control",
    )
    runner.add_argument(
        "--observation-interval", type=float, default=OUTPUT_INTERVAL
    )
    runner.add_argument("--timeout", type=float, default=900.0)
    runner.set_defaults(func=run)
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("runs", nargs="+", type=Path)
    analyzer.add_argument("--output", type=Path, required=True)
    analyzer.set_defaults(func=analyze)
    differ = commands.add_parser("compare")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.set_defaults(func=compare)
    summary = commands.add_parser("summarize")
    summary.add_argument("runs", nargs="+", type=Path)
    summary.add_argument("--output", type=Path, required=True)
    summary.set_defaults(func=summarize)
    args = parser.parse_args()
    result = args.func(args)
    if isinstance(result, int) and result:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
