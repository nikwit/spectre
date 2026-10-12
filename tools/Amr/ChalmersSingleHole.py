# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Stationary mass-one Schwarzschild control for Chalmers AMR in GH.

Uses Cartesian Kerr-Schild fields on static equiangular wedges, 1.8 <= r <= 6.
Diagnostics use coordinate-volume L2 norms with full symmetric-tensor
multiplicities, not invariant spacetime norms or proper-volume weighting.
"""

import argparse
import csv
import json
from itertools import combinations_with_replacement
from pathlib import Path

import h5py
import numpy as np
import yaml
from numpy.polynomial.legendre import legder, legval3d

from spectre.Domain import ElementId, ElementMap
from spectre.Informer import unit_test_src_path
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv

PAIRS = list(combinations_with_replacement(range(4), 2))
SUFFIXES = ["txyz"[b] + "txyz"[a] for a, b in PAIRS]
FIELDS = (
    [f"SpacetimeMetric_{s}" for s in SUFFIXES]
    + [f"Pi_{s}" for s in SUFFIXES]
    + [f"Phi_{d}{s}" for d in "xyz" for s in SUFFIXES]
)
MULTIPLICITY = np.array([1 if a == b else 2 for a, b in PAIRS] * 5)
ORDERS = (20, 28)
EXPERIMENT = "stationary-schwarzschild-v1"


def exact(xyz):
    """Analytic g_ab, Pi_ab and Phi_iab, in FIELDS order, for M=1, a=0."""
    xyz = np.asarray(xyz)
    r = np.sqrt(np.sum(xyz**2, axis=0))
    n = xyz / r
    ell = np.array([np.ones_like(r), *n])
    dl = np.zeros((3, 4, *r.shape))
    for i in range(3):
        for j in range(3):
            dl[i, j + 1] = ((i == j) - n[i] * n[j]) / r
    metric = np.array(
        [
            (-1 if a == 0 else 1) * (a == b) + 2 / r * ell[a] * ell[b]
            for a, b in PAIRS
        ]
    )
    phi = np.array(
        [
            [
                -2 / r**2 * n[i] * ell[a] * ell[b]
                + 2 / r * (dl[i, a] * ell[b] + ell[a] * dl[i, b])
                for a, b in PAIRS
            ]
            for i in range(3)
        ]
    )
    lapse = 1 / np.sqrt(1 + 2 / r)
    shift = 2 / (r + 2) * n
    pi = np.einsum("i...,ij...->j...", shift, phi) / lapse
    return np.concatenate([metric, pi, phi.reshape((30, *r.shape))])


def configuration(args):
    template = (
        Path(unit_test_src_path()).parent
        / "InputFiles"
        / "GeneralizedHarmonic"
        / "KerrSchild.yaml"
    )
    config = list(yaml.safe_load_all(template.read_text()))[1]
    config["Parallelization"]["ElementDistribution"] = "NumGridPoints"
    config["Evolution"].update(
        InitialTimeStep=0.0005,
        InitialSlabSize=0.0005,
        LocalTimeStepping="Off",
        LtsStepChoosers=None,
        TimeStepper={"DormandPrince5": None},
        VariableOrderAlgorithm={"GoalOrder": 5},
    )
    creator = config["DomainCreator"]["Sphere"]
    creator.update(
        InnerRadius=1.8,
        OuterRadius=6.0,
        RadialPartitioning=[3.0],
        RadialDistribution="Linear",
        InitialGridPoints=getattr(args, "grid_points", None) or args.points,
        InitialRefinement=[0, 0, getattr(args, "radial_level", 0)],
        Interior={"ExciseWithBoundaryCondition": "DemandOutgoingCharSpeeds"},
    )
    if getattr(args, "filter", "none") == "none":
        config["SpatialDiscretization"]["DiscontinuousGalerkin"][
            "Filtering"
        ] = [{"None": {"BlocksToFilter": "All"}}]
    elif args.filter != "hypercube":
        raise ValueError("Unknown filter")
    else:
        config["SpatialDiscretization"]["DiscontinuousGalerkin"]["Filtering"][
            0
        ]["Hypercube"]["HalfPower"] = getattr(args, "half_power", 64)
    criterion_variant = getattr(args, "criterion", "legacy")
    if (
        criterion_variant != "legacy"
        and getattr(args, "filter", "none") != "hypercube"
    ):
        raise ValueError(
            "Ignoring the highest mode requires the Hypercube filter"
        )
    criterion_name = {
        "legacy": "Chalmers",
        "filtered": "ChalmersIgnoringHighestMode",
        "parity": "ChalmersParityIgnoringHighestMode",
    }[criterion_variant]
    criteria = []
    if args.strategy == "split":
        criteria = [
            {
                "DriveToTargetRefinementLevels": {
                    "Target": [1, 0, 1],
                    "OscillationAtTarget": ["DoNothing"] * 3,
                }
            }
        ]
    elif args.strategy != "fixed":
        criteria = [
            {
                criterion_name: {
                    "VariablesToMonitor": ["SpacetimeMetric", "Pi", "Phi"],
                    "AbsoluteTarget": args.tolerance,
                    "RelativeTarget": 0.0,
                    "SmoothnessThreshold": (
                        0.5 if args.strategy == "hp05" else 1.0
                    ),
                    "PRefineOnly": args.strategy == "p",
                    "CoarseningFactor": getattr(args, "coarsening_factor", 0.0),
                    "HCoarsening": getattr(args, "h_coarsening", False),
                    "UseDirectionalEnergy": True,
                }
            }
        ]
    config["Amr"]["Criteria"] = criteria
    config["Amr"]["Policies"].update(
        AllowCoarsening=getattr(args, "coarsening_factor", 0.0) > 0.0
    )
    config["Amr"]["Policies"]["Limits"].update(
        RefinementLevel=[0, 2],
        NumPolynomialModes=[2, 20],
        ErrorBeyondLimits=True,
    )
    for interval in [args.end_time, args.observation_interval]:
        if interval <= 0 or not np.isclose(
            interval / 0.1, round(interval / 0.1)
        ):
            raise ValueError("End time and observations must align with AMR")
    times = np.round(
        np.arange(round(args.end_time / 0.1) + 1) * 0.1, 12
    ).tolist()
    observe = np.round(
        np.arange(round(args.end_time / args.observation_interval) + 1)
        * args.observation_interval,
        12,
    ).tolist()
    if observe[-1] != args.end_time:
        raise ValueError("End time must be divisible by observation interval")
    early_until = getattr(args, "early_observation_until", 0.0)
    if not 0 <= early_until <= args.end_time:
        raise ValueError("Early observation cutoff must lie within the run")
    # Add output on existing AMR boundaries without changing the time steps.
    observe = sorted(set(observe + [t for t in times if t <= early_until]))
    trigger = lambda values: {"Times": {"Specified": {"Values": values}}}
    config["PhaseChangeAndTriggers"] = (
        []
        if args.strategy == "fixed"
        else [
            {
                "Trigger": trigger(times[1:-1]),
                "PhaseChanges": [
                    "VisitAndReturn(EvaluateAmrCriteria)",
                    "VisitAndReturn(AdjustDomain)",
                    "VisitAndReturn(CheckDomain)",
                ],
            }
        ]
    )
    norms = config["EventsAndTriggersAtSlabs"][0]["Events"][1:]
    if getattr(args, "constraint_diagnostics", False):
        for event in norms:
            tensors = event["ObserveNorms"]["TensorsToObserve"]
            tensors.extend(
                {**tensor, "NormType": "Max"}
                for tensor in list(tensors)
                if "Constraint" in tensor["Name"]
            )
    config["EventsAndTriggersAtSlabs"] = [
        {
            "Trigger": trigger(observe),
            "Events": [
                *norms,
                {
                    "ObserveFields": {
                        "SubfileName": "Fields",
                        "VariablesToObserve": ["SpacetimeMetric", "Pi", "Phi"],
                        "InterpolateToMesh": "None",
                        "ProjectToMesh": "None",
                        "CoordinatesFloatingPointType": "Double",
                        "FloatingPointTypes": ["Double"],
                        "BlocksToObserve": "All",
                    }
                },
                "ApparentHorizon",
            ],
        },
        {
            "Trigger": trigger(times),
            "Events": [
                {
                    "ObserveAmrStats": {
                        "PrintToTerminal": False,
                        "ObservePerCore": False,
                    }
                }
            ],
        },
        {"Trigger": trigger([args.end_time]), "Events": ["Completion"]},
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
                            {"Cfl": {"SafetyFactor": args.cfl}},
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
    if getattr(args, "mesh_observations", False):
        # Lapse plus coordinates and element IDs/extents track every mesh
        # transition without saving all 50 GH fields at this cadence.
        config["EventsAndTriggersAtSlabs"][1]["Events"].append(
            {
                "ObserveFields": {
                    "SubfileName": "Meshes",
                    "VariablesToObserve": ["Lapse"],
                    "InterpolateToMesh": "None",
                    "ProjectToMesh": "None",
                    "CoordinatesFloatingPointType": "Double",
                    "FloatingPointTypes": ["Double"],
                    "BlocksToObserve": "All",
                }
            }
        )
    config["Observers"].update(
        VolumeFileName="Volume",
        ReductionFileName="Reductions",
        SurfaceFileName="Surfaces",
    )
    horizon = config["ApparentHorizons"]["ApparentHorizon"]
    horizon["InitialGuess"]["InitialShape"]["Sphere"]["Radius"] = 2.0
    horizon["Criteria"] = [
        {
            "Residual": {
                "MinResidual": 1e-8,
                "MaxResidual": 1e-6,
                "MinResolutionL": 6,
            }
        }
    ]
    horizon["FastFlow"]["TruncationTol"] = 1e-8
    horizon["Verbosity"] = "Silent"
    for name in ["ExcisionBoundary", "ControlSystemCharSpeedExcision"]:
        config["InterpolationTargets"][name]["Radius"] = 1.8
    return config


def elements(group):
    extents = np.asarray(group["total_extents"]).reshape(-1, 3)
    names = bytes(group["grid_names"][:]).decode().split(":")
    offset = 0
    for name, ext in zip(names, extents):
        size = int(np.prod(ext))
        selection = slice(offset, offset + size)
        fields = np.array(
            [
                np.asarray(group[field][selection]).reshape(ext, order="F")
                for field in FIELDS
            ]
        )
        xyz = np.array(
            [
                np.asarray(
                    group[f"InertialCoordinates_{d}"][selection]
                ).reshape(ext, order="F")
                for d in "xyz"
            ]
        )
        yield name, ext, xyz, fields
        offset += size


def element_diagnostics(modes, element_map, order):
    xi, _, weights = cube.quadrature(order)
    xyz, determinant, inverse = sphere.geometry(element_map, [xi] * 3)
    measure = weights * determinant
    numerical = np.array([cube.evaluate(c, order) for c in modes])
    squared = np.sum(measure * (numerical - exact(xyz)) ** 2, axis=(1, 2, 3))
    errors = [
        np.sum((squared * MULTIPLICITY)[s])
        for s in [slice(0, 10), slice(10, 20), slice(20, 50)]
    ]
    c3 = np.array(
        [sphere.physical_gradient(c, inverse, order) for c in modes[:10]]
    ) - numerical[20:].reshape(3, 10, *xyz.shape[1:]).swapaxes(0, 1)
    constraint = np.sum(
        measure * c3**2 * MULTIPLICITY[:10, None, None, None, None]
    )
    return np.array([*errors, constraint]), float(np.sum(measure))


def analyze_run(directory, allow_incomplete=False):
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata.get("experiment") != EXPERIMENT:
        raise ValueError("Not a stationary Schwarzschild experiment")
    result = {"run": directory.name, "processes": 1, **metadata}
    if metadata["returncode"] and not allow_incomplete:
        result["failure_log"] = str(directory / "Run.log")
        (directory / "Analysis.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result
    result["incomplete"] = metadata["returncode"] != 0
    domain = sphere.read_domain(directory)
    rows, meshes = [], []
    initial_nodal_error = initial_metric_error = initial_phi_constraint = (
        coordinate_error
    ) = 0.0
    with h5py.File(directory / "Volume0.h5") as handle:
        for group in sphere.groups(handle):
            t = float(group.attrs["observation_value"])
            sums = {order: np.zeros(4) for order in ORDERS}
            volume = 0.0
            for name, ext, xyz, fields in elements(group):
                element_id = ElementId[3](name)
                element_map = ElementMap(element_id, domain)
                mapped, _, inverse = sphere.geometry(
                    element_map, [legendre_grid(n)[0] for n in ext]
                )
                coordinate_error = max(
                    coordinate_error, float(np.max(abs(mapped - xyz)))
                )
                if t == 0:
                    initial_nodal_error = max(
                        initial_nodal_error,
                        float(np.max(abs(fields - exact(xyz)))),
                    )
                modes = [to_modal(v) for v in fields]
                if t == 0:
                    initial_metric_error = max(
                        initial_metric_error,
                        float(np.max(abs(fields[:10] - exact(xyz)[:10]))),
                    )
                    logical = np.meshgrid(
                        *[legendre_grid(n)[0] for n in ext], indexing="ij"
                    )
                    for k, mode in enumerate(modes[:10]):
                        gradient = np.einsum(
                            "ai...,a...->i...",
                            inverse,
                            np.array(
                                [
                                    legval3d(*logical, legder(mode, axis=d))
                                    for d in range(3)
                                ]
                            ),
                        )
                        initial_phi_constraint = max(
                            initial_phi_constraint,
                            float(
                                np.max(
                                    abs(
                                        gradient
                                        - fields[20:].reshape(3, 10, *ext)[:, k]
                                    )
                                )
                            ),
                        )
                for order in ORDERS:
                    contributions, element_volume = element_diagnostics(
                        modes, element_map, order
                    )
                    sums[order] += contributions
                volume += element_volume
                meshes.append(
                    {
                        "time": t,
                        "element": name,
                        "points": int(np.prod(ext)),
                        **{
                            f"degree_{d}": int(n - 1) for d, n in enumerate(ext)
                        },
                        **{
                            f"level_{d}": int(n)
                            for d, n in enumerate(element_id.refinement_levels)
                        },
                    }
                )
            low, high = [np.sqrt(sums[order]) for order in ORDERS]
            delta = float(np.max(abs(high - low)))
            if not np.all(np.isfinite(high)) or delta > max(
                1e-11, 1e-5 * max(high)
            ):
                raise RuntimeError(f"Unconverged diagnostics at {t}: {delta}")
            if (
                abs(volume - 4 * np.pi / 3 * (6**3 - 1.8**3)) > 1e-9
                or coordinate_error > 1e-12
            ):
                raise RuntimeError("Geometry verification failed")
            rows.append(
                {
                    "time": t,
                    "error_metric": float(high[0]),
                    "error_pi": float(high[1]),
                    "error_phi": float(high[2]),
                    "error_state": float(np.linalg.norm(high[:3])),
                    "three_index_constraint": float(high[3]),
                    "quadrature_change": delta,
                }
            )
    with h5py.File(directory / "Reductions.h5") as handle:
        for subfile in ["Errors", "Norms", "ApparentHorizon"]:
            data = handle[subfile + ".dat"]
            legend = [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in data.attrs["Legend"]
            ]
            write_csv(
                directory / (subfile + ".csv"),
                [dict(zip(legend, row)) for row in data[:]],
            )
            columns = dict(zip(legend, data[:].T))
            if not np.all(np.isfinite(data[:])) or not np.array_equal(
                columns["Time"], [r["time"] for r in rows]
            ):
                raise RuntimeError(f"Missing or invalid {subfile} observations")
            if subfile == "ApparentHorizon":
                result.update(
                    max_horizon_mass_error=float(
                        np.max(abs(columns["ChristodoulouMass"] - 1))
                    ),
                    max_horizon_area_relative_error=float(
                        np.max(abs(columns["Area"] / (16 * np.pi) - 1))
                    ),
                    max_dimensionless_spin=float(
                        np.max(columns["DimensionlessSpinMagnitude"])
                    ),
                )
            elif subfile == "Norms":
                result["native_constraint_maxima"] = {
                    key: float(np.max(value))
                    for key, value in columns.items()
                    if "Constraint" in key
                }
        amr = handle["Amr.dat"]
        stats, steps = amr[:], handle["TimeSteps.dat"][:]
        legend = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in amr.attrs["Legend"]
        ]
        write_csv(
            directory / "Amr.csv", [dict(zip(legend, row)) for row in stats]
        )
        deltas = np.diff(stats[:, 2])
        result.update(
            min_points=int(np.min(stats[:, 2])),
            max_points=int(np.max(stats[:, 2])),
            point_decrease_observations=int(np.sum(deltas < 0)),
            point_increase_observations=int(np.sum(deltas > 0)),
            amr_observation_count=len(stats),
            final_elements=int(stats[-1, 1]),
            final_points=int(stats[-1, 2]),
            evolution_seconds=float(steps[-1, -1] - steps[0, -1]),
        )
    result.update(
        last_observation_time=rows[-1]["time"],
        initial_nodal_error=initial_nodal_error,
        initial_metric_nodal_error=initial_metric_error,
        initial_phi_nodal_constraint=initial_phi_constraint,
        coordinate_error=coordinate_error,
        quadrature_orders=list(ORDERS),
        quadrature_change=max(r["quadrature_change"] for r in rows),
        initial_error=rows[0]["error_state"],
        final_error=rows[-1]["error_state"],
        max_error=max(r["error_state"] for r in rows),
        max_three_index_constraint=max(
            r["three_index_constraint"] for r in rows
        ),
        final_max_level=[
            max(
                r[f"level_{d}"] for r in meshes if r["time"] == rows[-1]["time"]
            )
            for d in range(3)
        ],
        final_max_degree=[
            max(
                r[f"degree_{d}"]
                for r in meshes
                if r["time"] == rows[-1]["time"]
            )
            for d in range(3)
        ],
    )
    # GH deliberately resets Pi and Phi to satisfy the discrete constraints.
    # Only the metric agrees with analytic nodal data to roundoff; Phi must
    # agree with the derivative of that interpolant, not the continuum field.
    if max(initial_metric_error, initial_phi_constraint) > 1e-12 or (
        not allow_incomplete and rows[-1]["time"] != metadata["end_time"]
    ):
        raise RuntimeError(
            "Initial data or final observation validation failed"
        )
    write_csv(directory / "History.csv", rows)
    write_csv(directory / "Mesh.csv", meshes)
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def compare(args):
    """Compare full GH states, integrating on a common partition if needed."""

    def history(directory):
        with (directory / "Mesh.csv").open() as handle:
            return sorted(
                tuple(sorted(row.items())) for row in csv.DictReader(handle)
            )

    same_history = history(args.first) == history(args.second)
    domain = sphere.read_domain(args.first)
    if domain != sphere.read_domain(args.second):
        raise ValueError("Domains differ")
    with (
        h5py.File(args.first / "Volume0.h5") as a,
        h5py.File(args.second / "Volume0.h5") as b,
    ):

        def observation(handle):
            groups = sphere.groups(handle)
            if getattr(args, "time", None) is None:
                return groups[-1]
            return next(
                g for g in groups if g.attrs["observation_value"] == args.time
            )

        group_a, group_b = observation(a), observation(b)
        if (
            group_a.attrs["observation_value"]
            != group_b.attrs["observation_value"]
        ):
            raise ValueError("Final times differ")
        fields_a = {
            name: (ext, xyz, fields)
            for name, ext, xyz, fields in elements(group_a)
        }
        fields_b = {
            name: (ext, xyz, fields)
            for name, ext, xyz, fields in elements(group_b)
        }
        same_mesh = fields_a.keys() == fields_b.keys() and all(
            np.array_equal(ext, fields_b[name][0])
            and np.allclose(xyz, fields_b[name][1], atol=1e-13, rtol=0)
            for name, (ext, xyz, _) in fields_a.items()
        )
        sums = np.zeros(2)
        if same_mesh:
            # Subtract nodal fields before the transform to retain exact zeros.
            for name, (_, _, fields) in fields_a.items():
                modes = [to_modal(v) for v in fields - fields_b[name][2]]
                for i, order in enumerate(ORDERS):
                    xi, _, weights = cube.quadrature(order)
                    _, determinant, _ = sphere.geometry(
                        ElementMap(ElementId[3](name), domain), [xi] * 3
                    )
                    sums[i] += sum(
                        w
                        * np.sum(
                            weights * determinant * cube.evaluate(c, order) ** 2
                        )
                        for w, c in zip(MULTIPLICITY, modes)
                    )
        else:
            modes_a = {
                name: [to_modal(v) for v in fields]
                for name, (_, _, fields) in fields_a.items()
            }
            modes_b = {
                name: [to_modal(v) for v in fields]
                for name, (_, _, fields) in fields_b.items()
            }
            for i, order in enumerate(ORDERS):
                differences = sphere.curved_difference(
                    modes_a, modes_b, domain, order
                )
                sums[i] = np.dot(MULTIPLICITY, differences**2)
    low, high = np.sqrt(sums)
    if abs(high - low) > max(1e-12, 1e-5 * high):
        raise RuntimeError("Unconverged GH state comparison")
    # Read the error at the selected observation, which need not be the last.
    with (args.second / "History.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    selected_time = getattr(args, "time", None)
    reference_error = float(
        rows[-1]["error_state"]
        if selected_time is None
        else next(
            row["error_state"]
            for row in rows
            if float(row["time"]) == selected_time
        )
    )
    result = {
        "first": str(args.first),
        "second": str(args.second),
        "same_recorded_meshes": same_history,
        "same_compared_mesh": same_mesh,
        "state_l2_difference": float(high),
        "fraction_of_reference_error": float(high / reference_error),
        "comparison_time": (
            selected_time
            if selected_time is not None
            else float(rows[-1]["time"])
        ),
        "quadrature_change": float(abs(high - low)),
        "quadrature_orders": list(ORDERS),
    }
    if selected_time is None:
        result["fraction_of_final_error"] = result[
            "fraction_of_reference_error"
        ]
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def select_initial_resolution(args):
    """Initialize each candidate natively and measure its complete GH error.

    Selection uses t=0 only, after GH has reset Pi/Phi from the discrete
    constraints. A short fixed-grid evolution also verifies the input. The
    selected mesh is then reinitialized analytically for the comparison runs.
    """
    if not np.isfinite(args.target) or args.target <= 0:
        raise ValueError("Initial error target must be finite and positive")
    args.output.mkdir(parents=True, exist_ok=False)
    candidates = []
    for points in args.candidates:
        probe = argparse.Namespace(
            output=args.output / f"points{points}",
            executable=args.executable,
            strategy="fixed",
            points=points,
            tolerance=1e-5,
            cfl=0.5,
            end_time=0.1,
            observation_interval=0.1,
            timeout=args.timeout,
            experiment=EXPERIMENT,
            processes=args.processes,
            coarsening_factor=0.0,
        )
        status = cube.run(
            probe,
            config=configuration(probe),
            analyzer=analyze_run,
            executable_name="EvolveGhSingleBlackHole",
        )
        if status:
            raise RuntimeError(f"Initial-resolution probe failed: {points}")
        result = json.loads((probe.output / "Analysis.json").read_text())
        candidates.append(
            {
                "points_per_direction": points,
                "initial_error": result["initial_error"],
                "run": str(probe.output),
            }
        )
        if result["initial_error"] <= args.target:
            selection = {
                "target": args.target,
                "selected_points": points,
                "candidates": candidates,
                "norm": "Full GH coordinate-volume L2 at t=0",
            }
            (args.output / "Selection.json").write_text(
                json.dumps(selection, indent=2) + "\n"
            )
            return selection
    (args.output / "Selection.json").write_text(
        json.dumps(
            {
                "target": args.target,
                "selected_points": None,
                "candidates": candidates,
            },
            indent=2,
        )
        + "\n"
    )
    raise RuntimeError("No candidate satisfies the initial error target")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    runner = commands.add_parser("run")
    runner.add_argument("--executable", type=Path, required=True)
    runner.add_argument("--output", type=Path, required=True)
    runner.add_argument(
        "--strategy",
        choices=["fixed", "split", "p", "hp05", "hp1"],
        required=True,
    )
    runner.add_argument("--points", type=int, default=8)
    runner.add_argument(
        "--grid-points",
        type=int,
        nargs=3,
        help="Override --points with angular, angular, radial extents.",
    )
    runner.add_argument("--tolerance", type=float, default=1e-4)
    runner.add_argument("--cfl", type=float, default=0.5)
    runner.add_argument("--end-time", type=float, default=2.0)
    runner.add_argument("--observation-interval", type=float, default=0.5)
    runner.add_argument(
        "--early-observation-until",
        type=float,
        default=0.0,
        help="Also save full diagnostics at every AMR check through this time.",
    )
    runner.add_argument(
        "--constraint-diagnostics",
        action="store_true",
        help="Also record the largest nodal GH constraint norms.",
    )
    runner.add_argument("--timeout", type=float, default=900)
    runner.add_argument("--processes", type=int, default=6)
    runner.add_argument("--coarsening-factor", type=float, default=0.0)
    runner.add_argument("--h-coarsening", action="store_true")
    runner.add_argument("--radial-level", type=int, default=0)
    runner.add_argument("--half-power", type=int, default=64)
    runner.add_argument(
        "--criterion",
        choices=["legacy", "filtered", "parity"],
        default="legacy",
    )
    runner.add_argument("--mesh-observations", action="store_true")
    runner.add_argument("--defer-analysis", action="store_true")
    runner.add_argument(
        "--filter", choices=["none", "hypercube"], default="none"
    )
    runner.set_defaults(
        experiment=EXPERIMENT,
        func=lambda a: cube.run(
            a,
            config=configuration(a),
            analyzer=(
                (lambda _: {"analysis": "deferred"})
                if a.defer_analysis
                else analyze_run
            ),
            executable_name="EvolveGhSingleBlackHole",
        ),
    )
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("run", type=Path)
    analyzer.add_argument("--allow-incomplete", action="store_true")
    analyzer.set_defaults(func=lambda a: analyze_run(a.run, a.allow_incomplete))
    differ = commands.add_parser("compare")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.add_argument("--time", type=float)
    differ.set_defaults(func=compare)
    selector = commands.add_parser("select-initial")
    selector.add_argument("--executable", type=Path, required=True)
    selector.add_argument("--output", type=Path, required=True)
    selector.add_argument("--target", type=float, default=1e-5)
    selector.add_argument(
        "--candidates", type=int, nargs="+", default=[8, 10, 12, 14, 16]
    )
    selector.add_argument("--processes", type=int, default=6)
    selector.add_argument("--timeout", type=float, default=900)
    selector.set_defaults(func=select_initial_resolution)
    args = parser.parse_args()
    result = args.func(args)
    if isinstance(result, int):
        if result:
            raise SystemExit(1)
    else:
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
