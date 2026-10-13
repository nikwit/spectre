# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Prepare controlled BBH AMR pilots from a rendered, working inspiral input.

This command writes inputs, not jobs. Generate them on the cluster so input
file paths resolve there. See ChalmersBbh.md for the staged experiment.
"""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import yaml

CRITERION = "SpecSpectral"
CASES = (
    "p-smoke",
    "restart-control",
    "restart-checkpoint",
    "fixed",
    "fixed-fine",
    "spec",
    "spec-tight",
)


def integer_leaves(value):
    if isinstance(value, dict):
        return [n for v in value.values() for n in integer_leaves(v)]
    if isinstance(value, list):
        return [n for v in value for n in integer_leaves(v)]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("Expected integer initial extents/refinement levels")
    return [value]


def add_points(value, increment):
    if isinstance(value, dict):
        return {k: add_points(v, increment) for k, v in value.items()}
    if isinstance(value, list):
        return [add_points(v, increment) for v in value]
    return value + increment


def resolve_input_paths(value, directory):
    """Only input paths are rebased; output paths remain local to each run."""
    if isinstance(value, dict):
        for k, v in value.items():
            if k in ("FileGlob", "H5Filename", "DataDirectory"):
                if isinstance(v, str) and v != "None":
                    value[k] = str((directory / v).resolve())
            else:
                resolve_input_paths(v, directory)
    elif isinstance(value, list):
        for v in value:
            resolve_input_paths(v, directory)


def trigger(times):
    return {"Times": {"Specified": {"Values": times}}}


def times(initial, duration, interval):
    count = round(duration / interval)
    if count < 1 or not math.isclose(count * interval, duration, abs_tol=1e-10):
        raise ValueError(
            "Durations must be multiples of observation/AMR intervals"
        )
    return [round(initial + n * interval, 12) for n in range(count + 1)]


def observe_fields(subfile, variables):
    return {
        "ObserveFields": {
            "SubfileName": subfile,
            "VariablesToObserve": variables,
            "InterpolateToMesh": "None",
            "ProjectToMesh": "None",
            "CoordinatesFloatingPointType": "Double",
            "FloatingPointTypes": ["Double"],
            "BlocksToObserve": "All",
        }
    }


def configuration(base, case, settings):
    """Return a fresh config; leave initial data, gauge, maps and LTS intact."""
    s = settings
    config = copy.deepcopy(base)
    if set(config["DomainCreator"]) != {"BinaryCompactObject"}:
        raise ValueError(
            "Only the rectangular BinaryCompactObject domain is supported"
        )
    domain = config["DomainCreator"]["BinaryCompactObject"]
    if domain["OuterShell"].get("UseSphericalHarmonics", False):
        raise ValueError("All blocks must use tensor-product Legendre bases")
    if not domain.get("TimeDependentMaps"):
        raise ValueError("Use a working BBH input with its moving maps")
    initial_points = integer_leaves(domain["InitialGridPoints"])
    initial_levels = integer_leaves(domain["InitialRefinement"])
    if min(initial_points) < s.min_points or max(initial_points) > s.max_points:
        raise ValueError(
            "Initial grid must lie inside the requested point limits"
        )
    if min(initial_levels) < s.min_level or max(initial_levels) > s.max_level:
        raise ValueError(
            "Initial h mesh must lie inside the requested h limits"
        )
    filters = config["SpatialDiscretization"]["DiscontinuousGalerkin"][
        "Filtering"
    ]
    if len(filters) != 1 or set(filters[0]) != {"Hypercube"}:
        raise ValueError("Expected one steep Hypercube filter on all blocks")
    filt = filters[0]["Hypercube"]
    if not (
        filt["Enable"]
        and filt["BlocksToFilter"] == "All"
        and filt["VolumeFilterOnSubstep"]
        and filt["HalfPower"] >= 256
    ):
        raise ValueError(
            "Expected the enabled steep substep filter on all blocks"
        )
    if config["Amr"].get("Criteria"):
        raise ValueError(
            "Start from a fixed-mesh BBH input with empty Amr.Criteria"
        )
    initial = float(config["Evolution"]["InitialTime"])
    slab = float(config["Evolution"]["InitialSlabSize"])
    for interval in (s.amr_interval, s.volume_interval):
        if not math.isclose(
            interval / slab, round(interval / slab), abs_tol=1e-9
        ):
            raise ValueError(
                "Diagnostic and AMR intervals must align with initial slabs"
            )
    # Restrict this first pilot to the established constant-slab inspiral phase.
    # Near-merger changes remain in the base input; finish before they activate.
    smoke = case.endswith("smoke") or case.startswith("restart-")
    duration = s.smoke_duration if smoke else s.duration
    end = initial + duration
    amr_times = times(initial, duration, s.amr_interval)
    volume_times = times(initial, duration, s.volume_interval)
    fixed = case.startswith("fixed")
    if case == "fixed-fine":
        domain["InitialGridPoints"] = add_points(
            domain["InitialGridPoints"], s.fine_increment
        )
        if max(integer_leaves(domain["InitialGridPoints"])) > s.max_points:
            raise ValueError(
                "Finer reference exceeds point/filter limit; raise it"
                " explicitly (<=20)"
            )
    tolerance_factor = s.tight_factor if case == "spec-tight" else 1.0
    target = s.spectral_target * tolerance_factor
    # SpEC's background plus two Gaussian profiles, expressed in the grid
    # frame. Centers remain attached to the objects under moving maps.
    y, z = domain["CenterOfMassOffset"]
    gaussians = []
    for name, other in (("ObjectA", "ObjectB"), ("ObjectB", "ObjectA")):
        gaussians.append(
            dict(
                Center=[domain[name]["XCoord"], y, z],
                Width=abs(domain[other]["XCoord"]),
                Amplitude=1.0,
                Target=target * s.near_hole_factor,
            )
        )
    if any(g["Width"] <= 0 for g in gaussians):
        raise ValueError("Gaussian widths require nonzero object X coordinates")
    spectral_target = dict(Background=target, Gaussians=gaussians)
    config["Amr"].update(
        Criteria=(
            None
            if fixed
            else [
                {
                    CRITERION: {
                        "RelativeTruncationTarget": spectral_target,
                        "ConstraintTargets": [
                            c * tolerance_factor for c in s.constraint_targets
                        ],
                        "ForceDropPileUpModes": True,
                        "MinimumModes": s.min_points,
                        "CallsBetweenChanges": 2,
                        "TimeBetweenChanges": s.time_between_changes,
                        "TruncationErrorAdjustmentLimit": -1.0,
                        "StartTime": initial,
                        "Verbose": not s.quiet_decisions,
                    }
                }
            ]
        ),
        EnableInBlocks="All",
        Verbosity="Quiet",
    )
    config["Amr"]["Policies"].update(
        Isotropy="Anisotropic",
        AllowCoarsening=True,
        EnforceTwoToOneBalanceInNormalDirection=True,
    )
    config["Amr"]["Policies"]["Limits"].update(
        RefinementLevel=[s.min_level, s.max_level],
        # SpECTRE's policy specifies polynomial degree, not point count.
        NumPolynomialModes=[s.min_points - 1, s.max_points - 1],
        ErrorBeyondLimits=True,
    )
    phases = []
    for entry in config.get("PhaseChangeAndTriggers") or []:
        changes = entry["PhaseChanges"]
        if any(
            isinstance(c, str) and ("AmrCriteria" in c or "AdjustDomain" in c)
            for c in changes
        ):
            raise ValueError("Base input already has an AMR phase trigger")
        changes = [
            c
            for c in changes
            if not (
                isinstance(c, dict) and "CheckpointAndExitAfterWallclock" in c
            )
        ]
        if changes:
            phases.append({**entry, "PhaseChanges": changes})
    phases.append(
        {
            "Trigger": "Always",
            "PhaseChanges": [
                {
                    "CheckpointAndExitAfterWallclock": {
                        "WallclockHours": s.checkpoint_hours
                    }
                }
            ],
        }
    )
    if not fixed:
        phases.append(
            {
                "Trigger": trigger(amr_times[1:-1]),
                "PhaseChanges": [
                    "VisitAndReturn(EvaluateAmrCriteria)",
                    "VisitAndReturn(AdjustDomain)",
                    "VisitAndReturn(CheckDomain)",
                ],
            }
        )
    # The two restart inputs are otherwise identical. A checkpoint is written
    # without exiting, providing both uninterrupted and restarted trajectories.
    if case == "restart-checkpoint":
        phases.append(
            {
                "Trigger": trigger([initial + s.restart_at]),
                "PhaseChanges": ["VisitAndReturn(WriteCheckpoint)"],
            }
        )
    config["PhaseChangeAndTriggers"] = phases
    events = []
    for entry in config.get("EventsAndTriggersAtSlabs") or []:
        # Replace the finite-duration Completion event, leaving safety checks,
        # horizon/excision observations and all dense control triggers intact.
        kept = [e for e in entry["Events"] if e != "Completion"]
        if kept:
            events.append({**entry, "Events": kept})
    events.extend(
        [
            {
                "Trigger": trigger(amr_times),
                "Events": [
                    {
                        "ObserveAmrStats": {
                            "PrintToTerminal": False,
                            "ObservePerCore": False,
                        }
                    },
                    {
                        "ObserveNorms": {
                            "SubfileName": "AmrConstraints",
                            "TensorsToObserve": [
                                {
                                    "Name": (
                                        f"PointwiseL2Norm({name}Constraint)"
                                    ),
                                    "NormType": norm,
                                    "Components": "Sum",
                                }
                                for name in ("Gauge", "TwoIndex", "ThreeIndex")
                                for norm in ("L2Norm", "Max")
                            ],
                        }
                    },
                ],
            },
            {
                "Trigger": trigger(volume_times),
                "Events": [
                    observe_fields(
                        "AmrFields",
                        [
                            "SpacetimeMetric",
                            "Pi",
                            "Phi",
                            "Lapse",
                            "GaugeH",
                            "PointwiseL2Norm(GaugeConstraint)",
                            "PointwiseL2Norm(TwoIndexConstraint)",
                            "PointwiseL2Norm(ThreeIndexConstraint)",
                        ],
                    )
                ],
            },
            {
                "Trigger": {
                    "TimeCompares": {
                        "Comparison": "GreaterThanOrEqualTo",
                        "Value": end,
                    }
                },
                "Events": ["Completion"],
            },
        ]
    )
    if s.mesh_observations:
        events.append(
            {
                "Trigger": trigger(amr_times),
                "Events": [observe_fields("AmrMeshes", ["Lapse"])],
            }
        )
    config["EventsAndTriggersAtSlabs"] = events
    return config


def prepare(args):
    if args.output.exists():
        raise ValueError("Output exists; choose a fresh directory")
    documents = list(yaml.safe_load_all(args.base_input.read_text()))
    if (
        len(documents) != 2
        or documents[0].get("Executable") != "EvolveGhBinaryBlackHole"
    ):
        raise ValueError(
            "Expected rendered EvolveGhBinaryBlackHole metadata + input"
        )
    metadata, base = documents
    resolve_input_paths(
        base, (args.input_run_dir or args.base_input.parent).resolve()
    )
    cases = {name: configuration(base, name, args) for name in CASES}
    # The two adaptive runs differ only in spectral and constraint targets.
    tighter = copy.deepcopy(cases["spec"])
    criterion = tighter["Amr"]["Criteria"][0][CRITERION]
    profile = criterion["RelativeTruncationTarget"]
    profile["Background"] *= args.tight_factor
    for gaussian in profile["Gaussians"]:
        gaussian["Target"] *= args.tight_factor
    criterion["ConstraintTargets"] = [
        c * args.tight_factor for c in criterion["ConstraintTargets"]
    ]
    actual = cases["spec-tight"]["Amr"]["Criteria"][0][CRITERION]
    assert math.isclose(
        profile["Background"], actual["RelativeTruncationTarget"]["Background"]
    )
    for expected, found in zip(
        profile["Gaussians"], actual["RelativeTruncationTarget"]["Gaussians"]
    ):
        assert math.isclose(expected["Target"], found["Target"])
        expected["Target"] = found["Target"]
    profile["Background"] = actual["RelativeTruncationTarget"]["Background"]
    assert tighter == cases["spec-tight"]
    restart = copy.deepcopy(cases["restart-checkpoint"])
    restart["PhaseChangeAndTriggers"].pop()
    assert restart == cases["restart-control"]
    args.output.mkdir(parents=True)
    outputs = {}
    for name, config in cases.items():
        # Never inherit pipeline continuation, Testing or expected-output rules.
        meta = {
            k: v
            for k, v in metadata.items()
            if k not in ("Next", "Testing", "ExpectedOutput")
        }
        meta["AmrExperiment"] = dict(
            case=name, base_input=str(args.base_input.resolve())
        )
        path = args.output / (name + ".yaml")
        path.write_text(yaml.safe_dump_all([meta, config], sort_keys=False))
        outputs[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = dict(
        base_input=str(args.base_input.resolve()),
        base_sha256=hashlib.sha256(args.base_input.read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        settings={
            k: str(v) if isinstance(v, Path) else v
            for k, v in vars(args).items()
        },
        inputs_sha256=outputs,
        matched_controller_inputs_verified=True,
        restart_inputs_differ_only_by_checkpoint=True,
        validation=(
            "not run; validate with the cluster-built executable before"
            " submission"
        ),
        note=(
            "Native RMS constraints are not physical-volume norms. Targets are"
            " provisional BBH values."
        ),
    )
    (args.output / "Manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(f"Prepared {len(cases)} inputs in {args.output}; no jobs submitted.")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("base_input", type=Path)
    p.add_argument("--input-run-dir", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--duration", type=float, default=200.0)
    p.add_argument("--smoke-duration", type=float, default=20.0)
    p.add_argument(
        "--restart-at", type=float, default=10.0, help="Time after InitialTime"
    )
    p.add_argument("--amr-interval", type=float, default=0.5)
    p.add_argument("--volume-interval", type=float, default=5.0)
    p.add_argument("--checkpoint-hours", type=float, required=True)
    p.add_argument("--min-points", type=int, default=7)
    p.add_argument("--max-points", type=int, default=20)
    p.add_argument("--min-level", type=int, default=1)
    p.add_argument("--max-level", type=int, default=4)
    p.add_argument("--fine-increment", type=int, default=2)
    p.add_argument("--spectral-target", type=float, default=1e-4)
    p.add_argument(
        "--near-hole-factor",
        type=float,
        default=0.01,
        help="Gaussian target/background ratio; center blend is geometric",
    )
    p.add_argument(
        "--tight-factor",
        type=float,
        default=0.25,
        help="Multiply all spectral and constraint targets in spec-tight",
    )
    p.add_argument("--time-between-changes", type=float, default=1.0)
    p.add_argument(
        "--constraint-targets", type=float, nargs=2, default=[0.002, 0.004]
    )
    p.add_argument("--quiet-decisions", action="store_true")
    p.add_argument(
        "--mesh-observations",
        action="store_true",
        help="Save every mesh; expensive for large BBH grids",
    )
    return p


def main():
    p = parser()
    args = p.parse_args()
    positive = [
        args.duration,
        args.smoke_duration,
        args.amr_interval,
        args.volume_interval,
        args.checkpoint_hours,
        args.near_hole_factor,
        args.tight_factor,
        args.spectral_target,
        *args.constraint_targets,
    ]
    if not all(math.isfinite(v) and v > 0 for v in positive):
        p.error("Durations, intervals and targets must be finite and positive")
    if not 0 < args.tight_factor < 1 or not 0 < args.near_hole_factor <= 1:
        p.error("Require 0 < tight-factor < 1 and 0 < near-hole-factor <= 1")
    if (
        not math.isfinite(args.time_between_changes)
        or args.time_between_changes < 0
    ):
        p.error("time-between-changes must be finite and nonnegative")
    if not 7 <= args.min_points <= args.max_points <= 20:
        p.error(
            "Require 7 <= min-points <= max-points <= 20 for the filter model"
        )
    if not 0 <= args.min_level <= args.max_level or args.fine_increment <= 0:
        p.error("Invalid refinement limits or fine increment")
    if not 0 < args.restart_at < args.smoke_duration or not math.isclose(
        args.restart_at / args.amr_interval,
        round(args.restart_at / args.amr_interval),
        abs_tol=1e-9,
    ):
        p.error("restart-at must lie on an interior AMR boundary")
    try:
        prepare(args)
    except ValueError as error:
        p.error(str(error))


if __name__ == "__main__":
    main()
