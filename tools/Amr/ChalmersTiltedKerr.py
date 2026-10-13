# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Tilted Kerr gauge relaxation with matched original/SpEC coarsening runs.

Fixed three-dimensional runs provide references; stationary Kerr-Schild
components are used only to check initial data, not as evolved exact fields.
The spec and hybrid-p policies hold h fixed and default to activation at t=0.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersGaugeRelaxation as relax
from spectre.tools.Amr import ChalmersSpecCoarsening as coarsening

EXPERIMENT = "tilted-kerr-damped-harmonic-relaxation-v1"


def configuration(args):
    p_only = args.policy in ("spec", "hybrid-p")
    args.strategy = (
        "fixed" if args.policy == "fixed" else ("p" if p_only else "hp05")
    )
    if p_only:
        args.h_coarsening = False
    args.constraint_feedback = "off" if args.policy == "fixed" else "tighten"
    config = relax.configuration(args)
    direction = np.asarray(args.spin_direction, dtype=float)
    if not np.all(np.isfinite(direction)) or np.linalg.norm(direction) == 0:
        raise ValueError("Spin direction must be finite and nonzero")
    if not np.isfinite(args.spin) or not 0 <= args.spin < 1:
        raise ValueError("Require a subextremal nonnegative spin magnitude")
    spin = args.spin * direction / np.linalg.norm(direction)
    config["InitialData"]["GeneralizedHarmonic(KerrSchild)"][
        "Spin"
    ] = spin.tolist()
    # Keep sufficient angular resolution for the tilted horizon's spin
    # diagnostic. Horizon controls and moving maps remain disabled.
    horizon = config["ApparentHorizons"]["ApparentHorizon"]
    horizon["Criteria"][0]["Residual"]["MinResolutionL"] = args.horizon_l
    horizon["InitialGuess"]["InitialL"] = max(8, args.horizon_l)
    config["ApparentHorizons"]["LMax"] = max(10, args.horizon_l)
    creator = config["DomainCreator"]["Sphere"]
    if args.angular_points is not None or args.radial_points is not None:
        creator["InitialGridPoints"] = [
            [
                *(
                    [args.angular_points[shell]] * 2
                    if args.angular_points is not None
                    else args.grid_points[:2]
                ),
                (
                    args.radial_points[shell]
                    if args.radial_points is not None
                    else args.grid_points[2]
                ),
            ]
            for shell in range(4)
            for _ in range(6)
        ]
    if args.radial_levels is not None or args.angular_levels is not None:
        creator["InitialRefinement"] = [
            [
                *(
                    [args.angular_levels[shell]] * 2
                    if args.angular_levels
                    else [0, 0]
                ),
                (
                    args.radial_levels[shell]
                    if args.radial_levels
                    else args.radial_level
                ),
            ]
            for shell in range(4)
            for _ in range(6)
        ]
    # The policy bounds polynomial degree: degree 6 gives 7 points.
    config["Amr"]["Policies"]["Limits"]["NumPolynomialModes"][0] = 6
    if args.policy != "fixed":
        old = config["Amr"]["Criteria"][0]["ChalmersWithConstraintFeedback"]
        config["Amr"]["Criteria"] = [
            {
                "ChalmersWithSpecCoarsening": {
                    "RefinementCriterion": old,
                    "RelativeTruncationTarget": args.spectral_target,
                    "ConstraintTargets": args.constraint_targets,
                    "EnableCoarsening": (
                        args.policy in ("ordinary", "pileup", "hybrid-p")
                    ),
                    "ForceDropPileUpModes": (
                        args.policy in ("pileup", "hybrid-p")
                    ),
                    "MinimumModes": 7,
                    "CallsBetweenChanges": args.calls_between_changes,
                    "TruncationErrorAdjustmentLimit": args.target_shift_limit,
                    "StartTime": args.start_time,
                    "Verbose": True,
                }
            }
        ]
        if args.policy == "spec":
            settings = config["Amr"]["Criteria"][0][
                "ChalmersWithSpecCoarsening"
            ]
            del settings["RefinementCriterion"]
            del settings["EnableCoarsening"]
            settings["ForceDropPileUpModes"] = True
            settings["TimeBetweenChanges"] = args.time_between_changes
            if args.target_profile is not None:
                settings["RelativeTruncationTarget"] = yaml.safe_load(
                    args.target_profile.read_text()
                )
            if getattr(args, "radial_minimum_modes", None) is not None:
                settings["MinimumModes"] = [
                    [7, 7, minimum]
                    for minimum in args.radial_minimum_modes
                    for _ in range(6)
                ]
            config["Amr"]["Criteria"] = [{"SpecSpectral": settings}]
    if args.checkpoint_at is not None:
        config["PhaseChangeAndTriggers"].append(
            {
                "Trigger": {
                    "Times": {"Specified": {"Values": [args.checkpoint_at]}}
                },
                "PhaseChanges": ["VisitAndReturn(WriteCheckpoint)"],
            }
        )
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--policy",
        choices=["original", "ordinary", "pileup", "fixed", "spec", "hybrid-p"],
        required=True,
    )
    parser.add_argument("--spin", type=float, default=0.5)
    parser.add_argument(
        "--spin-direction", type=float, nargs=3, default=[1.0, 2.0, 3.0]
    )
    parser.add_argument("--grid-points", type=int, nargs=3, default=[12, 12, 8])
    parser.add_argument(
        "--angular-points",
        type=int,
        nargs=4,
        help="Angular extents in the four radial shells, innermost first",
    )
    parser.add_argument(
        "--radial-points",
        type=int,
        nargs=4,
        help="Radial extents in the four radial shells, innermost first",
    )
    parser.add_argument("--radial-level", type=int, default=0)
    parser.add_argument("--radial-levels", type=int, nargs=4)
    parser.add_argument("--angular-levels", type=int, nargs=4)
    parser.add_argument(
        "--target-profile",
        type=Path,
        help="Spec policy grid-frame Background/Gaussians YAML",
    )
    parser.add_argument(
        "--radial-minimum-modes",
        type=int,
        nargs=4,
        help="Spec policy radial floors by shell; angular floors stay at 7",
    )
    parser.add_argument("--inner-radius", type=float, default=1.2)
    parser.add_argument("--horizon-l", type=int, default=10)
    parser.add_argument("--cfl", type=float, default=0.5)
    parser.add_argument("--end-time", type=float, default=5.0)
    parser.add_argument("--checkpoint-at", type=float)
    parser.add_argument("--start-time", type=float, default=None)
    parser.add_argument(
        "--calls-between-changes",
        type=int,
        default=2,
        help="P-change cooldown; ten calls equal 1M at this driver's cadence",
    )
    parser.add_argument(
        "--time-between-changes",
        type=float,
        help="Elapsed-time cooldown for spec (default 1M; otherwise zero)",
    )
    parser.add_argument("--processes", type=int, default=10)
    parser.add_argument(
        "--target-shift-limit",
        type=float,
        default=None,
        help="Stored log10 target limit; default -1 for spec/hybrid-p, else 0",
    )
    parser.add_argument("--timeout", type=float, default=86400.0)
    parser.add_argument("--spectral-target", type=float, default=1.0e-4)
    parser.add_argument(
        "--constraint-targets", type=float, nargs=2, default=[0.002, 0.004]
    )
    parser.add_argument("--write-only", action="store_true")
    args = parser.parse_args()
    if args.time_between_changes is None:
        args.time_between_changes = 1.0 if args.policy == "spec" else 0.0
    if args.policy != "spec" and args.time_between_changes != 0.0:
        parser.error("A nonzero --time-between-changes requires --policy spec")
    if args.checkpoint_at is not None and not (
        0 < args.checkpoint_at < args.end_time
        and np.isclose(args.checkpoint_at * 10, round(args.checkpoint_at * 10))
    ):
        parser.error(
            "--checkpoint-at must lie on an interior 0.1M AMR boundary"
        )
    if args.target_profile is not None and args.policy != "spec":
        parser.error("--target-profile requires --policy spec")
    if args.radial_minimum_modes is not None:
        if args.policy != "spec" or min(args.radial_minimum_modes) < 7:
            parser.error("Radial floors require --policy spec and values >= 7")
        initial = args.radial_points or [args.grid_points[2]] * 4
        if any(n < m for n, m in zip(initial, args.radial_minimum_modes)):
            parser.error("Initial radial points must meet each radial floor")
    if args.start_time is None:
        args.start_time = 0.0 if args.policy in ("spec", "hybrid-p") else 10.0
    if args.target_shift_limit is None:
        args.target_shift_limit = (
            -1.0 if args.policy in ("spec", "hybrid-p") else 0.0
        )
    if (
        not np.isfinite(args.time_between_changes)
        or args.time_between_changes < 0
    ):
        parser.error("--time-between-changes must be finite and nonnegative")
    if args.calls_between_changes < 1:
        parser.error("--calls-between-changes must be positive")
    if not np.isfinite(args.target_shift_limit) or args.target_shift_limit > 0:
        parser.error("--target-shift-limit must be finite and nonpositive")
    args.points = max(args.grid_points)
    args.tolerance = 1.0e-4
    args.coarsening_factor = 0.05
    args.h_coarsening = True
    args.observation_interval = min(1.0, args.end_time)
    args.early_observation_until = min(0.5, args.end_time)
    args.experiment = EXPERIMENT
    args.criterion = "parity"
    args.filter = "hypercube"
    args.half_power = 256
    args.constraint_diagnostics = True
    args.mesh_observations = True
    args.constraint_target = 3.0e-5
    args.minimum_target_factor = 0.1
    if args.policy in ("spec", "hybrid-p"):
        args.estimator_variant = (
            "spec-grouped-tensor"
            if args.policy == "spec"
            else "chalmers-directional-energy-with-spec-tensor"
        )
    config = configuration(args)
    if args.target_profile is not None:
        args.target_profile = str(args.target_profile.resolve())
    if args.write_only:
        args.output.mkdir(parents=True, exist_ok=False)
        (args.output / "Input.yaml").write_text(
            "Executable: EvolveGhSingleBlackHole\n---\n"
            + yaml.safe_dump(config, sort_keys=False)
        )
        return
    sources = {
        str(p): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (
            Path(__file__).resolve(),
            Path(relax.__file__).resolve(),
            Path(cube.__file__).resolve(),
            Path(coarsening.__file__).resolve(),
        )
    }
    try:
        status = cube.run(
            args,
            config=config,
            analyzer=(
                coarsening.analyze_log
                if args.policy != "fixed"
                else lambda _: {"analysis": "deferred 3D comparison"}
            ),
            executable_name="EvolveGhSingleBlackHole",
        )
    finally:
        if (args.output / "Run.json").exists():
            (args.output / "Sources.json").write_text(
                json.dumps(sources, indent=2) + "\n"
            )
    raise SystemExit(status)


if __name__ == "__main__":
    main()
