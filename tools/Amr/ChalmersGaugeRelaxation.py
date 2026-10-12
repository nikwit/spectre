# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Single-hole gauge relaxation with Chalmers AMR and BBH outer conditions.

The initial spatial geometry is mass-one Schwarzschild in Kerr-Schild
coordinates. GH initializes Pi/Phi consistently with damped harmonic gauge.
The coordinates then evolve, so the stationary Kerr-Schild fields are NOT an
analytic error reference. Compare evolved fields against converged fixed runs.
"""

import argparse
import hashlib
import json
from pathlib import Path

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSingleHole as gh

EXPERIMENT = "schwarzschild-damped-harmonic-relaxation-v1"


def configuration(args):
    config = gh.configuration(args)
    creator = config["DomainCreator"]["Sphere"]
    creator.update(
        InnerRadius=args.inner_radius,
        OuterRadius=20.0,
        RadialPartitioning=[2.0, 4.0, 8.0],
        RadialDistribution="Logarithmic",
        OuterBoundaryCondition={
            "ConstraintPreservingBjorhus": {
                "Type": "ConstraintPreservingPhysical",
                "IncomingWaveProfile": "None",
            }
        },
    )
    config["EvolutionSystem"]["GeneralizedHarmonic"]["GaugeCondition"] = {
        "DampedHarmonic": {
            "SpatialDecayWidth": 17.0152695482514,
            "Amplitudes": [1.0, 0.0, 1.0],
            "Exponents": [2, 2, 2],
        }
    }
    for name in ["ExcisionBoundary", "ControlSystemCharSpeedExcision"]:
        config["InterpolationTargets"][name]["Radius"] = args.inner_radius
    events = config["EventsAndTriggersAtSlabs"][0]["Events"]
    # Error(SpacetimeMetric/Pi/Phi) assumes the stationary analytic solution.
    # It is inappropriate once the coordinates change.
    events[:] = [
        event
        for event in events
        if not (
            isinstance(event, dict)
            and event.get("ObserveNorms", {}).get("SubfileName") == "Errors"
        )
    ]
    events.append("ExcisionBoundary")
    for event in events:
        if isinstance(event, dict) and "ObserveFields" in event:
            event["ObserveFields"]["VariablesToObserve"].extend(
                ["Lapse", "GaugeH"]
            )
    feedback = getattr(args, "constraint_feedback", "off")
    if feedback != "off":
        if args.strategy == "fixed":
            raise ValueError("Constraint feedback requires an adaptive run")
        criterion = config["Amr"]["Criteria"][0][
            "ChalmersParityIgnoringHighestMode"
        ]
        config["Amr"]["Criteria"] = [
            {
                "ChalmersWithConstraintFeedback": {
                    "Chalmers": criterion,
                    "ConstraintTarget": args.constraint_target,
                    "MinimumTargetFactor": (
                        args.minimum_target_factor
                        if feedback == "tighten"
                        else 1.0
                    ),
                    "EnableCoarseningVeto": feedback != "baseline",
                    "Verbose": True,
                }
            }
        ]
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--strategy", choices=["fixed", "p", "hp05"], required=True
    )
    parser.add_argument("--points", type=int, default=12)
    parser.add_argument("--grid-points", type=int, nargs=3)
    parser.add_argument("--radial-level", type=int, default=0)
    parser.add_argument("--inner-radius", type=float, default=1.2)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    parser.add_argument("--coarsening-factor", type=float, default=0.05)
    parser.add_argument("--h-coarsening", action="store_true")
    parser.add_argument(
        "--constraint-feedback",
        choices=["off", "baseline", "veto", "tighten"],
        default="off",
    )
    parser.add_argument("--constraint-target", type=float)
    parser.add_argument("--minimum-target-factor", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.5)
    parser.add_argument("--end-time", type=float, default=5.0)
    parser.add_argument("--observation-interval", type=float, default=1.0)
    parser.add_argument("--early-observation-until", type=float, default=0.5)
    parser.add_argument("--processes", type=int, default=10)
    parser.add_argument("--timeout", type=float, default=14400)
    parser.set_defaults(
        experiment=EXPERIMENT,
        criterion="parity",
        filter="hypercube",
        half_power=256,
        constraint_diagnostics=True,
        mesh_observations=True,
    )
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output directory already exists; preserve prior runs")
    if not 0 < args.inner_radius < 2:
        parser.error("Inner radius must lie inside the initial horizon")
    if args.constraint_feedback != "off" and (
        args.constraint_target is None or not args.constraint_target > 0
    ):
        parser.error("Feedback requires a positive --constraint-target")
    config = configuration(args)
    # Record the generator before launching; inputs and Run.json are written
    # by cube.run. In particular, preserve failures as well as successful runs.
    sources = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [Path(__file__).resolve(), Path(gh.__file__).resolve()]
    }
    try:
        status = cube.run(
            args,
            config=config,
            analyzer=lambda _: {
                "analysis": "deferred; use numerical reference"
            },
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
