# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Fixed spherical-symmetry reference for the identical gauge-relaxation problem.

This is a reference calculation, not a Chalmers AMR test. Angular derivatives
are supplied by the Cartoon spherical symmetry implementation. It must be
checked for radial/time convergence and against the full 3D evolutions.
"""

import argparse
import hashlib
import json
from pathlib import Path

import yaml

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr.ChalmersGaugeRelaxation import configuration


def config(args):
    result = configuration(args)
    sphere = result["DomainCreator"]["Sphere"]
    result["DomainCreator"] = {
        "CartoonSphere1D": {
            "InnerRadius": sphere["InnerRadius"],
            "OuterRadius": sphere["OuterRadius"],
            "InitialRadialRefinement": args.radial_levels or args.radial_level,
            "InitialNumberOfRadialGridPoints": args.points,
            "RadialPartitioning": sphere["RadialPartitioning"],
            "RadialDistributions": sphere["RadialDistribution"],
            "TimeDependence": "None",
            "InnerBoundaryCondition": "DemandOutgoingCharSpeeds",
            "OuterBoundaryCondition": sphere["OuterBoundaryCondition"],
        }
    }
    # No off-axis interpolation targets are used for the one-dimensional
    # reference. The three-dimensional runs retain horizon observations.
    events = result["EventsAndTriggersAtSlabs"][0]["Events"]
    events[:] = [e for e in events if isinstance(e, dict)]
    for bundle in result["EventsAndTriggersAtSlabs"]:
        bundle["Events"] = [
            e
            for e in bundle["Events"]
            if not isinstance(e, dict) or "ObserveAmrStats" not in e
        ]
    for key in [
        "Interpolator",
        "ApparentHorizons",
        "InterpolationTargets",
        "Cce",
        "ControlSystems",
        "EventsRunAtCleanup",
    ]:
        result.pop(key)
    result["Observers"].pop("SurfaceFileName")
    result["Amr"]["Policies"]["Limits"]["NumPolynomialModes"] = "Auto"
    result["Amr"]["Policies"]["Limits"]["RefinementLevel"] = "Auto"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--points", type=int, default=24)
    parser.add_argument("--radial-level", type=int, default=1)
    parser.add_argument("--radial-levels", type=int, nargs=4)
    parser.add_argument("--end-time", type=float, default=5)
    parser.add_argument("--observation-interval", type=float, default=1)
    parser.add_argument("--cfl", type=float, default=0.5)
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--write-only", action="store_true")
    parser.set_defaults(
        strategy="fixed",
        grid_points=None,
        inner_radius=1.2,
        tolerance=1e-5,
        coarsening_factor=0.05,
        h_coarsening=False,
        early_observation_until=0.0,
        timeout=1800,
        experiment="schwarzschild-damped-harmonic-cartoon-reference-v1",
        criterion="parity",
        filter="hypercube",
        half_power=256,
        constraint_diagnostics=True,
        mesh_observations=True,
    )
    args = parser.parse_args()
    if not 2 <= args.points <= 21:
        parser.error(
            "Legendre meshes in this build support 2 through 21 points"
        )
    args.early_observation_until = min(0.5, args.end_time)
    configuration = config(args)
    if args.write_only:
        args.output.write_text(
            "Executable: EvolveGhNoBlackHole3D\n---\n"
            + yaml.safe_dump(configuration, sort_keys=False)
        )
    else:
        source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        status = cube.run(
            args,
            config=configuration,
            analyzer=lambda _: {
                "analysis": "deferred spherical reconstruction"
            },
            executable_name="EvolveGhNoBlackHole3D",
        )
        (args.output / "Driver.sha256").write_text(source_hash + "\n")
        raise SystemExit(status)


if __name__ == "__main__":
    main()
