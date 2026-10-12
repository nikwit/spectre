# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Small scalar-wave controls for time stepping and p-coarsening constraints.

The fixed and adaptive controls share their initial mesh. Sibling joining is
disabled in this experiment to isolate the p-projection contribution.
"""

import argparse
from pathlib import Path

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersEvolution as pulse
from spectre.tools.Amr import ChalmersSphere as sphere


def configuration(args):
    base = argparse.Namespace(**vars(args))
    base.strategy = "fixed" if args.fixed_grid else "hp05"
    base.coarsening_factor = 0.0 if args.fixed_grid else 0.1
    base.h_coarsening = False
    if args.model == "pulse":
        config = pulse.input_file(
            base.strategy,
            args.tolerance,
            args.dt,
            args.end_time,
            True,
            coarsening_factor=base.coarsening_factor,
        )
        domain = config["DomainCreator"]["RotatedIntervals"]
        domain["InitialGridPoints"] = [[16, 16]]
        domain["InitialRefinement"] = [3]
        config["Amr"]["Policies"]["Limits"]["RefinementLevel"] = [3, 7]
    else:
        base.fixed_p = 13
        base.max_level = 3
        base.observation_interval = 0.2
        base.stepping = "cfl"
        config = sphere.configuration(base)
        config["DomainCreator"]["Sphere"]["InitialGridPoints"] = [14, 14, 18]
    if args.constraint_target is not None:
        if args.fixed_grid:
            raise ValueError("A fixed grid has no coarsening to guard")
        chalmers = config["Amr"]["Criteria"][0]["Chalmers"]
        config["Amr"]["Criteria"] = [
            {
                "ChalmersWithConstraintGuard": {
                    "Chalmers": chalmers,
                    "ConstraintChangeTarget": args.constraint_target,
                }
            }
        ]
    return config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=["pulse", "sphere"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--fixed-grid", action="store_true")
    parser.add_argument("--constraint-target", type=float)
    parser.add_argument("--tolerance", type=float, required=True)
    parser.add_argument("--end-time", type=float, required=True)
    parser.add_argument("--dt", type=float, default=0.0005)
    parser.add_argument("--cfl", type=float, default=0.25)
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    args.strategy = "fixed" if args.fixed_grid else "hp05"
    raise SystemExit(
        cube.run(
            args,
            config=configuration(args),
            analyzer=(
                pulse.analyze_run
                if args.model == "pulse"
                else sphere.analyze_run
            ),
            executable_name=(
                "EvolveScalarWave1D"
                if args.model == "pulse"
                else "EvolveScalarWave3D"
            ),
        )
    )


if __name__ == "__main__":
    main()
