# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Filtered Chalmers AMR for spherical-harmonic scalar waves on Kerr-Schild.

M=1, gamma1=0, gamma2=1. The very steep filter zeros the highest mode;
ChalmersIgnoringHighestMode fits the remaining coefficients. Fixed and
adaptive controls use the same domain, initial data and filtering.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.ChalmersEvolution import write_csv


def configuration(args):
    m = getattr(args, "m", 0)
    if abs(m) > args.ell:
        raise ValueError("Spherical-harmonic mode requires |m| <= ell")
    base = argparse.Namespace(**vars(args))
    base.amr_interval = 0.4
    base.orientation = "oblique"
    base.aggregation = "energy"
    base.endpoints_only = False
    base.max_level = getattr(args, "max_level", 4)
    base.stepping = "cfl"
    base.fixed_p = args.points - 1
    config = sphere.configuration(base)
    config["InitialData"] = {
        "PureSphericalHarmonic": {
            "Radius": 10.0,
            "Width": 2.0,
            "Mode": [args.ell, m],
        }
    }
    config["BackgroundSpacetime"] = {
        "KerrSchild": {
            "Mass": 1.0,
            "Spin": [0.0, 0.0, args.spin],
            "Center": [0.0, 0.0, 0.0],
            "Velocity": [0.0, 0.0, 0.0],
        }
    }
    domain = config["DomainCreator"]["Sphere"]
    domain.update(
        InnerRadius=1.7,
        OuterRadius=40.0,
        Interior={"ExciseWithBoundaryCondition": "DemandOutgoingCharSpeeds"},
        InitialGridPoints=[args.points, args.points, args.points],
        InitialRefinement=[0, 0, args.radial_level],
        RadialPartitioning=[3.0, 6.0, 10.0, 16.0, 24.0],
        OuterBoundaryCondition="ConstraintPreservingSphericalRadiation",
    )
    amr = config["Amr"]
    if args.strategy != "fixed":
        criterion = amr["Criteria"][0].pop("Chalmers")
        amr["Criteria"][0][
            (
                "Chalmers"
                if args.include_highest
                else "ChalmersIgnoringHighestMode"
            )
        ] = criterion
        if getattr(args, "separate_parity", False):
            if args.include_highest:
                raise ValueError(
                    "Parity evolution experiment requires the steep-filter"
                    " criterion"
                )
            amr["Criteria"][0] = {
                "ChalmersParityIgnoringHighestMode": criterion
            }
    # Limits currently interpret these bounds as polynomial degrees. The
    # criterion itself keeps at least six usable (unfiltered) modes.
    amr["Policies"]["Limits"]["NumPolynomialModes"] = [2, 19]
    config["SpatialDiscretization"]["DiscontinuousGalerkin"]["Filtering"] = [
        {
            "Hypercube": {
                "HalfPower": args.half_power,
                "Enable": True,
                "BlocksToFilter": "All",
                "VolumeFilterOnSubstep": True,
                "BoundaryCorrectionFilterOnSubstep": False,
                "VolumeFilterEveryNSteps": "None",
                "BoundaryCorrectionFilterEveryNSteps": "None",
            }
        }
    ]
    # Replace flat-space analytic error events by norms; errors are measured
    # independently against fixed-grid controls.
    for entry in config["EventsAndTriggersAtSlabs"]:
        for event in entry["Events"]:
            if not isinstance(event, dict):
                continue
            if "ObserveNorms" in event:
                event["ObserveNorms"] = {
                    "SubfileName": "Norms",
                    "TensorsToObserve": [
                        {
                            "Name": name,
                            "NormType": "L2Norm",
                            "Components": "Sum",
                        }
                        for name in [
                            "Psi",
                            "Pi",
                            "Phi",
                            "PointwiseL2Norm(OneIndexConstraint)",
                            "PointwiseL2Norm(TwoIndexConstraint)",
                        ]
                    ],
                }
            if "ObserveFields" in event:
                event["ObserveFields"]["VariablesToObserve"] = [
                    "Psi",
                    "Pi",
                    "Phi",
                    "OneIndexConstraint",
                    "TwoIndexConstraint",
                ]
    config["InterpolationTargets"] = {
        "SphericalSurface": {
            "LMax": 8,
            "Center": [0.0, 0.0, 0.0],
            "Radius": 12.0,
            "AngularOrdering": "Strahlkorper",
        }
    }
    config["Interpolator"] = {"Verbosity": "Silent"}
    return config


def analyze_run(directory):
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata["returncode"]:
        return cube.analyze_run(directory)
    result = {"run": directory.name, **metadata}
    with h5py.File(directory / "Reductions.h5") as f:
        stats = f["Amr.dat"][:]
        times = f["TimeSteps.dat"][:]
        norms = f["Norms.dat"][:]
        legend = [
            x.decode() if isinstance(x, bytes) else str(x)
            for x in f["Norms.dat"].attrs["Legend"]
        ]
        result.update(
            final_elements=int(stats[-1, 1]),
            final_points=int(stats[-1, 2]),
            minimum_elements=int(np.min(stats[:, 1])),
            minimum_points=int(np.min(stats[:, 2])),
            maximum_elements=int(np.max(stats[:, 1])),
            maximum_points=int(np.max(stats[:, 2])),
            initial_elements=int(stats[0, 1]),
            initial_points=int(stats[0, 2]),
            final_time=float(times[-1, 0]),
            steps=len(times) - 1,
            evolution_seconds=float(times[-1, -1] - times[0, -1]),
            norm_legend=legend,
            final_norms=norms[-1].tolist(),
        )
        write_csv(
            directory / "Norms.csv", [dict(zip(legend, row)) for row in norms]
        )
        write_csv(
            directory / "MeshHistory.csv",
            [
                {"time": row[0], "elements": int(row[1]), "points": int(row[2])}
                for row in stats
            ],
        )
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--ell", type=int, choices=[0, 2], default=0)
    parser.add_argument("--m", type=int, default=0)
    parser.add_argument("--spin", type=float, default=0.0)
    parser.add_argument(
        "--strategy", choices=["fixed", "p", "hp05"], default="hp05"
    )
    parser.add_argument("--points", type=int, default=12)
    parser.add_argument("--max-level", type=int, default=4)
    parser.add_argument("--radial-level", type=int, default=0)
    parser.add_argument("--half-power", type=int, default=256)
    parser.add_argument("--include-highest", action="store_true")
    parser.add_argument("--separate-parity", action="store_true")
    parser.add_argument("--coarsening-factor", type=float, default=0.1)
    parser.add_argument("--h-coarsening", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1.0e-4)
    parser.add_argument("--end-time", type=float, default=40.0)
    parser.add_argument("--observation-interval", type=float, default=2.0)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--cfl", type=float, default=0.25)
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    raise SystemExit(
        cube.run(
            args,
            config=configuration(args),
            analyzer=analyze_run,
            executable_name="EvolveCurvedScalarWaveKerrSchild3D",
        )
    )


if __name__ == "__main__":
    main()
