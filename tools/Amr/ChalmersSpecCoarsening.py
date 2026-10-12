# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Controlled single-hole tests of SpEC-style p-coarsening.

Keep the existing refinement estimator and physics. Switch p-coarsening and
the corresponding p-addition permission rules at a declared time. Record native
pre-AMR tensor spectra and both projected constraint norms; analyze evolved
fields against the spherical reference separately with
ChalmersGaugeRelaxationAnalysis.
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import yaml

from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersGaugeRelaxation as relax
from spectre.tools.Amr.ChalmersEvolution import write_csv


def audit_saved(source, output, times):
    """Apply the exact new native diagnostics to saved, post-AMR fields."""
    import h5py

    from spectre.DataStructures import DataVector
    from spectre.DataStructures.Tensor import Scalar, tnsr
    from spectre.Domain import ElementId, ElementMap
    from spectre.Evolution.Systems.GeneralizedHarmonic import (
        spec_amr_damped_harmonic_diagnostics,
    )
    from spectre.Spectral import Basis, Mesh3D, Quadrature, logical_coordinates
    from spectre.tools.Amr import ChalmersSingleHole as gh
    from spectre.tools.Amr import ChalmersSphere as sphere

    config = list(yaml.safe_load_all((source / "Input.yaml").read_text()))[1]
    system = config["EvolutionSystem"]["GeneralizedHarmonic"]
    gauge = system["GaugeCondition"]["DampedHarmonic"]
    damping = system["DampingFunctionGamma2"]["GaussianPlusConstant"]
    domain = sphere.read_domain(source)
    rows = []
    with h5py.File(source / "Volume0.h5") as handle:
        for group in sphere.groups(handle):
            time = float(group.attrs["observation_value"])
            if time not in times:
                continue
            for name, extents, xyz, fields in gh.elements(group):
                mesh = Mesh3D(extents, Basis.Legendre, Quadrature.GaussLobatto)
                element_map = ElementMap(ElementId[3](name), domain)
                logical = logical_coordinates(mesh)
                tensors = [
                    tnsr.aa[DataVector, 3](mesh.number_of_grid_points(), 0.0)
                    for _ in range(2)
                ]
                phi = tnsr.iaa[DataVector, 3](mesh.number_of_grid_points(), 0.0)
                for k, (a, b) in enumerate(gh.PAIRS):
                    for offset, tensor in zip([0, 10], tensors):
                        tensor[tensor.get_storage_index(a, b)] = DataVector(
                            fields[offset + k].ravel(order="F")
                        )
                    for d in range(3):
                        phi[phi.get_storage_index(d, a, b)] = DataVector(
                            fields[20 + 10 * d + k].ravel(order="F")
                        )
                radius_squared = np.sum(
                    (xyz - np.array(damping["Center"])[:, None, None, None])
                    ** 2,
                    axis=0,
                )
                gamma2 = damping["Constant"] + damping["Amplitude"] * np.exp(
                    -radius_squared / damping["Width"] ** 2
                )
                diagnostic = spec_amr_damped_harmonic_diagnostics(
                    *tensors,
                    phi,
                    mesh,
                    element_map.inv_jacobian(logical),
                    element_map(logical),
                    time,
                    Scalar[DataVector]([DataVector(gamma2.ravel(order="F"))]),
                    gauge["SpatialDecayWidth"],
                    gauge["Amplitudes"],
                    gauge["Exponents"],
                )
                rows.append(
                    dict(
                        time=time,
                        element=name,
                        extents=extents.tolist(),
                        **diagnostic,
                    )
                )
    if not rows:
        raise ValueError("No selected snapshots")
    output.mkdir(parents=True, exist_ok=True)
    (output / "SavedDiagnostics.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )
    summary = dict(
        elements=len(rows),
        source=str(source.resolve()),
        times=sorted(set(r["time"] for r in rows)),
        constraint_quantiles={
            name: (
                np.quantile(
                    [c for r in rows for c in r["constraints"][i]],
                    [0.5, 0.9, 0.99, 1.0],
                ).tolist()
            )
            for i, name in enumerate(["c2", "c3"])
        },
        maximum_minimum_pileup=max(
            min(r["spectra"][i][d]["pile_up"] for i in range(2))
            for r in rows
            for d in range(3)
        ),
        both_tensors_piled_up=sum(
            min(r["spectra"][i][d]["pile_up"] for i in range(2)) > 3
            for r in rows
            for d in range(3)
        ),
    )
    (output / "SavedDiagnosticsSummary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def configuration(args):
    config = relax.configuration(args)
    old = config["Amr"]["Criteria"][0]["ChalmersWithConstraintFeedback"]
    config["Amr"]["Criteria"] = [
        {
            "ChalmersWithSpecCoarsening": {
                "RefinementCriterion": old,
                "RelativeTruncationTarget": args.spectral_target,
                "ConstraintTargets": args.constraint_targets,
                "EnableCoarsening": args.policy != "observe",
                "ForceDropPileUpModes": args.policy == "pileup",
                "MinimumModes": 7,
                "CallsBetweenChanges": 2,
                "TruncationErrorAdjustmentLimit": 0.0,
                "StartTime": args.start_time,
                "Verbose": True,
            }
        }
    ]
    # The policy bounds polynomial degree: degree 6 gives 7 points.
    config["Amr"]["Policies"]["Limits"]["NumPolynomialModes"][0] = 6
    return config


def analyze_log(directory):
    pattern = re.compile(
        r"SpecCoarsening time=(\S+) element=(\S+) direction=(\d+) "
        r"extent=(\d+) metric=\(([^)]+)\) derivative=\(([^)]+)\) "
        r"constraints=\(([^)]+)\) before=(\w+) after=(\w+)"
        r"(?: target_shift=(\S+) effective_target=(\S+))?"
    )
    records = []
    for match in pattern.finditer((directory / "Run.log").read_text()):
        (
            time,
            element,
            direction,
            extent,
            metric,
            deriv,
            constraints,
            before,
            after,
            target_shift,
            effective_target,
        ) = match.groups()
        record = dict(
            time=float(time),
            element=element,
            direction=int(direction),
            extent=int(extent),
            before=before,
            after=after,
        )
        for name, values in [("metric", metric), ("derivative", deriv)]:
            record.update(
                zip(
                    [
                        name + "_" + key
                        for key in ("error", "next", "slope", "pileup")
                    ],
                    map(float, values.split(",")),
                )
            )
        record["c2"], record["c3"] = map(float, constraints.split(","))
        record["minimum_pileup"] = min(
            record["metric_pileup"], record["derivative_pileup"]
        )
        if target_shift is not None:
            record["target_shift"] = float(target_shift)
            record["effective_target"] = float(effective_target)
        records.append(record)
    if not records:
        raise ValueError("No native pre-AMR diagnostics")
    write_csv(directory / "SpecDecisions.csv", records)
    changed = [r for r in records if r["before"] != r["after"]]
    summary = dict(
        directional_evaluations=len(records),
        decisions=dict(Counter(r["after"] for r in records)),
        changed_decisions=len(changed),
        changes=dict(Counter(r["before"] + "->" + r["after"] for r in changed)),
        max_minimum_pileup=max(r["minimum_pileup"] for r in records),
        both_tensors_piled_up=sum(r["minimum_pileup"] > 3 for r in records),
        projected_constraint_quantiles={
            key: (
                np.quantile(
                    [r[key] for r in records], [0.5, 0.9, 0.99, 1.0]
                ).tolist()
            )
            for key in ["c2", "c3"]
        },
        first_change=min((r["time"] for r in changed), default=None),
    )
    (directory / "SpecAnalysis.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--policy",
        choices=["observe", "tolerance", "pileup"],
        default="observe",
    )
    parser.add_argument("--start-time", type=float, default=10.0)
    parser.add_argument("--end-time", type=float, default=20.0)
    parser.add_argument("--spectral-target", type=float, default=1.0e-4)
    parser.add_argument(
        "--constraint-targets", type=float, nargs=2, default=[3.0e-5, 3.0e-5]
    )
    parser.add_argument(
        "--constraint-feedback",
        choices=["baseline", "tighten"],
        default="tighten",
    )
    parser.add_argument("--processes", type=int, default=5)
    parser.add_argument("--write-only", action="store_true")
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--audit-source", type=Path)
    parser.add_argument(
        "--times", type=float, nargs="+", default=[10.0, 15.0, 20.0]
    )
    args = parser.parse_args()
    if args.audit_source is not None:
        print(
            json.dumps(
                audit_saved(args.audit_source, args.output, args.times),
                indent=2,
            )
        )
        return
    if args.analyze_only:
        print(json.dumps(analyze_log(args.output), indent=2))
        return
    if args.executable is None:
        parser.error("--executable is required for an evolution")
    # Match the established 20M gauge-relaxation comparison.
    args.strategy = "hp05"
    args.points = 12
    args.grid_points = [12, 12, 8]
    args.radial_level = 0
    args.inner_radius = 1.2
    args.tolerance = 1.0e-4
    args.coarsening_factor = 0.05
    args.h_coarsening = True
    args.cfl = 0.5
    args.observation_interval = min(1.0, args.end_time)
    args.early_observation_until = min(0.5, args.end_time)
    args.timeout = 18000
    args.experiment = relax.EXPERIMENT
    args.criterion = "parity"
    args.filter = "hypercube"
    args.half_power = 256
    args.constraint_diagnostics = True
    args.mesh_observations = True
    args.constraint_target = 3.0e-5
    args.minimum_target_factor = 0.1
    config = configuration(args)
    if args.write_only:
        args.output.mkdir(parents=True, exist_ok=False)
        (args.output / "Input.yaml").write_text(
            "Executable: EvolveGhSingleBlackHole\n---\n"
            + yaml.safe_dump(config, sort_keys=False)
        )
        return
    raise SystemExit(
        cube.run(
            args,
            config=config,
            analyzer=analyze_log,
            executable_name="EvolveGhSingleBlackHole",
        )
    )


if __name__ == "__main__":
    main()
