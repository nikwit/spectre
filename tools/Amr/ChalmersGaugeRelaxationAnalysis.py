# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Reference comparisons and independent constraints for gauge relaxation.

All integrals use coordinate volume and full symmetric-tensor multiplicities.
The Euclidean component norms are diagnostics, not invariant spacetime norms.
No stationary analytic GH solution is used as an evolved error reference.
"""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import yaml

from spectre.Domain import ElementId, ElementMap
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSingleHole as gh
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv
from spectre.tools.Amr.ChalmersMeshHistory import analyze as mesh_changes

ORDERS = (28, 36)


def initial_metric(xyz, prescription):
    """Initial Kerr metric from the native analytic solution's 3+1 fields."""
    from spectre.DataStructures import DataVector
    from spectre.DataStructures.Tensor import tnsr
    from spectre.PointwiseFunctions.AnalyticSolutions.GeneralRelativity import (
        KerrSchild,
    )

    if np.any(np.asarray(prescription["Velocity"]) != 0):
        raise ValueError("Initial-data check requires an unboosted hole")
    solution = KerrSchild(
        prescription["Mass"], prescription["Spin"], prescription["Center"]
    )
    data = solution.variables(
        tnsr.I[DataVector, 3](np.asarray(xyz).reshape(3, -1)),
        ["Lapse", "Shift", "SpatialMetric"],
    )
    spatial_tensor = data["SpatialMetric"]
    spatial = np.array(
        [
            [np.asarray(spatial_tensor.get(i, j)) for j in range(3)]
            for i in range(3)
        ]
    )
    shift = np.asarray(data["Shift"])
    lower_shift = np.einsum("ij...,j...->i...", spatial, shift)
    metric = np.empty((4, 4, xyz.reshape(3, -1).shape[1]))
    metric[1:, 1:] = spatial
    metric[0, 1:] = metric[1:, 0] = lower_shift
    metric[0, 0] = -np.asarray(data["Lapse"])[0] ** 2 + np.sum(
        shift * lower_shift, axis=0
    )
    return np.array([metric[a, b] for a, b in gh.PAIRS]).reshape(
        10, *xyz.shape[1:]
    )


def reductions(directory):
    result = {}
    with h5py.File(directory / "Reductions.h5") as handle:
        for name, data in handle.items():
            if not name.endswith(".dat"):
                continue
            if not np.all(np.isfinite(data[:])):
                raise ValueError(f"Nonfinite reductions: {name}")
            legend = [
                s.decode() if isinstance(s, bytes) else str(s)
                for s in data.attrs["Legend"]
            ]
            result[name[:-4]] = dict(zip(legend, data[:].T))
            write_csv(
                directory / (name[:-4] + ".csv"),
                [dict(zip(legend, row)) for row in data[:]],
            )
    return result


def inner_characteristic_speed(fields, xyz):
    """Smallest propagating speed, normal pointing out of the excised domain.

    For static spherical excision, n_i=-x_i/(r sqrt(gamma^jk x_j x_k/r^2)).
    Positive -beta.n-alpha means both lightlike characteristics leave the
    computational domain. The zero-speed metric characteristic is excluded.
    """
    shape = fields.shape[1:]
    metric = np.empty((4, 4, *shape))
    for value, (a, b) in zip(fields[:10], gh.PAIRS):
        metric[a, b] = metric[b, a] = value
    inverse_spatial = np.linalg.inv(
        np.moveaxis(metric[1:, 1:], (0, 1), (-2, -1))
    )
    inverse_spatial = np.moveaxis(inverse_spatial, (-2, -1), (0, 1))
    shift = np.einsum("ij...,j...->i...", inverse_spatial, metric[0, 1:])
    lapse = np.sqrt(-metric[0, 0] + np.sum(shift * metric[0, 1:], axis=0))
    radial = xyz / np.linalg.norm(xyz, axis=0)
    normal_factor = np.sqrt(
        np.einsum("i...,ij...,j...->...", radial, inverse_spatial, radial)
    )
    return np.sum(shift * radial, axis=0) / normal_factor - lapse


def analyze(directory, allow_incomplete=False):
    metadata = json.loads((directory / "Run.json").read_text())
    completed = metadata["returncode"] == 0
    if not completed and not allow_incomplete:
        raise ValueError("Analyze only completed evolutions")
    if metadata["experiment"] not in (
        "schwarzschild-damped-harmonic-relaxation-v1",
        "tilted-kerr-damped-harmonic-relaxation-v1",
    ):
        raise ValueError("Incorrect experiment")
    config = list(yaml.safe_load_all((directory / "Input.yaml").read_text()))[1]
    prescription = config["InitialData"]["GeneralizedHarmonic(KerrSchild)"]
    expected_spin = float(np.linalg.norm(prescription["Spin"]))
    creator = config["DomainCreator"]["Sphere"]
    if creator["TimeDependentMaps"] != "None":
        raise ValueError("These diagnostics require static coordinate maps")
    domain = sphere.read_domain(directory)
    expected_volume = (
        4
        * np.pi
        / 3
        * (creator["OuterRadius"] ** 3 - creator["InnerRadius"] ** 3)
    )
    rows, element_rows = [], []
    initial_metric_error = 0.0
    with h5py.File(directory / "Volume0.h5") as handle:
        for group in sphere.groups(handle):
            time = float(group.attrs["observation_value"])
            sums, volumes = np.zeros(2), np.zeros(2)
            inner_speed = np.inf
            coord_error = 0.0
            for name, ext, xyz, fields in gh.elements(group):
                if not np.all(np.isfinite(fields)) or not np.all(
                    np.isfinite(xyz)
                ):
                    raise ValueError(f"Nonfinite volume data at {time}: {name}")
                if time == 0:
                    initial_metric_error = max(
                        initial_metric_error,
                        float(
                            np.max(
                                abs(
                                    fields[:10]
                                    - initial_metric(xyz, prescription)
                                )
                            )
                        ),
                    )
                element_map = ElementMap(ElementId[3](name), domain)
                mapped, _, _ = sphere.geometry(
                    element_map, [legendre_grid(n)[0] for n in ext]
                )
                coord_error = max(coord_error, float(np.max(abs(mapped - xyz))))
                modes = [to_modal(v) for v in fields]
                element_constraints = []
                for i, order in enumerate(ORDERS):
                    xi, _, weights = cube.quadrature(order)
                    _, determinant, inverse = sphere.geometry(
                        element_map, [xi] * 3
                    )
                    measure = weights * determinant
                    numerical_phi = (
                        np.array([cube.evaluate(c, order) for c in modes[20:]])
                        .reshape(3, 10, order, order, order)
                        .swapaxes(0, 1)
                    )
                    constraint = (
                        np.array(
                            [
                                sphere.physical_gradient(c, inverse, order)
                                for c in modes[:10]
                            ]
                        )
                        - numerical_phi
                    )
                    contribution = np.sum(
                        measure
                        * constraint**2
                        * gh.MULTIPLICITY[:10, None, None, None, None]
                    )
                    sums[i] += contribution
                    element_constraints.append(float(np.sqrt(contribution)))
                    volumes[i] += np.sum(measure)
                if abs(element_constraints[1] - element_constraints[0]) > max(
                    1e-11, 1e-5 * element_constraints[1]
                ):
                    raise ValueError(
                        f"Unconverged element constraint at {time}: {name}"
                    )
                element_rows.append(
                    dict(
                        time=time,
                        element=name,
                        block=ElementId[3](name).block_id,
                        three_index_constraint=element_constraints[1],
                        quadrature_change=abs(
                            element_constraints[1] - element_constraints[0]
                        ),
                    )
                )
                inner = np.isclose(
                    np.linalg.norm(xyz, axis=0),
                    creator["InnerRadius"],
                    atol=1e-12,
                    rtol=0,
                )
                if np.any(inner):
                    inner_speed = min(
                        inner_speed,
                        float(
                            np.min(
                                inner_characteristic_speed(
                                    fields[:, inner], xyz[:, inner]
                                )
                            )
                        ),
                    )
            values = np.sqrt(sums)
            if not np.all(np.isfinite(values)):
                raise ValueError(f"Nonfinite integrated constraints at {time}")
            delta = float(abs(values[1] - values[0]))
            if delta > max(1e-11, 1e-5 * values[1]):
                raise ValueError(
                    f"Unconverged constraint quadrature at {time}: {delta}"
                )
            if (
                np.max(abs(volumes / expected_volume - 1)) > 1e-11
                or coord_error > 1e-11
            ):
                raise ValueError("Geometry verification failed")
            if not np.isfinite(inner_speed) or inner_speed <= 0:
                raise ValueError("Excision is not outflow")
            rows.append(
                dict(
                    time=time,
                    three_index_constraint=float(values[1]),
                    quadrature_change=delta,
                    coordinate_error=coord_error,
                    min_excision_light_speed=inner_speed,
                )
            )
            print(directory.name, rows[-1], flush=True)
    data = reductions(directory)
    with h5py.File(directory / "Surfaces.h5") as handle:
        radii = {}
        for group in handle["ApparentHorizon.vol"].values():
            if not isinstance(group, h5py.Group):
                continue
            radius = np.linalg.norm(
                [group[f"InertialCoordinates_{d}"][:] for d in "xyz"], axis=0
            )
            radii[float(group.attrs["observation_value"])] = (
                float(np.min(radius)),
                float(np.max(radius)),
            )
        for row in rows:
            (
                row["min_horizon_coordinate_radius"],
                row["max_horizon_coordinate_radius"],
            ) = radii[row["time"]]
    times = np.array([r["time"] for r in rows])
    expected_times = config["EventsAndTriggersAtSlabs"][0]["Trigger"]["Times"][
        "Specified"
    ]["Values"]
    if not completed:
        expected_times = [t for t in expected_times if t <= times[-1]]
    if not np.array_equal(times, expected_times) or (
        completed and times[-1] != metadata["end_time"]
    ):
        raise ValueError("Missing volume observations")
    for name in ("Norms", "ApparentHorizon"):
        if not np.array_equal(data[name]["Time"], times):
            raise ValueError(f"Missing reductions: {name}")
    if (
        initial_metric_error > 1e-12
        or max(
            data["Norms"][f"Max(PointwiseL2Norm({constraint}))"][0]
            for constraint in ("GaugeConstraint", "ThreeIndexConstraint")
        )
        > 1e-12
    ):
        raise ValueError(
            "Initial metric or discrete GH constraints are inconsistent"
        )
    amr, steps, horizon = (
        data[n] for n in ("Amr", "TimeSteps", "ApparentHorizon")
    )
    amr_times, nelem, points = (
        amr[n] for n in ("Time", "NumElements", "TotalNumPoints")
    )
    with h5py.File(directory / "Volume0.h5") as handle:
        meshes = {
            float(g.attrs["observation_value"]): (
                np.asarray(g["total_extents"]).reshape(-1, 3)
            )
            for g in handle["Meshes.vol"].values()
            if isinstance(g, h5py.Group)
        }
        for time, count, size in zip(amr_times, nelem, points):
            extents = meshes[time]
            if (
                len(extents) != count
                or np.sum(np.prod(extents, axis=1)) != size
            ):
                raise ValueError("Volume mesh and AMR statistics disagree")
    if not np.array_equal(
        amr_times,
        np.round(np.arange(round(amr_times[-1] / 0.1) + 1) * 0.1, 12),
    ):
        raise ValueError("Missing AMR observations")
    if (
        steps["Time"][0] != 0
        or (completed and steps["Time"][-1] != metadata["end_time"])
        or (completed and amr_times[-1] != metadata["end_time"])
        or not times[-1] <= amr_times[-1] <= steps["Time"][-1]
        or steps["Time"][-1] > metadata["end_time"]
        or np.any(np.diff(steps["Time"]) <= 0)
    ):
        raise ValueError("Incomplete time-step history")
    spin_vector = np.array(
        [horizon[f"DimensionfulSpinVector_{d}"] for d in "xyz"]
    )
    expected_vector = prescription["Mass"] ** 2 * np.asarray(
        prescription["Spin"]
    )
    spin_angle = None
    if expected_spin > 0:
        cosine = np.einsum("i,it->t", expected_vector, spin_vector) / (
            np.linalg.norm(expected_vector)
            * np.linalg.norm(spin_vector, axis=0)
        )
        spin_angle = float(
            np.max(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
        )
    result = dict(
        run=directory.name,
        **metadata,
        completed=completed,
        last_volume_time=float(times[-1]),
        last_mesh_time=float(amr_times[-1]),
        last_time_step_time=float(steps["Time"][-1]),
        quadrature_orders=list(ORDERS),
        initial_metric_nodal_error=initial_metric_error,
        max_three_index_constraint=max(
            r["three_index_constraint"] for r in rows
        ),
        final_three_index_constraint=rows[-1]["three_index_constraint"],
        max_quadrature_change=max(r["quadrature_change"] for r in rows),
        min_excision_light_speed=min(
            r["min_excision_light_speed"] for r in rows
        ),
        max_horizon_mass_error=float(
            np.max(abs(horizon["ChristodoulouMass"] - prescription["Mass"]))
        ),
        max_dimensionless_spin_error=float(
            np.max(abs(horizon["DimensionlessSpinMagnitude"] - expected_spin))
        ),
        max_horizon_spin_vector_error=float(
            np.max(
                np.linalg.norm(spin_vector - expected_vector[:, None], axis=0)
            )
        ),
        max_spin_direction_error_degrees=spin_angle,
        max_dimensionless_spin=float(
            np.max(horizon["DimensionlessSpinMagnitude"])
        ),
        initial_points=int(points[0]),
        final_points=int(points[-1]),
        min_points=int(min(points)),
        max_points=int(max(points)),
        final_elements=int(nelem[-1]),
        max_elements=int(max(nelem)),
        point_increase_observations=int(np.sum(np.diff(points) > 0)),
        point_decrease_observations=int(np.sum(np.diff(points) < 0)),
        mean_points=float(
            np.dot(points[:-1], np.diff(amr_times)) / amr_times[-1]
            if amr_times[-1] > 0
            else points[0]
        ),
        global_steps=len(steps["Time"]) - 1,
        grid_point_steps=float(
            np.sum(steps["Number of degrees of freedom"][:-1]) / len(gh.FIELDS)
        ),
        native_constraint_maxima={
            k: float(max(v))
            for k, v in data["Norms"].items()
            if "Constraint" in k
        },
        mesh_changes=mesh_changes(directory, 3, 0.1),
    )
    prefix = "" if completed else "Partial"
    write_csv(directory / f"{prefix}History.csv", rows)
    write_csv(directory / f"{prefix}ConstraintByElement.csv", element_rows)
    (directory / f"{prefix}Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def compare(first, second, output, selected_times=None, allow_incomplete=False):
    """Compare at common times on the union of the two element partitions."""
    if allow_incomplete and selected_times is None:
        raise ValueError("Select comparison times explicitly for partial runs")
    run_status = []
    for directory in (first, second):
        metadata = json.loads((directory / "Run.json").read_text())
        run_status.append(
            dict(
                directory=str(directory),
                returncode=metadata["returncode"],
                requested_end_time=metadata["end_time"],
            )
        )
        if metadata["returncode"] and not allow_incomplete:
            raise ValueError("Compare only completed evolutions")
    domain = sphere.read_domain(first)
    if domain != sphere.read_domain(second):
        raise ValueError("Domains differ")
    rows = []
    with (
        h5py.File(first / "Volume0.h5") as a,
        h5py.File(second / "Volume0.h5") as b,
    ):
        groups_a = {
            float(g.attrs["observation_value"]): g for g in sphere.groups(a)
        }
        groups_b = {
            float(g.attrs["observation_value"]): g for g in sphere.groups(b)
        }
        times = (
            sorted(groups_a.keys() & groups_b.keys())
            if selected_times is None
            else selected_times
        )
        if not times:
            raise ValueError("No common observation times")
        for time in times:
            modes_a, modes_b = [
                {
                    name: [to_modal(v) for v in fields]
                    for name, _, _, fields in gh.elements(g)
                }
                for g in (groups_a[time], groups_b[time])
            ]
            values = np.array(
                [
                    sphere.curved_difference(modes_a, modes_b, domain, order)
                    for order in ORDERS
                ]
            )
            norms = np.sqrt(np.sum(values**2 * gh.MULTIPLICITY, axis=1))
            if not np.all(np.isfinite(norms)):
                raise ValueError(f"Nonfinite reference difference at {time}")
            delta = float(abs(norms[1] - norms[0]))
            if delta > max(1e-11, 1e-5 * norms[1]):
                raise ValueError(
                    f"Unconverged reference comparison at {time}: {delta}"
                )
            rows.append(
                dict(
                    time=time,
                    state_l2_difference=float(norms[1]),
                    metric_l2_difference=float(
                        np.linalg.norm(
                            values[1, :10] * np.sqrt(gh.MULTIPLICITY[:10])
                        )
                    ),
                    pi_l2_difference=float(
                        np.linalg.norm(
                            values[1, 10:20] * np.sqrt(gh.MULTIPLICITY[10:20])
                        )
                    ),
                    phi_l2_difference=float(
                        np.linalg.norm(
                            values[1, 20:] * np.sqrt(gh.MULTIPLICITY[20:])
                        )
                    ),
                    quadrature_change=delta,
                )
            )
            print(first.name, second.name, rows[-1], flush=True)
    result = dict(
        first=str(first),
        second=str(second),
        run_status=run_status,
        quadrature_orders=list(ORDERS),
        observations=rows,
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("directory", type=Path)
    analyzer.add_argument("--allow-incomplete", action="store_true")
    differ = commands.add_parser("compare")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.add_argument("--times", type=float, nargs="+")
    differ.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    if args.command == "analyze":
        print(
            json.dumps(analyze(args.directory, args.allow_incomplete)),
            flush=True,
        )
    else:
        compare(
            args.first,
            args.second,
            args.output,
            args.times,
            args.allow_incomplete,
        )


if __name__ == "__main__":
    main()
