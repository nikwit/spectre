# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Native Chalmers experiments with a smooth pulse on curved cubed spheres.

The regular spherical wave has profile F(s)=exp(-((s-1.4)/0.4)^2).
Three shells cover 1 <= r <= 4 with six equiangular wedges each. Diagnostics
use the serialized, static element maps and independent Gauss quadrature.
"""

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import legder, legvander

import spectre.IO.H5 as spectre_h5
from spectre.DataStructures import DataVector
from spectre.DataStructures.Tensor import Frame, tnsr
from spectre.Domain import ElementId, ElementMap, deserialize_domain
from spectre.NumericalAlgorithms.LinearOperators import chalmers_energy_estimate
from spectre.Spectral import Basis, Mesh, Quadrature
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv

INNER_RADIUS = 1.0
OUTER_RADIUS = 4.0
WIDTH = 0.4
CENTER = 1.4
ORDERS = (24, 32)


def exact(xyz, t):
    """Return (Psi, Pi, Phi_i), with Psi=(F(r-t)-F(-r-t))/r, away from r=0."""
    r = np.sqrt(np.sum(np.asarray(xyz) ** 2, axis=0))

    def profile(s):
        u = s - CENTER
        f = np.exp(-((u / WIDTH) ** 2))
        return f, -2 * u / WIDTH**2 * f

    outgoing, dout = profile(r - t)
    incoming, din = profile(-r - t)
    psi = (outgoing - incoming) / r
    pi = (dout - din) / r
    phi = (dout + din - psi) / r**2 * np.asarray(xyz)
    return np.array([psi, pi, *phi])


def configuration(args):
    amr_interval = getattr(args, "amr_interval", cube.AMR_INTERVAL)
    for interval in [args.end_time, args.observation_interval]:
        if interval <= 0 or not np.isclose(
            interval / amr_interval,
            round(interval / amr_interval),
            rtol=0,
            atol=1e-12,
        ):
            raise ValueError(
                "End time and observations must align with the AMR schedule"
            )
    # Reuse the same physical AMR/observation schedule and time integrator.
    base = argparse.Namespace(**vars(args))
    base.orientation = "oblique"
    base.aggregation = "energy"
    base.endpoints_only = False
    config = cube.configuration(base)
    data = {
        "RegularSphericalWave": {
            "Profile": {
                "Gaussian": {
                    "Amplitude": 1.0,
                    "Width": WIDTH,
                    "Center": CENTER,
                }
            }
        }
    }
    config["InitialData"] = data
    boundary = {"DirichletAnalytic": {"AnalyticPrescription": data}}
    config["DomainCreator"] = {
        "Sphere": {
            "InnerRadius": INNER_RADIUS,
            "OuterRadius": OUTER_RADIUS,
            "Interior": {"ExciseWithBoundaryCondition": boundary},
            "InitialRefinement": [0, 0, 0],
            "InitialGridPoints": (
                [args.fixed_p + 1] * 3
                if args.strategy == "fixed"
                else {
                    "Shell0": [14, 14, 18],
                    "Shell1": [14, 14, 16],
                    "Shell2": [14, 14, 8],
                }
            ),
            "UseEquiangularMap": True,
            "EquatorialCompression": "None",
            "RadialPartitioning": [2.0, 3.0],
            "RadialDistribution": "Linear",
            "WhichWedges": "All",
            "TimeDependentMaps": "None",
            "OuterBoundaryCondition": boundary,
        }
    }
    return config


def logical_tensor(xyz):
    return tnsr.I[DataVector, 3, Frame.ElementLogical](
        np.array([x.ravel(order="F") for x in xyz])
    )


def geometry(element_map, axes):
    """Return coordinates, positive determinant and d xi^a / d x^i."""
    shape = tuple(len(x) for x in axes)
    logical = logical_tensor(np.meshgrid(*axes, indexing="ij"))
    xyz = np.asarray(element_map(logical)).reshape((3, *shape), order="F")
    jac = element_map.jacobian(logical)
    inverse = element_map.inv_jacobian(logical)
    jac = np.array([[jac.get(i, j) for j in range(3)] for i in range(3)])
    inverse = np.array(
        [[inverse.get(i, j) for j in range(3)] for i in range(3)]
    ).reshape((3, 3, *shape), order="F")
    determinant = np.linalg.det(np.moveaxis(jac, -1, 0)).reshape(
        shape, order="F"
    )
    if not np.all(determinant > 0):
        raise ValueError("Invalid or reversed element map")
    return xyz, determinant, inverse


def physical_gradient(modes, inverse, order):
    logical_gradient = np.array(
        [cube.evaluate(legder(modes, axis=d), order) for d in range(3)]
    )
    return np.einsum("ai...,a...->i...", inverse, logical_gradient)


def element_errors(modes, element_map, t, order, *, exact_solution=exact):
    xi, _, weights = cube.quadrature(order)
    xyz, determinant, inverse = geometry(element_map, [xi] * 3)
    measure = weights * determinant
    numerical = np.array([cube.evaluate(c, order) for c in modes])
    squared = np.sum(
        measure * (numerical - exact_solution(xyz, t)) ** 2, axis=(1, 2, 3)
    )
    constraint = physical_gradient(modes[0], inverse, order) - numerical[2:]
    return (
        squared,
        float(np.sum(measure * constraint**2)),
        float(np.sum(measure)),
    )


def read_domain(directory):
    with spectre_h5.H5File(str(directory / "Volume0.h5"), "r") as handle:
        return deserialize_domain[3](handle.get_vol("Fields").get_domain())


def groups(handle):
    return sorted(
        (g for g in handle["Fields.vol"].values() if isinstance(g, h5py.Group)),
        key=lambda g: g.attrs["observation_value"],
    )


def analyze_run(
    directory, *, exact_solution=exact, expected_volume=None, orders=ORDERS
):
    if len(orders) != 2 or not 1 <= orders[0] < orders[1]:
        raise ValueError("Specify two increasing quadrature orders")
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata["returncode"]:
        return cube.analyze_run(directory)
    result = {"run": directory.name, **metadata}
    domain = read_domain(directory)
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
    rows, meshes, fits = [], [], []
    coordinate_error = 0.0
    volume = (
        4 * np.pi / 3 * (OUTER_RADIUS**3 - INNER_RADIUS**3)
        if expected_volume is None
        else expected_volume
    )
    with h5py.File(directory / "Volume0.h5") as handle:
        for group in groups(handle):
            t = float(group.attrs["observation_value"])
            sums = {order: np.zeros(7) for order in orders}
            nelements, npoints = 0, 0
            for name, ext, _, xyz, fields in cube.elements(group):
                element_id = ElementId[3](name)
                element_map = ElementMap(element_id, domain)
                mapped, _, _ = geometry(
                    element_map, [legendre_grid(n)[0] for n in ext]
                )
                coordinate_error = max(
                    coordinate_error, float(np.max(abs(mapped - np.array(xyz))))
                )
                modes = [to_modal(v) for v in fields]
                for order in orders:
                    errors, constraint, measured_volume = element_errors(
                        modes,
                        element_map,
                        t,
                        order,
                        exact_solution=exact_solution,
                    )
                    sums[order] += np.array(
                        [*errors, constraint, measured_volume]
                    )
                nelements += 1
                npoints += int(np.prod(ext))
                meshes.append(
                    {
                        "time": t,
                        "element": name,
                        "block": element_id.block_id,
                        "physical_volume": measured_volume,
                        **{
                            f"level_{d}": level
                            for d, level in enumerate(
                                element_id.refinement_levels
                            )
                        },
                        **{
                            f"degree_{d}": int(n - 1) for d, n in enumerate(ext)
                        },
                    }
                )
                mesh = Mesh[3](ext, Basis.Legendre, Quadrature.GaussLobatto)
                for field, values in zip(cube.FIELDS, fields):
                    estimates = chalmers_energy_estimate(
                        DataVector(values.ravel(order="F")),
                        mesh,
                        metadata["tolerance"],
                    )
                    for d, fit in enumerate(estimates):
                        fits.append(
                            {
                                "time": t,
                                "element": name,
                                "field": field,
                                "direction": d,
                                **fit,
                            }
                        )
            high = np.sqrt(sums[orders[1]][:6])
            if not all(np.all(np.isfinite(values)) for values in sums.values()):
                raise RuntimeError(f"Nonfinite diagnostics at t={t}")
            delta = float(np.max(abs(high - np.sqrt(sums[orders[0]][:6]))))
            if delta > max(1e-12, 1e-5 * max(high)):
                raise RuntimeError(f"Unconverged quadrature at t={t}: {delta}")
            volume_error = abs(sums[orders[1]][6] - volume)
            if volume_error > 1e-10 * volume or coordinate_error > 1e-12:
                raise RuntimeError("Geometry verification failed")
            rows.append(
                {
                    "time": t,
                    "elements": nelements,
                    "points": npoints,
                    **{
                        f"error_{field}": float(e)
                        for field, e in zip(cube.FIELDS, high[:5])
                    },
                    "error_state": float(np.linalg.norm(high[:5])),
                    "constraint": float(high[5]),
                    "quadrature_change": delta,
                    "volume_error": float(volume_error),
                }
            )
    if not np.isclose(
        rows[-1]["time"], metadata["end_time"], atol=1e-12, rtol=0
    ):
        raise RuntimeError("Missing final observation")
    last_mesh = [r for r in meshes if r["time"] == rows[-1]["time"]]
    result.update(
        initial_error=rows[0]["error_state"],
        final_error=rows[-1]["error_state"],
        max_error=max(r["error_state"] for r in rows),
        final_constraint=rows[-1]["constraint"],
        max_constraint=max(r["constraint"] for r in rows),
        quadrature_change=max(r["quadrature_change"] for r in rows),
        volume_error=max(r["volume_error"] for r in rows),
        coordinate_error=coordinate_error,
        quadrature_orders=list(orders),
        final_max_degree=[
            max(r[f"degree_{d}"] for r in last_mesh) for d in range(3)
        ],
        final_max_level=[
            max(r[f"level_{d}"] for r in last_mesh) for d in range(3)
        ],
    )
    write_csv(directory / "History.csv", rows)
    write_csv(directory / "Mesh.csv", meshes)
    write_csv(directory / "Fits.csv", fits)
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def curved_difference(a, b, domain, order):
    """Integrate a polynomial difference on a common block-logical partition.

    Unlike affine boxes, the physical Jacobian varies, so the result requires
    a quadrature convergence check. Elements in different blocks never overlap.
    """
    xi, _, weights = cube.quadrature(order)
    if not a or not b:
        raise ValueError("Cannot compare empty meshes")
    components = len(next(iter(a.values())))
    if any(len(v) != components for v in [*a.values(), *b.values()]):
        raise ValueError("Field component counts differ")
    squared = np.zeros(components)
    covered = {}
    for name_a, modes_a in a.items():
        id_a = ElementId[3](name_a)
        bounds_a = [s.endpoints for s in id_a.segment_ids]
        for name_b, modes_b in b.items():
            id_b = ElementId[3](name_b)
            if id_a.block_id != id_b.block_id:
                continue
            bounds_b = [s.endpoints for s in id_b.segment_ids]
            overlap = [
                (max(a[0], b[0]), min(a[1], b[1]))
                for a, b in zip(bounds_a, bounds_b)
            ]
            if any(hi <= lo for lo, hi in overlap):
                continue
            scale = np.prod([(hi - lo) / 2 for lo, hi in overlap])
            covered[id_a.block_id] = covered.get(id_a.block_id, 0.0) + 8 * scale
            axes = [(lo + hi) / 2 + (hi - lo) / 2 * xi for lo, hi in overlap]
            block_map = ElementMap(ElementId[3](id_a.block_id), domain)
            _, determinant, _ = geometry(block_map, axes)

            def restrict(modes, bounds):
                result = modes
                for d, ((lo, hi), nodes) in enumerate(zip(bounds, axes)):
                    local = (2 * nodes - lo - hi) / (hi - lo)
                    matrix = legvander(local, modes.shape[d] - 1)
                    result = np.moveaxis(
                        np.tensordot(matrix, result, axes=(1, d)), 0, d
                    )
                return result

            for f, (ma, mb) in enumerate(zip(modes_a, modes_b)):
                delta = restrict(ma, bounds_a) - restrict(mb, bounds_b)
                squared[f] += scale * np.sum(weights * determinant * delta**2)
    blocks_a = {ElementId[3](n).block_id for n in a}
    blocks_b = {ElementId[3](n).block_id for n in b}
    if (
        blocks_a != blocks_b
        or set(covered) != blocks_a
        or any(abs(v - 8) > 1e-12 for v in covered.values())
    ):
        raise ValueError("Meshes do not cover the same complete blocks")
    return np.sqrt(squared)


def compare(args, *, orders=ORDERS):
    domains = [read_domain(d) for d in [args.first, args.second]]
    if domains[0] != domains[1]:
        raise ValueError("Domains differ")

    def final(directory):
        with h5py.File(directory / "Volume0.h5") as handle:
            group = groups(handle)[-1]
            return float(group.attrs["observation_value"]), {
                name: [to_modal(v) for v in fields]
                for name, _, _, _, fields in cube.elements(group)
            }

    ta, a = final(args.first)
    tb, b = final(args.second)
    if abs(ta - tb) > 1e-12:
        raise ValueError("Final times differ")
    low, high = [curved_difference(a, b, domains[0], order) for order in orders]
    change = float(np.max(abs(high - low)))
    if change > max(1e-12, 1e-5 * max(high)):
        raise RuntimeError("Difference quadrature not converged")
    result = {
        "first": str(args.first),
        "second": str(args.second),
        "same_h_mesh": a.keys() == b.keys(),
        "same_hp_mesh": (
            a.keys() == b.keys()
            and all(a[n][0].shape == b[n][0].shape for n in a)
        ),
        "field_l2_differences": high.tolist(),
        "state_l2_difference": float(np.linalg.norm(high)),
        "quadrature_change": change,
        "quadrature_orders": list(orders),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def summarize(args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.patches import Rectangle

    results = []
    for directory in args.runs:
        path = directory / "Analysis.json"
        record = (
            json.loads(path.read_text())
            if path.exists()
            else analyze_run(directory)
        )
        results.append({"directory": str(directory.resolve()), **record})
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "Summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    columns = [
        "run",
        "strategy",
        "tolerance",
        "fixed_p",
        "cfl",
        "returncode",
        "final_elements",
        "final_points",
        "initial_error",
        "max_error",
        "final_error",
        "max_constraint",
        "final_constraint",
        "evolution_seconds",
    ]
    write_csv(
        args.output / "RunSummary.csv",
        [{key: record.get(key, "") for key in columns} for record in results],
    )
    main = [
        r
        for r in results
        if r.get("returncode") == 0
        and r["stepping"] == "cfl"
        and r["cfl"] == 0.5
    ]
    styles = {
        "fixed": ("Fixed grids", "#555555"),
        "p": ("p only", "#2271b2"),
        "hp05": ("hp, sigma 0.5", "#009e73"),
        "hp1": ("hp, sigma 1", "#d55e00"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for strategy, (label, color) in styles.items():
        selected = sorted(
            [r for r in main if r["strategy"] == strategy],
            key=lambda r: r["max_error"],
        )
        for axis, cost in zip(axes, ["evolution_seconds", "final_points"]):
            axis.loglog(
                [r[cost] for r in selected],
                [r["max_error"] for r in selected],
                {"fixed": "^-", "p": "o-", "hp05": "s-", "hp1": "x--"}[
                    strategy
                ],
                label=label,
                color=color,
            )
    for axis, label in zip(
        axes, ["Evolution wall time [s]", "Final grid points"]
    ):
        axis.set(
            xlabel=label,
            ylabel="Maximum observed physical L2 error",
            title="Outgoing pulse on curved shells",
        )
        axis.set_xscale("linear")
        axis.set_xlim(left=0)
        axis.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
        axis.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:,.0f}"))
        axis.grid(alpha=0.25, which="both")
        axis.legend(fontsize=8)
    fig.savefig(args.output / "AccuracyCost.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(
        1, 3, figsize=(12, 4), constrained_layout=True, sharey=True
    )
    mesh_target = getattr(args, "mesh_target", 1e-6)
    fig.suptitle(f"Radial refinement at local target {mesh_target:g}")
    selected = [
        next(
            (
                r
                for r in main
                if r["strategy"] == strategy and r["tolerance"] == mesh_target
            ),
            None,
        )
        for strategy in ["p", "hp05", "hp1"]
    ]
    norm = matplotlib.colors.Normalize(vmin=5, vmax=20)
    cmap = matplotlib.colormaps["viridis"]
    for axis, record in zip(axes, selected):
        if record is None:
            axis.set_visible(False)
            continue
        directory = Path(record["directory"])
        domain = read_domain(directory)
        with (directory / "Mesh.csv").open() as handle:
            for row in csv.DictReader(handle):
                element_id = ElementId[3](row["element"])
                # Show the +z wedge of each shell. Angular refinements, if
                # present, would overlap here, so reject this visualization.
                if element_id.block_id % 6:
                    continue
                if any(element_id.refinement_levels[:2]):
                    raise ValueError(
                        "Radial-only mesh plot requires unrefined angles"
                    )
                xyz, _, _ = geometry(
                    ElementMap(element_id, domain),
                    [np.array([0.0]), np.array([0.0]), np.array([-1.0, 1.0])],
                )
                lo, hi = np.sqrt(np.sum(xyz**2, axis=0)).ravel()
                axis.add_patch(
                    Rectangle(
                        (lo, float(row["time"]) - 0.075),
                        hi - lo,
                        0.15,
                        facecolor=cmap(norm(int(row["degree_2"]))),
                        edgecolor="black",
                        linewidth=0.5,
                    )
                )
        times = np.linspace(0, record["end_time"], 100)
        axis.plot(
            CENTER + times, times, "r--", linewidth=1, label="Pulse center"
        )
        axis.set(
            xlim=(INNER_RADIUS, OUTER_RADIUS),
            ylim=(-0.12, record["end_time"] + 0.12),
            xlabel="Radius",
            title=styles[record["strategy"]][0],
        )
        axis.legend(fontsize=8)
    axes[0].set_ylabel("Observation time")
    axes[0].yaxis.set_major_locator(
        ticker.MultipleLocator(cube.OUTPUT_INTERVAL)
    )
    fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes,
        label="Radial polynomial degree",
    )
    fig.savefig(args.output / "RadialRefinement.png", dpi=180)
    plt.close(fig)

    matched = []
    for target in [1e-3, 1e-4, 1e-5]:
        for strategy in styles:
            candidates = [
                r
                for r in main
                if r["strategy"] == strategy and r["max_error"] <= target
            ]
            best = (
                min(candidates, key=lambda r: r["evolution_seconds"])
                if candidates
                else None
            )
            matched.append(
                {
                    "required_max_observed_error": target,
                    "strategy": strategy,
                    "run": best["run"] if best else "not attained",
                    "seconds": best["evolution_seconds"] if best else "",
                    "points": best["final_points"] if best else "",
                }
            )
    write_csv(args.output / "MatchedAccuracy.csv", matched)
    return {
        "runs": len(results),
        "complete": sum(r.get("returncode") == 0 for r in results),
        "output": str(args.output),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    runner = commands.add_parser("run")
    runner.add_argument("--executable", type=Path, required=True)
    runner.add_argument("--output", type=Path, required=True)
    runner.add_argument(
        "--strategy", choices=["fixed", "p", "hp1", "hp05"], required=True
    )
    runner.add_argument("--tolerance", type=float, default=1e-4)
    runner.add_argument("--fixed-p", type=int, default=13)
    runner.add_argument("--max-level", type=int, default=3)
    runner.add_argument("--end-time", type=float, default=1.6)
    runner.add_argument("--observation-interval", type=float, default=0.4)
    runner.add_argument("--dt", type=float, default=0.0005)
    runner.add_argument("--stepping", choices=["fixed", "cfl"], default="cfl")
    runner.add_argument("--cfl", type=float, default=0.5)
    runner.add_argument("--coarsening-factor", type=float, default=0.0)
    runner.add_argument("--h-coarsening", action="store_true")
    runner.add_argument("--processes", type=int, default=1)
    runner.add_argument("--timeout", type=float, default=900)
    runner.set_defaults(
        func=lambda a: cube.run(
            a, config=configuration(a), analyzer=analyze_run
        )
    )
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("run", type=Path)
    analyzer.set_defaults(func=lambda a: analyze_run(a.run))
    differ = commands.add_parser("compare")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.set_defaults(func=compare)
    summary = commands.add_parser("summarize")
    summary.add_argument("runs", nargs="+", type=Path)
    summary.add_argument("--output", type=Path, required=True)
    summary.add_argument("--mesh-target", type=float, default=1e-6)
    summary.set_defaults(func=summarize)
    args = parser.parse_args()
    result = args.func(args)
    if isinstance(result, int):
        if result:
            raise SystemExit(1)
    else:
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
