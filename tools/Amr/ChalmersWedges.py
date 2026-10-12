# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Oblique Gaussian pulse crossing angular interfaces of a cubed sphere.

The shell is 1 <= r <= 2, with six equiangular wedges. The pulse travels in
direction (1,2,3)/sqrt(14), with width 0.8 and initial phase center -2.8.
Use python-spectre; every native run writes to a new, independent directory.
"""

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import leggauss, legval3d

from spectre.Domain import ElementId, ElementMap
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv

INNER_RADIUS = 1.0
OUTER_RADIUS = 2.0
WIDTH = 0.8
CENTER = -2.8
NORMAL = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
INITIAL_EXTENTS = [[n, n, 14] for n in [10, 18, 14, 18, 14, 16]]
EXPERIMENT = "oblique-cubed-sphere-v1"


def exact(xyz, t):
    u = np.einsum("i,i...->...", NORMAL, xyz) - t - CENTER
    psi = np.exp(-((u / WIDTH) ** 2))
    derivative = -2 * u / WIDTH**2 * psi
    return np.array([psi, derivative, *(n * derivative for n in NORMAL)])


def configuration(args):
    config = sphere.configuration(args)
    data = {
        "PlaneWave": {
            "WaveVector": NORMAL.tolist(),
            "Center": (CENTER * NORMAL).tolist(),
            "Profile": {
                "Gaussian": {"Amplitude": 1.0, "Width": WIDTH, "Center": 0.0}
            },
        }
    }
    boundary = {"DirichletAnalytic": {"AnalyticPrescription": data}}
    config["InitialData"] = data
    creator = config["DomainCreator"]["Sphere"]
    creator.update(
        InnerRadius=INNER_RADIUS,
        OuterRadius=OUTER_RADIUS,
        RadialPartitioning=[],
        Interior={"ExciseWithBoundaryCondition": boundary},
        OuterBoundaryCondition=boundary,
        InitialGridPoints=(
            [args.fixed_p + 1] * 3
            if args.strategy == "fixed"
            else INITIAL_EXTENTS
        ),
    )
    return config


def analyze_run(directory):
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata.get("experiment") != EXPERIMENT:
        raise ValueError(
            "This analyzer requires the oblique cubed-sphere experiment"
        )
    result = sphere.analyze_run(
        directory,
        exact_solution=exact,
        expected_volume=4 * np.pi / 3 * (OUTER_RADIUS**3 - INNER_RADIUS**3),
        orders=(32, 40),
    )
    if metadata["returncode"]:
        result["startup_failure"] = (
            metadata["returncode"] == -11
            and not (directory / "Reductions.h5").exists()
            and not (directory / "Volume0.h5").exists()
        )
        if result["startup_failure"]:
            result["failure_reason"] = "Signal 11 before any evolution output"
        (directory / "Analysis.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result
    domain = sphere.read_domain(directory)
    rows = []
    with h5py.File(directory / "Volume0.h5") as handle:
        for group in sphere.groups(handle):
            t = float(group.attrs["observation_value"])
            modes = {
                name: [to_modal(v) for v in fields]
                for name, _, _, _, fields in cube.elements(group)
                if ElementId[3](name).block_id in [3, 4]
            }
            low, high = [
                interface_errors(modes, domain, t, order)
                for order in sphere.ORDERS
            ]
            delta = max(abs(low[key] - high[key]) for key in high)
            if delta > 1e-10:
                raise RuntimeError(
                    f"Interface quadrature not converged: {delta}"
                )
            rows.append({"time": t, **high, "quadrature_change": delta})
    result.update(
        max_interface_error=max(
            max(r["error_left"], r["error_right"]) for r in rows
        ),
        max_interface_jump=max(r["jump"] for r in rows),
        interface_quadrature_change=max(r["quadrature_change"] for r in rows),
    )
    write_csv(directory / "Interface.csv", rows)
    (directory / "Analysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def interface_fields(modes, domain, order):
    """Sample LowerY -> UpperX (x=-y>0) on a common face partition.

    These blocks share tangential logical directions 1 and 2. Direction 0
    is +1 on LowerY and -1 on UpperX. Integrating each element-face overlap
    separately resolves jumps at transverse h boundaries; whole-face Gauss
    quadrature would not reliably integrate these piecewise polynomials.
    """
    xi, w = leggauss(order)
    faces = {3: [], 4: []}
    for name, fields in modes.items():
        element_id = ElementId[3](name)
        block = element_id.block_id
        if block not in faces:
            continue
        bounds = np.array([s.endpoints for s in element_id.segment_ids])
        if (block == 3 and bounds[0, 1] == 1) or (
            block == 4 and bounds[0, 0] == -1
        ):
            faces[block].append((bounds, fields))
    coordinates, measures, left_values, right_values = [], [], [], []
    covered = 0.0
    for bounds_a, fields_a in faces[3]:
        for bounds_b, fields_b in faces[4]:
            overlap = [
                (max(a[0], b[0]), min(a[1], b[1]))
                for a, b in zip(bounds_a[1:], bounds_b[1:])
            ]
            if any(hi <= lo for lo, hi in overlap):
                continue
            patch_scale = np.prod([(hi - lo) / 2 for lo, hi in overlap])
            covered += 4 * patch_scale
            angle, radial = np.meshgrid(
                *[(lo + hi) / 2 + (hi - lo) / 2 * xi for lo, hi in overlap],
                indexing="ij",
            )
            radius = (INNER_RADIUS + OUTER_RADIUS) / 2 + (
                OUTER_RADIUS - INNER_RADIUS
            ) / 2 * radial
            zeta = np.tan(np.pi / 4 * angle)
            x = radius / np.sqrt(2 + zeta**2)
            xyz = np.array([x, -x, zeta * x])
            surface_jacobian = (
                radius
                * np.sqrt(2)
                / (2 + zeta**2)
                * (OUTER_RADIUS - INNER_RADIUS)
                / 2
                * np.pi
                / 4
                * (1 + zeta**2)
            )
            measures.append(
                (
                    patch_scale * w[:, None] * w[None, :] * surface_jacobian
                ).ravel(order="F")
            )
            coordinates.append(xyz.reshape(3, -1, order="F"))
            for block, face, bounds, fields, destination in [
                (3, 1.0, bounds_a, fields_a, left_values),
                (4, -1.0, bounds_b, fields_b, right_values),
            ]:
                logical = np.array([np.full_like(angle, face), angle, radial])
                mapped = np.asarray(
                    ElementMap(ElementId[3](block), domain)(
                        sphere.logical_tensor(logical)
                    )
                ).reshape(xyz.shape, order="F")
                if np.max(abs(mapped - xyz)) > 1e-12:
                    raise ValueError(
                        "Interface does not match the saved cubed-sphere map"
                    )
                logical = logical.reshape(3, -1, order="F")
                local = (2 * logical - bounds.sum(axis=1)[:, None]) / np.diff(
                    bounds, axis=1
                )
                destination.append(
                    np.array([legval3d(*local, c) for c in fields])
                )
    if abs(covered - 4) > 1e-12:
        raise ValueError(
            "Meshes do not cover the complete interface exactly once"
        )
    return (
        np.concatenate(coordinates, axis=1),
        np.concatenate(measures),
        [
            np.concatenate(left_values, axis=1),
            np.concatenate(right_values, axis=1),
        ],
    )


def interface_errors(modes, domain, t, order):
    xyz, weights, (left, right) = interface_fields(modes, domain, order)
    analytic = exact(xyz, t)
    normal = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)

    def flux(state):
        return float(
            np.sum(
                weights * state[1] * np.einsum("i,i...->...", normal, state[2:])
            )
        )

    def norm(values):
        return float(np.sqrt(np.sum(weights * values**2)))

    return {
        "error_left": norm(left - analytic),
        "error_right": norm(right - analytic),
        "jump": norm(left - right),
        "analytic_flux": flux(analytic),
        "flux_left": flux(left),
        "flux_right": flux(right),
    }


def summarize(args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import ticker

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
        "max_interface_error",
        "max_interface_jump",
        "evolution_seconds",
        "timing_eligible",
        "failure_reason",
    ]
    write_csv(
        args.output / "RunSummary.csv",
        [{key: r.get(key, "") for key in columns} for r in results],
    )
    main = [
        r
        for r in results
        if r.get("returncode") == 0
        and r["stepping"] == "cfl"
        and r["cfl"] == 0.5
    ]
    styles = {
        "fixed": ("Fixed grids", "#555555", "^-"),
        "p": ("p only", "#2271b2", "o-"),
        "hp05": ("hp, sigma 0.5", "#009e73", "s-"),
        "hp1": ("hp, sigma 1", "#d55e00", "x--"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for strategy, (label, color, marker) in styles.items():
        selected = sorted(
            [r for r in main if r["strategy"] == strategy],
            key=lambda r: r["max_error"],
        )
        for axis, cost in zip(axes, ["evolution_seconds", "final_points"]):
            samples = [
                r
                for r in selected
                if cost != "evolution_seconds" or r.get("timing_eligible", True)
            ]
            axis.semilogy(
                [r[cost] for r in samples],
                [r["max_error"] for r in samples],
                marker,
                label=label,
                color=color,
            )
    for axis, label in zip(
        axes, ["Evolution wall time [s]", "Final grid points"]
    ):
        axis.set(
            xlabel=label,
            ylabel="Maximum observed physical L2 error",
            title="Oblique pulse on curved wedges",
            xlim=(0, None),
        )
        axis.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
        axis.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:,.0f}"))
        axis.grid(alpha=0.25, which="both")
        axis.legend(fontsize=8)
    failed = sum(r.get("returncode", 0) != 0 for r in results)
    if failed:
        startup_failures = sum(r.get("startup_failure", False) for r in results)
        fig.suptitle(
            f"{startup_failures} startup failures; {failed - startup_failures}"
            " incomplete evolutions. Curves show completed runs."
        )
    fig.savefig(args.output / "AccuracyCost.png", dpi=180)
    plt.close(fig)

    def csv_rows(path):
        with path.open() as handle:
            return list(csv.DictReader(handle))

    selected = [
        next(
            (
                r
                for r in main
                if r["strategy"] == strategy
                and r["tolerance"] == args.mesh_target
            ),
            None,
        )
        for strategy in ["p", "hp05", "hp1"]
    ]
    fig, axes = plt.subplots(
        3, 1, figsize=(8, 8), constrained_layout=True, sharex=True
    )
    histories = []
    for record in selected:
        if record is None:
            continue
        directory = Path(record["directory"])
        mesh = csv_rows(directory / "Mesh.csv")
        history = csv_rows(directory / "History.csv")
        label, color, marker = styles[record["strategy"]]
        times, points, angular_h, radial_h, angular_p = [], [], [], [], []
        for entry in history:
            t = float(entry["time"])
            cells = [r for r in mesh if float(r["time"]) == t]
            volume = sum(float(r["physical_volume"]) for r in cells)
            hp = max(max(int(r["level_0"]), int(r["level_1"])) for r in cells)
            rp = max(int(r["level_2"]) for r in cells)
            degree = (
                sum(
                    float(r["physical_volume"])
                    * (int(r["degree_0"]) + int(r["degree_1"]))
                    / 2
                    for r in cells
                )
                / volume
            )
            times.append(t)
            points.append(int(entry["points"]))
            angular_h.append(hp)
            radial_h.append(rp)
            angular_p.append(degree)
            histories.append(
                {
                    "run": record["run"],
                    "time": t,
                    "points": points[-1],
                    "max_angular_level": hp,
                    "max_radial_level": rp,
                    "volume_weighted_angular_degree": degree,
                }
            )
        for axis, values in zip(axes, [points, angular_h, angular_p]):
            axis.plot(times, values, marker, label=label, color=color)
        axes[1].plot(times, radial_h, ":", color=color, alpha=0.6)
    axes[0].set(
        ylabel="Grid points",
        title=f"Adaptation at local target {args.mesh_target:g}",
    )
    axes[1].set_ylabel("Maximum h level\nsolid/marker: angular; dotted: radial")
    axes[1].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    axes[2].set(
        ylabel="Mean angular degree\n(weighted by physical volume)",
        xlabel="Time",
    )
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend()
    fig.savefig(args.output / "DirectionalRefinement.png", dpi=180)
    plt.close(fig)
    write_csv(args.output / "DirectionalHistory.csv", histories)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for record in selected:
        if record is None:
            continue
        rows = csv_rows(Path(record["directory"]) / "Interface.csv")
        t = [float(r["time"]) for r in rows]
        label, color, marker = styles[record["strategy"]]
        axes[1].semilogy(
            t,
            [float(r["jump"]) for r in rows],
            marker,
            color=color,
            label=label,
        )
        if record["strategy"] == "hp05":
            axes[0].plot(
                t,
                [float(r["analytic_flux"]) for r in rows],
                "k-",
                label="Exact",
            )
            axes[0].plot(
                t,
                [
                    (float(r["flux_left"]) + float(r["flux_right"])) / 2
                    for r in rows
                ],
                "s--",
                color=color,
                label="hp, sigma 0.5 (mean of traces)",
            )
    axes[0].set(
        xlabel="Time",
        ylabel="Integrated physical energy flux",
        title="Transport from LowerY to UpperX",
    )
    axes[1].set(
        xlabel="Time",
        ylabel="Physical L2 norm of the trace jump",
        title=f"Angular interface, target {args.mesh_target:g}",
    )
    for axis in axes:
        axis.legend(fontsize=8)
        axis.grid(alpha=0.25)
    fig.savefig(args.output / "InterfaceTransport.png", dpi=180)
    plt.close(fig)

    matched = []
    for target in [1e-3, 1e-4, 1e-5]:
        for strategy in styles:
            candidates = [
                r
                for r in main
                if r["strategy"] == strategy
                and r["max_error"] <= target
                and r.get("timing_eligible", True)
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
        "complete": len(results) - failed,
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
    runner.add_argument("--fixed-p", type=int, default=15)
    runner.add_argument("--max-level", type=int, default=3)
    runner.add_argument("--end-time", type=float, default=2.4)
    runner.add_argument("--observation-interval", type=float, default=0.4)
    runner.add_argument("--dt", type=float, default=0.0005)
    runner.add_argument("--stepping", choices=["fixed", "cfl"], default="cfl")
    runner.add_argument("--cfl", type=float, default=0.5)
    runner.add_argument("--timeout", type=float, default=900)
    runner.set_defaults(
        experiment=EXPERIMENT,
        func=lambda a: cube.run(
            a, config=configuration(a), analyzer=analyze_run
        ),
    )
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("run", type=Path)
    analyzer.set_defaults(func=lambda a: analyze_run(a.run))
    differ = commands.add_parser("compare")
    differ.add_argument("first", type=Path)
    differ.add_argument("second", type=Path)
    differ.add_argument("--output", type=Path, required=True)
    differ.set_defaults(func=lambda a: sphere.compare(a, orders=(32, 40)))
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
