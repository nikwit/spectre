# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Reconstruct spherical Cartesian GH tensors from a fixed Cartoon solution."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legvander

from spectre.Domain import ElementId, ElementMap
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSingleHole as gh
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv


def index(a, b):
    return gh.PAIRS.index(tuple(sorted((a, b))))


def coefficients(axis):
    """The 14 parity-even spherical coefficients of g, Pi and Phi."""
    result = []
    for offset in (0, 10):
        result.extend(
            [
                axis[offset + index(0, 0)],
                axis[offset + index(0, 1)],
                axis[offset + index(2, 2)],
                axis[offset + index(1, 1)] - axis[offset + index(2, 2)],
            ]
        )
    phi = lambda k, a, b: axis[20 + 10 * k + index(a, b)]
    a, b, d, e = phi(0, 0, 0), phi(1, 0, 2), phi(0, 2, 2), phi(1, 1, 2)
    result.extend([a, b, phi(0, 0, 1) - b, d, e, phi(0, 1, 1) - d - 2 * e])
    return np.array(result)


def reconstruct(c, n):
    """Lift the radial functions to Cartesian components in directions n."""
    result = []
    for offset in (0, 4):
        for a, b in gh.PAIRS:
            if a == b == 0:
                v = c[offset]
            elif a == 0:
                v = c[offset + 1] * n[b - 1]
            else:
                v = (
                    c[offset + 2] * (a == b)
                    + c[offset + 3] * n[a - 1] * n[b - 1]
                )
            result.append(v)
    for k in range(3):
        for a, b in gh.PAIRS:
            if a == b == 0:
                v = c[8] * n[k]
            elif a == 0:
                v = c[9] * (k == b - 1) + c[10] * n[k] * n[b - 1]
            else:
                i, j = a - 1, b - 1
                v = (
                    c[11] * n[k] * (i == j)
                    + c[12] * ((k == i) * n[j] + (k == j) * n[i])
                    + c[13] * n[k] * n[i] * n[j]
                )
            result.append(v)
    return np.array(result)


def read(directory, time):
    metadata = json.loads((directory / "Run.json").read_text())
    if (
        metadata["returncode"]
        or metadata["experiment"]
        != "schwarzschild-damped-harmonic-cartoon-reference-v1"
    ):
        raise ValueError("Need a successful spherical reference")
    cells = []
    with h5py.File(directory / "Volume0.h5") as f:
        group = next(
            g for g in sphere.groups(f) if g.attrs["observation_value"] == time
        )
        for name, ext, xyz, fields in gh.elements(group):
            assert tuple(ext[1:]) == (1, 1)
            r = xyz[0, :, 0, 0]
            if np.max(abs(xyz[1:])) > 1e-13 or not np.all(np.diff(r) > 0):
                raise ValueError("Reference must lie on the positive x axis")
            lo, hi = float(r[0]), float(r[-1])
            xi = (2 * np.log(r) - np.log(lo) - np.log(hi)) / np.log(hi / lo)
            np.testing.assert_allclose(
                xi, legendre_grid(len(r))[0], atol=1e-12, rtol=0
            )
            axis = fields[:, :, 0, 0]
            c = coefficients(axis)
            reconstructed = reconstruct(
                c,
                np.array([np.ones_like(r), np.zeros_like(r), np.zeros_like(r)]),
            )
            residual = float(np.max(abs(reconstructed - axis)))
            if not np.all(np.isfinite(axis)) or residual > 1e-12:
                raise ValueError(
                    f"Broken spherical symmetry: {name}, {residual}"
                )
            eid = ElementId[3](name)
            cells.append(
                dict(
                    name=name,
                    lo=lo,
                    hi=hi,
                    block=eid.block_id,
                    logical_bounds=eid.segment_ids[0].endpoints,
                    modes=np.array([to_modal(v) for v in c]),
                    symmetry_residual=residual,
                )
            )
    cells.sort(key=lambda c: c["lo"])
    np.testing.assert_allclose(
        [cells[0]["lo"], cells[-1]["hi"]], [1.2, 20.0], atol=1e-12, rtol=0
    )
    np.testing.assert_allclose(
        [c["hi"] for c in cells[:-1]],
        [c["lo"] for c in cells[1:]],
        atol=1e-12,
        rtol=0,
    )
    return cells


def evaluate(cell, r, derivative=False):
    xi = (2 * np.log(r) - np.log(cell["lo"]) - np.log(cell["hi"])) / np.log(
        cell["hi"] / cell["lo"]
    )
    if np.max(abs(xi)) > 1 + 1e-11:
        raise ValueError("Reference extrapolation requested")
    modes = cell["modes"]
    if derivative:
        modes = legder(modes, axis=1)
    values = np.moveaxis(legvander(xi, modes.shape[1] - 1) @ modes.T, -1, 0)
    return (
        values * (2 / (r * np.log(cell["hi"] / cell["lo"])))
        if derivative
        else values
    )


def reference_diagnostics(directory):
    with h5py.File(directory / "Volume0.h5") as f:
        times = [float(g.attrs["observation_value"]) for g in sphere.groups(f)]
    rows = []
    for time in times:
        cells = read(directory, time)
        integrals = []
        max_mass_error = 0.0
        for order in (64, 96):
            total = 0.0
            for cell in cells:
                xi, w = leggauss(order)
                width = np.log(cell["hi"] / cell["lo"])
                r = np.exp(
                    (np.log(cell["lo"]) + np.log(cell["hi"])) / 2
                    + width * xi / 2
                )
                c, dc = evaluate(cell, r), evaluate(cell, r, True)
                expected = c.copy()
                expected[8] = dc[0]
                expected[9] = c[1] / r
                expected[10] = dc[1] - c[1] / r
                expected[11] = dc[2]
                expected[12] = c[3] / r
                expected[13] = dc[3] - 2 * c[3] / r
                n = np.array(
                    [np.ones_like(r), np.zeros_like(r), np.zeros_like(r)]
                )
                difference = (reconstruct(expected, n) - reconstruct(c, n))[20:]
                measure = 4 * np.pi * r**2 * (r * width / 2) * w
                total += np.sum(
                    measure
                    * np.sum(difference**2 * gh.MULTIPLICITY[20:, None], axis=0)
                )
                # Spherical vacuum invariant: M=R/2 (1-g^ab d_a R d_b R),
                # where R=r sqrt(g_yy) on the positive x axis.
                A, B, C, D = c[:4]
                grr = C + D
                beta = B / grr
                lapse = np.sqrt(-A + B * beta)
                R = r * np.sqrt(C)
                drR = np.sqrt(C) + r * dc[2] / (2 * np.sqrt(C))
                dtR = r * (-lapse * c[6] + beta * c[11]) / (2 * np.sqrt(C))
                gradient_squared = (
                    grr * dtR**2 - 2 * B * dtR * drR + A * drR**2
                ) / (A * grr - B**2)
                mass = R / 2 * (1 - gradient_squared)
                if not np.all(np.isfinite(mass)):
                    raise ValueError("Nonfinite spherical mass diagnostic")
                max_mass_error = max(
                    max_mass_error, float(np.max(abs(mass - 1)))
                )
            integrals.append(float(np.sqrt(total)))
        if not np.all(np.isfinite(integrals)) or abs(
            integrals[1] - integrals[0]
        ) > max(1e-11, 1e-5 * integrals[1]):
            raise ValueError(
                "Reference constraint quadrature has not converged"
            )
        rows.append(
            dict(
                time=time,
                three_index_constraint=integrals[1],
                quadrature_change=abs(integrals[1] - integrals[0]),
                max_sampled_spherical_mass_error=max_mass_error,
                max_symmetry_residual=max(
                    c["symmetry_residual"] for c in cells
                ),
            )
        )
    write_csv(directory / "ReferenceHistory.csv", rows)
    return rows


def compare_references(first, second, output):
    with h5py.File(first / "Volume0.h5") as f:
        times = [float(g.attrs["observation_value"]) for g in sphere.groups(f)]
    rows = []
    for time in times:
        a, b = read(first, time), read(second, time)
        sums = []
        for order in (64, 96):
            xi, w = leggauss(order)
            total = 0.0
            for ca in a:
                for cb in b:
                    lo, hi = max(ca["lo"], cb["lo"]), min(ca["hi"], cb["hi"])
                    if hi - lo < 1e-12:
                        continue
                    width = np.log(hi / lo)
                    r = np.exp((np.log(lo) + np.log(hi)) / 2 + width * xi / 2)
                    n = np.array(
                        [np.ones_like(r), np.zeros_like(r), np.zeros_like(r)]
                    )
                    delta = reconstruct(evaluate(ca, r) - evaluate(cb, r), n)
                    total += np.sum(
                        4
                        * np.pi
                        * r**2
                        * (r * width / 2)
                        * w
                        * np.sum(delta**2 * gh.MULTIPLICITY[:, None], axis=0)
                    )
            sums.append(float(np.sqrt(total)))
        if not np.all(np.isfinite(sums)) or abs(sums[1] - sums[0]) > max(
            1e-11, 1e-5 * sums[1]
        ):
            raise ValueError(
                "Reference difference quadrature has not converged"
            )
        rows.append(
            dict(
                time=time,
                state_l2_difference=sums[1],
                quadrature_change=abs(sums[1] - sums[0]),
            )
        )
    result = dict(
        first=str(first),
        second=str(second),
        observations=rows,
        quadrature_orders=[64, 96],
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def tensor_evaluate(modes, axes):
    values = np.asarray(modes)
    for d, nodes in enumerate(axes, 1):
        matrix = legvander(nodes, values.shape[d] - 1)
        values = np.moveaxis(np.tensordot(matrix, values, axes=(1, d)), 0, d)
    return values


def compare_3d(directory, reference, output, selected_times=None):
    """Integrate on a common radial partition with exact spherical lifting.

    Both creators use identical logarithmic shells. The 3D angular element
    partition is retained; radial intersections with the Cartoon cells are
    integrated separately. Higher radial quadrature resolves the reference's
    radial polynomial order without unnecessarily enlarging angular grids.
    """
    meta = json.loads((directory / "Run.json").read_text())
    if (
        meta["returncode"]
        or meta["experiment"] != "schwarzschild-damped-harmonic-relaxation-v1"
    ):
        raise ValueError("Need a completed three-dimensional relaxation run")
    domain = sphere.read_domain(directory)
    rows = []
    orders = [(20, 40), (28, 48)]
    with h5py.File(directory / "Volume0.h5") as f:
        for group in sphere.groups(f):
            time = float(group.attrs["observation_value"])
            if selected_times is not None and time not in selected_times:
                continue
            cells = read(reference, time)
            squared = np.zeros((2, 50))
            volumes = np.zeros(2)
            for name, ext, xyz, fields in gh.elements(group):
                eid = ElementId[3](name)
                element_map = ElementMap(eid, domain)
                shell = eid.block_id // 6
                matching = [c for c in cells if c["block"] == shell]
                if not matching:
                    raise ValueError("No corresponding reference shell")
                lo, hi = eid.segment_ids[2].endpoints
                shell_lo, shell_hi = min(c["lo"] for c in matching), max(
                    c["hi"] for c in matching
                )
                boundaries, _, _ = sphere.geometry(
                    element_map,
                    [np.array([0.0]), np.array([0.0]), np.array([-1.0, 1.0])],
                )
                expected = np.exp(
                    (np.log(shell_lo) + np.log(shell_hi)) / 2
                    + np.log(shell_hi / shell_lo) * np.array([lo, hi]) / 2
                )
                np.testing.assert_allclose(
                    np.linalg.norm(boundaries, axis=0).ravel(),
                    expected,
                    rtol=0,
                    atol=1e-11,
                )
                modes = np.array([to_modal(v) for v in fields])
                for cell in matching:
                    a, b = max(lo, cell["logical_bounds"][0]), min(
                        hi, cell["logical_bounds"][1]
                    )
                    if b <= a:
                        continue
                    scale = (b - a) / (hi - lo)
                    for q, (na, nr) in enumerate(orders):
                        angular, wa = leggauss(na)
                        radial, wr = leggauss(nr)
                        block_radial = (a + b) / 2 + (b - a) * radial / 2
                        axes = [
                            angular,
                            angular,
                            (2 * block_radial - lo - hi) / (hi - lo),
                        ]
                        xyz, determinant, _ = sphere.geometry(element_map, axes)
                        r = np.linalg.norm(xyz, axis=0)
                        rline = r[0, 0, :]
                        if np.max(abs(r - rline)) > 1e-11:
                            raise ValueError("Radial map depends on angle")
                        c = evaluate(cell, rline)
                        c = np.broadcast_to(c[:, None, None, :], (14, *r.shape))
                        expected = reconstruct(c, xyz / r)
                        numerical = tensor_evaluate(modes, axes)
                        measure = (
                            scale
                            * determinant
                            * wa[:, None, None]
                            * wa[None, :, None]
                            * wr[None, None, :]
                        )
                        squared[q] += np.sum(
                            measure * (numerical - expected) ** 2,
                            axis=(1, 2, 3),
                        )
                        volumes[q] += np.sum(measure)
            expected_volume = 4 * np.pi / 3 * (20**3 - 1.2**3)
            if np.max(abs(volumes / expected_volume - 1)) > 1e-11:
                raise ValueError("Reference comparison volume does not match")
            component = np.sqrt(squared * gh.MULTIPLICITY)
            norm = np.linalg.norm(component, axis=1)
            delta = float(abs(norm[1] - norm[0]))
            if not np.all(np.isfinite(norm)) or delta > max(
                1e-11, 1e-5 * norm[1]
            ):
                raise ValueError(
                    f"Unconverged lifted-reference comparison at {time}:"
                    f" {delta}"
                )
            rows.append(
                dict(
                    time=time,
                    state_l2_difference=float(norm[1]),
                    metric_l2_difference=float(
                        np.linalg.norm(component[1, :10])
                    ),
                    pi_l2_difference=float(np.linalg.norm(component[1, 10:20])),
                    phi_l2_difference=float(np.linalg.norm(component[1, 20:])),
                    quadrature_change=delta,
                )
            )
            print(directory.name, reference.name, rows[-1], flush=True)
    if not rows or (
        selected_times is not None
        and set(row["time"] for row in rows) != set(selected_times)
    ):
        raise ValueError("Missing requested observations")
    result = dict(
        first=str(directory),
        second=str(reference),
        observations=rows,
        reference_kind=(
            "CartoonSphere1D, fixed mesh, spherical tensor reconstruction"
        ),
        quadrature_orders_angular_radial=orders,
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    analyzer = commands.add_parser("analyze")
    analyzer.add_argument("directory", type=Path)
    for command in ("compare-references", "compare-3d"):
        comparison = commands.add_parser(command)
        comparison.add_argument("first", type=Path)
        comparison.add_argument("second", type=Path)
        comparison.add_argument("--output", type=Path, required=True)
        if command == "compare-3d":
            comparison.add_argument("--times", type=float, nargs="+")
    args = parser.parse_args()
    if args.command == "analyze":
        print(json.dumps(reference_diagnostics(args.directory)))
    elif args.command == "compare-references":
        compare_references(args.first, args.second, args.output)
    else:
        compare_3d(args.first, args.second, args.output, args.times)


if __name__ == "__main__":
    main()
