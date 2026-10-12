# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Measure actual GH error and constraints after candidate p-coarsening.

Uses the saved 20M gauge-relaxation problem, its native projection and filter,
and both converged Cartoon references. No spectral error budget is used to
accept/reject the measured result. Norms use coordinate volume and Euclidean
tensor components with symmetric multiplicities, as in the preceding study.
"""

import argparse
import concurrent.futures
import json
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import legder, leggauss

from spectre.Domain import ElementId, ElementMap
from spectre.tools.Amr import ChalmersSingleHole as gh
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr import ChalmersSphericalReference as reference
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersKerrReversals import filter_parent, projection

VARIANTS = {
    "original": (),
    "refilter": (),
    "drop0": (0,),
    "drop1": (1,),
    "drop01": (0, 1),
}
ORDERS = ((20, 40), (28, 48))


def candidate(fields, directions):
    result = fields.copy()
    for d in directions:
        n = result.shape[d + 1]
        result = np.moveaxis(
            np.tensordot(projection(n, n - 1, 0), result, axes=(1, d + 1)),
            0,
            d + 1,
        )
    return filter_parent(result, 256)


def metric_tensor(fields):
    result = np.empty((4, 4, *fields.shape[1:]))
    for value, (a, b) in zip(fields[:10], gh.PAIRS):
        result[a, b] = result[b, a] = value
    return result


def gauge_constraint(fields, xyz):
    """H_a + Gamma_a with DH amplitudes (1,0,1), exponents (2,2,2)."""
    g = metric_tensor(fields)
    pi = metric_tensor(fields[10:20])
    phi = np.array(
        [metric_tensor(fields[20 + 10 * i : 30 + 10 * i]) for i in range(3)]
    )
    spatial = np.moveaxis(g[1:, 1:], (0, 1), (-2, -1))
    inv_spatial = np.moveaxis(np.linalg.inv(spatial), (-2, -1), (0, 1))
    beta = np.einsum("ij...,j...->i...", inv_spatial, g[0, 1:])
    alpha = np.sqrt(-g[0, 0] + np.sum(beta * g[0, 1:], axis=0))
    logfac = np.log(np.sqrt(np.linalg.det(spatial)) / alpha)
    mu = np.exp(-np.sum(xyz**2, axis=0) / 17.0152695482514**2) * logfac**2
    H = -mu / alpha * np.einsum("ai...,i...->a...", g[:, 1:], beta)
    H[0] -= alpha * mu * logfac
    inv = np.moveaxis(
        np.linalg.inv(np.moveaxis(g, (0, 1), (-2, -1))), (-2, -1), (0, 1)
    )
    d4g = np.concatenate(
        ((-alpha * pi + np.einsum("i...,iab...->ab...", beta, phi))[None], phi)
    )
    gamma = np.einsum("bc...,bca...->a...", inv, d4g) - 0.5 * np.einsum(
        "bc...,abc...->a...", inv, d4g
    )
    return H + gamma, H


def integrate_constraints(modes, element_map, order):
    x, w = leggauss(order)
    axes = [x] * 3
    xyz, det, inv = sphere.geometry(element_map, axes)
    values = reference.tensor_evaluate(modes, axes)
    logical = np.array(
        [
            reference.tensor_evaluate(legder(modes[:10], axis=d + 1), axes)
            for d in range(3)
        ]
    )
    C = np.einsum("di...,da...->ia...", inv, logical) - values[20:].reshape(
        (3, 10, order, order, order)
    )
    c_sq = np.sum(
        C**2 * gh.MULTIPLICITY[:10][None, :, None, None, None], axis=(0, 1)
    )
    G, _ = gauge_constraint(values, xyz)
    g_sq = np.sum(G**2, axis=0)
    measure = det * w[:, None, None] * w[None, :, None] * w[None, None, :]
    return np.array([np.sum(measure * c_sq), np.sum(measure * g_sq)]), np.sqrt(
        [np.max(c_sq), np.max(g_sq)]
    )


def integrate_errors(mode_sets, element_map, cells, orders):
    """Inner shell only: split radial integration at reference cell faces."""
    na, nr = orders
    xa, wa = leggauss(na)
    xr, wr = leggauss(nr)
    totals = {name: np.zeros(3) for name in mode_sets}
    for cell in cells:
        lo, hi = cell["logical_bounds"]
        axes = [xa, xa, (lo + hi) / 2 + (hi - lo) * xr / 2]
        xyz, det, _ = sphere.geometry(element_map, axes)
        radii = np.linalg.norm(xyz, axis=0)
        radial = reference.evaluate(cell, radii[0, 0, :])
        radial = np.broadcast_to(radial[:, None, None, :], (14, *radii.shape))
        expected = reference.reconstruct(radial, xyz / radii)
        measure = (
            (hi - lo)
            / 2
            * det
            * wa[:, None, None]
            * wa[None, :, None]
            * wr[None, None, :]
        )
        for name, modes in mode_sets.items():
            delta = reference.tensor_evaluate(modes, axes) - expected
            squared = (
                np.sum(measure * delta**2, axis=(1, 2, 3)) * gh.MULTIPLICITY
            )
            totals[name] += [
                np.sum(squared[:10]),
                np.sum(squared[10:20]),
                np.sum(squared[20:]),
            ]
    return totals


def snapshot(run, references, output, time):
    output.mkdir(parents=True, exist_ok=True)
    domain = sphere.read_domain(run)
    cells = [
        [c for c in reference.read(ref, time) if c["block"] == 0]
        for ref in references
    ]
    rows = []
    with h5py.File(run / "Volume0.h5") as f:
        group = next(
            g
            for g in sphere.groups(f)
            if float(g.attrs["observation_value"]) == time
        )
        offset = 0
        gauge_h_difference = 0.0
        for name, ext, xyz, fields in gh.elements(group):
            size = int(np.prod(ext))
            if ElementId[3](name).block_id >= 6:
                offset += size
                continue
            assert all(
                s.refinement_level == 0 for s in ElementId[3](name).segment_ids
            )
            _, H = gauge_constraint(fields, xyz)
            saved_H = np.array(
                [
                    group["GaugeH_" + a][offset : offset + size].reshape(
                        ext, order="F"
                    )
                    for a in ("t", "x", "y", "z")
                ]
            )
            gauge_h_difference = max(
                gauge_h_difference, float(np.max(abs(H - saved_H)))
            )
            offset += size
            element_map = ElementMap(ElementId[3](name), domain)
            variants = {
                key: (
                    fields
                    if key == "original"
                    else candidate(fields, directions)
                )
                for key, directions in VARIANTS.items()
            }
            modes = {
                key: np.array([to_modal(v) for v in vals])
                for key, vals in variants.items()
            }
            errors = [
                [
                    integrate_errors(modes, element_map, ref, order)
                    for order in ORDERS
                ]
                for ref in cells
            ]
            constraints = {
                key: [
                    integrate_constraints(mode, element_map, order)
                    for order in (20, 28)
                ]
                for key, mode in modes.items()
            }
            for key, values in variants.items():
                cs = [np.sqrt(c[0]) for c in constraints[key]]
                es = [[np.sqrt(order[key]) for order in ref] for ref in errors]
                for lo, hi in [cs, *es]:
                    np.testing.assert_allclose(lo, hi, atol=1e-10, rtol=2e-5)
                rows.append(
                    dict(
                        time=time,
                        element=name,
                        variant=key,
                        extents=list(values.shape[1:]),
                        points=int(np.prod(values.shape[1:])),
                        error_by_reference=[e[1].tolist() for e in es],
                        error_quadrature_changes=[
                            abs(e[1] - e[0]).tolist() for e in es
                        ],
                        constraints=cs[1].tolist(),
                        constraint_quadrature_change=abs(
                            cs[1] - cs[0]
                        ).tolist(),
                        constraint_maxima=constraints[key][1][1].tolist(),
                    )
                )
            print(time, name, "measured", flush=True)
    assert gauge_h_difference < 1e-12, gauge_h_difference
    original = {r["element"]: r for r in rows if r["variant"] == "original"}
    for row in rows:
        before = original[row["element"]]
        row["state_error_change_percent"] = 100 * (
            np.linalg.norm(row["error_by_reference"][1])
            / np.linalg.norm(before["error_by_reference"][1])
            - 1
        )
        row["constraint_change_percent"] = (
            100 * (np.array(row["constraints"]) / before["constraints"] - 1)
        ).tolist()
    result = dict(
        time=time,
        run=str(run.resolve()),
        reference_order=["coarse", "fine"],
        error_components=["metric", "pi", "phi"],
        constraint_components=["three_index", "gauge"],
        gauge_h_maximum_difference=gauge_h_difference,
        rows=rows,
    )
    (output / f"Snapshot{time:g}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--references", nargs=2, type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--times", nargs="+", type=float, default=[15, 20])
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers
    ) as pool:
        jobs = [
            pool.submit(snapshot, args.run, args.references, args.output, t)
            for t in args.times
        ]
        for job in jobs:
            result = job.result()
            print(result["time"], "complete", flush=True)


if __name__ == "__main__":
    main()
