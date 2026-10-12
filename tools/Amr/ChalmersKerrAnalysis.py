# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Independent proper-volume norms and reference comparisons for ChalmersKerr.

Uses analytic Kerr-Schild spatial metrics, map Jacobians, and two independent
Gauss quadratures. Waveforms are complex coordinate-sphere Y_lm coefficients
at r=12, with the Condon-Shortley convention and measure dOmega.
"""

import argparse
import json
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legvander

from spectre.DataStructures import DataVector
from spectre.DataStructures.Tensor import Frame, tnsr
from spectre.Domain import (
    ElementId,
    ElementMap,
    block_logical_coordinates,
    element_logical_coordinates,
)
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv

try:
    from scipy.special import sph_harm_y
except ImportError:
    from scipy.special import sph_harm

    def sph_harm_y(ell, m, theta, phi):
        return sph_harm(m, ell, phi, theta)


ORDERS = (24, 32)


def harmonic(ell, m, xyz):
    """Return normalized scalar Y_lm, including its complex phase."""
    if ell < 0 or abs(m) > ell:
        raise ValueError("Spherical-harmonic mode requires ell >= |m|")
    radius = np.linalg.norm(xyz, axis=0)
    if m == 0:
        # Preserve the real-valued axisymmetric diagnostics.
        return np.sqrt((2 * ell + 1) / (4 * np.pi)) * legval(
            xyz[2] / radius, [0] * ell + [1]
        )
    return sph_harm_y(
        ell,
        m,
        np.arctan2(np.hypot(xyz[0], xyz[1]), xyz[2]),
        np.arctan2(xyz[1], xyz[0]),
    )


def waveform_modes(metadata):
    """Select signal, conjugacy, and symmetry-forbidden control modes."""
    m = metadata.get("m", 0)
    if m == 0:
        return ((0, 0), (2, 0), (4, 0))
    ell = metadata["ell"]
    return (
        (ell, m),
        (ell + 2, m),
        (ell + 4, m),
        (ell, -m),
        (0, 0),
        (2, 0),
        (2, 2 if abs(m) != 2 else 1),
        (ell + 1, m),
        # Cubical wedges can couple azimuthal modes differing by four.
        (4, 3) if abs(m) == 1 else (6, 6),
    )


def waveform_columns(modes, values):
    result = {}
    for (ell, m), value in zip(modes, values):
        key = f"psi_l{ell}_m{m}"
        if m == 0:
            result[key] = float(np.real(value))
        else:
            result[key + "_real"] = float(np.real(value))
            result[key + "_imag"] = float(np.imag(value))
    return result


class PatchWaveforms:
    """Integrate coordinate-sphere harmonics on individual DG patches.

    This extractor is restricted to static spherical wedges with linear radial
    maps and logical axis 2 radial, as used by ChalmersKerr. Each intersecting
    element contributes its own trace, avoiding quadrature across DG seams.
    The measure is dOmega, not the Kerr proper-area element.
    """

    def __init__(self, domain, radius=12.0, modes=((0, 0), (2, 0), (4, 0))):
        if domain.is_time_dependent():
            raise ValueError("Patch waveforms require a static domain")
        self.domain = domain
        self.radius = radius
        self.modes = tuple(modes)
        if not self.modes or any(ell < abs(m) for ell, m in self.modes):
            raise ValueError("Invalid extraction modes")
        self.dtype = complex if any(m != 0 for _, m in self.modes) else float

    @lru_cache(maxsize=2048)
    def operator(self, name, extents, order):
        element_map = ElementMap(ElementId[3](name), self.domain)
        ends, _, _ = sphere.geometry(element_map, [[0.0], [0.0], [-1.0, 1.0]])
        radii = np.linalg.norm(ends.reshape(3, 2), axis=0)
        if not radii[0] < radii[1]:
            raise ValueError("Expected outward radial logical axis 2")
        if not radii[0] <= self.radius < radii[1]:
            return None
        zeta = (2 * self.radius - radii.sum()) / (radii[1] - radii[0])
        xi, weights = leggauss(order)
        xyz, det, inverse = sphere.geometry(element_map, [xi, xi, [zeta]])
        if not np.allclose(
            np.linalg.norm(xyz, axis=0), self.radius, rtol=0, atol=2e-11
        ):
            raise ValueError(
                "Extraction requires a spherical linear radial map"
            )
        # |dx/dxi cross dx/deta| = det(J) |grad(zeta)|.
        measure = (
            weights[:, None]
            * weights[None, :]
            * (det * np.linalg.norm(inverse[2], axis=0))[:, :, 0]
            / self.radius**2
        )
        angular = [legvander(xi, n - 1) for n in extents[:2]]
        matrices = np.array(
            [
                np.einsum(
                    "ia,jb,ij->ab",
                    *angular,
                    measure * np.conjugate(harmonic(ell, m, xyz[:, :, :, 0])),
                )
                for ell, m in self.modes
            ]
        )
        radial = legvander([zeta], extents[2] - 1)[0]
        return matrices, radial, float(measure.sum())

    def __call__(self, fields, order=32):
        waveform = np.zeros(len(self.modes), dtype=self.dtype)
        area = 0.0
        for name, modes in fields.items():
            psi = modes[0]
            operator = self.operator(name, psi.shape, order)
            if operator is None:
                continue
            angular, radial, solid_angle = operator
            waveform += np.einsum("lab,abc,c->l", angular, psi, radial)
            area += solid_angle
        if abs(area - 4 * np.pi) > 1e-10:
            raise ValueError(
                f"Incomplete extraction sphere: solid angle {area}"
            )
        return waveform


def metric(xyz, spin):
    """Return gamma^ij and sqrt(det(gamma)) for centered unit-mass Kerr."""
    x, y, z = xyz
    rho2 = x * x + y * y + z * z
    r2 = 0.5 * (
        rho2 - spin**2 + np.sqrt((rho2 - spin**2) ** 2 + 4 * spin**2 * z * z)
    )
    r = np.sqrt(r2)
    h = r**3 / (r2 * r2 + spin**2 * z * z)
    ell = np.array(
        [
            (r * x + spin * y) / (r2 + spin**2),
            (r * y - spin * x) / (r2 + spin**2),
            z / r,
        ]
    )
    inverse = (
        np.eye(3).reshape((3, 3) + (1,) * x.ndim)
        - (2 * h / (1 + 2 * h)) * ell[:, None] * ell[None, :]
    )
    return inverse, np.sqrt(1 + 2 * h)


def state_squared(fields, inverse):
    return (
        fields[0] ** 2
        + fields[1] ** 2
        + np.einsum("ij...,i...,j...->...", inverse, fields[2:], fields[2:])
    )


def polynomial_fields(group):
    return {
        name: [to_modal(v) for v in fields]
        for name, _, _, _, fields in cube.elements(group)
    }


def diagnose(directory, *, patch_waveforms=False):
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata["returncode"]:
        raise ValueError("Evolution failed")
    domain = sphere.read_domain(directory)
    modes_to_extract = waveform_modes(metadata)
    patch_extractor = (
        PatchWaveforms(domain, modes=modes_to_extract)
        if patch_waveforms
        else None
    )
    mu, w = leggauss(24)
    phi = 2 * np.pi * (np.arange(48) + 0.5) / 48
    mu, phi = np.meshgrid(mu, phi, indexing="ij")
    surface = np.array(
        [
            12 * np.sqrt(1 - mu**2) * np.cos(phi),
            12 * np.sqrt(1 - mu**2) * np.sin(phi),
            12 * mu,
        ]
    ).reshape(3, -1)
    surface_weights = np.repeat(w, 48) * 2 * np.pi / 48
    block_coords = block_logical_coordinates(
        domain,
        tnsr.I[DataVector, 3, Frame.Inertial](surface),
        time=0.0,
        functions_of_time={},
    )
    rows = []
    max_quad = 0.0
    with h5py.File(directory / "Volume0.h5") as handle:
        for group in sphere.groups(handle):
            t = float(group.attrs["observation_value"])
            fields = polynomial_fields(group)
            sums = {q: np.zeros(4) for q in ORDERS}
            for name, modes in fields.items():
                element_map = ElementMap(ElementId[3](name), domain)
                for q in ORDERS:
                    xi, _, weights = cube.quadrature(q)
                    xyz, det, jinv = sphere.geometry(element_map, [xi] * 3)
                    inverse, sqrtg = metric(xyz, metadata["spin"])
                    measure = weights * det * sqrtg
                    u = np.array([cube.evaluate(c, q) for c in modes])
                    c1 = sphere.physical_gradient(modes[0], jinv, q) - u[2:]
                    dphi = np.array(
                        [
                            sphere.physical_gradient(c, jinv, q)
                            for c in modes[2:]
                        ]
                    )
                    c2 = dphi.swapaxes(0, 1) - dphi
                    sums[q] += np.array(
                        [
                            np.sum(measure * state_squared(u, inverse)),
                            np.sum(
                                measure
                                * np.einsum(
                                    "ij...,i...,j...->...", inverse, c1, c1
                                )
                            ),
                            np.sum(
                                measure
                                * 0.5
                                * np.einsum(
                                    "ik...,jl...,ij...,kl...->...",
                                    inverse,
                                    inverse,
                                    c2,
                                    c2,
                                    optimize=True,
                                )
                            ),
                            np.sum(measure),
                        ]
                    )
            low, high = [np.sqrt(np.maximum(sums[q][:3], 0)) for q in ORDERS]
            change = float(np.max(abs(high - low)))
            if not np.allclose(low, high, atol=1.0e-10, rtol=1.0e-5):
                raise RuntimeError(f"Quadrature not converged: {change}")
            max_quad = max(max_quad, change)
            row = {
                "time": t,
                "state_l2": high[0],
                "constraint_l2": high[1],
                "curl_constraint_l2": high[2],
                "proper_volume": sums[ORDERS[-1]][3],
                "quadrature_change": change,
                "elements": len(fields),
                "points": sum(c[0].size for c in fields.values()),
            }
            if patch_extractor is not None:
                coarse, fine = [patch_extractor(fields, q) for q in ORDERS]
                change = float(np.max(abs(coarse - fine)))
                if change > 1e-11:
                    raise RuntimeError(
                        f"Waveform quadrature not converged: {change}"
                    )
                row["waveform_quadrature_change"] = change
                row.update(waveform_columns(modes_to_extract, fine))
            else:
                logical = element_logical_coordinates(
                    [ElementId[3](n) for n in fields], block_coords
                )
                surface_psi = np.full(surface.shape[1], np.nan)
                for eid, info in logical.items():
                    c = fields[str(eid)][0]
                    coords = np.asarray(info.element_logical_coords)
                    matrices = [
                        legvander(coords[d], c.shape[d] - 1) for d in range(3)
                    ]
                    surface_psi[info.offsets] = np.einsum(
                        "pi,pj,pk,ijk->p", *matrices, c
                    )
                if not np.all(np.isfinite(surface_psi)):
                    raise RuntimeError("Incomplete extraction sphere")
                values = [
                    np.sum(
                        surface_weights
                        * surface_psi
                        * np.conjugate(harmonic(ell, m, surface))
                    )
                    for ell, m in modes_to_extract
                ]
                row.update(waveform_columns(modes_to_extract, values))
            rows.append(row)
            print(directory.name, t, "C=", high[1], flush=True)
    write_csv(directory / "PhysicalHistory.csv", rows)
    result = {
        "run": directory.name,
        "maximum_quadrature_change": max_quad,
        "initial": rows[0],
        "final": rows[-1],
        "maximum_constraint": max(r["constraint_l2"] for r in rows),
        "maximum_curl_constraint": max(r["curl_constraint_l2"] for r in rows),
        "waveform_modes": modes_to_extract,
        "waveform_definition": "integral Psi conjugate(Y_lm) dOmega at r=12",
    }
    if patch_waveforms:
        result["waveform_method"] = (
            "Gauss quadrature on individual DG sphere patches"
        )
        result["maximum_waveform_quadrature_change"] = max(
            r["waveform_quadrature_change"] for r in rows
        )
    (directory / "PhysicalAnalysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def difference(a, b, domain, spin, order):
    xi, _, weights = cube.quadrature(order)
    squared = np.zeros(2)
    covered = {}
    for na, ma in a.items():
        ia = ElementId[3](na)
        ba = [s.endpoints for s in ia.segment_ids]
        for nb, mb in b.items():
            ib = ElementId[3](nb)
            if ia.block_id != ib.block_id:
                continue
            bb = [s.endpoints for s in ib.segment_ids]
            overlap = [
                (max(x[0], y[0]), min(x[1], y[1])) for x, y in zip(ba, bb)
            ]
            if any(hi <= lo for lo, hi in overlap):
                continue
            scale = np.prod([(hi - lo) / 2 for lo, hi in overlap])
            covered[ia.block_id] = covered.get(ia.block_id, 0) + 8 * scale
            axes = [(lo + hi) / 2 + (hi - lo) / 2 * xi for lo, hi in overlap]
            xyz, det, _ = sphere.geometry(
                ElementMap(ElementId[3](ia.block_id), domain), axes
            )
            inverse, sqrtg = metric(xyz, spin)

            def restrict(modes, bounds):
                result = modes
                for d, ((lo, hi), nodes) in enumerate(zip(bounds, axes)):
                    matrix = legvander(
                        (2 * nodes - lo - hi) / (hi - lo), modes.shape[d] - 1
                    )
                    result = np.moveaxis(
                        np.tensordot(matrix, result, axes=(1, d)), 0, d
                    )
                return result

            delta = np.array(
                [restrict(x, ba) - restrict(y, bb) for x, y in zip(ma, mb)]
            )
            measure = scale * weights * det * sqrtg
            squared += np.array(
                [
                    np.sum(measure * state_squared(delta, inverse)),
                    np.sum(measure * delta[0] ** 2),
                ]
            )
    blocks_a = {ElementId[3](n).block_id for n in a}
    blocks_b = {ElementId[3](n).block_id for n in b}
    if (
        blocks_a != blocks_b
        or set(covered) != blocks_a
        or any(abs(v - 8) > 1.0e-12 for v in covered.values())
    ):
        raise ValueError("Domain partition incomplete")
    return np.sqrt(squared)


def compare(first, second, output, *, time=None):
    domain = sphere.read_domain(first)
    if domain != sphere.read_domain(second):
        raise ValueError("Domains differ")
    metadata = [
        json.loads((d / "Run.json").read_text()) for d in [first, second]
    ]
    if any(metadata[0][k] != metadata[1][k] for k in ["spin", "ell"]) or (
        metadata[0].get("m", 0) != metadata[1].get("m", 0)
    ):
        raise ValueError("Different physical problems")
    with (
        h5py.File(first / "Volume0.h5") as fa,
        h5py.File(second / "Volume0.h5") as fb,
    ):

        def snapshot(handle):
            snapshots = sphere.groups(handle)
            if time is None:
                return snapshots[-1]
            matching = [
                g
                for g in snapshots
                if abs(g.attrs["observation_value"] - time) < 1.0e-10
            ]
            if len(matching) != 1:
                raise ValueError(f"Expected one snapshot at time {time}")
            return matching[0]

        a = snapshot(fa)
        b = snapshot(fb)
        if (
            abs(a.attrs["observation_value"] - b.attrs["observation_value"])
            > 1.0e-10
        ):
            raise ValueError("Times differ")
        comparison_time = float(a.attrs["observation_value"])
        a, b = polynomial_fields(a), polynomial_fields(b)
        low, high = [
            difference(a, b, domain, metadata[0]["spin"], q) for q in ORDERS
        ]
    change = float(max(abs(high - low)))
    if change > max(1.0e-11, 1.0e-5 * max(high)):
        raise RuntimeError("Difference quadrature not converged")
    result = {
        "first": str(first),
        "second": str(second),
        "time": comparison_time,
        "state_l2_difference": float(high[0]),
        "psi_l2_difference": float(high[1]),
        "quadrature_change": change,
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--compare", type=Path, metavar="OUTPUT_JSON")
    parser.add_argument("--time", type=float)
    parser.add_argument("--patch-waveforms", action="store_true")
    args = parser.parse_args()
    if args.compare:
        if len(args.runs) != 2:
            parser.error("--compare requires exactly two run directories")
        print(compare(*args.runs, args.compare, time=args.time))
        return
    for directory in args.runs:
        diagnose(directory, patch_waveforms=args.patch_waveforms)


if __name__ == "__main__":
    main()
