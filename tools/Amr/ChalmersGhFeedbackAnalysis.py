# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Directional GH constraints and spectral diagnostics for feedback trials.

The directional constraints match the controller's nodal RMS with full tensor
multiplicities. The integrated radial/tangential diagnostics use physical
coordinates and two coordinate-volume quadratures. Spectral diagnostics are
passive: the ordinary SpEC-style multi-fit slope is not symmetry aware, and
its raw value must not replace the parity-aware Chalmers decision.
Native logs describe the pre-AMR state. Volume observations at the same time
follow AMR, so changed elements need not have the same nodal constraint norm.
"""

import argparse
import gzip
import json
import re
from pathlib import Path

import h5py
import numpy as np
from numpy.polynomial.legendre import legder

from spectre.DataStructures import DataVector
from spectre.Domain import ElementId, ElementMap
from spectre.NumericalAlgorithms.LinearOperators import (
    chalmers_energy_estimate,
    convergence_rate_and_number_of_pile_up_modes,
    power_monitors,
)
from spectre.Spectral import Basis, Mesh3D, Quadrature
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSingleHole as gh
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import legendre_grid, to_modal
from spectre.tools.Amr.ChalmersEvolution import write_csv
from spectre.tools.Amr.ChalmersSphericalReference import tensor_evaluate


def directional_constraints(fields, modes, inverse, axes):
    logical = np.array(
        [
            tensor_evaluate(legder(modes[:10], axis=d + 1), axes)
            for d in range(3)
        ]
    )
    physical = np.einsum("di...,da...->ia...", inverse, logical)
    constraint = physical - fields[20:].reshape((3, 10, *fields.shape[1:]))
    jacobian = np.moveaxis(
        np.linalg.inv(np.moveaxis(inverse, (0, 1), (-2, -1))),
        (-2, -1),
        (0, 1),
    )
    projected = np.einsum("id...,ia...->da...", jacobian, constraint)
    squared = np.sum(
        projected**2 * gh.MULTIPLICITY[:10][None, :, None, None, None], axis=1
    ) / np.sum(jacobian**2, axis=0)
    return np.sqrt(np.mean(squared, axis=(1, 2, 3)))


def relative_tail_error(power):
    """SpECTRE's weighted tail diagnostic applied to unfiltered modes.

    Match its below-floor convention. Explicitly handle zero spectra and
    absent leading coefficients instead of taking log(0). These are passive
    diagnostics, not a replacement for the Chalmers absolute error estimate.
    """
    power = np.asarray(power)
    scale = np.max(power)
    if scale == 0:
        return 0.0
    cutoff = 100 * np.finfo(float).eps * scale
    if np.all(power[-2:] < cutoff):
        return cutoff * 0.01
    leading = max(power[:2])
    if leading <= cutoff:
        return None
    use = power >= cutoff
    weights = np.exp(-((len(power) - 1 - np.arange(len(power)) - 0.5) ** 2))
    return float(
        10 ** np.average(np.log10(power[use]), weights=weights[use]) / leading
    )


def finite(value):
    return float(value) if np.isfinite(value) else None


def spectral_diagnostics(fields, extents, target, factors):
    mesh = Mesh3D(extents, Basis.Legendre, Quadrature.GaussLobatto)
    records = []
    for field, values in zip(gh.FIELDS, fields):
        data = DataVector(values.ravel(order="F"))
        amplitude = float(np.max(abs(values)))
        negligible = bool(np.sqrt(8) * amplitude <= target * min(factors))
        estimates = (
            [dict(error=0, sigma=np.nan, max_log_fit_rms=0)] * 3
            if negligible
            else chalmers_energy_estimate(
                data, mesh, target * min(factors), 1, True
            )
        )
        for d, power in enumerate(power_monitors(data, mesh)):
            power = np.asarray(power)
            diagnostic = convergence_rate_and_number_of_pile_up_modes(
                DataVector(power), 1
            )
            relative = relative_tail_error(power[:-1])
            records.append(
                dict(
                    field=field,
                    direction=d,
                    amplitude=amplitude,
                    negligible=negligible,
                    power=power.tolist(),
                    chalmers_error=finite(estimates[d]["error"]),
                    chalmers_sigma=finite(estimates[d]["sigma"]),
                    chalmers_fit_rms=finite(estimates[d]["max_log_fit_rms"]),
                    effective_target=float(target * factors[d]),
                    spec_relative_tail=relative,
                    spec_scaled_tail=(
                        None if relative is None else relative * amplitude
                    ),
                    spec_log10_slope=finite(diagnostic["convergence_rate"]),
                    spec_pile_up=finite(diagnostic["number_of_pile_up_modes"]),
                )
            )
    return records


def analyze(directory, output=None, spectra=True, integrated=True):
    metadata = json.loads((directory / "Run.json").read_text())
    if metadata["returncode"]:
        raise ValueError("Need a completed evolution")
    output = directory if output is None else output
    output.mkdir(parents=True, exist_ok=True)
    target = metadata["tolerance"]
    constraint_target = metadata.get("constraint_target") or 1.0e-4
    feedback = metadata.get("constraint_feedback", "off")
    minimum_factor = (
        metadata.get("minimum_target_factor", 0.1)
        if feedback == "tighten"
        else 1.0
    )
    domain = sphere.read_domain(directory)
    rows, histories = [], []
    spectra_path = output / "FeedbackSpectra.json.gz"
    with (
        gzip.open(spectra_path, "wt") as spectral_file,
        h5py.File(directory / "Volume0.h5") as handle,
    ):
        for group in sphere.groups(handle):
            time = float(group.attrs["observation_value"])
            totals = np.zeros((2, 2))
            for name, ext, xyz, fields in gh.elements(group):
                modes = np.array([to_modal(v) for v in fields])
                eid = ElementId[3](name)
                element_map = ElementMap(eid, domain)
                axes = [legendre_grid(n)[0] for n in ext]
                _, _, inverse = sphere.geometry(element_map, axes)
                constraints = directional_constraints(
                    fields, modes, inverse, axes
                )
                if not np.all(np.isfinite(constraints)):
                    raise ValueError("Nonfinite constraints")
                factors = np.maximum(
                    minimum_factor,
                    np.minimum(
                        1.0, constraint_target / np.maximum(constraints, 1e-300)
                    ),
                )
                for d in range(3):
                    rows.append(
                        dict(
                            time=time,
                            element=name,
                            block=eid.block_id,
                            direction=d,
                            extent=int(ext[d]),
                            constraint=float(constraints[d]),
                            factor=float(factors[d]),
                        )
                    )
                if spectra:
                    for record in spectral_diagnostics(
                        fields, ext, target, factors
                    ):
                        record.update(time=time, element=name)
                        spectral_file.write(
                            json.dumps(record, allow_nan=False) + "\n"
                        )
                # The outer shell r=8..20 has blocks 18..23 in this study.
                if not integrated or eid.block_id < 18:
                    continue
                for q, order in enumerate((20, 28)):
                    xi, _, weights = cube.quadrature(order)
                    mapped, det, inv = sphere.geometry(element_map, [xi] * 3)
                    radial = mapped / np.linalg.norm(mapped, axis=0)
                    phi = tensor_evaluate(modes[20:], [xi] * 3).reshape(
                        3, 10, order, order, order
                    )
                    logical = np.array(
                        [
                            tensor_evaluate(
                                legder(modes[:10], axis=d + 1), [xi] * 3
                            )
                            for d in range(3)
                        ]
                    )
                    c = np.einsum("di...,da...->ia...", inv, logical) - phi
                    cr = np.sum(c * radial[:, None, ...], axis=0)
                    ct = c - cr[None, ...] * radial[:, None, ...]
                    measure = (
                        weights * det * gh.MULTIPLICITY[:10, None, None, None]
                    )
                    totals[q, 0] += np.sum(measure * cr**2)
                    totals[q, 1] += np.sum(measure[None, ...] * ct**2)
            norms = np.sqrt(totals)
            difference = float(np.max(abs(norms[1] - norms[0])))
            if difference > max(1e-11, 1e-5 * np.max(norms[1])):
                raise ValueError("Unconverged radial constraint quadrature")
            histories.append(
                dict(
                    time=time,
                    radial=float(norms[1, 0]),
                    tangential=float(norms[1, 1]),
                    quadrature_change=difference,
                )
            )
            print(
                directory.name,
                time,
                "feedback diagnostics complete",
                flush=True,
            )
    write_csv(output / "FeedbackByElement.csv", rows)
    write_csv(output / "OuterConstraintDirections.csv", histories)
    native = []
    pattern = re.compile(
        r"ChalmersFeedback time=(\S+) element=(\S+) constraints=\(([^)]*)\) "
        r"factors=\(([^)]*)\) before=\(([^)]*)\) after=\(([^)]*)\)"
    )
    for match in pattern.finditer((directory / "Run.log").read_text()):
        t, element, constraints, factors, before, after = match.groups()
        native.append(
            dict(
                time=float(t),
                element=element,
                constraints=[float(x) for x in constraints.split(",")],
                factors=[float(x) for x in factors.split(",")],
                before=before.split(","),
                after=after.split(","),
            )
        )
    result = dict(
        run=str(directory),
        constraint_target=constraint_target,
        mode=feedback,
        minimum_target_factor=minimum_factor,
        native_evaluations=native,
        outer_constraints=histories if integrated else [],
        passive_spectral_diagnostics=spectra,
        native_p_vetoes=sum(
            a.strip() == "DecreaseResolution" and b.strip() == "DoNothing"
            for row in native
            for a, b in zip(row["before"], row["after"])
        ),
        native_join_vetoes=sum(
            a.strip() == "Join" and b.strip() == "DoNothing"
            for row in native
            for a, b in zip(row["before"], row["after"])
        ),
    )
    (output / "FeedbackAnalysis.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-spectra", action="store_true")
    parser.add_argument("--no-integrated", action="store_true")
    args = parser.parse_args()
    analyze(
        args.directory, args.output, not args.no_spectra, not args.no_integrated
    )


if __name__ == "__main__":
    main()
