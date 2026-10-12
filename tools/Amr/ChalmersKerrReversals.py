# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Replay immediate split/join reversals from post-AMR volume snapshots.

The replay uses the native projection and estimator bindings. Snapshots must
cover every AMR round. The joined parent is checked before and after filtering.
The evolved parent is reconstructed from its children after the next split,
provided their extents retain its information. This distinguishes an immediate
estimator inconsistency from evolution or mesh-balance enforcement.
"""

import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np

from spectre.DataStructures import DataVector
from spectre.Domain import ElementId
from spectre.NumericalAlgorithms.LinearOperators import chalmers_energy_estimate
from spectre.Spectral import (
    Basis,
    Mesh,
    Quadrature,
    SegmentSize,
    exponential_filter,
    projection_matrix_child_to_parent,
    projection_matrix_parent_to_child,
)
from spectre.tools.Amr import ChalmersCube as cube
from spectre.tools.Amr import ChalmersSphere as sphere


def key(name):
    eid = ElementId[3](name)
    return (
        eid.block_id,
        eid.grid_index,
        tuple((s.refinement_level, s.index) for s in eid.segment_ids),
    )


def contains(parent, child):
    return parent[:2] == child[:2] and all(
        lc >= lp and ic >> (lc - lp) == ip
        for (lp, ip), (lc, ic) in zip(parent[2], child[2])
    )


def topology(old, new):
    added, removed = set(new) - set(old), set(old) - set(new)
    result = []
    for action, parents, children in (
        ("split", removed, added),
        ("join", added, removed),
    ):
        for parent in parents:
            family = sorted(c for c in children if contains(parent, c))
            if not family:
                continue
            axes = tuple(
                lc - lp for (lp, _), (lc, _) in zip(parent[2], family[0][2])
            )
            if len(family) != 2 ** sum(axes) or any(
                d not in (0, 1) for d in axes
            ):
                raise ValueError("Snapshots do not isolate one AMR round")
            result.append((action, parent, family, axes))
    return result


def mesh(extents):
    return Mesh[len(extents)](extents, Basis.Legendre, Quadrature.GaussLobatto)


@lru_cache(maxsize=None)
def projection(source, target, half):
    size = {
        0: SegmentSize.Full,
        -1: SegmentSize.LowerHalf,
        1: SegmentSize.UpperHalf,
    }[half]
    return np.asarray(
        projection_matrix_child_to_parent(
            mesh([source]), mesh([target]), size, False
        )
    )


def project_parent(children, axis):
    extents = tuple(
        max(a, b) for a, b in zip(children[0].shape[1:], children[1].shape[1:])
    )
    parent = np.zeros((5, *extents))
    for half, child in zip((-1, 1), children):
        values = child
        for d, (source, target) in enumerate(zip(child.shape[1:], extents)):
            matrix = projection(source, target, half if d == axis else 0)
            values = np.moveaxis(
                np.tensordot(matrix, values, axes=(1, d + 1)), 0, d + 1
            )
        parent += values
    return parent


def recover_parent(parent, family, fields, extents, filtered_modes=0):
    """Undo a split, optionally zero-filling modes suppressed before AMR.

    A one-point p decrease can be inverted to filter/roundoff accuracy when
    the preceding RK filter suppressed the old highest mode by exp(-36).
    """
    if any(
        any(
            n < p - filtered_modes for n, p in zip(fields[c].shape[1:], extents)
        )
        for c in family
    ):
        return None, {"reason": "simultaneous_p_coarsening"}
    assumes_filtered_modes = any(
        any(n < p for n, p in zip(fields[c].shape[1:], extents)) for c in family
    )
    result = np.zeros((5, *extents))
    for child in family:
        values = fields[child]
        for d, (source, target) in enumerate(zip(values.shape[1:], extents)):
            lp, _ = parent[2][d]
            lc, ic = child[2][d]
            half = 0 if lp == lc else 1 if ic % 2 else -1
            matrix = projection(source, target, half)
            values = np.moveaxis(
                np.tensordot(matrix, values, axes=(1, d + 1)), 0, d + 1
            )
        result += values
    residual = 0.0
    for child in family:
        values = result
        for d, (source, target) in enumerate(
            zip(extents, fields[child].shape[1:])
        ):
            lp, _ = parent[2][d]
            lc, ic = child[2][d]
            size = (
                SegmentSize.Full
                if lp == lc
                else SegmentSize.UpperHalf if ic % 2 else SegmentSize.LowerHalf
            )
            matrix = (
                projection(source, target, 0)
                if size == SegmentSize.Full
                else np.asarray(
                    projection_matrix_parent_to_child(
                        mesh([source]), mesh([target]), size
                    )
                )
            )
            values = np.moveaxis(
                np.tensordot(matrix, values, axes=(1, d + 1)), 0, d + 1
            )
        residual = max(residual, float(np.max(abs(values - fields[child]))))
    if residual > 1e-11:
        return None, {
            "reason": "split_projection_not_invertible",
            "maximum_child_residual": residual,
        }
    return result, {
        "maximum_child_residual": residual,
        "uses_suppressed_mode_assumption": assumes_filtered_modes,
        "filtered_modes": filtered_modes if assumes_filtered_modes else 0,
    }


def filter_parent(fields, half_power):
    result = fields
    for d, n in enumerate(fields.shape[1:]):
        matrix = np.asarray(exponential_filter(mesh([n]), 36.0, half_power))
        result = np.moveaxis(
            np.tensordot(matrix, result, axes=(1, d + 1)), 0, d + 1
        )
    return result


def refinement(fields, tolerance, separate_parity=False):
    estimates = []
    flags = np.zeros(3, dtype=int)
    for field, values in zip(cube.FIELDS, fields):
        if np.sqrt(8) * np.max(abs(values)) <= tolerance:
            estimates.append({"field": field, "negligible": True})
            continue
        fits = chalmers_energy_estimate(
            DataVector(values.ravel(order="F")),
            mesh(list(values.shape)),
            tolerance,
            1,
            separate_parity,
        )
        estimates.append({"field": field, "directions": fits})
        for d, fit in enumerate(fits):
            if fit["error"] > tolerance:
                if not np.isfinite(fit["sigma"]):
                    raise ValueError("Undetermined replay fit")
                flags[d] = max(flags[d], 1 if fit["sigma"] > 0.5 else 2)
    return {
        "flags": ["none" if f == 0 else "p" if f == 1 else "h" for f in flags],
        "fits": estimates,
    }


def analyze(directory):
    metadata = json.loads((directory / "Run.json").read_text())
    if (
        metadata["returncode"]
        or not metadata["h_coarsening"]
        or metadata["include_highest"]
    ):
        raise ValueError("Expected a successful filter-aware hp-coarsening run")

    def replay_refinement(fields, tolerance):
        return refinement(
            fields, tolerance, metadata.get("separate_parity", False)
        )

    tolerance = metadata["tolerance"]
    interval = metadata["observation_interval"]
    if abs(interval - 0.4) > 1e-12:
        raise ValueError("Need snapshots at every 0.4M AMR round")
    records = []
    opposite_events_on_different_axes = 0
    with h5py.File(directory / "Volume0.h5") as handle:
        groups = sphere.groups(handle)
        times = [float(g.attrs["observation_value"]) for g in groups]
        if not np.allclose(np.diff(times), interval, rtol=0, atol=1e-12):
            raise ValueError("Missing snapshot")
        frames = [
            {key(n): n for n in bytes(g["grid_names"][:]).decode().split(":")}
            for g in groups
        ]
        meshes = [
            dict(
                zip(
                    frame,
                    np.asarray(group["total_extents"]).reshape(-1, 3),
                )
            )
            for frame, group in zip(frames, groups)
        ]
        events = [topology(a, b) for a, b in zip(frames, frames[1:])]
        join_costs = []
        for step, changes in enumerate(events):
            for action, parent, family, axes in changes:
                if action != "join":
                    continue
                before = sum(int(np.prod(meshes[step][c])) for c in family)
                after = int(np.prod(meshes[step + 1][parent]))
                join_costs.append(
                    {
                        "time": times[step + 1],
                        "parent": repr(parent),
                        "child_points": before,
                        "parent_points": after,
                        "points_removed": before - after,
                    }
                )

        @lru_cache(maxsize=3)
        def fields(index):
            return {key(n): u for n, _, _, _, u in cube.elements(groups[index])}

        for step in range(len(events) - 1):
            following = {
                (action, parent): (family, axes)
                for action, parent, family, axes in events[step + 1]
            }
            for action, parent, family, axes in events[step]:
                opposite = "split" if action == "join" else "join"
                if (opposite, parent) not in following:
                    continue
                next_family, next_axes = following[opposite, parent]
                # Require a reversal along at least one of the same axes.
                if not any(a and b for a, b in zip(axes, next_axes)):
                    opposite_events_on_different_axes += 1
                    continue
                if action == "join":
                    joined = fields(step + 1)[parent]
                    raw = replay_refinement(joined, tolerance)
                    filtered = replay_refinement(
                        filter_parent(joined, metadata["half_power"]), tolerance
                    )
                    evolved, recovery = recover_parent(
                        parent,
                        next_family,
                        fields(step + 2),
                        joined.shape[1:],
                        filtered_modes=1,
                    )
                    later = (
                        replay_refinement(evolved, tolerance)
                        if evolved is not None
                        else None
                    )
                    split_axes = [d for d, v in enumerate(next_axes) if v]
                    if any(
                        f != "none" for f in raw["flags"] + filtered["flags"]
                    ):
                        cause = "joined_parent_requests_refinement"
                    elif evolved is None:
                        cause = "ambiguous_split_reconstruction"
                    elif any(later["flags"][d] == "h" for d in split_axes):
                        cause = "evolution_requests_split"
                    else:
                        cause = "split_not_requested_by_local_criterion"
                    detail = {
                        "joined_parent": raw,
                        "filtered_parent": filtered,
                        "evolved_parent": later,
                        "reconstruction": recovery,
                    }
                else:
                    original, recovery = recover_parent(
                        parent,
                        family,
                        fields(step + 1),
                        fields(step)[parent].shape[1:],
                        filtered_modes=1,
                    )
                    before = (
                        replay_refinement(original, tolerance)
                        if original is not None
                        else None
                    )
                    joined = fields(step + 2)[parent]
                    after = replay_refinement(
                        filter_parent(joined, metadata["half_power"]), tolerance
                    )
                    cause = (
                        "ambiguous_split_reconstruction"
                        if before is None
                        else (
                            "local_split_then_join"
                            if any(
                                before["flags"][d] == "h"
                                for d, v in enumerate(axes)
                                if v
                            )
                            else "nonlocal_split_then_join"
                        )
                    )
                    detail = {
                        "original_parent": before,
                        "filtered_rejoined_parent": after,
                        "reconstruction": recovery,
                    }
                records.append(
                    {
                        "first_decision_time": times[step + 1],
                        "second_decision_time": times[step + 2],
                        "sequence": action + "-" + opposite,
                        "parent": repr(parent),
                        "first_axes": axes,
                        "second_axes": next_axes,
                        "classification": cause,
                        **detail,
                    }
                )
    result = {
        "run": directory.name,
        "snapshot_order": "post-AMR, before the next evolution step",
        "tolerance": tolerance,
        "counts": dict(Counter(r["classification"] for r in records)),
        "same_axis_immediate_reversals": len(records),
        "opposite_events_on_different_axes": opposite_events_on_different_axes,
        "join_costs": {
            "join_events": len(join_costs),
            "total_points_removed_at_joins": sum(
                r["points_removed"] for r in join_costs
            ),
            "joins_without_point_reduction": [
                r for r in join_costs if r["points_removed"] <= 0
            ],
            "events": join_costs,
        },
        "reversals": records,
        "interpretation": (
            "A split absent from replayed local requests is consistent with"
            " neighbor balance enforcement; neighbor messages are not recorded."
            " A simultaneous one-point p decrease is reconstructed assuming the"
            " old highest mode was suppressed on the preceding RK substep; use"
            " of this assumption is recorded per event. Larger decreases remain"
            " ambiguous. Reconstructed fits near the coefficient cutoff can"
            " themselves change under roundoff; these are replayed requests,"
            " not recorded native criterion flags."
        ),
    }
    (directory / "Reversals.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(directory.name, result["counts"], flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    analyze(parser.parse_args().directory)
