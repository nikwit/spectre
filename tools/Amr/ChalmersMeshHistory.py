# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Count topology/degree changes with one AMR check between mesh snapshots.

The lightweight Meshes.vol stream must retain element IDs and extents at every
check. Counts describe realized mesh changes, including balance enforcement;
they do not count criterion requests overridden by other policies.
"""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

from spectre.Domain import ElementId


def analyze(directory, dim, interval):
    frames = []
    with h5py.File(directory / "Volume0.h5") as handle:
        groups = sorted(
            (
                g
                for g in handle["Meshes.vol"].values()
                if isinstance(g, h5py.Group)
            ),
            key=lambda g: g.attrs["observation_value"],
        )
        for group in groups:
            names = bytes(group["grid_names"][:]).decode().split(":")
            extents = np.asarray(group["total_extents"]).reshape((-1, dim))
            frame = {}
            for name, ext in zip(names, extents):
                eid = ElementId[dim](name)
                key = (
                    eid.block_id,
                    eid.grid_index,
                    tuple(
                        (s.refinement_level, s.index) for s in eid.segment_ids
                    ),
                )
                frame[key] = tuple(int(n) for n in ext)
            assert len(frame) == len(extents)
            frames.append((float(group.attrs["observation_value"]), frame))
    assert np.allclose(
        np.diff([t for t, f in frames]), interval, atol=1e-12, rtol=0
    )

    def contains(parent, child):
        if parent[:2] != child[:2]:
            return False
        return all(
            lc >= lp and ic >> (lc - lp) == ip
            for (lp, ip), (lc, ic) in zip(parent[2], child[2])
        )

    events = []
    p_signs = defaultdict(set)
    h_signs = defaultdict(set)
    last_p, last_h = {}, {}
    immediate_p = immediate_h = 0
    counts = dict(
        split_events=0,
        join_events=0,
        split_added_elements=0,
        join_removed_elements=0,
        p_increases=0,
        p_decreases=0,
    )
    for step, ((old_t, old), (t, new)) in enumerate(zip(frames, frames[1:])):
        added, removed = set(new) - set(old), set(old) - set(new)
        covered_add, covered_remove = set(), set()
        for action, parents, children in [
            ("split", removed, added),
            ("join", added, removed),
        ]:
            for parent in parents:
                family = {
                    child for child in children if contains(parent, child)
                }
                if not family:
                    continue
                differences = {
                    tuple(
                        lc - lp
                        for (lp, ip), (lc, ic) in zip(parent[2], child[2])
                    )
                    for child in family
                }
                assert len(differences) == 1, (t, parent, family)
                delta = next(iter(differences))
                assert all(d in (0, 1) for d in delta), (t, delta)
                assert len(family) == 2 ** sum(delta), (t, family)
                sign = 1 if action == "split" else -1
                hkey = (parent, delta)
                if hkey in last_h and last_h[hkey] == (step - 1, -sign):
                    immediate_h += 1
                last_h[hkey] = (step, sign)
                h_signs[hkey].add(sign)
                counts[f"{action}_events"] += 1
                counts[
                    (
                        "split_added_elements"
                        if action == "split"
                        else "join_removed_elements"
                    )
                ] += (len(family) - 1)
                events.append(
                    dict(
                        time=t,
                        action=action,
                        element=repr(parent),
                        axes=repr(delta),
                        children=len(family),
                    )
                )
                if action == "split":
                    covered_remove.add(parent)
                    covered_add.update(family)
                else:
                    covered_add.add(parent)
                    covered_remove.update(family)
        assert covered_add == added and covered_remove == removed, (
            t,
            added - covered_add,
            removed - covered_remove,
        )
        for key in old.keys() & new.keys():
            for d, (before, after) in enumerate(zip(old[key], new[key])):
                if before == after:
                    continue
                sign = 1 if after > before else -1
                pkey = (key, d)
                if pkey in last_p and last_p[pkey] == (step - 1, -sign):
                    immediate_p += 1
                last_p[pkey] = (step, sign)
                p_signs[pkey].add(sign)
                counts["p_increases" if sign == 1 else "p_decreases"] += 1
                events.append(
                    dict(
                        time=t,
                        action="p+" if sign == 1 else "p-",
                        element=repr(key),
                        axes=d,
                        children=0,
                    )
                )
    assert (
        len(frames[-1][1]) - len(frames[0][1])
        == counts["split_added_elements"] - counts["join_removed_elements"]
    )
    counts.update(
        snapshot_interval=interval,
        snapshots=len(frames),
        p_scope="directions on unchanged element IDs only",
        p_directions_with_both_signs=sum(len(x) == 2 for x in p_signs.values()),
        h_parent_axis_families_with_both_signs=sum(
            len(x) == 2 for x in h_signs.values()
        ),
        immediate_p_reversals=immediate_p,
        immediate_h_reversals=immediate_h,
    )
    (directory / "MeshChanges.json").write_text(
        json.dumps(counts, indent=2) + "\n"
    )
    if events:
        with (directory / "MeshChanges.csv").open("w") as f:
            writer = csv.DictWriter(f, fieldnames=list(events[0]))
            writer.writeheader()
            writer.writerows(events)
    print(directory.name, json.dumps(counts), flush=True)
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--dim", type=int, default=3)
    parser.add_argument("--interval", type=float, default=0.1)
    args = parser.parse_args()
    analyze(args.directory, args.dim, args.interval)
