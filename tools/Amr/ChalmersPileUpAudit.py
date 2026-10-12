# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Audit saved GH spectra for terminal plateaus, without changing evolution.

The SpEC pile-up count is retained as a diagnostic. An independent two-window
check distinguishes decaying-then-flat spectra from parity zeros and wholly
unresolved spectra. Its thresholds are exploratory, not a calibrated AMR rule.
All fits omit the highest filtered mode. Native decision comparisons are only
valid on elements unchanged by AMR at the observation time.
"""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

from spectre.DataStructures import DataVector
from spectre.NumericalAlgorithms.LinearOperators import (
    chalmers_p_coarsening_error,
    chalmers_parity_estimate,
    convergence_rate_and_number_of_pile_up_modes,
)
from spectre.Spectral import Basis, Mesh3D, Quadrature
from spectre.tools.Amr import ChalmersSingleHole as gh
from spectre.tools.Amr import ChalmersSphere as sphere
from spectre.tools.Amr.Chalmers import to_modal


def linear_fit(x, y):
    dx = x - np.mean(x)
    slope = float(np.dot(dx, y - np.mean(y)) / np.dot(dx, dx))
    residual = y - np.mean(y) - slope * dx
    return -slope, float(np.sqrt(np.mean(residual**2)))


def parity_diagnostics(spectrum, budget, parity):
    """Inspect one parity, retaining original polynomial degree spacing.

    spectrum has already had its filtered highest mode removed. Negligible
    parity and terminated-polynomial checks follow the Chalmers conventions.
    At least five significant consecutive parity samples are needed for the
    two-window test: a prefix of >=3 and a tail of 3 (one may overlap).
    Try starts 0, 1 and 2 in the parity sequence so exceptional low modes do
    not spoil a well-defined decay farther into the spectrum.
    """
    k = np.arange(parity, len(spectrum), 2)
    a = spectrum[k]
    varying = k > 0
    variation = float(
        np.linalg.norm(a[varying] * np.sqrt(2 / (2 * k[varying] + 1)))
    )
    floor = 100 * np.finfo(float).eps * np.max(spectrum)
    above = a > floor
    result = dict(
        parity=parity, variation=variation, samples=len(k), floor=float(floor)
    )
    if variation <= 0.01 * budget / np.sqrt(2):
        return dict(result, status="negligible_parity")
    if not np.any(above) or (len(k) >= 2 and not np.any(above[-2:])):
        return dict(result, status="terminated_or_roundoff_tail")
    if len(k) < 5 or not np.all(above):
        return dict(result, status="insufficient_significant_modes")
    log_a = np.log10(a)
    tail, tail_rms = linear_fit(k[-3:], log_a[-3:])
    fits = []
    for start in range(min(3, len(k) - 4)):
        slope, residual = linear_fit(k[start:-2], log_a[start:-2])
        drop = slope * (k[-3] - k[start])
        fits.append((residual, start, slope, float(drop)))
    usable = [fit for fit in fits if fit[2] > 0 and fit[3] >= 0.25]
    rms, start, prefix, prefix_drop = min(usable or fits)
    ratio = abs(tail) / prefix if prefix > 0 else None
    # A loose screening envelope is saved as well as the stricter candidate.
    # This prevents a negative result depending solely on one sharp threshold.
    loose = (
        prefix > 0
        and prefix_drop >= 0.25
        and rms <= 0.3
        and ratio <= 0.3
        and np.ptp(log_a[-3:]) <= 0.75
    )
    candidate = (
        loose
        and prefix_drop >= 0.5
        and rms <= 0.2
        and ratio <= 0.2
        and np.ptp(log_a[-3:]) <= 0.5
    )
    strict = candidate and rms <= 0.1 and ratio <= 0.1
    native = convergence_rate_and_number_of_pile_up_modes(DataVector(a), 0)
    return dict(
        result,
        status=(
            "terminal_plateau_candidate" if candidate else "no_terminal_plateau"
        ),
        loose_candidate=bool(loose),
        strict_candidate=bool(strict),
        prefix_log10_slope=prefix,
        prefix_fit_rms=rms,
        prefix_start_degree=int(k[start]),
        prefix_drop_decades=float(prefix_drop),
        tail_log10_slope=tail,
        tail_fit_rms=tail_rms,
        slope_ratio=ratio,
        tail_norm=float(np.linalg.norm(a[-3:] * np.sqrt(2 / (2 * k[-3:] + 1)))),
        # Native count is a count of parity samples, not original degrees.
        parity_multifit_log10_slope=float(native["convergence_rate"]) / 2,
        parity_pile_up_samples=float(native["number_of_pile_up_modes"]),
    )


def energy_spectra(modes):
    """Transverse Legendre L2 norms, matching directional_spectrum in C++."""
    w = [2 / (2 * np.arange(n) + 1) for n in modes.shape[1:]]
    squared = modes**2
    return [
        np.sqrt(np.einsum("fijk,j,k->fi", squared, w[1], w[2])),
        np.sqrt(np.einsum("fijk,i,k->fj", squared, w[0], w[2])),
        np.sqrt(np.einsum("fijk,i,j->fk", squared, w[0], w[1])),
    ]


def refinement(record):
    error = record["chalmers_error"]
    if error is not None and error <= record["effective_target"]:
        return "None"
    sigma = record["chalmers_sigma"]
    if sigma is None:
        return "Undetermined"
    return "IncreaseResolution" if sigma > 0.5 else "Split"


def aggregate(flags):
    return max(
        flags,
        key={
            "None": 0,
            "IncreaseResolution": 1,
            "Split": 2,
            "Undetermined": 3,
        }.__getitem__,
    )


def frame(group):
    names = bytes(group["grid_names"][:]).decode().split(":")
    extents = np.asarray(group["total_extents"]).reshape(-1, 3)
    return dict(zip(names, [tuple(map(int, n)) for n in extents]))


def self_test():
    # Genuine decay followed by a plateau; no mode grouping artefact.
    k = np.arange(18)
    plateau = np.exp(-0.8 * np.minimum(k, 12))
    assert (
        parity_diagnostics(plateau, 1e-6, 0)["status"]
        == "terminal_plateau_candidate"
    )
    exceptional_low_mode = plateau.copy()
    exceptional_low_mode[0] = 1000
    assert (
        parity_diagnostics(exceptional_low_mode, 1e-6, 0)["status"]
        == "terminal_plateau_candidate"
    )
    assert (
        parity_diagnostics(np.exp(-0.8 * k), 1e-6, 0)["status"]
        == "no_terminal_plateau"
    )
    assert (
        parity_diagnostics(np.ones(18), 1e-6, 0)["status"]
        == "no_terminal_plateau"
    )
    even = np.exp(-0.8 * k)
    even[1::2] = 1e-18
    assert parity_diagnostics(even, 1e-6, 1)["status"] == "negligible_parity"
    assert parity_diagnostics(even, 1e-6, 0)["status"] == "no_terminal_plateau"
    assert (
        parity_diagnostics(1e-15 * plateau, 1e-6, 0)["status"]
        == "negligible_parity"
    )
    assert (
        parity_diagnostics(np.exp(-np.arange(8)), 1e-6, 0)["status"]
        == "insufficient_significant_modes"
    )
    terminated = np.exp(-0.8 * k)
    terminated[10:] = 1e-18
    assert (
        parity_diagnostics(terminated, 1e-6, 0)["status"]
        == "terminated_or_roundoff_tail"
    )
    # Tensor-product Parseval identity for every directional energy spectrum.
    modes = np.random.default_rng(1729).normal(size=(2, 8, 9, 10))
    weights = [2 / (2 * np.arange(n) + 1) for n in modes.shape[1:]]
    norm = np.einsum("fijk,i,j,k->f", modes**2, *weights)
    for d, spectrum in enumerate(energy_spectra(modes)):
        np.testing.assert_allclose(
            np.sum(spectrum**2 * weights[d], axis=1), norm
        )
    return dict(analytic_cases=9, directional_parseval_checks=3, passed=True)


def analyze(directory, output):
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((directory / "Run.json").read_text())
    assert metadata["returncode"] == 0
    assert metadata["criterion"] == "parity" and metadata["tolerance"] > 0
    assert metadata["half_power"] == 256 and metadata["filter"] == "hypercube"
    native = json.loads((directory / "FeedbackAnalysis.json").read_text())[
        "native_evaluations"
    ]
    native_index = {(round(r["time"], 8), r["element"]): r for r in native}
    assert len(native_index) == len(native)
    records = defaultdict(list)
    with gzip.open(directory / "FeedbackSpectra.json.gz", "rt") as f:
        for line in f:
            record = json.loads(line)
            records[(round(record["time"], 8), record["element"])].append(
                record
            )
    counts, parity_counts, raw_classes = Counter(), Counter(), Counter()
    native_counts = Counter(a.strip() for row in native for a in row["after"])
    coverage = Counter()
    directions, examples = [], {}
    max_power_difference = max_error_difference = max_sigma_difference = 0.0
    with (
        h5py.File(directory / "Volume0.h5") as h5,
        gzip.open(output / "SpectraAudit.json.gz", "wt") as out,
    ):
        meshes = sorted(
            (g for g in h5["Meshes.vol"].values() if isinstance(g, h5py.Group)),
            key=lambda g: g.attrs["observation_value"],
        )
        old_frames = {}
        previous = {}
        for group in meshes:
            time = round(float(group.attrs["observation_value"]), 8)
            old_frames[time] = previous
            previous = frame(group)
        for group in sphere.groups(h5):
            time = round(float(group.attrs["observation_value"]), 8)
            for name, ext, xyz, fields in gh.elements(group):
                key = (time, name)
                saved = records.pop(key)
                assert len(saved) == 150
                by_key = {(r["field"], r["direction"]): r for r in saved}
                budget = min(r["effective_target"] for r in saved)
                modes = np.array([to_modal(v) for v in fields])
                spectra = energy_spectra(modes)
                flags = [[] for _ in range(3)]
                veto_flags = [[] for _ in range(3)]
                loose_veto_flags = [[] for _ in range(3)]
                raw_veto_flags = [[] for _ in range(3)]
                candidate_fields = [[] for _ in range(3)]
                for f, field in enumerate(gh.FIELDS):
                    for d in range(3):
                        record = by_key[(field, d)]
                        spectrum = spectra[d][f][:-1]
                        power = np.sqrt(
                            np.mean(
                                modes[f] ** 2,
                                axis=tuple(j for j in range(3) if j != d),
                            )
                        )
                        scale = max(np.max(power), np.finfo(float).tiny)
                        max_power_difference = max(
                            max_power_difference,
                            float(np.max(abs(power - record["power"])) / scale),
                        )
                        counts["spectra"] += 1
                        variation = float(
                            np.linalg.norm(
                                spectrum[1:]
                                * np.sqrt(
                                    2 / (2 * np.arange(1, len(spectrum)) + 1)
                                )
                            )
                        )
                        if record["negligible"] or variation <= budget:
                            replay = dict(error=0.0, sigma=np.nan)
                        else:
                            replay = chalmers_parity_estimate(
                                DataVector(spectrum), budget
                            )
                        old_error = record["chalmers_error"]
                        if old_error is not None:
                            assert np.isfinite(replay["error"])
                            max_error_difference = max(
                                max_error_difference,
                                abs(replay["error"] - old_error)
                                / max(record["effective_target"], old_error),
                            )
                        else:
                            assert not np.isfinite(replay["error"])
                        if record["chalmers_sigma"] is not None:
                            max_sigma_difference = max(
                                max_sigma_difference,
                                abs(replay["sigma"] - record["chalmers_sigma"]),
                            )
                        diagnostic = [
                            parity_diagnostics(spectrum, budget, p)
                            for p in (0, 1)
                        ]
                        for a in diagnostic:
                            parity_counts[a["status"]] += 1
                        raw = (record["spec_pile_up"] or 0) > 3
                        candidate = any(
                            a["status"] == "terminal_plateau_candidate"
                            for a in diagnostic
                        )
                        loose = any(
                            a.get("loose_candidate", False) for a in diagnostic
                        )
                        strict = any(
                            a.get("strict_candidate", False) for a in diagnostic
                        )
                        decision = refinement(record)
                        counts["refinement_" + decision] += 1
                        for label, active in [
                            ("raw", raw),
                            ("candidate", candidate),
                            ("loose", loose),
                            ("strict", strict),
                        ]:
                            if active:
                                counts[label + "_spectra"] += 1
                                counts[label + "_refinement_" + decision] += 1
                        if raw:
                            raw_classes[
                                (
                                    "direction_below_budget"
                                    if variation <= budget
                                    else "direction_above_budget"
                                )
                            ] += 1
                            raw_classes[
                                (
                                    "has_negligible_parity"
                                    if any(
                                        a["status"] == "negligible_parity"
                                        for a in diagnostic
                                    )
                                    else "both_parities_active"
                                )
                            ] += 1
                        # Deliberately permissive hypothetical veto: even one
                        # candidate parity can suppress this component's p
                        # request. This bounds the possible impact from above.
                        # Other components and every h request remain intact.
                        flags[d].append(decision)
                        veto_flags[d].append(
                            "None"
                            if candidate and decision == "IncreaseResolution"
                            else decision
                        )
                        loose_veto_flags[d].append(
                            "None"
                            if loose and decision == "IncreaseResolution"
                            else decision
                        )
                        raw_veto_flags[d].append(
                            "None"
                            if raw and decision == "IncreaseResolution"
                            else decision
                        )
                        if candidate:
                            candidate_fields[d].append(field)
                        row = dict(
                            time=time,
                            element=name,
                            field=field,
                            direction=d,
                            extent=int(ext[d]),
                            effective_target=record["effective_target"],
                            common_budget=budget,
                            variation=variation,
                            error=old_error,
                            sigma=record["chalmers_sigma"],
                            refinement=decision,
                            raw_pile_up=record["spec_pile_up"],
                            parities=diagnostic,
                        )
                        out.write(json.dumps(row, allow_nan=False) + "\n")
                        labels = []
                        if raw:
                            labels.append("raw_largest")
                        if candidate:
                            labels.append("plateau_largest_error")
                        if loose and decision != "None":
                            labels.append("loose_driving")
                        if decision == "IncreaseResolution":
                            labels.append("p_driving")
                        for label in labels:
                            score = (
                                record["spec_pile_up"]
                                if label == "raw_largest"
                                else (old_error or 0)
                                / record["effective_target"]
                            )
                            if (
                                label not in examples
                                or score > examples[label]["score"]
                            ):
                                examples[label] = dict(
                                    row,
                                    score=score,
                                    power=record["power"],
                                    energy=spectra[d][f].tolist(),
                                )
                native_row = native_index.get(key)
                exact = native_row is not None and old_frames.get(time, {}).get(
                    name
                ) == tuple(ext)
                for d in range(3):
                    original = aggregate(flags[d])
                    hypothetical = aggregate(veto_flags[d])
                    loose_hypothetical = aggregate(loose_veto_flags[d])
                    raw_hypothetical = aggregate(raw_veto_flags[d])
                    native_flag = (
                        native_row["after"][d].strip() if native_row else None
                    )
                    counts["snapshot_directions"] += 1
                    counts["snapshot_" + original] += 1
                    counts["candidate_changed_directions"] += (
                        original != hypothetical
                    )
                    counts["loose_changed_directions"] += (
                        original != loose_hypothetical
                    )
                    counts["raw_changed_directions"] += (
                        original != raw_hypothetical
                    )
                    if exact:
                        coverage["native_unchanged_directions"] += 1
                        coverage["native_unchanged_" + native_flag] += 1
                        expect = (
                            native_flag
                            if native_flag in ("Split", "IncreaseResolution")
                            else "None"
                        )
                        if original != expect:
                            raise ValueError(
                                (key, d, "replay mismatch", original, expect)
                            )
                        coverage["candidate_changed_native_directions"] += (
                            original != hypothetical
                        )
                    elif native_row:
                        coverage["native_changed_mesh_directions_excluded"] += 1
                    else:
                        coverage["no_matching_native_direction"] += 1
                    directions.append(
                        dict(
                            time=time,
                            element=name,
                            direction=d,
                            replay=original,
                            candidate_veto=hypothetical,
                            loose_veto=loose_hypothetical,
                            raw_veto=raw_hypothetical,
                            native_flag=native_flag,
                            native_exact=exact,
                            candidate_fields=candidate_fields[d],
                        )
                    )
            print(directory.name, time, "audited", flush=True)
    assert not records
    assert max_power_difference < 1e-10
    assert max_error_difference < 1e-7
    assert max_sigma_difference < 1e-7
    result = dict(
        run=str(directory.resolve()),
        counts=dict(counts),
        parity_counts=dict(parity_counts),
        raw_classes=dict(raw_classes),
        native_decisions=dict(native_counts),
        coverage=dict(coverage),
        validation=dict(
            max_power_relative_difference=max_power_difference,
            max_error_target_scaled_difference=max_error_difference,
            max_sigma_difference=max_sigma_difference,
        ),
        examples=examples,
        limitations=[
            "Fields sampled at 26 times; AMR occurs every 0.1M.",
            (
                "Saved fields follow AMR; changed elements cannot reconstruct"
                " native pre-AMR spectra."
            ),
            (
                "Counterfactual veto is diagnostic only; no coarsening is"
                " forced and no evolution is rerun."
            ),
            (
                "Fewer than five significant samples per active parity cannot"
                " establish a two-window plateau."
            ),
        ],
    )
    (output / "Audit.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    (output / "Directions.json").write_text(
        json.dumps(directions, indent=2) + "\n"
    )
    return result


def follow_up(directory, output):
    """Check plateau persistence, p-coarsening budgets and decision coverage."""
    candidates = defaultdict(list)
    time_counts, field_counts = Counter(), Counter()
    series = defaultdict(list)
    with gzip.open(output / "SpectraAudit.json.gz", "rt") as f:
        for line in f:
            row = json.loads(line)
            for parity in row["parities"]:
                if parity["status"] != "terminal_plateau_candidate":
                    continue
                key = (row["time"], row["element"], row["direction"])
                candidates[key].append(row["field"])
                time_counts[row["time"]] += 1
                field_counts[row["field"]] += 1
                series[
                    (
                        row["element"],
                        row["direction"],
                        row["field"],
                        parity["parity"],
                    )
                ].append((row["time"], row["extent"]))
    longest_same_mesh = 0
    for samples in series.values():
        current = 0
        previous = None
        for time, extent in sorted(samples):
            current = current + 1 if previous == (time - 1, extent) else 1
            longest_same_mesh = max(longest_same_mesh, current)
            previous = time, extent
    metadata = json.loads((directory / "Run.json").read_text())
    directions = json.loads((output / "Directions.json").read_text())
    by_key = {(r["time"], r["element"], r["direction"]): r for r in directions}
    times = sorted({r["time"] for r in directions})
    native = json.loads((directory / "FeedbackAnalysis.json").read_text())[
        "native_evaluations"
    ]
    brackets = []
    for row in native:
        for d, flag in enumerate(row["after"]):
            if flag.strip() != "IncreaseResolution":
                continue
            before = max(t for t in times if t < row["time"] - 1e-8)
            after = min(t for t in times if t >= row["time"] - 1e-8)
            # Nearest observations only: these are NOT the decision-time fields.
            a = by_key.get((before, row["element"], d))
            b = by_key.get((after, row["element"], d))
            brackets.append(
                dict(
                    time=row["time"],
                    element=row["element"],
                    direction=d,
                    before_time=before,
                    after_time=after,
                    before_fields=None if a is None else a["candidate_fields"],
                    after_fields=None if b is None else b["candidate_fields"],
                )
            )
    coarsening = []
    if candidates:
        feedback = {}
        with (directory / "FeedbackByElement.csv").open() as f:
            for r in csv.DictReader(f):
                feedback[
                    (
                        round(float(r["time"]), 8),
                        r["element"],
                        int(r["direction"]),
                    )
                ] = float(r["factor"])
        with h5py.File(directory / "Volume0.h5") as h5:
            for group in sphere.groups(h5):
                time = round(float(group.attrs["observation_value"]), 8)
                for name, ext, xyz, fields in gh.elements(group):
                    selected = [
                        d for d in range(3) if (time, name, d) in candidates
                    ]
                    if not selected:
                        continue
                    factors = [feedback[(time, name, d)] for d in range(3)]
                    common_budget = (
                        metadata["coarsening_factor"]
                        * metadata["tolerance"]
                        * min(factors)
                        / 3
                    )
                    mesh = Mesh3D(ext, Basis.Legendre, Quadrature.GaussLobatto)
                    errors = np.array(
                        [
                            chalmers_p_coarsening_error(
                                DataVector(v.ravel(order="F")),
                                mesh,
                                common_budget,
                                1,
                                True,
                            )
                            for v in fields
                        ]
                    )
                    for d in selected:
                        field = int(np.argmax(errors[:, d]))
                        budget = (
                            metadata["coarsening_factor"]
                            * metadata["tolerance"]
                            * factors[d]
                            / 3
                        )
                        ratio = float(errors[field, d] / budget)
                        coarsening.append(
                            dict(
                                time=time,
                                element=name,
                                direction=d,
                                maximum_error_over_budget=(
                                    ratio if np.isfinite(ratio) else None
                                ),
                                field=gh.FIELDS[field],
                                passes=bool(ratio <= 1),
                            )
                        )
    finite_ratios = [
        r["maximum_error_over_budget"]
        for r in coarsening
        if r["maximum_error_over_budget"] is not None
    ]
    result = dict(
        candidate_spectra_by_time=dict(time_counts),
        candidate_spectra_by_field=dict(field_counts),
        candidate_element_directions=len(candidates),
        longest_same_extent_consecutive_1M_observations=longest_same_mesh,
        candidate_directions_passing_existing_p_coarsening=sum(
            r["passes"] for r in coarsening
        ),
        minimum_coarsening_error_over_budget=min(finite_ratios, default=None),
        maximum_finite_coarsening_error_over_budget=max(
            finite_ratios, default=None
        ),
        native_p_requests=len(brackets),
        p_requests_with_plateau_at_preceding_observation=sum(
            bool(r["before_fields"]) for r in brackets
        ),
        largest_bracketing_interval=max(
            (r["after_time"] - r["before_time"] for r in brackets), default=None
        ),
        p_request_brackets=brackets,
        candidate_coarsening=coarsening,
    )
    (output / "FollowUp.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="*", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--follow-up-only", action="store_true")
    args = parser.parse_args()
    checks = self_test()
    print(json.dumps(checks), flush=True)
    if args.directories and args.output is None:
        parser.error("--output is required for an audit")
    results = []
    for d in args.directories:
        output = args.output / d.name
        result = (
            json.loads((output / "Audit.json").read_text())
            if args.follow_up_only
            else analyze(d, output)
        )
        extra = follow_up(d, output)
        result["follow_up"] = {
            k: v
            for k, v in extra.items()
            if k not in ("p_request_brackets", "candidate_coarsening")
        }
        results.append(result)
    if results:
        summary = dict(
            self_tests=checks,
            source_sha256=hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest(),
            runs=[
                {k: v for k, v in r.items() if k != "examples"} for r in results
            ],
        )
        (args.output / "Summary.json").write_text(
            json.dumps(summary, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
