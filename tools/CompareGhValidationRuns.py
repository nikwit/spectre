#!/usr/bin/env python3

# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Compare CPU and batched-GPU validation runs of the same GH evolution.

Compares (1) reduction norms/errors slab by slab and (2) volume data
pointwise per element and component. Exits nonzero if any comparison
exceeds its tolerance.
"""

import argparse
import sys

import h5py
import numpy as np


def read_dat(filename, subfile):
    with h5py.File(filename, "r") as f:
        dataset = f[subfile]
        legend = [
            s.decode() if isinstance(s, bytes) else s
            for s in dataset.attrs["Legend"]
        ]
        data = np.array(dataset)
    return legend, data


def normalize_column(name):
    """Map production and batched norm column names to a common form.

    Production ObserveNorms writes 'L2Norm(X)'; ObserveNormsBatched
    writes 'X_L2Norm'. Columns without a counterpart (NumberOfPoints,
    Volume) are skipped by the caller.
    """
    for norm in ("L2Norm", "Max", "Min"):
        prefix = norm + "("
        if name.startswith(prefix) and name.endswith(")"):
            return name[len(prefix) : -1] + "_" + norm
    return name


def compare_reductions(cpu_file, gpu_file, subfile, rtol):
    cpu_legend, cpu_data = read_dat(cpu_file, subfile)
    gpu_legend, gpu_data = read_dat(gpu_file, subfile)
    cpu_legend = [normalize_column(name) for name in cpu_legend]
    gpu_legend = [normalize_column(name) for name in gpu_legend]
    failures = 0
    if cpu_data.shape[0] != gpu_data.shape[0]:
        print(
            f"  {subfile}: row count differs "
            f"(cpu {cpu_data.shape[0]}, gpu {gpu_data.shape[0]})"
        )
        return 1
    # Match columns by name; the time column anchors row alignment.
    cpu_time = cpu_data[:, cpu_legend.index("Time")]
    gpu_time = gpu_data[:, gpu_legend.index("Time")]
    if not np.allclose(cpu_time, gpu_time, rtol=0, atol=1e-14):
        print(f"  {subfile}: observation times differ")
        return 1
    for name in cpu_legend:
        if name in ("Time", "NumberOfPoints", "Volume"):
            continue
        if name not in gpu_legend:
            print(f"  {subfile}: column '{name}' missing in GPU output")
            failures += 1
            continue
        cpu_col = cpu_data[:, cpu_legend.index(name)]
        gpu_col = gpu_data[:, gpu_legend.index(name)]
        scale = np.maximum(np.abs(cpu_col), np.abs(gpu_col))
        scale[scale == 0.0] = 1.0
        rel = np.abs(cpu_col - gpu_col) / scale
        worst = np.max(rel)
        status = "OK " if worst <= rtol else "FAIL"
        if worst > rtol:
            failures += 1
        print(f"  [{status}] {subfile} {name}: max rel diff {worst:.3e}")
    return failures


def compare_dat_pair(cpu_file, gpu_file, cpu_subfile, gpu_subfile, rtol, atol):
    """Compare corresponding generic reduction tables by column name."""
    cpu_legend, cpu_data = read_dat(cpu_file, cpu_subfile)
    gpu_legend, gpu_data = read_dat(gpu_file, gpu_subfile)
    failures = 0
    label = f"{cpu_subfile} : {gpu_subfile}"
    if cpu_data.shape[0] != gpu_data.shape[0]:
        print(
            f"  {label}: row count differs "
            f"(cpu {cpu_data.shape[0]}, gpu {gpu_data.shape[0]})"
        )
        return 1
    if "Time" not in cpu_legend or "Time" not in gpu_legend:
        print(f"  {label}: missing Time column")
        return 1
    cpu_time = cpu_data[:, cpu_legend.index("Time")]
    gpu_time = gpu_data[:, gpu_legend.index("Time")]
    if not np.allclose(cpu_time, gpu_time, rtol=0, atol=atol):
        print(f"  {label}: observation times differ")
        return 1
    for name in cpu_legend:
        if name == "Time":
            continue
        if name not in gpu_legend:
            print(f"  {label}: column '{name}' missing in GPU output")
            failures += 1
            continue
        cpu_col = cpu_data[:, cpu_legend.index(name)]
        gpu_col = gpu_data[:, gpu_legend.index(name)]
        difference = np.abs(cpu_col - gpu_col)
        allowed = atol + rtol * np.maximum(np.abs(cpu_col), np.abs(gpu_col))
        scaled_error = np.max(difference / allowed)
        worst_absolute = np.max(difference)
        status = "OK " if scaled_error <= 1.0 else "FAIL"
        if scaled_error > 1.0:
            failures += 1
        print(
            f"  [{status}] {label} {name}: max abs diff "
            f"{worst_absolute:.3e}, tolerance ratio {scaled_error:.3e}"
        )
    return failures


def read_volume_observations(filename, subfile):
    """Return {obs_value: {grid_name: {component: array}}}."""
    result = {}
    with h5py.File(filename, "r") as f:
        vol = f[subfile]
        for obs_group_name in vol:
            group = vol[obs_group_name]
            # Standard observer files also store domain metadata directly
            # under the volume subfile. Only ObservationId groups carry an
            # observation value.
            if (
                not isinstance(group, h5py.Group)
                or "observation_value" not in group.attrs
            ):
                continue
            obs_value = group.attrs["observation_value"]
            raw_names = (
                np.array(group["grid_names"]).tobytes().decode().strip("\x00")
            )
            # Grid names are ElementId strings like [B0,(L0I0,L0I0,L0I0)]
            # concatenated with ':' separators (no colon occurs inside a
            # name).
            grid_names = [s for s in raw_names.split(":") if s]
            extents = np.array(group["total_extents"])
            # Points per grid: product of extents per grid.
            dim = len(extents) // len(grid_names)
            points_per_grid = [
                int(np.prod(extents[i * dim : (i + 1) * dim]))
                for i in range(len(grid_names))
            ]
            offsets = np.concatenate(([0], np.cumsum(points_per_grid)))
            components = {}
            for key in group:
                if key in (
                    "BlockId",
                    "ElementId",
                    "grid_names",
                    "total_extents",
                    "bases",
                    "quadratures",
                    "domain",
                    "functions_of_time",
                    "connectivity",
                    "pole_connectivity",
                ):
                    continue
                components[key] = np.array(group[key])
            grids = {}
            for i, grid in enumerate(grid_names):
                grids[grid] = {
                    comp: data[offsets[i] : offsets[i + 1]]
                    for comp, data in components.items()
                }
            result[float(obs_value)] = grids
    return result


def compare_volumes(cpu_file, gpu_file, subfile, rtol, coord_atol):
    cpu_obs = read_volume_observations(cpu_file, subfile)
    gpu_obs = read_volume_observations(gpu_file, subfile)
    failures = 0
    cpu_times = sorted(cpu_obs)
    gpu_times = sorted(gpu_obs)
    if len(cpu_times) != len(gpu_times) or not np.allclose(
        cpu_times, gpu_times, rtol=0, atol=1e-14
    ):
        print(
            f"  volume observation times differ: cpu {cpu_times}, "
            f"gpu {gpu_times}"
        )
        return 1
    for time, gpu_time in zip(cpu_times, gpu_times):
        cpu_grids = cpu_obs[time]
        gpu_grids = gpu_obs[gpu_time]
        if set(cpu_grids) != set(gpu_grids):
            only_cpu = set(cpu_grids) - set(gpu_grids)
            only_gpu = set(gpu_grids) - set(cpu_grids)
            print(
                f"  t={time}: grid sets differ (only cpu: {only_cpu}, "
                f"only gpu: {only_gpu})"
            )
            failures += 1
            continue
        worst_by_component = {}
        for grid in cpu_grids:
            cpu_comps = cpu_grids[grid]
            gpu_comps = gpu_grids[grid]
            shared = set(cpu_comps) & set(gpu_comps)
            for comp in shared:
                cpu_vals = cpu_comps[comp]
                gpu_vals = gpu_comps[comp]
                if cpu_vals.shape != gpu_vals.shape:
                    print(
                        f"  t={time} {grid} {comp}: shape mismatch "
                        f"{cpu_vals.shape} vs {gpu_vals.shape}"
                    )
                    failures += 1
                    continue
                if comp.startswith("InertialCoordinates"):
                    # Coordinates must agree to roundoff; they anchor
                    # the point ordering.
                    diff = np.max(np.abs(cpu_vals - gpu_vals))
                    if diff > coord_atol:
                        print(
                            f"  t={time} {grid} {comp}: coordinate "
                            f"mismatch {diff:.3e} (points misordered?)"
                        )
                        failures += 1
                    continue
                scale = np.maximum(np.abs(cpu_vals), np.abs(gpu_vals))
                scale[scale == 0.0] = 1.0
                rel = np.max(np.abs(cpu_vals - gpu_vals) / scale)
                worst_by_component[comp] = max(
                    worst_by_component.get(comp, 0.0), rel
                )
        for comp in sorted(worst_by_component):
            worst = worst_by_component[comp]
            status = "OK " if worst <= rtol else "FAIL"
            if worst > rtol:
                failures += 1
            print(
                f"  [{status}] t={time} {comp}: max rel diff over all "
                f"grids {worst:.3e}"
            )
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu-reductions", required=True)
    parser.add_argument("--gpu-reductions", required=True)
    parser.add_argument("--cpu-volume", default=None)
    parser.add_argument("--gpu-volume", default=None)
    parser.add_argument("--volume-subfile", default="VolumeData.vol")
    parser.add_argument(
        "--rtol",
        type=float,
        default=1e-10,
        help=(
            "Relative tolerance for norms and fields; should be far "
            "below truncation error but above accumulated roundoff."
        ),
    )
    parser.add_argument("--coord-atol", type=float, default=1e-13)
    parser.add_argument(
        "--reduction-subfiles",
        nargs="+",
        default=("Norms.dat", "Errors.dat"),
        help="Reduction subfiles to compare (default: Norms.dat Errors.dat).",
    )
    parser.add_argument(
        "--dat-pair",
        action="append",
        default=[],
        metavar="CPU_SUBFILE=GPU_SUBFILE",
        help=(
            "Additional corresponding reduction tables to compare. May be "
            "specified multiple times, for example for horizon or control "
            "system diagnostics."
        ),
    )
    parser.add_argument(
        "--dat-atol",
        type=float,
        default=1e-13,
        help="Absolute tolerance for --dat-pair comparisons.",
    )
    args = parser.parse_args()

    failures = 0
    print("Reduction comparison:")
    for subfile in args.reduction_subfiles:
        failures += compare_reductions(
            args.cpu_reductions, args.gpu_reductions, subfile, args.rtol
        )
    if args.dat_pair:
        print("Additional reduction-table comparison:")
        for pair in args.dat_pair:
            if "=" not in pair:
                parser.error(
                    "--dat-pair must have the form CPU_SUBFILE=GPU_SUBFILE: "
                    f"{pair}"
                )
            cpu_subfile, gpu_subfile = pair.split("=", maxsplit=1)
            failures += compare_dat_pair(
                args.cpu_reductions,
                args.gpu_reductions,
                cpu_subfile,
                gpu_subfile,
                args.rtol,
                args.dat_atol,
            )
    if args.cpu_volume and args.gpu_volume:
        print("Volume comparison:")
        failures += compare_volumes(
            args.cpu_volume,
            args.gpu_volume,
            args.volume_subfile,
            args.rtol,
            args.coord_atol,
        )
    print("PASSED" if failures == 0 else f"FAILED ({failures} comparisons)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
