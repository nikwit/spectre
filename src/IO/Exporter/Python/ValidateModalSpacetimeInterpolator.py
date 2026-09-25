# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Validate a modal spacetime interpolator against held-out volume data."""

import csv
import glob
import hashlib
import json
import math
import resource
import sys
import time
from pathlib import Path

import click
import numpy as np

import spectre.IO.H5 as spectre_h5
from spectre.IO.Exporter import ModalSpacetimeInterpolator, Verbosity


def _subfile_path(name):
    return "/" + name.removeprefix("/")


def _resolve_files(files_or_globs):
    files = []
    for file_or_glob in files_or_globs:
        matches = glob.glob(file_or_glob)
        if not matches:
            raise click.ClickException(f"No files match '{file_or_glob}'.")
        files.extend(matches)
    return sorted(set(map(str, files)))


def _interior_collapsed_indices(extents, count, seed_key):
    """Choose deterministic, strictly interior grid points."""
    extents = tuple(map(int, extents))
    if any(extent <= 2 for extent in extents):
        raise ValueError(
            "Validation meshes need at least three points in every dimension "
            "to select points away from element boundaries."
        )
    interior_extents = tuple(extent - 2 for extent in extents)
    num_interior = math.prod(interior_extents)
    count = min(count, num_interior)
    seed = int.from_bytes(
        hashlib.sha256(seed_key.encode()).digest()[:8], "little"
    )
    rng = np.random.default_rng(seed)
    interior_collapsed = rng.choice(num_interior, size=count, replace=False)
    interior_indices = np.unravel_index(
        interior_collapsed, interior_extents, order="F"
    )
    indices = tuple(index + 1 for index in interior_indices)
    return np.ravel_multi_index(indices, extents, order="F")


def _max_rss_mib():
    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    divisor = 1024.0**2 if sys.platform == "darwin" else 1024.0
    return max_rss / divisor


def _summarize(reference, result):
    reference = np.asarray(reference)
    result = np.asarray(result)
    error = result - reference
    abs_error = np.abs(error)
    sum_squared_error = float(np.dot(error, error))
    sum_squared_reference = float(np.dot(reference, reference))
    return {
        "sample_count": len(error),
        "sum_squared_error": sum_squared_error,
        "sum_squared_reference": sum_squared_reference,
        "rms_abs_error": math.sqrt(sum_squared_error / len(error)),
        "relative_l2_error": (
            math.sqrt(sum_squared_error / sum_squared_reference)
            if sum_squared_reference > 0.0
            else math.nan
        ),
        "median_abs_error": float(np.quantile(abs_error, 0.5)),
        "p90_abs_error": float(np.quantile(abs_error, 0.9)),
        "p99_abs_error": float(np.quantile(abs_error, 0.99)),
        "max_abs_error": float(np.max(abs_error)),
        "worst_index": int(np.argmax(abs_error)),
    }


def _validation_metadata(filename, subfile_name):
    with spectre_h5.H5File(filename, "r") as h5file:
        volfile = h5file.get_vol(_subfile_path(subfile_name))
        if volfile.get_dimension() != 3:
            raise click.ClickException(
                "The validation utility currently supports 3D data only."
            )
        observation_ids = list(volfile.list_observation_ids())
        if not observation_ids:
            raise click.ClickException(
                f"No observations found in '{subfile_name}' of '{filename}'."
            )
        observation_times = [
            volfile.get_observation_value(obs_id) for obs_id in observation_ids
        ]
        components = list(volfile.list_tensor_components(observation_ids[0]))
    return observation_ids, observation_times, components


@click.command(name="validate-modal-spacetime-interpolator")
@click.argument("h5_files_or_globs", nargs=-1, required=True)
@click.option(
    "--training-subfile",
    "training_subfiles",
    multiple=True,
    default=("VeryCoarse", "Coarse", "Full"),
    show_default=True,
    help="Training subfiles, ordered from low-mode to final-mesh priority.",
)
@click.option(
    "--validation-subfile",
    default="Validation",
    show_default=True,
    help="Subfile containing held-out observations.",
)
@click.option(
    "--var",
    "tensor_components",
    multiple=True,
    help=(
        "Tensor component to validate. May be repeated. By default, validates "
        "all non-coordinate components in the validation subfile."
    ),
)
@click.option(
    "--samples-per-element",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of interior validation points sampled in each element.",
)
@click.option(
    "--seed",
    type=int,
    default=12345,
    show_default=True,
    help="Seed used for deterministic spatial sampling.",
)
@click.option(
    "--start-time",
    type=float,
    help="Optional lower time bound applied to the training data.",
)
@click.option(
    "--end-time",
    type=float,
    help="Optional upper time bound applied to the training data.",
)
@click.option(
    "--max-validation-observations",
    type=click.IntRange(min=1),
    help="Process only the first N validation observations for a smoke test.",
)
@click.option(
    "--verbosity",
    type=click.Choice(["silent", "quiet", "verbose", "debug"]),
    default="quiet",
    show_default=True,
)
@click.option(
    "--observation-batch-size",
    type=click.IntRange(min=1),
    default=16,
    show_default=True,
    help="Observations staged together for contiguous modal-history writes.",
)
@click.option(
    "--save-interpolator",
    type=click.Path(dir_okay=False, writable=True),
    help="Save immediately after construction, before starting validation.",
)
@click.option(
    "--load-interpolator",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Load a saved interpolator instead of reading training data.",
)
@click.option(
    "--build-only",
    is_flag=True,
    help="Construct and save without validation. Requires --save-interpolator.",
)
@click.option(
    "--output",
    "output_filename",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Output CSV file. Run metadata is written to a JSON sidecar.",
)
def validate_modal_spacetime_interpolator_command(
    h5_files_or_globs,
    training_subfiles,
    validation_subfile,
    tensor_components,
    samples_per_element,
    seed,
    start_time,
    end_time,
    max_validation_observations,
    verbosity,
    output_filename,
    observation_batch_size,
    save_interpolator,
    load_interpolator,
    build_only,
):
    """Compare modal spacetime interpolation with held-out observations.

    Pass all node-written H5 files to construct one interpolator. Save it with
    --save-interpolator before validation, or use --build-only to stop after
    saving. Later runs can use --load-interpolator with different validation
    samples or a subset of the saved components, without reading training data.
    Saved files are intended for reuse with the same SpECTRE build.
    """
    if load_interpolator and (save_interpolator or build_only):
        raise click.UsageError(
            "--load-interpolator cannot be combined with "
            "--save-interpolator or --build-only."
        )
    if load_interpolator and (start_time is not None or end_time is not None):
        raise click.UsageError(
            "Training time bounds cannot be changed when loading an"
            " interpolator."
        )
    if build_only and not save_interpolator:
        raise click.UsageError("--build-only requires --save-interpolator.")
    if save_interpolator:
        save_path = Path(save_interpolator)
        if save_path.exists() or Path(str(save_path) + ".partial").exists():
            raise click.ClickException(
                f"Refusing to overwrite '{save_path}' or its partial file."
            )
        save_path.parent.mkdir(parents=True, exist_ok=True)
    h5_files = _resolve_files(h5_files_or_globs)
    training_subfiles = list(map(_subfile_path, training_subfiles))
    validation_subfile = _subfile_path(validation_subfile)
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    (
        observation_ids,
        observation_times,
        available_components,
    ) = _validation_metadata(h5_files[0], validation_subfile)
    coordinate_components = [
        "InertialCoordinates_x",
        "InertialCoordinates_y",
        "InertialCoordinates_z",
    ]
    missing_coordinates = set(coordinate_components) - set(available_components)
    if missing_coordinates:
        raise click.ClickException(
            "The validation subfile is missing coordinate components: "
            + ", ".join(sorted(missing_coordinates))
        )
    interpolator = None
    construction_seconds = 0.0
    load_seconds = 0.0
    serialization_seconds = 0.0
    if load_interpolator:
        click.echo(
            f"Loading interpolator from {load_interpolator}...", err=True
        )
        load_start = time.perf_counter()
        interpolator = ModalSpacetimeInterpolator[3].load(load_interpolator)
        load_seconds = time.perf_counter() - load_start
        click.echo(f"Loaded in {load_seconds:.1f} s.", err=True)
        if not tensor_components:
            tensor_components = interpolator.tensor_components()
    if tensor_components:
        tensor_components = list(tensor_components)
    else:
        tensor_components = [
            component
            for component in available_components
            if component not in coordinate_components
        ]
    missing_components = set(tensor_components) - set(available_components)
    if missing_components:
        raise click.ClickException(
            "The validation subfile is missing requested components: "
            + ", ".join(sorted(missing_components))
        )
    if max_validation_observations is not None:
        observation_ids = observation_ids[:max_validation_observations]
        observation_times = observation_times[:max_validation_observations]

    for filename in h5_files[1:]:
        other_ids, other_times, other_components = _validation_metadata(
            filename, validation_subfile
        )
        if (
            other_ids[: len(observation_ids)] != observation_ids
            or other_times[: len(observation_times)] != observation_times
        ):
            raise click.ClickException(
                f"Validation observations in '{filename}' do not match "
                f"'{h5_files[0]}'."
            )
        missing = set(tensor_components + coordinate_components) - set(
            other_components
        )
        if missing:
            raise click.ClickException(
                f"Validation subfile in '{filename}' is missing: "
                + ", ".join(sorted(missing))
            )

    if interpolator is None:
        click.echo(
            (
                f"Constructing interpolator for {len(h5_files)} file(s), "
                f"{len(tensor_components)} component(s), and subfiles "
                f"{training_subfiles}..."
            ),
            err=True,
        )
        construction_start = time.perf_counter()
        interpolator = ModalSpacetimeInterpolator[3](
            h5_files,
            subfiles_in_priority_order=training_subfiles,
            tensor_components=tensor_components,
            start_time=start_time,
            end_time=end_time,
            verbosity=getattr(Verbosity, verbosity.capitalize()),
            observation_batch_size=observation_batch_size,
        )
        construction_seconds = time.perf_counter() - construction_start
        click.echo(
            (
                f"Interpolator constructed in {construction_seconds:.1f} s; "
                f"valid on {interpolator.time_bounds()}."
            ),
            err=True,
        )

    stored_components = list(interpolator.tensor_components())
    missing = set(tensor_components) - set(stored_components)
    if missing:
        raise click.ClickException(
            "The saved interpolator is missing requested components: "
            + ", ".join(sorted(missing))
        )
    component_indices = [stored_components.index(c) for c in tensor_components]
    if save_interpolator:
        click.echo(f"Saving interpolator to {save_interpolator}...", err=True)
        save_start = time.perf_counter()
        interpolator.save(save_interpolator)
        serialization_seconds = time.perf_counter() - save_start
        save_metadata = {
            "h5_files": h5_files,
            "training_subfiles": training_subfiles,
            "tensor_components": stored_components,
            "time_bounds": list(interpolator.time_bounds()),
            "construction_seconds": construction_seconds,
            "observation_batch_size": observation_batch_size,
            "serialization_seconds": serialization_seconds,
            "file_size_bytes": Path(save_interpolator).stat().st_size,
            "python_executable": sys.executable,
            "format": (
                "SpECTRE native PUP v1; use the same build and architecture"
            ),
        }
        Path(str(save_interpolator) + ".json").write_text(
            json.dumps(save_metadata, indent=2) + "\n"
        )
        click.echo(
            (
                f"Saved interpolator in {serialization_seconds:.1f} s "
                f"({save_metadata['file_size_bytes'] / 1024**3:.3f} GiB)."
            ),
            err=True,
        )
    if build_only:
        return

    time_bounds = interpolator.time_bounds()
    outside_times = [
        value
        for value in observation_times
        if not time_bounds[0] <= value <= time_bounds[1]
    ]
    if outside_times:
        raise click.ClickException(
            f"Validation times {outside_times} are outside the interpolator "
            f"time bounds {time_bounds}."
        )

    fieldnames = [
        "time",
        "component",
        "sample_count",
        "sum_squared_error",
        "sum_squared_reference",
        "rms_abs_error",
        "relative_l2_error",
        "median_abs_error",
        "p90_abs_error",
        "p99_abs_error",
        "max_abs_error",
        "worst_file",
        "worst_element",
        "worst_collapsed_index",
        "worst_x",
        "worst_y",
        "worst_z",
        "worst_reference",
        "worst_result",
    ]
    validation_start = time.perf_counter()
    with output_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        for observation_number, (observation_id, observation_time) in enumerate(
            zip(observation_ids, observation_times), start=1
        ):
            click.echo(
                (
                    f"Validation observation {observation_number}/"
                    f"{len(observation_ids)} at t={observation_time:.16g}..."
                ),
                err=True,
            )
            references = [[] for _ in tensor_components]
            results = [[] for _ in tensor_components]
            locations = []
            for file_number, filename in enumerate(h5_files, start=1):
                file_start = time.perf_counter()
                click.echo(
                    (
                        "  Reading validation file"
                        f" {file_number}/{len(h5_files)}: {filename}..."
                    ),
                    err=True,
                )
                with spectre_h5.H5File(filename, "r") as h5file:
                    volfile = h5file.get_vol(validation_subfile)
                    grid_names = list(volfile.get_grid_names(observation_id))
                    extents = list(volfile.get_extents(observation_id))
                    coordinates = [
                        np.array(
                            volfile.get_tensor_component(
                                observation_id, component
                            ).data,
                            copy=False,
                        )
                        for component in coordinate_components
                    ]
                    field_data = [
                        np.array(
                            volfile.get_tensor_component(
                                observation_id, component
                            ).data,
                            copy=False,
                        )
                        for component in tensor_components
                    ]
                    offset = 0
                    for element_number, (grid_name, grid_extents) in enumerate(
                        zip(grid_names, extents), start=1
                    ):
                        local_indices = _interior_collapsed_indices(
                            grid_extents,
                            samples_per_element,
                            f"{seed}:{grid_name}",
                        )
                        for local_index in local_indices:
                            global_index = offset + int(local_index)
                            target_point = np.array(
                                [
                                    coordinate[global_index]
                                    for coordinate in coordinates
                                ]
                            )
                            result = interpolator.interpolate_to_point(
                                target_point, time=observation_time
                            )
                            for component_index, component_data in enumerate(
                                field_data
                            ):
                                references[component_index].append(
                                    component_data[global_index]
                                )
                                results[component_index].append(
                                    result[component_indices[component_index]]
                                )
                            locations.append(
                                (
                                    filename,
                                    grid_name,
                                    int(local_index),
                                    *map(float, target_point),
                                )
                            )
                        offset += math.prod(grid_extents)
                        if element_number == len(grid_names) or (
                            verbosity in ("verbose", "debug")
                            and element_number % 250 == 0
                        ):
                            click.echo(
                                (
                                    "  Evaluated"
                                    f" {element_number}/{len(grid_names)} "
                                    "elements in file"
                                    f" {file_number}/{len(h5_files)}; "
                                    "elapsed "
                                    f"{time.perf_counter() - file_start:.1f}"
                                    " s."
                                ),
                                err=True,
                            )

            for component_index, component in enumerate(tensor_components):
                reference = np.asarray(references[component_index])
                result = np.asarray(results[component_index])
                summary = _summarize(reference, result)
                click.echo(
                    (
                        f"  {component}: RMS={summary['rms_abs_error']:.6e}, "
                        f"relative L2={summary['relative_l2_error']:.6e}, "
                        f"max={summary['max_abs_error']:.6e} "
                        f"({summary['sample_count']} samples)"
                    ),
                    err=True,
                )
                worst_index = summary.pop("worst_index")
                (
                    worst_file,
                    worst_element,
                    worst_collapsed_index,
                    worst_x,
                    worst_y,
                    worst_z,
                ) = locations[worst_index]
                writer.writerow(
                    {
                        "time": observation_time,
                        "component": component,
                        **summary,
                        "worst_file": worst_file,
                        "worst_element": worst_element,
                        "worst_collapsed_index": worst_collapsed_index,
                        "worst_x": worst_x,
                        "worst_y": worst_y,
                        "worst_z": worst_z,
                        "worst_reference": reference[worst_index],
                        "worst_result": result[worst_index],
                    }
                )
            output_file.flush()

    validation_seconds = time.perf_counter() - validation_start
    metadata = {
        "h5_files": h5_files,
        "training_subfiles": training_subfiles,
        "validation_subfile": validation_subfile,
        "tensor_components": tensor_components,
        "samples_per_element": samples_per_element,
        "seed": seed,
        "training_time_bounds": list(interpolator.time_bounds()),
        "validation_times": observation_times,
        "construction_seconds": construction_seconds,
        "observation_batch_size": (
            observation_batch_size if not load_interpolator else None
        ),
        "load_seconds": load_seconds,
        "serialization_seconds": serialization_seconds,
        "saved_interpolator": save_interpolator,
        "loaded_interpolator": load_interpolator,
        "validation_seconds": validation_seconds,
        "max_rss_mib": _max_rss_mib(),
        "csv_output": str(output_path),
    }
    metadata_path = output_path.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    click.echo(
        (
            f"Wrote {output_path} and {metadata_path}. Construction took "
            f"{construction_seconds:.1f} s, validation took "
            f"{validation_seconds:.1f} s, and peak RSS was "
            f"{metadata['max_rss_mib']:.1f} MiB."
        ),
        err=True,
    )


if __name__ == "__main__":
    validate_modal_spacetime_interpolator_command(
        help_option_names=["-h", "--help"]
    )
