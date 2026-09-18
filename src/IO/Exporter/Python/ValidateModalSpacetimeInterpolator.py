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
):
    """Compare modal spacetime interpolation with held-out observations.

    For cluster-scale data, invoke this command once per node-written H5 file
    in a job array. The CSV contains sums of squares and maxima that can be
    reduced across shards.
    """
    h5_files = _resolve_files(h5_files_or_globs)
    training_subfiles = list(map(_subfile_path, training_subfiles))
    validation_subfile = _subfile_path(validation_subfile)
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    observation_ids, observation_times, available_components = (
        _validation_metadata(h5_files[0], validation_subfile)
    )
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
    )
    construction_seconds = time.perf_counter() - construction_start
    click.echo(
        (
            f"Interpolator constructed in {construction_seconds:.1f} s; "
            f"valid on {interpolator.time_bounds()}."
        ),
        err=True,
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
            for filename in h5_files:
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
                    for grid_name, grid_extents in zip(grid_names, extents):
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
                                    result[component_index]
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

            for component_index, component in enumerate(tensor_components):
                reference = np.asarray(references[component_index])
                result = np.asarray(results[component_index])
                summary = _summarize(reference, result)
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
