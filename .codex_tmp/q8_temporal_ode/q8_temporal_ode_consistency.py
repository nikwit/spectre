#!/usr/bin/env python3
"""Test temporal ODE consistency of the q=8 matching coefficients.

The first-order matching coefficients are fitted independently on each of the
ten downloaded volume observations.  The seven-parameter second-order
subsystem (six symmetric spatial-strain rates and one time acceleration) is
then fitted and alternated with the first-order solve.

The resulting time series test

    qddot0 = (1 + qdot0) d(qdot0)/dT
    dot(L) = M_B (1 + qdot0) d(Lambda)/dT,

where T is the simulation observation time, Lambda=L/M_B, and overdots in the
matching formula denote derivatives with respect to the centered harmonic
time t in T=t+q0(t).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/q8-temporal-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/q8-temporal-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

SPECTRE_BUILD = Path(
    os.environ.get(
        "SPECTRE_BUILD", "/Users/niko/caltech/spectre/build-Release"
    )
)
sys.path.insert(0, str(SPECTRE_BUILD / "bin/python"))
sys.path.insert(0, str(SPECTRE_BUILD / "lib/python3.13/site-packages"))

AUDIT_ROOT = Path(
    os.environ.get(
        "GAUGE_MATCHING_AUDIT",
        "/Users/niko/caltech/worldtube-gravity/gauge-matching-audit",
    )
)
sys.path.insert(0, str(AUDIT_ROOT / "python"))

Q8_DIRECTORY = Path(
    os.environ.get(
        "Q8DIR", "/Users/niko/caltech/simulations/harmonic-sims/q8"
    )
)
sys.path.insert(0, str(Q8_DIRECTORY))

import spectre.IO.H5 as spectre_h5
from second_order_responses import (
    FULL_GROUP_SLICES,
    SYMMETRIC_BASIS,
    full_response_columns_at_points,
    known_quadratic_first_order_source,
)
from spectre.DataStructures import DataVector
from spectre.DataStructures.Tensor import Frame, tnsr
from spectre.IO.Exporter import ObservationId, interpolate_to_points
from spectre.SphericalHarmonics import Spherepack, SpherepackIterator
from worldtube_matching.responses import (
    first_order_response_columns,
    inverse_harmonic_schwarzschild,
)


VOLUME_FILES = [str(Q8_DIRECTORY / "subset.h5")]
REDUCTIONS_FILE = str(Q8_DIRECTORY / "BbhReductions.h5")
VOLUME_SUBFILE = "MetricAndScalars"
OUTPUT_DIRECTORY = Path(
    os.environ.get(
        "Q8_TEMPORAL_OUTPUT",
        "/Users/niko/caltech/spectre/.codex_tmp/q8_temporal_ode",
    )
)
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

SMALL_BLACK_HOLE_MASS = 1.0 / 9.0
SPHEREPACK_L_MAX = 16
CENTER_FINDING_RADIUS = 0.20
RADII = np.array([0.20, 0.24, 0.28, 0.34, 0.40, 0.48, 0.58, 0.68, 0.80])
NUMBER_OF_ALTERNATIONS = 4
VECTOR_DERIVATIVE_WEIGHT = 3.0
SECOND_ORDER_DERIVATIVE_WEIGHT = 1.0

INVERSE_METRIC_COMPONENT_NAMES = (
    "InverseSpacetimeMetric_tt",
    "InverseSpacetimeMetric_xt",
    "InverseSpacetimeMetric_yt",
    "InverseSpacetimeMetric_zt",
    "InverseSpacetimeMetric_xx",
    "InverseSpacetimeMetric_yx",
    "InverseSpacetimeMetric_yy",
    "InverseSpacetimeMetric_zx",
    "InverseSpacetimeMetric_zy",
    "InverseSpacetimeMetric_zz",
)
SPATIAL_COMPONENTS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 1),
    (1, 2),
    (2, 2),
)
DATA_SPATIAL_COMPONENTS = (
    (0, 0),
    (0, 1),
    (1, 1),
    (0, 2),
    (1, 2),
    (2, 2),
)
COVARIANT_COMPONENT_SUFFIXES = (
    "tt",
    "xt",
    "xx",
    "yt",
    "yx",
    "yy",
    "zt",
    "zx",
    "zy",
    "zz",
)
COVARIANT_DATA_COMPONENTS = (
    (0, 0),
    (0, 1),
    (1, 1),
    (0, 2),
    (1, 2),
    (2, 2),
    (0, 3),
    (1, 3),
    (2, 3),
    (3, 3),
)
PHI_COMPONENT_NAMES = tuple(
    f"Phi_{derivative}{suffix}"
    for derivative in ("x", "y", "z")
    for suffix in COVARIANT_COMPONENT_SUFFIXES
)
VECTOR_PARAMETER_INDICES = np.arange(1, 7)
CLOCK_STRAIN_PARAMETER_INDICES = np.r_[0, np.arange(7, 13)]
VECTOR_BLOCKS = (("GTT", 1), ("GTI", 0))
CLOCK_STRAIN_BLOCKS = (("GTT", 0), ("GIJ", 2))
SECOND_ORDER_INDICES = np.r_[
    np.arange(
        FULL_GROUP_SLICES["symmetric spatial-strain rate"].start,
        FULL_GROUP_SLICES["symmetric spatial-strain rate"].stop,
    ),
    np.arange(
        FULL_GROUP_SLICES["time acceleration"].start,
        FULL_GROUP_SLICES["time acceleration"].stop,
    ),
]


def spherepack_grid(ell_max):
    spherepack = Spherepack(ell_max, ell_max)
    theta, phi = map(np.asarray, spherepack.theta_phi_points)
    directions = np.column_stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        )
    )
    return spherepack, directions


SPHEREPACK, DIRECTIONS = spherepack_grid(SPHEREPACK_L_MAX)


def spherical_harmonic_coefficients(values, ell_max):
    spectral = np.asarray(
        SPHEREPACK.phys_to_spec(DataVector(np.asarray(values)))
    )
    iterator = SpherepackIterator(SPHEREPACK_L_MAX, SPHEREPACK_L_MAX)
    indices = []
    for ell in range(ell_max + 1):
        for m in range(-ell, ell + 1):
            iterator.set(ell, m)
            indices.append(iterator())
    return spectral[indices]


def ylm_to_cartesian_matrix():
    nx, ny, nz = DIRECTIONS.T
    basis = np.column_stack(
        (
            np.ones_like(nx),
            nx,
            ny,
            nz,
            nx * ny,
            nx * nz,
            ny * nz,
            nx * nx - ny * ny,
            2.0 * nz * nz - nx * nx - ny * ny,
        )
    )
    transform = np.column_stack(
        [
            spherical_harmonic_coefficients(basis[:, column], 2)
            for column in range(basis.shape[1])
        ]
    )
    return np.linalg.inv(transform)


YLM_TO_CARTESIAN = ylm_to_cartesian_matrix()


def cartesian_dipole(values):
    coefficients = (
        YLM_TO_CARTESIAN @ spherical_harmonic_coefficients(values, 2)
    )
    return coefficients[1:4]


def volume_observations():
    with spectre_h5.H5File(VOLUME_FILES[0], "r") as h5file:
        volume = h5file.get_vol("/" + VOLUME_SUBFILE)
        return sorted(
            (
                volume.get_observation_value(observation_id),
                observation_id,
            )
            for observation_id in volume.list_observation_ids()
        )


def horizon_b_center(time):
    with spectre_h5.H5File(REDUCTIONS_FILE, "r") as h5file:
        centers = np.asarray(
            h5file.get_dat(
                "/ApparentHorizons/ControlSystemAhB_Centers"
            ).get_data()
        )
    return np.array(
        [
            np.interp(time, centers[:, 0], centers[:, component])
            for component in range(4, 7)
        ]
    )


def interpolate_components(observation_id, points, component_names):
    points = np.asarray(points)
    target_points = tnsr.I[DataVector, 3, Frame.Inertial](
        points.T.copy()
    )
    return np.asarray(
        interpolate_to_points(
            VOLUME_FILES,
            subfile_name=VOLUME_SUBFILE,
            observation=ObservationId(observation_id),
            tensor_components=list(component_names),
            target_points=target_points,
            error_on_missing_points=True,
        )
    )


def iterate_curvature_center(
    observation_id,
    initial_center,
    radius=CENTER_FINDING_RADIUS,
    maximum_iterations=10,
    tolerance=1.0e-11,
):
    center = np.asarray(initial_center, dtype=float).copy()
    correction = np.full(3, np.nan)
    for iteration in range(maximum_iterations):
        points = center[None, :] + radius * DIRECTIONS
        kretschmann = interpolate_components(
            observation_id, points, ("GaussBonnetScalar",)
        )[0]
        areal_radius = (
            48.0 * SMALL_BLACK_HOLE_MASS**2 / kretschmann
        ) ** (1.0 / 6.0)
        inferred_harmonic_radius = (
            areal_radius - SMALL_BLACK_HOLE_MASS
        )
        residual = radius**2 - inferred_harmonic_radius**2
        correction = cartesian_dipole(residual) / (2.0 * radius)
        center += correction
        if np.linalg.norm(correction) < tolerance:
            return center, iteration + 1, correction
    return center, maximum_iterations, correction


def interpolate_inverse_metric(observation_id, center):
    points = np.concatenate(
        [center[None, :] + radius * DIRECTIONS for radius in RADII]
    )
    values = interpolate_components(
        observation_id, points, INVERSE_METRIC_COMPONENT_NAMES
    )
    return values.reshape(
        len(INVERSE_METRIC_COMPONENT_NAMES),
        len(RADII),
        len(DIRECTIONS),
    )


def unpack_inverse_metric(field_data):
    gtt = field_data[0]
    gti = np.moveaxis(field_data[1:4], 0, -1)
    gij = np.zeros((len(gtt), 3, 3))
    for values, (i, j) in zip(
        field_data[4:], DATA_SPATIAL_COMPONENTS
    ):
        gij[:, i, j] = gij[:, j, i] = values
    return gtt, gti, gij


def inverse_metric_matrix(field_data):
    gtt, gti, gij = unpack_inverse_metric(field_data)
    matrix = np.zeros((len(gtt), 4, 4))
    matrix[:, 0, 0] = gtt
    matrix[:, 0, 1:] = matrix[:, 1:, 0] = gti
    matrix[:, 1:, 1:] = gij
    return matrix


def metric_tuple_from_matrix(matrix):
    return matrix[:, 0, 0], matrix[:, 0, 1:], matrix[:, 1:, 1:]


def interpolate_radial_derivative(
    observation_id, center, inverse_metric_data
):
    points = np.concatenate(
        [center[None, :] + radius * DIRECTIONS for radius in RADII]
    )
    values = interpolate_components(
        observation_id, points, PHI_COMPONENT_NAMES
    ).reshape(
        3,
        len(COVARIANT_DATA_COMPONENTS),
        len(RADII),
        len(DIRECTIONS),
    )
    phi = np.zeros((len(RADII), len(DIRECTIONS), 3, 4, 4))
    for derivative in range(3):
        for component_values, (a, b) in zip(
            values[derivative], COVARIANT_DATA_COMPONENTS
        ):
            phi[:, :, derivative, a, b] = component_values
            phi[:, :, derivative, b, a] = component_values
    inverse_matrices = np.asarray(
        [
            inverse_metric_matrix(inverse_metric_data[:, radius_index])
            for radius_index in range(len(RADII))
        ]
    )
    spatial_derivative = -np.einsum(
        "rpac,rpbd,rpicd->rpiab",
        inverse_matrices,
        inverse_matrices,
        phi,
    )
    return np.einsum("pi,rpiab->rpab", DIRECTIONS, spatial_derivative)


def ell_block(values, ell):
    coefficients = spherical_harmonic_coefficients(values, ell)
    return coefficients[ell**2 : (ell + 1) ** 2]


def metric_stf_block(metric, block):
    gtt, gti, gij = metric
    metric_component, ell = block
    if metric_component == "GTT":
        return ell_block(gtt, ell)
    if metric_component == "GTI":
        return np.concatenate(
            [ell_block(gti[:, component], ell) for component in range(3)]
        )
    if metric_component == "GIJ":
        return np.concatenate(
            [ell_block(gij[:, i, j], ell) for i, j in SPATIAL_COMPONENTS]
        )
    raise ValueError(block)


def projected_design(columns, blocks, parameter_indices):
    return np.concatenate(
        [
            np.column_stack(
                [
                    metric_stf_block(columns[index], block)
                    for index in parameter_indices
                ]
            )
            for block in blocks
        ],
        axis=0,
    )


def projected_data(metric, blocks):
    return np.concatenate(
        [metric_stf_block(metric, block) for block in blocks]
    )


def metric_values_by_point(metric):
    gtt, gti, gij = metric
    return np.column_stack(
        [gtt, gti, *[gij[:, i, j] for i, j in SPATIAL_COMPONENTS]]
    )


def scaled_least_squares(matrix, data):
    scales = np.linalg.norm(matrix, axis=0)
    if np.any(scales == 0.0):
        raise RuntimeError("Zero response column")
    scaled_solution, _, rank, singular_values = np.linalg.lstsq(
        matrix / scales, data, rcond=1.0e-10
    )
    return (
        scaled_solution / scales,
        rank,
        singular_values[0] / singular_values[rank - 1],
    )


def metric_linear_combination(columns, parameters):
    return tuple(
        sum(
            parameter * column[component]
            for parameter, column in zip(parameters, columns)
        )
        for component in range(3)
    )


def subtract_metric(left, *right):
    return tuple(
        left_component
        - sum(item[component] for item in right)
        for component, left_component in enumerate(left)
    )


def add_metrics(*metrics):
    return tuple(
        sum(metric[component] for metric in metrics)
        for component in range(3)
    )


def symmetric_matrix_from_six(values):
    matrix = np.zeros((3, 3))
    for value, (i, j) in zip(values, SPATIAL_COMPONENTS):
        matrix[i, j] = matrix[j, i] = value
    return matrix


def five_point_radial_derivative(function, radius, relative_step=1.0e-4):
    step = relative_step * radius
    samples = [
        function((radius + offset * step) * DIRECTIONS)
        for offset in (-2.0, -1.0, 1.0, 2.0)
    ]
    return tuple(
        (
            samples[0][component]
            - 8.0 * samples[1][component]
            + 8.0 * samples[2][component]
            - samples[3][component]
        )
        / (12.0 * step)
        for component in range(3)
    )


def response_radial_derivatives(response_function, radius):
    step = 1.0e-4 * radius
    samples = [
        response_function((radius + offset * step) * DIRECTIONS)
        for offset in (-2.0, -1.0, 1.0, 2.0)
    ]
    return tuple(
        tuple(
            (
                samples[0][parameter][component]
                - 8.0 * samples[1][parameter][component]
                + 8.0 * samples[2][parameter][component]
                - samples[3][parameter][component]
            )
            / (12.0 * step)
            for component in range(3)
        )
        for parameter in range(len(samples[0]))
    )


def prepare_responses():
    backgrounds = []
    background_derivatives = []
    first_order_columns = []
    first_order_derivatives = []
    second_order_columns = []
    second_order_derivatives = []
    for radius in RADII:
        points = radius * DIRECTIONS
        backgrounds.append(
            inverse_harmonic_schwarzschild(
                points, SMALL_BLACK_HOLE_MASS
            )
        )
        background_derivatives.append(
            tuple(
                radius * component
                for component in five_point_radial_derivative(
                    lambda local_points: inverse_harmonic_schwarzschild(
                        local_points, SMALL_BLACK_HOLE_MASS
                    ),
                    radius,
                )
            )
        )
        first_columns = first_order_response_columns(
            points, SMALL_BLACK_HOLE_MASS
        )
        first_order_columns.append(first_columns)
        first_order_derivatives.append(
            tuple(
                tuple(radius * component for component in column)
                for column in response_radial_derivatives(
                    lambda local_points: first_order_response_columns(
                        local_points, SMALL_BLACK_HOLE_MASS
                    ),
                    radius,
                )
            )
        )
        second_columns = full_response_columns_at_points(
            points, SMALL_BLACK_HOLE_MASS
        )
        second_order_columns.append(second_columns)
        second_order_derivatives.append(
            tuple(
                tuple(radius * component for component in column)
                for column in response_radial_derivatives(
                    lambda local_points: full_response_columns_at_points(
                        local_points, SMALL_BLACK_HOLE_MASS
                    ),
                    radius,
                )
            )
        )
    return (
        backgrounds,
        background_derivatives,
        first_order_columns,
        first_order_derivatives,
        second_order_columns,
        second_order_derivatives,
    )


(
    BACKGROUNDS,
    BACKGROUND_DERIVATIVES,
    FIRST_ORDER_COLUMNS,
    FIRST_ORDER_DERIVATIVES,
    SECOND_ORDER_COLUMNS,
    SECOND_ORDER_DERIVATIVES,
) = prepare_responses()


def quadratic_value_and_derivative(first_order_parameters, radius_index):
    radius = RADII[radius_index]
    beta = first_order_parameters[1:4]
    strain = symmetric_matrix_from_six(first_order_parameters[7:13])
    value = known_quadratic_first_order_source(
        radius * DIRECTIONS, beta, strain, SMALL_BLACK_HOLE_MASS
    )
    derivative = five_point_radial_derivative(
        lambda points: known_quadratic_first_order_source(
            points, beta, strain, SMALL_BLACK_HOLE_MASS
        ),
        radius,
    )
    return value, tuple(radius * component for component in derivative)


def fit_first_order(
    raw_value,
    raw_derivative,
    radius_index,
    correction_value=None,
    correction_derivative=None,
):
    if correction_value is None:
        correction_value = tuple(
            np.zeros_like(component) for component in raw_value
        )
        correction_derivative = tuple(
            np.zeros_like(component) for component in raw_derivative
        )
    result = np.zeros(13)
    sectors = (
        (
            VECTOR_BLOCKS,
            VECTOR_PARAMETER_INDICES,
            VECTOR_DERIVATIVE_WEIGHT,
        ),
        (CLOCK_STRAIN_BLOCKS, CLOCK_STRAIN_PARAMETER_INDICES, 0.0),
    )
    for blocks, indices, derivative_weight in sectors:
        corrected_value = subtract_metric(raw_value, correction_value)
        matrix_parts = [
            projected_design(
                FIRST_ORDER_COLUMNS[radius_index], blocks, indices
            )
        ]
        data_parts = [projected_data(corrected_value, blocks)]
        if derivative_weight > 0.0:
            corrected_derivative = subtract_metric(
                raw_derivative, correction_derivative
            )
            matrix_parts.append(
                derivative_weight
                * projected_design(
                    FIRST_ORDER_DERIVATIVES[radius_index],
                    blocks,
                    indices,
                )
            )
            data_parts.append(
                derivative_weight
                * projected_data(corrected_derivative, blocks)
            )
        estimate, rank, _ = scaled_least_squares(
            np.concatenate(matrix_parts),
            np.concatenate(data_parts),
        )
        if rank != len(indices):
            raise RuntimeError("First-order block lost rank")
        result[indices] = estimate
    return result


def fit_second_order(
    first_order_parameters, raw_value, raw_derivative, radius_index
):
    first_value = metric_linear_combination(
        FIRST_ORDER_COLUMNS[radius_index], first_order_parameters
    )
    first_derivative = metric_linear_combination(
        FIRST_ORDER_DERIVATIVES[radius_index], first_order_parameters
    )
    quadratic_value, quadratic_derivative = quadratic_value_and_derivative(
        first_order_parameters, radius_index
    )
    candidate_value = subtract_metric(
        raw_value, first_value, quadratic_value
    )
    candidate_derivative = subtract_metric(
        raw_derivative, first_derivative, quadratic_derivative
    )
    selected_columns = tuple(
        SECOND_ORDER_COLUMNS[radius_index][index]
        for index in SECOND_ORDER_INDICES
    )
    selected_derivatives = tuple(
        SECOND_ORDER_DERIVATIVES[radius_index][index]
        for index in SECOND_ORDER_INDICES
    )
    value_matrix = np.column_stack(
        [metric_values_by_point(column).ravel() for column in selected_columns]
    )
    derivative_matrix = np.column_stack(
        [
            metric_values_by_point(column).ravel()
            for column in selected_derivatives
        ]
    )
    estimate, rank, _ = scaled_least_squares(
        np.concatenate(
            (
                value_matrix,
                SECOND_ORDER_DERIVATIVE_WEIGHT * derivative_matrix,
            )
        ),
        np.concatenate(
            (
                metric_values_by_point(candidate_value).ravel(),
                SECOND_ORDER_DERIVATIVE_WEIGHT
                * metric_values_by_point(candidate_derivative).ravel(),
            )
        ),
    )
    if rank != len(SECOND_ORDER_INDICES):
        raise RuntimeError("Second-order block lost rank")
    return estimate


def fit_observation(observation_id, time):
    ah_center = horizon_b_center(time)
    center, center_iterations, final_correction = iterate_curvature_center(
        observation_id, ah_center
    )
    inverse_data = interpolate_inverse_metric(observation_id, center)
    numerical_derivatives = interpolate_radial_derivative(
        observation_id, center, inverse_data
    )

    number_of_radii = len(RADII)
    first_history = np.empty(
        (number_of_radii, NUMBER_OF_ALTERNATIONS + 1, 13)
    )
    second_history = np.empty(
        (number_of_radii, NUMBER_OF_ALTERNATIONS, len(SECOND_ORDER_INDICES))
    )
    residual_history = np.empty(
        (number_of_radii, NUMBER_OF_ALTERNATIONS)
    )

    for radius_index, radius in enumerate(RADII):
        numerical_metric = unpack_inverse_metric(
            inverse_data[:, radius_index]
        )
        raw_value = subtract_metric(
            numerical_metric, BACKGROUNDS[radius_index]
        )
        numerical_derivative = metric_tuple_from_matrix(
            numerical_derivatives[radius_index]
        )
        raw_derivative = tuple(
            radius * component
            for component in numerical_derivative
        )
        raw_derivative = subtract_metric(
            raw_derivative, BACKGROUND_DERIVATIVES[radius_index]
        )

        first = fit_first_order(
            raw_value, raw_derivative, radius_index
        )
        first_history[radius_index, 0] = first
        for iteration in range(NUMBER_OF_ALTERNATIONS):
            second = fit_second_order(
                first, raw_value, raw_derivative, radius_index
            )
            second_history[radius_index, iteration] = second

            quadratic_value, quadratic_derivative = (
                quadratic_value_and_derivative(first, radius_index)
            )
            second_full = np.zeros(43)
            second_full[SECOND_ORDER_INDICES] = second
            second_value = metric_linear_combination(
                SECOND_ORDER_COLUMNS[radius_index], second_full
            )
            second_derivative = metric_linear_combination(
                SECOND_ORDER_DERIVATIVES[radius_index], second_full
            )
            first_prediction = metric_linear_combination(
                FIRST_ORDER_COLUMNS[radius_index], first
            )
            first_derivative_prediction = metric_linear_combination(
                FIRST_ORDER_DERIVATIVES[radius_index], first
            )
            value_residual = subtract_metric(
                raw_value,
                first_prediction,
                quadratic_value,
                second_value,
            )
            derivative_residual = subtract_metric(
                raw_derivative,
                first_derivative_prediction,
                quadratic_derivative,
                second_derivative,
            )
            residual_history[radius_index, iteration] = np.linalg.norm(
                np.concatenate(
                    (
                        metric_values_by_point(value_residual).ravel(),
                        metric_values_by_point(derivative_residual).ravel(),
                    )
                )
            ) / np.linalg.norm(
                np.concatenate(
                    (
                        metric_values_by_point(raw_value).ravel(),
                        metric_values_by_point(raw_derivative).ravel(),
                    )
                )
            )

            correction_value = add_metrics(
                quadratic_value, second_value
            )
            correction_derivative = add_metrics(
                quadratic_derivative, second_derivative
            )
            first = fit_first_order(
                raw_value,
                raw_derivative,
                radius_index,
                correction_value,
                correction_derivative,
            )
            first_history[radius_index, iteration + 1] = first

    print(
        f"t={time:.1f}: center iterations={center_iterations}, "
        f"|last dz|={np.linalg.norm(final_correction):.2e}, "
        f"center={center}"
    )
    return (
        ah_center,
        center,
        first_history,
        second_history,
        residual_history,
    )


def savgol_time_derivative(values, times):
    """Cubic local-polynomial derivative on the uniform 5M grid."""
    values = np.asarray(values)
    delta = float(np.median(np.diff(times)))
    return savgol_filter(
        values,
        window_length=7,
        polyorder=3,
        deriv=1,
        delta=delta,
        axis=0,
        mode="interp",
    )


def cumulative_trapezoid_samples(values, times):
    increments = (
        0.5
        * (values[1:] + values[:-1])
        * np.diff(times).reshape((-1,) + (1,) * (values.ndim - 1))
    )
    result = np.zeros_like(values)
    result[1:] = np.cumsum(increments, axis=0)
    return result


def temporal_diagnostics(times, first, second):
    """Return local and integrated ODE diagnostics.

    first has shape (time, radius, stage, parameter), where stage 0 through
    NUMBER_OF_ALTERNATIONS-1 are paired with the second-order solve at the
    same stage.  second has shape (time, radius, stage, seven).
    """
    number_of_stages = second.shape[2]
    number_of_radii = first.shape[1]
    local_clock = np.empty((number_of_stages, number_of_radii))
    local_clock_uncorrected = np.empty_like(local_clock)
    local_strain = np.empty_like(local_clock)
    local_strain_uncorrected = np.empty_like(local_clock)
    integral_clock = np.empty_like(local_clock)
    integral_strain = np.empty_like(local_clock)

    predicted_acceleration = np.empty(
        (len(times), number_of_radii, number_of_stages)
    )
    fitted_acceleration = np.empty_like(predicted_acceleration)
    predicted_strain_rate = np.empty(
        (len(times), number_of_radii, number_of_stages, 3, 3)
    )
    fitted_strain_rate = np.empty_like(predicted_strain_rate)

    for stage in range(number_of_stages):
        for radius_index in range(number_of_radii):
            qdot = first[:, radius_index, stage, 0]
            strain = np.asarray(
                [
                    symmetric_matrix_from_six(values)
                    for values in first[:, radius_index, stage, 7:13]
                ]
            )
            acceleration = second[:, radius_index, stage, -1]
            strain_rate = np.einsum(
                "ta,aij->tij",
                second[:, radius_index, stage, :6],
                SYMMETRIC_BASIS,
            )

            dqdot_dT = savgol_time_derivative(qdot, times)
            dstrain_dT = savgol_time_derivative(strain, times)
            clock_prediction = (1.0 + qdot) * dqdot_dT
            strain_prediction = (
                SMALL_BLACK_HOLE_MASS
                * (1.0 + qdot)[:, None, None]
                * dstrain_dT
            )
            predicted_acceleration[:, radius_index, stage] = (
                clock_prediction
            )
            fitted_acceleration[:, radius_index, stage] = acceleration
            predicted_strain_rate[
                :, radius_index, stage
            ] = strain_prediction
            fitted_strain_rate[:, radius_index, stage] = strain_rate

            interior = slice(2, -2)
            local_clock[stage, radius_index] = np.linalg.norm(
                acceleration[interior] - clock_prediction[interior]
            ) / max(
                np.linalg.norm(acceleration[interior]),
                np.linalg.norm(clock_prediction[interior]),
                1.0e-300,
            )
            local_clock_uncorrected[stage, radius_index] = np.linalg.norm(
                acceleration[interior] - dqdot_dT[interior]
            ) / max(
                np.linalg.norm(acceleration[interior]),
                np.linalg.norm(dqdot_dT[interior]),
                1.0e-300,
            )
            local_strain[stage, radius_index] = np.linalg.norm(
                strain_rate[interior] - strain_prediction[interior]
            ) / max(
                np.linalg.norm(strain_rate[interior]),
                np.linalg.norm(strain_prediction[interior]),
                1.0e-300,
            )
            uncorrected_strain_prediction = (
                SMALL_BLACK_HOLE_MASS * dstrain_dT
            )
            local_strain_uncorrected[stage, radius_index] = np.linalg.norm(
                strain_rate[interior]
                - uncorrected_strain_prediction[interior]
            ) / max(
                np.linalg.norm(strain_rate[interior]),
                np.linalg.norm(
                    uncorrected_strain_prediction[interior]
                ),
                1.0e-300,
            )

            # dt/dT=1/(1+qdot).  Integrate fitted harmonic-time rates over
            # simulation time and compare changes, removing the arbitrary
            # initial value.
            clock_integrand = acceleration / (1.0 + qdot)
            reconstructed_qdot = (
                qdot[0]
                + cumulative_trapezoid_samples(
                    clock_integrand[:, None], times
                )[:, 0]
            )
            strain_integrand = (
                strain_rate
                / (
                    SMALL_BLACK_HOLE_MASS
                    * (1.0 + qdot)[:, None, None]
                )
            )
            reconstructed_strain = (
                strain[0]
                + cumulative_trapezoid_samples(
                    strain_integrand, times
                )
            )
            integral_clock[stage, radius_index] = np.linalg.norm(
                (qdot - qdot[0])
                - (reconstructed_qdot - reconstructed_qdot[0])
            ) / max(
                np.linalg.norm(qdot - qdot[0]),
                np.linalg.norm(
                    reconstructed_qdot - reconstructed_qdot[0]
                ),
                1.0e-300,
            )
            integral_strain[stage, radius_index] = np.linalg.norm(
                (strain - strain[0])
                - (reconstructed_strain - reconstructed_strain[0])
            ) / max(
                np.linalg.norm(strain - strain[0]),
                np.linalg.norm(
                    reconstructed_strain - reconstructed_strain[0]
                ),
                1.0e-300,
            )
    return {
        "local_clock": local_clock,
        "local_clock_uncorrected": local_clock_uncorrected,
        "local_strain": local_strain,
        "local_strain_uncorrected": local_strain_uncorrected,
        "integral_clock": integral_clock,
        "integral_strain": integral_strain,
        "predicted_acceleration": predicted_acceleration,
        "fitted_acceleration": fitted_acceleration,
        "predicted_strain_rate": predicted_strain_rate,
        "fitted_strain_rate": fitted_strain_rate,
    }


def make_plots(times, first, diagnostics):
    stages = np.arange(NUMBER_OF_ALTERNATIONS)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, NUMBER_OF_ALTERNATIONS))

    amplitude_ratio_clock = np.empty(
        (NUMBER_OF_ALTERNATIONS, len(RADII)))
    amplitude_ratio_strain = np.empty_like(amplitude_ratio_clock)
    alignment_clock = np.empty_like(amplitude_ratio_clock)
    alignment_strain = np.empty_like(amplitude_ratio_clock)
    interior = slice(2, -2)
    for stage in stages:
        for radius_index in range(len(RADII)):
            fitted_clock = diagnostics["fitted_acceleration"][
                interior, radius_index, stage]
            predicted_clock = diagnostics["predicted_acceleration"][
                interior, radius_index, stage]
            fitted_strain = diagnostics["fitted_strain_rate"][
                interior, radius_index, stage]
            predicted_strain = diagnostics["predicted_strain_rate"][
                interior, radius_index, stage]
            amplitude_ratio_clock[stage, radius_index] = (
                np.linalg.norm(fitted_clock)
                / np.linalg.norm(predicted_clock))
            amplitude_ratio_strain[stage, radius_index] = (
                np.linalg.norm(fitted_strain)
                / np.linalg.norm(predicted_strain))
            alignment_clock[stage, radius_index] = (
                np.sum(fitted_clock * predicted_clock)
                / (
                    np.linalg.norm(fitted_clock)
                    * np.linalg.norm(predicted_clock)
                )
            )
            alignment_strain[stage, radius_index] = (
                np.sum(fitted_strain * predicted_strain)
                / (
                    np.linalg.norm(fitted_strain)
                    * np.linalg.norm(predicted_strain)
                )
            )

    fig, axes = plt.subplots(
        2, 2, figsize=(13.5, 9.0), constrained_layout=True)
    for stage, color in zip(stages, colors):
        label = f"alternation {stage}"
        axes[0, 0].loglog(
            RADII / SMALL_BLACK_HOLE_MASS,
            amplitude_ratio_clock[stage],
            "o-",
            color=color,
            label=label,
        )
        axes[0, 1].loglog(
            RADII / SMALL_BLACK_HOLE_MASS,
            amplitude_ratio_strain[stage],
            "o-",
            color=color,
            label=label,
        )
        axes[1, 0].plot(
            RADII / SMALL_BLACK_HOLE_MASS,
            alignment_clock[stage],
            "o-",
            color=color,
            label=label,
        )
        axes[1, 1].plot(
            RADII / SMALL_BLACK_HOLE_MASS,
            alignment_strain[stage],
            "o-",
            color=color,
            label=label,
        )
    axes[0, 0].set_title(
        r"$\|\ddot q^0_{\rm fit}\|/"
        r"\|(1+\dot q^0)d\dot q^0/dT\|$")
    axes[0, 1].set_title(
        r"$\|\dot L_{\rm fit}\|/"
        r"\|M_B(1+\dot q^0)d\Lambda/dT\|$")
    axes[1, 0].set_title("Clock-rate alignment cosine")
    axes[1, 1].set_title("Strain-rate Frobenius alignment cosine")
    for axis in axes.flat:
        axis.set_xlabel(r"$R/M_B$")
    axes[0, 0].set_ylabel("amplitude ratio")
    axes[0, 1].set_ylabel("amplitude ratio")
    axes[1, 0].set_ylabel("cosine")
    axes[1, 1].set_ylabel("cosine")
    axes[1, 0].set_ylim(-1.05, 1.05)
    axes[1, 1].set_ylim(-1.05, 1.05)
    axes[1, 0].axhline(0.0, color="0.5", lw=1)
    axes[1, 1].axhline(0.0, color="0.5", lw=1)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle(
        "Temporal ODE consistency across worldtube radius "
        "(central six observations)")
    fig.savefig(
        OUTPUT_DIRECTORY / "q8_temporal_ode_consistency_by_radius.png",
        dpi=180,
    )
    fig.savefig(
        OUTPUT_DIRECTORY / "q8_temporal_ode_consistency_by_radius.pdf"
    )
    plt.close(fig)

    best_stage, best_radius = np.unravel_index(
        np.argmin(
            diagnostics["local_clock"]
            + diagnostics["local_strain"]
        ),
        diagnostics["local_clock"].shape,
    )
    qdot = first[:, best_radius, best_stage, 0]
    fitted_acceleration = diagnostics["fitted_acceleration"][
        :, best_radius, best_stage
    ]
    predicted_acceleration = diagnostics["predicted_acceleration"][
        :, best_radius, best_stage
    ]
    fitted_strain_rate = diagnostics["fitted_strain_rate"][
        :, best_radius, best_stage
    ]
    predicted_strain_rate = diagnostics["predicted_strain_rate"][
        :, best_radius, best_stage
    ]

    fig, axes = plt.subplots(
        2, 2, figsize=(15.0, 9.0), constrained_layout=True)
    axes[0, 0].plot(times, qdot, "o-")
    axes[0, 0].set_title(r"First-order clock rate $\dot q^0$")
    axes[0, 1].semilogy(
        times,
        np.abs(fitted_acceleration),
        "o-",
        label=r"$|\ddot q^0|$ from second-order fit",
    )
    axes[0, 1].semilogy(
        times,
        np.abs(predicted_acceleration),
        "s--",
        label=r"$|(1+\dot q^0)d\dot q^0/dT|$",
    )
    axes[0, 1].set_title("Time-acceleration magnitude")
    axes[0, 1].legend()
    axes[1, 0].semilogy(
        times,
        np.linalg.norm(fitted_strain_rate, axis=(1, 2)),
        "o-",
        label="second-order fit",
    )
    axes[1, 0].semilogy(
        times,
        np.linalg.norm(predicted_strain_rate, axis=(1, 2)),
        "s--",
        label=r"$M_B(1+\dot q^0)d\Lambda/dT$",
    )
    axes[1, 0].set_title("Strain-rate Frobenius norm")
    axes[1, 0].legend()
    axes[1, 1].semilogy(
        times,
        np.linalg.norm(
            fitted_strain_rate - predicted_strain_rate, axis=(1, 2)
        ),
        "o-",
        label="difference",
    )
    axes[1, 1].semilogy(
        times,
        np.linalg.norm(fitted_strain_rate, axis=(1, 2)),
        "--",
        label="fitted norm",
    )
    axes[1, 1].set_title("Pointwise strain-rate discrepancy")
    axes[1, 1].legend()
    for axis in axes.flat:
        axis.set_xlabel(r"simulation time $T/M_{\rm tot}$")
    fig.suptitle(
        f"Best local ODE comparison: R={RADII[best_radius]:.2f} "
        f"({RADII[best_radius]/SMALL_BLACK_HOLE_MASS:.2f} M_B), "
        f"alternation {best_stage}"
    )
    fig.savefig(
        OUTPUT_DIRECTORY / "q8_temporal_ode_best_timeseries.png", dpi=180
    )
    fig.savefig(
        OUTPUT_DIRECTORY / "q8_temporal_ode_best_timeseries.pdf"
    )
    plt.close(fig)
    return best_stage, best_radius


def main():
    observations = volume_observations()
    observations = [
        item for item in observations if 949.0 <= item[0] <= 996.0
    ]
    times = np.asarray([item[0] for item in observations])
    print(
        f"Analyzing {len(observations)} observations: "
        f"{times[0]:.1f} through {times[-1]:.1f}"
    )
    print(
        f"Spherepack l_max={SPHEREPACK_L_MAX}, "
        f"{len(DIRECTIONS)} angular points, radii={RADII}"
    )

    ah_centers = []
    curvature_centers = []
    first_histories = []
    second_histories = []
    residual_histories = []
    for time, observation_id in observations:
        (
            ah_center,
            curvature_center,
            first_history,
            second_history,
            residual_history,
        ) = fit_observation(observation_id, time)
        ah_centers.append(ah_center)
        curvature_centers.append(curvature_center)
        first_histories.append(first_history)
        second_histories.append(second_history)
        residual_histories.append(residual_history)

    first = np.asarray(first_histories)
    second = np.asarray(second_histories)
    residual = np.asarray(residual_histories)
    diagnostics = temporal_diagnostics(times, first, second)
    best_stage, best_radius = make_plots(times, first, diagnostics)

    np.savez(
        OUTPUT_DIRECTORY / "q8_temporal_ode_consistency.npz",
        times=times,
        radii=RADII,
        small_black_hole_mass=SMALL_BLACK_HOLE_MASS,
        ah_centers=np.asarray(ah_centers),
        curvature_centers=np.asarray(curvature_centers),
        first_order_history=first,
        second_order_history=second,
        surface_residual_history=residual,
        **diagnostics,
    )

    print("\nBest combined local mismatch:")
    print(
        f"  radius={RADII[best_radius]:.2f} "
        f"({RADII[best_radius]/SMALL_BLACK_HOLE_MASS:.2f} M_B), "
        f"alternation={best_stage}"
    )
    print(
        f"  clock local={diagnostics['local_clock'][best_stage,best_radius]:.6g}"
    )
    print(
        f"  strain local={diagnostics['local_strain'][best_stage,best_radius]:.6g}"
    )
    print(
        f"  clock integral={diagnostics['integral_clock'][best_stage,best_radius]:.6g}"
    )
    print(
        f"  strain integral={diagnostics['integral_strain'][best_stage,best_radius]:.6g}"
    )
    print("\nTime-coordinate correction effect at this point:")
    print(
        "  clock: "
        f"{diagnostics['local_clock_uncorrected'][best_stage,best_radius]:.6g}"
        " -> "
        f"{diagnostics['local_clock'][best_stage,best_radius]:.6g}"
    )
    print(
        "  strain: "
        f"{diagnostics['local_strain_uncorrected'][best_stage,best_radius]:.6g}"
        " -> "
        f"{diagnostics['local_strain'][best_stage,best_radius]:.6g}"
    )
    print("\nMinimum mismatch over radius for each alternation:")
    for stage in range(NUMBER_OF_ALTERNATIONS):
        clock_index = int(np.argmin(diagnostics["local_clock"][stage]))
        strain_index = int(np.argmin(diagnostics["local_strain"][stage]))
        print(
            f"  {stage}: clock={diagnostics['local_clock'][stage,clock_index]:.4g}"
            f" at R={RADII[clock_index]:.2f}; "
            f"strain={diagnostics['local_strain'][stage,strain_index]:.4g}"
            f" at R={RADII[strain_index]:.2f}"
        )


if __name__ == "__main__":
    main()
