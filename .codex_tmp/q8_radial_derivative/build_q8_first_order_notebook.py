#!/usr/bin/env python
"""Build the focused q=8 first-order matching notebook."""

from __future__ import annotations

from pathlib import Path

import nbformat


OUTPUT_DIRECTORY = Path(__file__).resolve().parent


def markdown(source: str):
    return nbformat.v4.new_markdown_cell(source.strip() + "\n")


def code(source: str):
    return nbformat.v4.new_code_cell(source.strip() + "\n")


def main() -> None:
    notebook = nbformat.v4.new_notebook()
    notebook.metadata.kernelspec = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook.metadata.language_info = {"name": "python"}
    notebook.cells = [
        markdown(
            r"""
# q=8 first-order matching and held-out closure

This notebook analyzes one full-resolution q=8 volume-data slice at
$t=950M_{\rm tot}$.  It addresses one question:

> Do different minimal first-order metric/STF blocks recover compatible
> coordinate-map coefficients and predict the blocks that were not fitted?

The small-black-hole mass is fixed to its initial-data value
$M_B=1/9$.  A single common center is determined first from the dipole of the
Gauss--Bonnet scalar, which equals the Kretschmann scalar in vacuum.  That
center is then held fixed for every first-order estimator.

The response functions are imported from the audited q=4 package.  In
particular, the spatial-strain response is evaluated at fixed final/inertial
points and includes transport of the spatially varying Schwarzschild
background.
"""
        ),
        code(
            r"""
import os
import sys
from itertools import combinations
from pathlib import Path

SPECTRE_BUILD = Path(os.environ.get(
    "SPECTRE_BUILD", "/Users/niko/caltech/spectre/build-Release"))
sys.path.insert(0, str(SPECTRE_BUILD / "bin/python"))
sys.path.insert(
    0, str(SPECTRE_BUILD / "lib/python3.13/site-packages"))

AUDIT_PACKAGE_ROOT = Path(os.environ.get(
    "GAUGE_MATCHING_AUDIT",
    "/Users/niko/caltech/worldtube-gravity/gauge-matching-audit",
))
sys.path.insert(0, str(AUDIT_PACKAGE_ROOT / "python"))

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/q8-matching-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/q8-matching-cache")

import matplotlib
if "ipykernel" in sys.modules:
    matplotlib.use("module://matplotlib_inline.backend_inline")
import matplotlib.pyplot as plt
import numpy as np

import spectre.IO.H5 as spectre_h5
from spectre.DataStructures import DataVector
from spectre.DataStructures.Tensor import Frame, tnsr
from spectre.IO.Exporter import ObservationId, interpolate_to_points
from spectre.SphericalHarmonics import Spherepack, SpherepackIterator
from worldtube_matching.responses import (
    first_order_response_columns,
    inverse_harmonic_schwarzschild,
)

Q8_DIRECTORY = Path(os.environ.get(
    "Q8DIR", "/Users/niko/caltech/simulations/harmonic-sims/q8"))
VOLUME_FILES = [str(Q8_DIRECTORY / "subset.h5")]
REDUCTIONS_FILE = str(Q8_DIRECTORY / "BbhReductions.h5")
VOLUME_SUBFILE = "MetricAndScalars"
OUTPUT_DIRECTORY = Path(os.environ.get("Q8_OUTPUT_DIR", Q8_DIRECTORY))
OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

REQUESTED_TIME = 950.0
SMALL_BLACK_HOLE_MASS = 1.0 / 9.0
SPHEREPACK_L_MAX = 16
CENTER_FINDING_RADIUS = 0.20
RADIUS_TEST_VALUES = np.array([
    0.20, 0.24, 0.28, 0.34, 0.40,
    0.48, 0.58, 0.68, 0.80,
    0.96, 1.16, 1.36, 1.60,
])
RADIUS_OVER_SMALL_MASS = RADIUS_TEST_VALUES / SMALL_BLACK_HOLE_MASS

plt.rcParams.update({
    "figure.figsize": (10.0, 4.8),
    "figure.max_open_warning": 0,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "font.size": 11,
})
"""
        ),
        markdown(
            r"""
## 1. Angular projection conventions

SpECTRE's `Spherepack` performs the physical-to-spectral transform.  We retain
ordinary real spherical-harmonic modes and, where useful for the center,
convert through $\ell=2$ to the Cartesian polynomial basis

$$
1,\quad n_x,\ n_y,\ n_z,\quad
n_xn_y,\ n_xn_z,\ n_yn_z,\quad
n_x^2-n_y^2,\quad 2n_z^2-n_x^2-n_y^2.
$$

The first-order closure calculation itself keeps complete real-$Y_{\ell m}$
blocks.  It therefore does not depend on manually identifying individual
Cartesian STF components.
"""
        ),
        code(
            r"""
def spherepack_grid(ell_max):
    spherepack = Spherepack(ell_max, ell_max)
    theta, phi = map(np.asarray, spherepack.theta_phi_points)
    directions = np.column_stack((
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ))
    return spherepack, directions


def spherical_harmonic_modes(ell_max):
    return tuple(
        (ell, m)
        for ell in range(ell_max + 1)
        for m in range(-ell, ell + 1)
    )


def spherical_harmonic_coefficients(
        values, spherepack, ell_max):
    spectral = np.asarray(
        spherepack.phys_to_spec(DataVector(np.asarray(values))))
    iterator = SpherepackIterator(
        SPHEREPACK_L_MAX, SPHEREPACK_L_MAX)
    indices = []
    for ell, m in spherical_harmonic_modes(ell_max):
        iterator.set(ell, m)
        indices.append(iterator())
    return spectral[indices]


def ylm_to_cartesian_matrix(spherepack, directions):
    nx, ny, nz = directions.T
    cartesian_basis = np.column_stack((
        np.ones_like(nx), nx, ny, nz,
        nx * ny, nx * nz, ny * nz,
        nx * nx - ny * ny,
        2.0 * nz * nz - nx * nx - ny * ny,
    ))
    cartesian_to_ylm = np.column_stack([
        spherical_harmonic_coefficients(
            cartesian_basis[:, column], spherepack, ell_max=2)
        for column in range(cartesian_basis.shape[1])
    ])
    return np.linalg.inv(cartesian_to_ylm)


def cartesian_multipoles(values):
    ylm = spherical_harmonic_coefficients(
        values, SPHEREPACK, ell_max=2)
    coefficients = YLM_TO_CARTESIAN @ ylm
    monopole = coefficients[0]
    dipole = coefficients[1:4].copy()
    quadrupole = np.zeros((3, 3))
    quadrupole[0, 1] = quadrupole[1, 0] = coefficients[4]
    quadrupole[0, 2] = quadrupole[2, 0] = coefficients[5]
    quadrupole[1, 2] = quadrupole[2, 1] = coefficients[6]
    quadrupole[0, 0] += coefficients[7] - coefficients[8]
    quadrupole[1, 1] -= coefficients[7] + coefficients[8]
    quadrupole[2, 2] += 2.0 * coefficients[8]
    return monopole, dipole, quadrupole


SPHEREPACK, DIRECTIONS = spherepack_grid(SPHEREPACK_L_MAX)
YLM_TO_CARTESIAN = ylm_to_cartesian_matrix(
    SPHEREPACK, DIRECTIONS)

print(
    f"Spherepack l_max={SPHEREPACK_L_MAX}: "
    f"{len(DIRECTIONS)} collocation points")
"""
        ),
        markdown(
            r"""
## 2. Select the first downloaded slice

The volume file contains ten observations from $950M_{\rm tot}$ through
$995M_{\rm tot}$ in steps of $5M_{\rm tot}$.  We select the observation nearest
$950M_{\rm tot}$ and use the inertial apparent-horizon center only as the
initial guess for the curvature-center iteration.
"""
        ),
        code(
            r"""
def volume_observations():
    with spectre_h5.H5File(VOLUME_FILES[0], "r") as h5file:
        volume = h5file.get_vol("/" + VOLUME_SUBFILE)
        return sorted([
            (
                volume.get_observation_value(observation_id),
                observation_id,
            )
            for observation_id in volume.list_observation_ids()
        ])


def horizon_b_center(time):
    with spectre_h5.H5File(REDUCTIONS_FILE, "r") as h5file:
        centers = np.asarray(h5file.get_dat(
            "/ApparentHorizons/ControlSystemAhB_Centers"
        ).get_data())
    return np.array([
        np.interp(time, centers[:, 0], centers[:, component])
        for component in range(4, 7)
    ])


def horizon_b_mass(time):
    with spectre_h5.H5File(REDUCTIONS_FILE, "r") as h5file:
        horizon = h5file.get_dat("/ObservationAhB")
        legend = list(horizon.get_legend())
        data = np.asarray(horizon.get_data())
    mass_column = legend.index("ChristodoulouMass")
    return np.interp(time, data[:, 0], data[:, mass_column])


OBSERVATIONS = volume_observations()
observation_times = np.array([time for time, _ in OBSERVATIONS])
observation_index = int(np.argmin(
    np.abs(observation_times - REQUESTED_TIME)))
TEST_TIME, TEST_OBSERVATION_ID = OBSERVATIONS[observation_index]
AH_B_CENTER = horizon_b_center(TEST_TIME)
AH_B_MASS = horizon_b_mass(TEST_TIME)

print(f"selected time:                 {TEST_TIME:.12f}")
print(f"observation id:                {TEST_OBSERVATION_ID}")
print("AH-B inertial center:          ", AH_B_CENTER)
print(f"fixed initial-data mass:       {SMALL_BLACK_HOLE_MASS:.12f}")
print(f"measured AH-B mass at slice:   {AH_B_MASS:.12f}")
print(
    "worldtube radii R/M_B:       ",
    np.array2string(RADIUS_OVER_SMALL_MASS, precision=3),
)
"""
        ),
        markdown(
            r"""
## 3. One common Kretschmann-dipole center

For Schwarzschild,

$$
\mathcal K=\frac{48M_B^2}{r^6},
\qquad
\rho=r-M_B,
$$

where $r$ is the areal radius and $\rho$ the harmonic radius.  We therefore
infer a Schwarzschild harmonic radius from the numerical Gauss--Bonnet scalar,

$$
\rho_{\mathcal K}
=\left(\frac{48M_B^2}{\mathcal K}\right)^{1/6}-M_B.
$$

On a coordinate sphere of radius $R$ centered at the current guess, define

$$
F(n^i)=R^2-\rho_{\mathcal K}^2.
$$

If the physical center is displaced by $\delta z^i$, then to leading order
$F=2R\,\delta z_i n^i+\cdots$.  Thus the Cartesian dipole gives the Newton
correction

$$
\delta z_i=\frac{F_i^{\ell=1}}{2R}.
$$

We resample after every correction at $R=0.2M_{\rm tot}$ and use the
converged result as one common center for all fitting radii.
"""
        ),
        code(
            r"""
def interpolate_components(
        observation_id, points, component_names):
    points = np.asarray(points)
    target_points = tnsr.I[DataVector, 3, Frame.Inertial](
        points.T.copy())
    return np.asarray(interpolate_to_points(
        VOLUME_FILES,
        subfile_name=VOLUME_SUBFILE,
        observation=ObservationId(observation_id),
        tensor_components=list(component_names),
        target_points=target_points,
        error_on_missing_points=True,
    ))


def interpolate_gauss_bonnet(
        observation_id, centers, radii):
    centers = np.asarray(centers)
    radii = np.broadcast_to(np.asarray(radii), (len(centers),))
    points = np.concatenate([
        center[None, :] + radius * DIRECTIONS
        for center, radius in zip(centers, radii)
    ])
    values = interpolate_components(
        observation_id, points, ["GaussBonnetScalar"])[0]
    return values.reshape(len(centers), len(DIRECTIONS))


def kretschmann_harmonic_radius(
        kretschmann, mass=SMALL_BLACK_HOLE_MASS):
    areal_radius = (
        48.0 * mass**2 / np.asarray(kretschmann)
    ) ** (1.0 / 6.0)
    return areal_radius - mass


def curvature_center_correction(kretschmann, radius):
    inferred_radius = kretschmann_harmonic_radius(kretschmann)
    radial_residual = radius**2 - inferred_radius**2
    _, dipole, _ = cartesian_multipoles(radial_residual)
    return dipole / (2.0 * radius)


def iterate_curvature_center(
        initial_center, radius, maximum_iterations=12,
        tolerance=1.0e-12):
    center = np.array(initial_center, dtype=float)
    history = []
    for _ in range(maximum_iterations):
        kretschmann = interpolate_gauss_bonnet(
            TEST_OBSERVATION_ID, [center], [radius])[0]
        correction = curvature_center_correction(
            kretschmann, radius)
        center += correction
        history.append((center.copy(), correction.copy()))
        if np.linalg.norm(correction) < tolerance:
            break
    return center, history


CURVATURE_CENTER, CENTER_ITERATION_HISTORY = (
    iterate_curvature_center(
        AH_B_CENTER, CENTER_FINDING_RADIUS))
center_positions = np.asarray([
    entry[0] for entry in CENTER_ITERATION_HISTORY])
center_corrections = np.asarray([
    entry[1] for entry in CENTER_ITERATION_HISTORY])

# Diagnose the remaining one-step correction on every fitting radius without
# allowing the metric fit itself to change centers.
kretschmann_by_radius = interpolate_gauss_bonnet(
    TEST_OBSERVATION_ID,
    np.repeat(CURVATURE_CENTER[None, :], len(RADIUS_TEST_VALUES), axis=0),
    RADIUS_TEST_VALUES,
)
CENTER_CORRECTIONS_BY_RADIUS = np.asarray([
    curvature_center_correction(values, radius)
    for values, radius in zip(
        kretschmann_by_radius, RADIUS_TEST_VALUES)
])

print("AH center:             ", AH_B_CENTER)
print("curvature center:      ", CURVATURE_CENTER)
print("curvature - AH center: ", CURVATURE_CENTER - AH_B_CENTER)
print(
    "last center correction:",
    center_corrections[-1],
    "norm =", np.linalg.norm(center_corrections[-1]),
)

fig_center, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
axes[0].semilogy(
    np.arange(1, len(center_corrections) + 1),
    np.linalg.norm(center_corrections, axis=1),
    "o-",
)
axes[0].set_xlabel("outer center iteration")
axes[0].set_ylabel("correction norm")
axes[0].set_title(
    rf"Kretschmann center at $R={CENTER_FINDING_RADIUS:g}M_{{\rm tot}}$")

axes[1].loglog(
    RADIUS_OVER_SMALL_MASS,
    np.linalg.norm(CENTER_CORRECTIONS_BY_RADIUS, axis=1),
    "o-",
)
axes[1].set_xlabel(r"$R/M_B$")
axes[1].set_ylabel("one-step correction norm")
axes[1].set_title("Residual dipole about the common center")
fig_center.tight_layout()
fig_center.savefig(
    OUTPUT_DIRECTORY / "q8_t950_curvature_center.pdf")
plt.show()
"""
        ),
        markdown(
            r"""
## 4. First-order response blocks

The thirteen unknown first-order affine-map coefficients are

$$
p_{(1)}
=\left(
\dot q^0,\,
\beta_i,\,
\dot q^i,\,
\sigma_{ij}
\right),
\qquad \sigma_{ij}=\sigma_{ji}.
$$

After subtracting harmonic Schwarzschild at the common curvature center,

$$
d_{\rm block}
=G_{\rm numerical,block}-G_{\rm Schwarzschild,block}
=R_{\rm block}\,p_{(1)}+\text{unmodeled terms}.
$$

Parity and rotational structure separate the response into two independent
sectors:

| Sector | Unknowns | Candidate metric/STF blocks |
|---|---|---|
| Vector | $\beta_i,\dot q^i$ | $G^{TT}_{\ell=1}$, $G^{Ti}_{\ell=0,2}$, $G^{ij}_{\ell=1}$ |
| Clock/strain | $\dot q^0,\sigma_{ij}$ | $G^{TT}_{\ell=0,2}$, $G^{Ti}_{\ell=1,3}$, $G^{ij}_{\ell=0,2,4}$ |

A *minimal estimator* is an inclusion-minimal subset of candidate blocks with
full column rank.  It is fitted using only its selected blocks.  Its held-out
closure residual is

$$
\mathcal C
=\frac{\lVert d_{\rm omitted}-R_{\rm omitted}\hat p\rVert_2}
{\lVert d_{\rm omitted}\rVert_2}.
$$

This tests whether parameters inferred from one set of metric projections
predict independent projections that were not used in the fit.
"""
        ),
        code(
            r"""
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
    (0, 0), (0, 1), (0, 2),
    (1, 1), (1, 2), (2, 2),
)
DATA_SPATIAL_COMPONENTS = (
    (0, 0), (0, 1), (1, 1),
    (0, 2), (1, 2), (2, 2),
)
FIRST_ORDER_PARAMETER_NAMES = (
    "qdot0",
    "beta_x", "beta_y", "beta_z",
    "qdot_x", "qdot_y", "qdot_z",
    "sigma_xx", "sigma_xy", "sigma_xz",
    "sigma_yy", "sigma_yz", "sigma_zz",
)


def interpolate_inverse_metric(
        observation_id, center, radii):
    points = np.concatenate([
        center[None, :] + radius * DIRECTIONS
        for radius in radii
    ])
    values = interpolate_components(
        observation_id,
        points,
        INVERSE_METRIC_COMPONENT_NAMES,
    )
    return values.reshape(
        len(INVERSE_METRIC_COMPONENT_NAMES),
        len(radii),
        len(DIRECTIONS),
    )


def unpack_inverse_metric(field_data):
    gtt = field_data[0]
    gti = np.moveaxis(field_data[1:4], 0, -1)
    gij = np.zeros((len(gtt), 3, 3))
    for values, (i, j) in zip(
            field_data[4:], DATA_SPATIAL_COMPONENTS):
        gij[:, i, j] = values
        gij[:, j, i] = values
    return gtt, gti, gij


def ell_block(values, ell):
    coefficients = spherical_harmonic_coefficients(
        values, SPHEREPACK, ell_max=ell)
    return coefficients[ell**2:(ell + 1)**2]


def metric_stf_block(gtt, gti, gij, block):
    metric_component, ell = block
    if metric_component == "GTT":
        return ell_block(gtt, ell)
    if metric_component == "GTI":
        return np.concatenate([
            ell_block(gti[:, component], ell)
            for component in range(3)
        ])
    if metric_component == "GIJ":
        return np.concatenate([
            ell_block(gij[:, i, j], ell)
            for i, j in SPATIAL_COMPONENTS
        ])
    raise ValueError(f"Unknown metric component {metric_component}")


FIRST_ORDER_BLOCKS = (
    ("GTT", 0), ("GTT", 1), ("GTT", 2),
    ("GTI", 0), ("GTI", 1), ("GTI", 2), ("GTI", 3),
    ("GIJ", 0), ("GIJ", 1), ("GIJ", 2), ("GIJ", 4),
)

FIRST_ORDER_METRIC_DATA = interpolate_inverse_metric(
    TEST_OBSERVATION_ID,
    CURVATURE_CENTER,
    RADIUS_TEST_VALUES,
)

BLOCK_DATA = {block: [] for block in FIRST_ORDER_BLOCKS}
BLOCK_RESPONSE = {block: [] for block in FIRST_ORDER_BLOCKS}

for radius_index, radius in enumerate(RADIUS_TEST_VALUES):
    numerical_metric = unpack_inverse_metric(
        FIRST_ORDER_METRIC_DATA[:, radius_index])
    local_points = radius * DIRECTIONS
    background_metric = inverse_harmonic_schwarzschild(
        local_points, SMALL_BLACK_HOLE_MASS)
    metric_residual = tuple(
        numerical - background
        for numerical, background in zip(
            numerical_metric, background_metric)
    )
    response_columns = first_order_response_columns(
        local_points, SMALL_BLACK_HOLE_MASS)

    for block in FIRST_ORDER_BLOCKS:
        BLOCK_DATA[block].append(
            metric_stf_block(*metric_residual, block))
        BLOCK_RESPONSE[block].append(np.column_stack([
            metric_stf_block(*response, block)
            for response in response_columns
        ]))

for block in FIRST_ORDER_BLOCKS:
    BLOCK_DATA[block] = np.asarray(BLOCK_DATA[block])
    BLOCK_RESPONSE[block] = np.asarray(BLOCK_RESPONSE[block])

print("first-order parameter order:")
print(FIRST_ORDER_PARAMETER_NAMES)
"""
        ),
        code(
            r"""
VECTOR_PARAMETER_INDICES = np.arange(1, 7)
CLOCK_STRAIN_PARAMETER_INDICES = np.r_[
    0, np.arange(7, 13)]

VECTOR_BLOCKS = (
    ("GTT", 1), ("GTI", 0),
    ("GTI", 2), ("GIJ", 1),
)
CLOCK_STRAIN_BLOCKS = (
    ("GTT", 0), ("GTT", 2),
    ("GTI", 1), ("GTI", 3),
    ("GIJ", 0), ("GIJ", 2), ("GIJ", 4),
)


def block_name(block):
    return f"{block[0]} l={block[1]}"


def estimator_label(selected_blocks):
    return " + ".join(map(block_name, selected_blocks))


def stacked_system(
        radius_index, parameter_indices, selected_blocks,
        block_data=BLOCK_DATA, block_response=BLOCK_RESPONSE):
    matrix = np.concatenate([
        block_response[block][radius_index][
            :, parameter_indices]
        for block in selected_blocks
    ], axis=0)
    data = np.concatenate([
        block_data[block][radius_index]
        for block in selected_blocks
    ])
    return matrix, data


def numerical_rank(matrix, relative_tolerance=1.0e-10):
    singular_values = np.linalg.svd(
        matrix, compute_uv=False)
    if len(singular_values) == 0 or singular_values[0] == 0.0:
        return 0
    return np.count_nonzero(
        singular_values
        > relative_tolerance * singular_values[0])


def minimal_full_rank_block_sets(
        parameter_indices, candidate_blocks):
    number_of_unknowns = len(parameter_indices)
    full_rank_sets = []
    for number_of_blocks in range(
            1, len(candidate_blocks) + 1):
        for selected_blocks in combinations(
                candidate_blocks, number_of_blocks):
            full_rank_at_every_radius = all(
                numerical_rank(stacked_system(
                    radius_index,
                    parameter_indices,
                    selected_blocks,
                )[0]) == number_of_unknowns
                for radius_index in range(
                    len(RADIUS_TEST_VALUES))
            )
            if not full_rank_at_every_radius:
                continue
            if any(
                set(previous).issubset(selected_blocks)
                for previous in full_rank_sets
            ):
                continue
            full_rank_sets.append(selected_blocks)
    return tuple(full_rank_sets)


def evaluate_estimators(
        parameter_indices, candidate_blocks,
        estimator_block_sets,
        block_data=BLOCK_DATA, block_response=BLOCK_RESPONSE):
    number_of_estimators = len(estimator_block_sets)
    number_of_radii = len(RADIUS_TEST_VALUES)
    number_of_parameters = len(parameter_indices)

    estimates = np.empty((
        number_of_estimators,
        number_of_radii,
        number_of_parameters,
    ))
    condition_numbers = np.empty((
        number_of_estimators, number_of_radii))
    fitted_block_residuals = np.empty((
        number_of_estimators, number_of_radii))
    held_out_residuals = np.empty((
        number_of_estimators, number_of_radii))
    held_out_rms = np.empty((
        number_of_estimators, number_of_radii))

    for estimator_index, selected_blocks in enumerate(
            estimator_block_sets):
        omitted_blocks = tuple(
            block for block in candidate_blocks
            if block not in selected_blocks)
        for radius_index in range(number_of_radii):
            matrix, data = stacked_system(
                radius_index,
                parameter_indices,
                selected_blocks,
                block_data,
                block_response,
            )
            estimates[
                estimator_index, radius_index
            ], *_ = np.linalg.lstsq(
                matrix, data, rcond=None)
            fitted_residual = (
                data
                - matrix
                @ estimates[estimator_index, radius_index]
            )
            condition_numbers[
                estimator_index, radius_index
            ] = np.linalg.cond(matrix)
            fitted_block_residuals[
                estimator_index, radius_index
            ] = (
                np.linalg.norm(fitted_residual)
                / max(np.linalg.norm(data), 1.0e-300)
            )

            omitted_matrix, omitted_data = stacked_system(
                radius_index,
                parameter_indices,
                omitted_blocks,
                block_data,
                block_response,
            )
            omitted_residual = (
                omitted_data
                - omitted_matrix
                @ estimates[estimator_index, radius_index]
            )
            held_out_residuals[
                estimator_index, radius_index
            ] = (
                np.linalg.norm(omitted_residual)
                / max(np.linalg.norm(omitted_data), 1.0e-300)
            )
            held_out_rms[
                estimator_index, radius_index
            ] = (
                np.linalg.norm(omitted_residual)
                / np.sqrt(len(omitted_residual))
            )

    return {
        "estimates": estimates,
        "condition_numbers": condition_numbers,
        "fitted_block_residuals": fitted_block_residuals,
        "held_out_residuals": held_out_residuals,
        "held_out_rms": held_out_rms,
    }


VECTOR_ESTIMATOR_BLOCKS = minimal_full_rank_block_sets(
    VECTOR_PARAMETER_INDICES, VECTOR_BLOCKS)
CLOCK_STRAIN_ESTIMATOR_BLOCKS = (
    minimal_full_rank_block_sets(
        CLOCK_STRAIN_PARAMETER_INDICES,
        CLOCK_STRAIN_BLOCKS,
    )
)

print("\nVector-sector minimal estimators")
for index, blocks in enumerate(VECTOR_ESTIMATOR_BLOCKS):
    print(f"  V{index + 1}: {estimator_label(blocks)}")
print("\nClock/strain-sector minimal estimators")
for index, blocks in enumerate(
        CLOCK_STRAIN_ESTIMATOR_BLOCKS):
    print(f"  C{index + 1}: {estimator_label(blocks)}")
"""
        ),
        markdown(
            r"""
### Analytic radial sensitivity of the first-order terms

Before perturbing fitted coefficients numerically, we can inspect the response
functions themselves.  Let $\rho$ be harmonic radius, $r=\rho+M_B$, and
define

$$
F=1+\frac{2M_B}{r},\qquad
H=\frac{4M_B^2}{r^2},\qquad
P=\frac{\rho^2}{r^2},\qquad
L=-\frac{M_B^2}{r^2},\qquad
T=G^{TT}_{(0)}.
$$

For a Cartesian component $a$, the clock, boost, and velocity responses are

$$
\begin{aligned}
\partial_{\dot q^0}G^{TT}&=2T,
&
\partial_{\dot q^0}G^{Ti}&=Hn_i,
\\
\partial_{\beta_a}G^{TT}&=2Hn_a,
&
\partial_{\beta_a}G^{Ti}
&=P\delta_{ia}+L n_i n_a,
\\
\partial_{\dot q^a}G^{Ti}
&=-F(1+H)\delta_{ia},
&
\partial_{\dot q^a}G^{ij}
&=H(\delta_{ia}n_j+\delta_{ja}n_i).
\end{aligned}
$$

Writing $\sigma_{nn}=\sigma_{ab}n_an_b$, the strain response at fixed final
harmonic coordinates is

$$
\begin{aligned}
\partial_\sigma G^{TT}
&=-\rho\,\sigma_{nn}T',
\\
\partial_\sigma G^{Ti}
&=\sigma_{nn}(H-\rho H')n_i,
\\
\partial_\sigma G^{ij}
&=2P\sigma_{ij}
-\rho\,\sigma_{nn}P'\delta_{ij}
\sigma_{nn}(2L-\rho L')n_i n_j.
\end{aligned}
$$

These formulas make an important point: every denominator is a power of
$r=\rho+M_B$, not $\rho$.  In these horizon-penetrating harmonic coordinates,
the first-order response columns do **not** contain a pole at
$\rho=M_B$ or as $\rho$ decreases toward zero.  Indeed,

$$
F\to3,\quad H\to4,\quad
P=O((\rho/M_B)^2),\quad L\to-1
\qquad(\rho\to0).
$$

The clock, boost, and velocity responses therefore remain finite.  Most
strain responses remain finite as well, but the strain contribution to
$G^{TT}$ vanishes linearly:

$$
\partial_\sigma G^{TT}=O(\rho/M_B).
$$

The expected block behavior is:

| Unknown group | Nonzero STF blocks | Small-$\rho$ behavior | Large-$\rho$ behavior and likely weakness |
|---|---|---|---|
| $\dot q^0$ | $G^{TT}_{0}$, $G^{Ti}_{1}$ | both finite | $G^{TT}_{0}=O(1)$; $G^{Ti}_{1}=O((M_B/\rho)^2)$ |
| $\beta_i$ | $G^{TT}_{1}$, $G^{Ti}_{0,2}$ | all finite | $G^{TT}_{1}$ and $G^{Ti}_{2}$ decay as $(M_B/\rho)^2$; only the $G^{Ti}_{0}$ combination remains $O(1)$ |
| $\dot q^i$ | $G^{Ti}_{0}$, $G^{ij}_{1}$ | both finite | $G^{Ti}_{0}=O(1)$ while $G^{ij}_{1}=O((M_B/\rho)^2)$ |
| $\sigma_{ij}$ | $G^{TT}_{0,2}$, $G^{Ti}_{1,3}$, $G^{ij}_{0,2,4}$ | $G^{TT}_{0,2}$ lose strain sensitivity as $O(\rho/M_B)$; the other blocks remain finite | $G^{ij}_{0}=O(1)$, $G^{TT}_{0,2}$ and $G^{ij}_{2}$ decay roughly as $M_B/\rho$, and $G^{Ti}_{1,3}$ and $G^{ij}_{4}$ decay as $(M_B/\rho)^2$ |

This predicts no generic small-radius explosion of the first-order metric.
It does predict loss of *absolute* strain sensitivity in
$G^{TT}_{\ell=0,2}$ at sufficiently small radius.  At large radius it
predicts a more general loss of information: boost and velocity become
difficult to separate, and clock information in $G^{Ti}_{1}$ becomes weak.

Consequently, a fractional closure error for a radially weak block can become
large even when its absolute residual is small.  The numerical analysis below
therefore keeps separate:

1. response-column strength;
2. response-matrix singular values;
3. fitted-parameter ambiguity;
4. the resulting error in the held-out interior metric.
"""
        ),
        markdown(
            r"""
### Synthetic round-trip test

Before using simulation data, generate every block from one known set of all
thirteen coefficients and pass those data through exactly the same rank
selection, fitting, and held-out prediction code.  Every minimal estimator
must recover its coefficients and all omitted blocks to roundoff.
"""
        ),
        code(
            r"""
SYNTHETIC_PARAMETERS = np.array([
    0.037,
    -0.12, 0.08, 0.025,
    -0.17, 0.045, -0.015,
    0.031, -0.012, 0.009,
    -0.021, 0.014, 0.006,
])
SYNTHETIC_BLOCK_DATA = {
    block: np.einsum(
        "rkp,p->rk",
        BLOCK_RESPONSE[block],
        SYNTHETIC_PARAMETERS,
    )
    for block in FIRST_ORDER_BLOCKS
}

synthetic_vector = evaluate_estimators(
    VECTOR_PARAMETER_INDICES,
    VECTOR_BLOCKS,
    VECTOR_ESTIMATOR_BLOCKS,
    block_data=SYNTHETIC_BLOCK_DATA,
)
synthetic_clock = evaluate_estimators(
    CLOCK_STRAIN_PARAMETER_INDICES,
    CLOCK_STRAIN_BLOCKS,
    CLOCK_STRAIN_ESTIMATOR_BLOCKS,
    block_data=SYNTHETIC_BLOCK_DATA,
)

vector_truth = SYNTHETIC_PARAMETERS[
    VECTOR_PARAMETER_INDICES]
clock_truth = SYNTHETIC_PARAMETERS[
    CLOCK_STRAIN_PARAMETER_INDICES]
maximum_coefficient_error = max(
    np.max(np.abs(
        synthetic_vector["estimates"]
        - vector_truth[None, None, :]
    )),
    np.max(np.abs(
        synthetic_clock["estimates"]
        - clock_truth[None, None, :]
    )),
)
maximum_closure_error = max(
    synthetic_vector["held_out_residuals"].max(),
    synthetic_clock["held_out_residuals"].max(),
)
assert maximum_coefficient_error < 1.0e-10
assert maximum_closure_error < 1.0e-10

print(
    "maximum absolute coefficient error:",
    f"{maximum_coefficient_error:.3e}",
)
print(
    "maximum fractional closure residual:",
    f"{maximum_closure_error:.3e}",
)
"""
        ),
        markdown(
            r"""
## 5. Numerical minimal-estimator comparison

No estimator is selected in advance.  The following cells show:

1. fractional held-out closure;
2. absolute held-out RMS residual;
3. residual on the blocks that were fitted;
4. condition number of each minimal response matrix;
5. the recovered coefficient groups as functions of worldtube radius.

All radius scans use $R/M_B$ on logarithmic axes.
"""
        ),
        code(
            r"""
VECTOR_RESULTS = evaluate_estimators(
    VECTOR_PARAMETER_INDICES,
    VECTOR_BLOCKS,
    VECTOR_ESTIMATOR_BLOCKS,
)
CLOCK_STRAIN_RESULTS = evaluate_estimators(
    CLOCK_STRAIN_PARAMETER_INDICES,
    CLOCK_STRAIN_BLOCKS,
    CLOCK_STRAIN_ESTIMATOR_BLOCKS,
)

VECTOR_LABELS = tuple(
    f"V{index + 1}: {estimator_label(blocks)}"
    for index, blocks in enumerate(VECTOR_ESTIMATOR_BLOCKS)
)
CLOCK_STRAIN_LABELS = tuple(
    f"C{index + 1}: {estimator_label(blocks)}"
    for index, blocks in enumerate(
        CLOCK_STRAIN_ESTIMATOR_BLOCKS)
)

fig_closure, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.0))
for estimator_index, label in enumerate(VECTOR_LABELS):
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_RESULTS["held_out_residuals"][
            estimator_index],
        "o-",
        label=label,
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_RESULTS["held_out_rms"][
            estimator_index],
        "o-",
        label=label,
    )
for estimator_index, label in enumerate(
        CLOCK_STRAIN_LABELS):
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_STRAIN_RESULTS["held_out_residuals"][
            estimator_index],
        "s-",
        label=label,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_STRAIN_RESULTS["held_out_rms"][
            estimator_index],
        "s-",
        label=label,
    )

axes[0, 0].set_title("Vector-sector fractional closure")
axes[0, 1].set_title("Clock/strain fractional closure")
axes[1, 0].set_title("Vector-sector absolute closure")
axes[1, 1].set_title("Clock/strain absolute closure")
for axis in axes[0]:
    axis.set_ylabel("held-out residual / held-out data")
for axis in axes[1]:
    axis.set_ylabel("held-out residual RMS")
for axis in axes.flat:
    axis.set_xlabel(r"$R/M_B$")
    axis.legend(fontsize=7)
fig_closure.suptitle(
    rf"q=8 first-order held-out closure at "
    rf"$t={TEST_TIME:g}M_{{\rm tot}}$")
fig_closure.tight_layout()
fig_closure.savefig(
    OUTPUT_DIRECTORY / "q8_t950_first_order_closure.pdf")
plt.show()
"""
        ),
        code(
            r"""
fig_fit_quality, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.0))
for estimator_index, label in enumerate(VECTOR_LABELS):
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_RESULTS["fitted_block_residuals"][
            estimator_index],
        "o-",
        label=label,
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_RESULTS["condition_numbers"][
            estimator_index],
        "o-",
        label=label,
    )
for estimator_index, label in enumerate(
        CLOCK_STRAIN_LABELS):
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_STRAIN_RESULTS["fitted_block_residuals"][
            estimator_index],
        "s-",
        label=label,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_STRAIN_RESULTS["condition_numbers"][
            estimator_index],
        "s-",
        label=label,
    )
axes[0, 0].set_title("Vector fitted-block residual")
axes[0, 1].set_title("Clock/strain fitted-block residual")
axes[1, 0].set_title("Vector response condition number")
axes[1, 1].set_title("Clock/strain response condition number")
for axis in axes[0]:
    axis.set_ylabel("fit residual / fitted data")
for axis in axes[1]:
    axis.set_ylabel("condition number")
for axis in axes.flat:
    axis.set_xlabel(r"$R/M_B$")
    axis.legend(fontsize=7)
fig_fit_quality.tight_layout()
fig_fit_quality.savefig(
    OUTPUT_DIRECTORY / "q8_t950_first_order_fit_quality.pdf")
plt.show()
"""
        ),
        code(
            r"""
fig_parameters, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.0))

for estimator_index, label in enumerate(VECTOR_LABELS):
    estimates = VECTOR_RESULTS["estimates"][estimator_index]
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.linalg.norm(estimates[:, :3], axis=1),
        "o-",
        label=label,
    )
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.linalg.norm(estimates[:, 3:6], axis=1),
        "o-",
        label=label,
    )

for estimator_index, label in enumerate(
        CLOCK_STRAIN_LABELS):
    estimates = CLOCK_STRAIN_RESULTS["estimates"][
        estimator_index]
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.abs(estimates[:, 0]),
        "s-",
        label=label,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.linalg.norm(estimates[:, 1:], axis=1),
        "s-",
        label=label,
    )

axes[0, 0].set_title(r"$\|\beta_i\|$")
axes[0, 1].set_title(r"$\|\dot q^i\|$")
axes[1, 0].set_title(r"$|\dot q^0|$")
axes[1, 1].set_title(r"$\|\sigma_{ij}\|$")
for axis in axes.flat:
    axis.set_xlabel(r"$R/M_B$")
    axis.set_ylabel("fitted coefficient norm")
    axis.legend(fontsize=7)
fig_parameters.suptitle(
    "First-order estimates from every minimal block choice")
fig_parameters.tight_layout()
fig_parameters.savefig(
    OUTPUT_DIRECTORY / "q8_t950_first_order_parameters.pdf")
plt.show()
"""
        ),
        markdown(
            r"""
## 6. Out-of-sample inverse-metric reconstruction

Surface closure is an intermediate test.  The quantity ultimately needed by
the worldtube scheme is the metric extrapolated into the excised region.  We
therefore draw one fixed set of random validation points strictly inside the
smallest fitting sphere, $R=0.2M_{\rm tot}$.  The numerical domain is excised
at smaller radii, so the validation points occupy three volume-uniform shells,

$$
0.105 \leq \rho < 0.135,\qquad
0.135 \leq \rho < 0.165,\qquad
0.165 \leq \rho < 0.195,
$$

with 256 points per shell.  None of these 768 points participates in a
surface fit.

Every vector estimator is combined with every clock/strain estimator.  This
gives 25 complete first-order parameter vectors at each fitting radius.  For
each one we evaluate

$$
E_{\rm total}
=\frac{\lVert G_{\rm numerical}-G_{\rm model}\rVert_2}
{\lVert G_{\rm numerical}\rVert_2},
\qquad
E_{\rm pert}
=\frac{\lVert G_{\rm numerical}-G_{\rm model}\rVert_2}
{\lVert G_{\rm numerical}-G_{(0)}\rVert_2}.
$$

Here $G_{(0)}$ is harmonic Schwarzschild.  The first normalization measures
the error in the full inverse metric.  The second asks what fraction of the
non-Schwarzschild content remains.  The norm contains all ten independent
inverse-metric components at all validation points.
"""
        ),
        code(
            r"""
INTERIOR_RADIAL_EDGES = np.array([
    0.105, 0.135, 0.165, 0.195,
])
INTERIOR_POINTS_PER_SHELL = 256
INTERIOR_RANDOM_SEED = 20260728
interior_rng = np.random.default_rng(INTERIOR_RANDOM_SEED)

interior_local_points = []
interior_shell_indices = []
for shell_index, (inner_radius, outer_radius) in enumerate(zip(
        INTERIOR_RADIAL_EDGES[:-1],
        INTERIOR_RADIAL_EDGES[1:])):
    directions = interior_rng.normal(
        size=(INTERIOR_POINTS_PER_SHELL, 3))
    directions /= np.linalg.norm(
        directions, axis=1)[:, None]
    radii = (
        inner_radius**3
        + interior_rng.random(INTERIOR_POINTS_PER_SHELL)
        * (outer_radius**3 - inner_radius**3)
    ) ** (1.0 / 3.0)
    interior_local_points.append(
        radii[:, None] * directions)
    interior_shell_indices.append(np.full(
        INTERIOR_POINTS_PER_SHELL,
        shell_index,
        dtype=int,
    ))

INTERIOR_LOCAL_POINTS = np.concatenate(
    interior_local_points)
INTERIOR_SHELL_INDEX = np.concatenate(
    interior_shell_indices)
INTERIOR_RADII = np.linalg.norm(
    INTERIOR_LOCAL_POINTS, axis=1)
INTERIOR_INERTIAL_POINTS = (
    CURVATURE_CENTER[None, :] + INTERIOR_LOCAL_POINTS)

interior_raw_metric = interpolate_components(
    TEST_OBSERVATION_ID,
    INTERIOR_INERTIAL_POINTS,
    INVERSE_METRIC_COMPONENT_NAMES,
)
INTERIOR_NUMERICAL_METRIC = unpack_inverse_metric(
    interior_raw_metric)
INTERIOR_BACKGROUND_METRIC = inverse_harmonic_schwarzschild(
    INTERIOR_LOCAL_POINTS, SMALL_BLACK_HOLE_MASS)
INTERIOR_RESPONSE_COLUMNS = first_order_response_columns(
    INTERIOR_LOCAL_POINTS, SMALL_BLACK_HOLE_MASS)


def metric_values_by_point(metric):
    # Pack the ten independent components at every point.
    gtt, gti, gij = metric
    return np.column_stack([
        gtt,
        gti,
        *[gij[:, i, j] for i, j in SPATIAL_COMPONENTS],
    ])


INTERIOR_NUMERICAL_VALUES = metric_values_by_point(
    INTERIOR_NUMERICAL_METRIC)
INTERIOR_BACKGROUND_VALUES = metric_values_by_point(
    INTERIOR_BACKGROUND_METRIC)
INTERIOR_RESPONSE_DESIGN = np.stack([
    metric_values_by_point(column)
    for column in INTERIOR_RESPONSE_COLUMNS
], axis=-1)
INTERIOR_PERTURBATION_VALUES = (
    INTERIOR_NUMERICAL_VALUES
    - INTERIOR_BACKGROUND_VALUES)

number_of_vector_estimators = len(
    VECTOR_ESTIMATOR_BLOCKS)
number_of_clock_estimators = len(
    CLOCK_STRAIN_ESTIMATOR_BLOCKS)
number_of_radii = len(RADIUS_TEST_VALUES)
number_of_shells = len(INTERIOR_RADIAL_EDGES) - 1

INTERIOR_FIRST_ORDER_PARAMETERS = np.zeros((
    number_of_vector_estimators,
    number_of_clock_estimators,
    number_of_radii,
    len(FIRST_ORDER_PARAMETER_NAMES),
))
INTERIOR_TOTAL_ERROR = np.empty((
    number_of_vector_estimators,
    number_of_clock_estimators,
    number_of_radii,
))
INTERIOR_PERTURBATION_RESIDUAL = np.empty_like(
    INTERIOR_TOTAL_ERROR)
INTERIOR_SHELL_ERROR = np.empty((
    number_of_vector_estimators,
    number_of_clock_estimators,
    number_of_radii,
    number_of_shells,
))

numerical_norm = np.linalg.norm(
    INTERIOR_NUMERICAL_VALUES)
perturbation_norm = np.linalg.norm(
    INTERIOR_PERTURBATION_VALUES)
INTERIOR_ZEROTH_ORDER_ERROR = (
    perturbation_norm / numerical_norm)

for vector_index in range(number_of_vector_estimators):
    for clock_index in range(number_of_clock_estimators):
        for radius_index in range(number_of_radii):
            parameters = INTERIOR_FIRST_ORDER_PARAMETERS[
                vector_index, clock_index, radius_index]
            parameters[VECTOR_PARAMETER_INDICES] = (
                VECTOR_RESULTS["estimates"][
                    vector_index, radius_index])
            parameters[CLOCK_STRAIN_PARAMETER_INDICES] = (
                CLOCK_STRAIN_RESULTS["estimates"][
                    clock_index, radius_index])

            model_values = (
                INTERIOR_BACKGROUND_VALUES
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_RESPONSE_DESIGN,
                    parameters,
                )
            )
            residual = (
                INTERIOR_NUMERICAL_VALUES - model_values)
            INTERIOR_TOTAL_ERROR[
                vector_index, clock_index, radius_index
            ] = np.linalg.norm(residual) / numerical_norm
            INTERIOR_PERTURBATION_RESIDUAL[
                vector_index, clock_index, radius_index
            ] = np.linalg.norm(residual) / perturbation_norm

            for shell_index in range(number_of_shells):
                in_shell = (
                    INTERIOR_SHELL_INDEX == shell_index)
                INTERIOR_SHELL_ERROR[
                    vector_index,
                    clock_index,
                    radius_index,
                    shell_index,
                ] = (
                    np.linalg.norm(residual[in_shell])
                    / np.linalg.norm(
                        INTERIOR_NUMERICAL_VALUES[in_shell])
                )

median_perturbation_by_pair = np.median(
    INTERIOR_PERTURBATION_RESIDUAL, axis=2)
SELECTED_INTERIOR_PAIR = np.unravel_index(
    np.argmin(median_perturbation_by_pair),
    median_perturbation_by_pair.shape,
)
selected_vector_index, selected_clock_index = (
    SELECTED_INTERIOR_PAIR)
selected_pair_errors = INTERIOR_TOTAL_ERROR[
    selected_vector_index, selected_clock_index]
SELECTED_INTERIOR_RADIUS_INDEX = int(np.argmin(
    selected_pair_errors))
GLOBAL_INTERIOR_MINIMUM_INDEX = np.unravel_index(
    np.argmin(INTERIOR_TOTAL_ERROR),
    INTERIOR_TOTAL_ERROR.shape,
)
(
    global_vector_index,
    global_clock_index,
    global_radius_index,
) = GLOBAL_INTERIOR_MINIMUM_INDEX

print(
    f"held-out interior points: {len(INTERIOR_LOCAL_POINTS)}")
print(
    "held-out radial range:",
    f"[{INTERIOR_RADII.min():.6f}, "
    f"{INTERIOR_RADII.max():.6f}] M_tot",
)
print(
    "zeroth-order total relative RMS:",
    f"{INTERIOR_ZEROTH_ORDER_ERROR:.6e}",
)
print(
    "lowest-median first-order pair:",
    f"V{selected_vector_index + 1} + "
    f"C{selected_clock_index + 1}",
)
print(
    "best total relative RMS for that pair:",
    f"{selected_pair_errors[SELECTED_INTERIOR_RADIUS_INDEX]:.6e}",
    "at R =",
    f"{RADIUS_TEST_VALUES[SELECTED_INTERIOR_RADIUS_INDEX]:.3f} M_tot",
)
print(
    "smallest individual first-order RMS:",
    f"{INTERIOR_TOTAL_ERROR[GLOBAL_INTERIOR_MINIMUM_INDEX]:.6e}",
    "from",
    f"V{global_vector_index + 1} + C{global_clock_index + 1}",
    "at R =",
    f"{RADIUS_TEST_VALUES[global_radius_index]:.3f} M_tot",
)
"""
        ),
        code(
            r"""
fig_interior, axes = plt.subplots(
    2, 2, figsize=(14.0, 10.0))

all_total_errors = INTERIOR_TOTAL_ERROR.reshape(
    -1, number_of_radii)
all_perturbation_residuals = (
    INTERIOR_PERTURBATION_RESIDUAL.reshape(
        -1, number_of_radii))

for values in all_total_errors:
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        values,
        color="0.75",
        linewidth=0.8,
        alpha=0.7,
    )
for values in all_perturbation_residuals:
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        values,
        color="0.75",
        linewidth=0.8,
        alpha=0.7,
    )

selected_label = (
    f"V{selected_vector_index + 1} + "
    f"C{selected_clock_index + 1} "
    "(lowest median)")
axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    selected_pair_errors,
    "o-",
    linewidth=2.0,
    label=selected_label,
)
axes[0, 0].axhline(
    INTERIOR_ZEROTH_ORDER_ERROR,
    color="black",
    linestyle="--",
    label="zeroth order",
)
axes[0, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    INTERIOR_PERTURBATION_RESIDUAL[
        selected_vector_index, selected_clock_index],
    "o-",
    linewidth=2.0,
    label=selected_label,
)
axes[0, 1].axhline(
    1.0,
    color="black",
    linestyle="--",
    label="zeroth order",
)

axes[0, 0].set_title(
    "All ten inverse-metric components")
axes[0, 0].set_ylabel(
    r"$\|G_{\rm num}-G_{\rm model}\|/"
    r"\|G_{\rm num}\|$")
axes[0, 1].set_title(
    "Residual non-Schwarzschild content")
axes[0, 1].set_ylabel(
    r"$\|G_{\rm num}-G_{\rm model}\|/"
    r"\|G_{\rm num}-G_{(0)}\|$")
for axis in axes[0]:
    axis.set_xlabel(r"$R/M_B$")
    axis.legend(fontsize=8)

minimum_error_by_pair = np.min(
    INTERIOR_TOTAL_ERROR, axis=2)
best_radius_by_pair = RADIUS_TEST_VALUES[
    np.argmin(INTERIOR_TOTAL_ERROR, axis=2)]

image_error = axes[1, 0].imshow(
    minimum_error_by_pair,
    origin="lower",
    aspect="auto",
)
image_radius = axes[1, 1].imshow(
    best_radius_by_pair,
    origin="lower",
    aspect="auto",
)
for vector_index in range(number_of_vector_estimators):
    for clock_index in range(number_of_clock_estimators):
        axes[1, 0].text(
            clock_index,
            vector_index,
            f"{minimum_error_by_pair[vector_index, clock_index]:.3f}",
            ha="center",
            va="center",
            color="white",
            fontsize=8,
        )
        axes[1, 1].text(
            clock_index,
            vector_index,
            f"{best_radius_by_pair[vector_index, clock_index]:.2f}",
            ha="center",
            va="center",
            color="white",
            fontsize=8,
        )

for axis in axes[1]:
    axis.set_xticks(
        np.arange(number_of_clock_estimators),
        [f"C{index + 1}"
         for index in range(number_of_clock_estimators)],
    )
    axis.set_yticks(
        np.arange(number_of_vector_estimators),
        [f"V{index + 1}"
         for index in range(number_of_vector_estimators)],
    )
    axis.set_xlabel("clock/strain estimator")
    axis.set_ylabel("vector estimator")
axes[1, 0].set_title(
    "Minimum total RMS over fitting radius")
axes[1, 1].set_title(
    r"Fitting radius $R/M_{\rm tot}$ at minimum")
fig_interior.colorbar(
    image_error, ax=axes[1, 0], shrink=0.8)
fig_interior.colorbar(
    image_radius, ax=axes[1, 1], shrink=0.8)

fig_interior.suptitle(
    rf"Held-out interior inverse metric at "
    rf"$t={TEST_TIME:g}M_{{\rm tot}}$ "
    rf"({len(INTERIOR_LOCAL_POINTS)} fixed random points)")
fig_interior.tight_layout()
fig_interior.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_interior_residual.pdf")
plt.show()


fig_interior_shells, axis = plt.subplots(
    1, 1, figsize=(8.5, 5.5))
for shell_index, (inner_radius, outer_radius) in enumerate(zip(
        INTERIOR_RADIAL_EDGES[:-1],
        INTERIOR_RADIAL_EDGES[1:])):
    axis.loglog(
        RADIUS_OVER_SMALL_MASS,
        INTERIOR_SHELL_ERROR[
            selected_vector_index,
            selected_clock_index,
            :,
            shell_index,
        ],
        "o-",
        label=(
            fr"${inner_radius:.3f}\leq\rho"
            fr"<{outer_radius:.3f}$"
        ),
    )
axis.set_xlabel(r"$R/M_B$")
axis.set_ylabel("ten-component relative RMS")
axis.set_title(
    selected_label + ": error by interior shell")
axis.legend()
fig_interior_shells.tight_layout()
fig_interior_shells.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_interior_shells.pdf")
plt.show()
"""
        ),
        markdown(
            r"""
## 7. Radius-dependent first-order sensitivity

The response model is linear in all thirteen first-order parameters, so its
Jacobian is exactly the set of response columns already used in the fit.  No
finite differencing is required.

We distinguish four quantities that are sometimes all called
“sensitivity”:

1. **Raw surface response:** the RMS change in the worldtube metric per unit
   change in one parameter.
2. **Absolute identifiability:** the smallest singular value of a block or
   estimator response matrix.  A small value means that some parameter
   combination produces little measurable surface signal.
3. **Empirical parameter ambiguity:** the spread among the 25 algebraically
   viable estimator pairs at the same radius.
4. **Interior impact:** the metric change on the fixed random validation
   points produced by that ambiguity.

For a fitted surface system $d_R=A_Rp+h_R$,

$$
\delta p=A_R^+h_R,\qquad
\delta G_{\rm int}=J_{\rm int}A_R^+h_R.
$$

We therefore also compute the operator norm of
$J_{\rm int}A_R^+$.  It measures how strongly unit RMS contamination in the
selected surface blocks can be amplified into the held-out interior metric.
"""
        ),
        code(
            r"""
number_of_parameters = len(FIRST_ORDER_PARAMETER_NAMES)

# Raw response strength on each fitting sphere.
SURFACE_PARAMETER_RESPONSE_RMS = np.empty((
    number_of_parameters, number_of_radii))
for radius_index, radius in enumerate(RADIUS_TEST_VALUES):
    surface_columns = first_order_response_columns(
        radius * DIRECTIONS, SMALL_BLACK_HOLE_MASS)
    for parameter_index, column in enumerate(surface_columns):
        values = metric_values_by_point(column)
        SURFACE_PARAMETER_RESPONSE_RMS[
            parameter_index, radius_index
        ] = np.linalg.norm(values) / np.sqrt(values.size)


# Absolute sensitivity of each individual STF block, restricted to the
# parameter directions to which that block responds.
BLOCK_SMALLEST_NONZERO_SINGULAR_VALUE = np.empty((
    len(FIRST_ORDER_BLOCKS), number_of_radii))
BLOCK_ACTIVE_PARAMETER_COUNT = np.empty(
    len(FIRST_ORDER_BLOCKS), dtype=int)
for block_index, block in enumerate(FIRST_ORDER_BLOCKS):
    representative_matrix = BLOCK_RESPONSE[block][0]
    column_norms = np.linalg.norm(
        representative_matrix, axis=0)
    active_parameters = np.flatnonzero(
        column_norms
        > 1.0e-12 * max(column_norms.max(), 1.0))
    BLOCK_ACTIVE_PARAMETER_COUNT[block_index] = len(
        active_parameters)
    for radius_index in range(number_of_radii):
        matrix = BLOCK_RESPONSE[block][radius_index][
            :, active_parameters]
        singular_values = np.linalg.svd(
            matrix, compute_uv=False)
        nonzero = singular_values[
            singular_values
            > 1.0e-10 * singular_values[0]]
        BLOCK_SMALLEST_NONZERO_SINGULAR_VALUE[
            block_index, radius_index
        ] = (
            nonzero[-1] / np.sqrt(matrix.shape[0])
            if len(nonzero) else 0.0
        )


def estimator_stability(
        parameter_indices, estimator_block_sets):
    smallest_singular_values = np.empty((
        len(estimator_block_sets), number_of_radii))
    interior_amplification = np.empty_like(
        smallest_singular_values)
    interior_design = INTERIOR_RESPONSE_DESIGN[
        :, :, parameter_indices].reshape(
            -1, len(parameter_indices))

    for estimator_index, selected_blocks in enumerate(
            estimator_block_sets):
        for radius_index in range(number_of_radii):
            matrix, _ = stacked_system(
                radius_index,
                parameter_indices,
                selected_blocks,
            )
            singular_values = np.linalg.svd(
                matrix, compute_uv=False)
            smallest_singular_values[
                estimator_index, radius_index
            ] = (
                singular_values[-1]
                / np.sqrt(matrix.shape[0])
            )
            # If the surface contamination has unit RMS, its Euclidean norm
            # is sqrt(number of fitted data).  Divide the predicted interior
            # norm by the numerical metric norm.
            interior_amplification[
                estimator_index, radius_index
            ] = (
                np.linalg.norm(
                    interior_design @ np.linalg.pinv(matrix),
                    ord=2,
                )
                * np.sqrt(matrix.shape[0])
                / numerical_norm
            )
    return smallest_singular_values, interior_amplification


(
    VECTOR_SMALLEST_SINGULAR_VALUE,
    VECTOR_INTERIOR_AMPLIFICATION,
) = estimator_stability(
    VECTOR_PARAMETER_INDICES,
    VECTOR_ESTIMATOR_BLOCKS,
)
(
    CLOCK_SMALLEST_SINGULAR_VALUE,
    CLOCK_INTERIOR_AMPLIFICATION,
) = estimator_stability(
    CLOCK_STRAIN_PARAMETER_INDICES,
    CLOCK_STRAIN_ESTIMATOR_BLOCKS,
)


# Empirical ambiguity among every algebraically viable complete first-order
# estimator pair.
PARAMETER_ENSEMBLE = INTERIOR_FIRST_ORDER_PARAMETERS.reshape(
    -1, number_of_radii, number_of_parameters)
PARAMETER_ENSEMBLE_MEAN = np.mean(
    PARAMETER_ENSEMBLE, axis=0)
PARAMETER_ENSEMBLE_SPREAD = np.std(
    PARAMETER_ENSEMBLE, axis=0)

interior_column_norm = np.array([
    np.linalg.norm(INTERIOR_RESPONSE_DESIGN[:, :, index])
    for index in range(number_of_parameters)
])
INTERIOR_PARAMETER_CONTRIBUTION = (
    np.abs(PARAMETER_ENSEMBLE_MEAN).T
    * interior_column_norm[:, None]
    / numerical_norm
)
INTERIOR_PARAMETER_UNCERTAINTY_IMPACT = (
    PARAMETER_ENSEMBLE_SPREAD.T
    * interior_column_norm[:, None]
    / numerical_norm
)

FIRST_ORDER_GROUPS = (
    ("clock rate", np.array([0])),
    ("boost", np.arange(1, 4)),
    ("velocity", np.arange(4, 7)),
    ("strain", np.arange(7, 13)),
)
INTERIOR_GROUP_PREDICTION_SPREAD = np.empty((
    len(FIRST_ORDER_GROUPS), number_of_radii))
INTERIOR_TOTAL_PREDICTION_SPREAD = np.empty(
    number_of_radii)

for radius_index in range(number_of_radii):
    parameter_deviation = (
        PARAMETER_ENSEMBLE[:, radius_index]
        - PARAMETER_ENSEMBLE_MEAN[radius_index])
    total_metric_deviation = np.einsum(
        "ncp,ep->enc",
        INTERIOR_RESPONSE_DESIGN,
        parameter_deviation,
    )
    INTERIOR_TOTAL_PREDICTION_SPREAD[radius_index] = (
        np.sqrt(np.mean(np.sum(
            total_metric_deviation**2, axis=(1, 2))))
        / numerical_norm
    )

    for group_index, (_, parameter_indices) in enumerate(
            FIRST_ORDER_GROUPS):
        group_metric_deviation = np.einsum(
            "ncg,eg->enc",
            INTERIOR_RESPONSE_DESIGN[
                :, :, parameter_indices],
            parameter_deviation[:, parameter_indices],
        )
        INTERIOR_GROUP_PREDICTION_SPREAD[
            group_index, radius_index
        ] = (
            np.sqrt(np.mean(np.sum(
                group_metric_deviation**2,
                axis=(1, 2),
            )))
            / numerical_norm
        )

print(
    "small-radius / minimum raw response ratio by parameter:")
minimum_radius_index = 0
for parameter_index, name in enumerate(
        FIRST_ORDER_PARAMETER_NAMES):
    values = SURFACE_PARAMETER_RESPONSE_RMS[
        parameter_index]
    print(
        f"  {name:10s}: "
        f"{values[minimum_radius_index] / values.min():.6e}"
    )

print("\nLargest individual uncertainty impact at each radius:")
for radius_index, radius in enumerate(RADIUS_TEST_VALUES):
    parameter_index = int(np.argmax(
        INTERIOR_PARAMETER_UNCERTAINTY_IMPACT[:, radius_index]))
    largest_impact = INTERIOR_PARAMETER_UNCERTAINTY_IMPACT[
        parameter_index, radius_index]
    print(
        f"  R={radius:.3f}: "
        f"{FIRST_ORDER_PARAMETER_NAMES[parameter_index]} "
        f"({largest_impact:.6e})"
    )
"""
        ),
        code(
            r"""
def heatmap(
        axis, values, row_labels, title, colorbar_label,
        floor=1.0e-16):
    log_values = np.log10(np.maximum(values, floor))
    image = axis.imshow(
        log_values,
        origin="lower",
        aspect="auto",
    )
    axis.set_yticks(
        np.arange(len(row_labels)), row_labels)
    tick_indices = np.arange(
        0, number_of_radii, 2)
    axis.set_xticks(
        tick_indices,
        [f"{RADIUS_OVER_SMALL_MASS[index]:g}"
         for index in tick_indices],
        rotation=35,
        ha="right",
    )
    axis.set_xlabel(r"$R/M_B$")
    axis.set_title(title)
    colorbar = axis.figure.colorbar(
        image, ax=axis, shrink=0.85)
    colorbar.set_label(colorbar_label)
    return image


fig_response_sensitivity, axes = plt.subplots(
    2, 1, figsize=(13.5, 12.0))
heatmap(
    axes[0],
    SURFACE_PARAMETER_RESPONSE_RMS,
    FIRST_ORDER_PARAMETER_NAMES,
    "Raw surface response per unit parameter",
    r"$\log_{10}$ metric RMS",
)
heatmap(
    axes[1],
    BLOCK_SMALLEST_NONZERO_SINGULAR_VALUE,
    [block_name(block) for block in FIRST_ORDER_BLOCKS],
    "Smallest nonzero singular value of each STF block",
    r"$\log_{10}(s_{\min}/\sqrt{N_{\rm block}})$",
)
fig_response_sensitivity.suptitle(
    "Analytic first-order response strength versus radius")
fig_response_sensitivity.tight_layout()
fig_response_sensitivity.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_response_sensitivity.pdf")
plt.show()


fig_estimator_stability, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))
for estimator_index, label in enumerate(VECTOR_LABELS):
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_SMALLEST_SINGULAR_VALUE[estimator_index],
        "o-",
        label=label,
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        VECTOR_INTERIOR_AMPLIFICATION[estimator_index],
        "o-",
        label=label,
    )
for estimator_index, label in enumerate(
        CLOCK_STRAIN_LABELS):
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_SMALLEST_SINGULAR_VALUE[estimator_index],
        "s-",
        label=label,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        CLOCK_INTERIOR_AMPLIFICATION[estimator_index],
        "s-",
        label=label,
    )
axes[0, 0].set_title(
    "Vector absolute identifiability")
axes[0, 1].set_title(
    "Clock/strain absolute identifiability")
axes[1, 0].set_title(
    "Vector surface-to-interior amplification")
axes[1, 1].set_title(
    "Clock/strain surface-to-interior amplification")
for axis in axes[0]:
    axis.set_ylabel(
        r"$s_{\min}/\sqrt{N_{\rm fitted}}$")
for axis in axes[1]:
    axis.set_ylabel(
        "interior relative error / unit surface RMS")
for axis in axes.flat:
    axis.set_xlabel(r"$R/M_B$")
    axis.legend(fontsize=7)
fig_estimator_stability.tight_layout()
fig_estimator_stability.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_estimator_sensitivity.pdf")
plt.show()


fig_uncertainty, axes = plt.subplots(
    2, 2, figsize=(14.0, 10.5))
heatmap(
    axes[0, 0],
    PARAMETER_ENSEMBLE_SPREAD.T,
    FIRST_ORDER_PARAMETER_NAMES,
    "Empirical parameter spread",
    r"$\log_{10}\sigma(p_A)$",
)
heatmap(
    axes[0, 1],
    INTERIOR_PARAMETER_CONTRIBUTION,
    FIRST_ORDER_PARAMETER_NAMES,
    "Mean fitted contribution to interior metric",
    r"$\log_{10}(\|J_A\bar p_A\|/\|G_{\rm num}\|)$",
)
heatmap(
    axes[1, 0],
    INTERIOR_PARAMETER_UNCERTAINTY_IMPACT,
    FIRST_ORDER_PARAMETER_NAMES,
    "Interior impact of parameter ambiguity",
    r"$\log_{10}(\|J_A\|\sigma_A/\|G_{\rm num}\|)$",
)
for group_index, (group_name, _) in enumerate(
        FIRST_ORDER_GROUPS):
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        INTERIOR_GROUP_PREDICTION_SPREAD[group_index],
        "o-",
        label=group_name,
    )
axes[1, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    INTERIOR_TOTAL_PREDICTION_SPREAD,
    "k--",
    linewidth=2.0,
    label="all parameters",
)
axes[1, 1].set_xlabel(r"$R/M_B$")
axes[1, 1].set_ylabel(
    "RMS estimator-induced interior spread")
axes[1, 1].set_title(
    "Propagated model-choice ambiguity")
axes[1, 1].legend()
fig_uncertainty.suptitle(
    "First-order coefficient sensitivity and interior impact")
fig_uncertainty.tight_layout()
fig_uncertainty.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_uncertainty_propagation.pdf")
plt.show()
"""
        ),
        markdown(
            r"""
### Interpretation for choosing first-order blocks

The analytic and numerical tests agree that the small-radius branch is **not**
caused by a divergent first-order response:

- every raw parameter response changes by only an order-unity factor over the
  scan;
- every minimal estimator has better absolute identifiability and smaller
  surface-to-interior amplification at the smallest radii than at large
  radius;
- the strong large-radius amplification is exactly where the analytic table
  predicts that the separating blocks decay.

The inward growth instead comes from disagreement among otherwise
well-conditioned clock/strain fits.  At $R=0.2M_{\rm tot}$ the inferred clock
rates are

$$
\dot q^0=(0.0747,\ 0.1057,\ 0.0767,\ 0.2611,\ 0.0811)
$$

for C1 through C5.  C4 is the conspicuous outlier.  The pairwise structure
localizes the problem:

- C4 and C5 share $G^{Ti}_{1}$ but replace $G^{ij}_{0}$ by $G^{ij}_{2}$;
- C2 and C3 share $G^{TT}_{0}$ and make the same replacement.

Both comparisons point to systematic non-first-order content in
$G^{ij}_{0}$ at small radius, rather than poor linear conditioning.  This
block retains an $O(1)$ strain response inward, so higher-order monopole
content can be absorbed efficiently into the fitted trace strain and then
feed back into the inferred clock rate.

For the present slice, the cleanest conservative hierarchy is therefore:

1. use V1, $G^{TT}_{1}+G^{Ti}_{0}$, for the vector sector; it has the largest
   small-radius singular value and the smallest interior amplification;
2. use C3, $G^{TT}_{0}+G^{ij}_{2}$, for the clock/strain sector;
3. use $G^{Ti}_{2}+G^{ij}_{1}$ as independent vector closure and
   $G^{Ti}_{1}$ as independent clock/strain closure;
4. retain $G^{ij}_{0}$ as a diagnostic absolute residual, but do not use it
   to determine the small-radius clock/strain parameters until its
   higher-order monopole content is modeled;
5. avoid interpreting a large *fractional* closure error in a block whose
   first-order response is radially weak.  Plot its absolute residual
   alongside the fractional one.

This is a data-driven block choice for one slice, not yet a universal
prescription.  It should be checked on later q=8 slices and against the q=4
data before being promoted to the matching algorithm.
"""
        ),
        markdown(
            r"""
## 8. Provisional V1+C3 hierarchy

We now display the provisional hierarchy directly, without averaging it with
the other estimators:

$$
\begin{aligned}
\text{vector fit V1}:&\quad
G^{TT}_{1}+G^{Ti}_{0},
\\
\text{clock/strain fit C3}:&\quad
G^{TT}_{0}+G^{ij}_{2}.
\end{aligned}
$$

The primary independent closure blocks are

$$
G^{Ti}_{2}+G^{ij}_{1}
\quad\text{and}\quad
G^{Ti}_{1}.
$$

The spatial monopole $G^{ij}_{0}$ is shown only through its absolute RMS
residual.  It is not included in the fitted system or in the recommended
fractional closure score.
"""
        ),
        code(
            r"""
PROVISIONAL_VECTOR_ESTIMATOR_INDEX = 0  # V1
PROVISIONAL_CLOCK_ESTIMATOR_INDEX = 2   # C3

PROVISIONAL_PARAMETERS = np.zeros((
    number_of_radii, number_of_parameters))
PROVISIONAL_PARAMETERS[:, VECTOR_PARAMETER_INDICES] = (
    VECTOR_RESULTS["estimates"][
        PROVISIONAL_VECTOR_ESTIMATOR_INDEX])
PROVISIONAL_PARAMETERS[
    :, CLOCK_STRAIN_PARAMETER_INDICES] = (
        CLOCK_STRAIN_RESULTS["estimates"][
            PROVISIONAL_CLOCK_ESTIMATOR_INDEX])


def provisional_closure(selected_blocks):
    fractional = np.empty(number_of_radii)
    absolute_rms = np.empty(number_of_radii)
    for radius_index in range(number_of_radii):
        data = np.concatenate([
            BLOCK_DATA[block][radius_index]
            for block in selected_blocks
        ])
        matrix = np.concatenate([
            BLOCK_RESPONSE[block][radius_index]
            for block in selected_blocks
        ], axis=0)
        residual = (
            data
            - matrix @ PROVISIONAL_PARAMETERS[radius_index])
        fractional[radius_index] = (
            np.linalg.norm(residual)
            / max(np.linalg.norm(data), 1.0e-300)
        )
        absolute_rms[radius_index] = (
            np.linalg.norm(residual)
            / np.sqrt(len(residual))
        )
    return fractional, absolute_rms


PROVISIONAL_VECTOR_CLOSURE_BLOCKS = (
    ("GTI", 2), ("GIJ", 1))
(
    PROVISIONAL_VECTOR_CLOSURE_FRACTION,
    PROVISIONAL_VECTOR_CLOSURE_RMS,
) = provisional_closure(
    PROVISIONAL_VECTOR_CLOSURE_BLOCKS)
(
    PROVISIONAL_GTI2_CLOSURE_FRACTION,
    PROVISIONAL_GTI2_CLOSURE_RMS,
) = provisional_closure((("GTI", 2),))
(
    PROVISIONAL_GIJ1_CLOSURE_FRACTION,
    PROVISIONAL_GIJ1_CLOSURE_RMS,
) = provisional_closure((("GIJ", 1),))
(
    PROVISIONAL_CLOCK_CLOSURE_FRACTION,
    PROVISIONAL_CLOCK_CLOSURE_RMS,
) = provisional_closure((("GTI", 1),))
(
    PROVISIONAL_GIJ0_DIAGNOSTIC_FRACTION,
    PROVISIONAL_GIJ0_DIAGNOSTIC_RMS,
) = provisional_closure((("GIJ", 0),))


fig_provisional_values, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))

axes[0, 0].semilogx(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_PARAMETERS[:, 0],
    "o-",
    label=r"$\dot q^0$",
)
for component, label in enumerate(("x", "y", "z")):
    axes[0, 1].semilogx(
        RADIUS_OVER_SMALL_MASS,
        PROVISIONAL_PARAMETERS[:, 1 + component],
        "o-",
        label=fr"$\beta_{label}$",
    )
    axes[1, 0].semilogx(
        RADIUS_OVER_SMALL_MASS,
        PROVISIONAL_PARAMETERS[:, 4 + component],
        "o-",
        label=fr"$\dot q^{label}$",
    )

strain_labels = (
    r"$\sigma_{xx}$", r"$\sigma_{xy}$",
    r"$\sigma_{xz}$", r"$\sigma_{yy}$",
    r"$\sigma_{yz}$", r"$\sigma_{zz}$",
)
for component, label in enumerate(strain_labels):
    axes[1, 1].semilogx(
        RADIUS_OVER_SMALL_MASS,
        PROVISIONAL_PARAMETERS[:, 7 + component],
        "o-",
        label=label,
    )

axes[0, 0].set_title("Clock rate from C3")
axes[0, 1].set_title("Boost from V1")
axes[1, 0].set_title("Spatial velocity from V1")
axes[1, 1].set_title("Spatial strain from C3")
for axis in axes.flat:
    axis.axhline(0.0, color="0.5", linewidth=0.8)
    axis.set_xlabel(r"$R/M_B$")
    axis.set_ylabel("fitted coefficient")
    axis.legend(fontsize=8)
fig_provisional_values.suptitle(
    "Signed first-order coefficients in the provisional V1+C3 hierarchy")
fig_provisional_values.tight_layout()
fig_provisional_values.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_provisional_values.pdf")
plt.show()


fig_provisional_closure, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))
vector_closure_curves = (
    (
        "combined",
        PROVISIONAL_VECTOR_CLOSURE_FRACTION,
        PROVISIONAL_VECTOR_CLOSURE_RMS,
    ),
    (
        r"$G^{Ti}_{2}$",
        PROVISIONAL_GTI2_CLOSURE_FRACTION,
        PROVISIONAL_GTI2_CLOSURE_RMS,
    ),
    (
        r"$G^{ij}_{1}$",
        PROVISIONAL_GIJ1_CLOSURE_FRACTION,
        PROVISIONAL_GIJ1_CLOSURE_RMS,
    ),
)
for label, fractional, absolute_rms in vector_closure_curves:
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        fractional,
        "o-",
        label=label,
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        absolute_rms,
        "o-",
        label=label,
    )

axes[0, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_CLOCK_CLOSURE_FRACTION,
    "s-",
    label=r"$G^{Ti}_{1}$ primary closure",
)
axes[1, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_CLOCK_CLOSURE_RMS,
    "s-",
    label=r"$G^{Ti}_{1}$ primary closure",
)
axes[1, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_GIJ0_DIAGNOSTIC_RMS,
    "D--",
    label=r"$G^{ij}_{0}$ diagnostic only",
)

axes[0, 0].set_title("Vector fractional closure")
axes[1, 0].set_title("Vector absolute closure")
axes[0, 1].set_title("Clock/strain fractional closure")
axes[1, 1].set_title(
    "Clock/strain absolute residuals")
for axis in axes[0]:
    axis.set_ylabel(
        "closure residual / numerical block")
for axis in axes[1]:
    axis.set_ylabel("closure residual RMS")
for axis in axes.flat:
    axis.set_xlabel(r"$R/M_B$")
    axis.legend(fontsize=8)
fig_provisional_closure.suptitle(
    "Independent closures for the provisional V1+C3 hierarchy")
fig_provisional_closure.tight_layout()
fig_provisional_closure.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_provisional_closure.pdf")
plt.show()


def print_provisional_minimum(name, values):
    minimum_index = int(np.argmin(values))
    print(
        f"{name:<32s}: {values[minimum_index]:.6e} "
        f"at R={RADIUS_TEST_VALUES[minimum_index]:.3f} M_tot"
    )


print("\nProvisional hierarchy closure minima")
print_provisional_minimum(
    "vector combined fractional",
    PROVISIONAL_VECTOR_CLOSURE_FRACTION,
)
print_provisional_minimum(
    "clock GTI l=1 fractional",
    PROVISIONAL_CLOCK_CLOSURE_FRACTION,
)
print_provisional_minimum(
    "GIJ l=0 diagnostic RMS",
    PROVISIONAL_GIJ0_DIAGNOSTIC_RMS,
)
"""
        ),
        markdown(
            r"""
### Reading the provisional curves

The fitted coefficients are smooth over the small and intermediate radii.
The sharp drift in the boost and velocity components begins only at the
larger radii, in the same region where their separating response blocks lose
absolute sensitivity.

The closure plots also show why absolute and fractional residuals must be
kept together:

- the combined vector fractional closure is smallest at the innermost sphere;
- the fractional $G^{Ti}_{2}$ error rises rapidly with radius even while its
  absolute residual decreases, because the numerical $\ell=2$ block itself
  is becoming small;
- the $G^{Ti}_{1}$ clock/strain fractional closure has a convex minimum, but
  its absolute residual decreases throughout the scan.  Its large-radius
  fractional upturn is therefore predominantly loss of signal in the
  denominator, not growth of the mismatch;
- the absolute $G^{ij}_{0}$ diagnostic has a genuine convex shape and remains
  large at the smallest radii, consistent with the monopole-contamination
  diagnosis above.

Thus the most meaningful provisional closure is not one scalar fractional
score.  It is the collection of blockwise absolute residuals together with a
fractional score only where the corresponding first-order block retains
appreciable signal.
"""
        ),
        markdown(
            r"""
### Held-out interior metric for V1+C3

Finally, evaluate only the provisional hierarchy on the same 768 random
points inside the smallest worldtube.  These points were not used to choose
or fit any surface block.  We show the full ten-component error, the residual
non-Schwarzschild content, the three validation shells, and the pointwise
error distribution at the best V1+C3 fitting radius.
"""
        ),
        code(
            r"""
PROVISIONAL_INTERIOR_TOTAL_ERROR = INTERIOR_TOTAL_ERROR[
    PROVISIONAL_VECTOR_ESTIMATOR_INDEX,
    PROVISIONAL_CLOCK_ESTIMATOR_INDEX,
]
PROVISIONAL_INTERIOR_PERTURBATION_RESIDUAL = (
    INTERIOR_PERTURBATION_RESIDUAL[
        PROVISIONAL_VECTOR_ESTIMATOR_INDEX,
        PROVISIONAL_CLOCK_ESTIMATOR_INDEX,
    ]
)
PROVISIONAL_INTERIOR_SHELL_ERROR = INTERIOR_SHELL_ERROR[
    PROVISIONAL_VECTOR_ESTIMATOR_INDEX,
    PROVISIONAL_CLOCK_ESTIMATOR_INDEX,
]
PROVISIONAL_BEST_INTERIOR_RADIUS_INDEX = int(np.argmin(
    PROVISIONAL_INTERIOR_TOTAL_ERROR))
PROVISIONAL_BEST_INTERIOR_RADIUS = RADIUS_TEST_VALUES[
    PROVISIONAL_BEST_INTERIOR_RADIUS_INDEX]
PROVISIONAL_BEST_INTERIOR_ERROR = (
    PROVISIONAL_INTERIOR_TOTAL_ERROR[
        PROVISIONAL_BEST_INTERIOR_RADIUS_INDEX])
PROVISIONAL_INTERIOR_IMPROVEMENT = (
    INTERIOR_ZEROTH_ORDER_ERROR
    / PROVISIONAL_BEST_INTERIOR_ERROR)

provisional_best_model_values = (
    INTERIOR_BACKGROUND_VALUES
    + np.einsum(
        "ncp,p->nc",
        INTERIOR_RESPONSE_DESIGN,
        PROVISIONAL_PARAMETERS[
            PROVISIONAL_BEST_INTERIOR_RADIUS_INDEX],
    )
)
provisional_best_pointwise_residual = (
    INTERIOR_NUMERICAL_VALUES
    - provisional_best_model_values)
zeroth_pointwise_residual = (
    INTERIOR_NUMERICAL_VALUES
    - INTERIOR_BACKGROUND_VALUES)
pointwise_numerical_norm = np.linalg.norm(
    INTERIOR_NUMERICAL_VALUES, axis=1)
PROVISIONAL_POINTWISE_RELATIVE_ERROR = (
    np.linalg.norm(
        provisional_best_pointwise_residual, axis=1)
    / pointwise_numerical_norm
)
ZEROTH_POINTWISE_RELATIVE_ERROR = (
    np.linalg.norm(
        zeroth_pointwise_residual, axis=1)
    / pointwise_numerical_norm
)

pointwise_bin_edges = np.linspace(
    INTERIOR_RADII.min(),
    INTERIOR_RADII.max(),
    13,
)
pointwise_bin_centers = 0.5 * (
    pointwise_bin_edges[:-1] + pointwise_bin_edges[1:])
provisional_pointwise_median = np.full(
    len(pointwise_bin_centers), np.nan)
zeroth_pointwise_median = np.full(
    len(pointwise_bin_centers), np.nan)
for bin_index, (lower, upper) in enumerate(zip(
        pointwise_bin_edges[:-1],
        pointwise_bin_edges[1:])):
    in_bin = (
        (INTERIOR_RADII >= lower)
        & (INTERIOR_RADII < upper)
    )
    if np.any(in_bin):
        provisional_pointwise_median[bin_index] = np.median(
            PROVISIONAL_POINTWISE_RELATIVE_ERROR[in_bin])
        zeroth_pointwise_median[bin_index] = np.median(
            ZEROTH_POINTWISE_RELATIVE_ERROR[in_bin])


fig_provisional_interior, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))

axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_INTERIOR_TOTAL_ERROR,
    "o-",
    label="V1+C3",
)
axes[0, 0].axhline(
    INTERIOR_ZEROTH_ORDER_ERROR,
    color="black",
    linestyle="--",
    label="zeroth order",
)
axes[0, 0].set_title(
    "All ten inverse-metric components")
axes[0, 0].set_ylabel(
    r"$\|G_{\rm num}-G_{\rm model}\|/"
    r"\|G_{\rm num}\|$")

axes[0, 1].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_INTERIOR_PERTURBATION_RESIDUAL,
    "o-",
    label="V1+C3",
)
axes[0, 1].axhline(
    1.0,
    color="black",
    linestyle="--",
    label="zeroth order",
)
axes[0, 1].set_title(
    "Residual non-Schwarzschild content")
axes[0, 1].set_ylabel(
    r"$\|G_{\rm num}-G_{\rm model}\|/"
    r"\|G_{\rm num}-G_{(0)}\|$")

for shell_index, (inner_radius, outer_radius) in enumerate(zip(
        INTERIOR_RADIAL_EDGES[:-1],
        INTERIOR_RADIAL_EDGES[1:])):
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        PROVISIONAL_INTERIOR_SHELL_ERROR[:, shell_index],
        "o-",
        label=(
            fr"${inner_radius:.3f}\leq\rho"
            fr"<{outer_radius:.3f}$"
        ),
    )
axes[1, 0].set_title("Error by validation shell")
axes[1, 0].set_ylabel(
    "ten-component relative RMS")

axes[1, 1].scatter(
    INTERIOR_RADII / SMALL_BLACK_HOLE_MASS,
    ZEROTH_POINTWISE_RELATIVE_ERROR,
    s=7,
    alpha=0.12,
    color="0.4",
)
axes[1, 1].scatter(
    INTERIOR_RADII / SMALL_BLACK_HOLE_MASS,
    PROVISIONAL_POINTWISE_RELATIVE_ERROR,
    s=7,
    alpha=0.12,
    color="C0",
)
axes[1, 1].semilogy(
    pointwise_bin_centers / SMALL_BLACK_HOLE_MASS,
    zeroth_pointwise_median,
    "o--",
    color="0.2",
    linewidth=2.0,
    label="zeroth-order median",
)
axes[1, 1].semilogy(
    pointwise_bin_centers / SMALL_BLACK_HOLE_MASS,
    provisional_pointwise_median,
    "o-",
    color="C0",
    linewidth=2.0,
    label="V1+C3 median",
)
axes[1, 1].set_title(
    "Pointwise error at best fitting radius")
axes[1, 1].set_xlabel(r"validation radius $\rho/M_B$")
axes[1, 1].set_ylabel(
    "pointwise ten-component relative error")

for axis in axes[0]:
    axis.set_xlabel(r"fitting radius $R/M_B$")
for axis in axes[:, 0]:
    axis.legend(fontsize=8)
axes[0, 1].legend(fontsize=8)
axes[1, 1].legend(fontsize=8)

fig_provisional_interior.suptitle(
    rf"V1+C3 held-out interior inverse metric at "
    rf"$t={TEST_TIME:g}M_{{\rm tot}}$; "
    rf"best fit $R={PROVISIONAL_BEST_INTERIOR_RADIUS:g}"
    rf"M_{{\rm tot}}$")
fig_provisional_interior.tight_layout()
fig_provisional_interior.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_provisional_interior.pdf")
plt.show()

print(
    "\nProvisional V1+C3 held-out interior minimum:",
    f"{PROVISIONAL_BEST_INTERIOR_ERROR:.6e}",
    "at R =",
    f"{PROVISIONAL_BEST_INTERIOR_RADIUS:.3f} M_tot",
)
print(
    "improvement over zeroth order:",
    f"{PROVISIONAL_INTERIOR_IMPROVEMENT:.3f}x",
)
"""
        ),
        markdown(
            r"""
### Rendered V1+C3 interior-validation result

The preview below is stored beside this notebook, so it remains visible
without rerunning the preceding cells.

![V1+C3 held-out interior inverse-metric validation](q8_t950_first_order_provisional_interior.png)
"""
        ),
        markdown(
            r"""
## 9. Spatial radial derivative as a strictly held-out closure

The provisional V1+C3 parameters above were determined from values of the
inverse metric only.  We now test those same parameters against an independent
observable,

$$
\mathcal D G^{ab}
\equiv R n^i\partial_i G^{ab}.
$$

No derivative data enter the fit.  The evolved first-order reduction variable
is

$$
\Phi_{iab}=\partial_i g_{ab},
$$

so the derivative of the inverse metric follows algebraically:

$$
\partial_iG^{ab}
=-G^{ac}G^{bd}\Phi_{icd}.
$$

The factor of \(R\) makes \(\mathcal D\) dimensionless and prevents a trivial
overall \(1/R\) from dominating comparisons across the radius scan.  We apply
the same STF projections as for the metric values.  The analytic Schwarzschild
background and every first-order response column are differentiated along a
fixed angular ray with a centered five-point stencil.  This numerical
differentiation acts only on smooth analytic functions; it does not
finite-difference the simulation data.

Before using the result, we check the component ordering and inverse-metric
identity by comparing the \(\Phi\)-based derivative with a centered radial
finite difference of independently interpolated \(G^{ab}\).  That check is a
diagnostic only and does not enter either the fit or the closure.
"""
        ),
        code(
            r"""
COVARIANT_COMPONENT_SUFFIXES = (
    "tt", "xt", "xx", "yt", "yx",
    "yy", "zt", "zx", "zy", "zz",
)
COVARIANT_DATA_COMPONENTS = (
    (0, 0), (0, 1), (1, 1), (0, 2), (1, 2),
    (2, 2), (0, 3), (1, 3), (2, 3), (3, 3),
)
PHI_COMPONENT_NAMES = tuple(
    f"Phi_{derivative}{suffix}"
    for derivative in ("x", "y", "z")
    for suffix in COVARIANT_COMPONENT_SUFFIXES
)


def inverse_metric_matrix(field_data):
    gtt, gti, gij = unpack_inverse_metric(field_data)
    matrix = np.zeros((len(gtt), 4, 4))
    matrix[:, 0, 0] = gtt
    matrix[:, 0, 1:] = gti
    matrix[:, 1:, 0] = gti
    matrix[:, 1:, 1:] = gij
    return matrix


def metric_tuple_from_matrix(matrix):
    return (
        matrix[:, 0, 0],
        matrix[:, 0, 1:],
        matrix[:, 1:, 1:],
    )


def interpolate_phi(observation_id, center, radii):
    points = np.concatenate([
        center[None, :] + radius * DIRECTIONS
        for radius in radii
    ])
    values = interpolate_components(
        observation_id, points, PHI_COMPONENT_NAMES)
    values = values.reshape(
        3,
        len(COVARIANT_DATA_COMPONENTS),
        len(radii),
        len(DIRECTIONS),
    )
    phi = np.zeros((
        len(radii), len(DIRECTIONS), 3, 4, 4))
    for derivative in range(3):
        for component_values, (a, b) in zip(
                values[derivative],
                COVARIANT_DATA_COMPONENTS):
            phi[:, :, derivative, a, b] = component_values
            phi[:, :, derivative, b, a] = component_values
    return phi


PHI_COVARIANT = interpolate_phi(
    TEST_OBSERVATION_ID,
    CURVATURE_CENTER,
    RADIUS_TEST_VALUES,
)
INVERSE_METRIC_MATRICES = np.asarray([
    inverse_metric_matrix(
        FIRST_ORDER_METRIC_DATA[:, radius_index])
    for radius_index in range(len(RADIUS_TEST_VALUES))
])

# d_i G^{ab} = -G^{ac} G^{bd} Phi_{icd}
INVERSE_METRIC_SPATIAL_DERIVATIVE = -np.einsum(
    "rpac,rpbd,rpicd->rpiab",
    INVERSE_METRIC_MATRICES,
    INVERSE_METRIC_MATRICES,
    PHI_COVARIANT,
)
NUMERICAL_RADIAL_DERIVATIVE = np.einsum(
    "pi,rpiab->rpab",
    DIRECTIONS,
    INVERSE_METRIC_SPATIAL_DERIVATIVE,
)


# Independent implementation check: finite-difference only the interpolated
# inverse metric, not the closure model.
PHI_CHECK_STEP = 2.0e-4
metric_plus = interpolate_inverse_metric(
    TEST_OBSERVATION_ID,
    CURVATURE_CENTER,
    RADIUS_TEST_VALUES + PHI_CHECK_STEP,
)
metric_minus = interpolate_inverse_metric(
    TEST_OBSERVATION_ID,
    CURVATURE_CENTER,
    RADIUS_TEST_VALUES - PHI_CHECK_STEP,
)
finite_difference_radial_derivative = np.asarray([
    (
        inverse_metric_matrix(metric_plus[:, radius_index])
        - inverse_metric_matrix(metric_minus[:, radius_index])
    ) / (2.0 * PHI_CHECK_STEP)
    for radius_index in range(len(RADIUS_TEST_VALUES))
])
PHI_FINITE_DIFFERENCE_RELATIVE_ERROR = np.asarray([
    np.linalg.norm(
        NUMERICAL_RADIAL_DERIVATIVE[radius_index]
        - finite_difference_radial_derivative[radius_index])
    / np.linalg.norm(
        finite_difference_radial_derivative[radius_index])
    for radius_index in range(len(RADIUS_TEST_VALUES))
])

print(
    "Phi-to-inverse radial derivative check:",
    f"min={PHI_FINITE_DIFFERENCE_RELATIVE_ERROR.min():.3e},",
    f"max={PHI_FINITE_DIFFERENCE_RELATIVE_ERROR.max():.3e}",
)
"""
        ),
        code(
            r"""
def five_point_radial_derivative(
        metric_function, radius, relative_step=1.0e-4):
    step = relative_step * radius
    samples = [
        metric_function((radius + offset * step) * DIRECTIONS)
        for offset in (-2.0, -1.0, 1.0, 2.0)
    ]
    return tuple(
        (
            samples[0][component]
            - 8.0 * samples[1][component]
            + 8.0 * samples[2][component]
            - samples[3][component]
        ) / (12.0 * step)
        for component in range(3)
    )


def analytic_background_radial_derivative(
        radius, relative_step=1.0e-4):
    return five_point_radial_derivative(
        lambda points: inverse_harmonic_schwarzschild(
            points, SMALL_BLACK_HOLE_MASS),
        radius,
        relative_step,
    )


def analytic_response_radial_derivatives(
        radius, relative_step=1.0e-4):
    step = relative_step * radius
    samples = [
        first_order_response_columns(
            (radius + offset * step) * DIRECTIONS,
            SMALL_BLACK_HOLE_MASS,
        )
        for offset in (-2.0, -1.0, 1.0, 2.0)
    ]
    return tuple(
        tuple(
            (
                samples[0][parameter][component]
                - 8.0 * samples[1][parameter][component]
                + 8.0 * samples[2][parameter][component]
                - samples[3][parameter][component]
            ) / (12.0 * step)
            for component in range(3)
        )
        for parameter in range(len(FIRST_ORDER_PARAMETER_NAMES))
    )


DERIVATIVE_BLOCK_DATA = {
    block: [] for block in FIRST_ORDER_BLOCKS}
DERIVATIVE_BLOCK_RESPONSE = {
    block: [] for block in FIRST_ORDER_BLOCKS}
ANALYTIC_DERIVATIVE_STENCIL_ERROR = []

for radius_index, radius in enumerate(RADIUS_TEST_VALUES):
    numerical_derivative = metric_tuple_from_matrix(
        NUMERICAL_RADIAL_DERIVATIVE[radius_index])
    background_derivative = (
        analytic_background_radial_derivative(radius))
    response_derivatives = (
        analytic_response_radial_derivatives(radius))

    # Verify that halving the analytic stencil changes neither background nor
    # response derivatives appreciably.
    background_half_step = (
        analytic_background_radial_derivative(
            radius, relative_step=5.0e-5))
    response_half_step = (
        analytic_response_radial_derivatives(
            radius, relative_step=5.0e-5))
    full_step_values = np.concatenate([
        metric_values_by_point(background_derivative).ravel(),
        *[
            metric_values_by_point(column).ravel()
            for column in response_derivatives
        ],
    ])
    half_step_values = np.concatenate([
        metric_values_by_point(background_half_step).ravel(),
        *[
            metric_values_by_point(column).ravel()
            for column in response_half_step
        ],
    ])
    ANALYTIC_DERIVATIVE_STENCIL_ERROR.append(
        np.linalg.norm(full_step_values - half_step_values)
        / np.linalg.norm(half_step_values)
    )

    derivative_residual = tuple(
        radius * (numerical - background)
        for numerical, background in zip(
            numerical_derivative, background_derivative)
    )
    dimensionless_response_derivatives = tuple(
        tuple(radius * component for component in column)
        for column in response_derivatives
    )

    for block in FIRST_ORDER_BLOCKS:
        DERIVATIVE_BLOCK_DATA[block].append(
            metric_stf_block(*derivative_residual, block))
        DERIVATIVE_BLOCK_RESPONSE[block].append(
            np.column_stack([
                metric_stf_block(*response, block)
                for response in dimensionless_response_derivatives
            ])
        )

for block in FIRST_ORDER_BLOCKS:
    DERIVATIVE_BLOCK_DATA[block] = np.asarray(
        DERIVATIVE_BLOCK_DATA[block])
    DERIVATIVE_BLOCK_RESPONSE[block] = np.asarray(
        DERIVATIVE_BLOCK_RESPONSE[block])
ANALYTIC_DERIVATIVE_STENCIL_ERROR = np.asarray(
    ANALYTIC_DERIVATIVE_STENCIL_ERROR)

print(
    "analytic five-point stencil self-check:",
    f"max={ANALYTIC_DERIVATIVE_STENCIL_ERROR.max():.3e}",
)
"""
        ),
        code(
            r"""
def closure_from_blocks(data_blocks, response_blocks, selected_blocks):
    fractional = np.empty(number_of_radii)
    absolute_rms = np.empty(number_of_radii)
    for radius_index in range(number_of_radii):
        data = np.concatenate([
            data_blocks[block][radius_index]
            for block in selected_blocks
        ])
        matrix = np.concatenate([
            response_blocks[block][radius_index]
            for block in selected_blocks
        ], axis=0)
        residual = (
            data
            - matrix @ PROVISIONAL_PARAMETERS[radius_index])
        fractional[radius_index] = (
            np.linalg.norm(residual)
            / max(np.linalg.norm(data), 1.0e-300)
        )
        absolute_rms[radius_index] = (
            np.linalg.norm(residual)
            / np.sqrt(len(residual))
        )
    return fractional, absolute_rms


DERIVATIVE_CLOSURES = {}
derivative_block_groups = {
    "V1 fitted blocks": (("GTT", 1), ("GTI", 0)),
    "C3 fitted blocks": (("GTT", 0), ("GIJ", 2)),
    "vector closure": (("GTI", 2), ("GIJ", 1)),
    r"$G^{Ti}_{2}$": (("GTI", 2),),
    r"$G^{ij}_{1}$": (("GIJ", 1),),
    r"$G^{Ti}_{1}$": (("GTI", 1),),
    r"$G^{ij}_{0}$ diagnostic": (("GIJ", 0),),
}
for label, blocks in derivative_block_groups.items():
    DERIVATIVE_CLOSURES[label] = closure_from_blocks(
        DERIVATIVE_BLOCK_DATA,
        DERIVATIVE_BLOCK_RESPONSE,
        blocks,
    )


fig_derivative, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))

axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_VECTOR_CLOSURE_FRACTION,
    "o--",
    color="C0",
    label="metric values: vector",
)
axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    DERIVATIVE_CLOSURES["vector closure"][0],
    "o-",
    color="C0",
    label=r"$R\partial_R G$: vector",
)
axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_CLOCK_CLOSURE_FRACTION,
    "s--",
    color="C1",
    label=r"metric values: $G^{Ti}_{1}$",
)
axes[0, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    DERIVATIVE_CLOSURES[r"$G^{Ti}_{1}$"][0],
    "s-",
    color="C1",
    label=r"$R\partial_R G$: $G^{Ti}_{1}$",
)
axes[0, 0].set_title(
    "Primary held-out fractional closures")
axes[0, 0].set_ylabel(
    r"$\|d-\mathcal{R}\hat p\|/\|d\|$")

for label in (
        "V1 fitted blocks", "C3 fitted blocks",
        "vector closure", r"$G^{Ti}_{1}$"):
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        DERIVATIVE_CLOSURES[label][0],
        "o-",
        label=label,
    )
axes[0, 1].set_title(
    r"All independent $R\partial_R G$ tests")
axes[0, 1].set_ylabel("fractional derivative residual")

axes[1, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_VECTOR_CLOSURE_RMS,
    "o--",
    color="C0",
    label="metric values: vector",
)
axes[1, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    DERIVATIVE_CLOSURES["vector closure"][1],
    "o-",
    color="C0",
    label=r"$R\partial_R G$: vector",
)
axes[1, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    PROVISIONAL_CLOCK_CLOSURE_RMS,
    "s--",
    color="C1",
    label=r"metric values: $G^{Ti}_{1}$",
)
axes[1, 0].loglog(
    RADIUS_OVER_SMALL_MASS,
    DERIVATIVE_CLOSURES[r"$G^{Ti}_{1}$"][1],
    "s-",
    color="C1",
    label=r"$R\partial_R G$: $G^{Ti}_{1}$",
)
axes[1, 0].set_title(
    "Primary held-out absolute closures")
axes[1, 0].set_ylabel("absolute RMS")

for label in (
        r"$G^{Ti}_{2}$", r"$G^{ij}_{1}$",
        r"$G^{Ti}_{1}$", r"$G^{ij}_{0}$ diagnostic"):
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        DERIVATIVE_CLOSURES[label][1],
        "o-",
        label=label,
    )
axes[1, 1].set_title(
    r"Blockwise $R\partial_R G$ absolute residual")
axes[1, 1].set_ylabel("absolute RMS")

for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=8)
fig_derivative.suptitle(
    r"Metric-only V1+C3 fit tested against held-out radial derivatives")
fig_derivative.tight_layout()
fig_derivative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_radial_derivative_closure.pdf")
fig_derivative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_radial_derivative_closure.png",
    dpi=180,
)
plt.show()

print("\nHeld-out derivative minima")
for label, (fractional, absolute_rms) in (
        DERIVATIVE_CLOSURES.items()):
    best_index = int(np.argmin(fractional))
    print(
        f"{label:<28s} "
        f"fraction={fractional[best_index]:.6e}, "
        f"R={RADIUS_TEST_VALUES[best_index]:.3f}, "
        f"absolute RMS={absolute_rms[best_index]:.6e}"
    )
"""
        ),
        markdown(
            r"""
### How to read this test

The derivative curves are not additional fitting residuals.  They answer the
stronger question: do coefficients inferred from the metric values reproduce
the radial variation of independent metric/STF blocks?

- The V1 and C3 curves test radial closure even in the blocks used for the
  value fit.
- The vector and \(G^{Ti}_{1}\) curves test blocks omitted from the fit in
  both component/STF structure and radial dependence.
- \(G^{ij}_{0}\) remains an absolute-only diagnostic because its small
  first-order signal makes a fractional score misleading.

Agreement of the \(\Phi\) and finite-difference derivatives is a prerequisite
for interpreting these curves.  A small analytic-stencil self-error only
checks the differentiation of the exact response functions.

![Held-out radial-derivative closure](q8_t950_first_order_radial_derivative_closure.png)
"""
        ),
        code(
            r"""
def print_best_closure(
        sector_name, labels, results):
    print(f"\n{sector_name}")
    for estimator_index, label in enumerate(labels):
        values = results["held_out_residuals"][
            estimator_index]
        best_index = int(np.argmin(values))
        print(
            f"{label:<42s} "
            f"best closure={values[best_index]:.6e} "
            f"at R/M_B={RADIUS_OVER_SMALL_MASS[best_index]:.3f}; "
            f"range=[{values.min():.6e}, {values.max():.6e}]"
        )


print_best_closure(
    "Vector-sector closure",
    VECTOR_LABELS,
    VECTOR_RESULTS,
)
print_best_closure(
    "Clock/strain-sector closure",
    CLOCK_STRAIN_LABELS,
    CLOCK_STRAIN_RESULTS,
)

np.savez_compressed(
    OUTPUT_DIRECTORY / "q8_t950_first_order_results.npz",
    test_time=TEST_TIME,
    observation_id=TEST_OBSERVATION_ID,
    mass=SMALL_BLACK_HOLE_MASS,
    ah_center=AH_B_CENTER,
    curvature_center=CURVATURE_CENTER,
    radii=RADIUS_TEST_VALUES,
    radius_over_small_mass=RADIUS_OVER_SMALL_MASS,
    center_corrections_by_radius=CENTER_CORRECTIONS_BY_RADIUS,
    vector_labels=np.asarray(VECTOR_LABELS),
    clock_strain_labels=np.asarray(CLOCK_STRAIN_LABELS),
    vector_estimates=VECTOR_RESULTS["estimates"],
    clock_strain_estimates=CLOCK_STRAIN_RESULTS["estimates"],
    vector_condition_numbers=VECTOR_RESULTS["condition_numbers"],
    clock_strain_condition_numbers=(
        CLOCK_STRAIN_RESULTS["condition_numbers"]),
    vector_fit_residuals=VECTOR_RESULTS["fitted_block_residuals"],
    clock_strain_fit_residuals=(
        CLOCK_STRAIN_RESULTS["fitted_block_residuals"]),
    vector_held_out_residuals=(
        VECTOR_RESULTS["held_out_residuals"]),
    clock_strain_held_out_residuals=(
        CLOCK_STRAIN_RESULTS["held_out_residuals"]),
    vector_held_out_rms=VECTOR_RESULTS["held_out_rms"],
    clock_strain_held_out_rms=(
        CLOCK_STRAIN_RESULTS["held_out_rms"]),
    interior_random_seed=INTERIOR_RANDOM_SEED,
    interior_radial_edges=INTERIOR_RADIAL_EDGES,
    interior_local_points=INTERIOR_LOCAL_POINTS,
    interior_shell_index=INTERIOR_SHELL_INDEX,
    interior_zeroth_order_error=(
        INTERIOR_ZEROTH_ORDER_ERROR),
    interior_first_order_parameters=(
        INTERIOR_FIRST_ORDER_PARAMETERS),
    interior_total_error=INTERIOR_TOTAL_ERROR,
    interior_perturbation_residual=(
        INTERIOR_PERTURBATION_RESIDUAL),
    interior_shell_error=INTERIOR_SHELL_ERROR,
    selected_interior_pair=np.asarray(
        SELECTED_INTERIOR_PAIR),
    global_interior_minimum_index=np.asarray(
        GLOBAL_INTERIOR_MINIMUM_INDEX),
    surface_parameter_response_rms=(
        SURFACE_PARAMETER_RESPONSE_RMS),
    block_smallest_nonzero_singular_value=(
        BLOCK_SMALLEST_NONZERO_SINGULAR_VALUE),
    block_active_parameter_count=(
        BLOCK_ACTIVE_PARAMETER_COUNT),
    vector_smallest_singular_value=(
        VECTOR_SMALLEST_SINGULAR_VALUE),
    clock_smallest_singular_value=(
        CLOCK_SMALLEST_SINGULAR_VALUE),
    vector_interior_amplification=(
        VECTOR_INTERIOR_AMPLIFICATION),
    clock_interior_amplification=(
        CLOCK_INTERIOR_AMPLIFICATION),
    parameter_ensemble_mean=PARAMETER_ENSEMBLE_MEAN,
    parameter_ensemble_spread=(
        PARAMETER_ENSEMBLE_SPREAD),
    interior_parameter_contribution=(
        INTERIOR_PARAMETER_CONTRIBUTION),
    interior_parameter_uncertainty_impact=(
        INTERIOR_PARAMETER_UNCERTAINTY_IMPACT),
    interior_group_prediction_spread=(
        INTERIOR_GROUP_PREDICTION_SPREAD),
    interior_total_prediction_spread=(
        INTERIOR_TOTAL_PREDICTION_SPREAD),
    provisional_parameters=PROVISIONAL_PARAMETERS,
    provisional_vector_closure_fraction=(
        PROVISIONAL_VECTOR_CLOSURE_FRACTION),
    provisional_vector_closure_rms=(
        PROVISIONAL_VECTOR_CLOSURE_RMS),
    provisional_gti2_closure_fraction=(
        PROVISIONAL_GTI2_CLOSURE_FRACTION),
    provisional_gti2_closure_rms=(
        PROVISIONAL_GTI2_CLOSURE_RMS),
    provisional_gij1_closure_fraction=(
        PROVISIONAL_GIJ1_CLOSURE_FRACTION),
    provisional_gij1_closure_rms=(
        PROVISIONAL_GIJ1_CLOSURE_RMS),
    provisional_clock_closure_fraction=(
        PROVISIONAL_CLOCK_CLOSURE_FRACTION),
    provisional_clock_closure_rms=(
        PROVISIONAL_CLOCK_CLOSURE_RMS),
    provisional_gij0_diagnostic_fraction=(
        PROVISIONAL_GIJ0_DIAGNOSTIC_FRACTION),
    provisional_gij0_diagnostic_rms=(
        PROVISIONAL_GIJ0_DIAGNOSTIC_RMS),
    provisional_interior_total_error=(
        PROVISIONAL_INTERIOR_TOTAL_ERROR),
    provisional_interior_perturbation_residual=(
        PROVISIONAL_INTERIOR_PERTURBATION_RESIDUAL),
    provisional_interior_shell_error=(
        PROVISIONAL_INTERIOR_SHELL_ERROR),
    provisional_best_interior_radius_index=(
        PROVISIONAL_BEST_INTERIOR_RADIUS_INDEX),
    provisional_pointwise_relative_error=(
        PROVISIONAL_POINTWISE_RELATIVE_ERROR),
    zeroth_pointwise_relative_error=(
        ZEROTH_POINTWISE_RELATIVE_ERROR),
    phi_finite_difference_relative_error=(
        PHI_FINITE_DIFFERENCE_RELATIVE_ERROR),
    analytic_derivative_stencil_error=(
        ANALYTIC_DERIVATIVE_STENCIL_ERROR),
    derivative_v1_fraction=(
        DERIVATIVE_CLOSURES["V1 fitted blocks"][0]),
    derivative_v1_rms=(
        DERIVATIVE_CLOSURES["V1 fitted blocks"][1]),
    derivative_c3_fraction=(
        DERIVATIVE_CLOSURES["C3 fitted blocks"][0]),
    derivative_c3_rms=(
        DERIVATIVE_CLOSURES["C3 fitted blocks"][1]),
    derivative_vector_closure_fraction=(
        DERIVATIVE_CLOSURES["vector closure"][0]),
    derivative_vector_closure_rms=(
        DERIVATIVE_CLOSURES["vector closure"][1]),
    derivative_clock_closure_fraction=(
        DERIVATIVE_CLOSURES[r"$G^{Ti}_{1}$"][0]),
    derivative_clock_closure_rms=(
        DERIVATIVE_CLOSURES[r"$G^{Ti}_{1}$"][1]),
    derivative_gij0_diagnostic_fraction=(
        DERIVATIVE_CLOSURES[
            r"$G^{ij}_{0}$ diagnostic"][0]),
    derivative_gij0_diagnostic_rms=(
        DERIVATIVE_CLOSURES[
            r"$G^{ij}_{0}$ diagnostic"][1]),
)
print(
    "\nSaved compact results to",
    OUTPUT_DIRECTORY / "q8_t950_first_order_results.npz",
)
"""
        ),
        markdown(
            r"""
The notebook does not promote one block choice to the matching prescription.
The highlighted interior pair is selected only as a diagnostic: it is the
pair with the lowest median interior perturbation residual over this radius
scan.  A final choice should use surface closure, radius convergence,
conditioning, and interior prediction together, and should be checked on
additional time slices.
"""
        ),
    ]

    output = OUTPUT_DIRECTORY / "q8_first_order_matching.ipynb"
    nbformat.write(notebook, output)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
