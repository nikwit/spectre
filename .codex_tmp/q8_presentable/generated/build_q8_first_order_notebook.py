#!/usr/bin/env python
"""Build the q=8 hierarchical harmonic-worldtube matching notebook."""

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
# q=8 hierarchical harmonic-worldtube matching

This notebook is the executable companion to
[`q8_harmonic_worldtube_matching_report.pdf`](q8_harmonic_worldtube_matching_report.pdf).
It analyzes one full-resolution q=8 volume-data slice at
$T=950M_{\rm tot}$ and develops the matching calculation from first-order
minimal blocks through selected second-order subsystems and alternating
first-/second-order solves.

The central question throughout is:

> Do coefficients inferred from a restricted set of metric/STF blocks predict
> independent blocks, radial derivatives, and inverse-metric values inside the
> worldtube?

The small-black-hole mass is fixed to its initial-data value $M_B=1/9$.
The bookkeeping parameter $\epsilon$ labels perturbative order and is set to
one only when the complete truncated metric is evaluated.  It is not a
numerical velocity.

A common operational center is obtained from the Kretschmann dipole and then
held fixed across all fitting radii.  This is not an assertion that the
Kretschmann center, apparent-horizon areal center, and physical singularity
position coincide at finite perturbative order.

The response functions are imported from the audited q=4 package.  In
particular, the spatial-strain response is evaluated at fixed final/inertial
points and includes transport of the spatially varying Schwarzschild
background.
"""
        ),
        markdown(
            r"""
## How to navigate this notebook

The calculation is organized into four phases:

1. **Setup and first order (Sections 1--8).** Construct the common center,
   enumerate every minimal first-order estimator, test omitted STF blocks,
   propagate estimator ambiguity into the interior metric, and select the
   provisional V1+C3 hierarchy.
2. **Spatial derivatives (Sections 9--10).** Use the evolved
   $\Phi_{iab}$ field first as held-out radial closure and then as an
   additional equation in an overdetermined value-plus-derivative fit.
3. **Second order (Sections 11--14).** Subtract the fixed first-order model
   and known quadratic source, test the complete 43-column response space,
   isolate recoverable STF subsystems, and add radial derivatives.
4. **Cross-order coupling (Sections 15--16).** Alternate the selected
   first- and second-order solves using block Gauss--Seidel updates and test
   every iteration against independent interior points.

Four residuals appear repeatedly:

| name | data used | role |
|---|---|---|
| **training residual** | blocks included in a solve | verifies that the selected equations can be fitted |
| **surface closure** | component/STF blocks omitted from the solve | tests angular and component consistency |
| **radial closure** | $R n^i\partial_iG^{ab}$ from $\Phi_{iab}$ | tests an independent radial observable |
| **interior error** | 768 fixed random points inside the smallest sphere | primary out-of-sample metric test |

The notebook distinguishes three statuses:

- **algebraic validation:** manufactured data are recovered to roundoff;
- **numerical validation:** a fitted subsystem predicts held-out simulation
  data;
- **provisional interpretation:** a trend is useful but has not yet been
  demonstrated across enough radii, times, and simulations.

The single-slice notebook does not contain the later ten-slice time-derivative
analysis.  A summary and links to those companion calculations are included
after Section 16.
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

Q8_DIRECTORY = Path(os.environ.get(
    "Q8DIR", "/Users/niko/caltech/simulations/harmonic-sims/q8"))
sys.path.insert(0, str(Q8_DIRECTORY))

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
from second_order_responses import (
    FULL_GROUP_SLICES,
    FULL_PARAMETER_NAMES,
    STF_BASIS,
    STAGE1_INDICES,
    STAGE1_GROUP_SLICES,
    STAGE1_PARAMETER_NAMES,
    STAGE2_INDICES,
    STAGE3_INDICES,
    electric_response_at_direction,
    full_response_columns_at_points,
    known_quadratic_first_order_source,
    matrix_to_metric_tuple,
    spatial_quadrupole_response_at_direction,
    stage1_response_columns_at_points,
    strain_rate_response_at_direction,
    time_quadrupole_response_at_direction,
)

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
## Perturbative hierarchy used below

After centering, write the inverse metric as

$$
G^{ab}
=G_{(0)}^{ab}
+\epsilon A_1^{ab}p_1
+\epsilon^2\left[
Q_{11}^{ab}(p_1,p_1)+A_2^{ab}p_2
\right]
+O(\epsilon^3).
$$

At each order, all lower-order contributions are moved to the data side.
The new coefficients then enter linearly:

$$
d_{(n),B}=A_{(n),BA}p_{(n)}^A+e_{(n),B}.
$$

The source subscript denotes its **total** perturbative order.  For example,
$Q_{11}(p_1,p_1)$ is constructed from first-order coefficients but belongs to
the second-order source because it is quadratic.

The angular STF decomposition block-diagonalizes the exact response by
rotation and parity.  It does not guarantee statistical independence on
finite-mass-ratio data: omitted higher orders can project into the same
blocks.  This is why every solve below is paired with omitted-block and
interior validation.

The hierarchy is linear only after lower orders are held fixed.  The known
quadratic source and the later first-/second-order Gauss--Seidel iteration are
the two places where cross-order nonlinearity enters explicitly.

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
        markdown(
            r"""
## 10. Joint value-and-radial-derivative fits

We now allow the derivative data to participate in the fit.  For one sector
and one fitting radius, the augmented system is

$$
\begin{pmatrix}
\mathcal R_B\\
w_D\,R\partial_R\mathcal R_B
\end{pmatrix}
\widehat p
\simeq
\begin{pmatrix}
d_B\\
w_D\,R\partial_Rd_B
\end{pmatrix}.
$$

Here \(B\) is exactly the block already used by the provisional hierarchy:
V1 in the vector sector and C3 in the clock/strain sector.  Thus this test
adds radial information without adding a different component or STF
projection.  The unknowns remain the same six vector and seven clock/strain
coefficients, so every nonzero \(w_D\) produces an overdetermined linear
least-squares system.

The dimensionless weight \(w_D\) controls the relative influence of metric
values and derivatives.  We scan it rather than selecting it from the fitted
residual.  For every weight we retain as held-out validation:

- the value and derivative of \(G^{Ti}_{2}+G^{ij}_{1}\) in the vector sector;
- the value and derivative of \(G^{Ti}_{1}\) in the clock/strain sector;
- all random interior inverse-metric values.

Consequently, any improvement shown below cannot be obtained merely by
driving the newly included derivative equations toward zero.
"""
        ),
        code(
            r"""
JOINT_DERIVATIVE_WEIGHTS = np.array([
    0.0, 0.1, 0.3, 1.0, 3.0, 10.0,
])
JOINT_WEIGHT_LABELS = tuple(
    "metric only" if weight == 0.0
    else rf"$w_D={weight:g}$"
    for weight in JOINT_DERIVATIVE_WEIGHTS
)
JOINT_VECTOR_FIT_BLOCKS = (
    ("GTT", 1), ("GTI", 0))
JOINT_CLOCK_FIT_BLOCKS = (
    ("GTT", 0), ("GIJ", 2))


def stacked_block_system(
        data_blocks, response_blocks, blocks,
        radius_index, parameter_indices):
    data = np.concatenate([
        data_blocks[block][radius_index]
        for block in blocks
    ])
    matrix = np.concatenate([
        response_blocks[block][radius_index][
            :, parameter_indices]
        for block in blocks
    ], axis=0)
    return data, matrix


def solve_joint_sector(
        fit_blocks, parameter_indices, derivative_weight,
        radius_index):
    value_data, value_matrix = stacked_block_system(
        BLOCK_DATA,
        BLOCK_RESPONSE,
        fit_blocks,
        radius_index,
        parameter_indices,
    )
    derivative_data, derivative_matrix = stacked_block_system(
        DERIVATIVE_BLOCK_DATA,
        DERIVATIVE_BLOCK_RESPONSE,
        fit_blocks,
        radius_index,
        parameter_indices,
    )
    augmented_data = np.concatenate((
        value_data,
        derivative_weight * derivative_data,
    ))
    augmented_matrix = np.concatenate((
        value_matrix,
        derivative_weight * derivative_matrix,
    ), axis=0)
    estimate, _, rank, singular_values = np.linalg.lstsq(
        augmented_matrix, augmented_data, rcond=None)
    if rank != len(parameter_indices):
        raise RuntimeError(
            "Joint value/derivative system lost rank: "
            f"{rank} != {len(parameter_indices)}")
    condition_number = (
        singular_values[0] / singular_values[-1])
    return estimate, condition_number


JOINT_PARAMETERS = np.zeros((
    len(JOINT_DERIVATIVE_WEIGHTS),
    number_of_radii,
    number_of_parameters,
))
JOINT_VECTOR_CONDITION_NUMBER = np.zeros((
    len(JOINT_DERIVATIVE_WEIGHTS), number_of_radii))
JOINT_CLOCK_CONDITION_NUMBER = np.zeros_like(
    JOINT_VECTOR_CONDITION_NUMBER)

for weight_index, derivative_weight in enumerate(
        JOINT_DERIVATIVE_WEIGHTS):
    for radius_index in range(number_of_radii):
        (
            JOINT_PARAMETERS[
                weight_index,
                radius_index,
                VECTOR_PARAMETER_INDICES],
            JOINT_VECTOR_CONDITION_NUMBER[
                weight_index, radius_index],
        ) = solve_joint_sector(
            JOINT_VECTOR_FIT_BLOCKS,
            VECTOR_PARAMETER_INDICES,
            derivative_weight,
            radius_index,
        )
        (
            JOINT_PARAMETERS[
                weight_index,
                radius_index,
                CLOCK_STRAIN_PARAMETER_INDICES],
            JOINT_CLOCK_CONDITION_NUMBER[
                weight_index, radius_index],
        ) = solve_joint_sector(
            JOINT_CLOCK_FIT_BLOCKS,
            CLOCK_STRAIN_PARAMETER_INDICES,
            derivative_weight,
            radius_index,
        )

np.testing.assert_allclose(
    JOINT_PARAMETERS[0],
    PROVISIONAL_PARAMETERS,
    rtol=2.0e-12,
    atol=2.0e-13,
)
maximum_metric_only_difference = np.max(np.abs(
    JOINT_PARAMETERS[0] - PROVISIONAL_PARAMETERS))
print(
    "metric-only row reproduces V1+C3 to",
    f"{maximum_metric_only_difference:.3e}",
)
"""
        ),
        code(
            r"""
def closure_for_parameter_grid(
        data_blocks, response_blocks, selected_blocks,
        parameter_grid):
    fractional = np.empty(parameter_grid.shape[:2])
    absolute_rms = np.empty(parameter_grid.shape[:2])
    for weight_index in range(parameter_grid.shape[0]):
        for radius_index in range(parameter_grid.shape[1]):
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
                - matrix @ parameter_grid[
                    weight_index, radius_index])
            fractional[weight_index, radius_index] = (
                np.linalg.norm(residual)
                / max(np.linalg.norm(data), 1.0e-300)
            )
            absolute_rms[weight_index, radius_index] = (
                np.linalg.norm(residual)
                / np.sqrt(len(residual))
            )
    return fractional, absolute_rms


(
    JOINT_VALUE_VECTOR_CLOSURE_FRACTION,
    JOINT_VALUE_VECTOR_CLOSURE_RMS,
) = closure_for_parameter_grid(
    BLOCK_DATA,
    BLOCK_RESPONSE,
    PROVISIONAL_VECTOR_CLOSURE_BLOCKS,
    JOINT_PARAMETERS,
)
(
    JOINT_DERIVATIVE_VECTOR_CLOSURE_FRACTION,
    JOINT_DERIVATIVE_VECTOR_CLOSURE_RMS,
) = closure_for_parameter_grid(
    DERIVATIVE_BLOCK_DATA,
    DERIVATIVE_BLOCK_RESPONSE,
    PROVISIONAL_VECTOR_CLOSURE_BLOCKS,
    JOINT_PARAMETERS,
)
(
    JOINT_VALUE_CLOCK_CLOSURE_FRACTION,
    JOINT_VALUE_CLOCK_CLOSURE_RMS,
) = closure_for_parameter_grid(
    BLOCK_DATA,
    BLOCK_RESPONSE,
    (("GTI", 1),),
    JOINT_PARAMETERS,
)
(
    JOINT_DERIVATIVE_CLOCK_CLOSURE_FRACTION,
    JOINT_DERIVATIVE_CLOCK_CLOSURE_RMS,
) = closure_for_parameter_grid(
    DERIVATIVE_BLOCK_DATA,
    DERIVATIVE_BLOCK_RESPONSE,
    (("GTI", 1),),
    JOINT_PARAMETERS,
)

JOINT_INTERIOR_TOTAL_ERROR = np.empty((
    len(JOINT_DERIVATIVE_WEIGHTS), number_of_radii))
JOINT_INTERIOR_PERTURBATION_RESIDUAL = np.empty_like(
    JOINT_INTERIOR_TOTAL_ERROR)
for weight_index in range(
        len(JOINT_DERIVATIVE_WEIGHTS)):
    for radius_index in range(number_of_radii):
        model_values = (
            INTERIOR_BACKGROUND_VALUES
            + np.einsum(
                "ncp,p->nc",
                INTERIOR_RESPONSE_DESIGN,
                JOINT_PARAMETERS[
                    weight_index, radius_index],
            )
        )
        residual = (
            INTERIOR_NUMERICAL_VALUES - model_values)
        JOINT_INTERIOR_TOTAL_ERROR[
            weight_index, radius_index] = (
                np.linalg.norm(residual)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )
        JOINT_INTERIOR_PERTURBATION_RESIDUAL[
            weight_index, radius_index] = (
                np.linalg.norm(residual)
                / np.linalg.norm(
                    INTERIOR_PERTURBATION_VALUES)
            )

np.testing.assert_allclose(
    JOINT_INTERIOR_TOTAL_ERROR[0],
    PROVISIONAL_INTERIOR_TOTAL_ERROR,
    rtol=2.0e-12,
    atol=2.0e-13,
)

metric_only_parameter_norm = np.linalg.norm(
    JOINT_PARAMETERS[0], axis=-1)
JOINT_PARAMETER_DISPLACEMENT = (
    np.linalg.norm(
        JOINT_PARAMETERS
        - JOINT_PARAMETERS[0][None, :, :],
        axis=-1,
    )
    / np.maximum(
        metric_only_parameter_norm[None, :],
        1.0e-300,
    )
)

# The two sectors are orthogonal, so they need not share a derivative weight.
# Assemble every vector-weight / clock-weight combination without refitting.
HYBRID_PARAMETERS = np.zeros((
    len(JOINT_DERIVATIVE_WEIGHTS),
    len(JOINT_DERIVATIVE_WEIGHTS),
    number_of_radii,
    number_of_parameters,
))
for vector_weight_index in range(
        len(JOINT_DERIVATIVE_WEIGHTS)):
    for clock_weight_index in range(
            len(JOINT_DERIVATIVE_WEIGHTS)):
        HYBRID_PARAMETERS[
            vector_weight_index,
            clock_weight_index,
            :,
            VECTOR_PARAMETER_INDICES,
        ] = JOINT_PARAMETERS[
            vector_weight_index,
            :,
            VECTOR_PARAMETER_INDICES,
        ]
        HYBRID_PARAMETERS[
            vector_weight_index,
            clock_weight_index,
            :,
            CLOCK_STRAIN_PARAMETER_INDICES,
        ] = JOINT_PARAMETERS[
            clock_weight_index,
            :,
            CLOCK_STRAIN_PARAMETER_INDICES,
        ]

HYBRID_INTERIOR_TOTAL_ERROR = np.empty((
    len(JOINT_DERIVATIVE_WEIGHTS),
    len(JOINT_DERIVATIVE_WEIGHTS),
    number_of_radii,
))
for vector_weight_index in range(
        len(JOINT_DERIVATIVE_WEIGHTS)):
    for clock_weight_index in range(
            len(JOINT_DERIVATIVE_WEIGHTS)):
        for radius_index in range(number_of_radii):
            model_values = (
                INTERIOR_BACKGROUND_VALUES
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_RESPONSE_DESIGN,
                    HYBRID_PARAMETERS[
                        vector_weight_index,
                        clock_weight_index,
                        radius_index,
                    ],
                )
            )
            HYBRID_INTERIOR_TOTAL_ERROR[
                vector_weight_index,
                clock_weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - model_values)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )

HYBRID_MINIMUM_INTERIOR_ERROR = np.min(
    HYBRID_INTERIOR_TOTAL_ERROR, axis=-1)
HYBRID_MINIMUM_RADIUS_INDEX = np.argmin(
    HYBRID_INTERIOR_TOTAL_ERROR, axis=-1)
HYBRID_MINIMUM_RADIUS = RADIUS_TEST_VALUES[
    HYBRID_MINIMUM_RADIUS_INDEX]
HYBRID_GLOBAL_INDEX = np.unravel_index(
    np.argmin(HYBRID_INTERIOR_TOTAL_ERROR),
    HYBRID_INTERIOR_TOTAL_ERROR.shape,
)
"""
        ),
        code(
            r"""
fig_joint, axes = plt.subplots(
    2, 3, figsize=(17.0, 9.5))

for weight_index, label in enumerate(
        JOINT_WEIGHT_LABELS):
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_INTERIOR_TOTAL_ERROR[weight_index],
        "o-",
        label=label,
    )
    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_VALUE_VECTOR_CLOSURE_FRACTION[
            weight_index],
        "o-",
        label=label,
    )
    axes[0, 2].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_DERIVATIVE_VECTOR_CLOSURE_FRACTION[
            weight_index],
        "o-",
        label=label,
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_PARAMETER_DISPLACEMENT[weight_index],
        "o-",
        label=label,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_VALUE_CLOCK_CLOSURE_FRACTION[
            weight_index],
        "o-",
        label=label,
    )
    axes[1, 2].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_DERIVATIVE_CLOCK_CLOSURE_FRACTION[
            weight_index],
        "o-",
        label=label,
    )

axes[0, 0].axhline(
    INTERIOR_ZEROTH_ORDER_ERROR,
    color="black",
    linestyle="--",
    label="zeroth order",
)
axes[0, 0].set_title("Held-out interior inverse metric")
axes[0, 0].set_ylabel("ten-component relative RMS")
axes[0, 1].set_title(
    "Held-out vector value closure")
axes[0, 1].set_ylabel("fractional residual")
axes[0, 2].set_title(
    "Held-out vector derivative closure")
axes[0, 2].set_ylabel("fractional residual")
axes[1, 0].set_title(
    "Change from metric-only parameters")
axes[1, 0].set_ylabel(
    r"$\|\hat p(w_D)-\hat p(0)\|/\|\hat p(0)\|$")
axes[1, 1].set_title(
    r"Held-out $G^{Ti}_{1}$ value closure")
axes[1, 1].set_ylabel("fractional residual")
axes[1, 2].set_title(
    r"Held-out $G^{Ti}_{1}$ derivative closure")
axes[1, 2].set_ylabel("fractional residual")

for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=7)
fig_joint.suptitle(
    "Joint V1+C3 value-and-radial-derivative fit")
fig_joint.tight_layout()
fig_joint.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_fit.pdf")
fig_joint.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_fit.png",
    dpi=180,
)
plt.show()


fig_joint_absolute, axes = plt.subplots(
    1, 2, figsize=(13.5, 4.8))
for weight_index, label in enumerate(
        JOINT_WEIGHT_LABELS):
    axes[0].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_DERIVATIVE_VECTOR_CLOSURE_RMS[
            weight_index],
        "o-",
        label=label,
    )
    axes[1].loglog(
        RADIUS_OVER_SMALL_MASS,
        JOINT_DERIVATIVE_CLOCK_CLOSURE_RMS[
            weight_index],
        "o-",
        label=label,
    )
axes[0].set_title(
    "Held-out vector derivative")
axes[1].set_title(
    r"Held-out $G^{Ti}_{1}$ derivative")
for axis in axes:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.set_ylabel("absolute RMS")
    axis.legend(fontsize=8)
fig_joint_absolute.suptitle(
    "Absolute derivative closure for the joint fit")
fig_joint_absolute.tight_layout()
fig_joint_absolute.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_absolute.pdf")
fig_joint_absolute.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_absolute.png",
    dpi=180,
)
plt.show()


fig_hybrid, axes = plt.subplots(
    1, 2, figsize=(13.5, 5.2))
weight_tick_labels = [
    f"{weight:g}" for weight in JOINT_DERIVATIVE_WEIGHTS]

image_error = axes[0].imshow(
    HYBRID_MINIMUM_INTERIOR_ERROR,
    origin="lower",
    aspect="auto",
)
axes[0].set_title(
    "Best held-out interior error")
axes[0].set_xlabel(r"clock/strain weight $w_D^{\rm clock}$")
axes[0].set_ylabel(r"vector weight $w_D^{\rm vector}$")
fig_hybrid.colorbar(
    image_error, ax=axes[0], label="relative RMS")

image_radius = axes[1].imshow(
    HYBRID_MINIMUM_RADIUS,
    origin="lower",
    aspect="auto",
)
axes[1].set_title(
    "Fitting radius of that minimum")
axes[1].set_xlabel(r"clock/strain weight $w_D^{\rm clock}$")
axes[1].set_ylabel(r"vector weight $w_D^{\rm vector}$")
fig_hybrid.colorbar(
    image_radius, ax=axes[1],
    label=r"$R/M_{\rm tot}$")

for axis in axes:
    axis.set_xticks(
        np.arange(len(JOINT_DERIVATIVE_WEIGHTS)),
        weight_tick_labels,
    )
    axis.set_yticks(
        np.arange(len(JOINT_DERIVATIVE_WEIGHTS)),
        weight_tick_labels,
    )
for row in range(len(JOINT_DERIVATIVE_WEIGHTS)):
    for column in range(len(JOINT_DERIVATIVE_WEIGHTS)):
        axes[0].text(
            column, row,
            f"{HYBRID_MINIMUM_INTERIOR_ERROR[row, column]:.4f}",
            ha="center", va="center",
            color="white"
            if HYBRID_MINIMUM_INTERIOR_ERROR[row, column]
            > np.median(HYBRID_MINIMUM_INTERIOR_ERROR)
            else "black",
            fontsize=8,
        )
        axes[1].text(
            column, row,
            f"{HYBRID_MINIMUM_RADIUS[row, column]:.2f}",
            ha="center", va="center",
            color="white"
            if HYBRID_MINIMUM_RADIUS[row, column]
            > np.median(HYBRID_MINIMUM_RADIUS)
            else "black",
            fontsize=8,
        )
fig_hybrid.suptitle(
    "Independent derivative weights for the two sectors")
fig_hybrid.tight_layout()
fig_hybrid.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_hybrid.pdf")
fig_hybrid.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_first_order_joint_derivative_hybrid.png",
    dpi=180,
)
plt.show()
"""
        ),
        code(
            r"""
print("\nJoint-fit validation summary")
for weight_index, (weight, label) in enumerate(zip(
        JOINT_DERIVATIVE_WEIGHTS,
        JOINT_WEIGHT_LABELS)):
    interior_index = int(np.argmin(
        JOINT_INTERIOR_TOTAL_ERROR[weight_index]))
    vector_value_index = int(np.argmin(
        JOINT_VALUE_VECTOR_CLOSURE_FRACTION[
            weight_index]))
    vector_derivative_index = int(np.argmin(
        JOINT_DERIVATIVE_VECTOR_CLOSURE_FRACTION[
            weight_index]))
    clock_value_index = int(np.argmin(
        JOINT_VALUE_CLOCK_CLOSURE_FRACTION[
            weight_index]))
    clock_derivative_index = int(np.argmin(
        JOINT_DERIVATIVE_CLOCK_CLOSURE_FRACTION[
            weight_index]))
    print(f"\n{label}")
    print(
        "  interior:",
        f"{JOINT_INTERIOR_TOTAL_ERROR[weight_index, interior_index]:.6e}",
        f"at R={RADIUS_TEST_VALUES[interior_index]:.3f}",
    )
    print(
        "  vector value / derivative:",
        f"{JOINT_VALUE_VECTOR_CLOSURE_FRACTION[weight_index, vector_value_index]:.6e}",
        "/",
        f"{JOINT_DERIVATIVE_VECTOR_CLOSURE_FRACTION[weight_index, vector_derivative_index]:.6e}",
    )
    print(
        "  clock value / derivative:",
        f"{JOINT_VALUE_CLOCK_CLOSURE_FRACTION[weight_index, clock_value_index]:.6e}",
        "/",
        f"{JOINT_DERIVATIVE_CLOCK_CLOSURE_FRACTION[weight_index, clock_derivative_index]:.6e}",
    )

(
    best_vector_weight_index,
    best_clock_weight_index,
    best_hybrid_radius_index,
) = HYBRID_GLOBAL_INDEX
print("\nBest sector-weight combination by held-out interior error")
print(
    "  vector weight:",
    JOINT_DERIVATIVE_WEIGHTS[
        best_vector_weight_index],
)
print(
    "  clock weight:",
    JOINT_DERIVATIVE_WEIGHTS[
        best_clock_weight_index],
)
print(
    "  fitting radius:",
    RADIUS_TEST_VALUES[
        best_hybrid_radius_index],
)
print(
    "  interior relative RMS:",
    HYBRID_INTERIOR_TOTAL_ERROR[
        HYBRID_GLOBAL_INDEX],
)
"""
        ),
        markdown(
            r"""
### Interpretation of the weight scan

The derivative weight should not be chosen from the residual of the augmented
V1+C3 equations, since increasing \(w_D\) directly changes that objective.
The useful evidence is instead the behavior of blocks and interior points
that never entered the joint fit.

A credible improvement would therefore have to satisfy several conditions at
once:

1. reduce, or at least preserve, the held-out interior metric error;
2. improve derivative closure without destroying value closure;
3. remain stable over a range of nearby radii and nearby \(w_D\);
4. avoid a large displacement of the inferred parameters.

The two plots retain absolute derivative residuals alongside fractional
closures, because a small fractional denominator can otherwise make a
radially decaying block appear to become less accurate.

![Joint value-and-derivative validation](q8_t950_first_order_joint_derivative_fit.png)

![Absolute joint derivative closure](q8_t950_first_order_joint_derivative_absolute.png)

Because the vector and clock/strain systems are exactly orthogonal at first
order, the final heat map also permits different \(w_D\) values in the two
sectors.  This distinguishes a failure of radial derivatives in general from
a failure confined to one response block.

![Independent sector weights](q8_t950_first_order_joint_derivative_hybrid.png)

For this \(t=950M_{\rm tot}\) slice the distinction is clear.  Increasing the
vector derivative weight improves both omitted vector closures.  In contrast,
adding the C3 derivative generally moves the clock/strain parameters away
from the omitted \(G^{Ti}_{1}\) value and derivative.  The best entry in the
finite scan therefore uses

$$
w_D^{\rm vector}=3,\qquad w_D^{\rm clock}=0,
$$

at \(R=0.2M_{\rm tot}\).  Its held-out interior error is \(0.0356993\),
compared with \(0.0360554\) for the metric-only V1+C3 fit.  The improvement is
real but only about one percent, so this is evidence that the spatial
derivative can modestly regularize the vector sector—not yet evidence for a
universal optimal weight.  On this slice the clock/strain derivative should
remain a closure diagnostic rather than a fitted equation.
"""
        ),
        markdown(
            r"""
## 11. Controlled second-order start: the spatial 29-variable block

We now keep two first-order baselines in parallel:

1. metric-only V1+C3;
2. derivative-assisted V1 with \(w_D^{\rm vector}=3\), together with
   metric-only C3.

At every radius we construct the candidate second-order residual

$$
r_{(2)}^{ab}
=G_{\rm num}^{ab}
-G_{(0)}^{ab}
-\mathcal R_{(1)}^{ab}\widehat p_{(1)}
-Q_{11}^{ab}[\widehat p_{(1)},\widehat p_{(1)}].
$$

The known quadratic source used here is the audited
\(\beta^2+\beta\Lambda+\Lambda^2\) contribution.  It is evaluated with the
baseline first-order time gradient \(\beta_i\) and symmetric spatial strain
\(\Lambda_{ij}\).  Since \(\epsilon\) is an order-counting device, its numerical
value is one throughout this comparison.

We begin only with the 29-variable spatial block visible in \(G^{ij}\):

| group | number |
|---|---:|
| spatial acceleration \(\ddot q^i\) | 3 |
| symmetric strain rate \(\dot\Lambda_{(ij)}\) | 6 |
| physical electric STF tensor \(\mathcal E_{ij}\) | 5 |
| spatial harmonic quadrupoles \(Q^i{}_{\langle jk\rangle}\) | 15 |

Two solves are compared:

* **frozen first order:** fit these 29 coefficients to \(G^{ij}\);
* **first-order correction:** fit the same coefficients together with
  \(\delta\dot q^i\) and \(\delta\Lambda_{ij}\), the nine first-order
  corrections that are actually visible in \(G^{ij}\).

Thus the second solve implements

$$
r_{(2)}^{ij}
\simeq
\mathcal R_{(1)}^{ij}\delta p_{(1)}
+\mathcal R_{(2),A}^{ij}p_{(2)}^A.
$$

It is still a linear overdetermined problem.  The correction is diagnostic:
if it becomes large or destroys radius convergence, the apparent
second-order coefficients are likely repairing an imperfect first-order
baseline.

This section is deliberately not a full second-order reconstruction.  The
nine temporal-even and five magnetic coefficients are not yet fitted, so
\(G^{TT}\) and \(G^{Ti}\) are not final closure blocks at this stage.
"""
        ),
        code(
            r"""
SECOND_ORDER_BASELINE_NAMES = (
    "metric-only V1+C3",
    "derivative-assisted V1 + metric C3",
)
vector_weight_three_index = int(np.flatnonzero(
    np.isclose(JOINT_DERIVATIVE_WEIGHTS, 3.0))[0])
metric_only_weight_index = int(np.flatnonzero(
    np.isclose(JOINT_DERIVATIVE_WEIGHTS, 0.0))[0])
SECOND_ORDER_FIRST_ORDER_BASELINES = np.stack((
    PROVISIONAL_PARAMETERS,
    HYBRID_PARAMETERS[
        vector_weight_three_index,
        metric_only_weight_index,
    ],
))
SECOND_ORDER_SOLVE_NAMES = (
    "frozen first order",
    r"with $\delta p_{(1)}$",
)
VISIBLE_FIRST_ORDER_CORRECTION_INDICES = np.r_[
    np.arange(4, 7), np.arange(7, 13)]


def symmetric_matrix_from_six(values):
    matrix = np.zeros((3, 3))
    for value, (i, j) in zip(
            values, SPATIAL_COMPONENTS):
        matrix[i, j] = matrix[j, i] = value
    return matrix


def metric_linear_combination(columns, parameters):
    return tuple(
        sum(
            parameter * column[component]
            for parameter, column in zip(
                parameters, columns)
        )
        for component in range(3)
    )


def subtract_metric(left, *right_metrics):
    return tuple(
        left_component
        - sum(
            right[component]
            for right in right_metrics
        )
        for component, left_component in enumerate(left)
    )


def gij_values(metric):
    return np.column_stack([
        metric[2][:, i, j]
        for i, j in SPATIAL_COMPONENTS
    ]).ravel()


def scaled_least_squares(matrix, data):
    column_scales = np.linalg.norm(matrix, axis=0)
    if np.any(column_scales == 0.0):
        raise RuntimeError(
            "A fitted response column is identically zero")
    scaled_matrix = matrix / column_scales
    scaled_solution, _, rank, singular_values = (
        np.linalg.lstsq(
            scaled_matrix, data, rcond=1.0e-10)
    )
    solution = scaled_solution / column_scales
    condition_number = (
        singular_values[0] / singular_values[rank - 1])
    return solution, rank, condition_number


def response_design(columns, value_function):
    return np.column_stack([
        value_function(column) for column in columns
    ])


number_of_second_order_baselines = len(
    SECOND_ORDER_BASELINE_NAMES)
number_of_second_order_solves = len(
    SECOND_ORDER_SOLVE_NAMES)
number_of_stage1_parameters = len(
    STAGE1_PARAMETER_NAMES)

SURFACE_STAGE1_COLUMNS = []
SURFACE_STAGE1_GIJ_DESIGN = []
SURFACE_VISIBLE_FIRST_ORDER_GIJ_DESIGN = []
for radius in RADIUS_TEST_VALUES:
    points = radius * DIRECTIONS
    stage1_columns = stage1_response_columns_at_points(
        points, SMALL_BLACK_HOLE_MASS)
    first_order_columns = first_order_response_columns(
        points, SMALL_BLACK_HOLE_MASS)
    SURFACE_STAGE1_COLUMNS.append(stage1_columns)
    SURFACE_STAGE1_GIJ_DESIGN.append(
        response_design(stage1_columns, gij_values))
    SURFACE_VISIBLE_FIRST_ORDER_GIJ_DESIGN.append(
        response_design(
            tuple(
                first_order_columns[index]
                for index
                in VISIBLE_FIRST_ORDER_CORRECTION_INDICES
            ),
            gij_values,
        )
    )
"""
        ),
        code(
            r"""
SECOND_ORDER_CANDIDATE_RESIDUALS = [
    [None for _ in RADIUS_TEST_VALUES]
    for _ in SECOND_ORDER_BASELINE_NAMES
]
SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION = np.empty((
    number_of_second_order_baselines, number_of_radii))
SECOND_ORDER_QUADRATIC_RESIDUAL_FRACTION = np.empty_like(
    SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION)
SECOND_ORDER_QUADRATIC_SOURCE_FRACTION = np.empty_like(
    SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION)
SECOND_ORDER_STAGE1_PARAMETERS = np.zeros((
    number_of_second_order_baselines,
    number_of_second_order_solves,
    number_of_radii,
    number_of_stage1_parameters,
))
SECOND_ORDER_FIRST_ORDER_CORRECTIONS = np.zeros((
    number_of_second_order_baselines,
    number_of_second_order_solves,
    number_of_radii,
    number_of_parameters,
))
SECOND_ORDER_STAGE1_GIJ_RESIDUAL = np.empty((
    number_of_second_order_baselines,
    number_of_second_order_solves,
    number_of_radii,
))
SECOND_ORDER_STAGE1_CONDITION_NUMBER = np.empty_like(
    SECOND_ORDER_STAGE1_GIJ_RESIDUAL)
SECOND_ORDER_STAGE1_RANK = np.empty(
    SECOND_ORDER_STAGE1_GIJ_RESIDUAL.shape,
    dtype=int,
)
SECOND_ORDER_FIRST_ORDER_CORRECTION_FRACTION = np.zeros_like(
    SECOND_ORDER_STAGE1_GIJ_RESIDUAL)

for baseline_index in range(
        number_of_second_order_baselines):
    for radius_index, radius in enumerate(
            RADIUS_TEST_VALUES):
        points = radius * DIRECTIONS
        parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index]
        )
        numerical_metric = unpack_inverse_metric(
            FIRST_ORDER_METRIC_DATA[:, radius_index])
        background_metric = inverse_harmonic_schwarzschild(
            points, SMALL_BLACK_HOLE_MASS)
        raw_residual = subtract_metric(
            numerical_metric, background_metric)
        first_order_columns = first_order_response_columns(
            points, SMALL_BLACK_HOLE_MASS)
        first_order_prediction = metric_linear_combination(
            first_order_columns, parameters)
        first_order_residual = subtract_metric(
            raw_residual, first_order_prediction)
        beta = parameters[1:4]
        spatial_strain = symmetric_matrix_from_six(
            parameters[7:13])
        quadratic_source = (
            known_quadratic_first_order_source(
                points,
                beta,
                spatial_strain,
                SMALL_BLACK_HOLE_MASS,
            )
        )
        candidate_residual = subtract_metric(
            first_order_residual, quadratic_source)
        SECOND_ORDER_CANDIDATE_RESIDUALS[
            baseline_index][radius_index] = (
                candidate_residual)

        raw_norm = np.linalg.norm(
            metric_values_by_point(raw_residual))
        SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    metric_values_by_point(
                        first_order_residual))
                / raw_norm
            )
        SECOND_ORDER_QUADRATIC_RESIDUAL_FRACTION[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    metric_values_by_point(
                        candidate_residual))
                / raw_norm
            )
        SECOND_ORDER_QUADRATIC_SOURCE_FRACTION[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    metric_values_by_point(
                        quadratic_source))
                / raw_norm
            )

        stage1_matrix = (
            SURFACE_STAGE1_GIJ_DESIGN[radius_index])
        first_order_correction_matrix = (
            SURFACE_VISIBLE_FIRST_ORDER_GIJ_DESIGN[
                radius_index])
        candidate_gij = gij_values(candidate_residual)

        frozen_solution, frozen_rank, frozen_condition = (
            scaled_least_squares(
                stage1_matrix, candidate_gij))
        SECOND_ORDER_STAGE1_PARAMETERS[
            baseline_index, 0, radius_index] = (
                frozen_solution)
        SECOND_ORDER_STAGE1_RANK[
            baseline_index, 0, radius_index] = (
                frozen_rank)
        SECOND_ORDER_STAGE1_CONDITION_NUMBER[
            baseline_index, 0, radius_index] = (
                frozen_condition)

        augmented_matrix = np.column_stack((
            first_order_correction_matrix,
            stage1_matrix,
        ))
        (
            augmented_solution,
            augmented_rank,
            augmented_condition,
        ) = scaled_least_squares(
            augmented_matrix, candidate_gij)
        visible_correction = augmented_solution[:len(
            VISIBLE_FIRST_ORDER_CORRECTION_INDICES)]
        SECOND_ORDER_FIRST_ORDER_CORRECTIONS[
            baseline_index,
            1,
            radius_index,
            VISIBLE_FIRST_ORDER_CORRECTION_INDICES,
        ] = visible_correction
        SECOND_ORDER_STAGE1_PARAMETERS[
            baseline_index, 1, radius_index] = (
                augmented_solution[len(
                    VISIBLE_FIRST_ORDER_CORRECTION_INDICES):]
            )
        SECOND_ORDER_STAGE1_RANK[
            baseline_index, 1, radius_index] = (
                augmented_rank)
        SECOND_ORDER_STAGE1_CONDITION_NUMBER[
            baseline_index, 1, radius_index] = (
                augmented_condition)
        SECOND_ORDER_FIRST_ORDER_CORRECTION_FRACTION[
            baseline_index, 1, radius_index] = (
                np.linalg.norm(visible_correction)
                / max(
                    np.linalg.norm(parameters[
                        VISIBLE_FIRST_ORDER_CORRECTION_INDICES]),
                    1.0e-300,
                )
            )

        for solve_index in range(
                number_of_second_order_solves):
            correction = (
                SECOND_ORDER_FIRST_ORDER_CORRECTIONS[
                    baseline_index,
                    solve_index,
                    radius_index,
                    VISIBLE_FIRST_ORDER_CORRECTION_INDICES,
                ]
            )
            stage1_parameters = (
                SECOND_ORDER_STAGE1_PARAMETERS[
                    baseline_index,
                    solve_index,
                    radius_index,
                ]
            )
            gij_remainder = (
                candidate_gij
                - first_order_correction_matrix @ correction
                - stage1_matrix @ stage1_parameters
            )
            SECOND_ORDER_STAGE1_GIJ_RESIDUAL[
                baseline_index,
                solve_index,
                radius_index,
            ] = (
                np.linalg.norm(gij_remainder)
                / max(
                    np.linalg.norm(candidate_gij),
                    1.0e-300,
                )
            )

print("Second-order Stage-1 rank ranges")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    for solve_index, solve_name in enumerate(
            SECOND_ORDER_SOLVE_NAMES):
        print(
            f"  {baseline_name}; {solve_name}: "
            f"rank "
            f"{SECOND_ORDER_STAGE1_RANK[baseline_index, solve_index].min()}"
            f"--"
            f"{SECOND_ORDER_STAGE1_RANK[baseline_index, solve_index].max()}, "
            f"condition "
            f"{SECOND_ORDER_STAGE1_CONDITION_NUMBER[baseline_index, solve_index].min():.3e}"
            f"--"
            f"{SECOND_ORDER_STAGE1_CONDITION_NUMBER[baseline_index, solve_index].max():.3e}"
        )
"""
        ),
        markdown(
            r"""
### Manufactured recovery before interpreting the numerical fit

The frozen system must recover 29 arbitrary coefficients, and the augmented
system must recover nine arbitrary first-order corrections plus the same 29
second-order coefficients.  We test both systems on every radius using the
actual Spherepack grid and response matrices.  This is an implementation and
rank test; it does not test the physical adequacy of the truncated model.
"""
        ),
        code(
            r"""
stage1_manufactured_rng = np.random.default_rng(20260728)
MANUFACTURED_STAGE1_PARAMETERS = (
    stage1_manufactured_rng.normal(
        size=number_of_stage1_parameters))
MANUFACTURED_FIRST_ORDER_CORRECTIONS = (
    stage1_manufactured_rng.normal(
        size=len(
            VISIBLE_FIRST_ORDER_CORRECTION_INDICES)))
MANUFACTURED_FROZEN_ERROR = []
MANUFACTURED_AUGMENTED_ERROR = []

for radius_index in range(number_of_radii):
    stage1_matrix = (
        SURFACE_STAGE1_GIJ_DESIGN[radius_index])
    correction_matrix = (
        SURFACE_VISIBLE_FIRST_ORDER_GIJ_DESIGN[
            radius_index])

    frozen_data = (
        stage1_matrix @ MANUFACTURED_STAGE1_PARAMETERS)
    frozen_solution, _, _ = scaled_least_squares(
        stage1_matrix, frozen_data)
    MANUFACTURED_FROZEN_ERROR.append(
        np.max(np.abs(
            frozen_solution
            - MANUFACTURED_STAGE1_PARAMETERS))
    )

    augmented_matrix = np.column_stack((
        correction_matrix, stage1_matrix))
    prescribed_augmented = np.concatenate((
        MANUFACTURED_FIRST_ORDER_CORRECTIONS,
        MANUFACTURED_STAGE1_PARAMETERS,
    ))
    augmented_data = (
        augmented_matrix @ prescribed_augmented)
    augmented_solution, _, _ = scaled_least_squares(
        augmented_matrix, augmented_data)
    MANUFACTURED_AUGMENTED_ERROR.append(
        np.max(np.abs(
            augmented_solution
            - prescribed_augmented))
    )

MANUFACTURED_FROZEN_ERROR = np.asarray(
    MANUFACTURED_FROZEN_ERROR)
MANUFACTURED_AUGMENTED_ERROR = np.asarray(
    MANUFACTURED_AUGMENTED_ERROR)
print(
    "maximum manufactured frozen error:",
    f"{MANUFACTURED_FROZEN_ERROR.max():.3e}",
)
print(
    "maximum manufactured augmented error:",
    f"{MANUFACTURED_AUGMENTED_ERROR.max():.3e}",
)
assert MANUFACTURED_FROZEN_ERROR.max() < 2.0e-10
assert MANUFACTURED_AUGMENTED_ERROR.max() < 2.0e-9
"""
        ),
        code(
            r"""
# Evaluate the partial Stage-1 model at the same random interior points used
# for first-order validation.  The temporal and magnetic second-order
# responses are intentionally absent.
INTERIOR_STAGE1_COLUMNS = stage1_response_columns_at_points(
    INTERIOR_LOCAL_POINTS, SMALL_BLACK_HOLE_MASS)
INTERIOR_STAGE1_DESIGN = np.stack([
    metric_values_by_point(column)
    for column in INTERIOR_STAGE1_COLUMNS
], axis=-1)

SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR = np.empty((
    number_of_second_order_baselines, number_of_radii))
SECOND_ORDER_STAGE1_INTERIOR_ERROR = np.empty((
    number_of_second_order_baselines,
    number_of_second_order_solves,
    number_of_radii,
))

for baseline_index in range(
        number_of_second_order_baselines):
    for radius_index in range(number_of_radii):
        baseline_parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index])
        first_order_model = (
            INTERIOR_BACKGROUND_VALUES
            + np.einsum(
                "ncp,p->nc",
                INTERIOR_RESPONSE_DESIGN,
                baseline_parameters,
            )
        )
        SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - first_order_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )

        beta = baseline_parameters[1:4]
        spatial_strain = symmetric_matrix_from_six(
            baseline_parameters[7:13])
        interior_quadratic_source = (
            known_quadratic_first_order_source(
                INTERIOR_LOCAL_POINTS,
                beta,
                spatial_strain,
                SMALL_BLACK_HOLE_MASS,
            )
        )
        interior_quadratic_values = metric_values_by_point(
            interior_quadratic_source)

        for solve_index in range(
                number_of_second_order_solves):
            corrected_first_order_parameters = (
                baseline_parameters
                + SECOND_ORDER_FIRST_ORDER_CORRECTIONS[
                    baseline_index,
                    solve_index,
                    radius_index,
                ]
            )
            model_values = (
                INTERIOR_BACKGROUND_VALUES
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_RESPONSE_DESIGN,
                    corrected_first_order_parameters,
                )
                + interior_quadratic_values
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_STAGE1_DESIGN,
                    SECOND_ORDER_STAGE1_PARAMETERS[
                        baseline_index,
                        solve_index,
                        radius_index,
                    ],
                )
            )
            SECOND_ORDER_STAGE1_INTERIOR_ERROR[
                baseline_index,
                solve_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - model_values)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )
"""
        ),
        code(
            r"""
fig_second_order_stage1, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))
baseline_colors = ("C0", "C1")
solve_markers = ("o", "s")

for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION[
            baseline_index],
        "o-",
        color=baseline_colors[baseline_index],
        label=baseline_name + r": after $G_{(1)}$",
    )
    axes[0, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_QUADRATIC_RESIDUAL_FRACTION[
            baseline_index],
        "s--",
        color=baseline_colors[baseline_index],
        label=baseline_name + r": after $G_{(1)}+Q_{11}$",
    )

    for solve_index, solve_name in enumerate(
            SECOND_ORDER_SOLVE_NAMES):
        axes[0, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_STAGE1_GIJ_RESIDUAL[
                baseline_index, solve_index],
            marker=solve_markers[solve_index],
            linestyle=(
                "-" if solve_index == 0 else "--"),
            color=baseline_colors[baseline_index],
            label=baseline_name + "; " + solve_name,
        )
        axes[1, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_STAGE1_INTERIOR_ERROR[
                baseline_index, solve_index],
            marker=solve_markers[solve_index],
            linestyle=(
                "-" if solve_index == 0 else "--"),
            color=baseline_colors[baseline_index],
            label=baseline_name + "; " + solve_name,
        )

    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_FIRST_ORDER_CORRECTION_FRACTION[
            baseline_index, 1],
        "o-",
        color=baseline_colors[baseline_index],
        label=baseline_name,
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR[
            baseline_index],
        ":",
        linewidth=2.0,
        color=baseline_colors[baseline_index],
        label=baseline_name + ": first order",
    )

axes[0, 0].set_title(
    "Candidate second-order residual")
axes[0, 0].set_ylabel(
    r"norm divided by $\|G_{\rm num}-G_{(0)}\|$")
axes[0, 1].set_title(
    r"Stage-1 fit residual in $G^{ij}$")
axes[0, 1].set_ylabel("fractional residual")
axes[1, 0].set_title(
    "First-order correction admitted by Stage 1")
axes[1, 0].set_ylabel(
    r"$\|\delta p_{(1)}^{G^{ij}}\|/"
    r"\|p_{(1)}^{G^{ij}}\|$")
axes[1, 1].set_title(
    "Held-out interior metric: partial Stage-1 model")
axes[1, 1].set_ylabel(
    "ten-component relative RMS")
for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=7)
fig_second_order_stage1.suptitle(
    "First controlled q=8 second-order block")
fig_second_order_stage1.tight_layout()
fig_second_order_stage1.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_stage1_validation.pdf")
fig_second_order_stage1.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_stage1_validation.png",
    dpi=180,
)
plt.show()


fig_second_order_groups, axes = plt.subplots(
    2, 2, figsize=(13.5, 9.0))
for axis, (group_name, group_slice) in zip(
        axes.flat, STAGE1_GROUP_SLICES.items()):
    for baseline_index, baseline_name in enumerate(
            SECOND_ORDER_BASELINE_NAMES):
        for solve_index, solve_name in enumerate(
                SECOND_ORDER_SOLVE_NAMES):
            axis.loglog(
                RADIUS_OVER_SMALL_MASS,
                np.linalg.norm(
                    SECOND_ORDER_STAGE1_PARAMETERS[
                        baseline_index,
                        solve_index,
                        :,
                        group_slice,
                    ],
                    axis=1,
                ),
                marker=solve_markers[solve_index],
                linestyle=(
                    "-" if solve_index == 0 else "--"),
                color=baseline_colors[baseline_index],
                label=baseline_name + "; " + solve_name,
            )
    axis.set_title(group_name)
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.set_ylabel("coefficient norm")
    axis.legend(fontsize=6)
fig_second_order_groups.suptitle(
    "Stage-1 coefficient groups and first-order leakage")
fig_second_order_groups.tight_layout()
fig_second_order_groups.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_stage1_parameters.pdf")
fig_second_order_groups.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_stage1_parameters.png",
    dpi=180,
)
plt.show()


print("\nSecond-order Stage-1 summary")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    print(f"\n{baseline_name}")
    for solve_index, solve_name in enumerate(
            SECOND_ORDER_SOLVE_NAMES):
        best_index = int(np.argmin(
            SECOND_ORDER_STAGE1_INTERIOR_ERROR[
                baseline_index, solve_index]))
        print(
            f"  {solve_name:<24s}: "
            f"interior="
            f"{SECOND_ORDER_STAGE1_INTERIOR_ERROR[baseline_index, solve_index, best_index]:.6e} "
            f"at R={RADIUS_TEST_VALUES[best_index]:.3f}; "
            f"GIJ fit="
            f"{SECOND_ORDER_STAGE1_GIJ_RESIDUAL[baseline_index, solve_index, best_index]:.6e}"
        )
"""
        ),
        markdown(
            r"""
### Reading the first second-order block

The \(G^{ij}\) residual can become very small because Stage 1 contains 29
independent columns.  That alone is not evidence of successful physical
recovery.  The stronger diagnostics are:

- convergence of each coefficient-group norm with radius;
- agreement between the two first-order baselines;
- the size of the admitted \(\delta p_{(1)}\);
- prediction at the random interior points.

The interior curve is intentionally labeled a *partial* Stage-1 model.  It
contains the known quadratic source and 29 fitted spatial coefficients, but
not the remaining temporal-even or magnetic second-order responses.  It is
therefore useful for detecting catastrophic overfitting or clear improvement,
not yet for declaring a complete second-order boundary model.

![Second-order Stage-1 validation](q8_t950_second_order_stage1_validation.png)

![Second-order Stage-1 parameter groups](q8_t950_second_order_stage1_parameters.png)

On this slice the frozen 29-column system is numerically healthy: its
normalized condition number is only \(2.4\)--\(2.9\), and manufactured
coefficients are recovered to \(9\times10^{-12}\).  The derivative-assisted
first-order baseline gives the better partial Stage-1 interior result,
\(0.04108\) at \(R=0.8M_{\rm tot}\), compared with \(0.04546\) for the
metric-only baseline.  Neither improves on the best first-order-only result
\(0.03570\).  This is not yet a failure of second order, because 14
second-order responses remain absent and may restore component cancellations.

The unconstrained first-order-correction solve is much less credible.  Although
it is algebraically full rank, its condition number grows as high as
\(2\times10^4\); it admits corrections comparable to or larger than the
original spatial first-order coefficients and worsens the best interior error
to about \(0.11\).  We should therefore keep the first-order baseline frozen
for the next temporal-even and magnetic stages, or introduce a quantitatively
motivated prior before allowing \(\delta p_{(1)}\) to vary.
"""
        ),
        markdown(
            r"""
## 12. Complete second-order triangular solve

We now add the fourteen responses omitted from Stage 1, without refitting any
first-order coefficient.  The 43 second-order unknowns are ordered as

\[
\underbrace{a^i,\ \dot b^{(ij)},\ {\cal E}_{ij},\
             Q^i{}_{\langle jk\rangle}}_{\text{29 spatial coefficients}},
\qquad
\underbrace{a^T,\ \dot b^T{}_i,\
             Q^T{}_{\langle ij\rangle}}_{\text{9 temporal-even coefficients}},
\qquad
\underbrace{{\cal B}_{ij}}_{\text{5 magnetic STF coefficients}} .
\]

The response structure makes a triangular component solve possible:

1. fit the 29 spatial coefficients to all six independent components of
   \(G^{ij}\);
2. subtract their response and fit the nine temporal-even coefficients to
   \(G^{TT}\);
3. subtract that response and fit the five magnetic coefficients to all three
   components of \(G^{Ti}\).

Each line is still an overdetermined linear least-squares problem.  Later
responses can affect components used in earlier stages only through structures
that vanish in those projection blocks.  We verify this operationally with a
manufactured 43-parameter recovery test before interpreting the simulation.
"""
        ),
        code(
            r"""
def gtt_values(metric):
    return np.asarray(metric[0]).ravel()


def gti_values(metric):
    return np.asarray(metric[1]).ravel()


number_of_full_second_order_parameters = len(
    FULL_PARAMETER_NAMES)
assert number_of_full_second_order_parameters == 43
assert np.array_equal(STAGE1_INDICES, np.arange(29))
assert np.array_equal(STAGE2_INDICES, np.arange(29, 38))
assert np.array_equal(STAGE3_INDICES, np.arange(38, 43))

SURFACE_FULL_SECOND_ORDER_COLUMNS = []
SURFACE_STAGE2_GTT_DESIGN = []
SURFACE_STAGE3_GTI_DESIGN = []
for radius_index, radius in enumerate(
        RADIUS_TEST_VALUES):
    columns = full_response_columns_at_points(
        radius * DIRECTIONS, SMALL_BLACK_HOLE_MASS)
    SURFACE_FULL_SECOND_ORDER_COLUMNS.append(columns)
    SURFACE_STAGE2_GTT_DESIGN.append(
        response_design(
            tuple(columns[index] for index in STAGE2_INDICES),
            gtt_values,
        )
    )
    SURFACE_STAGE3_GTI_DESIGN.append(
        response_design(
            tuple(columns[index] for index in STAGE3_INDICES),
            gti_values,
        )
    )

    # Guard against transcription drift between the original Stage-1
    # implementation and the complete response generator.
    for old_column, new_column in zip(
            SURFACE_STAGE1_COLUMNS[radius_index],
            columns[:len(STAGE1_PARAMETER_NAMES)]):
        for old_component, new_component in zip(
                old_column, new_column):
            np.testing.assert_allclose(
                old_component,
                new_component,
                rtol=2.0e-13,
                atol=2.0e-13,
            )
"""
        ),
        code(
            r"""
COMPLETE_SECOND_ORDER_STAGE_NAMES = (
    r"after $G_{(1)}+Q_{11}$",
    "after 29 spatial responses",
    "after 9 temporal-even responses",
    "after 5 magnetic responses",
)
number_of_complete_stages = len(
    COMPLETE_SECOND_ORDER_STAGE_NAMES)

COMPLETE_SECOND_ORDER_PARAMETERS = np.zeros((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_full_second_order_parameters,
))
COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_complete_stages,
))
COMPLETE_SECOND_ORDER_STAGE_RAW_FRACTION = np.empty_like(
    COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION)
COMPLETE_SECOND_ORDER_COMPONENT_FRACTION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    3,
))
COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
))
COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE = np.empty_like(
    COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE)
COMPLETE_SECOND_ORDER_STAGE2_RANK = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
), dtype=int)
COMPLETE_SECOND_ORDER_STAGE3_RANK = np.empty_like(
    COMPLETE_SECOND_ORDER_STAGE2_RANK)
COMPLETE_SECOND_ORDER_STAGE2_CONDITION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
))
COMPLETE_SECOND_ORDER_STAGE3_CONDITION = np.empty_like(
    COMPLETE_SECOND_ORDER_STAGE2_CONDITION)

for baseline_index in range(
        number_of_second_order_baselines):
    for radius_index, radius in enumerate(
            RADIUS_TEST_VALUES):
        columns = SURFACE_FULL_SECOND_ORDER_COLUMNS[
            radius_index]
        candidate = SECOND_ORDER_CANDIDATE_RESIDUALS[
            baseline_index][radius_index]
        candidate_norm = max(
            np.linalg.norm(metric_values_by_point(candidate)),
            1.0e-300,
        )

        points = radius * DIRECTIONS
        numerical_metric = unpack_inverse_metric(
            FIRST_ORDER_METRIC_DATA[:, radius_index])
        background_metric = inverse_harmonic_schwarzschild(
            points, SMALL_BLACK_HOLE_MASS)
        raw_residual = subtract_metric(
            numerical_metric, background_metric)
        raw_norm = max(
            np.linalg.norm(
                metric_values_by_point(raw_residual)),
            1.0e-300,
        )

        parameters = COMPLETE_SECOND_ORDER_PARAMETERS[
            baseline_index, radius_index]
        parameters[STAGE1_INDICES] = (
            SECOND_ORDER_STAGE1_PARAMETERS[
                baseline_index, 0, radius_index]
        )
        stage1_prediction = metric_linear_combination(
            tuple(columns[index] for index in STAGE1_INDICES),
            parameters[STAGE1_INDICES],
        )
        after_stage1 = subtract_metric(
            candidate, stage1_prediction)

        (
            temporal_solution,
            temporal_rank,
            temporal_condition,
        ) = scaled_least_squares(
            SURFACE_STAGE2_GTT_DESIGN[radius_index],
            gtt_values(after_stage1),
        )
        parameters[STAGE2_INDICES] = temporal_solution
        COMPLETE_SECOND_ORDER_STAGE2_RANK[
            baseline_index, radius_index] = temporal_rank
        COMPLETE_SECOND_ORDER_STAGE2_CONDITION[
            baseline_index, radius_index] = temporal_condition
        stage2_prediction = metric_linear_combination(
            tuple(columns[index] for index in STAGE2_INDICES),
            temporal_solution,
        )
        after_stage2 = subtract_metric(
            after_stage1, stage2_prediction)

        (
            magnetic_solution,
            magnetic_rank,
            magnetic_condition,
        ) = scaled_least_squares(
            SURFACE_STAGE3_GTI_DESIGN[radius_index],
            gti_values(after_stage2),
        )
        parameters[STAGE3_INDICES] = magnetic_solution
        COMPLETE_SECOND_ORDER_STAGE3_RANK[
            baseline_index, radius_index] = magnetic_rank
        COMPLETE_SECOND_ORDER_STAGE3_CONDITION[
            baseline_index, radius_index] = magnetic_condition
        stage3_prediction = metric_linear_combination(
            tuple(columns[index] for index in STAGE3_INDICES),
            magnetic_solution,
        )
        after_stage3 = subtract_metric(
            after_stage2, stage3_prediction)
        COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE[
            baseline_index, radius_index] = (
                np.linalg.norm(after_stage2[0])
                / max(
                    np.linalg.norm(after_stage1[0]),
                    1.0e-300,
                )
            )
        COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE[
            baseline_index, radius_index] = (
                np.linalg.norm(after_stage3[1])
                / max(
                    np.linalg.norm(after_stage2[1]),
                    1.0e-300,
                )
            )

        stage_metrics = (
            candidate,
            after_stage1,
            after_stage2,
            after_stage3,
        )
        for stage_index, stage_metric in enumerate(
                stage_metrics):
            stage_norm = np.linalg.norm(
                metric_values_by_point(stage_metric))
            COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION[
                baseline_index,
                radius_index,
                stage_index,
            ] = stage_norm / candidate_norm
            COMPLETE_SECOND_ORDER_STAGE_RAW_FRACTION[
                baseline_index,
                radius_index,
                stage_index,
            ] = stage_norm / raw_norm

        for component_index, component in enumerate(
                after_stage3):
            COMPLETE_SECOND_ORDER_COMPONENT_FRACTION[
                baseline_index,
                radius_index,
                component_index,
            ] = (
                np.linalg.norm(component)
                / max(
                    np.linalg.norm(candidate[component_index]),
                    1.0e-300,
                )
            )

print("Complete second-order triangular ranks and conditions")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    print(f"\n{baseline_name}")
    print(
        "  temporal-even GTT block: "
        f"rank {COMPLETE_SECOND_ORDER_STAGE2_RANK[baseline_index].min()}"
        f"--{COMPLETE_SECOND_ORDER_STAGE2_RANK[baseline_index].max()}, "
        "condition "
        f"{COMPLETE_SECOND_ORDER_STAGE2_CONDITION[baseline_index].min():.3e}"
        "--"
        f"{COMPLETE_SECOND_ORDER_STAGE2_CONDITION[baseline_index].max():.3e}"
    )
    print(
        "  magnetic GTI block:      "
        f"rank {COMPLETE_SECOND_ORDER_STAGE3_RANK[baseline_index].min()}"
        f"--{COMPLETE_SECOND_ORDER_STAGE3_RANK[baseline_index].max()}, "
        "condition "
        f"{COMPLETE_SECOND_ORDER_STAGE3_CONDITION[baseline_index].min():.3e}"
        "--"
        f"{COMPLETE_SECOND_ORDER_STAGE3_CONDITION[baseline_index].max():.3e}"
    )
    print(
        "  temporal GTT closure:    "
        f"{COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE[baseline_index].min():.3e}"
        "--"
        f"{COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE[baseline_index].max():.3e}"
    )
    print(
        "  magnetic GTI closure:    "
        f"{COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE[baseline_index].min():.3e}"
        "--"
        f"{COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE[baseline_index].max():.3e}"
    )
"""
        ),
        markdown(
            r"""
### Manufactured recovery of the complete hierarchy

We prescribe all 43 coefficients, construct their exact response on each
worldtube, and apply precisely the same three solves.  This is stronger than a
rank printout: it checks that subtracting each fitted response leaves the
correct data for the next block and that the component ordering is consistent.
"""
        ),
        code(
            r"""
complete_manufactured_rng = np.random.default_rng(20260729)
MANUFACTURED_COMPLETE_PARAMETERS = (
    complete_manufactured_rng.normal(
        size=number_of_full_second_order_parameters))
MANUFACTURED_COMPLETE_PARAMETER_ERROR = np.empty(
    number_of_radii)
MANUFACTURED_COMPLETE_METRIC_CLOSURE = np.empty(
    number_of_radii)

for radius_index in range(number_of_radii):
    columns = SURFACE_FULL_SECOND_ORDER_COLUMNS[
        radius_index]
    synthetic_metric = metric_linear_combination(
        columns, MANUFACTURED_COMPLETE_PARAMETERS)
    recovered = np.zeros(
        number_of_full_second_order_parameters)

    recovered[STAGE1_INDICES], _, _ = (
        scaled_least_squares(
            SURFACE_STAGE1_GIJ_DESIGN[radius_index],
            gij_values(synthetic_metric),
        )
    )
    after_spatial = subtract_metric(
        synthetic_metric,
        metric_linear_combination(
            tuple(columns[index] for index in STAGE1_INDICES),
            recovered[STAGE1_INDICES],
        ),
    )
    recovered[STAGE2_INDICES], _, _ = (
        scaled_least_squares(
            SURFACE_STAGE2_GTT_DESIGN[radius_index],
            gtt_values(after_spatial),
        )
    )
    after_temporal = subtract_metric(
        after_spatial,
        metric_linear_combination(
            tuple(columns[index] for index in STAGE2_INDICES),
            recovered[STAGE2_INDICES],
        ),
    )
    recovered[STAGE3_INDICES], _, _ = (
        scaled_least_squares(
            SURFACE_STAGE3_GTI_DESIGN[radius_index],
            gti_values(after_temporal),
        )
    )
    reconstructed = metric_linear_combination(
        columns, recovered)

    MANUFACTURED_COMPLETE_PARAMETER_ERROR[radius_index] = (
        np.max(np.abs(
            recovered
            - MANUFACTURED_COMPLETE_PARAMETERS))
    )
    MANUFACTURED_COMPLETE_METRIC_CLOSURE[radius_index] = (
        np.linalg.norm(metric_values_by_point(
            subtract_metric(
                synthetic_metric, reconstructed)))
        / np.linalg.norm(
            metric_values_by_point(synthetic_metric))
    )

print(
    "maximum manufactured 43-parameter error:",
    f"{MANUFACTURED_COMPLETE_PARAMETER_ERROR.max():.3e}",
)
print(
    "maximum manufactured metric closure:",
    f"{MANUFACTURED_COMPLETE_METRIC_CLOSURE.max():.3e}",
)
assert MANUFACTURED_COMPLETE_PARAMETER_ERROR.max() < 3.0e-10
assert MANUFACTURED_COMPLETE_METRIC_CLOSURE.max() < 3.0e-12
"""
        ),
        code(
            r"""
# Complete held-out interior prediction.  The random points and numerical
# target are exactly the same as in the zeroth- and first-order comparisons.
INTERIOR_FULL_SECOND_ORDER_COLUMNS = (
    full_response_columns_at_points(
        INTERIOR_LOCAL_POINTS,
        SMALL_BLACK_HOLE_MASS,
    )
)
INTERIOR_FULL_SECOND_ORDER_DESIGN = np.stack([
    metric_values_by_point(column)
    for column in INTERIOR_FULL_SECOND_ORDER_COLUMNS
], axis=-1)

COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
))
COMPLETE_SECOND_ORDER_INTERIOR_ERROR = np.empty_like(
    COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR)
COMPLETE_SECOND_ORDER_INTERIOR_PERTURBATION_RESIDUAL = (
    np.empty_like(
        COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR)
)
COMPLETE_SECOND_ORDER_POINTWISE_RELATIVE_ERROR = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    len(INTERIOR_LOCAL_POINTS),
))

for baseline_index in range(
        number_of_second_order_baselines):
    for radius_index in range(number_of_radii):
        first_order_parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index])
        first_order_values = np.einsum(
            "ncp,p->nc",
            INTERIOR_RESPONSE_DESIGN,
            first_order_parameters,
        )
        beta = first_order_parameters[1:4]
        spatial_strain = symmetric_matrix_from_six(
            first_order_parameters[7:13])
        quadratic_values = metric_values_by_point(
            known_quadratic_first_order_source(
                INTERIOR_LOCAL_POINTS,
                beta,
                spatial_strain,
                SMALL_BLACK_HOLE_MASS,
            )
        )
        quadratic_model = (
            INTERIOR_BACKGROUND_VALUES
            + first_order_values
            + quadratic_values
        )
        COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - quadratic_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )

        complete_model = (
            quadratic_model
            + np.einsum(
                "ncp,p->nc",
                INTERIOR_FULL_SECOND_ORDER_DESIGN,
                COMPLETE_SECOND_ORDER_PARAMETERS[
                    baseline_index, radius_index],
            )
        )
        pointwise_error = np.linalg.norm(
            INTERIOR_NUMERICAL_VALUES
            - complete_model,
            axis=1,
        ) / np.maximum(
            np.linalg.norm(
                INTERIOR_NUMERICAL_VALUES, axis=1),
            1.0e-300,
        )
        COMPLETE_SECOND_ORDER_POINTWISE_RELATIVE_ERROR[
            baseline_index, radius_index] = pointwise_error
        COMPLETE_SECOND_ORDER_INTERIOR_ERROR[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - complete_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )
        COMPLETE_SECOND_ORDER_INTERIOR_PERTURBATION_RESIDUAL[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - complete_model)
                / np.linalg.norm(
                    INTERIOR_PERTURBATION_VALUES)
            )
"""
        ),
        code(
            r"""
fig_second_order_complete, axes = plt.subplots(
    2, 2, figsize=(14.5, 10.0))
stage_markers = ("o", "s", "^", "D")
component_labels = (
    r"$G^{TT}$", r"$G^{Ti}$", r"$G^{ij}$")

for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    for stage_index, stage_name in enumerate(
            COMPLETE_SECOND_ORDER_STAGE_NAMES):
        axes[0, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION[
                baseline_index, :, stage_index],
            marker=stage_markers[stage_index],
            linestyle=(
                "-" if baseline_index == 0 else "--"),
            color=f"C{stage_index}",
            label=baseline_name + "; " + stage_name,
        )

    for component_index, component_label in enumerate(
            component_labels):
        axes[0, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            COMPLETE_SECOND_ORDER_COMPONENT_FRACTION[
                baseline_index, :, component_index],
            marker=stage_markers[component_index],
            linestyle=(
                "-" if baseline_index == 0 else "--"),
            color=f"C{component_index}",
            label=baseline_name + "; " + component_label,
        )

    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR[
            baseline_index],
        "o:",
        color=baseline_colors[baseline_index],
        label=baseline_name + ": first order",
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR[
            baseline_index],
        "s--",
        color=baseline_colors[baseline_index],
        label=baseline_name + r": $+Q_{11}$",
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        SECOND_ORDER_STAGE1_INTERIOR_ERROR[
            baseline_index, 0],
        "^-.",
        color=baseline_colors[baseline_index],
        label=baseline_name + ": +spatial",
    )
    axes[1, 0].loglog(
        RADIUS_OVER_SMALL_MASS,
        COMPLETE_SECOND_ORDER_INTERIOR_ERROR[
            baseline_index],
        "D-",
        color=baseline_colors[baseline_index],
        label=baseline_name + ": complete",
    )

    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE[
            baseline_index],
        "o-",
        color=baseline_colors[baseline_index],
        label=baseline_name + r": temporal-even $G^{TT}$",
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE[
            baseline_index],
        "s--",
        color=baseline_colors[baseline_index],
        label=baseline_name + r": magnetic $G^{Ti}$",
    )

axes[0, 0].set_title(
    "Surface residual after each triangular stage")
axes[0, 0].set_ylabel(
    r"norm divided by candidate $\epsilon^2$ norm")
axes[0, 1].set_title(
    "Final surface residual by metric component")
axes[0, 1].set_ylabel(
    "norm divided by candidate component norm")
axes[1, 0].set_title(
    "Held-out random interior metric")
axes[1, 0].set_ylabel(
    "ten-component relative RMS")
axes[1, 1].set_title(
    "Residual in each newly fitted component block")
axes[1, 1].set_ylabel(
    "post-fit / pre-fit component norm")
for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=6)
fig_second_order_complete.suptitle(
    "Complete frozen-first-order q=8 second-order hierarchy")
fig_second_order_complete.tight_layout()
fig_second_order_complete.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_complete_validation.pdf")
fig_second_order_complete.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_complete_validation.png",
    dpi=180,
)
plt.show()


fig_second_order_new_groups, axes = plt.subplots(
    2, 2, figsize=(13.5, 9.0))
new_group_names = (
    "time acceleration",
    "time-gradient rate",
    "time harmonic quadrupole",
    "physical magnetic STF",
)
for axis, group_name in zip(
        axes.flat, new_group_names):
    group_slice = FULL_GROUP_SLICES[group_name]
    for baseline_index, baseline_name in enumerate(
            SECOND_ORDER_BASELINE_NAMES):
        values = np.linalg.norm(
            COMPLETE_SECOND_ORDER_PARAMETERS[
                baseline_index, :, group_slice],
            axis=1,
        )
        axis.loglog(
            RADIUS_OVER_SMALL_MASS,
            np.maximum(values, 1.0e-300),
            "o-" if baseline_index == 0 else "s--",
            color=baseline_colors[baseline_index],
            label=baseline_name,
        )
    axis.set_title(group_name)
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.set_ylabel("coefficient norm")
    axis.legend(fontsize=7)
fig_second_order_new_groups.suptitle(
    "New temporal-even and magnetic second-order coefficients")
fig_second_order_new_groups.tight_layout()
fig_second_order_new_groups.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_complete_parameters.pdf")
fig_second_order_new_groups.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_complete_parameters.png",
    dpi=180,
)
plt.show()


print("\nComplete second-order held-out interior summary")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    best_index = int(np.argmin(
        COMPLETE_SECOND_ORDER_INTERIOR_ERROR[
            baseline_index]))
    first_order_best = np.min(
        SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR[
            baseline_index])
    print(f"\n{baseline_name}")
    print(
        "  best complete interior error = "
        f"{COMPLETE_SECOND_ORDER_INTERIOR_ERROR[baseline_index, best_index]:.6e} "
        f"at R={RADIUS_TEST_VALUES[best_index]:.3f} "
        f"(R/M_B={RADIUS_OVER_SMALL_MASS[best_index]:.3f})"
    )
    print(
        "  best first-order interior error = "
        f"{first_order_best:.6e}"
    )
    print(
        "  final candidate-normalized surface residual = "
        f"{COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION[baseline_index, best_index, -1]:.6e}"
    )
"""
        ),
        markdown(
            r"""
### How to read the complete fit

The manufactured test establishes that the triangular algebra itself recovers
all 43 coefficients and closes the metric when the data are generated by the
same model.  It does **not** establish that the numerical binary metric is
described accurately by the truncation on this slice.

The decisive physical diagnostic remains the random interior prediction.  A
small fitted surface residual accompanied by a worse interior residual is
evidence that the second-order columns are absorbing omitted higher-order
structure or numerical error at the fitting sphere.  For that reason no
coefficient is promoted to a boundary prescription solely because its own
block closes.

That failure mode is what occurs on this slice.  The temporal \(G^{TT}\)
block closes to \(4\times10^{-4}\)--\(4\times10^{-3}\) of its pre-fit norm,
but its predicted \(G^{Ti}\) contribution is not supported by the numerical
residual.  The magnetic fit removes essentially none of that mismatch: its
post-fit/pre-fit \(G^{Ti}\) norm is unity to the displayed precision.  Thus
the final \(G^{TT}\) and \(G^{ij}\) closures are small while the \(G^{Ti}\)
residual is tens to hundreds of times its candidate norm.

The held-out interior check reaches the same conclusion.  The complete
hierarchy bottoms out near \(9.7\%\)--\(10.1\%\), whereas the corresponding
first-order models reach \(3.57\%\)--\(3.61\%\).  Since the manufactured
43-parameter test recovers coefficients to \(8\times10^{-12}\) and the two
new matrices have condition numbers \(1.53\) and \(1.00\), this is not a
rank or numerical-conditioning failure.  It is a cross-component closure
failure of this strictly triangular second-order estimator on the q=8
\(t=950M\) data.
"""
        ),
        markdown(
            r"""
## 13. Which second-order blocks carry recoverable information?

The failed complete triangular solve does not imply that every second-order
response is unusable.  We now repeat the more discriminating q=4 analysis.
There are three different questions:

1. **Response sensitivity:** can a physical coefficient group make a
   measurable change to the metric?
2. **Fit utility:** when that group is fitted by itself, does it reduce the
   surface candidate residual and improve the held-out interior metric?
3. **Cross-block consistency:** do coefficients recovered from one STF block
   correctly predict another block that was not fitted?

A small in-block residual answers none of the other questions.  Throughout
this section the first-order parameters remain frozen and the known
quadratic source \(Q_{11}\) is subtracted before fitting.
"""
        ),
        markdown(
            r"""
### 13.1 Irreducible STF subsystems

We begin with the same three targeted systems used for q=4.

The polar-\(J=2\) system combines
\[
  G^{TT}_{\ell=2},\qquad G^{Ti}_{\ell=3},\qquad
  G^{ij}_{\ell=4}
\]
to recover the STF spatial-strain rate, physical electric STF tensor, and
time-coordinate quadrupole: fifteen coefficients in total.

The spatial system fits the three accelerations and fifteen spatial
quadrupole coefficients using \(G^{ij}_{\ell=1}\), then predicts the held-out
\(G^{ij}_{\ell=3,5}\) blocks.

Finally, the fully symmetric \(J=3\) part of the spatial quadrupole has seven
coefficients that can be recovered independently from each of
\(G^{TT}_{\ell=3}\), \(G^{Ti}_{\ell=4}\), and \(G^{ij}_{\ell=5}\).  Agreement
between these independent estimators is a particularly clean closure test.
"""
        ),
        code(
            r"""
from itertools import permutations

POLAR_J2_BLOCKS = (
    ("GTT", 2), ("GTI", 3), ("GIJ", 4))
SPATIAL_ESTIMATOR_BLOCKS = (("GIJ", 1),)
SPATIAL_CLOSURE_BLOCKS = (
    ("GIJ", 3), ("GIJ", 5))
SPATIAL_J3_BLOCKS = (
    ("GTT", 3), ("GTI", 4), ("GIJ", 5))


def response_field_on_sphere(radius, builder):
    matrices = np.asarray([
        builder(direction) for direction in DIRECTIONS
    ])
    return matrix_to_metric_tuple(matrices)


def polar_j2_response_columns(radius):
    columns = []
    for basis in STF_BASIS:
        columns.append(response_field_on_sphere(
            radius,
            lambda direction, basis=basis:
            strain_rate_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            ),
        ))
    for basis in STF_BASIS:
        columns.append(response_field_on_sphere(
            radius,
            lambda direction, basis=basis:
            electric_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            ),
        ))
    for basis in STF_BASIS:
        columns.append(response_field_on_sphere(
            radius,
            lambda direction, basis=basis:
            time_quadrupole_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            ),
        ))
    return tuple(columns)


def projected_design(
        columns, blocks, parameter_indices=None):
    if parameter_indices is None:
        parameter_indices = np.arange(len(columns))
    return np.concatenate([
        np.column_stack([
            metric_stf_block(*columns[index], block)
            for index in parameter_indices
        ])
        for block in blocks
    ], axis=0)


def projected_data(metric, blocks):
    return np.concatenate([
        metric_stf_block(*metric, block)
        for block in blocks
    ])


def projected_residual_fraction(
        columns, estimate, metric, block,
        parameter_indices=None):
    matrix = projected_design(
        columns, (block,), parameter_indices)
    data = projected_data(metric, (block,))
    return (
        np.linalg.norm(data - matrix @ estimate)
        / max(np.linalg.norm(data), 1.0e-300)
    )


def orthonormal_stf_rank3_basis():
    symmetric_basis = []
    for i in range(3):
        for j in range(i, 3):
            for k in range(j, 3):
                tensor = np.zeros((3, 3, 3))
                index_permutations = sorted(
                    set(permutations((i, j, k))))
                for indices in index_permutations:
                    tensor[indices] = (
                        1.0 / np.sqrt(
                            len(index_permutations)))
                symmetric_basis.append(tensor)
    symmetric_basis = np.asarray(symmetric_basis)
    trace_constraints = np.column_stack([
        np.einsum(
            "aiik->ak", symmetric_basis)[:, component]
        for component in range(3)
    ]).T
    _, singular_values, right_vectors = np.linalg.svd(
        trace_constraints)
    rank = np.count_nonzero(
        singular_values
        > 1.0e-12 * singular_values[0])
    basis = np.einsum(
        "ab,bijk->aijk",
        right_vectors[rank:],
        symmetric_basis,
    )
    assert basis.shape == (7, 3, 3, 3)
    assert (
        np.max(np.abs(
            np.einsum("aiik->ak", basis)))
        < 1.0e-12
    )
    return basis


SPATIAL_J3_BASIS = orthonormal_stf_rank3_basis()


def spatial_j3_response_columns(radius):
    return tuple(
        response_field_on_sphere(
            radius,
            lambda direction, basis=basis:
            spatial_quadrupole_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            ),
        )
        for basis in SPATIAL_J3_BASIS
    )
"""
        ),
        code(
            r"""
POLAR_J2_ESTIMATES = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    15,
))
POLAR_J2_BLOCK_RESIDUALS = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    len(POLAR_J2_BLOCKS),
))
POLAR_J2_RANK = np.empty(
    number_of_radii, dtype=int)
POLAR_J2_CONDITION = np.empty(number_of_radii)

SPATIAL_ORTHOGONAL_ESTIMATES = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    18,
))
SPATIAL_ORTHOGONAL_IN_SAMPLE = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
))
SPATIAL_ORTHOGONAL_CLOSURE = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    len(SPATIAL_CLOSURE_BLOCKS),
))
SPATIAL_ORTHOGONAL_RANK = np.empty(
    number_of_radii, dtype=int)
SPATIAL_ORTHOGONAL_CONDITION = np.empty(
    number_of_radii)

SPATIAL_J3_ESTIMATES = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    len(SPATIAL_J3_BLOCKS),
    7,
))
SPATIAL_J3_CLOSURE = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    len(SPATIAL_J3_BLOCKS),
    len(SPATIAL_J3_BLOCKS),
))
SPATIAL_J3_RANK = np.empty((
    number_of_radii,
    len(SPATIAL_J3_BLOCKS),
), dtype=int)
SPATIAL_J3_CONDITION = np.empty((
    number_of_radii,
    len(SPATIAL_J3_BLOCKS),
))

spatial_parameter_indices = np.r_[
    np.arange(0, 3), np.arange(14, 29)]

for radius_index, radius in enumerate(
        RADIUS_TEST_VALUES):
    polar_columns = polar_j2_response_columns(radius)
    polar_matrix = projected_design(
        polar_columns, POLAR_J2_BLOCKS)
    _, polar_rank, polar_condition = (
        scaled_least_squares(
            polar_matrix,
            np.zeros(polar_matrix.shape[0]),
        )
    )
    POLAR_J2_RANK[radius_index] = polar_rank
    POLAR_J2_CONDITION[radius_index] = (
        polar_condition)

    complete_columns = (
        SURFACE_FULL_SECOND_ORDER_COLUMNS[
            radius_index])
    spatial_matrix = projected_design(
        complete_columns,
        SPATIAL_ESTIMATOR_BLOCKS,
        spatial_parameter_indices,
    )
    _, spatial_rank, spatial_condition = (
        scaled_least_squares(
            spatial_matrix,
            np.zeros(spatial_matrix.shape[0]),
        )
    )
    SPATIAL_ORTHOGONAL_RANK[
        radius_index] = spatial_rank
    SPATIAL_ORTHOGONAL_CONDITION[
        radius_index] = spatial_condition

    j3_columns = spatial_j3_response_columns(radius)
    for estimator_index, estimator_block in enumerate(
            SPATIAL_J3_BLOCKS):
        j3_matrix = projected_design(
            j3_columns, (estimator_block,))
        _, j3_rank, j3_condition = (
            scaled_least_squares(
                j3_matrix,
                np.zeros(j3_matrix.shape[0]),
            )
        )
        SPATIAL_J3_RANK[
            radius_index, estimator_index] = j3_rank
        SPATIAL_J3_CONDITION[
            radius_index, estimator_index] = (
                j3_condition)

    for baseline_index in range(
            number_of_second_order_baselines):
        candidate = SECOND_ORDER_CANDIDATE_RESIDUALS[
            baseline_index][radius_index]

        polar_data = projected_data(
            candidate, POLAR_J2_BLOCKS)
        polar_estimate, _, _ = scaled_least_squares(
            polar_matrix, polar_data)
        POLAR_J2_ESTIMATES[
            baseline_index, radius_index] = (
                polar_estimate)
        for block_index, block in enumerate(
                POLAR_J2_BLOCKS):
            POLAR_J2_BLOCK_RESIDUALS[
                baseline_index,
                radius_index,
                block_index,
            ] = projected_residual_fraction(
                polar_columns,
                polar_estimate,
                candidate,
                block,
            )

        spatial_data = projected_data(
            candidate, SPATIAL_ESTIMATOR_BLOCKS)
        spatial_estimate, _, _ = scaled_least_squares(
            spatial_matrix, spatial_data)
        SPATIAL_ORTHOGONAL_ESTIMATES[
            baseline_index, radius_index] = (
                spatial_estimate)
        SPATIAL_ORTHOGONAL_IN_SAMPLE[
            baseline_index, radius_index] = (
                np.linalg.norm(
                    spatial_data
                    - spatial_matrix @ spatial_estimate)
                / max(
                    np.linalg.norm(spatial_data),
                    1.0e-300,
                )
            )
        for block_index, block in enumerate(
                SPATIAL_CLOSURE_BLOCKS):
            closure_matrix = projected_design(
                complete_columns,
                (block,),
                spatial_parameter_indices,
            )
            closure_data = projected_data(
                candidate, (block,))
            SPATIAL_ORTHOGONAL_CLOSURE[
                baseline_index,
                radius_index,
                block_index,
            ] = (
                np.linalg.norm(
                    closure_data
                    - closure_matrix @ spatial_estimate)
                / max(
                    np.linalg.norm(closure_data),
                    1.0e-300,
                )
            )

        for estimator_index, estimator_block in enumerate(
                SPATIAL_J3_BLOCKS):
            j3_matrix = projected_design(
                j3_columns, (estimator_block,))
            j3_data = projected_data(
                candidate, (estimator_block,))
            j3_estimate, _, _ = scaled_least_squares(
                j3_matrix, j3_data)
            SPATIAL_J3_ESTIMATES[
                baseline_index,
                radius_index,
                estimator_index,
            ] = j3_estimate
            for prediction_index, prediction_block in enumerate(
                    SPATIAL_J3_BLOCKS):
                prediction_matrix = projected_design(
                    j3_columns, (prediction_block,))
                prediction_data = projected_data(
                    candidate, (prediction_block,))
                SPATIAL_J3_CLOSURE[
                    baseline_index,
                    radius_index,
                    estimator_index,
                    prediction_index,
                ] = (
                    np.linalg.norm(
                        prediction_data
                        - prediction_matrix @ j3_estimate)
                    / max(
                        np.linalg.norm(prediction_data),
                        1.0e-300,
                    )
                )

assert np.all(POLAR_J2_RANK == 15)
assert np.all(SPATIAL_ORTHOGONAL_RANK == 18)
assert np.all(SPATIAL_J3_RANK == 7)
print(
    "targeted condition ranges:",
    f"polar J2 {POLAR_J2_CONDITION.min():.3e}"
    f"--{POLAR_J2_CONDITION.max():.3e};",
    "spatial J1 "
    f"{SPATIAL_ORTHOGONAL_CONDITION.min():.3e}"
    f"--{SPATIAL_ORTHOGONAL_CONDITION.max():.3e};",
    "spatial J3 "
    f"{SPATIAL_J3_CONDITION.min():.3e}"
    f"--{SPATIAL_J3_CONDITION.max():.3e}",
)
"""
        ),
        code(
            r"""
fig_targeted_blocks, axes = plt.subplots(
    2, 2, figsize=(14.5, 10.0))
block_markers = ("o", "s", "^")

for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    for block_index, block in enumerate(
            POLAR_J2_BLOCKS):
        axes[0, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            POLAR_J2_BLOCK_RESIDUALS[
                baseline_index, :, block_index],
            marker=block_markers[block_index],
            linestyle=(
                "-" if baseline_index == 0 else "--"),
            color=f"C{block_index}",
            label=baseline_name + "; " + block_name(block),
        )

    axes[0, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        SPATIAL_ORTHOGONAL_IN_SAMPLE[
            baseline_index],
        "o-" if baseline_index == 0 else "o--",
        color=f"C{baseline_index}",
        label=baseline_name + r"; fit $G^{ij}_{\ell=1}$",
    )
    for block_index, block in enumerate(
            SPATIAL_CLOSURE_BLOCKS):
        axes[0, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SPATIAL_ORTHOGONAL_CLOSURE[
                baseline_index, :, block_index],
            marker=block_markers[block_index + 1],
            linestyle=(
                "-" if baseline_index == 0 else "--"),
            color=f"C{block_index + 2}",
            label=baseline_name + "; predict "
            + block_name(block),
        )

    polar_group_slices = (
        slice(0, 5), slice(5, 10), slice(10, 15))
    polar_group_labels = (
        "STF strain rate", "electric STF",
        "time quadrupole")
    for group_index, (
            group_slice, group_label) in enumerate(zip(
                polar_group_slices,
                polar_group_labels)):
        axes[1, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            np.linalg.norm(
                POLAR_J2_ESTIMATES[
                    baseline_index, :, group_slice],
                axis=1,
            ),
            marker=block_markers[group_index],
            linestyle=(
                "-" if baseline_index == 0 else "--"),
            color=f"C{group_index}",
            label=baseline_name + "; " + group_label,
        )

    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.linalg.norm(
            SPATIAL_ORTHOGONAL_ESTIMATES[
                baseline_index, :, :3],
            axis=1,
        ),
        "o-" if baseline_index == 0 else "o--",
        color="C0",
        label=baseline_name + "; acceleration",
    )
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        np.linalg.norm(
            SPATIAL_ORTHOGONAL_ESTIMATES[
                baseline_index, :, 3:],
            axis=1,
        ),
        "s-" if baseline_index == 0 else "s--",
        color="C1",
        label=baseline_name + "; spatial quadrupole",
    )

axes[0, 0].set_title(
    r"Polar-$J=2$ in-block residual")
axes[0, 0].set_ylabel(
    "fractional projected residual")
axes[0, 1].set_title(
    r"Spatial $G^{ij}_{\ell=1}$ fit and held-out closure")
axes[0, 1].set_ylabel(
    "fractional projected residual")
axes[1, 0].set_title(
    r"Polar-$J=2$ coefficient-group norms")
axes[1, 0].set_ylabel("coefficient norm")
axes[1, 1].set_title(
    "Spatial coefficient-group norms")
axes[1, 1].set_ylabel("coefficient norm")
for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=6)
fig_targeted_blocks.suptitle(
    "q=8 targeted second-order STF subsystems")
fig_targeted_blocks.tight_layout()
fig_targeted_blocks.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_targeted_blocks.pdf")
fig_targeted_blocks.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_targeted_blocks.png",
    dpi=180,
)
plt.show()


fig_j3_closure, axes = plt.subplots(
    1, 2, figsize=(15.0, 5.5), sharey=True)
for baseline_index, (
        baseline_name, axis) in enumerate(zip(
            SECOND_ORDER_BASELINE_NAMES, axes)):
    for estimator_index, estimator_block in enumerate(
            SPATIAL_J3_BLOCKS):
        for prediction_index, prediction_block in enumerate(
                SPATIAL_J3_BLOCKS):
            if estimator_index == prediction_index:
                continue
            axis.loglog(
                RADIUS_OVER_SMALL_MASS,
                SPATIAL_J3_CLOSURE[
                    baseline_index,
                    :,
                    estimator_index,
                    prediction_index,
                ],
                marker="o",
                label=(
                    "fit " + block_name(estimator_block)
                    + ", predict "
                    + block_name(prediction_block)
                ),
            )
    axis.set_title(baseline_name)
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.set_ylabel(
        "held-out fractional projected residual")
    axis.legend(fontsize=6)
fig_j3_closure.suptitle(
    r"Cross-component closure of spatial $J=3$ quadrupole")
fig_j3_closure.tight_layout()
fig_j3_closure.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_spatial_j3_closure.pdf")
fig_j3_closure.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_spatial_j3_closure.png",
    dpi=180,
)
plt.show()
"""
        ),
        markdown(
            r"""
### 13.2 Independent group fits and direct sensitivity

The STF tests expose consistency, but they do not directly measure how much
each physical group changes the complete inverse metric.  For every radius and
first-order baseline we therefore perform a separate global least-squares fit
for each of the eight second-order groups.

For each group we record:

- the reduction of the full ten-component surface candidate norm;
- the change in held-out interior error relative to the
  first-order-plus-\(Q_{11}\) model;
- the norm of the group's fitted interior contribution relative to the
  numerical non-Schwarzschild perturbation.

The groups are not combined in this test, so a useful group cannot borrow
coefficients from another group to hide a closure failure.
"""
        ),
        code(
            r"""
SECOND_ORDER_SENSITIVITY_GROUP_NAMES = tuple(
    FULL_GROUP_SLICES.keys())
SECOND_ORDER_SENSITIVITY_GROUP_SLICES = tuple(
    FULL_GROUP_SLICES.values())
number_of_sensitivity_groups = len(
    SECOND_ORDER_SENSITIVITY_GROUP_NAMES)

GROUP_ONLY_PARAMETERS = np.zeros((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_sensitivity_groups,
    number_of_full_second_order_parameters,
))
GROUP_ONLY_SURFACE_RESIDUAL_FRACTION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_sensitivity_groups,
))
GROUP_ONLY_SURFACE_REDUCTION = np.empty_like(
    GROUP_ONLY_SURFACE_RESIDUAL_FRACTION)
GROUP_ONLY_INTERIOR_ERROR = np.empty_like(
    GROUP_ONLY_SURFACE_RESIDUAL_FRACTION)
GROUP_ONLY_INTERIOR_ERROR_RATIO = np.empty_like(
    GROUP_ONLY_SURFACE_RESIDUAL_FRACTION)
GROUP_ONLY_INTERIOR_CONTRIBUTION = np.empty_like(
    GROUP_ONLY_SURFACE_RESIDUAL_FRACTION)
GROUP_ONLY_RANK = np.empty((
    number_of_radii,
    number_of_sensitivity_groups,
), dtype=int)
GROUP_ONLY_CONDITION = np.empty((
    number_of_radii,
    number_of_sensitivity_groups,
))

for radius_index in range(number_of_radii):
    columns = SURFACE_FULL_SECOND_ORDER_COLUMNS[
        radius_index]
    for group_index, group_slice in enumerate(
            SECOND_ORDER_SENSITIVITY_GROUP_SLICES):
        group_indices = np.arange(
            group_slice.start, group_slice.stop)
        surface_matrix = response_design(
            tuple(columns[index] for index in group_indices),
            lambda metric:
            metric_values_by_point(metric).ravel(),
        )
        _, group_rank, group_condition = (
            scaled_least_squares(
                surface_matrix,
                np.zeros(surface_matrix.shape[0]),
            )
        )
        GROUP_ONLY_RANK[
            radius_index, group_index] = group_rank
        GROUP_ONLY_CONDITION[
            radius_index, group_index] = group_condition

        interior_group_design = (
            INTERIOR_FULL_SECOND_ORDER_DESIGN[
                :, :, group_indices]
        )
        for baseline_index in range(
                number_of_second_order_baselines):
            candidate_values = metric_values_by_point(
                SECOND_ORDER_CANDIDATE_RESIDUALS[
                    baseline_index][radius_index]
            )
            estimate, _, _ = scaled_least_squares(
                surface_matrix,
                candidate_values.ravel(),
            )
            GROUP_ONLY_PARAMETERS[
                baseline_index,
                radius_index,
                group_index,
                group_indices,
            ] = estimate
            surface_remainder = (
                candidate_values.ravel()
                - surface_matrix @ estimate
            )
            surface_fraction = (
                np.linalg.norm(surface_remainder)
                / max(
                    np.linalg.norm(candidate_values),
                    1.0e-300,
                )
            )
            GROUP_ONLY_SURFACE_RESIDUAL_FRACTION[
                baseline_index,
                radius_index,
                group_index,
            ] = surface_fraction
            GROUP_ONLY_SURFACE_REDUCTION[
                baseline_index,
                radius_index,
                group_index,
            ] = max(0.0, 1.0 - surface_fraction)

            first_order_parameters = (
                SECOND_ORDER_FIRST_ORDER_BASELINES[
                    baseline_index, radius_index])
            first_order_values = np.einsum(
                "ncp,p->nc",
                INTERIOR_RESPONSE_DESIGN,
                first_order_parameters,
            )
            quadratic_values = metric_values_by_point(
                known_quadratic_first_order_source(
                    INTERIOR_LOCAL_POINTS,
                    first_order_parameters[1:4],
                    symmetric_matrix_from_six(
                        first_order_parameters[7:13]),
                    SMALL_BLACK_HOLE_MASS,
                )
            )
            quadratic_model = (
                INTERIOR_BACKGROUND_VALUES
                + first_order_values
                + quadratic_values
            )
            group_contribution = np.einsum(
                "ncp,p->nc",
                interior_group_design,
                estimate,
            )
            group_model = (
                quadratic_model + group_contribution)
            group_error = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - group_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )
            GROUP_ONLY_INTERIOR_ERROR[
                baseline_index,
                radius_index,
                group_index,
            ] = group_error
            GROUP_ONLY_INTERIOR_ERROR_RATIO[
                baseline_index,
                radius_index,
                group_index,
            ] = (
                group_error
                / COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR[
                    baseline_index, radius_index]
            )
            GROUP_ONLY_INTERIOR_CONTRIBUTION[
                baseline_index,
                radius_index,
                group_index,
            ] = (
                np.linalg.norm(group_contribution)
                / np.linalg.norm(
                    INTERIOR_PERTURBATION_VALUES)
            )
"""
        ),
        code(
            r"""
fig_group_only, axes = plt.subplots(
    2, 2, figsize=(15.0, 10.0))
group_colors = [
    f"C{index}" for index in range(
        number_of_sensitivity_groups)]

for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    for group_index, group_name in enumerate(
            SECOND_ORDER_SENSITIVITY_GROUP_NAMES):
        linestyle = (
            "-" if baseline_index == 0 else "--")
        label = baseline_name + "; " + group_name
        axes[0, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            np.maximum(
                GROUP_ONLY_SURFACE_REDUCTION[
                    baseline_index, :, group_index],
                1.0e-12,
            ),
            linestyle=linestyle,
            marker="o",
            color=group_colors[group_index],
            label=label,
        )
        axes[0, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            GROUP_ONLY_INTERIOR_ERROR_RATIO[
                baseline_index, :, group_index],
            linestyle=linestyle,
            marker="o",
            color=group_colors[group_index],
            label=label,
        )
        axes[1, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            np.maximum(
                GROUP_ONLY_INTERIOR_CONTRIBUTION[
                    baseline_index, :, group_index],
                1.0e-15,
            ),
            linestyle=linestyle,
            marker="o",
            color=group_colors[group_index],
            label=label,
        )

for group_index, group_name in enumerate(
        SECOND_ORDER_SENSITIVITY_GROUP_NAMES):
    axes[1, 1].loglog(
        RADIUS_OVER_SMALL_MASS,
        GROUP_ONLY_CONDITION[:, group_index],
        "o-",
        color=group_colors[group_index],
        label=group_name,
    )

axes[0, 0].set_title(
    "Surface residual reduction from group-only fit")
axes[0, 0].set_ylabel(
    r"$1-\|r_{\rm after}\|/\|r_{\rm candidate}\|$")
axes[0, 1].set_title(
    r"Held-out interior error relative to $G_{(1)}+Q_{11}$")
axes[0, 1].set_ylabel(
    r"$E_{\rm group}/E_{G_{(1)}+Q_{11}}$")
axes[0, 1].axhline(
    1.0, color="black", linestyle=":", linewidth=1.0)
axes[1, 0].set_title(
    "Fitted group contribution inside worldtube")
axes[1, 0].set_ylabel(
    r"$\|\Delta G_g\|/\|G_{\rm num}-G_{(0)}\|$")
axes[1, 1].set_title(
    "Normalized group-only condition number")
axes[1, 1].set_ylabel("condition number")
for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=5, ncol=2)
fig_group_only.suptitle(
    "Independent second-order group sensitivity")
fig_group_only.tight_layout()
fig_group_only.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_group_sensitivity.pdf")
fig_group_only.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_group_sensitivity.png",
    dpi=180,
)
plt.show()


print("\nBest independent group fits")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    print(f"\n{baseline_name}")
    print(
        "group                         "
        "best surface reduction  "
        "best interior ratio (radius)"
    )
    for group_index, group_name in enumerate(
            SECOND_ORDER_SENSITIVITY_GROUP_NAMES):
        surface_best = np.max(
            GROUP_ONLY_SURFACE_REDUCTION[
                baseline_index, :, group_index])
        interior_best_index = int(np.argmin(
            GROUP_ONLY_INTERIOR_ERROR_RATIO[
                baseline_index, :, group_index]))
        print(
            f"{group_name:29s} "
            f"{surface_best: .4e}          "
            f"{GROUP_ONLY_INTERIOR_ERROR_RATIO[baseline_index, interior_best_index, group_index]: .4e} "
            f"(R={RADIUS_TEST_VALUES[interior_best_index]:.2f})"
        )
"""
        ),
        markdown(
            r"""
### 13.3 Small fixed group combinations

The individual fits identify spatial-strain rate and time acceleration as the
only groups with a several-percent interior effect.  We now fit a short,
predeclared list of combinations.  This is not a search over the held-out
points: every combination is chosen from the physical group structure before
the interior errors are evaluated.

The comparison includes the two groups separately and together, two modest
extensions, the complete polar-\(J=2\) family, and a globally balanced fit of
all 43 columns.  Unlike the triangular solve, each combination minimizes the
full ten-component surface residual simultaneously.
"""
        ),
        code(
            r"""
SECOND_ORDER_SUBSET_DEFINITIONS = (
    ("strain rate", (1,)),
    ("time acceleration", (4,)),
    ("strain rate + time acceleration", (1, 4)),
    (
        "strain rate + time acceleration + electric",
        (1, 2, 4),
    ),
    (
        "strain rate + time acceleration + time quadrupole",
        (1, 4, 6),
    ),
    (
        "polar J2 groups",
        (1, 2, 6),
    ),
    (
        "all 43, globally balanced",
        tuple(range(number_of_sensitivity_groups)),
    ),
)
SECOND_ORDER_SUBSET_NAMES = tuple(
    definition[0]
    for definition in SECOND_ORDER_SUBSET_DEFINITIONS)
SECOND_ORDER_SUBSET_PARAMETER_INDICES = tuple(
    np.concatenate([
        np.arange(
            SECOND_ORDER_SENSITIVITY_GROUP_SLICES[
                group_index].start,
            SECOND_ORDER_SENSITIVITY_GROUP_SLICES[
                group_index].stop,
        )
        for group_index in group_indices
    ])
    for _, group_indices
    in SECOND_ORDER_SUBSET_DEFINITIONS
)
number_of_second_order_subsets = len(
    SECOND_ORDER_SUBSET_NAMES)

SECOND_ORDER_SUBSET_PARAMETERS = np.zeros((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_second_order_subsets,
    number_of_full_second_order_parameters,
))
SECOND_ORDER_SUBSET_SURFACE_FRACTION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_second_order_subsets,
))
SECOND_ORDER_SUBSET_INTERIOR_ERROR = np.empty_like(
    SECOND_ORDER_SUBSET_SURFACE_FRACTION)
SECOND_ORDER_SUBSET_INTERIOR_RATIO = np.empty_like(
    SECOND_ORDER_SUBSET_SURFACE_FRACTION)
SECOND_ORDER_SUBSET_RANK = np.empty((
    number_of_radii,
    number_of_second_order_subsets,
), dtype=int)
SECOND_ORDER_SUBSET_CONDITION = np.empty((
    number_of_radii,
    number_of_second_order_subsets,
))

for radius_index in range(number_of_radii):
    columns = SURFACE_FULL_SECOND_ORDER_COLUMNS[
        radius_index]
    for subset_index, parameter_indices in enumerate(
            SECOND_ORDER_SUBSET_PARAMETER_INDICES):
        surface_matrix = response_design(
            tuple(
                columns[index]
                for index in parameter_indices),
            lambda metric:
            metric_values_by_point(metric).ravel(),
        )
        for baseline_index in range(
                number_of_second_order_baselines):
            candidate_values = metric_values_by_point(
                SECOND_ORDER_CANDIDATE_RESIDUALS[
                    baseline_index][radius_index]
            ).ravel()
            estimate, rank, condition = (
                scaled_least_squares(
                    surface_matrix, candidate_values))
            SECOND_ORDER_SUBSET_PARAMETERS[
                baseline_index,
                radius_index,
                subset_index,
                parameter_indices,
            ] = estimate
            SECOND_ORDER_SUBSET_RANK[
                radius_index, subset_index] = rank
            SECOND_ORDER_SUBSET_CONDITION[
                radius_index, subset_index] = condition
            SECOND_ORDER_SUBSET_SURFACE_FRACTION[
                baseline_index,
                radius_index,
                subset_index,
            ] = (
                np.linalg.norm(
                    candidate_values
                    - surface_matrix @ estimate)
                / max(
                    np.linalg.norm(candidate_values),
                    1.0e-300,
                )
            )

            first_order_parameters = (
                SECOND_ORDER_FIRST_ORDER_BASELINES[
                    baseline_index, radius_index])
            quadratic_model = (
                INTERIOR_BACKGROUND_VALUES
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_RESPONSE_DESIGN,
                    first_order_parameters,
                )
                + metric_values_by_point(
                    known_quadratic_first_order_source(
                        INTERIOR_LOCAL_POINTS,
                        first_order_parameters[1:4],
                        symmetric_matrix_from_six(
                            first_order_parameters[7:13]),
                        SMALL_BLACK_HOLE_MASS,
                    )
                )
            )
            subset_prediction = np.einsum(
                "ncp,p->nc",
                INTERIOR_FULL_SECOND_ORDER_DESIGN[
                    :, :, parameter_indices],
                estimate,
            )
            subset_model = (
                quadratic_model + subset_prediction)
            subset_error = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - subset_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )
            SECOND_ORDER_SUBSET_INTERIOR_ERROR[
                baseline_index,
                radius_index,
                subset_index,
            ] = subset_error
            SECOND_ORDER_SUBSET_INTERIOR_RATIO[
                baseline_index,
                radius_index,
                subset_index,
            ] = (
                subset_error
                / COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR[
                    baseline_index, radius_index]
            )
"""
        ),
        code(
            r"""
fig_subset_comparison, axes = plt.subplots(
    1, 2, figsize=(15.0, 6.0))
subset_colors = [
    f"C{index}" for index in range(
        number_of_second_order_subsets)]
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    for subset_index, subset_name in enumerate(
            SECOND_ORDER_SUBSET_NAMES):
        linestyle = (
            "-" if baseline_index == 0 else "--")
        label = baseline_name + "; " + subset_name
        axes[0].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_SUBSET_SURFACE_FRACTION[
                baseline_index, :, subset_index],
            marker="o",
            linestyle=linestyle,
            color=subset_colors[subset_index],
            label=label,
        )
        axes[1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_SUBSET_INTERIOR_RATIO[
                baseline_index, :, subset_index],
            marker="o",
            linestyle=linestyle,
            color=subset_colors[subset_index],
            label=label,
        )
axes[0].set_title(
    "Full surface candidate residual")
axes[0].set_ylabel(
    r"$\|r_{\rm after}\|/\|r_{\rm candidate}\|$")
axes[1].set_title(
    r"Held-out interior error relative to $G_{(1)}+Q_{11}$")
axes[1].set_ylabel(
    r"$E_{\rm subset}/E_{G_{(1)}+Q_{11}}$")
axes[1].axhline(
    1.0, color="black", linestyle=":",
    linewidth=1.0)
for axis in axes:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=5, ncol=2)
fig_subset_comparison.suptitle(
    "Small globally balanced second-order subsets")
fig_subset_comparison.tight_layout()
fig_subset_comparison.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_subset_comparison.pdf")
fig_subset_comparison.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_subset_comparison.png",
    dpi=180,
)
plt.show()


print("\nFixed second-order subset comparison")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    print(f"\n{baseline_name}")
    for subset_index, subset_name in enumerate(
            SECOND_ORDER_SUBSET_NAMES):
        best_index = int(np.argmin(
            SECOND_ORDER_SUBSET_INTERIOR_RATIO[
                baseline_index, :, subset_index]))
        print(
            f"{subset_name:52s} "
            f"best ratio="
            f"{SECOND_ORDER_SUBSET_INTERIOR_RATIO[baseline_index, best_index, subset_index]:.6e} "
            f"at R={RADIUS_TEST_VALUES[best_index]:.2f}; "
            f"surface="
            f"{SECOND_ORDER_SUBSET_SURFACE_FRACTION[baseline_index, best_index, subset_index]:.6e}"
        )
"""
        ),
        markdown(
            r"""
### 13.4 Direct ablation of the complete triangular fit

The independent fits ask what each group can do alone.  The complementary
ablation test keeps the complete triangular coefficients fixed and removes one
group without refitting the others.  The contribution norm measures
sensitivity; the error ratio
\[
  E_{-g}/E_{\rm full}
\]
measures utility.  Values above one mean the group improves the interior
prediction, values near one mean it is irrelevant, and values below one mean
that removing the group makes the prediction better.
"""
        ),
        code(
            r"""
COMPLETE_ABLATION_CONTRIBUTION = np.empty((
    number_of_second_order_baselines,
    number_of_radii,
    number_of_sensitivity_groups,
))
COMPLETE_ABLATION_ERROR_RATIO = np.empty_like(
    COMPLETE_ABLATION_CONTRIBUTION)

for baseline_index in range(
        number_of_second_order_baselines):
    for radius_index in range(number_of_radii):
        first_order_parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index])
        first_order_values = np.einsum(
            "ncp,p->nc",
            INTERIOR_RESPONSE_DESIGN,
            first_order_parameters,
        )
        quadratic_values = metric_values_by_point(
            known_quadratic_first_order_source(
                INTERIOR_LOCAL_POINTS,
                first_order_parameters[1:4],
                symmetric_matrix_from_six(
                    first_order_parameters[7:13]),
                SMALL_BLACK_HOLE_MASS,
            )
        )
        complete_second_order_values = np.einsum(
            "ncp,p->nc",
            INTERIOR_FULL_SECOND_ORDER_DESIGN,
            COMPLETE_SECOND_ORDER_PARAMETERS[
                baseline_index, radius_index],
        )
        full_model = (
            INTERIOR_BACKGROUND_VALUES
            + first_order_values
            + quadratic_values
            + complete_second_order_values
        )
        full_error = np.linalg.norm(
            INTERIOR_NUMERICAL_VALUES - full_model)

        for group_index, group_slice in enumerate(
                SECOND_ORDER_SENSITIVITY_GROUP_SLICES):
            group_indices = np.arange(
                group_slice.start, group_slice.stop)
            group_contribution = np.einsum(
                "ncp,p->nc",
                INTERIOR_FULL_SECOND_ORDER_DESIGN[
                    :, :, group_indices],
                COMPLETE_SECOND_ORDER_PARAMETERS[
                    baseline_index,
                    radius_index,
                    group_indices,
                ],
            )
            ablated_model = (
                full_model - group_contribution)
            COMPLETE_ABLATION_CONTRIBUTION[
                baseline_index,
                radius_index,
                group_index,
            ] = (
                np.linalg.norm(group_contribution)
                / np.linalg.norm(
                    INTERIOR_PERTURBATION_VALUES)
            )
            COMPLETE_ABLATION_ERROR_RATIO[
                baseline_index,
                radius_index,
                group_index,
            ] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - ablated_model)
                / full_error
            )
"""
        ),
        code(
            r"""
fig_ablation_scan, axes = plt.subplots(
    2, 2, figsize=(15.0, 10.0))
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    contribution_axis = axes[0, baseline_index]
    error_axis = axes[1, baseline_index]
    for group_index, group_name in enumerate(
            SECOND_ORDER_SENSITIVITY_GROUP_NAMES):
        contribution_axis.loglog(
            RADIUS_OVER_SMALL_MASS,
            np.maximum(
                COMPLETE_ABLATION_CONTRIBUTION[
                    baseline_index, :, group_index],
                1.0e-15,
            ),
            "o-",
            color=group_colors[group_index],
            label=group_name,
        )
        error_axis.loglog(
            RADIUS_OVER_SMALL_MASS,
            COMPLETE_ABLATION_ERROR_RATIO[
                baseline_index, :, group_index],
            "o-",
            color=group_colors[group_index],
            label=group_name,
        )
    contribution_axis.set_title(
        baseline_name + ": metric sensitivity")
    contribution_axis.set_ylabel(
        r"$\|\Delta G_g\|/\|G_{\rm num}-G_{(0)}\|$")
    error_axis.set_title(
        baseline_name + ": removal test")
    error_axis.set_ylabel(r"$E_{-g}/E_{\rm full}$")
    error_axis.axhline(
        1.0, color="black", linestyle=":",
        linewidth=1.0)
    for axis in (contribution_axis, error_axis):
        axis.set_xlabel(r"fitting radius $R/M_B$")
        axis.legend(fontsize=6)
fig_ablation_scan.suptitle(
    "Direct physical-group ablation of complete triangular fit")
fig_ablation_scan.tight_layout()
fig_ablation_scan.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_group_ablation.pdf")
fig_ablation_scan.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_group_ablation.png",
    dpi=180,
)
plt.show()
"""
        ),
        markdown(
            r"""
### Reading the sensitivity analysis

The three diagnostics must be read together.  A large contribution norm means
the metric is sensitive to the fitted coefficient, not that the coefficient
is accurate.  A strong independent surface reduction with an interior error
ratio above one indicates overfitting at the worldtube.  Conversely, a group
with a modest surface effect but a stable interior improvement may be a better
candidate for the matching subsystem.

The irreducible cross-block plots are the strictest recovery test: a block is
credible only if its coefficients are reasonably stable with radius and
predict multipoles that were not included in the fit.

On the q=8 slice, the sensitivity is concentrated but the result is more
promising than the triangular solve suggested:

- fitted independently, spatial-strain rate lowers the best interior error by
  about \(7.8\%\)--\(7.9\%\), and time acceleration by about
  \(7.1\%\)--\(7.2\%\);
- fitting those two groups jointly lowers it by \(10.6\%\) for the metric-only
  baseline and by as much as \(15.5\%\) for the derivative-assisted baseline;
- electric STF, time quadrupole, acceleration, time-gradient rate, and spatial
  quadrupole each change the best group-only interior error by less than one
  percent; the magnetic STF response is numerically irrelevant at this
  accuracy;
- adding electric STF or time quadrupole to the strain-rate/time-acceleration
  pair does not improve its best small-radius result.

The polar-\(J=2\) subsystem is the strongest targeted block.  At the smallest
radius its held-in \(G^{TT}_{\ell=2}\) block closes algebraically, while the
\(G^{Ti}_{\ell=3}\) and \(G^{ij}_{\ell=4}\) residuals are approximately
\(10\%\) and \(5\%\), respectively.  Its held-out interior prediction reaches
\(3.13\times10^{-2}\), compared with \(3.57\times10^{-2}\) for the best
first-order baseline.

By contrast, the spatial \(G^{ij}_{\ell=1}\) subsystem is not validated.  It
fits its square estimator block to roundoff, but the held-out
\(\ell=3,5\) residuals remain near or above unity.  The three independent
spatial-\(J=3\) estimators also fail to give uniformly good cross-component
closure.  These spatial quadrupole coefficients should therefore not yet be
used as individually recovered quantities.

A globally balanced 43-column fit reaches the lowest interior error,
\(3.01\times10^{-2}\), without the catastrophic \(G^{Ti}\) amplification of
the triangular solve.  This shows that second-order response directions do
contain useful aggregate information, but it does not establish that all 43
coefficients are separately resolved.  The conservative provisional
subsystem is therefore the spatial-strain rate plus time acceleration; the
polar-\(J=2\) system is the next block worth testing with radial derivatives.
"""
        ),
        markdown(
            r"""
## 14. Joint value-and-radial-derivative fits at second order

We now add
\[
  {\cal D}G^{ab}\equiv R n^i\partial_iG^{ab}
\]
to the two provisional second-order systems:

1. symmetric spatial-strain rate plus time acceleration;
2. the fifteen-parameter polar-\(J=2\) subsystem.

No new unknown is introduced.  The same coefficient must multiply a response
\(\mathcal R_A(R)\) in the metric equation and
\(R\partial_R\mathcal R_A(R)\) in the derivative equation.  The numerical
derivative comes directly from the evolved \(\Phi_{iab}\).  Only the analytic
background, first-order responses, quadratic source, and second-order
responses are differentiated with the converged five-point radial stencil.

For derivative weight \(w_D\), the joint system is
\[
\begin{pmatrix}
 A\\ w_D A_D
\end{pmatrix}p
\simeq
\begin{pmatrix}
 d\\ w_D d_D
\end{pmatrix}.
\]
We scan the same weights used at first order.  Metric values at the fixed
random interior points remain completely held out.
"""
        ),
        code(
            r"""
SECOND_ORDER_DERIVATIVE_WEIGHTS = (
    JOINT_DERIVATIVE_WEIGHTS.copy())
SECOND_ORDER_DERIVATIVE_WEIGHT_LABELS = (
    JOINT_WEIGHT_LABELS)


def analytic_full_second_order_response_derivatives(
        radius, relative_step=1.0e-4):
    step = relative_step * radius
    samples = [
        full_response_columns_at_points(
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
        for parameter in range(
            number_of_full_second_order_parameters)
    )


def analytic_polar_j2_response_derivatives(
        radius, relative_step=1.0e-4):
    step = relative_step * radius
    samples = [
        polar_j2_response_columns(
            radius + offset * step)
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
        for parameter in range(15)
    )


def polar_j2_response_columns_at_points(points):
    points = np.asarray(points)
    radii = np.linalg.norm(points, axis=1)
    directions = points / radii[:, None]
    builders = []
    for basis in STF_BASIS:
        builders.append(
            lambda radius, direction, basis=basis:
            strain_rate_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            )
        )
    for basis in STF_BASIS:
        builders.append(
            lambda radius, direction, basis=basis:
            electric_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            )
        )
    for basis in STF_BASIS:
        builders.append(
            lambda radius, direction, basis=basis:
            time_quadrupole_response_at_direction(
                SMALL_BLACK_HOLE_MASS,
                radius,
                direction,
                basis,
            )
        )
    return tuple(
        matrix_to_metric_tuple(np.asarray([
            builder(radius, direction)
            for radius, direction in zip(
                radii, directions)
        ]))
        for builder in builders
    )


SECOND_ORDER_DERIVATIVE_CANDIDATES = [
    [None for _ in RADIUS_TEST_VALUES]
    for _ in SECOND_ORDER_BASELINE_NAMES
]
SURFACE_FULL_SECOND_ORDER_DERIVATIVES = []
SURFACE_POLAR_J2_DERIVATIVES = []
SECOND_ORDER_DERIVATIVE_STENCIL_ERROR = np.empty(
    number_of_radii)

for radius_index, radius in enumerate(
        RADIUS_TEST_VALUES):
    numerical_derivative = metric_tuple_from_matrix(
        NUMERICAL_RADIAL_DERIVATIVE[radius_index])
    background_derivative = (
        analytic_background_radial_derivative(radius))
    first_order_derivatives = (
        analytic_response_radial_derivatives(radius))
    full_second_order_derivatives = (
        analytic_full_second_order_response_derivatives(
            radius))
    polar_derivatives = (
        analytic_polar_j2_response_derivatives(radius))
    SURFACE_FULL_SECOND_ORDER_DERIVATIVES.append(
        tuple(
            tuple(radius * component for component in column)
            for column in full_second_order_derivatives
        )
    )
    SURFACE_POLAR_J2_DERIVATIVES.append(
        tuple(
            tuple(radius * component for component in column)
            for column in polar_derivatives
        )
    )

    half_step = (
        analytic_full_second_order_response_derivatives(
            radius, relative_step=5.0e-5))
    full_values = np.concatenate([
        metric_values_by_point(column).ravel()
        for column in full_second_order_derivatives
    ])
    half_values = np.concatenate([
        metric_values_by_point(column).ravel()
        for column in half_step
    ])
    SECOND_ORDER_DERIVATIVE_STENCIL_ERROR[
        radius_index] = (
            np.linalg.norm(full_values - half_values)
            / np.linalg.norm(half_values)
        )

    for baseline_index in range(
            number_of_second_order_baselines):
        first_order_parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index])
        first_order_derivative_prediction = (
            metric_linear_combination(
                first_order_derivatives,
                first_order_parameters,
            )
        )
        quadratic_derivative = (
            five_point_radial_derivative(
                lambda points:
                known_quadratic_first_order_source(
                    points,
                    first_order_parameters[1:4],
                    symmetric_matrix_from_six(
                        first_order_parameters[7:13]),
                    SMALL_BLACK_HOLE_MASS,
                ),
                radius,
            )
        )
        SECOND_ORDER_DERIVATIVE_CANDIDATES[
            baseline_index][radius_index] = tuple(
                radius * (
                    numerical
                    - background
                    - first_order
                    - quadratic
                )
                for numerical, background, first_order, quadratic
                in zip(
                    numerical_derivative,
                    background_derivative,
                    first_order_derivative_prediction,
                    quadratic_derivative,
                )
        )

print(
    "second-order response derivative stencil check:",
    f"max={SECOND_ORDER_DERIVATIVE_STENCIL_ERROR.max():.3e}",
)
assert (
    SECOND_ORDER_DERIVATIVE_STENCIL_ERROR.max()
    < 2.0e-9
)
"""
        ),
        code(
            r"""
PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES = np.r_[
    np.arange(
        FULL_GROUP_SLICES[
            "symmetric spatial-strain rate"].start,
        FULL_GROUP_SLICES[
            "symmetric spatial-strain rate"].stop,
    ),
    np.arange(
        FULL_GROUP_SLICES[
            "time acceleration"].start,
        FULL_GROUP_SLICES[
            "time acceleration"].stop,
    ),
]
number_of_second_order_weights = len(
    SECOND_ORDER_DERIVATIVE_WEIGHTS)

SECOND_ORDER_JOINT_SUBSET_PARAMETERS = np.zeros((
    number_of_second_order_baselines,
    number_of_second_order_weights,
    number_of_radii,
    number_of_full_second_order_parameters,
))
SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL = np.empty((
    number_of_second_order_baselines,
    number_of_second_order_weights,
    number_of_radii,
))
SECOND_ORDER_JOINT_SUBSET_DERIVATIVE_RESIDUAL = (
    np.empty_like(
        SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL)
)
SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR = np.empty_like(
    SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL)
SECOND_ORDER_JOINT_SUBSET_CONDITION = np.empty_like(
    SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL)

SECOND_ORDER_JOINT_POLAR_PARAMETERS = np.empty((
    number_of_second_order_baselines,
    number_of_second_order_weights,
    number_of_radii,
    15,
))
SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL = np.empty((
    number_of_second_order_baselines,
    number_of_second_order_weights,
    number_of_radii,
))
SECOND_ORDER_JOINT_POLAR_DERIVATIVE_RESIDUAL = (
    np.empty_like(
        SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL)
)
SECOND_ORDER_JOINT_POLAR_INTERIOR_ERROR = np.empty_like(
    SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL)
SECOND_ORDER_JOINT_POLAR_CONDITION = np.empty_like(
    SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL)

INTERIOR_POLAR_J2_COLUMNS = (
    polar_j2_response_columns_at_points(
        INTERIOR_LOCAL_POINTS))
INTERIOR_POLAR_J2_DESIGN = np.stack([
    metric_values_by_point(column)
    for column in INTERIOR_POLAR_J2_COLUMNS
], axis=-1)

for radius_index, radius in enumerate(
        RADIUS_TEST_VALUES):
    full_columns = SURFACE_FULL_SECOND_ORDER_COLUMNS[
        radius_index]
    full_derivatives = (
        SURFACE_FULL_SECOND_ORDER_DERIVATIVES[
            radius_index])
    subset_value_matrix = response_design(
        tuple(
            full_columns[index]
            for index in
            PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES),
        lambda metric:
        metric_values_by_point(metric).ravel(),
    )
    subset_derivative_matrix = response_design(
        tuple(
            full_derivatives[index]
            for index in
            PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES),
        lambda metric:
        metric_values_by_point(metric).ravel(),
    )

    polar_columns = polar_j2_response_columns(radius)
    polar_derivatives = (
        SURFACE_POLAR_J2_DERIVATIVES[radius_index])
    polar_value_matrix = projected_design(
        polar_columns, POLAR_J2_BLOCKS)
    polar_derivative_matrix = projected_design(
        polar_derivatives, POLAR_J2_BLOCKS)

    for baseline_index in range(
            number_of_second_order_baselines):
        candidate = SECOND_ORDER_CANDIDATE_RESIDUALS[
            baseline_index][radius_index]
        derivative_candidate = (
            SECOND_ORDER_DERIVATIVE_CANDIDATES[
                baseline_index][radius_index])
        subset_value_data = (
            metric_values_by_point(candidate).ravel())
        subset_derivative_data = (
            metric_values_by_point(
                derivative_candidate).ravel())
        polar_value_data = projected_data(
            candidate, POLAR_J2_BLOCKS)
        polar_derivative_data = projected_data(
            derivative_candidate, POLAR_J2_BLOCKS)

        first_order_parameters = (
            SECOND_ORDER_FIRST_ORDER_BASELINES[
                baseline_index, radius_index])
        quadratic_model = (
            INTERIOR_BACKGROUND_VALUES
            + np.einsum(
                "ncp,p->nc",
                INTERIOR_RESPONSE_DESIGN,
                first_order_parameters,
            )
            + metric_values_by_point(
                known_quadratic_first_order_source(
                    INTERIOR_LOCAL_POINTS,
                    first_order_parameters[1:4],
                    symmetric_matrix_from_six(
                        first_order_parameters[7:13]),
                    SMALL_BLACK_HOLE_MASS,
                )
            )
        )

        for weight_index, derivative_weight in enumerate(
                SECOND_ORDER_DERIVATIVE_WEIGHTS):
            subset_matrix = np.concatenate((
                subset_value_matrix,
                derivative_weight
                * subset_derivative_matrix,
            ), axis=0)
            subset_data = np.concatenate((
                subset_value_data,
                derivative_weight
                * subset_derivative_data,
            ))
            subset_estimate, subset_rank, subset_condition = (
                scaled_least_squares(
                    subset_matrix, subset_data))
            assert subset_rank == len(
                PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES)
            SECOND_ORDER_JOINT_SUBSET_PARAMETERS[
                baseline_index,
                weight_index,
                radius_index,
                PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES,
            ] = subset_estimate
            SECOND_ORDER_JOINT_SUBSET_CONDITION[
                baseline_index,
                weight_index,
                radius_index,
            ] = subset_condition
            SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    subset_value_data
                    - subset_value_matrix
                    @ subset_estimate)
                / np.linalg.norm(subset_value_data)
            )
            SECOND_ORDER_JOINT_SUBSET_DERIVATIVE_RESIDUAL[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    subset_derivative_data
                    - subset_derivative_matrix
                    @ subset_estimate)
                / np.linalg.norm(
                    subset_derivative_data)
            )
            subset_model = (
                quadratic_model
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_FULL_SECOND_ORDER_DESIGN[
                        :,
                        :,
                        PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES,
                    ],
                    subset_estimate,
                )
            )
            SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - subset_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )

            polar_matrix = np.concatenate((
                polar_value_matrix,
                derivative_weight
                * polar_derivative_matrix,
            ), axis=0)
            polar_data = np.concatenate((
                polar_value_data,
                derivative_weight
                * polar_derivative_data,
            ))
            polar_estimate, polar_rank, polar_condition = (
                scaled_least_squares(
                    polar_matrix, polar_data))
            assert polar_rank == 15
            SECOND_ORDER_JOINT_POLAR_PARAMETERS[
                baseline_index,
                weight_index,
                radius_index,
            ] = polar_estimate
            SECOND_ORDER_JOINT_POLAR_CONDITION[
                baseline_index,
                weight_index,
                radius_index,
            ] = polar_condition
            SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    polar_value_data
                    - polar_value_matrix @ polar_estimate)
                / np.linalg.norm(polar_value_data)
            )
            SECOND_ORDER_JOINT_POLAR_DERIVATIVE_RESIDUAL[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    polar_derivative_data
                    - polar_derivative_matrix @ polar_estimate)
                / np.linalg.norm(
                    polar_derivative_data)
            )
            polar_model = (
                quadratic_model
                + np.einsum(
                    "ncp,p->nc",
                    INTERIOR_POLAR_J2_DESIGN,
                    polar_estimate,
                )
            )
            SECOND_ORDER_JOINT_POLAR_INTERIOR_ERROR[
                baseline_index,
                weight_index,
                radius_index,
            ] = (
                np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES
                    - polar_model)
                / np.linalg.norm(
                    INTERIOR_NUMERICAL_VALUES)
            )

# Weight zero must reproduce the corresponding value-only fits.
subset_definition_index = (
    SECOND_ORDER_SUBSET_NAMES.index(
        "strain rate + time acceleration"))
np.testing.assert_allclose(
    SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR[:, 0],
    SECOND_ORDER_SUBSET_INTERIOR_ERROR[
        :, :, subset_definition_index],
    rtol=3.0e-12,
    atol=3.0e-13,
)
np.testing.assert_allclose(
    SECOND_ORDER_JOINT_POLAR_PARAMETERS[:, 0],
    POLAR_J2_ESTIMATES,
    rtol=3.0e-11,
    atol=3.0e-12,
)
"""
        ),
        code(
            r"""
fig_second_order_joint_derivative, axes = plt.subplots(
    2, 2, figsize=(15.0, 10.0))
weight_colors = [
    f"C{index}" for index in range(
        number_of_second_order_weights)]

for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    linestyle = (
        "-" if baseline_index == 0 else "--")
    for weight_index, weight_label in enumerate(
            SECOND_ORDER_DERIVATIVE_WEIGHT_LABELS):
        label = baseline_name + "; " + weight_label
        axes[0, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR[
                baseline_index, weight_index],
            marker="o",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            label=label,
        )
        axes[0, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_POLAR_INTERIOR_ERROR[
                baseline_index, weight_index],
            marker="o",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            label=label,
        )
        axes[1, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL[
                baseline_index, weight_index],
            marker="o",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            label=label + ": value",
        )
        axes[1, 0].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_SUBSET_DERIVATIVE_RESIDUAL[
                baseline_index, weight_index],
            marker="s",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            alpha=0.65,
            label=label + ": derivative",
        )
        axes[1, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL[
                baseline_index, weight_index],
            marker="o",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            label=label + ": value",
        )
        axes[1, 1].loglog(
            RADIUS_OVER_SMALL_MASS,
            SECOND_ORDER_JOINT_POLAR_DERIVATIVE_RESIDUAL[
                baseline_index, weight_index],
            marker="s",
            linestyle=linestyle,
            color=weight_colors[weight_index],
            alpha=0.65,
            label=label + ": derivative",
        )

axes[0, 0].set_title(
    "Strain rate + time acceleration: interior metric")
axes[0, 0].set_ylabel(
    "ten-component relative RMS")
axes[0, 1].set_title(
    r"Polar-$J=2$: interior metric")
axes[0, 1].set_ylabel(
    "ten-component relative RMS")
axes[1, 0].set_title(
    "Subset surface value/derivative residual")
axes[1, 0].set_ylabel("fractional residual")
axes[1, 1].set_title(
    r"Polar-$J=2$ value/derivative residual")
axes[1, 1].set_ylabel("fractional residual")
for axis in axes.flat:
    axis.set_xlabel(r"fitting radius $R/M_B$")
    axis.legend(fontsize=4.5, ncol=2)
fig_second_order_joint_derivative.suptitle(
    "Joint value-and-radial-derivative second-order fits")
fig_second_order_joint_derivative.tight_layout()
fig_second_order_joint_derivative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_joint_radial_derivative.pdf")
fig_second_order_joint_derivative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_second_order_joint_radial_derivative.png",
    dpi=180,
)
plt.show()


print("\nSecond-order joint radial-derivative summary")
for baseline_index, baseline_name in enumerate(
        SECOND_ORDER_BASELINE_NAMES):
    print(f"\n{baseline_name}")
    for family_name, errors in (
        (
            "strain rate + time acceleration",
            SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR[
                baseline_index],
        ),
        (
            "polar J2",
            SECOND_ORDER_JOINT_POLAR_INTERIOR_ERROR[
                baseline_index],
        ),
    ):
        best = np.unravel_index(
            np.argmin(errors), errors.shape)
        print(
            f"  {family_name:36s}: "
            f"interior={errors[best]:.6e}, "
            f"weight={SECOND_ORDER_DERIVATIVE_WEIGHTS[best[0]]:g}, "
            f"R={RADIUS_TEST_VALUES[best[1]]:.3f}"
        )
"""
        ),
        markdown(
            r"""
### Reading the second-order derivative fit

The derivative should be judged by whether it improves the held-out interior
metric over a finite range of weights and radii, not by whether it makes the
derivative residual alone small.  A very large weight can trade away the
metric-value fit and is not automatically preferable.

The weight-zero curves exactly reproduce the corresponding value-only
systems.  This provides an end-to-end check of the response ordering before
interpreting any derivative-assisted change.

The derivative data are useful, but in different ways for the two systems.
For spatial-strain rate plus time acceleration, \(w_D=1\) moves the preferred
radius from \(0.20M\) to \(0.24M\) and lowers the interior error from
\(3.27\times10^{-2}\) to \(3.19\times10^{-2}\) for the metric-only
first-order baseline.  With the derivative-assisted first-order baseline it
lowers the error from \(3.24\times10^{-2}\) to
\(3.17\times10^{-2}\).  This is a modest but genuine improvement at held-out
points; weights \(3\) and \(10\) begin to trade away too much metric-value
information.

The targeted polar-\(J=2\) system benefits more dramatically.  Its
value-only projected-block fit extrapolates poorly, with a best interior error
near \(4.1\times10^{-2}\).  Adding only \(w_D=0.1\) reduces this to
\(3.35\times10^{-2}\) and \(3.32\times10^{-2}\) for the two first-order
baselines.  Larger derivative weights improve the normalized matrix condition
number—from as large as \(62\) in the value-only system to below \(10\)—but
do not further improve the interior metric.  Better conditioning is therefore
not by itself the correct weight-selection criterion.

The analytic second-order derivative columns are numerically converged: halving
the radial stencil changes them by at most \(4.4\times10^{-12}\).  Together
with the earlier \(\Phi\)-versus-interpolation check, this makes it unlikely
that the improvement is a finite-difference artifact.

Thus the radial derivative should be retained.  The most conservative current
choice is the seven-parameter strain-rate/time-acceleration system with
\(w_D\) near one.  The polar-\(J=2\) derivative system is also informative,
but still does not outperform the globally balanced value-only polar-group
fit; its cross-block closure should be examined before promoting all fifteen
coefficients.
"""
        ),
        markdown(
            r"""
## 15. Exploratory multi-radius alternation (superseded below)

The preceding fits freeze the first-order coefficients before solving for the
second-order coefficients.  That is useful for diagnosis, but it leaves a
possible leakage mechanism: an imperfect first-order estimate can project
onto the second-order response space, and a fitted second-order contribution
can in turn change the first-order blocks.

This first experiment uses several radii in one solve.  It is retained as a
diagnostic of cross-order coupling, but it is *not* the matching prescription.
Section 16 replaces it with the radius-by-radius construction used throughout
the rest of the notebook.

Here we test an alternating, block-Gauss--Seidel iteration.  Write the
background-subtracted numerical metric and its dimensionless radial derivative
as
\[
\begin{aligned}
d &=
 A_1p_1+Q_{11}(p_1,p_1)+A_2p_2+\mathcal O(\epsilon^3),\\
d_D &=
 A_{1D}p_1+Q_{11,D}(p_1,p_1)+A_{2D}p_2
 +\mathcal O(\epsilon^3).
\end{aligned}
\]
At iteration \(k\), we first hold \(p_1^{(k)}\) fixed and solve
\[
 p_{2,\star}^{(k+1)}
 =
 F_2\!\left[
 d-A_1p_1^{(k)}-Q_{11}(p_1^{(k)},p_1^{(k)})
 \right].
\]
We then subtract both the known quadratic source and the new second-order
prediction before reapplying the *original* first-order estimators,
\[
 p_{1,\star}^{(k+1)}
 =
 F_1\!\left[
 d-Q_{11}(p_1^{(k)},p_1^{(k)})
       -A_2p_2^{(k+1)}
 \right].
\]
Thus the iteration does not introduce a free thirteen-component
\(\delta p_1\).  Every first-order update must still be inferred through the
same \(V1+C3\) blocks used and validated above.

Both block updates are relaxed,
\[
 p_j^{(k+1)}
 = (1-\lambda)p_j^{(k)}
   +\lambda p_{j,\star}^{(k+1)},\qquad j=1,2,
\]
and we compare \(\lambda=1/4,1/2,1\).  The first-order vector estimator uses
the radial derivative with \(w_D=3\), while the clock/strain estimator remains
metric-only.  The second-order solve contains only the seven provisional
coefficients (six symmetric spatial-strain-rate components and the time
acceleration) and uses values plus \(R\partial_R\) with \(w_D=1\).

Most importantly, one common \(p_1,p_2\) is fitted simultaneously over
\(R=0.20\)--\(0.40M\).  The different radii are extra equations, not separate
estimates.  Random interior points and all radii larger than \(0.40M\) remain
held out.
"""
        ),
        code(
            r"""
ITERATIVE_TRAINING_RADIUS_INDICES = np.flatnonzero(
    (RADIUS_TEST_VALUES >= 0.20)
    & (RADIUS_TEST_VALUES <= 0.40)
)
ITERATIVE_HELD_OUT_RADIUS_INDICES = np.flatnonzero(
    RADIUS_TEST_VALUES > 0.40)
ITERATIVE_DAMPING_VALUES = np.array([0.25, 0.5, 1.0])
ITERATIVE_MAX_ITERATIONS = 30
ITERATIVE_VECTOR_DERIVATIVE_WEIGHT = 3.0
ITERATIVE_SECOND_ORDER_DERIVATIVE_WEIGHT = 1.0


def add_metrics(*metrics):
    return tuple(
        sum(metric[component] for metric in metrics)
        for component in range(3)
    )


def zero_metric_like(metric):
    return tuple(np.zeros_like(component) for component in metric)


ITERATIVE_RAW_VALUES = []
ITERATIVE_RAW_DERIVATIVES = []
ITERATIVE_FIRST_ORDER_COLUMNS = []
ITERATIVE_FIRST_ORDER_DERIVATIVES = []

for radius_index, radius in enumerate(RADIUS_TEST_VALUES):
    numerical_metric = unpack_inverse_metric(
        FIRST_ORDER_METRIC_DATA[:, radius_index])
    background_metric = inverse_harmonic_schwarzschild(
        radius * DIRECTIONS, SMALL_BLACK_HOLE_MASS)
    ITERATIVE_RAW_VALUES.append(
        subtract_metric(numerical_metric, background_metric))

    numerical_derivative = metric_tuple_from_matrix(
        NUMERICAL_RADIAL_DERIVATIVE[radius_index])
    background_derivative = (
        analytic_background_radial_derivative(radius))
    ITERATIVE_RAW_DERIVATIVES.append(tuple(
        radius * (numerical - background)
        for numerical, background in zip(
            numerical_derivative, background_derivative)
    ))

    ITERATIVE_FIRST_ORDER_COLUMNS.append(
        first_order_response_columns(
            radius * DIRECTIONS, SMALL_BLACK_HOLE_MASS))
    ITERATIVE_FIRST_ORDER_DERIVATIVES.append(tuple(
        tuple(radius * component for component in column)
        for column in analytic_response_radial_derivatives(radius)
    ))


def quadratic_value_and_derivative(
        first_order_parameters, radius_index):
    radius = RADIUS_TEST_VALUES[radius_index]
    beta = first_order_parameters[1:4]
    spatial_strain = symmetric_matrix_from_six(
        first_order_parameters[7:13])
    value = known_quadratic_first_order_source(
        radius * DIRECTIONS,
        beta,
        spatial_strain,
        SMALL_BLACK_HOLE_MASS,
    )
    derivative = five_point_radial_derivative(
        lambda points:
        known_quadratic_first_order_source(
            points,
            beta,
            spatial_strain,
            SMALL_BLACK_HOLE_MASS,
        ),
        radius,
    )
    return value, tuple(
        radius * component for component in derivative)


def second_order_value_and_derivative(
        second_order_parameters, radius_index):
    return (
        metric_linear_combination(
            SURFACE_FULL_SECOND_ORDER_COLUMNS[radius_index],
            second_order_parameters,
        ),
        metric_linear_combination(
            SURFACE_FULL_SECOND_ORDER_DERIVATIVES[radius_index],
            second_order_parameters,
        ),
    )


def fit_iterative_first_order(
        raw_values, raw_derivatives,
        correction_values, correction_derivatives,
        radius_indices=ITERATIVE_TRAINING_RADIUS_INDICES):
    result = np.zeros(len(FIRST_ORDER_PARAMETER_NAMES))
    sectors = (
        (
            JOINT_VECTOR_FIT_BLOCKS,
            VECTOR_PARAMETER_INDICES,
            ITERATIVE_VECTOR_DERIVATIVE_WEIGHT,
        ),
        (
            JOINT_CLOCK_FIT_BLOCKS,
            CLOCK_STRAIN_PARAMETER_INDICES,
            0.0,
        ),
    )
    for blocks, parameter_indices, derivative_weight in sectors:
        matrix_rows = []
        data_rows = []
        for radius_index in radius_indices:
            corrected_value = subtract_metric(
                raw_values[radius_index],
                correction_values[radius_index],
            )
            matrix_rows.append(projected_design(
                ITERATIVE_FIRST_ORDER_COLUMNS[radius_index],
                blocks,
                parameter_indices,
            ))
            data_rows.append(projected_data(
                corrected_value, blocks))
            if derivative_weight > 0.0:
                corrected_derivative = subtract_metric(
                    raw_derivatives[radius_index],
                    correction_derivatives[radius_index],
                )
                matrix_rows.append(
                    derivative_weight * projected_design(
                        ITERATIVE_FIRST_ORDER_DERIVATIVES[
                            radius_index],
                        blocks,
                        parameter_indices,
                    )
                )
                data_rows.append(
                    derivative_weight * projected_data(
                        corrected_derivative, blocks)
                )
        estimate, rank, _ = scaled_least_squares(
            np.concatenate(matrix_rows, axis=0),
            np.concatenate(data_rows),
        )
        if rank != len(parameter_indices):
            raise RuntimeError(
                "Iterative first-order block lost rank: "
                f"{rank} != {len(parameter_indices)}")
        result[parameter_indices] = estimate
    return result


def fit_iterative_second_order(
        first_order_parameters,
        raw_values, raw_derivatives,
        radius_indices=ITERATIVE_TRAINING_RADIUS_INDICES):
    matrix_rows = []
    data_rows = []
    for radius_index in radius_indices:
        first_order_value = metric_linear_combination(
            ITERATIVE_FIRST_ORDER_COLUMNS[radius_index],
            first_order_parameters,
        )
        first_order_derivative = metric_linear_combination(
            ITERATIVE_FIRST_ORDER_DERIVATIVES[radius_index],
            first_order_parameters,
        )
        quadratic_value, quadratic_derivative = (
            quadratic_value_and_derivative(
                first_order_parameters, radius_index))
        candidate_value = subtract_metric(
            raw_values[radius_index],
            first_order_value,
            quadratic_value,
        )
        candidate_derivative = subtract_metric(
            raw_derivatives[radius_index],
            first_order_derivative,
            quadratic_derivative,
        )
        matrix_rows.append(response_design(
            tuple(
                SURFACE_FULL_SECOND_ORDER_COLUMNS[
                    radius_index][index]
                for index in
                PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES
            ),
            lambda metric:
            metric_values_by_point(metric).ravel(),
        ))
        data_rows.append(
            metric_values_by_point(
                candidate_value).ravel())
        matrix_rows.append(
            ITERATIVE_SECOND_ORDER_DERIVATIVE_WEIGHT
            * response_design(
                tuple(
                    SURFACE_FULL_SECOND_ORDER_DERIVATIVES[
                        radius_index][index]
                    for index in
                    PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES
                ),
                lambda metric:
                metric_values_by_point(metric).ravel(),
            )
        )
        data_rows.append(
            ITERATIVE_SECOND_ORDER_DERIVATIVE_WEIGHT
            * metric_values_by_point(
                candidate_derivative).ravel()
        )
    estimate, rank, _ = scaled_least_squares(
        np.concatenate(matrix_rows, axis=0),
        np.concatenate(data_rows),
    )
    if rank != len(
            PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES):
        raise RuntimeError(
            "Iterative second-order block lost rank: "
            f"{rank} != "
            f"{len(PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES)}")
    result = np.zeros(number_of_full_second_order_parameters)
    result[
        PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES] = estimate
    return result


def iterative_training_residual(
        first_order_parameters,
        second_order_parameters,
        raw_values, raw_derivatives,
        radius_indices=ITERATIVE_TRAINING_RADIUS_INDICES):
    residual_parts = []
    data_parts = []
    for radius_index in radius_indices:
        first_order_value = metric_linear_combination(
            ITERATIVE_FIRST_ORDER_COLUMNS[radius_index],
            first_order_parameters,
        )
        first_order_derivative = metric_linear_combination(
            ITERATIVE_FIRST_ORDER_DERIVATIVES[radius_index],
            first_order_parameters,
        )
        quadratic_value, quadratic_derivative = (
            quadratic_value_and_derivative(
                first_order_parameters, radius_index))
        second_order_value, second_order_derivative = (
            second_order_value_and_derivative(
                second_order_parameters, radius_index))
        residual_parts.extend((
            metric_values_by_point(subtract_metric(
                raw_values[radius_index],
                first_order_value,
                quadratic_value,
                second_order_value,
            )).ravel(),
            ITERATIVE_SECOND_ORDER_DERIVATIVE_WEIGHT
            * metric_values_by_point(subtract_metric(
                raw_derivatives[radius_index],
                first_order_derivative,
                quadratic_derivative,
                second_order_derivative,
            )).ravel(),
        ))
        data_parts.extend((
            metric_values_by_point(
                raw_values[radius_index]).ravel(),
            ITERATIVE_SECOND_ORDER_DERIVATIVE_WEIGHT
            * metric_values_by_point(
                raw_derivatives[radius_index]).ravel(),
        ))
    return (
        np.linalg.norm(np.concatenate(residual_parts))
        / np.linalg.norm(np.concatenate(data_parts))
    )


def iterative_interior_error(
        first_order_parameters, second_order_parameters):
    first_order_values = np.einsum(
        "ncp,p->nc",
        INTERIOR_RESPONSE_DESIGN,
        first_order_parameters,
    )
    quadratic_values = metric_values_by_point(
        known_quadratic_first_order_source(
            INTERIOR_LOCAL_POINTS,
            first_order_parameters[1:4],
            symmetric_matrix_from_six(
                first_order_parameters[7:13]),
            SMALL_BLACK_HOLE_MASS,
        )
    )
    second_order_values = np.einsum(
        "ncp,p->nc",
        INTERIOR_FULL_SECOND_ORDER_DESIGN,
        second_order_parameters,
    )
    model = (
        INTERIOR_BACKGROUND_VALUES
        + first_order_values
        + quadratic_values
        + second_order_values
    )
    return (
        np.linalg.norm(INTERIOR_NUMERICAL_VALUES - model)
        / np.linalg.norm(INTERIOR_NUMERICAL_VALUES)
    )


def run_alternating_iteration(
        raw_values, raw_derivatives, damping,
        number_of_iterations=ITERATIVE_MAX_ITERATIONS,
        interior_error_function=None,
        true_first_order=None, true_second_order=None,
        radius_indices=ITERATIVE_TRAINING_RADIUS_INDICES):
    zero_value_corrections = [
        zero_metric_like(metric) for metric in raw_values]
    zero_derivative_corrections = [
        zero_metric_like(metric) for metric in raw_derivatives]
    first_order = fit_iterative_first_order(
        raw_values,
        raw_derivatives,
        zero_value_corrections,
        zero_derivative_corrections,
        radius_indices,
    )
    second_order = np.zeros(
        number_of_full_second_order_parameters)

    history = {
        "first_order": [first_order.copy()],
        "second_order": [second_order.copy()],
        "training_residual": [
            iterative_training_residual(
                first_order, second_order,
                raw_values, raw_derivatives,
                radius_indices)],
        "first_order_update": [np.nan],
        "second_order_update": [np.nan],
        "interior_error": [
            np.nan if interior_error_function is None
            else interior_error_function(
                first_order, second_order)],
        "first_order_error": [
            np.nan if true_first_order is None
            else np.linalg.norm(
                first_order - true_first_order)
            / max(np.linalg.norm(true_first_order), 1.0e-300)],
        "second_order_error": [
            np.nan if true_second_order is None
            else np.linalg.norm(
                second_order - true_second_order)
            / max(np.linalg.norm(true_second_order), 1.0e-300)],
    }

    for _ in range(number_of_iterations):
        second_order_target = fit_iterative_second_order(
            first_order, raw_values, raw_derivatives,
            radius_indices)
        next_second_order = (
            (1.0 - damping) * second_order
            + damping * second_order_target
        )

        correction_values = []
        correction_derivatives = []
        for radius_index in range(number_of_radii):
            quadratic_value, quadratic_derivative = (
                quadratic_value_and_derivative(
                    first_order, radius_index))
            second_order_value, second_order_derivative = (
                second_order_value_and_derivative(
                    next_second_order, radius_index))
            correction_values.append(add_metrics(
                quadratic_value, second_order_value))
            correction_derivatives.append(add_metrics(
                quadratic_derivative, second_order_derivative))

        first_order_target = fit_iterative_first_order(
            raw_values,
            raw_derivatives,
            correction_values,
            correction_derivatives,
            radius_indices,
        )
        next_first_order = (
            (1.0 - damping) * first_order
            + damping * first_order_target
        )

        history["first_order_update"].append(
            np.linalg.norm(next_first_order - first_order)
            / max(np.linalg.norm(next_first_order), 1.0e-300))
        history["second_order_update"].append(
            np.linalg.norm(next_second_order - second_order)
            / max(np.linalg.norm(next_second_order), 1.0e-300))
        first_order = next_first_order
        second_order = next_second_order
        history["first_order"].append(first_order.copy())
        history["second_order"].append(second_order.copy())
        history["training_residual"].append(
            iterative_training_residual(
                first_order, second_order,
                raw_values, raw_derivatives,
                radius_indices))
        history["interior_error"].append(
            np.nan if interior_error_function is None
            else interior_error_function(
                first_order, second_order))
        history["first_order_error"].append(
            np.nan if true_first_order is None
            else np.linalg.norm(
                first_order - true_first_order)
            / max(np.linalg.norm(true_first_order), 1.0e-300))
        history["second_order_error"].append(
            np.nan if true_second_order is None
            else np.linalg.norm(
                second_order - true_second_order)
            / max(np.linalg.norm(true_second_order), 1.0e-300))

    return {
        name: np.asarray(values)
        for name, values in history.items()
    }
"""
        ),
        markdown(
            r"""
### Manufactured fixed-point test

Before applying the iteration to simulation data, we manufacture data from
exactly the same nonlinear truncated model.  The chosen coefficients have
realistic magnitudes taken from the current q=8 fits, but the test is otherwise
noise-free.  Successful recovery verifies more than the individual linear
solvers: it checks the signs, ordering, quadratic subtraction, radial
derivatives, and alternating update as one coupled calculation.
"""
        ),
        code(
            r"""
MANUFACTURED_ITERATIVE_FIRST_ORDER = np.mean(
    SECOND_ORDER_FIRST_ORDER_BASELINES[
        1, ITERATIVE_TRAINING_RADIUS_INDICES],
    axis=0,
)
weight_one_index = int(np.flatnonzero(np.isclose(
    SECOND_ORDER_DERIVATIVE_WEIGHTS, 1.0))[0])
MANUFACTURED_ITERATIVE_SECOND_ORDER = np.mean(
    SECOND_ORDER_JOINT_SUBSET_PARAMETERS[
        1,
        weight_one_index,
        ITERATIVE_TRAINING_RADIUS_INDICES,
    ],
    axis=0,
)

MANUFACTURED_ITERATIVE_RAW_VALUES = []
MANUFACTURED_ITERATIVE_RAW_DERIVATIVES = []
for radius_index in range(number_of_radii):
    first_order_value = metric_linear_combination(
        ITERATIVE_FIRST_ORDER_COLUMNS[radius_index],
        MANUFACTURED_ITERATIVE_FIRST_ORDER,
    )
    first_order_derivative = metric_linear_combination(
        ITERATIVE_FIRST_ORDER_DERIVATIVES[radius_index],
        MANUFACTURED_ITERATIVE_FIRST_ORDER,
    )
    quadratic_value, quadratic_derivative = (
        quadratic_value_and_derivative(
            MANUFACTURED_ITERATIVE_FIRST_ORDER,
            radius_index,
        )
    )
    second_order_value, second_order_derivative = (
        second_order_value_and_derivative(
            MANUFACTURED_ITERATIVE_SECOND_ORDER,
            radius_index,
        )
    )
    MANUFACTURED_ITERATIVE_RAW_VALUES.append(add_metrics(
        first_order_value,
        quadratic_value,
        second_order_value,
    ))
    MANUFACTURED_ITERATIVE_RAW_DERIVATIVES.append(add_metrics(
        first_order_derivative,
        quadratic_derivative,
        second_order_derivative,
    ))

# Verify the exact fixed-point identity separately from the convergence rate
# of an iteration started far from that fixed point.
MANUFACTURED_FIXED_SECOND_ORDER = fit_iterative_second_order(
    MANUFACTURED_ITERATIVE_FIRST_ORDER,
    MANUFACTURED_ITERATIVE_RAW_VALUES,
    MANUFACTURED_ITERATIVE_RAW_DERIVATIVES,
)
manufactured_fixed_correction_values = []
manufactured_fixed_correction_derivatives = []
for radius_index in range(number_of_radii):
    quadratic_value, quadratic_derivative = (
        quadratic_value_and_derivative(
            MANUFACTURED_ITERATIVE_FIRST_ORDER,
            radius_index,
        )
    )
    second_order_value, second_order_derivative = (
        second_order_value_and_derivative(
            MANUFACTURED_FIXED_SECOND_ORDER,
            radius_index,
        )
    )
    manufactured_fixed_correction_values.append(add_metrics(
        quadratic_value, second_order_value))
    manufactured_fixed_correction_derivatives.append(add_metrics(
        quadratic_derivative, second_order_derivative))
MANUFACTURED_FIXED_FIRST_ORDER = fit_iterative_first_order(
    MANUFACTURED_ITERATIVE_RAW_VALUES,
    MANUFACTURED_ITERATIVE_RAW_DERIVATIVES,
    manufactured_fixed_correction_values,
    manufactured_fixed_correction_derivatives,
)
MANUFACTURED_FIXED_FIRST_ORDER_ERROR = (
    np.linalg.norm(
        MANUFACTURED_FIXED_FIRST_ORDER
        - MANUFACTURED_ITERATIVE_FIRST_ORDER)
    / np.linalg.norm(MANUFACTURED_ITERATIVE_FIRST_ORDER)
)
MANUFACTURED_FIXED_SECOND_ORDER_ERROR = (
    np.linalg.norm(
        MANUFACTURED_FIXED_SECOND_ORDER
        - MANUFACTURED_ITERATIVE_SECOND_ORDER)
    / np.linalg.norm(MANUFACTURED_ITERATIVE_SECOND_ORDER)
)
print(
    "manufactured fixed-point identity:",
    f"p1 error={MANUFACTURED_FIXED_FIRST_ORDER_ERROR:.3e},",
    f"p2 error={MANUFACTURED_FIXED_SECOND_ORDER_ERROR:.3e}",
)
assert MANUFACTURED_FIXED_FIRST_ORDER_ERROR < 2.0e-11
assert MANUFACTURED_FIXED_SECOND_ORDER_ERROR < 2.0e-11

MANUFACTURED_ITERATION_RESULTS = [
    run_alternating_iteration(
        MANUFACTURED_ITERATIVE_RAW_VALUES,
        MANUFACTURED_ITERATIVE_RAW_DERIVATIVES,
        damping,
        true_first_order=MANUFACTURED_ITERATIVE_FIRST_ORDER,
        true_second_order=MANUFACTURED_ITERATIVE_SECOND_ORDER,
    )
    for damping in ITERATIVE_DAMPING_VALUES
]

for damping, results in zip(
        ITERATIVE_DAMPING_VALUES,
        MANUFACTURED_ITERATION_RESULTS):
    print(
        f"manufactured lambda={damping:g}: "
        f"p1 error={results['first_order_error'][-1]:.3e}, "
        f"p2 error={results['second_order_error'][-1]:.3e}, "
        f"residual={results['training_residual'][-1]:.3e}"
    )

manufactured_best_index = int(np.argmin([
    results["training_residual"][-1]
    for results in MANUFACTURED_ITERATION_RESULTS
]))
assert (
    MANUFACTURED_ITERATION_RESULTS[
        manufactured_best_index]["training_residual"][-1]
    < MANUFACTURED_ITERATION_RESULTS[
        manufactured_best_index]["training_residual"][0]
)
"""
        ),
        markdown(
            r"""
### Alternation on the q=8 slice

We now run precisely the same iteration on the numerical data.  Convergence
of the iteration is distinct from accuracy of the truncated model.  A stable
fixed point shows that the two selected blocks can be made mutually
consistent; only the held-out interior metric and held-out radii tell us
whether that fixed point predicts the simulation.
"""
        ),
        code(
            r"""
REAL_ITERATION_RESULTS = [
    run_alternating_iteration(
        ITERATIVE_RAW_VALUES,
        ITERATIVE_RAW_DERIVATIVES,
        damping,
        interior_error_function=iterative_interior_error,
    )
    for damping in ITERATIVE_DAMPING_VALUES
]

REAL_ITERATION_HELD_OUT_RESIDUAL = np.empty((
    len(ITERATIVE_DAMPING_VALUES),
    ITERATIVE_MAX_ITERATIONS + 1,
))
for damping_index, results in enumerate(
        REAL_ITERATION_RESULTS):
    for iteration in range(ITERATIVE_MAX_ITERATIONS + 1):
        REAL_ITERATION_HELD_OUT_RESIDUAL[
            damping_index, iteration] = (
                iterative_training_residual(
                    results["first_order"][iteration],
                    results["second_order"][iteration],
                    ITERATIVE_RAW_VALUES,
                    ITERATIVE_RAW_DERIVATIVES,
                    ITERATIVE_HELD_OUT_RADIUS_INDICES,
                )
            )

fig_iterative, axes = plt.subplots(
    2, 2, figsize=(13.5, 9.0))
iteration_numbers = np.arange(
    ITERATIVE_MAX_ITERATIONS + 1)
for damping, manufactured, real, held_out in zip(
        ITERATIVE_DAMPING_VALUES,
        MANUFACTURED_ITERATION_RESULTS,
        REAL_ITERATION_RESULTS,
        REAL_ITERATION_HELD_OUT_RESIDUAL,
):
    label = rf"$\lambda={damping:g}$"
    axes[0, 0].semilogy(
        iteration_numbers,
        manufactured["first_order_error"],
        "-",
        label=label + r": $p_1$",
    )
    axes[0, 0].semilogy(
        iteration_numbers,
        manufactured["second_order_error"],
        "--",
        label=label + r": $p_2$",
    )
    axes[0, 1].semilogy(
        iteration_numbers,
        real["training_residual"],
        "o-",
        markevery=3,
        label=label,
    )
    axes[1, 0].semilogy(
        iteration_numbers,
        real["interior_error"],
        "o-",
        markevery=3,
        label=label,
    )
    axes[1, 1].semilogy(
        iteration_numbers,
        held_out,
        "o-",
        markevery=3,
        label=label,
    )

axes[0, 0].set_title("Manufactured parameter recovery")
axes[0, 0].set_ylabel("relative parameter error")
axes[0, 1].set_title(
    r"Real-data training residual, $0.20\leq R/M\leq0.40$")
axes[0, 1].set_ylabel("joint value/derivative residual")
axes[1, 0].set_title("Held-out random interior inverse metric")
axes[1, 0].set_ylabel("ten-component relative RMS")
axes[1, 1].set_title(
    r"Held-out radii, $R/M>0.40$")
axes[1, 1].set_ylabel("joint value/derivative residual")
for axis in axes.flat:
    axis.set_xlabel("alternating iteration")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=7)
fig_iterative.suptitle(
    "Alternating first- and second-order matching")
fig_iterative.tight_layout()
fig_iterative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_alternating_first_second_order.pdf")
fig_iterative.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_alternating_first_second_order.png",
    dpi=180,
)
plt.show()


print("\nAlternating real-data summary")
for damping, results, held_out in zip(
        ITERATIVE_DAMPING_VALUES,
        REAL_ITERATION_RESULTS,
        REAL_ITERATION_HELD_OUT_RESIDUAL):
    best_interior_iteration = int(np.argmin(
        results["interior_error"]))
    print(
        f"lambda={damping:g}: "
        f"final training={results['training_residual'][-1]:.6e}, "
        f"final held-out radii={held_out[-1]:.6e}, "
        f"best interior={results['interior_error'][best_interior_iteration]:.6e} "
        f"at iteration {best_interior_iteration}, "
        f"final interior={results['interior_error'][-1]:.6e}"
    )
"""
        ),
        markdown(
            r"""
### Interpretation of the iteration

The manufactured test determines whether the alternating map is a valid
solver for the truncated equations.  The real-data curves answer a different
question: whether iterating toward that fixed point improves predictions that
were not used in either block solve.

If the training residual converges while the interior or large-radius
residual worsens, the iteration is successfully reallocating truncation error
between orders but is overfitting the selected blocks.  Conversely, a stable
minimum shared by several damping factors is evidence that the first- and
second-order estimates are mutually consistent.  This separation is why the
iteration number is treated as a validation axis rather than automatically
running to convergence.

The manufactured fixed-point identity is satisfied to
\(4\times10^{-14}\) in \(p_1\) and \(10^{-15}\) in \(p_2\), so the two
block solvers and both subtractions are algebraically consistent.  Starting
from the naive first-order estimate is nevertheless slow: even the undamped
iteration retains \(1.2\%\) relative error in \(p_1\) and \(19\%\) in \(p_2\)
after thirty sweeps.  The first- and second-order response spaces are
therefore strongly coupled rather than approximately block diagonal.

For the numerical slice, all three damping choices initially lower the
training residual.  Their minima occur after approximately the same amount of
effective iteration: iteration \(8\) for \(\lambda=1/4\), \(4\) for
\(\lambda=1/2\), and \(2\) for \(\lambda=1\), at a residual near \(0.146\).
Thereafter the residual rises.  The held-out-radii residual is even more
restrictive: it is minimized at iteration \(2\), \(1\), and \(0\),
respectively, and subsequently grows substantially.

The random interior metric behaves differently.  Its error decreases
monotonically from \(3.72\times10^{-2}\) to
\(2.96,2.89,2.80\times10^{-2}\) for
\(\lambda=1/4,1/2,1\) after thirty iterations.  The smallest value is better
than the approximately \(3.01\times10^{-2}\) globally balanced value-only
second-order fit found above, but it coincides with worse extrapolation to
larger worldtubes and neither parameter block has converged.

So the experiment supports the *possibility* of iterative decontamination,
but not unrestricted iteration to a fixed point.  At present it is best
viewed as a regularization path: an iteration count of order one to four gives
a reproducible compromise between the small-radius equations and nearby
held-out radii, while later sweeps continue improving this particular
interior point set by reallocating unresolved higher-order content into
\(p_1\) and \(p_2\).  Additional time slices or derivative observables are
needed to decide which validation target is physically predictive.
"""
        ),
        markdown(
            r"""
## 16. Radius-by-radius alternating matching

We now perform the matching independently on each worldtube.  For a fixed
radius \(R\),
\[
\begin{aligned}
d(R) &= A_1(R)p_1(R)
      +Q_{11}[p_1(R),p_1(R)]+A_2(R)p_2(R),\\
d_D(R) &= A_{1D}(R)p_1(R)
      +Q_{11,D}[p_1(R),p_1(R)]+A_{2D}(R)p_2(R).
\end{aligned}
\]
No information from any other radius enters either least-squares solve.  The
system is nevertheless overdetermined because the ten metric components at
all angular points and their radial derivatives provide many more equations
than the thirteen first-order and seven selected second-order coefficients.

At every radius from \(0.20M\) through \(0.80M\) we initialize
\(p_1^{(0)}(R)\) with the same
derivative-assisted \(V1+C3\) estimator, alternate the two undamped block
solves for thirty iterations, and retain the complete trajectory.  The
value-plus-derivative residual on that sphere is training information.  The
metric at the common random points inside the \(R=0.20M\) sphere is not used
by either solve and is the primary held-out diagnostic.
"""
        ),
        code(
            r"""
RADIUS_BY_RADIUS_DAMPING = 1.0
RADIUS_BY_RADIUS_MAX_RADIUS = 0.8
RADIUS_BY_RADIUS_RADIUS_INDICES = np.flatnonzero(
    RADIUS_TEST_VALUES <= RADIUS_BY_RADIUS_MAX_RADIUS)
RADIUS_BY_RADIUS_RADII = RADIUS_TEST_VALUES[
    RADIUS_BY_RADIUS_RADIUS_INDICES]
RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS = (
    RADIUS_OVER_SMALL_MASS[
        RADIUS_BY_RADIUS_RADIUS_INDICES])
RADIUS_BY_RADIUS_ITERATION_RESULTS = [
    run_alternating_iteration(
        ITERATIVE_RAW_VALUES,
        ITERATIVE_RAW_DERIVATIVES,
        RADIUS_BY_RADIUS_DAMPING,
        interior_error_function=iterative_interior_error,
        radius_indices=np.asarray([radius_index]),
    )
    for radius_index in RADIUS_BY_RADIUS_RADIUS_INDICES
]

RADIUS_BY_RADIUS_FIRST_ORDER = np.asarray([
    result["first_order"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])
RADIUS_BY_RADIUS_SECOND_ORDER = np.asarray([
    result["second_order"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])
RADIUS_BY_RADIUS_SURFACE_RESIDUAL = np.asarray([
    result["training_residual"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])
RADIUS_BY_RADIUS_INTERIOR_ERROR = np.asarray([
    result["interior_error"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])
RADIUS_BY_RADIUS_FIRST_ORDER_UPDATE = np.asarray([
    result["first_order_update"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])
RADIUS_BY_RADIUS_SECOND_ORDER_UPDATE = np.asarray([
    result["second_order_update"]
    for result in RADIUS_BY_RADIUS_ITERATION_RESULTS
])

# The known manufactured solution must be an exact fixed point at every
# individual radius, independently of how fast an iteration approaches it.
RADIUS_BY_RADIUS_MANUFACTURED_P1_ERROR = np.empty(
    len(RADIUS_BY_RADIUS_RADIUS_INDICES))
RADIUS_BY_RADIUS_MANUFACTURED_P2_ERROR = np.empty(
    len(RADIUS_BY_RADIUS_RADIUS_INDICES))
for local_radius_index, radius_index in enumerate(
        RADIUS_BY_RADIUS_RADIUS_INDICES):
    radius_indices = np.asarray([radius_index])
    fixed_second_order = fit_iterative_second_order(
        MANUFACTURED_ITERATIVE_FIRST_ORDER,
        MANUFACTURED_ITERATIVE_RAW_VALUES,
        MANUFACTURED_ITERATIVE_RAW_DERIVATIVES,
        radius_indices,
    )
    correction_values = [
        zero_metric_like(metric)
        for metric in MANUFACTURED_ITERATIVE_RAW_VALUES
    ]
    correction_derivatives = [
        zero_metric_like(metric)
        for metric in MANUFACTURED_ITERATIVE_RAW_DERIVATIVES
    ]
    quadratic_value, quadratic_derivative = (
        quadratic_value_and_derivative(
            MANUFACTURED_ITERATIVE_FIRST_ORDER,
            radius_index,
        )
    )
    second_order_value, second_order_derivative = (
        second_order_value_and_derivative(
            fixed_second_order, radius_index)
    )
    correction_values[radius_index] = add_metrics(
        quadratic_value, second_order_value)
    correction_derivatives[radius_index] = add_metrics(
        quadratic_derivative, second_order_derivative)
    fixed_first_order = fit_iterative_first_order(
        MANUFACTURED_ITERATIVE_RAW_VALUES,
        MANUFACTURED_ITERATIVE_RAW_DERIVATIVES,
        correction_values,
        correction_derivatives,
        radius_indices,
    )
    RADIUS_BY_RADIUS_MANUFACTURED_P1_ERROR[
        local_radius_index] = (
        np.linalg.norm(
            fixed_first_order
            - MANUFACTURED_ITERATIVE_FIRST_ORDER)
        / np.linalg.norm(MANUFACTURED_ITERATIVE_FIRST_ORDER)
    )
    RADIUS_BY_RADIUS_MANUFACTURED_P2_ERROR[
        local_radius_index] = (
        np.linalg.norm(
            fixed_second_order
            - MANUFACTURED_ITERATIVE_SECOND_ORDER)
        / np.linalg.norm(MANUFACTURED_ITERATIVE_SECOND_ORDER)
    )

print(
    "maximum single-radius manufactured fixed-point errors:",
    f"p1={RADIUS_BY_RADIUS_MANUFACTURED_P1_ERROR.max():.3e},",
    f"p2={RADIUS_BY_RADIUS_MANUFACTURED_P2_ERROR.max():.3e}",
)
assert RADIUS_BY_RADIUS_MANUFACTURED_P1_ERROR.max() < 2.0e-10
assert RADIUS_BY_RADIUS_MANUFACTURED_P2_ERROR.max() < 2.0e-10
"""
        ),
        code(
            r"""
RADIUS_BY_RADIUS_PLOT_ITERATIONS = np.array([
    0, 1, 2, 4, 8, 15, ITERATIVE_MAX_ITERATIONS,
])
fig_radius_iteration, axes = plt.subplots(
    2, 2, figsize=(14.0, 9.5))

for iteration in RADIUS_BY_RADIUS_PLOT_ITERATIONS:
    label = (
        "initial first order" if iteration == 0
        else f"iteration {iteration}")
    axes[0, 0].loglog(
        RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS,
        RADIUS_BY_RADIUS_SURFACE_RESIDUAL[:, iteration],
        "o-",
        label=label,
    )
    axes[0, 1].loglog(
        RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS,
        RADIUS_BY_RADIUS_INTERIOR_ERROR[:, iteration],
        "o-",
        label=label,
    )

for local_radius_index, radius in enumerate(
        RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS):
    axes[1, 0].plot(
        np.arange(ITERATIVE_MAX_ITERATIONS + 1),
        RADIUS_BY_RADIUS_SURFACE_RESIDUAL[
            local_radius_index],
        "o-",
        markevery=3,
        label=rf"$R/M_B={radius:.2f}$",
    )
    axes[1, 1].plot(
        np.arange(ITERATIVE_MAX_ITERATIONS + 1),
        RADIUS_BY_RADIUS_INTERIOR_ERROR[
            local_radius_index],
        "o-",
        markevery=3,
        label=rf"$R/M_B={radius:.2f}$",
    )

axes[0, 0].set_title(
    "One-worldtube value/derivative residual")
axes[0, 0].set_ylabel("fractional training residual")
axes[0, 1].set_title(
    "Held-out random interior inverse metric")
axes[0, 1].set_ylabel("ten-component relative RMS")
axes[1, 0].set_title(
    "Surface residual along each local iteration")
axes[1, 0].set_ylabel("fractional training residual")
axes[1, 1].set_title(
    "Interior prediction along each local iteration")
axes[1, 1].set_ylabel("ten-component relative RMS")
for axis in axes.flat:
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(fontsize=7)
axes[0, 0].set_xlabel(r"fitting radius $R/M_B$")
axes[0, 1].set_xlabel(r"fitting radius $R/M_B$")
axes[1, 0].set_xlabel("alternating iteration")
axes[1, 1].set_xlabel("alternating iteration")
fig_radius_iteration.suptitle(
    "Independent alternating solve at each worldtube radius")
fig_radius_iteration.tight_layout()
fig_radius_iteration.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_radius_by_radius_iteration.pdf")
fig_radius_iteration.savefig(
    OUTPUT_DIRECTORY
    / "q8_t950_radius_by_radius_iteration.png",
    dpi=180,
)
plt.show()


RADIUS_BY_RADIUS_BEST_SURFACE_ITERATION = np.argmin(
    RADIUS_BY_RADIUS_SURFACE_RESIDUAL, axis=1)
RADIUS_BY_RADIUS_BEST_INTERIOR_ITERATION = np.argmin(
    RADIUS_BY_RADIUS_INTERIOR_ERROR, axis=1)
print("\nRadius-by-radius alternating summary")
print(
    " R/M_B   best surface (iteration, residual)"
    "   best interior (iteration, error)"
)
for local_radius_index, radius in enumerate(
        RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS):
    surface_iteration = int(
        RADIUS_BY_RADIUS_BEST_SURFACE_ITERATION[
            local_radius_index])
    interior_iteration = int(
        RADIUS_BY_RADIUS_BEST_INTERIOR_ITERATION[
            local_radius_index])
    print(
        f" {radius:5.2f}"
        f"   ({surface_iteration:2d}, "
        f"{RADIUS_BY_RADIUS_SURFACE_RESIDUAL[
            local_radius_index, surface_iteration]:.6e})"
        f"   ({interior_iteration:2d}, "
        f"{RADIUS_BY_RADIUS_INTERIOR_ERROR[
            local_radius_index, interior_iteration]:.6e})"
    )
"""
        ),
        markdown(
            r"""
### Reading the radius-by-radius iteration

The surface curve measures how closely the selected truncated model satisfies
the equations on the worldtube that determined it.  The interior curve asks
whether the same coefficients predict independent metric values closer to the
black hole.  Agreement between their preferred iteration counts is stronger
evidence than improvement of either diagnostic alone.

Parameter convergence is not required for a useful early-stopped model.
However, a rapidly growing \(\|p_2(R)\|\), strong radius dependence, or
continued surface improvement accompanied by worsening interior prediction
indicates that the iteration is assigning truncation error to poorly
identified coefficients.

The manufactured fixed-point identities hold independently at all nine
fitted radii, with maximum relative errors \(9.4\times10^{-14}\) in \(p_1\) and
\(3.2\times10^{-15}\) in \(p_2\).  Thus the single-radius systems are
algebraically consistent and retain full rank.

The useful behavior is concentrated at the smaller radii.  At \(R=0.24M\),
the held-out interior error falls from \(3.65\times10^{-2}\) to
\(2.98\times10^{-2}\) and reaches its minimum after three iterations.
The neighboring \(R=0.28M\), \(0.34M\), and \(0.40M\) fits prefer roughly
three to five iterations and give minimum errors
\(3.01,3.12,3.23\times10^{-2}\), respectively.  At the smallest
\(R=0.20M\), the improvement is slower but remains monotone through iteration
thirty, reaching \(2.98\times10^{-2}\).

The surface residual usually reaches its minimum earlier—after one to three
iterations—than the interior metric.  This is direct evidence that iterating
to optimize the worldtube equations alone is not the correct stopping rule.
The analysis is deliberately capped at \(R=0.8M\), before the large-radius
instability seen in the exploratory scan can compress the plotting scale.
Within this range, the radius-by-radius result supports a small number of
alternating decontamination steps, with the interior metric used to choose
the iteration count.  It does not support unrestricted iteration at every
radius.
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
    joint_derivative_weights=JOINT_DERIVATIVE_WEIGHTS,
    joint_parameters=JOINT_PARAMETERS,
    joint_vector_condition_number=(
        JOINT_VECTOR_CONDITION_NUMBER),
    joint_clock_condition_number=(
        JOINT_CLOCK_CONDITION_NUMBER),
    joint_parameter_displacement=(
        JOINT_PARAMETER_DISPLACEMENT),
    joint_value_vector_closure_fraction=(
        JOINT_VALUE_VECTOR_CLOSURE_FRACTION),
    joint_value_vector_closure_rms=(
        JOINT_VALUE_VECTOR_CLOSURE_RMS),
    joint_derivative_vector_closure_fraction=(
        JOINT_DERIVATIVE_VECTOR_CLOSURE_FRACTION),
    joint_derivative_vector_closure_rms=(
        JOINT_DERIVATIVE_VECTOR_CLOSURE_RMS),
    joint_value_clock_closure_fraction=(
        JOINT_VALUE_CLOCK_CLOSURE_FRACTION),
    joint_value_clock_closure_rms=(
        JOINT_VALUE_CLOCK_CLOSURE_RMS),
    joint_derivative_clock_closure_fraction=(
        JOINT_DERIVATIVE_CLOCK_CLOSURE_FRACTION),
    joint_derivative_clock_closure_rms=(
        JOINT_DERIVATIVE_CLOCK_CLOSURE_RMS),
    joint_interior_total_error=(
        JOINT_INTERIOR_TOTAL_ERROR),
    joint_interior_perturbation_residual=(
        JOINT_INTERIOR_PERTURBATION_RESIDUAL),
    hybrid_parameters=HYBRID_PARAMETERS,
    hybrid_interior_total_error=(
        HYBRID_INTERIOR_TOTAL_ERROR),
    hybrid_minimum_interior_error=(
        HYBRID_MINIMUM_INTERIOR_ERROR),
    hybrid_minimum_radius=(
        HYBRID_MINIMUM_RADIUS),
    hybrid_global_index=np.asarray(
        HYBRID_GLOBAL_INDEX),
    second_order_baseline_names=np.asarray(
        SECOND_ORDER_BASELINE_NAMES),
    second_order_first_order_baselines=(
        SECOND_ORDER_FIRST_ORDER_BASELINES),
    second_order_first_order_residual_fraction=(
        SECOND_ORDER_FIRST_ORDER_RESIDUAL_FRACTION),
    second_order_quadratic_residual_fraction=(
        SECOND_ORDER_QUADRATIC_RESIDUAL_FRACTION),
    second_order_quadratic_source_fraction=(
        SECOND_ORDER_QUADRATIC_SOURCE_FRACTION),
    second_order_stage1_parameter_names=np.asarray(
        STAGE1_PARAMETER_NAMES),
    second_order_stage1_parameters=(
        SECOND_ORDER_STAGE1_PARAMETERS),
    second_order_first_order_corrections=(
        SECOND_ORDER_FIRST_ORDER_CORRECTIONS),
    second_order_stage1_gij_residual=(
        SECOND_ORDER_STAGE1_GIJ_RESIDUAL),
    second_order_stage1_condition_number=(
        SECOND_ORDER_STAGE1_CONDITION_NUMBER),
    second_order_stage1_rank=(
        SECOND_ORDER_STAGE1_RANK),
    second_order_first_order_correction_fraction=(
        SECOND_ORDER_FIRST_ORDER_CORRECTION_FRACTION),
    second_order_first_order_interior_error=(
        SECOND_ORDER_FIRST_ORDER_INTERIOR_ERROR),
    second_order_stage1_interior_error=(
        SECOND_ORDER_STAGE1_INTERIOR_ERROR),
    manufactured_stage1_frozen_error=(
        MANUFACTURED_FROZEN_ERROR),
    manufactured_stage1_augmented_error=(
        MANUFACTURED_AUGMENTED_ERROR),
    complete_second_order_parameter_names=np.asarray(
        FULL_PARAMETER_NAMES),
    complete_second_order_parameters=(
        COMPLETE_SECOND_ORDER_PARAMETERS),
    complete_second_order_stage_names=np.asarray(
        COMPLETE_SECOND_ORDER_STAGE_NAMES),
    complete_second_order_stage_candidate_fraction=(
        COMPLETE_SECOND_ORDER_STAGE_CANDIDATE_FRACTION),
    complete_second_order_stage_raw_fraction=(
        COMPLETE_SECOND_ORDER_STAGE_RAW_FRACTION),
    complete_second_order_component_fraction=(
        COMPLETE_SECOND_ORDER_COMPONENT_FRACTION),
    complete_second_order_stage2_gtt_closure=(
        COMPLETE_SECOND_ORDER_STAGE2_GTT_CLOSURE),
    complete_second_order_stage3_gti_closure=(
        COMPLETE_SECOND_ORDER_STAGE3_GTI_CLOSURE),
    complete_second_order_stage2_rank=(
        COMPLETE_SECOND_ORDER_STAGE2_RANK),
    complete_second_order_stage3_rank=(
        COMPLETE_SECOND_ORDER_STAGE3_RANK),
    complete_second_order_stage2_condition=(
        COMPLETE_SECOND_ORDER_STAGE2_CONDITION),
    complete_second_order_stage3_condition=(
        COMPLETE_SECOND_ORDER_STAGE3_CONDITION),
    complete_second_order_quadratic_interior_error=(
        COMPLETE_SECOND_ORDER_QUADRATIC_INTERIOR_ERROR),
    complete_second_order_interior_error=(
        COMPLETE_SECOND_ORDER_INTERIOR_ERROR),
    complete_second_order_interior_perturbation_residual=(
        COMPLETE_SECOND_ORDER_INTERIOR_PERTURBATION_RESIDUAL),
    complete_second_order_pointwise_relative_error=(
        COMPLETE_SECOND_ORDER_POINTWISE_RELATIVE_ERROR),
    manufactured_complete_parameter_error=(
        MANUFACTURED_COMPLETE_PARAMETER_ERROR),
    manufactured_complete_metric_closure=(
        MANUFACTURED_COMPLETE_METRIC_CLOSURE),
    polar_j2_blocks=np.asarray(
        POLAR_J2_BLOCKS),
    polar_j2_estimates=POLAR_J2_ESTIMATES,
    polar_j2_block_residuals=(
        POLAR_J2_BLOCK_RESIDUALS),
    polar_j2_rank=POLAR_J2_RANK,
    polar_j2_condition=POLAR_J2_CONDITION,
    spatial_orthogonal_estimates=(
        SPATIAL_ORTHOGONAL_ESTIMATES),
    spatial_orthogonal_in_sample=(
        SPATIAL_ORTHOGONAL_IN_SAMPLE),
    spatial_orthogonal_closure=(
        SPATIAL_ORTHOGONAL_CLOSURE),
    spatial_orthogonal_rank=(
        SPATIAL_ORTHOGONAL_RANK),
    spatial_orthogonal_condition=(
        SPATIAL_ORTHOGONAL_CONDITION),
    spatial_j3_estimates=SPATIAL_J3_ESTIMATES,
    spatial_j3_closure=SPATIAL_J3_CLOSURE,
    spatial_j3_rank=SPATIAL_J3_RANK,
    spatial_j3_condition=SPATIAL_J3_CONDITION,
    second_order_sensitivity_group_names=np.asarray(
        SECOND_ORDER_SENSITIVITY_GROUP_NAMES),
    group_only_parameters=GROUP_ONLY_PARAMETERS,
    group_only_surface_residual_fraction=(
        GROUP_ONLY_SURFACE_RESIDUAL_FRACTION),
    group_only_surface_reduction=(
        GROUP_ONLY_SURFACE_REDUCTION),
    group_only_interior_error=(
        GROUP_ONLY_INTERIOR_ERROR),
    group_only_interior_error_ratio=(
        GROUP_ONLY_INTERIOR_ERROR_RATIO),
    group_only_interior_contribution=(
        GROUP_ONLY_INTERIOR_CONTRIBUTION),
    group_only_rank=GROUP_ONLY_RANK,
    group_only_condition=GROUP_ONLY_CONDITION,
    second_order_subset_names=np.asarray(
        SECOND_ORDER_SUBSET_NAMES),
    second_order_subset_parameters=(
        SECOND_ORDER_SUBSET_PARAMETERS),
    second_order_subset_surface_fraction=(
        SECOND_ORDER_SUBSET_SURFACE_FRACTION),
    second_order_subset_interior_error=(
        SECOND_ORDER_SUBSET_INTERIOR_ERROR),
    second_order_subset_interior_ratio=(
        SECOND_ORDER_SUBSET_INTERIOR_RATIO),
    second_order_subset_rank=(
        SECOND_ORDER_SUBSET_RANK),
    second_order_subset_condition=(
        SECOND_ORDER_SUBSET_CONDITION),
    complete_ablation_contribution=(
        COMPLETE_ABLATION_CONTRIBUTION),
    complete_ablation_error_ratio=(
        COMPLETE_ABLATION_ERROR_RATIO),
    second_order_derivative_weights=(
        SECOND_ORDER_DERIVATIVE_WEIGHTS),
    second_order_derivative_stencil_error=(
        SECOND_ORDER_DERIVATIVE_STENCIL_ERROR),
    second_order_joint_subset_parameter_indices=(
        PROVISIONAL_SECOND_ORDER_PARAMETER_INDICES),
    second_order_joint_subset_parameters=(
        SECOND_ORDER_JOINT_SUBSET_PARAMETERS),
    second_order_joint_subset_value_residual=(
        SECOND_ORDER_JOINT_SUBSET_VALUE_RESIDUAL),
    second_order_joint_subset_derivative_residual=(
        SECOND_ORDER_JOINT_SUBSET_DERIVATIVE_RESIDUAL),
    second_order_joint_subset_interior_error=(
        SECOND_ORDER_JOINT_SUBSET_INTERIOR_ERROR),
    second_order_joint_subset_condition=(
        SECOND_ORDER_JOINT_SUBSET_CONDITION),
    second_order_joint_polar_parameters=(
        SECOND_ORDER_JOINT_POLAR_PARAMETERS),
    second_order_joint_polar_value_residual=(
        SECOND_ORDER_JOINT_POLAR_VALUE_RESIDUAL),
    second_order_joint_polar_derivative_residual=(
        SECOND_ORDER_JOINT_POLAR_DERIVATIVE_RESIDUAL),
    second_order_joint_polar_interior_error=(
        SECOND_ORDER_JOINT_POLAR_INTERIOR_ERROR),
    second_order_joint_polar_condition=(
        SECOND_ORDER_JOINT_POLAR_CONDITION),
    iterative_training_radius_indices=(
        ITERATIVE_TRAINING_RADIUS_INDICES),
    iterative_held_out_radius_indices=(
        ITERATIVE_HELD_OUT_RADIUS_INDICES),
    iterative_damping_values=ITERATIVE_DAMPING_VALUES,
    iterative_first_order_parameters=np.asarray([
        results["first_order"]
        for results in REAL_ITERATION_RESULTS
    ]),
    iterative_second_order_parameters=np.asarray([
        results["second_order"]
        for results in REAL_ITERATION_RESULTS
    ]),
    iterative_training_residual=np.asarray([
        results["training_residual"]
        for results in REAL_ITERATION_RESULTS
    ]),
    iterative_held_out_residual=(
        REAL_ITERATION_HELD_OUT_RESIDUAL),
    iterative_interior_error=np.asarray([
        results["interior_error"]
        for results in REAL_ITERATION_RESULTS
    ]),
    iterative_first_order_update=np.asarray([
        results["first_order_update"]
        for results in REAL_ITERATION_RESULTS
    ]),
    iterative_second_order_update=np.asarray([
        results["second_order_update"]
        for results in REAL_ITERATION_RESULTS
    ]),
    manufactured_iterative_first_order=(
        MANUFACTURED_ITERATIVE_FIRST_ORDER),
    manufactured_iterative_second_order=(
        MANUFACTURED_ITERATIVE_SECOND_ORDER),
    manufactured_iterative_first_order_error=np.asarray([
        results["first_order_error"]
        for results in MANUFACTURED_ITERATION_RESULTS
    ]),
    manufactured_iterative_second_order_error=np.asarray([
        results["second_order_error"]
        for results in MANUFACTURED_ITERATION_RESULTS
    ]),
    manufactured_iterative_training_residual=np.asarray([
        results["training_residual"]
        for results in MANUFACTURED_ITERATION_RESULTS
    ]),
    radius_by_radius_damping=RADIUS_BY_RADIUS_DAMPING,
    radius_by_radius_max_radius=(
        RADIUS_BY_RADIUS_MAX_RADIUS),
    radius_by_radius_radius_indices=(
        RADIUS_BY_RADIUS_RADIUS_INDICES),
    radius_by_radius_radii=RADIUS_BY_RADIUS_RADII,
    radius_by_radius_radius_over_small_mass=(
        RADIUS_BY_RADIUS_RADIUS_OVER_SMALL_MASS),
    radius_by_radius_first_order=(
        RADIUS_BY_RADIUS_FIRST_ORDER),
    radius_by_radius_second_order=(
        RADIUS_BY_RADIUS_SECOND_ORDER),
    radius_by_radius_surface_residual=(
        RADIUS_BY_RADIUS_SURFACE_RESIDUAL),
    radius_by_radius_interior_error=(
        RADIUS_BY_RADIUS_INTERIOR_ERROR),
    radius_by_radius_first_order_update=(
        RADIUS_BY_RADIUS_FIRST_ORDER_UPDATE),
    radius_by_radius_second_order_update=(
        RADIUS_BY_RADIUS_SECOND_ORDER_UPDATE),
    radius_by_radius_best_surface_iteration=(
        RADIUS_BY_RADIUS_BEST_SURFACE_ITERATION),
    radius_by_radius_best_interior_iteration=(
        RADIUS_BY_RADIUS_BEST_INTERIOR_ITERATION),
    radius_by_radius_manufactured_p1_error=(
        RADIUS_BY_RADIUS_MANUFACTURED_P1_ERROR),
    radius_by_radius_manufactured_p2_error=(
        RADIUS_BY_RADIUS_MANUFACTURED_P2_ERROR),
)
print(
    "\nSaved compact results to",
    OUTPUT_DIRECTORY / "q8_t950_first_order_results.npz",
)
"""
        ),
        markdown(
            r"""
## 17. Companion time-derivative analysis

The single-slice analysis above treats the second-order time jets as
instantaneous algebraic unknowns.  Ten downloaded slices from
$T=950M_{\rm tot}$ through $995M_{\rm tot}$ provide two stronger dynamical
tests, stored under [`temporal_ode_consistency/`](temporal_ode_consistency/).

First, the evolved generalized-harmonic fields give

$$
\partial_Tg_{ab}
=\beta^i\Phi_{iab}-\alpha\Pi_{ab},
\qquad
\partial_TG^{ab}
=-G^{ac}G^{bd}\partial_Tg_{cd}.
$$

At fixed local coordinates on a moving sphere,

$$
D_TG^{ab}
=\partial_TG^{ab}
+v_{\rm center}^i\partial_iG^{ab}.
$$

Fitting the first-order response to this derivative measures the rates of all
thirteen first-order coefficients.  At $R=0.20M_{\rm tot}$, the directly
measured clock acceleration and strain-rate scales are approximately

$$
\operatorname{RMS}|\ddot q^0|=3.95\times10^{-5},
\qquad
\operatorname{RMS}|M_B\dot\Lambda_{ij}|=1.25\times10^{-5}.
$$

The corresponding freely fitted seven-parameter second-order values are
$4.08\times10^{-2}$ and $7.16\times10^{-3}$.  They are therefore effective
residual-reduction directions, not physical measurements of the time
derivatives.  The direct derivative also shows that the omitted vector rates
$\dot b_i$ and $\ddot q^i$ are larger than the clock/strain rates and must be
included in any dynamical second-order subsystem.

Second, the fixed-inertial derivative retains the translation response,

$$
\partial_TG^{ab}
=-\partial_iG_{(0)}^{ab}\frac{dq^i}{dT}
+\sum_A\mathcal R_A^{ab}\frac{dp_A}{dT}.
$$

A joint full-rank fit therefore gives a metric-only measurement of the
spatial center velocity.  This is distinct from the comoving derivative,
which cancels translation and measures $\ddot q^i$.  There is no analogous
absolute $q^0$ translation measurement because the Schwarzschild background
is stationary.

## Current status

The following statements are supported by the calculations in this notebook
and its time-derivative companion:

- The corrected first-order response and implemented second-order response
  ordering pass manufactured tests to roundoff.
- V1 is the most reliable current vector estimator; radial derivatives
  modestly improve it.
- C3 is the most conservative current clock/strain estimator, but its radial
  derivative is presently more useful as closure than as a fitted equation.
- The complete triangular 43-parameter second-order solve is algebraically
  correct but fails cross-component and interior validation on this slice.
- The seven-parameter strain-rate/time-acceleration subsystem and the
  polar-$J=2$ subsystem contain the clearest second-order metric information.
- Direct time derivatives show that freely fitted time jets must be replaced
  or constrained by their dynamical measurements.
- A few radius-by-radius Gauss--Seidel sweeps can reduce interior error, but
  unrestricted iteration is not supported; the surface and interior
  diagnostics prefer different stopping times.

No block is promoted to a production matching prescription solely because it
has a small training residual.  A final subsystem must remain stable across
radii and times and must predict omitted STF blocks, spatial derivatives, and
interior metric values.
"""
        ),
    ]

    output = OUTPUT_DIRECTORY / "q8_first_order_matching.ipynb"
    nbformat.write(notebook, output)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
