# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Tidal moments from the Coulomb channel with exact type-D kinematics.

The order-two target of eq:second-order-psi0 needs the quadrupole moments
(E, B) of the tide.  Reading them from the leaving mode Psi4 closes an
unstable loop when the target is imposed online (the injected entering mode
returns as Psi4 and is re-emitted amplified by the slice's boost factor).
This module reads them from the Coulomb channel instead.

On a single sphere the electric tide is degenerate with a displacement of
the sphere in the Coulomb scalar alone: Psi2^K = -M/r^3 + (1/2) E_nn and an
unknown areal-radius field r(n) absorb each other.  The normal derivative
breaks the degeneracy because a static tide does not change it at first
order (profile e_L = 1 in Table tab:responses).  For a type-D background
the proper radial derivative of the Coulomb scalar is

    dPsi2/dl = B(r) = sqrt(1 - 2M/r) 3M / r^4

in the rest frame of the hole, and the slice measures its component along
the sphere normal s.  The frame registration relates the two exactly at
every point: with the tangent member (Gamma, w, r_hat) and the invariant
rapidity eta, the rest-frame radial unit vector has slice components
sinh(eta) Gamma w + cosh(eta) r_hat, so

    d_s Psi2^K = B(r) [cosh(eta) (r_hat . s) + sinh(eta) Gamma (w . s)].

Solving this for r pointwise gives the areal radius of every face point
free of the tide, whatever the offset, shape and motion of the coordinate
sphere; nothing is linearized in the kinematics.  The tide is then the
pointwise excess Psi2^K + M/r^3, fitted to the Psi2 slot of the tidal model
(all channels, boosted to the slice) for the ten real moments, with l = 0
and l = 1 nuisance terms absorbing a mass mismatch and a residual dipole.
Psi0 and Psi4 enter only through the invariants, quadratically.

The map r -> B(r) has its turning point at r = 9M/4, so the decode needs the
sphere outside about 2.4M; its sensitivity is 1/|(M/r)/f - 4|.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import profiles
import psi0
import scalars
import tides

TURNING_POINT_OVER_MASS = 2.25


@dataclass(frozen=True)
class CoulombDecode:
    """The decoded tide of one sphere.

    ``radius`` is the areal radius per point from the normal derivative,
    ``components`` the five complex direct components (H11, H22, H12, H13,
    H23) of E + iB in the convention of :func:`psi0.direct_tide_scalar_columns`,
    ``nuisance`` the fitted l = 0 and l = 1 Coulomb offsets (complex),
    ``relative_residual`` the weighted relative residual of the Coulomb
    excess fit, ``valid`` whether every point stayed on the outer branch of
    B(r) and the radius solve converged.
    """

    radius: np.ndarray
    components: np.ndarray
    nuisance: np.ndarray
    relative_residual: float
    valid: bool
    newton_residual: float


def background_radial_derivative(radius: np.ndarray, mass: float) -> np.ndarray:
    """B(r) = sqrt(1 - 2M/r) 3M/r^4, the proper radial derivative of the
    Coulomb scalar of a hole at rest."""
    radius = np.asarray(radius, dtype=float)
    return np.sqrt(1.0 - 2.0 * mass / radius) * 3.0 * mass / radius**4


def normal_derivative_factor(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    lorentz_factor: np.ndarray,
    rapidity: np.ndarray,
) -> np.ndarray:
    """cosh(eta) (r_hat . s) + sinh(eta) Gamma (w . s) for the tangent member
    in adapted-triad components (row 0 is the sphere normal s)."""
    radial_direction = np.asarray(radial_direction, dtype=float)
    transverse_velocity = np.asarray(transverse_velocity, dtype=float)
    return (
        np.cosh(rapidity) * radial_direction[..., 0]
        + np.sinh(rapidity)
        * np.asarray(lorentz_factor)
        * transverse_velocity[..., 0]
    )


def radius_from_normal_derivative(
    d_s_coulomb_real: np.ndarray,
    factor: np.ndarray,
    initial_radius: np.ndarray,
    mass: float,
    *,
    iterations: int = 40,
) -> tuple[np.ndarray, bool, float]:
    """Solve B(r) |factor| = |d_s Re Psi2^K| for r on the outer branch by
    Newton's method from ``initial_radius``.  Returns (radius, valid,
    max relative Newton residual)."""
    magnitude = np.abs(np.asarray(d_s_coulomb_real, dtype=float)) / np.abs(
        factor
    )
    floor = (TURNING_POINT_OVER_MASS + 0.01) * mass
    radius = np.maximum(np.asarray(initial_radius, dtype=float).copy(), floor)
    for _ in range(iterations):
        f = 1.0 - 2.0 * mass / radius
        b = np.sqrt(f) * 3.0 * mass / radius**4
        db = b * ((mass / radius) / f - 4.0) / radius
        step = (b - magnitude) / db
        step = np.clip(step, -0.2 * radius, 0.2 * radius)
        radius = np.maximum(radius - step, floor)
    residual = float(
        np.max(
            np.abs(background_radial_derivative(radius, mass) - magnitude)
            / magnitude
        )
    )
    valid = bool(np.all(radius > floor * (1.0 + 1.0e-12)) and residual < 1.0e-8)
    return radius, valid, residual


def coulomb_tide_columns(
    radial_direction: np.ndarray,
    radius: np.ndarray,
    mass: float,
) -> np.ndarray:
    """The response of the Coulomb scalar -3J/I to the ten real unit
    moments: five electric direct components then five magnetic ones, shape
    ``(10, ...)``.

    -3J/I is a tetrad invariant, so the tide's contribution to it is the
    rest-frame value (1/2) Q(R_hat, R_hat) = (1/2) e_L H(R_hat, R_hat) along
    the rest-frame radial direction, with no boost mixing whatsoever; the
    boost only rearranges the scalars between the slots of a tetrad, not
    the invariant.  ``radial_direction`` is the measured radial direction in
    the Cholesky frame, which the tidal model identifies with the rest-frame
    radial direction (deviation D5), and ``radius`` enters only through the
    profile e_L (= 1 for the quadrupole)."""
    zero = np.zeros((3, 3))
    r_hat = np.asarray(radial_direction, dtype=float)
    columns = []
    for kind in ("electric", "magnetic"):
        for component in tides.DIRECT_STF_COMPONENT_TENSORS:
            electric = component if kind == "electric" else zero
            magnetic = component if kind == "magnetic" else zero
            tide_rest = profiles.quadrupole_tide_tensor(
                electric, magnetic, r_hat, radius, mass, transverse_only=False
            )
            columns.append(
                0.5 * np.einsum("...i,...ij,...j->...", r_hat, tide_rest, r_hat)
            )
    return np.stack(columns)


def decode_tidal_moments_from_coulomb(
    coulomb: np.ndarray,
    d_s_coulomb_real: np.ndarray,
    member_radial_direction: np.ndarray,
    member_transverse_velocity: np.ndarray,
    member_lorentz_factor: np.ndarray,
    rapidity: np.ndarray,
    measured_radius: np.ndarray,
    radial_cholesky: np.ndarray,
    transverse_cholesky: np.ndarray,
    mass: float,
    point_weights: np.ndarray | None = None,
) -> CoulombDecode:
    """Decode the ten quadrupole moments of one sphere from the Coulomb
    scalar ``coulomb`` and the real part of its derivative along the sphere
    normal ``d_s_coulomb_real`` (both pointwise), given the registered frame:
    the tangent member in adapted-triad components, the same vectors in the
    Cholesky frame, the invariant rapidity and the measured radius as the
    starting point of the radius solve."""
    coulomb = np.asarray(coulomb, dtype=complex)
    n_point = coulomb.shape[-1]
    factor = normal_derivative_factor(
        member_radial_direction,
        member_transverse_velocity,
        member_lorentz_factor,
        rapidity,
    )
    radius, valid, newton_residual = radius_from_normal_derivative(
        d_s_coulomb_real, factor, measured_radius, mass
    )
    excess = coulomb + mass / radius**3

    columns = coulomb_tide_columns(
        radial_direction=radial_cholesky, radius=radius, mass=mass
    ).reshape(10, n_point)
    n_hat = np.asarray(radial_cholesky).reshape(n_point, 3)
    design = np.zeros((2 * n_point, 18))
    design[:n_point, :10] = columns.real.T
    design[n_point:, :10] = columns.imag.T
    design[:n_point, 10] = 1.0
    design[n_point:, 11] = 1.0
    design[:n_point, 12:15] = n_hat
    design[n_point:, 15:18] = n_hat
    target = np.concatenate([excess.real.reshape(-1), excess.imag.reshape(-1)])
    weights = (
        np.ones(n_point)
        if point_weights is None
        else np.asarray(point_weights, dtype=float)
    )
    root = np.sqrt(np.concatenate([weights, weights]))
    solution, *_ = np.linalg.lstsq(
        design * root[:, None], target * root, rcond=None
    )
    fitted = design @ solution
    residual = float(
        np.linalg.norm((fitted - target) * root)
        / max(np.linalg.norm(target * root), 1.0e-300)
    )
    return CoulombDecode(
        radius=radius.reshape(coulomb.shape),
        components=solution[:5] + 1j * solution[5:10],
        nuisance=np.concatenate(
            [
                [solution[10] + 1j * solution[11]],
                solution[12:15] + 1j * solution[15:18],
            ]
        ),
        relative_residual=residual,
        valid=valid,
        newton_residual=newton_residual,
    )
