# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Leading order: the Kinnersley frame from curvature (Sec. 3 of np_matching.tex).

The measured invariants hand over the Kinnersley Coulomb scalar and the
background radius without knowing the Kinnersley frame; the Psi1 and Psi2
equations of the transformed type-D scalars then fix the two null rotations
connecting the Kinnersley frame to the NR tetrad, and the remaining three
equations are O(eps^2) constraints (the tide meters).  All functions are
pure NumPy, vectorized over leading point axes.

Sign convention: Psi2^K is the Kinnersley Coulomb scalar itself, so for
Schwarzschild Psi2^K = -M/r^3 < 0.  (The older pipeline's "coulomb" is
-Psi2^K = +M/r^3.)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scalars


@dataclass(frozen=True)
class TypeDRotation:
    """The Sec.-3 null rotations fixed by the Psi1 and Psi2 equations.

    ``a_bar`` and ``b`` are the effective parameters with the undetermined
    type-III factor absorbed (Sec. 3: the type-D model cannot fix type III);
    only ``x = a_bar * b`` and ``predicted_psi`` are convention-free.
    ``predicted_psi`` are the five type-D scalars in the NR tetrad; its
    slots 1 and 2 match the numerical scalars by construction (roundoff),
    and slots 0, 3, 4 differ from them at O(eps^2) -- the constraints.
    """

    a_bar: np.ndarray
    b: np.ndarray
    x: np.ndarray
    predicted_psi: np.ndarray


def coulomb_scalar(
    invariant_i: np.ndarray, invariant_j: np.ndarray
) -> np.ndarray:
    """The Kinnersley Coulomb scalar Psi2^K = -3J/I, eq:typeD-IJ."""
    return (
        -3.0
        * np.asarray(invariant_j, dtype=complex)
        / np.asarray(invariant_i, dtype=complex)
    )


def coulomb_scalar_alternatives(
    invariant_i: np.ndarray,
    invariant_j: np.ndarray,
    reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """The other two Coulomb estimators of eq:typeD-IJ's remark.

    Returns (-sqrt(I/3), (-J)^(1/3)) with the square-root sign and cube-root
    branch chosen nearest ``reference`` (normally the -3J/I value).  The
    note states all three agree at O(eps^2); their measured spread is a
    consistency diagnostic.
    """
    reference = np.asarray(reference, dtype=complex)
    root = np.sqrt(np.asarray(invariant_i, dtype=complex) / 3.0)
    candidates = np.stack((-root, root))
    pick = np.argmin(np.abs(candidates - reference), axis=0)
    alternative_sqrt = np.take_along_axis(candidates, pick[None], axis=0)[0]

    principal = np.power(-np.asarray(invariant_j, dtype=complex), 1.0 / 3.0)
    turn = np.exp(2j * np.pi / 3.0)
    candidates = np.stack((principal, principal * turn, principal * turn**2))
    pick = np.argmin(np.abs(candidates - reference), axis=0)
    alternative_cbrt = np.take_along_axis(candidates, pick[None], axis=0)[0]
    return alternative_sqrt, alternative_cbrt


def background_radius(
    invariant_i: np.ndarray, invariant_j: np.ndarray, mass: float
) -> np.ndarray:
    """The background areal radius r = (M I / 3J)^(1/3), eq:radius.

    Returned complex (principal branch): the real part is the radius
    estimate and the imaginary part is an O(eps^2) diagnostic (it measures
    the odd tide through Im 3J/I).
    """
    argument = (
        mass
        * np.asarray(invariant_i, dtype=complex)
        / (3.0 * np.asarray(invariant_j, dtype=complex))
    )
    return np.power(argument, 1.0 / 3.0)


def kinnersley_scalars(coulomb: np.ndarray) -> np.ndarray:
    """The Kinnersley-frame scalars (0, 0, Psi2^K, 0, 0), eq:kinnersley-scalars."""
    coulomb = np.asarray(coulomb, dtype=complex)
    psi = np.zeros(coulomb.shape + (5,), dtype=complex)
    psi[..., 2] = coulomb
    return psi


def type_d_scalars(
    coulomb: np.ndarray,
    a_bar: np.ndarray | complex,
    b: np.ndarray | complex,
) -> np.ndarray:
    """Type-D scalars in a generic tetrad: the Psi'' system of Sec. 3.

    Implemented literally as the note's composition: type II (parameter b,
    fixing k), then type I (parameter a, fixing l), applied to the
    Kinnersley scalars.  The type-III transformation drops out
    (eq:kinnersley-scalars is invariant under it).
    """
    return scalars.type_i(
        scalars.type_ii(kinnersley_scalars(coulomb), b), a_bar
    )


def solve_type_d_rotation(
    psi: np.ndarray, coulomb: np.ndarray
) -> TypeDRotation:
    """Fix (a_bar, b) from the Psi1 and Psi2 equations of Sec. 3.

    The Psi2'' equation gives the quadratic of Sec. 3,

        a_bar b = (-6 +- sqrt(36 - 24 (1 - Psi2''/Psi2^K))) / 12 ,

    whose near-zero root aligns the ingoing/outgoing null legs with those of
    the NR tetrad (the other root swaps them).  The Psi1'' equation,
    Psi1'' = 3 b (1 + 2 a_bar b) Psi2^K, then gives b, and a_bar = x/b.
    The remaining equations (slots 0, 3, 4 of ``predicted_psi``) are the
    O(eps^2) constraints.
    """
    psi = np.asarray(psi, dtype=complex)
    coulomb = np.asarray(coulomb, dtype=complex)
    ratio = psi[..., 2] / coulomb
    discriminant = np.sqrt(36.0 - 24.0 * (1.0 - ratio))
    roots = np.stack(
        ((-6.0 + discriminant) / 12.0, (-6.0 - discriminant) / 12.0)
    )
    pick = np.argmin(np.abs(roots), axis=0)
    x = np.take_along_axis(roots, pick[None], axis=0)[0]
    if np.any(np.abs(x) > 0.5):
        raise RuntimeError(
            "The aligning root a_bar*b should be near zero; a magnitude "
            "above 0.5 means the tetrad is longitudinally exchanged and the "
            "root selection is unreliable here."
        )

    b = psi[..., 1] / (3.0 * (1.0 + 2.0 * x) * coulomb)
    aligned = np.abs(b) == 0.0
    if np.any(aligned & (np.abs(x) > 0.0)):
        raise RuntimeError(
            "Psi1 vanishes while a_bar*b does not; the Psi1/Psi2 equation "
            "pair cannot fix the rotation at such points."
        )
    a_bar = np.where(aligned, 0.0 + 0.0j, x / np.where(aligned, 1.0, b))
    return TypeDRotation(
        a_bar=a_bar,
        b=b,
        x=x,
        predicted_psi=type_d_scalars(coulomb, a_bar, b),
    )


def pull_back(
    psi: np.ndarray,
    a_bar: np.ndarray | complex,
    b: np.ndarray | complex,
) -> np.ndarray:
    """Apply the inverse of the Kinnersley rotations (Sec. 6.1, first step).

    The forward map is type II (b) followed by type I (a_bar); each family
    is an additive group, so the exact inverse is type I (-a_bar) followed
    by type II (-b).  On exact type-D scalars this lands on
    eq:kinnersley-scalars; on data it exposes the longitudinal residuals
    eps^2 Psi1^(2), eps^2 Psi3^(2) of Sec. 6.1.
    """
    return scalars.type_ii(
        scalars.type_i(psi, -np.asarray(a_bar)), -np.asarray(b)
    )


def psi0_leading(coulomb: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The leading-order boundary value Psi0^NR = 6 b^2 Psi2^K,
    eq:psi0-leading.  Equal to slot 0 of the solve's ``predicted_psi``."""
    return (
        6.0
        * np.asarray(b, dtype=complex) ** 2
        * np.asarray(coulomb, dtype=complex)
    )
