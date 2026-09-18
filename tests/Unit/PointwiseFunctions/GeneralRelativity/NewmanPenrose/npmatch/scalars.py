# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Newman--Penrose scalars on the worldtube: Secs. 1--2 of np_matching.tex.

The five complex Weyl scalars always occupy the LAST axis, indexed
A = 0..4, and every function is vectorized over arbitrary leading point
axes.  Each public function cites the equation label of
``notes/np_matching.tex`` that it implements.  Nothing in this module
touches simulation data.

Frame conventions.  Spatial components are taken in an orthonormal triad.
The worldtube-adapted triad has rows (s, theta_hat, phi_hat): s is the unit
outward normal of the excision sphere within the slice and the tangent legs
come from the gauge-coordinate angular directions, so the NR tetrad of
eq:surface-tetrad is

    l = (n + s)/sqrt(2),  k = (n - s)/sqrt(2),
    m = (theta_hat + i phi_hat)/sqrt(2).

The dyad convention fixes the type-III spin chi that the note leaves free
(harmless in w^-, eq:incoming-weyl); see deviations.md entry D2.
"""

from __future__ import annotations

import numpy as np


def adapted_tetrad_rotation(
    gamma: np.ndarray, directions: np.ndarray
) -> np.ndarray:
    """Rows ``(s, theta_hat, phi_hat)`` of the surface triad.

    ``gamma[..., 3, 3]`` is the covariant spatial metric in simulation
    coordinates and ``directions[..., 3]`` the Euclidean unit normal of the
    coordinate sphere.  The rows are components in the Cholesky-orthonormal
    triad of gamma (gamma = L L^T): the sphere normal is the covector
    proportional to ``directions``, whose orthonormal components are
    ``L^{-1} d``, normalized.  The tangent legs start from the standard
    gauge-coordinate spherical vectors ``partial_theta x`` and
    ``partial_phi x``.  Their coordinate-vector components are carried to
    the orthonormal triad with ``L^T`` and orthonormalized there with the
    orientation ``s cross theta_hat = phi_hat``.  At a coordinate pole the
    azimuth is fixed to zero as an explicit angular-patch convention.
    """
    gamma = np.asarray(gamma, dtype=float)
    directions = np.broadcast_to(
        np.asarray(directions, dtype=float), gamma.shape[:-2] + (3,)
    )
    cholesky = np.linalg.cholesky(gamma)
    cholesky_inverse = np.linalg.inv(cholesky)
    radial = np.einsum("...ij,...j->...i", cholesky_inverse, directions)
    radial = radial / np.linalg.norm(radial, axis=-1)[..., None]

    cylindrical_radius = np.linalg.norm(directions[..., :2], axis=-1)
    safe_radius = np.maximum(cylindrical_radius, 1.0e-300)
    cos_phi = np.where(
        cylindrical_radius > 1.0e-14,
        directions[..., 0] / safe_radius,
        1.0,
    )
    sin_phi = np.where(
        cylindrical_radius > 1.0e-14,
        directions[..., 1] / safe_radius,
        0.0,
    )
    theta_coordinate = np.stack(
        (
            directions[..., 2] * cos_phi,
            directions[..., 2] * sin_phi,
            -cylindrical_radius,
        ),
        axis=-1,
    )
    phi_coordinate = np.stack(
        (-sin_phi, cos_phi, np.zeros_like(cos_phi)), axis=-1
    )

    theta = np.einsum("...ji,...j->...i", cholesky, theta_coordinate)
    phi = np.einsum("...ji,...j->...i", cholesky, phi_coordinate)
    theta -= radial * np.einsum("...i,...i->...", radial, theta)[..., None]
    theta /= np.linalg.norm(theta, axis=-1)[..., None]
    phi -= radial * np.einsum("...i,...i->...", radial, phi)[..., None]
    phi -= theta * np.einsum("...i,...i->...", theta, phi)[..., None]
    phi /= np.linalg.norm(phi, axis=-1)[..., None]
    orientation = np.einsum("...i,...i->...", np.cross(radial, theta), phi)
    phi = np.where((orientation < 0.0)[..., None], -phi, phi)
    return np.stack((radial, theta, phi), axis=-2)


def rotate_symmetric(tensor: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    """Components of a rank-2 tensor in the triad whose basis vectors are
    the rows of ``rotation``."""
    return np.einsum("...ik,...kl,...jl->...ij", rotation, tensor, rotation)


def rotate_vector(vector: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    """Components of a vector in the triad whose basis vectors are the rows
    of ``rotation``."""
    return np.einsum("...ij,...j->...i", rotation, vector)


def coordinate_vector_from_adapted(
    vector: np.ndarray, adapted_rotation: np.ndarray, gamma: np.ndarray
) -> np.ndarray:
    """Simulation-coordinate components of an adapted-triad vector.

    A coordinate vector v has Cholesky-orthonormal components L^T v
    (gamma = L L^T) and adapted components R L^T v, so the inverse map is
    v_coord = L^{-T} R^T v_adapted.
    """
    cholesky_transpose = np.linalg.cholesky(
        np.asarray(gamma, dtype=float)
    ).swapaxes(-1, -2)
    in_cholesky_triad = np.einsum(
        "...ji,...j->...i", np.asarray(adapted_rotation, dtype=float), vector
    )
    return np.linalg.solve(cholesky_transpose, in_cholesky_triad[..., None])[
        ..., 0
    ]


def psi_from_q(q: np.ndarray) -> np.ndarray:
    """Weyl scalars from the tidal tensor, eq:scalars-from-Q.

    ``q[..., 3, 3]`` holds the components of Q = E + iB in the adapted
    triad (rows s, t1, t2), and m = (t1 + i t2)/sqrt(2).  Explicitly,

        Psi_0 = Q(m, m),          Psi_1 = -Q(s, m)/sqrt(2),
        Psi_2 = Q(s, s)/2,        Psi_3 = +Q(s, mbar)/sqrt(2),
        Psi_4 = Q(mbar, mbar).
    """
    q = np.asarray(q, dtype=complex)
    psi = np.empty(q.shape[:-2] + (5,), dtype=complex)
    transverse_difference = 0.5 * (q[..., 1, 1] - q[..., 2, 2])
    psi[..., 0] = transverse_difference + 1j * q[..., 1, 2]
    psi[..., 1] = -0.5 * (q[..., 0, 1] + 1j * q[..., 0, 2])
    psi[..., 2] = 0.5 * q[..., 0, 0]
    psi[..., 3] = 0.5 * (q[..., 0, 1] - 1j * q[..., 0, 2])
    psi[..., 4] = transverse_difference - 1j * q[..., 1, 2]
    return psi


def q_from_psi(psi: np.ndarray) -> np.ndarray:
    """Tidal tensor from the Weyl scalars: the inverse of eq:scalars-from-Q.

    Returns the symmetric trace-free ``Q`` in the adapted triad; five
    complex scalars are exactly the ten real components (Sec. 1).
    """
    psi = np.asarray(psi, dtype=complex)
    p0, p1, p2, p3, p4 = np.moveaxis(psi, -1, 0)
    q = np.empty(psi.shape[:-1] + (3, 3), dtype=complex)
    q[..., 0, 0] = 2.0 * p2
    q[..., 0, 1] = q[..., 1, 0] = p3 - p1
    q[..., 0, 2] = q[..., 2, 0] = 1j * (p1 + p3)
    q[..., 1, 1] = -p2 + 0.5 * (p0 + p4)
    q[..., 1, 2] = q[..., 2, 1] = 0.5j * (p4 - p0)
    q[..., 2, 2] = -p2 - 0.5 * (p0 + p4)
    return q


def w_minus(psi0: np.ndarray, adapted_rotation: np.ndarray) -> np.ndarray:
    """The incoming characteristic field of eq:incoming-weyl,

        w^-_ij = 2 ( conj(Psi0) m_i m_j + Psi0 mbar_i mbar_j ),

    returned as real components in the Cholesky-orthonormal triad, with
    m = (t1 + i t2)/sqrt(2) built from the rows of ``adapted_rotation``.

    The field is real, symmetric, trace-free, and tangent to the cut.
    It is invariant under the dyad spin m -> e^{i chi} m: Psi0 rotates
    oppositely to mbar_i mbar_j (eq:type-III-scalars), so w^- depends on
    the dyad only through the tangent plane it spans -- two Psi0's can
    be compared through w^- without sharing a dyad convention, provided
    they refer to the same cut normal s.  Pointwise
    |w^-|_F^2 = 8 |Psi0|^2, so relative errors in w^- and Psi0 coincide.
    """
    psi0 = np.asarray(psi0, dtype=complex)
    rotation = np.asarray(adapted_rotation, dtype=float)
    m = (rotation[..., 1, :] + 1j * rotation[..., 2, :]) / np.sqrt(2.0)
    mm = np.einsum("...i,...j->...ij", m, m)
    return np.real(
        2.0
        * (
            np.conj(psi0)[..., None, None] * mm
            + psi0[..., None, None] * np.conj(mm)
        )
    )


def w_plus(psi4: np.ndarray, adapted_rotation: np.ndarray) -> np.ndarray:
    """The outgoing characteristic field, the spin-independent partner of
    :func:`w_minus`,

        w^+_ij = 2 (Psi4 m_i m_j + conj(Psi4) mbar_i mbar_j).

    The result is returned as real components in the Cholesky-orthonormal
    triad.  Under a dyad spin ``m -> exp(i chi) m`` the scalar transforms as
    ``Psi4 -> exp(-2 i chi) Psi4``, so the tensor is unchanged.  Pointwise
    ``|w^+|_F^2 = 8 |Psi4|^2``.
    """
    psi4 = np.asarray(psi4, dtype=complex)
    rotation = np.asarray(adapted_rotation, dtype=float)
    m = (rotation[..., 1, :] + 1j * rotation[..., 2, :]) / np.sqrt(2.0)
    mm = np.einsum("...i,...j->...ij", m, m)
    return np.real(
        2.0
        * (
            psi4[..., None, None] * mm
            + np.conj(psi4)[..., None, None] * np.conj(mm)
        )
    )


def w_plus_screen(psi4: np.ndarray) -> np.ndarray:
    """Components of ``w^+`` in the local screen basis of ``psi4``.

    The first axis is the radial leg and the last two axes span the screen.
    A point-dependent dyad spin rotates these components orthogonally, so a
    Frobenius-norm least-squares problem built from data and model columns in
    the same local screen basis is independent of that arbitrary spin.
    """
    psi4 = np.asarray(psi4, dtype=complex)
    out = np.zeros(psi4.shape + (3, 3), dtype=float)
    out[..., 1, 1] = 2.0 * psi4.real
    out[..., 1, 2] = out[..., 2, 1] = -2.0 * psi4.imag
    out[..., 2, 2] = -2.0 * psi4.real
    return out


def type_i(psi: np.ndarray, a_bar: np.ndarray | complex) -> np.ndarray:
    """Type-I null rotation of the scalars, eq:type-I-scalars.

    The parameter is written ``a_bar`` because the scalars transform with
    the conjugate of the tetrad parameter a in eq:typeITransformation.
    """
    psi = np.asarray(psi, dtype=complex)
    a = np.asarray(a_bar, dtype=complex)
    p0, p1, p2, p3, p4 = np.moveaxis(psi, -1, 0)
    out = np.empty(np.broadcast_shapes(p0.shape, a.shape) + (5,), dtype=complex)
    out[..., 0] = p0
    out[..., 1] = p1 + a * p0
    out[..., 2] = p2 + 2.0 * a * p1 + a**2 * p0
    out[..., 3] = p3 + 3.0 * a * p2 + 3.0 * a**2 * p1 + a**3 * p0
    out[..., 4] = (
        p4 + 4.0 * a * p3 + 6.0 * a**2 * p2 + 4.0 * a**3 * p1 + a**4 * p0
    )
    return out


def type_ii(psi: np.ndarray, b: np.ndarray | complex) -> np.ndarray:
    """Type-II null rotation of the scalars, eq:type-II-scalars."""
    psi = np.asarray(psi, dtype=complex)
    b = np.asarray(b, dtype=complex)
    p0, p1, p2, p3, p4 = np.moveaxis(psi, -1, 0)
    out = np.empty(np.broadcast_shapes(p0.shape, b.shape) + (5,), dtype=complex)
    out[..., 4] = p4
    out[..., 3] = p3 + b * p4
    out[..., 2] = p2 + 2.0 * b * p3 + b**2 * p4
    out[..., 1] = p1 + 3.0 * b * p2 + 3.0 * b**2 * p3 + b**3 * p4
    out[..., 0] = (
        p0 + 4.0 * b * p1 + 6.0 * b**2 * p2 + 4.0 * b**3 * p3 + b**4 * p4
    )
    return out


def type_iii(
    psi: np.ndarray,
    eta: np.ndarray | float = 0.0,
    chi: np.ndarray | float = 0.0,
) -> np.ndarray:
    """Type-III boost/spin of the scalars, eq:type-III-scalars:
    Psi_A -> exp((2 - A)(eta + i chi)) Psi_A."""
    psi = np.asarray(psi, dtype=complex)
    exponent = np.asarray(eta, dtype=complex) + 1j * np.asarray(
        chi, dtype=complex
    )
    weights = 2.0 - np.arange(5.0)
    return psi * np.exp(weights * exponent[..., None])


def invariants(psi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The curvature invariants I and J from the scalars, eq:invariants.

    I = Psi0 Psi4 - 4 Psi1 Psi3 + 3 Psi2^2 and J is the determinant of the
    Weyl-scalar Hankel matrix, expanded explicitly.  Both are unchanged
    under every tetrad transformation of Sec. 2.
    """
    psi = np.asarray(psi, dtype=complex)
    p0, p1, p2, p3, p4 = np.moveaxis(psi, -1, 0)
    invariant_i = p0 * p4 - 4.0 * p1 * p3 + 3.0 * p2**2
    invariant_j = (
        p0 * p2 * p4 - p4 * p1**2 - p0 * p3**2 + 2.0 * p1 * p2 * p3 - p2**3
    )
    return invariant_i, invariant_j


def invariants_from_q(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """I and J from the tidal tensor, eq:invariants-from-Q.

    ``q[..., 3, 3]`` holds components in any orthonormal triad (equivalently
    the mixed components Q^i_j):  I = Q_ij Q^ij / 2  and  J = -det(Q^i_j)/2.
    """
    q = np.asarray(q, dtype=complex)
    invariant_i = 0.5 * np.einsum("...ij,...ji->...", q, q)
    invariant_j = -0.5 * np.linalg.det(q)
    return invariant_i, invariant_j
