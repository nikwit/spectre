# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Fitting the tides (Sec. 5 of np_matching.tex) -- and finding N1.

The measured invariant combination 3J/I carries the quadrupole moments
through eq:coulomb-quadrupole,

    3J/I = M/r^3 - (1/2)(E + iB)_nn + O(eps^3).

This module implements Sec. 5 literally: the radial derivative along the
measured radial direction with the measured radius as Jacobian
(eq:radial-chain-rule, eq:radial-derivative), the Kretschmann
combination eq:kretschmann-radius,

    R~^2 = -(1/9) ( r d_r Re I + 3 Re I ),

the pointwise moment fields E_nn = 2 (R~ - Re 3J/I) and
B_nn = -2 Im 3J/I, and their least-squares fit over the boundary points
against the l=2 contractions r_hat<i r_hat j> (eq:design; pointwise, no
harmonic projection).

FINDING N1 (deviations.md; verified analytically, in the closure tests,
and on q8): the electric sector of this construction is DEGENERATE.
Through linear order in the tide,

    Re I = 3 A^2 - (3/4) E_nn^2      with  A := 3J/I,

so I, J -- and therefore the radius field r = (M/A)^(1/3) of eq:radius
and every spatial derivative of any of them -- are functions of the
single measured field A.  In the ratio of eq:radial-derivative the
gradients of A cancel between numerator and Jacobian, and the
Kretschmann field evaluates to R~ = A + O(eps^4) identically: the
contaminated Coulomb again, not the de-tided background.  E_nn as
defined above vanishes at working order and the electric moments are
unobservable in the invariant channel; eq:kretschmann-radius de-tides
only if the Jacobian d_i r comes from OUTSIDE the invariants (a
coordinate map or a metric model).  The magnetic sector is unaffected:
Im 3J/I has no radius contamination and determines B cleanly.

All derivative inputs are element-spectral gradients from the slice
cache; everything here is single-sphere.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations

import kinnersley
import numpy as np


def stf_basis() -> np.ndarray:
    """Frobenius-orthonormal basis of real symmetric trace-free tensors.

    The convention matches the older pipeline
    (scalar_radius_tide_iteration.STF_BASIS) so fitted coefficient vectors
    are directly comparable: diag(1,-1,0)/sqrt(2), diag(1,1,-2)/sqrt(6),
    then the symmetrized (01), (02), (12) pairs.
    """
    basis = [
        np.diag([1.0, -1.0, 0.0]) / np.sqrt(2.0),
        np.diag([1.0, 1.0, -2.0]) / np.sqrt(6.0),
    ]
    for first, second in ((0, 1), (0, 2), (1, 2)):
        tensor = np.zeros((3, 3))
        tensor[first, second] = tensor[second, first] = 1.0 / np.sqrt(2.0)
        basis.append(tensor)
    return np.asarray(basis)


STF_BASIS = stf_basis()


def stf_from_components(components: np.ndarray) -> np.ndarray:
    """Build an STF tensor from its five direct Cartesian components.

    The last axis of ``components`` is ordered as
    ``(X11, X22, X12, X13, X23)`` in one common orthonormal Cartesian
    frame.  Symmetry supplies the lower triangle and trace-freeness fixes
    ``X33 = -X11 - X22``.  Complex input represents ``E + i B`` directly.
    """
    components = np.asarray(components)
    if components.shape[-1] != 5:
        raise ValueError("expected components (X11, X22, X12, X13, X23)")
    tensor = np.zeros(components.shape[:-1] + (3, 3), dtype=components.dtype)
    tensor[..., 0, 0] = components[..., 0]
    tensor[..., 1, 1] = components[..., 1]
    tensor[..., 2, 2] = -components[..., 0] - components[..., 1]
    tensor[..., 0, 1] = tensor[..., 1, 0] = components[..., 2]
    tensor[..., 0, 2] = tensor[..., 2, 0] = components[..., 3]
    tensor[..., 1, 2] = tensor[..., 2, 1] = components[..., 4]
    return tensor


def stf_components(tensor: np.ndarray) -> np.ndarray:
    """Return ``(X11, X22, X12, X13, X23)`` from an STF tensor."""
    tensor = np.asarray(tensor)
    if tensor.shape[-2:] != (3, 3):
        raise ValueError("expected a tensor with final shape (3, 3)")
    return np.stack(
        (
            tensor[..., 0, 0],
            tensor[..., 1, 1],
            tensor[..., 0, 1],
            tensor[..., 0, 2],
            tensor[..., 1, 2],
        ),
        axis=-1,
    )


# Unit variations of the five stored components.  These are not an
# auxiliary orthonormal STF basis: column A is literally the derivative
# with respect to X11, X22, X12, X13, or X23, respectively.
DIRECT_STF_COMPONENT_TENSORS = stf_from_components(np.eye(5))


def rank3_stf_from_components(components: np.ndarray) -> np.ndarray:
    """Rank-three STF tensor from (111,112,113,122,123,222,223).

    Symmetry fills permutations; the three trace conditions supply
    133=-111-122, 233=-112-222, and 333=-113-223. Complex components
    can represent E_ijk + (4i/3) B_ijk for the transverse octupole.
    """
    components = np.asarray(components)
    if components.shape[-1] != 7:
        raise ValueError("expected seven independent octupole components")
    tensor = np.zeros(components.shape[:-1] + (3, 3, 3), dtype=components.dtype)
    indices = (
        (0, 0, 0),
        (0, 0, 1),
        (0, 0, 2),
        (0, 1, 1),
        (0, 1, 2),
        (1, 1, 1),
        (1, 1, 2),
        (0, 2, 2),
        (1, 2, 2),
        (2, 2, 2),
    )
    values = [components[..., i] for i in range(7)] + [
        -components[..., 0] - components[..., 3],
        -components[..., 1] - components[..., 5],
        -components[..., 2] - components[..., 6],
    ]
    for index, value in zip(indices, values):
        for permutation in set(permutations(index)):
            tensor[(...,) + permutation] = value
    return tensor


DIRECT_OCTUPOLE_COMPONENT_TENSORS = rank3_stf_from_components(np.eye(7))


def rank3_stf_basis() -> np.ndarray:
    """Frobenius-orthonormal basis of the seven real rank-3 STF tensors.

    Deterministic Gram-Schmidt over the same seed set as the older
    pipeline (pv_octupole_responses.rank3_stf_basis), so fitted octupole
    coefficient vectors are directly comparable with the old audits.
    """
    seeds = (
        (0, 0, 0),
        (0, 0, 1),
        (0, 0, 2),
        (0, 1, 1),
        (0, 1, 2),
        (1, 1, 1),
        (1, 1, 2),
    )
    identity = np.eye(3)
    basis = []
    for indices in seeds:
        symmetric = np.zeros((3, 3, 3))
        for permutation in set(permutations(indices)):
            symmetric[permutation] = 1.0
        trace = np.einsum("iik->k", symmetric)
        candidate = (
            symmetric
            - (
                np.einsum("ij,k->ijk", identity, trace)
                + np.einsum("ik,j->ijk", identity, trace)
                + np.einsum("jk,i->ijk", identity, trace)
            )
            / 5.0
        )
        for previous in basis:
            candidate -= np.einsum("ijk,ijk", candidate, previous) * previous
        norm = np.linalg.norm(candidate)
        if norm <= 1.0e-12:
            raise RuntimeError("rank-3 STF seed set is linearly dependent")
        basis.append(candidate / norm)
    result = np.asarray(basis)
    if np.max(np.abs(np.einsum("aiik->ak", result))) > 1.0e-13:
        raise RuntimeError("rank-3 STF projection failed")
    return result


STF_BASIS_RANK3 = rank3_stf_basis()


def stf_angular_functions(direction: np.ndarray) -> np.ndarray:
    """The l=2 design functions r_hat<i r_hat j> of eq:design, contracted
    with the STF basis: S_a = basis_a,ij r_hat^i r_hat^j."""
    direction = np.asarray(direction, dtype=float)
    return np.einsum("...i,aij,...j->...a", direction, STF_BASIS, direction)


def radius_and_gradient(
    invariant_i: np.ndarray,
    invariant_j: np.ndarray,
    grad_invariant_i: np.ndarray,
    grad_invariant_j: np.ndarray,
    mass: float,
) -> tuple[np.ndarray, np.ndarray]:
    """The measured radius field and its coordinate gradient.

    Sec. 5 uses that r is itself a measured scalar field through
    eq:radius; the chain rule on 3 ln r = ln(M I / 3J) gives

        d_i r = (r/3) ( d_i I / I - d_i J / J ).
    """
    invariant_i = np.asarray(invariant_i, dtype=complex)
    invariant_j = np.asarray(invariant_j, dtype=complex)
    radius = kinnersley.background_radius(invariant_i, invariant_j, mass)
    grad_radius = (radius[..., None] / 3.0) * (
        np.asarray(grad_invariant_i, dtype=complex) / invariant_i[..., None]
        - np.asarray(grad_invariant_j, dtype=complex) / invariant_j[..., None]
    )
    return radius, grad_radius


def radial_derivative(
    gradient: np.ndarray,
    grad_radius: np.ndarray,
    radial_direction: np.ndarray,
) -> np.ndarray:
    """eq:radial-derivative: d_r f = (r_hat^i d_i f) / (r_hat^j d_j r).

    ``gradient`` and ``grad_radius`` are coordinate gradients (covector
    components) and ``radial_direction`` the measured radial direction in
    coordinate vector components; the denominator supplies the Jacobian
    between the simulation coordinates and the background radius, so no
    coordinate map is needed.
    """
    numerator = np.einsum(
        "...i,...i->...", np.asarray(radial_direction), np.asarray(gradient)
    )
    denominator = np.einsum(
        "...i,...i->...",
        np.asarray(radial_direction),
        np.asarray(grad_radius),
    )
    return numerator / denominator


def detided_coulomb(
    invariant_i: np.ndarray,
    dr_re_invariant_i: np.ndarray,
    radius: np.ndarray,
) -> np.ndarray:
    """eq:kretschmann-radius: R~^2 = -(1/9)(r d_r Re I + 3 Re I).

    The note's derivation cancels the R*E cross term with the TRUE radius
    as Jacobian.  With the measured radius of eq:radius -- the only one
    available -- the result is R~ = 3J/I + O(eps^4) instead (finding N1 in
    the module docstring): the returned field is the contaminated Coulomb,
    kept here as the literal Sec.-5 construction and as the measured
    demonstration of the degeneracy.
    """
    square = (
        -(
            np.asarray(radius, dtype=float) * np.asarray(dr_re_invariant_i)
            + 3.0 * np.real(invariant_i)
        )
        / 9.0
    )
    if np.any(square <= 0.0):
        raise RuntimeError("the de-tided Coulomb square is not positive")
    return np.sqrt(square)


@dataclass(frozen=True)
class QuadrupoleFit:
    """One sphere's moment fit from the invariant channel (eq:design).

    ``electric``/``magnetic`` are the five STF_BASIS coefficients of
    E_ij and B_ij; the ``*_nn`` fields are the measured pointwise moment
    fields the fit consumes; the residuals are relative L2 over the
    sphere.  ``corrected_coulomb`` and ``corrected_radius`` implement
    eq:radius-correction with the fitted moments (step 2 of the Sec. 6.3
    iteration); ``detided_coulomb`` is the Kretschmann field, so the
    difference between the two is a consistency diagnostic of the whole
    section.
    """

    electric: np.ndarray
    magnetic: np.ndarray
    electric_nn: np.ndarray
    magnetic_nn: np.ndarray
    electric_residual: float
    magnetic_residual: float
    condition_number: float
    detided_coulomb: np.ndarray
    corrected_coulomb: np.ndarray
    corrected_radius: np.ndarray


def fit_quadrupole_moments(
    invariant_i: np.ndarray,
    invariant_j: np.ndarray,
    dr_re_invariant_i: np.ndarray,
    radial_direction: np.ndarray,
    mass: float,
) -> QuadrupoleFit:
    """Fit the ten quadrupole moments on one sphere (Sec. 5, literal).

    ``radial_direction`` holds the measured radial direction in the same
    orthonormal triad in which the moment tensors are to be expressed;
    ``dr_re_invariant_i`` is the measured radial derivative of Re I from
    eq:radial-derivative.  Using the leading-order r_hat here is exact at
    working order (the eps^2 correction of the radial vector first enters
    at eps^4, eq:design).

    Because of finding N1 the ``electric`` result carries no tide signal
    (it measures the O(eps^4) leftover of the degeneracy and is reported
    as that diagnostic); ``magnetic`` is a genuine determination.
    """
    invariant_i = np.asarray(invariant_i, dtype=complex)
    invariant_j = np.asarray(invariant_j, dtype=complex)
    coulomb_measured = 3.0 * invariant_j / invariant_i
    radius = kinnersley.background_radius(invariant_i, invariant_j, mass).real
    background = detided_coulomb(invariant_i, dr_re_invariant_i, radius)

    electric_nn = 2.0 * (background - coulomb_measured.real)
    magnetic_nn = -2.0 * coulomb_measured.imag

    design = stf_angular_functions(radial_direction).reshape(-1, 5)
    singular_values = np.linalg.svd(design, compute_uv=False)
    electric, *_ = np.linalg.lstsq(design, electric_nn.reshape(-1), rcond=None)
    magnetic, *_ = np.linalg.lstsq(design, magnetic_nn.reshape(-1), rcond=None)
    electric_residual = float(
        np.linalg.norm(design @ electric - electric_nn.reshape(-1))
        / max(np.linalg.norm(electric_nn), 1.0e-300)
    )
    magnetic_residual = float(
        np.linalg.norm(design @ magnetic - magnetic_nn.reshape(-1))
        / max(np.linalg.norm(magnetic_nn), 1.0e-300)
    )

    # eq:radius-correction: remove the fitted tide from the invariant
    # Coulomb (step 2 of the Sec. 6.3 pass).
    fitted_nn = design @ (electric + 1j * magnetic)
    corrected_coulomb = coulomb_measured + 0.5 * fitted_nn.reshape(
        coulomb_measured.shape
    )
    corrected_radius = np.power(mass / corrected_coulomb, 1.0 / 3.0)
    return QuadrupoleFit(
        electric=electric,
        magnetic=magnetic,
        electric_nn=electric_nn,
        magnetic_nn=magnetic_nn,
        electric_residual=electric_residual,
        magnetic_residual=magnetic_residual,
        condition_number=float(singular_values[0] / singular_values[-1]),
        detided_coulomb=background,
        corrected_coulomb=corrected_coulomb,
        corrected_radius=corrected_radius,
    )
