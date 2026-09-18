# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Measuring the radial direction and the rest frame (Sec. 4 of np_matching.tex).

The type-D solve of Sec. 3 fixes the timelike--radial plane at each
boundary point but leaves a type-III boost within it.  This module
constructs the principal null pair in the NR tetrad (eq:kinnersley-uN),
selects the boost-family member whose radial leg is tangent to the
numerical slice (eq:boost-family, eq:transverse-member) -- defining the
measured radial direction r_hat and transverse velocity w -- and fixes
the missing radial velocity component by the sphere-wide least-squares
fit of eq:velocity-matching, from which the pointwise type-III rapidity
follows.

Components: spatial vectors live in the worldtube-adapted orthonormal
triad of scalars.adapted_tetrad_rotation (rows s, t1, t2); 4-vectors
carry (n, s, t1, t2) components, with the slice normal n as time axis and
Minkowski metric diag(-1, 1, 1, 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_SQRT2 = np.sqrt(2.0)


@dataclass(frozen=True)
class TangentBoostMember:
    """The boost-family member of eq:transverse-member.

    ``radial_direction`` is the measured radial direction r_hat and
    ``transverse_velocity`` the Eulerian velocity w of the member whose
    radial leg is tangent to the slice; ``lorentz_factor`` is
    Gamma_p = (1 - w.w)^(-1/2) = t^0.  The ``max_*`` fields are roundoff
    diagnostics of identities the construction satisfies exactly:
    the null legs stay real, g(l, k) = -1 is preserved by every rotation,
    the radial leg is unit, and w is orthogonal to r_hat (the note's
    tangency argument below eq:transverse-member).
    """

    radial_direction: np.ndarray
    transverse_velocity: np.ndarray
    lorentz_factor: np.ndarray
    max_imaginary_part: float
    max_null_product_error: float
    max_radial_norm_error: float
    max_transverse_radial_overlap: float


@dataclass(frozen=True)
class VelocityFit:
    """One sphere-wide velocity fit (eq:velocity-matching) and its
    pointwise rapidity reconstruction.

    ``coordinate_velocity`` is the single V^i of the simulation gauge;
    ``relative_residual`` is the fit residual the note names as the direct
    test that all boundary points share one rest frame.  ``local_velocity``
    is u_p(V) of eq:local-velocity in adapted components, ``rapidity`` the
    pointwise eta_p, and ``reconstruction_error`` the pointwise distance
    |u_p - (w + tanh(eta)/Gamma r_hat)| against eq:family-velocities.
    """

    coordinate_velocity: np.ndarray
    relative_residual: float
    singular_values: np.ndarray
    condition_number: float
    local_velocity: np.ndarray
    rapidity: np.ndarray
    reconstruction_error: np.ndarray


def _base_nr_tetrad(shape: tuple) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The NR tetrad of eq:surface-tetrad in adapted components."""
    ell = np.zeros(shape + (4,), dtype=complex)
    kay = np.zeros(shape + (4,), dtype=complex)
    emm = np.zeros(shape + (4,), dtype=complex)
    ell[..., 0] = ell[..., 1] = 1.0 / _SQRT2
    kay[..., 0] = 1.0 / _SQRT2
    kay[..., 1] = -1.0 / _SQRT2
    emm[..., 2] = 1.0 / _SQRT2
    emm[..., 3] = 1j / _SQRT2
    return ell, kay, emm


def _type_i_tetrad(ell, kay, emm, a):
    """eq:typeITransformation with tetrad parameter a (l is fixed)."""
    a = a[..., None]
    return (
        ell,
        kay + np.conj(a) * emm + a * np.conj(emm) + np.abs(a) ** 2 * ell,
        emm + a * ell,
    )


def _type_ii_tetrad(ell, kay, emm, b):
    """eq:typeIITransformation with tetrad parameter b (k is fixed)."""
    b = b[..., None]
    return (
        ell + np.conj(b) * emm + b * np.conj(emm) + np.abs(b) ** 2 * kay,
        kay,
        emm + b * kay,
    )


def principal_null_pair(
    a_bar: np.ndarray, b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float]:
    """The principal null pair in NR adapted components, up to type III.

    Sec. 3 writes the Kinnersley-to-NR map as type III, then type II (b),
    then type I (a); inverting it leg by leg, the Kinnersley pair in NR
    components is type I (-a) followed by type II (-b) applied to the NR
    tetrad, with the unknown type-III factor only rescaling the two legs
    -- it cancels in the tangent-member selection downstream.  The
    scalar-side solve returns a_bar = conj(a), so the tetrad parameter is
    its conjugate.  Returns (l, k, max imaginary part), the legs being
    exactly real by the conjugate pairing of the rotation terms.
    """
    a_bar = np.asarray(a_bar, dtype=complex)
    b = np.asarray(b, dtype=complex)
    shape = np.broadcast_shapes(a_bar.shape, b.shape)
    ell, kay, emm = _base_nr_tetrad(shape)
    tetrad_a = np.broadcast_to(np.conj(a_bar), shape)
    ell, kay, emm = _type_i_tetrad(ell, kay, emm, -tetrad_a)
    ell, kay, emm = _type_ii_tetrad(ell, kay, emm, -np.broadcast_to(b, shape))
    max_imaginary = float(
        max(np.max(np.abs(ell.imag)), np.max(np.abs(kay.imag)))
    )
    return ell.real, kay.real, max_imaginary


def _minkowski_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return -left[..., 0] * right[..., 0] + np.einsum(
        "...i,...i->...", left[..., 1:], right[..., 1:]
    )


def tangent_boost_member(
    a_bar: np.ndarray, b: np.ndarray
) -> TangentBoostMember:
    """Select the boost-family member tangent to the slice.

    From the principal pair, t = (l + k)/sqrt(2) and r = (l - k)/sqrt(2)
    (eq:kinnersley-uN); the type-III freedom is the boost family
    eq:boost-family.  Rescaling the null legs so their time components
    agree is a type-III transformation and makes r^0 = 0 -- exactly the
    tangency condition g(r(eta), n) = 0 -- so the member
    (t_0p, r_0p) of eq:transverse-member is read off directly:
    Gamma_p = t^0, w_p = t_spatial/t^0, r_hat = r_spatial.
    """
    ell, kay, max_imaginary = principal_null_pair(a_bar, b)
    if np.any(ell[..., 0] <= 0.0) or np.any(kay[..., 0] <= 0.0):
        raise RuntimeError("principal null legs are not future-directed")
    boost = np.sqrt(kay[..., 0] / ell[..., 0])
    ell = boost[..., None] * ell
    kay = kay / boost[..., None]
    t_leg = (ell + kay) / _SQRT2
    r_leg = (ell - kay) / _SQRT2

    lorentz_factor = t_leg[..., 0]
    transverse_velocity = t_leg[..., 1:] / lorentz_factor[..., None]
    radial = r_leg[..., 1:]
    radial_norm = np.linalg.norm(radial, axis=-1)
    radial_direction = radial / radial_norm[..., None]
    if np.any(
        np.einsum("...i,...i->...", transverse_velocity, transverse_velocity)
        >= 1.0
    ):
        raise RuntimeError("decoded transverse velocity is not subluminal")
    return TangentBoostMember(
        radial_direction=radial_direction,
        transverse_velocity=transverse_velocity,
        lorentz_factor=lorentz_factor,
        max_imaginary_part=max_imaginary,
        max_null_product_error=float(
            np.max(np.abs(_minkowski_product(ell, kay) + 1.0))
        ),
        max_radial_norm_error=float(np.max(np.abs(radial_norm - 1.0))),
        max_transverse_radial_overlap=float(
            np.max(
                np.abs(
                    np.einsum(
                        "...i,...i->...",
                        transverse_velocity,
                        radial_direction,
                    )
                )
            )
        ),
    )


def fit_coordinate_velocity(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    adapted_rotation: np.ndarray,
    spatial_metric: np.ndarray,
    lapse: np.ndarray,
    shift: np.ndarray,
) -> VelocityFit:
    """Fit one coordinate velocity V to a sphere (eq:velocity-matching).

    The Eulerian observer measures u_p(V) = (V + beta)/alpha
    (eq:local-velocity); in adapted components u_p = R L^T (V + beta) /
    alpha with gamma = L L^T and R the adapted rotation.  Consistency with
    eq:family-velocities requires P u_p(V) = w_p at every point, an
    overdetermined linear system for the three components of V solved by
    least squares; its residual directly tests the one-rest-frame
    assumption.  The remaining radial component gives the pointwise
    type-III rapidity, tanh(eta_p) = Gamma_p (u_p(V) - w_p) . r_hat.
    """
    radial_direction = np.asarray(radial_direction, dtype=float)
    transverse_velocity = np.asarray(transverse_velocity, dtype=float)
    adapted_rotation = np.asarray(adapted_rotation, dtype=float)
    spatial_metric = np.asarray(spatial_metric, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    point_shape = radial_direction.shape[:-1]
    if transverse_velocity.shape != point_shape + (3,):
        raise ValueError("transverse_velocity does not match the grid")
    if adapted_rotation.shape != point_shape + (3, 3):
        raise ValueError("adapted_rotation does not match the grid")
    if spatial_metric.shape != point_shape + (3, 3):
        raise ValueError("spatial_metric does not match the grid")
    if lapse.shape != point_shape or shift.shape != point_shape + (3,):
        raise ValueError("lapse or shift does not match the grid")

    cholesky_transpose = np.linalg.cholesky(spatial_metric).swapaxes(-1, -2)
    to_adapted = (
        np.einsum("...ij,...jk->...ik", adapted_rotation, cholesky_transpose)
        / lapse[..., None, None]
    )
    projector = np.eye(3) - np.einsum(
        "...i,...j->...ij", radial_direction, radial_direction
    )
    design_pointwise = np.einsum("...ij,...jk->...ik", projector, to_adapted)
    target_pointwise = transverse_velocity - np.einsum(
        "...ij,...j->...i", design_pointwise, shift
    )
    design = design_pointwise.reshape(-1, 3)
    target = target_pointwise.reshape(-1)
    coordinate_velocity, *_ = np.linalg.lstsq(design, target, rcond=None)
    singular_values = np.linalg.svd(design, compute_uv=False)
    relative_residual = float(
        np.linalg.norm(design @ coordinate_velocity - target)
        / max(np.linalg.norm(target), 1.0e-300)
    )

    local_velocity = np.einsum(
        "...ij,...j->...i", to_adapted, coordinate_velocity + shift
    )
    transverse_speed_squared = np.einsum(
        "...i,...i->...", transverse_velocity, transverse_velocity
    )
    lorentz_factor = 1.0 / np.sqrt(1.0 - transverse_speed_squared)
    tanh_rapidity = lorentz_factor * np.einsum(
        "...i,...i->...",
        local_velocity - transverse_velocity,
        radial_direction,
    )
    if np.any(np.abs(tanh_rapidity) >= 1.0):
        raise RuntimeError("inferred radial rapidity is not physical")
    rapidity = np.arctanh(tanh_rapidity)
    reconstructed = (
        transverse_velocity
        + (np.tanh(rapidity) / lorentz_factor)[..., None] * radial_direction
    )
    reconstruction_error = np.linalg.norm(
        reconstructed - local_velocity, axis=-1
    )
    return VelocityFit(
        coordinate_velocity=coordinate_velocity,
        relative_residual=relative_residual,
        singular_values=singular_values,
        condition_number=float(singular_values[0] / singular_values[-1]),
        local_velocity=local_velocity,
        rapidity=rapidity,
        reconstruction_error=reconstruction_error,
    )


@dataclass(frozen=True)
class SliceTilt:
    """Model-time tilt across one excision sphere (Sec. slicing-caveat).

    ``delta_t`` is the per-point model time offset t - t_0 (zero-mean:
    the constant of integration is the choice of t_0, absorbed into the
    moments); ``gradient`` the measured covector of
    eq:slice-tilt-gradient; ``fit_residual`` the relative residual of the
    tangential potential fit -- the direct test that the measured
    gradient field is integrable over the sphere at this order.
    """

    delta_t: np.ndarray
    gradient: np.ndarray
    fit_residual: float


def slice_tilt(
    positions: np.ndarray,
    coordinate_velocity: np.ndarray,
    lapse: np.ndarray,
    shift: np.ndarray,
    spatial_metric: np.ndarray,
    root_f: np.ndarray,
) -> SliceTilt:
    """Integrate eq:slice-tilt-gradient over one sphere.

    The model-time gradient along any slice-tangent displacement is
    measured, d_i t = -Gamma_u (gamma u_p)_i / sqrt(f) with
    u_p = (V + beta)/alpha (eq:local-velocity), so delta_t is the
    potential of a known covector field.  It is recovered by a
    least-squares potential fit (linear plus quadratic monomials in the
    coordinate position; delta_t is position-linear at leading order,
    the quadratic terms absorb the sphere-scale variation of the
    coefficient) matching the sphere-tangential gradient components.
    ``positions`` are coordinate positions of the sphere points relative
    to any fixed origin.
    """
    positions = np.asarray(positions, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    spatial_metric = np.asarray(spatial_metric, dtype=float)
    root_f = np.asarray(root_f, dtype=float)
    local_velocity = (
        np.asarray(coordinate_velocity, dtype=float) + shift
    ) / lapse[..., None]
    velocity_lower = np.einsum(
        "...ij,...j->...i", spatial_metric, local_velocity
    )
    speed_squared = np.einsum("...i,...i->...", velocity_lower, local_velocity)
    lorentz = 1.0 / np.sqrt(1.0 - speed_squared)
    gradient = -(lorentz / root_f)[..., None] * velocity_lower

    # Tangential potential fit: only displacements on the coordinate
    # sphere are sampled, so only the projected gradient is constrained.
    center = np.mean(positions.reshape(-1, 3), axis=0)
    relative = positions - center
    normal = relative / np.linalg.norm(relative, axis=-1)[..., None]
    projector = np.eye(3) - np.einsum("...i,...j->...ij", normal, normal)
    pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    basis_values = [relative[..., i] for i in range(3)]
    basis_gradients = [
        np.broadcast_to(np.eye(3)[i], relative.shape).copy() for i in range(3)
    ]
    for i, j in pairs:
        basis_values.append(relative[..., i] * relative[..., j])
        grad = np.zeros_like(relative)
        grad[..., i] += relative[..., j]
        grad[..., j] += relative[..., i]
        basis_gradients.append(grad)
    tangential = lambda field: np.einsum(  # noqa: E731
        "...ij,...j->...i", projector, field
    )
    design = np.stack(
        [tangential(grad).reshape(-1) for grad in basis_gradients], axis=-1
    )
    data = tangential(gradient).reshape(-1)
    coefficients, *_ = np.linalg.lstsq(design, data, rcond=None)
    fit_residual = float(
        np.linalg.norm(design @ coefficients - data)
        / max(np.linalg.norm(data), 1.0e-300)
    )
    delta_t = np.einsum("a...,a->...", np.stack(basis_values), coefficients)
    delta_t = delta_t - np.mean(delta_t)
    return SliceTilt(
        delta_t=delta_t, gradient=gradient, fit_residual=fit_residual
    )


@dataclass(frozen=True)
class InvariantRapidity:
    """Pointwise rapidity from the curvature gradient (note Sec. 3.2).

    ``rapidity`` is eta_p of eq:invariant-boost-condition in the
    boost-family parametrization anchored at the tangent member;
    ``normal_gradient`` is n.grad(K) and ``radial_gradient`` the adapted
    r_hat.DK entering the ratio; ``max_abs_tanh`` is the sanity bound
    (must stay below one for a physical boost).
    """

    rapidity: np.ndarray
    normal_gradient: np.ndarray
    radial_gradient: np.ndarray
    max_abs_tanh: float


def invariant_rapidity(
    member: TangentBoostMember,
    adapted_rotation: np.ndarray,
    spatial_metric: np.ndarray,
    lapse: np.ndarray,
    shift: np.ndarray,
    spatial_gradient: np.ndarray,
    time_derivative: np.ndarray,
) -> InvariantRapidity:
    """The invariant boost of note Sec. 3.2 (eq:invariant-boost-condition).

    The rest-frame member of the boost family annihilates the spacetime
    gradient of the Kretschmann field:  with
    t(eta) = cosh(eta) t_0p + sinh(eta) r_0p anchored at the tangent
    member (t_0p = Gamma (n + w), r_0p = (0, r_hat)),

        t(eta_p) . grad K = 0
        =>  tanh(eta_p) = - Gamma (n.grad K + w . DK) / (r_hat . DK),

    the tangent-member form of eq:boost-from-gradient.  Inputs:
    ``spatial_gradient`` is D_i K in inertial coordinate (covector)
    components -- for K = 16 Re I this is 16 Re grad_invariant_i of the
    slice cache; ``time_derivative`` is dt K at fixed inertial
    coordinates, supplied by the on-shell single-slice evaluation (the
    vacuum-ADM synthetic step) or by a time series;
    n.grad K = (dt K - beta^i D_i K)/alpha.  Covector components in the
    adapted triad are L^{-1} D K rotated by ``adapted_rotation``
    (gamma = L L^T), so that contractions with the member's triad
    vectors are Euclidean.
    """
    spatial_gradient = np.asarray(spatial_gradient, dtype=float)
    time_derivative = np.asarray(time_derivative, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    normal_gradient = (
        time_derivative - np.einsum("...i,...i->...", shift, spatial_gradient)
    ) / lapse

    cholesky = np.linalg.cholesky(np.asarray(spatial_metric, dtype=float))
    orthonormal = np.linalg.solve(cholesky, spatial_gradient[..., None])[..., 0]
    adapted = np.einsum(
        "...ij,...j->...i",
        np.asarray(adapted_rotation, dtype=float),
        orthonormal,
    )
    radial_gradient = np.einsum(
        "...i,...i->...", member.radial_direction, adapted
    )
    tangential = np.einsum(
        "...i,...i->...", member.transverse_velocity, adapted
    )
    tanh_rapidity = (
        -member.lorentz_factor
        * (normal_gradient + tangential)
        / radial_gradient
    )
    max_abs_tanh = float(np.max(np.abs(tanh_rapidity)))
    if max_abs_tanh >= 1.0:
        raise RuntimeError("invariant-boost rapidity is not physical")
    return InvariantRapidity(
        rapidity=np.arctanh(tanh_rapidity),
        normal_gradient=normal_gradient,
        radial_gradient=radial_gradient,
        max_abs_tanh=max_abs_tanh,
    )
