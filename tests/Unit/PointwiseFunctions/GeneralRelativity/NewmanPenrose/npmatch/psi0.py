# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Constructing Psi0 in the NR frame (Sec. 6.2 of np_matching.tex),
with the transverse-channel moment fit that finding N1 makes necessary.

The boundary value eq:psi0-nr is

    Psi0^NR = e^{2(eta+i chi)} Psi0^model + 6 b^2 Psi2^QK
              + b^4 e^{-2(eta+i chi)} Psi4^QK + O(eps^3):

the Coulomb slot keeps the measured Psi2^QK = -3J/I with its tidal
shift, and the tide enters only through its transverse channel
(eq:second-order-psi0) -- the longitudinal and vector channels are
already carried by the measured Coulomb and the measured frame
rotations, so a target that added them would double count.

Implementation note (deviation D5): the Lorentz map from the hole rest
frame to the numerical triad is applied as the SO(3,C) action on the
self-dual Q -- radial boost (r_hat, eta) followed by the transverse
boost (w), the tensorial image of the Sec.-4 registration -- rather
than through the note's type-III/II/I scalar ladder.  The two agree at
working order (the map is THE Lorentz transformation); the tensor route
carries the e^{2 eta} Psi0^model and b^4 e^{-2 eta} Psi4^QK terms of
eq:psi0-nr exactly and needs no spin bookkeeping (chi cancels between
the model contraction and the target's dyad).

Moment determination (deviation D4, justified by finding N1; note
appendix app:transverse-fit): the invariant channel cannot see the
electric moments, so E (and, as a cross-check, B) is fitted from the
quasi-Kinnersley Psi4 of eq:qk-psi4 -- after the Sec.-6.1 corrections
the frame has spent the tide's longitudinal and vector channels, and
the transverse pair is all that survives.  The design columns are the
full unit-moment tides differenced through the whole measurement chain
(solve, pull-back, Newton corrections), eq:transverse-design, so the
type-III registration and the frame's own response to each moment are
supplied automatically.  Numerical Psi0 never enters the fitted rows;
it reaches the fit only through the invariants in the Coulomb (the
conditional caveat of the old ledgers).
"""

from __future__ import annotations

from dataclasses import dataclass

import kinnersley
import numpy as np
import profiles
import qk
import scalars
import tides


@dataclass(frozen=True)
class TransverseMomentFit:
    """Moments from the quasi-Kinnersley Psi4 (deviation D4)."""

    electric: np.ndarray
    magnetic: np.ndarray
    newton_steps: int
    relative_residual: float
    condition_number: float


@dataclass(frozen=True)
class ThirdOrderMomentFit:
    """Joint eps^2 + eps^3 moments from the quasi-Kinnersley Psi4.

    All 34 parameters are determined in ONE least-squares problem: the
    eps^3 channels are not orthogonal to the quadrupole columns (the
    induction sector leaks into the magnetic quadrupole, the tilt term
    into the octupole), so stacking a second fit on frozen quadrupoles
    would bias exactly the moments being added.  Block iteration on the
    same objective converges to this joint solution.
    """

    electric: np.ndarray
    magnetic: np.ndarray
    electric_octupole: np.ndarray
    magnetic_octupole: np.ndarray
    electric_dot: np.ndarray
    magnetic_dot: np.ndarray
    newton_steps: int
    relative_residual: float
    condition_number: float

    def vector(self) -> np.ndarray:
        """The 34 coefficients in column order."""
        return np.concatenate(
            (
                self.electric,
                self.magnetic,
                self.electric_octupole,
                self.magnetic_octupole,
                self.electric_dot,
                self.magnetic_dot,
            )
        )


def _cross_matrix(axis: np.ndarray) -> np.ndarray:
    matrix = np.zeros(axis.shape[:-1] + (3, 3))
    matrix[..., 0, 1] = -axis[..., 2]
    matrix[..., 0, 2] = axis[..., 1]
    matrix[..., 1, 0] = axis[..., 2]
    matrix[..., 1, 2] = -axis[..., 0]
    matrix[..., 2, 0] = -axis[..., 1]
    matrix[..., 2, 1] = axis[..., 0]
    return matrix


def self_dual_boost(axis: np.ndarray, rapidity: np.ndarray) -> np.ndarray:
    """SO(3,C) representation of a pure boost acting on Q = E + iB.

    Convention: a boost transverse to a real direction N sends it to
    cosh(eta) (N + i velocity x N) -- the convention validated on q8 by
    the old pipeline's trajectory registration.
    """
    axis = np.asarray(axis, dtype=float)
    rapidity = np.asarray(rapidity, dtype=float)
    parallel = np.einsum("...i,...j->...ij", axis, axis)
    return (
        parallel
        + np.cosh(rapidity)[..., None, None] * (np.eye(3) - parallel)
        + 1j * np.sinh(rapidity)[..., None, None] * _cross_matrix(axis)
    )


def rest_to_slice_map(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    rapidity: np.ndarray,
) -> np.ndarray:
    """The SO(3,C) image of the Sec.-4 registration (deviation D5).

    Radial boost by the fitted rapidity along r_hat, then the transverse
    boost with velocity w: components of a self-dual tensor in the hole
    rest frame are carried to the numerical orthonormal triad.  All
    inputs are in one fixed orthonormal triad (the Cholesky triad here).
    """
    transverse_velocity = np.asarray(transverse_velocity, dtype=float)
    speed = np.linalg.norm(transverse_velocity, axis=-1)
    # At a zero of the transverse speed the boost is the identity and any
    # regular axis represents it.
    safe_speed = np.maximum(speed, 1.0e-300)
    axis = np.where(
        (speed > 1.0e-14)[..., None],
        transverse_velocity / safe_speed[..., None],
        np.asarray(radial_direction, dtype=float),
    )
    transverse_map = self_dual_boost(axis, np.arctanh(speed))
    radial_map = self_dual_boost(radial_direction, np.asarray(rapidity))
    return np.einsum("...ij,...jk->...ik", transverse_map, radial_map)


def tide_scalar_columns(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    rapidity: np.ndarray,
    radius: np.ndarray,
    adapted_rotation: np.ndarray,
    mass: float,
    *,
    transverse_only: bool,
) -> np.ndarray:
    """NR-frame scalar columns of the ten unit-moment tides.

    For each STF basis element the tide tensor of eq:tidal-tensor
    (profiles of Table tab:responses) is built in the rest frame, boosted
    with the Sec.-4 registration, rotated to the adapted triad, and read
    off as the five Weyl scalars.  Directions and velocities are
    Cholesky-triad components; ``radius`` is the measured background
    radius.  Returns shape (10, ..., 5): electric basis columns first,
    then magnetic.

    The two variants serve different roles.  The FULL columns
    (``transverse_only=False``) are the physical unit-moment
    perturbations fed through the measurement chain in the fit design
    eq:transverse-design -- the chain redistributes their longitudinal
    and vector channels into the frame, and the columns must contain
    them for that response to be captured.  The TRANSVERSE columns are
    the Psi0-target slot of eq:psi0-nr: there the vector channel enters
    through the measured corrections b^(2) and the longitudinal channel
    through the measured Coulomb, so a full-tide target slot would
    double count them.
    """
    boost = rest_to_slice_map(radial_direction, transverse_velocity, rapidity)
    rest_tensors = _quadrupole_rest_tensors(
        radial_direction, radius, mass, transverse_only=transverse_only
    )
    return _scalar_columns(rest_tensors, boost, adapted_rotation)


def direct_tide_scalar_columns(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    rapidity: np.ndarray,
    radius: np.ndarray,
    adapted_rotation: np.ndarray,
    mass: float,
    *,
    transverse_only: bool,
) -> np.ndarray:
    """Five scalar columns for the direct complex components of ``E+iB``.

    Column order is ``(H11, H22, H12, H13, H23)`` with
    ``H = E + i B`` and ``H33 = -H11 - H22``.  Each returned column is the
    response to a real unit value of the corresponding component; complex
    least-squares coefficients supply its electric and magnetic parts.

    This direct parameterization is valid for the transverse response used
    in Eq. (44), because ``e_T = b_T = f`` and the self-dual Lorentz action
    is complex-linear.  It deliberately does not introduce the normalized
    ``tides.STF_BASIS`` used by the older analysis pipeline.
    """
    if not transverse_only:
        raise ValueError(
            "direct complex components are implemented only for the "
            "transverse Eq. (44) response"
        )
    boost = rest_to_slice_map(radial_direction, transverse_velocity, rapidity)
    zero = np.zeros((3, 3))
    rest_tensors = [
        profiles.quadrupole_tide_tensor(
            component,
            zero,
            radial_direction,
            radius,
            mass,
            transverse_only=True,
        )
        for component in tides.DIRECT_STF_COMPONENT_TENSORS
    ]
    return _scalar_columns(rest_tensors, boost, adapted_rotation)


def _scalar_columns(rest_tensors, boost, adapted_rotation) -> np.ndarray:
    """Boost each rest-frame tide tensor to the slice (deviation D5) and
    read off the five scalars in the adapted triad."""
    columns = []
    for tide_rest in rest_tensors:
        tide_slice = np.einsum(
            "...ik,...kl,...jl->...ij", boost, tide_rest, boost
        )
        columns.append(
            scalars.psi_from_q(
                scalars.rotate_symmetric(tide_slice, adapted_rotation)
            )
        )
    return np.stack(columns)


def _quadrupole_rest_tensors(
    radial_direction, radius, mass, *, transverse_only
) -> list:
    zero = np.zeros((3, 3))
    tensors = []
    for kind in ("electric", "magnetic"):
        for basis in tides.STF_BASIS:
            electric = basis if kind == "electric" else zero
            magnetic = basis if kind == "magnetic" else zero
            tensors.append(
                profiles.quadrupole_tide_tensor(
                    electric,
                    magnetic,
                    radial_direction,
                    radius,
                    mass,
                    transverse_only=transverse_only,
                )
            )
    return tensors


def third_order_tide_scalar_columns(
    radial_direction: np.ndarray,
    transverse_velocity: np.ndarray,
    rapidity: np.ndarray,
    radius: np.ndarray,
    adapted_rotation: np.ndarray,
    mass: float,
    *,
    transverse_only: bool,
    slice_tilt: np.ndarray,
    calibration: tuple[float, float] = profiles.HORIZON_CALIBRATION,
) -> np.ndarray:
    """The 34 unit-moment scalar columns of the eps^3 rung.

    Column order: 5 electric + 5 magnetic quadrupoles (as in
    :func:`tide_scalar_columns`), 7 electric + 7 magnetic octupoles
    (rank-3 STF basis), 5 electric-dot + 5 magnetic-dot moments.  Each
    dotted column carries all three structures the appendix assigns to
    one dotted moment: the flat induction sector, the near-zone M-terms
    (in the fixed ``calibration``), and the slice-tilt term
    eq:slice-tilt -- the measured delta_t times the quadrupole
    structure.  ``slice_tilt`` is the per-point delta_t from
    restframe.slice_tilt (its constant is absorbed into the undotted
    moments).
    """
    boost = rest_to_slice_map(radial_direction, transverse_velocity, rapidity)
    tilt = np.asarray(slice_tilt, dtype=float)[..., None, None]
    zero2 = np.zeros((3, 3))
    zero3 = np.zeros((3, 3, 3))
    rest_tensors = _quadrupole_rest_tensors(
        radial_direction, radius, mass, transverse_only=transverse_only
    )
    for kind in ("electric", "magnetic"):
        for basis in tides.STF_BASIS_RANK3:
            electric = basis if kind == "electric" else zero3
            magnetic = basis if kind == "magnetic" else zero3
            rest_tensors.append(
                profiles.octupole_tide_tensor(
                    electric,
                    magnetic,
                    radial_direction,
                    radius,
                    mass,
                    transverse_only=transverse_only,
                )
            )
    for kind in ("electric", "magnetic"):
        for basis in tides.STF_BASIS:
            electric = basis if kind == "electric" else zero2
            magnetic = basis if kind == "magnetic" else zero2
            rest_tensors.append(
                profiles.induction_tide_tensor(
                    electric,
                    magnetic,
                    radial_direction,
                    radius,
                    mass,
                    transverse_only=transverse_only,
                )
                + profiles.nearzone_tide_tensor(
                    electric,
                    magnetic,
                    radial_direction,
                    radius,
                    mass,
                    calibration=calibration,
                    transverse_only=transverse_only,
                )
                + tilt
                * profiles.quadrupole_tide_tensor(
                    electric,
                    magnetic,
                    radial_direction,
                    radius,
                    mass,
                    transverse_only=transverse_only,
                )
            )
    return _scalar_columns(rest_tensors, boost, adapted_rotation)


def quasi_kinnersley_psi4(
    psi: np.ndarray, *, newton_steps: int = 2
) -> np.ndarray:
    """Psi4 of the quasi-Kinnersley frame: the measured field of eq:qk-psi4.

    The full measurement chain -- invariants, type-D solve, pull-back,
    Newton corrections (eq:qk-construction) -- applied to the scalars.
    Numerical Psi0 enters only through the invariants in the Coulomb.
    """
    invariant_i, invariant_j = scalars.invariants(psi)
    coulomb = kinnersley.coulomb_scalar(invariant_i, invariant_j)
    fit = kinnersley.solve_type_d_rotation(psi, coulomb)
    frame = qk.quasi_kinnersley(psi, fit, coulomb, steps=newton_steps)
    return frame.psi[..., 4]


def fit_transverse_moments(
    psi: np.ndarray,
    tide_columns: np.ndarray,
    *,
    amplitude_fraction: float = 1.0e-3,
    newton_steps: int = 2,
) -> TransverseMomentFit:
    """Fit the ten moments from the quasi-Kinnersley Psi4
    (eq:transverse-design; deviation D4).

    Each design column is the response of eq:qk-psi4 to one full
    unit-moment tide, differenced through the whole chain, so the
    type-III registration and the reabsorption of the moment's
    longitudinal and vector imprints into the frame are included
    exactly.  Numerical Psi0 is absent from the fitted rows and
    retained as a prediction diagnostic (it still conditions the
    invariants in the Coulomb).
    """
    psi = np.asarray(psi, dtype=complex)
    tide_columns = np.asarray(tide_columns, dtype=complex)
    if tide_columns.shape != (10,) + psi.shape:
        raise ValueError("tide_columns do not match the scalar grid")
    solution, relative_residual, condition_number = _fit_moment_columns(
        psi,
        tide_columns,
        amplitude_fraction=amplitude_fraction,
        newton_steps=newton_steps,
    )
    return TransverseMomentFit(
        electric=solution[:5],
        magnetic=solution[5:],
        newton_steps=newton_steps,
        relative_residual=relative_residual,
        condition_number=condition_number,
    )


@dataclass(frozen=True)
class JointMomentFit:
    """Joint solution of an arbitrary column subset (order as built)."""

    coefficients: np.ndarray
    newton_steps: int
    relative_residual: float
    condition_number: float


def fit_joint_moments(
    psi: np.ndarray,
    tide_columns: np.ndarray,
    *,
    amplitude_fraction: float = 1.0e-3,
    newton_steps: int = 2,
    fixed_trailing: np.ndarray | None = None,
) -> JointMomentFit:
    """Joint fit of any leading subset of the eps^3 column set (e.g. the
    24 quadrupole + octupole columns of the production target, deviation
    D6) -- same mechanism as :func:`fit_transverse_moments`.

    ``fixed_trailing`` holds the LAST k columns' coefficients fixed at
    the given values (e.g. dotted moments supplied by the time series,
    deviation D6/finding N3): their known response is moved to the data
    side and only the leading columns are solved for.  The returned
    coefficient vector always covers all columns (fitted then fixed);
    the condition number refers to the solved block.
    """
    psi = np.asarray(psi, dtype=complex)
    tide_columns = np.asarray(tide_columns, dtype=complex)
    if tide_columns.shape[1:] != psi.shape:
        raise ValueError("tide_columns do not match the scalar grid")
    design, data = _design_and_data(
        psi,
        tide_columns,
        amplitude_fraction=amplitude_fraction,
        newton_steps=newton_steps,
    )
    if fixed_trailing is not None:
        fixed_trailing = np.asarray(fixed_trailing, dtype=float)
        n_fixed = len(fixed_trailing)
        if not 0 < n_fixed < design.shape[1]:
            raise ValueError("fixed_trailing does not fit the column set")
        data = data - design[:, -n_fixed:] @ fixed_trailing
        design = design[:, :-n_fixed]
    solution, relative_residual, condition_number = _solve_design(design, data)
    if fixed_trailing is not None:
        solution = np.concatenate((solution, fixed_trailing))
    return JointMomentFit(
        coefficients=solution,
        newton_steps=newton_steps,
        relative_residual=relative_residual,
        condition_number=condition_number,
    )


def fit_third_order_moments(
    psi: np.ndarray,
    tide_columns: np.ndarray,
    *,
    amplitude_fraction: float = 1.0e-3,
    newton_steps: int = 2,
) -> ThirdOrderMomentFit:
    """Joint fit of the 34 eps^2 + eps^3 moments from the
    quasi-Kinnersley Psi4 (columns from
    :func:`third_order_tide_scalar_columns`); one linear least-squares
    problem, same mechanism as :func:`fit_transverse_moments`."""
    psi = np.asarray(psi, dtype=complex)
    tide_columns = np.asarray(tide_columns, dtype=complex)
    if tide_columns.shape != (34,) + psi.shape:
        raise ValueError("tide_columns do not match the scalar grid")
    solution, relative_residual, condition_number = _fit_moment_columns(
        psi,
        tide_columns,
        amplitude_fraction=amplitude_fraction,
        newton_steps=newton_steps,
    )
    return ThirdOrderMomentFit(
        electric=solution[:5],
        magnetic=solution[5:10],
        electric_octupole=solution[10:17],
        magnetic_octupole=solution[17:24],
        electric_dot=solution[24:29],
        magnetic_dot=solution[29:34],
        newton_steps=newton_steps,
        relative_residual=relative_residual,
        condition_number=condition_number,
    )


def _design_and_data(
    psi: np.ndarray,
    tide_columns: np.ndarray,
    *,
    amplitude_fraction: float,
    newton_steps: int,
):
    """eq:transverse-design for an arbitrary column set: difference
    eq:qk-psi4 through the whole measurement chain per column and stack
    the real system."""
    base = quasi_kinnersley_psi4(psi, newton_steps=newton_steps)
    scale = float(np.median(np.abs(psi[..., 2])))

    design_columns = []
    for column in tide_columns:
        amplitude = amplitude_fraction * scale / np.max(np.abs(column))
        plus = quasi_kinnersley_psi4(
            psi + amplitude * column, newton_steps=newton_steps
        )
        minus = quasi_kinnersley_psi4(
            psi - amplitude * column, newton_steps=newton_steps
        )
        design_columns.append((plus - minus) / (2.0 * amplitude))
    data_complex = base.reshape(-1)
    design_complex = np.stack(
        [column.reshape(-1) for column in design_columns], axis=-1
    )
    design = np.concatenate((design_complex.real, design_complex.imag))
    data = np.concatenate((data_complex.real, data_complex.imag))
    return design, data


def _solve_design(design: np.ndarray, data: np.ndarray):
    column_norm = np.linalg.norm(design, axis=0)
    if np.any(column_norm <= 0.0):
        raise RuntimeError("transverse moment design has an empty column")
    solution, *_ = np.linalg.lstsq(design / column_norm, data, rcond=None)
    solution = solution / column_norm
    singular_values = np.linalg.svd(design / column_norm, compute_uv=False)
    relative_residual = float(
        np.linalg.norm(design @ solution - data) / np.linalg.norm(data)
    )
    return (
        solution,
        relative_residual,
        float(singular_values[0] / singular_values[-1]),
    )


def _fit_moment_columns(
    psi: np.ndarray,
    tide_columns: np.ndarray,
    *,
    amplitude_fraction: float,
    newton_steps: int,
):
    design, data = _design_and_data(
        psi,
        tide_columns,
        amplitude_fraction=amplitude_fraction,
        newton_steps=newton_steps,
    )
    return _solve_design(design, data)


def psi0_target_ladder(
    coulomb: np.ndarray,
    b_total: np.ndarray,
    tide_scalars: np.ndarray,
) -> np.ndarray:
    """The boundary value assembled as the note writes it, eq:psi0-nr.

    ``coulomb`` is the measured Psi2^QK = -3J/I (kept WITH its tidal
    shift, eq:qk-orders), ``b_total`` = b^(0) + b^(2) the summed type-II
    parameter of eq:qk-rotation-expansion, and ``tide_scalars`` the
    NR-frame scalars of the fitted TRANSVERSE tide (its Psi0 slot carries
    the e^{2 eta} Psi0^model and b^4 e^{-2 eta} Psi4^QK terms exactly).

    Finding N2 (measured in the closure test): summing the rotation
    parameters, as eq:qk-rotation-expansion does, mis-composes the
    non-commuting type-I/type-II factors at O(a^(0) b^(2)) = O(eps^2 v),
    which leaves a tide-level error ~ 10 percent of the transverse
    content at v ~ 0.25.  Kept as the note-literal diagnostic; the
    production target is :func:`psi0_target` below.
    """
    return (
        kinnersley.psi0_leading(coulomb, b_total)
        + np.asarray(tide_scalars, dtype=complex)[..., 0]
    )


def psi0_target(
    coulomb: np.ndarray,
    fit: kinnersley.TypeDRotation,
    a_bar_correction: np.ndarray,
    b_correction: np.ndarray,
    tide_scalars: np.ndarray,
) -> np.ndarray:
    """The boundary value eq:psi0-nr with exact rotation composition.

    Kinematic slot: the Kinnersley scalars (0, 0, Psi2^K, 0, 0) with the
    measured Coulomb (kept WITH its tidal shift, eq:qk-orders) are
    carried to the NR frame by the exact inverse of the chain the data
    defined -- the inverse Newton corrections of Sec. 6.1 followed by the
    leading Sec.-3 rotations.  Composing the rotations instead of adding
    their parameters avoids the O(eps^2 v) error of finding N2; the
    unknown type-III factor drops because the input has only the
    weight-zero Psi2 slot.  The measured corrections carry the tide's
    vector channel, the measured Coulomb its longitudinal channel.

    Tide slot: the Psi0 component of ``tide_scalars`` -- the boosted
    TRANSVERSE model tide (tensor route, deviation D5), which carries the
    e^{2 eta} Psi0^model and b^4 e^{-2 eta} Psi4^QK terms of eq:psi0-nr
    exactly.
    """
    kinematic = kinnersley.kinnersley_scalars(coulomb)
    kinematic = scalars.type_i(kinematic, -np.asarray(a_bar_correction))
    kinematic = scalars.type_ii(kinematic, -np.asarray(b_correction))
    kinematic = scalars.type_ii(kinematic, fit.b)
    kinematic = scalars.type_i(kinematic, fit.a_bar)
    return kinematic[..., 0] + np.asarray(tide_scalars, dtype=complex)[..., 0]
