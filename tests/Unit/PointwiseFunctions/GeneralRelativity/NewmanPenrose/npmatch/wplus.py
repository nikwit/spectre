# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Minimal second-order matching directly through ``Psi4``.

This module is deliberately narrower than :mod:`psi0`:

* the unknown is the five complex Cartesian components of ``E + i B``;
* the NR-to-quasi-Kinnersley map is stored and applied as an ordered product;
* the transverse Eq. (44) columns are projected through that fixed map, with
  no finite-difference derivative;
* the fitted quantity is the complex outgoing scalar ``Psi4`` itself;
* a fixed-point driver can run with numerical Psi0 removed from every input.

The components are ordered ``(H11, H22, H12, H13, H23)`` in the existing
Cholesky-orthonormal rest-frame convention, with ``H = E + i B`` and
``H33 = -H11 - H22``.  A local tangent-dyad spin multiplies the data and
every design column by the same phase, so the complex least-squares solution
is independent of that spin without constructing ``w+``.
"""

from __future__ import annotations

from dataclasses import dataclass

import kinnersley
import numpy as np
import psi0
import restframe
import scalars


@dataclass(frozen=True)
class OrderedQKMap:
    """The measured NR-to-QK scalar map as an ordered Lorentz product."""

    coulomb: np.ndarray
    leading: kinnersley.TypeDRotation
    correction_steps: tuple[tuple[np.ndarray, np.ndarray], ...]
    psi_qk: np.ndarray
    residual_history: tuple[float, ...]

    def apply(self, fields: np.ndarray) -> np.ndarray:
        """Carry one or more scalar five-vectors from NR to QK."""
        out = scalars.type_i(fields, -self.leading.a_bar)
        out = scalars.type_ii(out, -self.leading.b)
        for a_bar, b in self.correction_steps:
            out = scalars.type_ii(out, b)
            out = scalars.type_i(out, a_bar)
        return out

    def inverse(self, fields: np.ndarray) -> np.ndarray:
        """Carry one or more scalar five-vectors from QK back to NR."""
        out = np.asarray(fields, dtype=complex)
        for a_bar, b in reversed(self.correction_steps):
            out = scalars.type_i(out, -a_bar)
            out = scalars.type_ii(out, -b)
        out = scalars.type_ii(out, self.leading.b)
        out = scalars.type_i(out, self.leading.a_bar)
        return out


def _longitudinal_residual(fields: np.ndarray) -> float:
    scale = np.maximum(np.abs(fields[..., 2]), 1.0e-300)
    return float(
        np.max(
            np.maximum(np.abs(fields[..., 1]), np.abs(fields[..., 3])) / scale
        )
    )


def measure_qk_map(psi: np.ndarray, *, newton_steps: int = 1) -> OrderedQKMap:
    """Measure the QK frame while retaining every correction in order."""
    psi = np.asarray(psi, dtype=complex)
    invariant_i, invariant_j = scalars.invariants(psi)
    coulomb = kinnersley.coulomb_scalar(invariant_i, invariant_j)
    leading = kinnersley.solve_type_d_rotation(psi, coulomb)
    pulled = kinnersley.pull_back(psi, leading.a_bar, leading.b)
    background = -coulomb
    history = [_longitudinal_residual(pulled)]
    steps = []
    for _ in range(newton_steps):
        b_step = pulled[..., 1] / (3.0 * background)
        a_step = pulled[..., 3] / (3.0 * background)
        pulled = scalars.type_ii(pulled, b_step)
        pulled = scalars.type_i(pulled, a_step)
        steps.append((a_step, b_step))
        history.append(_longitudinal_residual(pulled))
    return OrderedQKMap(
        coulomb=coulomb,
        leading=leading,
        correction_steps=tuple(steps),
        psi_qk=pulled,
        residual_history=tuple(history),
    )


@dataclass(frozen=True)
class DirectPsi4Fit:
    """Five complex components of ``E + i B`` fitted from ``Psi4``."""

    components: np.ndarray
    relative_residual: float
    condition_number: float

    @property
    def electric(self) -> np.ndarray:
        return self.components.real

    @property
    def magnetic(self) -> np.ndarray:
        return self.components.imag

    def vector(self) -> np.ndarray:
        return np.concatenate((self.electric, self.magnetic))


def _solve_complex_design(
    data_psi4: np.ndarray,
    column_psi4: np.ndarray,
    *,
    point_weights: np.ndarray | None = None,
    point_mask: np.ndarray | None = None,
) -> DirectPsi4Fit:
    """Solve Eq. (44) for five complex Cartesian components.

    ``point_weights`` (optional) are quadrature weights: rows of the
    design and data are scaled by their square root, turning the plain
    collocation sum into the continuum L2(S^2) inner product.  Because
    the columns are band-limited, this weighted solve is identical to a
    modal least squares on any mode window containing the columns'
    support.  ``point_mask`` (optional, boolean) restricts the fit to a
    subset of the sphere points.
    """
    data = np.asarray(data_psi4, dtype=complex).reshape(-1)
    columns = np.asarray(column_psi4, dtype=complex)
    design = np.stack([column.reshape(-1) for column in columns], axis=-1)
    if point_mask is not None:
        mask = np.asarray(point_mask, dtype=bool).reshape(-1)
        data = data[mask]
        design = design[mask]
    if point_weights is not None:
        root = np.sqrt(np.asarray(point_weights, dtype=float).reshape(-1))
        if point_mask is not None:
            root = root[mask]
        data = data * root
        design = design * root[:, None]
    norms = np.linalg.norm(design, axis=0)
    if np.any(norms <= 0.0):
        raise RuntimeError("Psi4 design contains an empty component column")
    normalized = design / norms
    scaled_solution, *_ = np.linalg.lstsq(normalized, data, rcond=None)
    solution = scaled_solution / norms
    singular_values = np.linalg.svd(normalized, compute_uv=False)
    denominator = max(np.linalg.norm(data), 1.0e-300)
    residual = float(np.linalg.norm(design @ solution - data) / denominator)
    return DirectPsi4Fit(
        components=solution,
        relative_residual=residual,
        condition_number=float(singular_values[0] / singular_values[-1]),
    )


def fit_psi4(
    measured_qk_psi4: np.ndarray,
    projected_columns: np.ndarray,
    *,
    point_weights: np.ndarray | None = None,
    point_mask: np.ndarray | None = None,
) -> DirectPsi4Fit:
    """Fit the five complex components of ``E + i B`` from QK ``Psi4``.

    ``projected_columns`` has shape ``(5, ..., 5)`` and has already been
    carried through a fixed :class:`OrderedQKMap`.  Under a pointwise dyad
    spin the data and all five Psi4 columns acquire the same phase, leaving
    the complex least-squares problem unchanged.
    """
    projected_columns = np.asarray(projected_columns, dtype=complex)
    if projected_columns.shape[0] != 5 or projected_columns.shape[-1] != 5:
        raise ValueError("expected five projected scalar columns")
    return _solve_complex_design(
        measured_qk_psi4,
        projected_columns[..., 4],
        point_weights=point_weights,
        point_mask=point_mask,
    )


@dataclass(frozen=True)
class SecondOrderEvaluation:
    """One evaluation of the minimal second-order construction."""

    frame: OrderedQKMap
    fits: tuple[DirectPsi4Fit, ...]
    psi0_target: np.ndarray
    measured_radius: np.ndarray
    coordinate_velocities: np.ndarray
    rapidity: np.ndarray
    rapidity_fit: np.ndarray


@dataclass(frozen=True)
class FrameRegistration:
    """Shared leading frame for the second- and third-order models."""

    frame: OrderedQKMap
    member: restframe.TangentBoostMember
    measured_radius: np.ndarray
    radial_cholesky: np.ndarray
    transverse_cholesky: np.ndarray
    rapidity: np.ndarray
    velocity_fits: tuple[restframe.VelocityFit, ...]

    def column_arguments(self, slice_data) -> tuple:
        return (
            self.radial_cholesky,
            self.transverse_cholesky,
            self.rapidity,
            self.measured_radius,
            slice_data.adapted_rotation,
            slice_data.mass,
        )


def register_frame(
    slice_data,
    psi_for_frame: np.ndarray,
    *,
    newton_steps: int = 1,
    boost: str = "fit",
    dt_kretschmann: np.ndarray | None = None,
) -> FrameRegistration:
    """Measure the frame independently of the tidal response order.

    ``boost`` selects the type-III rapidity: ``"fit"`` is the sphere-wide
    velocity fit (eq:velocity-matching); ``"gradient"`` is the invariant
    curvature-gradient boost of note Sec. 3.2
    (eq:invariant-boost-condition), which is pointwise and consumes no
    trajectory or rigid-frame assumption.  The gradient boost needs
    ``dt_kretschmann``, dt of K = 16 Re I at fixed inertial coordinates
    on the ladder (shape ``(n_shell, n_point)``), from the on-shell
    single-slice evaluation; the spatial gradient comes from the slice
    cache itself (16 Re grad_invariant_i).
    """
    psi_for_frame = np.asarray(psi_for_frame, dtype=complex)
    frame = measure_qk_map(psi_for_frame, newton_steps=newton_steps)
    member = restframe.tangent_boost_member(
        frame.leading.a_bar, frame.leading.b
    )
    velocity_fits = tuple(
        restframe.fit_coordinate_velocity(
            member.radial_direction[shell],
            member.transverse_velocity[shell],
            slice_data.adapted_rotation[shell],
            slice_data.gamma[shell],
            slice_data.lapse[shell],
            slice_data.shift[shell],
        )
        for shell in range(len(slice_data.radii))
    )
    measured_radius = np.power(-slice_data.mass / frame.coulomb, 1.0 / 3.0).real
    radial_cholesky = np.einsum(
        "...ji,...j->...i",
        slice_data.adapted_rotation,
        member.radial_direction,
    )
    transverse_cholesky = np.einsum(
        "...ji,...j->...i",
        slice_data.adapted_rotation,
        member.transverse_velocity,
    )
    rapidity_fit = np.stack([fit.rapidity for fit in velocity_fits])
    if boost == "fit":
        rapidity = rapidity_fit
    elif boost == "gradient":
        if dt_kretschmann is None:
            raise ValueError("the gradient boost needs dt_kretschmann")
        rapidity = restframe.invariant_rapidity(
            member,
            slice_data.adapted_rotation,
            slice_data.gamma,
            slice_data.lapse,
            slice_data.shift,
            16.0 * slice_data.grad_invariant_i.real,
            dt_kretschmann,
        ).rapidity
    else:
        raise ValueError(f"unknown boost mode: {boost}")
    return FrameRegistration(
        frame,
        member,
        measured_radius,
        radial_cholesky,
        transverse_cholesky,
        rapidity,
        velocity_fits,
    )


def evaluate_second_order(
    slice_data,
    psi_for_frame: np.ndarray,
    *,
    newton_steps: int = 1,
    boost: str = "fit",
    dt_kretschmann: np.ndarray | None = None,
    fit_point_weights: np.ndarray | None = None,
    fit_point_mask: np.ndarray | None = None,
) -> SecondOrderEvaluation:
    """Measure the frame, fit QK ``Psi4``, and construct an NR Psi0 target."""
    registered = register_frame(
        slice_data,
        psi_for_frame,
        newton_steps=newton_steps,
        boost=boost,
        dt_kretschmann=dt_kretschmann,
    )
    frame = registered.frame
    direct_columns = psi0.direct_tide_scalar_columns(
        *registered.column_arguments(slice_data), transverse_only=True
    )
    projected_columns = frame.apply(direct_columns)
    fits = tuple(
        fit_psi4(
            frame.psi_qk[shell, ..., 4],
            projected_columns[:, shell],
            point_weights=fit_point_weights,
            point_mask=fit_point_mask,
        )
        for shell in range(len(slice_data.radii))
    )
    coefficients = np.stack([fit.components for fit in fits])
    fitted_tide_nr = np.einsum("sa,as...->s...", coefficients, direct_columns)
    qk_coulomb = kinnersley.kinnersley_scalars(frame.coulomb)
    kinematic_nr = frame.inverse(qk_coulomb)
    target = kinematic_nr[..., 0] + fitted_tide_nr[..., 0]
    return SecondOrderEvaluation(
        frame=frame,
        fits=fits,
        psi0_target=target,
        measured_radius=registered.measured_radius,
        coordinate_velocities=np.stack(
            [fit.coordinate_velocity for fit in registered.velocity_fits]
        ),
        rapidity=registered.rapidity,
        rapidity_fit=np.stack(
            [fit.rapidity for fit in registered.velocity_fits]
        ),
    )


@dataclass(frozen=True)
class BlindFixedPoint:
    """A converged evaluation that never consumed numerical Psi0."""

    evaluation: SecondOrderEvaluation
    psi0: np.ndarray
    relative_change_history: tuple[float, ...]
    shell_change_history: tuple[np.ndarray, ...]
    converged: bool


def blind_fixed_point(
    slice_data,
    *,
    tolerance: float = 1.0e-10,
    max_iterations: int = 12,
    damping: float = 1.0,
    newton_steps: int = 1,
    boost: str = "fit",
    dt_kretschmann: np.ndarray | None = None,
    fit_point_weights: np.ndarray | None = None,
    fit_point_mask: np.ndarray | None = None,
) -> BlindFixedPoint:
    """Close the boundary construction with numerical Psi0 removed.

    The known slots Psi1..Psi4 are copied from the numerical slice only for
    this offline study.  Slot zero starts at zero and is replaced solely by
    the previous model target.  Convergence is measured with the physical
    ``w-`` Frobenius norm, equivalently the complex Psi0 norm.

    With ``boost="gradient"`` the only numerical-Psi0 touchpoint of the
    construction is the invariant field K = 16 Re I inside the boost
    (through its Psi0*Psi4 term, an O(eps^4) relative contribution).
    """
    if not 0.0 < damping <= 1.0:
        raise ValueError("damping must lie in (0, 1]")
    working = np.asarray(slice_data.psi, dtype=complex).copy()
    working[..., 0] = 0.0
    history = []
    shell_history = []
    converged = False
    for _ in range(max_iterations):
        evaluation = evaluate_second_order(
            slice_data,
            working,
            newton_steps=newton_steps,
            boost=boost,
            dt_kretschmann=dt_kretschmann,
            fit_point_weights=fit_point_weights,
            fit_point_mask=fit_point_mask,
        )
        proposed = evaluation.psi0_target
        old = working[..., 0]
        shell_change = np.linalg.norm(proposed - old, axis=1) / np.maximum(
            np.linalg.norm(proposed, axis=1), 1.0e-300
        )
        change = float(np.max(shell_change))
        history.append(change)
        shell_history.append(shell_change)
        working[..., 0] = (1.0 - damping) * old + damping * proposed
        if change < tolerance:
            converged = True
            break
    final = evaluate_second_order(
        slice_data,
        working,
        newton_steps=newton_steps,
        boost=boost,
        dt_kretschmann=dt_kretschmann,
        fit_point_weights=fit_point_weights,
        fit_point_mask=fit_point_mask,
    )
    return BlindFixedPoint(
        evaluation=final,
        psi0=final.psi0_target,
        relative_change_history=tuple(history),
        shell_change_history=tuple(shell_history),
        converged=converged,
    )
