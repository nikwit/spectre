# Distributed under the MIT License.
# See LICENSE.txt for details.

"""The quasi-Kinnersley frame (Sec. 6.1 of np_matching.tex).

The type-D solve of Sec. 3 imposed the Kinnersley conditions through the
Psi1 and Psi2 equations; the tide violates the remaining conditions at
O(eps^2).  Pulling the numerical scalars back with the leading rotations
exposes the longitudinal residuals eps^2 Psi1^(2), eps^2 Psi3^(2), and
the transversality conditions Psi1 = Psi3 = 0 are then restored by the
correction rotations of eq:qk-construction,

    b^(2) = Psi1^(2) / (3R),      a_bar^(2) = Psi3^(2) / (3R),

with R = M/r^3 = -Psi2^K from eq:typeD-IJ.  The closed-form step solves
the linearized conditions; reapplying it removes the reappearing
residuals order by order (the note states one step suffices at
O(eps^2) -- the returned history measures that claim).
"""

from __future__ import annotations

from dataclasses import dataclass

import kinnersley
import numpy as np
import scalars


@dataclass(frozen=True)
class QuasiKinnersleyFrame:
    """The pulled-back scalars, nulled to transversality.

    ``psi`` are the scalars in the quasi-Kinnersley frame (Psi1 = Psi3 =
    0 to the achieved tolerance); ``a_bar_correction`` and
    ``b_correction`` are the accumulated eps^2 rotations of
    eq:qk-rotation-expansion; ``residual_history[n]`` is the maximum of
    |Psi1|, |Psi3| over |Psi2| after n Newton steps (entry 0 is the raw
    pulled-back residual the note calls eps^2 Psi1^(2), Psi3^(2)).
    """

    psi: np.ndarray
    a_bar_correction: np.ndarray
    b_correction: np.ndarray
    residual_history: tuple[float, ...]


def _longitudinal_residual(psi: np.ndarray) -> float:
    return float(
        np.max(
            np.maximum(np.abs(psi[..., 1]), np.abs(psi[..., 3]))
            / np.abs(psi[..., 2])
        )
    )


def quasi_kinnersley(
    psi: np.ndarray,
    fit: kinnersley.TypeDRotation,
    coulomb: np.ndarray,
    *,
    steps: int = 1,
) -> QuasiKinnersleyFrame:
    """Construct the quasi-Kinnersley frame from measured data (Sec. 6.1).

    ``psi`` are the numerical NR-tetrad scalars, ``fit`` the Sec.-3 solve,
    and ``coulomb`` the measured Psi2^K.  Every quantity on the right of
    eq:qk-construction is measured data.  ``steps`` Newton passes are
    applied; the residual after each is recorded.
    """
    background = -np.asarray(coulomb, dtype=complex)  # R = -Psi2^K
    pulled = kinnersley.pull_back(psi, fit.a_bar, fit.b)
    a_bar_correction = np.zeros(pulled.shape[:-1], dtype=complex)
    b_correction = np.zeros(pulled.shape[:-1], dtype=complex)
    history = [_longitudinal_residual(pulled)]
    for _ in range(steps):
        b_step = pulled[..., 1] / (3.0 * background)
        a_step = pulled[..., 3] / (3.0 * background)
        pulled = scalars.type_i(scalars.type_ii(pulled, b_step), a_step)
        b_correction = b_correction + b_step
        a_bar_correction = a_bar_correction + a_step
        history.append(_longitudinal_residual(pulled))
    return QuasiKinnersleyFrame(
        psi=pulled,
        a_bar_correction=a_bar_correction,
        b_correction=b_correction,
        residual_history=tuple(history),
    )
