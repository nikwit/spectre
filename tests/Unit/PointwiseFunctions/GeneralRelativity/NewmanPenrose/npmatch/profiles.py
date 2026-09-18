# Distributed under the MIT License.
# See LICENSE.txt for details.

"""The tidal-tensor radial profiles (Appendix of np_matching.tex).

Table `tab:responses`: every r_hat-irreducible channel of every tidal
structure carries its own radial profile, exact in M/r at first order in
the moments.  ``quadrupole_tide_tensor`` assembles the eps^2 model tidal
tensor eq:tidal-tensor; ``octupole_tide_tensor``,
``induction_tide_tensor``, and ``nearzone_tide_tensor`` assemble the
three eps^3 structures of eq:app-Qflat and eq:app-nearzone.  Everything
is built from constant STF moment tensors, the measured radial
direction, and the measured radius; no metric evaluation and no
coordinate map appear -- the profiles are closed-form.

The near-zone sector carries the calibration constants (e1, b1) of the
time-differentiated moments (Table tab:radial-functions): shifting e1 by
delta is identical to redefining E -> E + delta M Edot, i.e. the shift
reproduces the quadrupole structure exactly (unit-tested).  One
convention must therefore be fixed and used everywhere; the default is
the horizon calibration.
"""

from __future__ import annotations

import numpy as np

# Calibrations (e1, b1) of the time-differentiated moments, note appendix.
HORIZON_CALIBRATION = (-92.0 / 15.0, -76.0 / 15.0)
POST_NEWTONIAN_CALIBRATION = (-52.0 / 15.0, -29.0 / 10.0)


def quadrupole_profiles(
    radius: np.ndarray, mass: float
) -> dict[str, np.ndarray]:
    """The six quadrupole radial profiles of Table tab:responses."""
    radius = np.asarray(radius, dtype=float)
    compactness = mass / radius
    f = 1.0 - 2.0 * compactness
    if np.any(f <= 0.0):
        raise ValueError("the profile radius lies inside the horizon")
    root_f = np.sqrt(f)
    return {
        "e_L": np.ones_like(f),
        "e_V": root_f * (1.0 + 2.0 * compactness),
        "e_T": f,
        "b_L": np.ones_like(f),
        "b_V": root_f,
        "b_T": f,
    }


def octupole_profiles(radius: np.ndarray, mass: float) -> dict[str, np.ndarray]:
    """The six octupole radial profiles of Table tab:responses."""
    radius = np.asarray(radius, dtype=float)
    compactness = mass / radius
    f = 1.0 - 2.0 * compactness
    if np.any(f <= 0.0):
        raise ValueError("the profile radius lies inside the horizon")
    root_f = np.sqrt(f)
    return {
        "e_L": 1.0 - 5.0 * compactness / 3.0,
        "e_V": root_f * (1.0 - compactness / 2.0 - compactness**2),
        "e_T": f * (1.0 - compactness),
        "b_L": 1.0 - 5.0 * compactness / 3.0,
        "b_V": root_f * (1.0 - 5.0 * compactness / 4.0),
        "b_T": f * (1.0 - compactness),
    }


def induction_profiles(
    radius: np.ndarray, mass: float
) -> dict[str, np.ndarray]:
    """The induction-sector radial profiles of Table tab:responses.

    The structure has no longitudinal part; the Edot row multiplies the
    imaginary (magnetic-type) pattern and the Bdot row the real one.
    """
    radius = np.asarray(radius, dtype=float)
    compactness = mass / radius
    f = 1.0 - 2.0 * compactness
    if np.any(f <= 0.0):
        raise ValueError("the profile radius lies inside the horizon")
    root_f = np.sqrt(f)
    return {
        "edot_V": (1.0 + 3.0 * compactness - 6.0 * compactness**3) / root_f,
        "edot_T": 1.0 - compactness / (4.0 * f) + 3.0 * compactness * f / 4.0,
        "bdot_V": (1.0 + 1.5 * compactness) / root_f,
        "bdot_T": 1.0 - compactness / (4.0 * f),
    }


def nearzone_profiles(
    radius: np.ndarray,
    mass: float,
    calibration: tuple[float, float] = HORIZON_CALIBRATION,
) -> dict[str, np.ndarray]:
    """The near-zone w-profiles of eq:app-nearzone-w (all proportional to M).

    ``calibration`` = (e1, b1) selects the time-differentiated-moment
    convention; the profiles carry it through e1 + 2 ln f and
    b1 + 2 ln f.
    """
    radius = np.asarray(radius, dtype=float)
    e1, b1 = calibration
    x = radius / mass
    f = 1.0 - 2.0 / x
    if np.any(f <= 0.0):
        raise ValueError("the profile radius lies inside the horizon")
    root_f = np.sqrt(f)
    log_e = e1 + 2.0 * np.log(f)
    log_b = b1 + 2.0 * np.log(f)
    return {
        "w_L": (
            mass
            * (
                log_e
                + 52.0 / 15.0
                + 4.0 / x
                + 4.0 / x**2
                + 16.0 / (3.0 * x**3)
                + 8.0 / x**4
            )
        ),
        "w_V": (
            mass
            / root_f
            * (
                (1.0 - 4.0 / x**2) * log_e
                + 52.0 / 15.0
                + 4.0 / x
                - 148.0 / (15.0 * x**2)
                - 32.0 / (3.0 * x**3)
                - 8.0 / x**4
            )
        ),
        "w_T": (
            mass
            * (
                f * log_e
                + (
                    52.0 / 15.0
                    - 148.0 / (15.0 * x)
                    + 28.0 / (15.0 * x**2)
                    + 16.0 / (3.0 * x**3)
                    + 8.0 / (3.0 * x**4)
                )
                / f
            )
        ),
        "wt_L": (
            mass
            * (
                log_b
                + 29.0 / 10.0
                + 4.0 / x
                + 4.0 / x**2
                + 16.0 / (3.0 * x**3)
                + 8.0 / x**4
            )
        ),
        "wt_V": (
            mass
            / root_f
            * (
                f * log_b
                + 29.0 / 10.0
                - 9.0 / (5.0 * x)
                - 4.0 / x**2
                - 8.0 / (3.0 * x**3)
                - 8.0 / (3.0 * x**4)
                + 16.0 / (3.0 * x**5)
            )
        ),
        "wt_T": (
            mass
            * (
                f * log_b
                + (
                    29.0 / 10.0
                    - 38.0 / (5.0 * x)
                    - 2.0 / (5.0 * x**2)
                    + 16.0 / (3.0 * x**3)
                    + 8.0 / (3.0 * x**4)
                )
                / f
            )
        ),
    }


def _irreducible_parts(tensor: np.ndarray, direction: np.ndarray):
    """The r_hat-irreducible pieces of eq:tide-irreducible.

    Returns (scalar, vector, transverse) with scalar = T:dd, vector =
    P T d, and transverse = P T P + (1/2) P scalar.
    """
    tensor = np.asarray(tensor, dtype=float)
    direction = np.asarray(direction, dtype=float)
    contracted = np.einsum("...ij,...j->...i", tensor, direction)
    scalar = np.einsum("...i,...i->...", direction, contracted)
    vector = contracted - scalar[..., None] * direction
    projector = np.eye(3) - np.einsum("...i,...j->...ij", direction, direction)
    transverse = (
        np.einsum("...ik,...kl,...lj->...ij", projector, tensor, projector)
        + 0.5 * scalar[..., None, None] * projector
    )
    return scalar, vector, transverse


def _assemble(
    scalar,
    vector,
    transverse,
    direction,
    profile_l,
    profile_v,
    profile_t,
    *,
    transverse_only,
):
    """eq:tidal-tensor for one parity sector."""
    if transverse_only:
        return profile_t[..., None, None] * transverse
    longitudinal_pattern = (
        np.einsum("...i,...j->...ij", direction, direction) - np.eye(3) / 3.0
    )
    vector_pattern = np.einsum(
        "...i,...j->...ij", direction, vector
    ) + np.einsum("...i,...j->...ij", vector, direction)
    return (
        1.5 * (profile_l * scalar)[..., None, None] * longitudinal_pattern
        + profile_v[..., None, None] * vector_pattern
        + profile_t[..., None, None] * transverse
    )


def quadrupole_tide_tensor(
    electric: np.ndarray,
    magnetic: np.ndarray,
    direction: np.ndarray,
    radius: np.ndarray,
    mass: float,
    *,
    transverse_only: bool = False,
) -> np.ndarray:
    """The model tidal tensor eq:tidal-tensor at first order in the moments.

    ``electric`` and ``magnetic`` are constant real STF tensors in the
    same orthonormal triad as ``direction`` (the measured radial
    direction); ``radius`` is the measured background radius.  With
    ``transverse_only`` the longitudinal and vector channels are omitted:
    that is the piece entering the Psi0/Psi4 slots of the aligned frame
    (eq:second-order-psi0) -- the L and V channels are carried by the
    measured Coulomb shift and the measured frame rotations instead, so
    including them in a boundary target would double count.

    In the far-field limit the full tensor reduces to E + iB
    (eq:app-Qflat, quadrupole sector).
    """
    profiles = quadrupole_profiles(radius, mass)
    e_scalar, e_vector, e_transverse = _irreducible_parts(electric, direction)
    b_scalar, b_vector, b_transverse = _irreducible_parts(magnetic, direction)
    electric_part = _assemble(
        e_scalar,
        e_vector,
        e_transverse,
        direction,
        profiles["e_L"],
        profiles["e_V"],
        profiles["e_T"],
        transverse_only=transverse_only,
    )
    magnetic_part = _assemble(
        b_scalar,
        b_vector,
        b_transverse,
        direction,
        profiles["b_L"],
        profiles["b_V"],
        profiles["b_T"],
        transverse_only=transverse_only,
    )
    return electric_part + 1j * magnetic_part


LEVI_CIVITA = np.zeros((3, 3, 3))
for _perm, _sign in (
    ((0, 1, 2), 1.0),
    ((1, 2, 0), 1.0),
    ((2, 0, 1), 1.0),
    ((0, 2, 1), -1.0),
    ((2, 1, 0), -1.0),
    ((1, 0, 2), -1.0),
):
    LEVI_CIVITA[_perm] = _sign


def octupole_tide_tensor(
    electric: np.ndarray,
    magnetic: np.ndarray,
    direction: np.ndarray,
    radius: np.ndarray,
    mass: float,
    *,
    transverse_only: bool = False,
) -> np.ndarray:
    """The eps^3 octupole structure of eq:app-Qflat with the octupole
    profiles of Table tab:responses.

    ``electric``/``magnetic`` are constant real rank-3 STF tensors; in
    the far field the structure is r (E_ijk + (4/3) i B_ijk) r_hat^k.
    """
    direction = np.asarray(direction, dtype=float)
    radius = np.asarray(radius, dtype=float)
    profiles = octupole_profiles(radius, mass)
    contracted_electric = np.einsum(
        "ijk,...k->...ij", np.asarray(electric, dtype=float), direction
    )
    contracted_magnetic = np.einsum(
        "ijk,...k->...ij", np.asarray(magnetic, dtype=float), direction
    )
    electric_part = _assemble(
        *_irreducible_parts(contracted_electric, direction),
        direction,
        profiles["e_L"],
        profiles["e_V"],
        profiles["e_T"],
        transverse_only=transverse_only,
    )
    magnetic_part = _assemble(
        *_irreducible_parts(contracted_magnetic, direction),
        direction,
        profiles["b_L"],
        profiles["b_V"],
        profiles["b_T"],
        transverse_only=transverse_only,
    )
    return radius[..., None, None] * (
        electric_part + 4.0j / 3.0 * magnetic_part
    )


def _induction_structure(
    moment: np.ndarray, direction: np.ndarray
) -> np.ndarray:
    """The symmetric structure eps_{kl(i} X_{j)}^k r_hat^l of eq:app-Qflat.

    Trace-free (X symmetric) and with no longitudinal part
    (eps_{kli} r_hat^l r_hat^i = 0)."""
    raw = np.einsum(
        "kli,...jk,...l->...ij",
        LEVI_CIVITA,
        np.asarray(moment, dtype=float),
        np.asarray(direction, dtype=float),
    )
    return 0.5 * (raw + raw.swapaxes(-2, -1))


def induction_tide_tensor(
    electric_dot: np.ndarray,
    magnetic_dot: np.ndarray,
    direction: np.ndarray,
    radius: np.ndarray,
    mass: float,
    *,
    transverse_only: bool = False,
) -> np.ndarray:
    """The flat-space induction sector of eq:app-Qflat with the
    induction profiles of Table tab:responses.

    In the far field the structure is
    (2r/3) eps_{kl(i} [i Edot_{j)}^k - Bdot_{j)}^k] r_hat^l -- a changing
    electric quadrupole sources a magnetic tidal pattern and vice versa.
    There is no longitudinal channel.
    """
    direction = np.asarray(direction, dtype=float)
    radius = np.asarray(radius, dtype=float)
    profiles = induction_profiles(radius, mass)
    zero_profile = np.zeros_like(profiles["edot_V"])
    electric_part = _assemble(
        *_irreducible_parts(
            _induction_structure(electric_dot, direction), direction
        ),
        direction,
        zero_profile,
        profiles["edot_V"],
        profiles["edot_T"],
        transverse_only=transverse_only,
    )
    magnetic_part = _assemble(
        *_irreducible_parts(
            _induction_structure(magnetic_dot, direction), direction
        ),
        direction,
        zero_profile,
        profiles["bdot_V"],
        profiles["bdot_T"],
        transverse_only=transverse_only,
    )
    return (
        (2.0 / 3.0)
        * radius[..., None, None]
        * (1j * electric_part - magnetic_part)
    )


def nearzone_tide_tensor(
    electric_dot: np.ndarray,
    magnetic_dot: np.ndarray,
    direction: np.ndarray,
    radius: np.ndarray,
    mass: float,
    *,
    calibration: tuple[float, float] = HORIZON_CALIBRATION,
    transverse_only: bool = False,
) -> np.ndarray:
    """The near-zone same-parity sector eq:app-nearzone.

    Proportional to M with no flat-space counterpart: the w-profiles of
    eq:app-nearzone-w multiply the L/V/T parts of the dotted moments
    with the SAME parity as the undotted quadrupole.  Carries the
    calibration constants; recalibrating is identical to the moment
    redefinition E -> E + delta_e1 M Edot (unit-tested).
    """
    profiles = nearzone_profiles(radius, mass, calibration)
    electric_part = _assemble(
        *_irreducible_parts(np.asarray(electric_dot, dtype=float), direction),
        direction,
        profiles["w_L"],
        profiles["w_V"],
        profiles["w_T"],
        transverse_only=transverse_only,
    )
    magnetic_part = _assemble(
        *_irreducible_parts(np.asarray(magnetic_dot, dtype=float), direction),
        direction,
        profiles["wt_L"],
        profiles["wt_V"],
        profiles["wt_T"],
        transverse_only=transverse_only,
    )
    return electric_part + 1j * magnetic_part
