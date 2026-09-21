# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Adapter between pypp and the vendored `npmatch` reference implementation.

The reference lives unchanged in the `npmatch/` subdirectory (copied from the
worldtube NP-matching study). pypp calls the pointwise functions below once
per grid point with plain NumPy arrays of the tensor's index shape, complex
where the C++ side is complex; the sphere-wide functions receive lists of
one-dimensional arrays, one per tensor component, and reassemble them.
"""

import os
import sys
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "npmatch"))

import coulomb  # noqa: E402
import kinnersley  # noqa: E402
import manufactured_wplus  # noqa: E402
import profiles  # noqa: E402
import psi0  # noqa: E402
import restframe  # noqa: E402
import scalars  # noqa: E402
import wplus  # noqa: E402

# --- Tetrad ---------------------------------------------------------------


def cholesky_factor(metric):
    return np.linalg.cholesky(metric)


def inverse_lower_triangular(lower):
    return np.linalg.inv(lower)


def adapted_triad(metric, directions):
    return scalars.adapted_tetrad_rotation(metric, directions)


def rotate_symmetric(tensor, rotation):
    return scalars.rotate_symmetric(tensor, rotation)


def weyl_scalars_from_tidal_tensor(q):
    return scalars.psi_from_q(q)


def tidal_tensor_from_weyl_scalars(psi):
    return scalars.q_from_psi(psi)


def weyl_scalars_from_electric_magnetic(electric, magnetic, gamma, directions):
    # The derivation of data.load_slice: Cholesky triad, hygiene projection,
    # adapted triad, scalars
    cholesky_inverse = np.linalg.inv(np.linalg.cholesky(gamma))
    q_raw = cholesky_inverse @ (electric + 1j * magnetic) @ cholesky_inverse.T
    q_symmetric = 0.5 * (q_raw + q_raw.T)
    q = q_symmetric - np.trace(q_symmetric) * np.eye(3) / 3.0
    rotation = scalars.adapted_tetrad_rotation(gamma, directions)
    return scalars.psi_from_q(scalars.rotate_symmetric(q, rotation))


def incoming_weyl_field(psi0_value, rotation):
    return scalars.w_minus(psi0_value, rotation)


def orthonormal_to_coordinate_covariant(orthonormal, cholesky):
    return np.einsum("ia,jb,ab->ij", cholesky, cholesky, orthonormal)


# --- Null rotations and invariants -----------------------------------------


def type_i(psi, a_bar):
    return scalars.type_i(psi, a_bar)


def type_ii(psi, b):
    return scalars.type_ii(psi, b)


def type_iii(psi, eta, chi):
    return scalars.type_iii(psi, eta, chi)


def invariant_i(psi):
    return scalars.invariants(psi)[0]


def invariant_j(psi):
    return scalars.invariants(psi)[1]


# --- Type D ----------------------------------------------------------------


def coulomb_scalar(invariant_i_value, invariant_j_value):
    return kinnersley.coulomb_scalar(invariant_i_value, invariant_j_value)


def background_radius(invariant_i_value, invariant_j_value, mass):
    return kinnersley.background_radius(
        invariant_i_value, invariant_j_value, mass
    )


def kinnersley_scalars(coulomb):
    return kinnersley.kinnersley_scalars(coulomb)


def type_d_scalars(coulomb, a_bar, b):
    return kinnersley.type_d_scalars(coulomb, a_bar, b)


def solve_type_d_a_bar(psi, coulomb):
    return kinnersley.solve_type_d_rotation(psi, coulomb).a_bar


def solve_type_d_b(psi, coulomb):
    return kinnersley.solve_type_d_rotation(psi, coulomb).b


def solve_type_d_x(psi, coulomb):
    return kinnersley.solve_type_d_rotation(psi, coulomb).x


def solve_type_d_predicted_psi(psi, coulomb):
    return kinnersley.solve_type_d_rotation(psi, coulomb).predicted_psi


def pull_back(psi, a_bar, b):
    return kinnersley.pull_back(psi, a_bar, b)


def psi0_leading(coulomb, b):
    return kinnersley.psi0_leading(coulomb, b)


# --- Rest frame --------------------------------------------------------------


def principal_null_outgoing(a_bar, b):
    return restframe.principal_null_pair(a_bar, b)[0]


def principal_null_incoming(a_bar, b):
    return restframe.principal_null_pair(a_bar, b)[1]


def tangent_radial_direction(a_bar, b):
    return restframe.tangent_boost_member(a_bar, b).radial_direction


def tangent_transverse_velocity(a_bar, b):
    return restframe.tangent_boost_member(a_bar, b).transverse_velocity


def tangent_lorentz_factor(a_bar, b):
    return restframe.tangent_boost_member(a_bar, b).lorentz_factor


def invariant_rapidity(
    radial_direction,
    transverse_velocity,
    lorentz_factor,
    rotation,
    gamma,
    lapse,
    shift,
    spatial_gradient,
    time_derivative,
):
    member = restframe.TangentBoostMember(
        radial_direction=np.asarray(radial_direction),
        transverse_velocity=np.asarray(transverse_velocity),
        lorentz_factor=np.asarray(lorentz_factor),
        max_imaginary_part=0.0,
        max_null_product_error=0.0,
        max_radial_norm_error=0.0,
        max_transverse_radial_overlap=0.0,
    )
    return restframe.invariant_rapidity(
        member, rotation, gamma, lapse, shift, spatial_gradient, time_derivative
    ).rapidity


# --- Tidal response ----------------------------------------------------------


def quadrupole_profile(radius, mass, which):
    return profiles.quadrupole_profiles(radius, mass)[which]


def quadrupole_tide_tensor(
    electric, magnetic, direction, radius, mass, transverse_only
):
    return profiles.quadrupole_tide_tensor(
        electric,
        magnetic,
        direction,
        radius,
        mass,
        transverse_only=bool(transverse_only),
    )


def self_dual_boost(axis, rapidity):
    return psi0.self_dual_boost(axis, rapidity)


def rest_to_slice_map(radial_direction, transverse_velocity, rapidity):
    return psi0.rest_to_slice_map(
        radial_direction, transverse_velocity, rapidity
    )


def direct_tide_scalar_column(
    radial_direction,
    transverse_velocity,
    rapidity,
    radius,
    rotation,
    mass,
    component,
):
    return psi0.direct_tide_scalar_columns(
        radial_direction,
        transverse_velocity,
        rapidity,
        radius,
        rotation,
        mass,
        transverse_only=True,
    )[int(component)]


# --- Sphere-wide functions -----------------------------------------------------
#
# Real face data travel as a list of 25 one-dimensional arrays in the order
# gamma (9, row major), rotation (9, row major), lapse (1), shift (3),
# directions (3); the Weyl scalars as a list of five complex arrays.


def _unpack_reals(reals):
    reals = [np.asarray(r, dtype=float) for r in reals]
    n = reals[0].size
    gamma = np.stack(reals[0:9], axis=-1).reshape(n, 3, 3)
    rotation = np.stack(reals[9:18], axis=-1).reshape(n, 3, 3)
    lapse = reals[18]
    shift = np.stack(reals[19:22], axis=-1)
    directions = np.stack(reals[22:25], axis=-1)
    return gamma, rotation, lapse, shift, directions


def _pack_reals(gamma, rotation, lapse, shift, directions):
    n = lapse.size
    return (
        [gamma.reshape(n, 9)[:, k] for k in range(9)]
        + [rotation.reshape(n, 9)[:, k] for k in range(9)]
        + [lapse]
        + [shift[:, k] for k in range(3)]
        + [directions[:, k] for k in range(3)]
    )


def _unpack_psi(psi_list):
    return np.stack([np.asarray(p, dtype=complex) for p in psi_list], axis=-1)


def _slice_data(psi_list, reals, mass):
    gamma, rotation, lapse, shift, _ = _unpack_reals(reals)
    return SimpleNamespace(
        mass=float(mass),
        radii=np.asarray([1.0]),
        gamma=gamma[None],
        lapse=lapse[None],
        shift=shift[None],
        adapted_rotation=rotation[None],
        psi=_unpack_psi(psi_list)[None],
    )


def fit_psi4(data, column_list, weights):
    columns = np.stack([np.asarray(c, dtype=complex) for c in column_list])
    fit = wplus._solve_complex_design(
        np.asarray(data, dtype=complex),
        columns,
        point_weights=None if len(weights) == 0 else np.asarray(weights),
    )
    return list(fit.components.real) + list(fit.components.imag)


def _evaluate_second_order(psi_list, reals, mass, rapidity):
    """wplus.evaluate_second_order without Newton corrections and with the
    type-III rapidity supplied instead of fitted: the leading type-D map,
    the tangent member, the direct columns pulled back through the same map,
    the complex fit of the pulled-back Psi4, and the exact-composition
    target."""
    gamma, rotation, lapse, shift, _ = _unpack_reals(reals)
    psi = _unpack_psi(psi_list)
    frame = wplus.measure_qk_map(psi, newton_steps=0)
    member = restframe.tangent_boost_member(
        frame.leading.a_bar, frame.leading.b
    )
    measured_radius = np.power(-mass / frame.coulomb, 1.0 / 3.0).real
    radial_cholesky = np.einsum(
        "...ji,...j->...i", rotation, member.radial_direction
    )
    transverse_cholesky = np.einsum(
        "...ji,...j->...i", rotation, member.transverse_velocity
    )
    direct_columns = psi0.direct_tide_scalar_columns(
        radial_cholesky,
        transverse_cholesky,
        np.asarray(rapidity),
        measured_radius,
        rotation,
        mass,
        transverse_only=True,
    )
    projected = frame.apply(direct_columns)
    fit = wplus.fit_psi4(frame.psi_qk[..., 4], projected)
    fitted_tide = np.einsum("a,a...->...", fit.components, direct_columns)
    kinematic = frame.inverse(kinnersley.kinnersley_scalars(frame.coulomb))
    return SimpleNamespace(
        components=fit.components,
        psi0_target=kinematic[..., 0] + fitted_tide[..., 0],
        measured_radius=measured_radius,
    )


def evaluate_second_order_psi0_target(psi_list, reals, mass, rapidity):
    return _evaluate_second_order(psi_list, reals, mass, rapidity).psi0_target


def evaluate_second_order_components(psi_list, reals, mass, rapidity):
    components = _evaluate_second_order(
        psi_list, reals, mass, rapidity
    ).components
    return list(components.real) + list(components.imag)


def register_frame_radius(psi_list, reals, mass):
    return _evaluate_second_order(
        psi_list, reals, mass, np.zeros(len(reals[0]))
    ).measured_radius


# --- Manufactured slice of the minimal study ------------------------------------


def manufactured_psi(seed, n):
    data, _ = manufactured_wplus._manufactured_slice(seed=int(seed), n=int(n))
    return [data.psi[0][:, a] for a in range(5)]


def manufactured_reals(seed, n):
    data, _ = manufactured_wplus._manufactured_slice(seed=int(seed), n=int(n))
    # The manufactured slice is flat with unit lapse and zero shift; the
    # sphere directions are the first row of the adapted triad because the
    # metric is flat.
    directions = data.adapted_rotation[0][:, 0, :]
    return _pack_reals(
        data.gamma[0],
        data.adapted_rotation[0],
        data.lapse[0],
        data.shift[0],
        directions,
    )


def manufactured_rapidity(seed, n):
    """The exact type-III rapidity of the manufactured boost at each point,
    recomputed from the construction in manufactured_wplus."""
    data, _ = manufactured_wplus._manufactured_slice(seed=int(seed), n=int(n))
    direction = data.adapted_rotation[0][:, 0, :]
    velocity = (
        0.22 * np.asarray([0.5, -0.8, 0.3]) / np.linalg.norm([0.5, -0.8, 0.3])
    )
    radial_speed = direction @ velocity
    transverse_velocity = velocity[None, :] - radial_speed[:, None] * direction
    transverse_lorentz = 1.0 / np.sqrt(
        1.0 - np.einsum("ni,ni->n", transverse_velocity, transverse_velocity)
    )
    return np.arctanh(transverse_lorentz * radial_speed)


def manufactured_truth(seed, n):
    _, truth = manufactured_wplus._manufactured_slice(seed=int(seed), n=int(n))
    return [float(t) for t in truth]


def manufactured_mass():
    return float(manufactured_wplus.MASS)


# --- Coulomb-channel decode ----------------------------------------------------


def _coulomb_registration(psi_list, reals, mass):
    """The frame registration of gr::np::register_frame: leading type-D map,
    tangent member, measured radius and the member's vectors in the
    Cholesky frame."""
    _, rotation, _, _, _ = _unpack_reals(reals)
    psi = _unpack_psi(psi_list)
    frame = wplus.measure_qk_map(psi, newton_steps=0)
    member = restframe.tangent_boost_member(
        frame.leading.a_bar, frame.leading.b
    )
    measured_radius = np.power(-mass / frame.coulomb, 1.0 / 3.0).real
    radial_cholesky = np.einsum(
        "...ji,...j->...i", rotation, member.radial_direction
    )
    transverse_cholesky = np.einsum(
        "...ji,...j->...i", rotation, member.transverse_velocity
    )
    return SimpleNamespace(
        coulomb=frame.coulomb,
        member=member,
        measured_radius=measured_radius,
        radial_cholesky=radial_cholesky,
        transverse_cholesky=transverse_cholesky,
    )


def coulomb_background_radial_derivative(radius, mass):
    return coulomb.background_radial_derivative(np.asarray(radius), mass)


def coulomb_normal_derivative_factor(psi_list, reals, mass, rapidity):
    reg = _coulomb_registration(psi_list, reals, mass)
    return coulomb.normal_derivative_factor(
        reg.member.radial_direction,
        reg.member.transverse_velocity,
        reg.member.lorentz_factor,
        np.asarray(rapidity),
    )


def coulomb_radius_from_normal_derivative(d_s_coulomb, factor, initial, mass):
    radius, _, _ = coulomb.radius_from_normal_derivative(
        np.asarray(d_s_coulomb), np.asarray(factor), np.asarray(initial), mass
    )
    return radius


def coulomb_radius_solve_valid(d_s_coulomb, factor, initial, mass):
    _, valid, _ = coulomb.radius_from_normal_derivative(
        np.asarray(d_s_coulomb), np.asarray(factor), np.asarray(initial), mass
    )
    return bool(valid)


def coulomb_tide_column(psi_list, reals, mass, radius, index):
    reg = _coulomb_registration(psi_list, reals, mass)
    return coulomb.coulomb_tide_columns(
        reg.radial_cholesky, np.asarray(radius), mass
    )[int(index)]


def _coulomb_decode(psi_list, reals, mass, rapidity, d_s_coulomb, weights):
    reg = _coulomb_registration(psi_list, reals, mass)
    return coulomb.decode_tidal_moments_from_coulomb(
        reg.coulomb,
        np.asarray(d_s_coulomb),
        reg.member.radial_direction,
        reg.member.transverse_velocity,
        reg.member.lorentz_factor,
        np.asarray(rapidity),
        reg.measured_radius,
        reg.radial_cholesky,
        reg.transverse_cholesky,
        mass,
        None if len(weights) == 0 else np.asarray(weights),
    )


def coulomb_decode_components(psi_list, reals, mass, rapidity, d_s, weights):
    decode = _coulomb_decode(psi_list, reals, mass, rapidity, d_s, weights)
    return list(decode.components.real) + list(decode.components.imag)


def coulomb_decode_radius(psi_list, reals, mass, rapidity, d_s, weights):
    return _coulomb_decode(psi_list, reals, mass, rapidity, d_s, weights).radius


def coulomb_decode_residual(psi_list, reals, mass, rapidity, d_s, weights):
    return float(
        _coulomb_decode(
            psi_list, reals, mass, rapidity, d_s, weights
        ).relative_residual
    )
