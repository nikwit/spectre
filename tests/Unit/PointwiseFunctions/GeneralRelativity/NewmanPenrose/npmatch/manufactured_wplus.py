# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Closure tests for the minimal direct-Psi4 matching study."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import profiles  # noqa: E402
import psi0  # noqa: E402
import scalars  # noqa: E402
import tides  # noqa: E402
import wplus  # noqa: E402

MASS = 1.0 / 9.0


def _random_stf(rng, scale):
    raw = scale * rng.normal(size=(3, 3))
    symmetric = 0.5 * (raw + raw.T)
    return symmetric - np.trace(symmetric) * np.eye(3) / 3.0


def _manufactured_slice(seed=71, n=240):
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=(n, 3))
    direction /= np.linalg.norm(direction, axis=-1)[:, None]
    gamma = np.broadcast_to(np.eye(3), (n, 3, 3)).copy()
    rotation = scalars.adapted_tetrad_rotation(gamma, direction)
    radius = 0.5
    background = MASS / radius**3
    electric = _random_stf(rng, 1.0e-3 * background)
    magnetic = _random_stf(rng, 0.5e-3 * background)
    velocity = (
        0.22 * np.asarray([0.5, -0.8, 0.3]) / np.linalg.norm([0.5, -0.8, 0.3])
    )
    q_rest = background * (
        np.eye(3) - 3.0 * direction[:, :, None] * direction[:, None, :]
    ) + profiles.quadrupole_tide_tensor(
        electric, magnetic, direction, np.full(n, radius), MASS
    )
    radial_speed = direction @ velocity
    transverse_velocity = velocity[None, :] - radial_speed[:, None] * direction
    transverse_lorentz = 1.0 / np.sqrt(
        1.0 - np.einsum("ni,ni->n", transverse_velocity, transverse_velocity)
    )
    rapidity = np.arctanh(transverse_lorentz * radial_speed)
    boost = psi0.rest_to_slice_map(direction, transverse_velocity, rapidity)
    q_slice = np.einsum("nik,nkl,njl->nij", boost, q_rest, boost)
    psi = scalars.psi_from_q(scalars.rotate_symmetric(q_slice, rotation))
    data = SimpleNamespace(
        mass=MASS,
        radii=np.asarray([radius]),
        gamma=gamma[None],
        lapse=np.ones((1, n)),
        shift=np.zeros((1, n, 3)),
        adapted_rotation=rotation[None],
        psi=psi[None],
    )
    truth_complex = tides.stf_components(electric) + 1j * tides.stf_components(
        magnetic
    )
    truth = np.concatenate((truth_complex.real, truth_complex.imag))
    return data, truth


def _spin_slice(data, chi):
    rotation = data.adapted_rotation.copy()
    tangent_1 = rotation[..., 1, :].copy()
    tangent_2 = rotation[..., 2, :].copy()
    rotation[..., 1, :] = (
        np.cos(chi)[..., None] * tangent_1 + np.sin(chi)[..., None] * tangent_2
    )
    rotation[..., 2, :] = (
        -np.sin(chi)[..., None] * tangent_1 + np.cos(chi)[..., None] * tangent_2
    )
    return SimpleNamespace(
        mass=data.mass,
        radii=data.radii,
        gamma=data.gamma,
        lapse=data.lapse,
        shift=data.shift,
        adapted_rotation=rotation,
        psi=scalars.type_iii(data.psi, chi=-chi),
    )


def test_ordered_qk_map_roundtrip():
    rng = np.random.default_rng(72)
    arbitrary = rng.normal(size=(32, 5)) + 1j * rng.normal(size=(32, 5))
    coulomb = -(0.5 + rng.random(32))
    a_bar = 0.1 * (rng.normal(size=32) + 1j * rng.normal(size=32))
    b = 0.1 * (rng.normal(size=32) + 1j * rng.normal(size=32))
    base = np.zeros((32, 5), dtype=complex)
    base[..., 2] = coulomb
    base = scalars.type_i(scalars.type_ii(base, b), a_bar)
    frame = wplus.measure_qk_map(base, newton_steps=2)
    assert np.allclose(
        frame.inverse(frame.apply(arbitrary)), arbitrary, rtol=1e-10, atol=1e-11
    )


def test_direct_psi4_fit_recovers_components_and_ignores_dyad_phase():
    """Eq. (44) is one complex least-squares problem in direct components."""
    rng = np.random.default_rng(74)
    components = rng.normal(size=5) + 1j * rng.normal(size=5)
    design = rng.normal(size=(180, 5)) + 1j * rng.normal(size=(180, 5))
    projected_columns = np.zeros((5, 180, 5), dtype=complex)
    projected_columns[..., 4] = design.T
    data_psi4 = design @ components

    fit = wplus.fit_psi4(data_psi4, projected_columns)
    assert np.allclose(fit.components, components, rtol=1e-12, atol=1e-12)

    phase = np.exp(2j * rng.uniform(-np.pi, np.pi, size=180))
    spun_columns = projected_columns * phase[None, :, None]
    spun_fit = wplus.fit_psi4(data_psi4 * phase, spun_columns)
    assert np.allclose(spun_fit.components, components, rtol=1e-12, atol=1e-12)


def test_blind_psi4_closure_and_spin_independence():
    data, truth = _manufactured_slice()
    blind = wplus.blind_fixed_point(data, tolerance=1.0e-11, max_iterations=12)
    assert blind.converged
    recovered = blind.evaluation.fits[0].vector()
    assert np.linalg.norm(recovered - truth) / np.linalg.norm(truth) < 0.05
    target_error = np.linalg.norm(
        blind.psi0 - data.psi[..., 0]
    ) / np.linalg.norm(data.psi[..., 0])
    assert target_error < 5.0e-4

    rng = np.random.default_rng(73)
    corrupted_psi = data.psi.copy()
    corrupted_psi[..., 0] = 100.0 * (
        rng.normal(size=corrupted_psi.shape[:-1])
        + 1j * rng.normal(size=corrupted_psi.shape[:-1])
    )
    corrupted_data = SimpleNamespace(**{**vars(data), "psi": corrupted_psi})
    corrupted = wplus.blind_fixed_point(
        corrupted_data, tolerance=1.0e-11, max_iterations=12
    )
    assert np.allclose(corrupted.psi0, blind.psi0, rtol=0.0, atol=0.0)
    assert np.allclose(
        corrupted.evaluation.fits[0].vector(), recovered, rtol=0.0, atol=0.0
    )

    chi = rng.uniform(-np.pi, np.pi, size=data.psi.shape[:-1])
    spun_data = _spin_slice(data, chi)
    spun = wplus.blind_fixed_point(
        spun_data, tolerance=1.0e-11, max_iterations=12
    )
    assert spun.converged
    spun_recovered = spun.evaluation.fits[0].vector()
    assert np.allclose(spun_recovered, recovered, rtol=2e-8, atol=2e-11)
    field = scalars.w_minus(blind.psi0, data.adapted_rotation)
    spun_field = scalars.w_minus(spun.psi0, spun_data.adapted_rotation)
    assert np.allclose(spun_field, field, rtol=2e-8, atol=2e-11)


ALL_TESTS = (
    test_ordered_qk_map_roundtrip,
    test_direct_psi4_fit_recovers_components_and_ignores_dyad_phase,
    test_blind_psi4_closure_and_spin_independence,
)


if __name__ == "__main__":
    for test in ALL_TESTS:
        test()
        print(f"{test.__name__}: ok")
    print(f"{len(ALL_TESTS)} tests passed")
