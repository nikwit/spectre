# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def _full_quadrupole(psi_monopole, psi_quadrupole, psi0, wt_radius):
    # the STF quadrupole plus the trace reconstructed from the boundary
    # monopole and the evolved constant coefficient
    return psi_quadrupole + np.eye(3) * (psi_monopole - psi0) / wt_radius**2


def dt_psi0(
    psi_monopole,
    psi_dipole,
    psi_quadrupole,
    dt_psi_dipole,
    inverse_spacetime_metric,
    trace_christoffel,
    evolved_vars,
    wt_radius,
    particle_velocity,
    particle_acceleration,
):
    return evolved_vars[1]


def dt2_psi0(
    psi_monopole,
    psi_dipole,
    psi_quadrupole,
    dt_psi_dipole,
    inverse_spacetime_metric,
    trace_christoffel,
    evolved_vars,
    wt_radius,
    particle_velocity,
    particle_acceleration,
):
    psi0 = evolved_vars[0]
    dt_psi0 = evolved_vars[1]
    v = particle_velocity
    a = particle_acceleration
    g = inverse_spacetime_metric
    psi_ij = _full_quadrupole(psi_monopole, psi_quadrupole, psi0, wt_radius)
    # time derivative of the dipole coefficient: the projection of the
    # time-derivative field corrected for the motion of the expansion center
    psi_dot = dt_psi_dipole + 2.0 * psi_ij @ v

    result = trace_christoffel[0] * dt_psi0
    result += np.dot(2.0 * g[0, 0] * v - 2.0 * g[0, 1:], psi_dot)
    result += np.dot(
        g[0, 0] * a - v * trace_christoffel[0] + trace_christoffel[1:],
        psi_dipole,
    )
    result -= np.einsum(
        "ij,ij",
        2.0 * g[0, 0] * np.outer(v, v)
        - 4.0 * np.outer(g[0, 1:], v)
        + 2.0 * g[1:, 1:],
        psi_ij,
    )
    return result / g[0, 0]
