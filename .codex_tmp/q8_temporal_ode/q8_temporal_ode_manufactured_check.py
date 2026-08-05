#!/usr/bin/env python3
"""Manufactured check for the q=8 temporal ODE analysis."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ANALYSIS = (
    Path(__file__).resolve().parent / "q8_temporal_ode_consistency.py"
)
spec = importlib.util.spec_from_file_location("q8_temporal_ode", ANALYSIS)
ode = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(ode)

times = np.arange(950.0, 1000.0, 5.0)
radius_index = 3
first_order_exact = []
second_order_exact = []
first_order_recovered = []
second_order_recovered = []

for time in times:
    shifted_time = time - 972.5
    first_order = np.zeros(13)
    first_order[0] = (
        0.08 - 4.0e-5 * shifted_time
        + 2.0e-8 * shifted_time**2
    )
    first_order[1:4] = (0.002, -0.001, 0.0003)
    first_order[4:7] = (-0.01, 0.015, 0.0)
    strain = np.array([
        [
            0.01 + 2.0e-4 * shifted_time
            + 1.0e-7 * shifted_time**2,
            0.003 - 8.0e-5 * shifted_time,
            0.0,
        ],
        [
            0.003 - 8.0e-5 * shifted_time,
            -0.006 + 1.0e-4 * shifted_time
            - 2.0e-7 * shifted_time**2,
            0.0,
        ],
        [0.0, 0.0, 0.002 + 3.0e-5 * shifted_time],
    ])
    first_order[7:13] = (
        strain[0, 0],
        strain[0, 1],
        strain[0, 2],
        strain[1, 1],
        strain[1, 2],
        strain[2, 2],
    )

    clock_derivative = -4.0e-5 + 4.0e-8 * shifted_time
    strain_derivative = np.array([
        [
            2.0e-4 + 2.0e-7 * shifted_time,
            -8.0e-5,
            0.0,
        ],
        [
            -8.0e-5,
            1.0e-4 - 4.0e-7 * shifted_time,
            0.0,
        ],
        [0.0, 0.0, 3.0e-5],
    ])
    strain_rate = (
        ode.SMALL_BLACK_HOLE_MASS
        * (1.0 + first_order[0])
        * strain_derivative
    )
    second_order = np.r_[
        np.einsum(
            "aij,ij->a", ode.SYMMETRIC_BASIS, strain_rate),
        (1.0 + first_order[0]) * clock_derivative,
    ]

    quadratic_value, quadratic_derivative = (
        ode.quadratic_value_and_derivative(
            first_order, radius_index))
    full_second_order = np.zeros(43)
    full_second_order[ode.SECOND_ORDER_INDICES] = second_order
    second_order_value = ode.metric_linear_combination(
        ode.SECOND_ORDER_COLUMNS[radius_index],
        full_second_order,
    )
    second_order_derivative = ode.metric_linear_combination(
        ode.SECOND_ORDER_DERIVATIVES[radius_index],
        full_second_order,
    )
    raw_value = ode.add_metrics(
        ode.metric_linear_combination(
            ode.FIRST_ORDER_COLUMNS[radius_index], first_order),
        quadratic_value,
        second_order_value,
    )
    raw_derivative = ode.add_metrics(
        ode.metric_linear_combination(
            ode.FIRST_ORDER_DERIVATIVES[radius_index],
            first_order,
        ),
        quadratic_derivative,
        second_order_derivative,
    )

    recovered_second_order = ode.fit_second_order(
        first_order, raw_value, raw_derivative, radius_index)
    recovered_first_order = ode.fit_first_order(
        raw_value,
        raw_derivative,
        radius_index,
        ode.add_metrics(quadratic_value, second_order_value),
        ode.add_metrics(
            quadratic_derivative, second_order_derivative),
    )
    first_order_exact.append(first_order)
    second_order_exact.append(second_order)
    first_order_recovered.append(recovered_first_order)
    second_order_recovered.append(recovered_second_order)

first_order_exact = np.asarray(first_order_exact)
second_order_exact = np.asarray(second_order_exact)
first_order_recovered = np.asarray(first_order_recovered)
second_order_recovered = np.asarray(second_order_recovered)

first_for_diagnostic = np.zeros((len(times), 1, 5, 13))
second_for_diagnostic = np.zeros((len(times), 1, 4, 7))
first_for_diagnostic[:, :, :4] = (
    first_order_recovered[:, None, None])
second_for_diagnostic[:] = (
    second_order_recovered[:, None, None])
diagnostics = ode.temporal_diagnostics(
    times, first_for_diagnostic, second_for_diagnostic)

print(
    "first-order block relative error:",
    np.linalg.norm(first_order_recovered - first_order_exact)
    / np.linalg.norm(first_order_exact),
)
print(
    "second-order block relative error:",
    np.linalg.norm(second_order_recovered - second_order_exact)
    / np.linalg.norm(second_order_exact),
)
print("clock local mismatch:", diagnostics["local_clock"][:, 0])
print("strain local mismatch:", diagnostics["local_strain"][:, 0])
print("clock integral mismatch:", diagnostics["integral_clock"][:, 0])
print("strain integral mismatch:", diagnostics["integral_strain"][:, 0])
