# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Offline symmetry-aware alternatives to the experimental Chalmers fit.

Both methods retain a measured high-mode guard and omit deliberately filtered
modes before calling these functions. Neither bounds unsampled structure.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import hyp2f1


@dataclass(frozen=True)
class Estimate:
    error: float
    sigma: float
    residual: float
    insufficient: int = 0
    ignored: int = 0


def fit(indices, amplitudes, floor):
    selected = amplitudes > floor
    if np.count_nonzero(selected) < 3:
        return np.nan, np.nan, np.nan
    x = indices[selected]
    y = np.log(amplitudes[selected] / np.max(amplitudes))
    slope, intercept = np.polyfit(x, y, 1)
    return (
        -float(slope),
        float(intercept + np.log(np.max(amplitudes))),
        float(np.sqrt(np.mean((y - intercept - slope * x) ** 2))),
    )


def tail(first, stride, sigma, log_amplitude):
    """Squared Legendre tail starting at an explicitly specified degree."""
    b = (2 * first + 1) / (2 * stride)
    return (
        np.exp(2 * (log_amplitude - sigma * first))
        * 2
        / (2 * first + 1)
        * hyp2f1(1, b, b + 1, np.exp(-2 * stride * sigma))
    )


def separate_parity(coefficients, budget=0.0):
    """Fit both parities; budget at most 1% for negligible resolved variation.

    k=0 is omitted only from the negligibility test, since it cannot require
    directional refinement. Two absent terminal modes *of the same parity*
    denote a terminated polynomial. Each active parity guards its last mode.
    """
    a = np.abs(np.asarray(coefficients, dtype=float))
    if (
        a.ndim != 1
        or len(a) < 2
        or not np.all(np.isfinite(a))
        or not np.isfinite(budget)
        or budget < 0
    ):
        raise ValueError("Expected finite coefficients and nonnegative budget")
    floor = 100 * np.finfo(float).eps * np.max(a)
    errors, slopes, residuals = [], [], []
    insufficient = ignored = 0
    for parity in (0, 1):
        k = np.arange(parity, len(a), 2)
        values = a[k]
        varying = k != 0
        norm = np.linalg.norm(
            values[varying] * np.sqrt(2 / (2 * k[varying] + 1))
        )
        if budget > 0 and norm <= 0.01 * budget / np.sqrt(2):
            errors.append(norm)
            ignored += 1
            continue
        significant = np.flatnonzero(values > floor)
        if not len(significant) or significant[-1] < len(k) - 2:
            errors.append(0.0)
            continue
        sigma, intercept, residual = fit(k, values, floor)
        if np.isnan(sigma):
            insufficient += 1
            errors.append(np.inf)
            continue
        slopes.append(sigma)
        residuals.append(residual)
        errors.append(
            np.sqrt(
                2 / (2 * k[-1] + 1) * values[-1] ** 2
                + tail(k[-1] + 2, 2, sigma, intercept)
            )
            if sigma > 0
            else np.inf
        )
    return Estimate(
        float(np.linalg.norm(errors)),
        min(slopes) if slopes and not insufficient else np.nan,
        max(residuals, default=0.0),
        insufficient,
        ignored,
    )


def paired(coefficients, budget=0.0):
    """Fit complete adjacent amplitude pairs against original degree spacing.

    Pair energy uses the larger of the two Legendre weights as a conservative
    envelope. An incomplete final pair is excluded from the slope, but its
    measured coefficient is retained in the error. The extrapolated envelope
    is still only an indicator, particularly for unequal parity decay rates.
    """
    a = np.abs(np.asarray(coefficients, dtype=float))
    if (
        a.ndim != 1
        or len(a) < 2
        or not np.all(np.isfinite(a))
        or not np.isfinite(budget)
        or budget < 0
    ):
        raise ValueError("Expected finite coefficients and nonnegative budget")
    k = np.arange(len(a))
    norm = np.linalg.norm(a[1:] * np.sqrt(2 / (2 * k[1:] + 1)))
    if budget > 0 and norm <= 0.01 * budget:
        return Estimate(float(norm), np.nan, 0.0, ignored=1)
    n = len(a) // 2
    values = np.hypot(a[: 2 * n : 2], a[1 : 2 * n : 2])
    x = 2 * np.arange(n)
    floor = 100 * np.finfo(float).eps * np.max(a)
    significant = np.flatnonzero(values > floor)
    if not len(significant) or (
        significant[-1] < n - 2 and np.all(a[2 * n :] <= floor)
    ):
        return Estimate(0.0, np.nan, 0.0)
    sigma, intercept, residual = fit(x, values, floor)
    if np.isnan(sigma):
        return Estimate(np.inf, np.nan, np.nan, insufficient=1)
    if sigma <= 0:
        return Estimate(np.inf, sigma, residual)
    guard = np.sum(a[2 * (n - 1) :] ** 2 * 2 / (2 * k[2 * (n - 1) :] + 1))
    return Estimate(
        float(np.sqrt(guard + tail(2 * n, 2, sigma, intercept))),
        sigma,
        residual,
    )
