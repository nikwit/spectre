# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Experimental Chalmers et al. hp estimator for Legendre tensor products.

This is an offline research prototype, not an evolution AMR criterion.
See https://doi.org/10.1016/j.jcpx.2019.100023, section 2.4.6. We fit each
coefficient sequence separately, not an RMS power monitor. The paper does not
specify a fit window or handling of zero coefficients. Here we fit log
amplitudes
over all coefficients above a roundoff floor, retaining their original indices.
These choices, and the treatment of numerically absent tails, must be calibrated
before using this estimator in an evolution.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.special import hyp2f1

from spectre.Spectral import (
    Basis,
    Mesh1D,
    Quadrature,
    collocation_points,
    nodal_to_modal_matrix,
)


@dataclass(frozen=True)
class SequenceEstimate:
    """Error is a reference-interval L2 estimate, not a physical/global bound."""

    error: float
    sigma: float
    log_fit_rms: float
    used_modes: int
    status: str


@dataclass(frozen=True)
class DirectionEstimate:
    error: float
    sigma: float
    sequences: tuple[SequenceEstimate, ...]


def estimate_sequence(coefficients, *, coefficient_floor=None):
    """Estimate one sequence in the standard P_k(1)=1 Legendre normalization.

    The error squared is 2*c_p**2/(2*p+1) plus the extrapolated squared L2 tail.
    Natural logarithms give the paper's sigma convention. At least three
    significant coefficients are needed to fit. Symmetry zeros are skipped, not
    replaced by tiny numbers or assigned consecutive indices. Two absent final
    modes are treated as a numerically terminated polynomial; this is NOT proof
    that the underlying function is resolved (aliasing can fool this test).

    A nondecaying fit has an infinite extrapolated tail. An insufficient fit
    also
    returns infinite error, with NaN sigma, to make the uncertainty explicit.
    """
    coeffs = np.asarray(coefficients, dtype=float)
    if coeffs.ndim != 1 or len(coeffs) < 2:
        raise ValueError(
            "Expected a one-dimensional sequence of at least 2 modes"
        )
    if not np.all(np.isfinite(coeffs)):
        raise ValueError("Modal coefficients must be finite")
    amplitudes = np.abs(coeffs)
    scale = np.max(amplitudes)
    if coefficient_floor is None:
        coefficient_floor = 100 * np.finfo(float).eps * scale
    if not np.isfinite(coefficient_floor) or coefficient_floor < 0:
        raise ValueError("The coefficient floor must be finite and nonnegative")
    indices = np.flatnonzero(amplitudes > coefficient_floor)
    if len(indices) == 0:
        return SequenceEstimate(0.0, np.nan, np.nan, 0, "below_floor")
    if indices[-1] < len(coeffs) - 2:
        return SequenceEstimate(
            0.0, np.nan, np.nan, len(indices), "terminated_tail"
        )
    if len(indices) < 3:
        return SequenceEstimate(
            np.inf, np.nan, np.nan, len(indices), "insufficient_modes"
        )

    log_amplitudes = np.log(amplitudes[indices] / scale)
    centered_indices = indices - np.mean(indices)
    slope = np.dot(
        centered_indices, log_amplitudes - np.mean(log_amplitudes)
    ) / np.dot(centered_indices, centered_indices)
    intercept = np.mean(log_amplitudes) - slope * np.mean(indices)
    sigma = -float(slope)
    residual = float(
        np.sqrt(np.mean((log_amplitudes - intercept - slope * indices) ** 2))
    )
    if sigma <= 0:
        return SequenceEstimate(
            np.inf, sigma, residual, len(indices), "nondecaying"
        )

    # Sum the infinite tail analytically. For q=exp(-2*sigma),
    # sum_{k=p+1}^inf 2*A^2*q^k/(2*k+1)
    # = 2*a_{p+1}^2/(2*p+3) * 2F1(1, p+3/2; p+5/2; q).
    degree = len(coeffs) - 1
    q = np.exp(-2 * sigma)
    first_omitted = scale * np.exp(intercept - sigma * (degree + 1))
    tail_factor = (
        2 / (2 * degree + 3) * hyp2f1(1, degree + 1.5, degree + 2.5, q)
    )
    error = np.hypot(
        np.sqrt(2 / (2 * degree + 1)) * amplitudes[-1],
        first_omitted * np.sqrt(tail_factor),
    )
    return SequenceEstimate(float(error), sigma, residual, len(indices), "fit")


def estimate_tensor(coefficients):
    """Return directional max errors and min decay rates in one to three D.

    Array axes are logical directions; axis d indexes the degree in direction d.
    A single global coefficient floor prevents roundoff-only transverse
    sequences from controlling the minimum. Numerically terminated sequences
    do not contribute a slope. Any active sequence with insufficient data makes
    the directional slope undetermined. This extends the paper's 2D recipe to
    3D; it does not supply a global multidimensional error bound.
    """
    coeffs = np.asarray(coefficients, dtype=float)
    if not 1 <= coeffs.ndim <= 3 or min(coeffs.shape) < 2:
        raise ValueError(
            "Expected a 1D, 2D, or 3D tensor with >=2 modes per axis"
        )
    if not np.all(np.isfinite(coeffs)):
        raise ValueError("Modal coefficients must be finite")
    floor = 100 * np.finfo(float).eps * np.max(np.abs(coeffs))
    result = []
    for axis, extent in enumerate(coeffs.shape):
        lines = np.moveaxis(coeffs, axis, 0).reshape(extent, -1)
        estimates = tuple(
            estimate_sequence(line, coefficient_floor=floor) for line in lines.T
        )
        slopes = [e.sigma for e in estimates if np.isfinite(e.sigma)]
        sigma = min(slopes) if slopes else np.nan
        if any(e.status == "insufficient_modes" for e in estimates):
            sigma = np.nan
        result.append(
            DirectionEstimate(max(e.error for e in estimates), sigma, estimates)
        )
    return tuple(result)


def refinement_choice(estimate, tolerance, *, threshold=1.0):
    """Return 'hold', 'p', 'h', or 'undetermined'; no silent fallback for NaNs."""
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("Tolerance must be finite and nonnegative")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("Threshold must be finite and positive")
    if estimate.error <= tolerance:
        return "hold"
    if not np.isfinite(estimate.sigma):
        return "undetermined"
    return "p" if estimate.sigma > threshold else "h"


@lru_cache(maxsize=None)
def legendre_grid(extent):
    """SpECTRE's LGL nodes and nodal-to-modal transform (cached, read-only)."""
    mesh = Mesh1D(extent, Basis.Legendre, Quadrature.GaussLobatto)
    nodes = np.array(collocation_points(mesh))
    transform = np.array(nodal_to_modal_matrix(mesh))
    nodes.flags.writeable = False
    transform.flags.writeable = False
    return nodes, transform


def to_modal(nodal_values):
    """Apply SpECTRE's transforms along each axis of a shaped nodal array."""
    coeffs = np.asarray(nodal_values, dtype=float)
    if not 1 <= coeffs.ndim <= 3:
        raise ValueError("Expected one to three logical axes")
    for axis, extent in enumerate(coeffs.shape):
        transform = legendre_grid(extent)[1]
        coeffs = np.moveaxis(
            np.tensordot(transform, coeffs, axes=(1, axis)), 0, axis
        )
    return coeffs


def sample_modal(function, extents, bounds=None):
    """Independently sample an analytic function on a box and transform it."""
    if bounds is None:
        bounds = [(-1.0, 1.0)] * len(extents)
    if len(extents) != len(bounds):
        raise ValueError("One bound pair is required per dimension")
    coordinates = []
    for extent, (lower, upper) in zip(extents, bounds):
        if not np.isfinite([lower, upper]).all() or upper <= lower:
            raise ValueError("Bounds must be finite and increasing")
        coordinates.append(
            lower + (legendre_grid(extent)[0] + 1) * (upper - lower) / 2
        )
    points = np.meshgrid(*coordinates, indexing="ij")
    values = np.broadcast_to(function(*points), tuple(extents))
    return to_modal(values)
