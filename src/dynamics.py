"""Exact periodic one-dimensional cosmological particle dynamics.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit(parallel=True, cache=True)
def rank_ordered_particle_force(positions: np.ndarray, box_length: float) -> np.ndarray:
    """Return the exact zero-softening periodic particle force.

    For each realization, the force is the displacement of every sorted
    particle from an equally spaced lattice, with the mean displacement removed.
    Removing that mean makes the expression independent of where the periodic
    cut is placed.  The returned array follows the original particle labels.
    """

    realization_count, particle_count = positions.shape
    spacing = box_length / particle_count
    forces = np.empty_like(positions)

    for realization in prange(realization_count):
        order = np.argsort(positions[realization])
        mean_offset = 0.0
        for rank in range(particle_count):
            label = order[rank]
            mean_offset += (
                positions[realization, label] - (rank + 0.5) * spacing
            )
        mean_offset /= particle_count

        for rank in range(particle_count):
            label = order[rank]
            forces[realization, label] = (
                positions[realization, label]
                - (rank + 0.5) * spacing
                - mean_offset
            )
    return forces


def direct_periodic_particle_force(
    positions: np.ndarray, box_length: float
) -> np.ndarray:
    """Slow pairwise reference for :func:`rank_ordered_particle_force`.

    This routine is intentionally quadratic and is only for validation.  Each
    periodic pair contributes the exact sawtooth Green-function force.
    """

    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2:
        raise ValueError("positions must have shape (realizations, particles)")
    particle_count = positions.shape[1]
    result = np.zeros_like(positions)
    for first in range(particle_count):
        difference = positions[:, first, None] - positions
        difference = (difference + 0.5 * box_length) % box_length - 0.5 * box_length
        result[:, first] = np.sum(
            difference / particle_count
            - box_length * np.sign(difference) / (2.0 * particle_count),
            axis=1,
        )
    return result


def expansion_rate(
    scale_factor: np.ndarray | float,
    omega_matter: float,
    omega_lambda: float,
) -> np.ndarray:
    """Return ``H(a) / H0`` for the configured flat matter-lambda model."""

    scale_factor = np.asarray(scale_factor, dtype=np.float64)
    return np.sqrt(omega_matter * scale_factor**-3 + omega_lambda)


def scale_factor_edges(
    lower: float, upper: float, maximum_delta_sqrt_scale_factor: float
) -> np.ndarray:
    """Make exact-endpoint steps uniform in square root of scale factor."""

    if upper <= lower:
        raise ValueError("upper scale factor must exceed lower scale factor")
    if maximum_delta_sqrt_scale_factor <= 0.0:
        raise ValueError("maximum step must be positive")
    square_root_lower = np.sqrt(lower)
    square_root_upper = np.sqrt(upper)
    step_count = int(
        np.ceil(
            (square_root_upper - square_root_lower)
            / maximum_delta_sqrt_scale_factor
        )
    )
    return np.linspace(square_root_lower, square_root_upper, step_count + 1) ** 2


def _gauss_legendre_integral(
    lower: np.ndarray,
    upper: np.ndarray,
    omega_matter: float,
    omega_lambda: float,
    coefficient: str,
) -> np.ndarray:
    """Integrate one time-dependent leapfrog coefficient."""

    nodes, weights = np.polynomial.legendre.leggauss(4)
    midpoint = 0.5 * (lower + upper)
    half_width = 0.5 * (upper - lower)
    samples = midpoint[:, None] + half_width[:, None] * nodes[None, :]
    expansion = expansion_rate(samples, omega_matter, omega_lambda)
    if coefficient == "drift":
        values = 1.0 / (samples**3 * expansion)
    elif coefficient == "kick":
        values = 1.5 * omega_matter / (samples**2 * expansion)
    else:
        raise ValueError(f"unknown coefficient: {coefficient}")
    return half_width * np.sum(values * weights[None, :], axis=1)


def leapfrog_coefficients(
    edges: np.ndarray, omega_matter: float, omega_lambda: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return first-kick, drift, and second-kick coefficients per step."""

    lower = np.asarray(edges[:-1], dtype=np.float64)
    upper = np.asarray(edges[1:], dtype=np.float64)
    square_root_midpoint = 0.5 * (np.sqrt(lower) + np.sqrt(upper))
    midpoint = square_root_midpoint**2
    first_kick = _gauss_legendre_integral(
        lower, midpoint, omega_matter, omega_lambda, "kick"
    )
    drift = _gauss_legendre_integral(
        lower, upper, omega_matter, omega_lambda, "drift"
    )
    second_kick = _gauss_legendre_integral(
        midpoint, upper, omega_matter, omega_lambda, "kick"
    )
    return first_kick, drift, second_kick


def advance_to_scale_factor(
    positions: np.ndarray,
    momenta: np.ndarray,
    forces: np.ndarray,
    lower_scale_factor: float,
    upper_scale_factor: float,
    maximum_delta_sqrt_scale_factor: float,
    omega_matter: float,
    omega_lambda: float,
    box_length: float,
) -> tuple[np.ndarray, int]:
    """Advance all realizations in place with kick-drift-kick steps."""

    edges = scale_factor_edges(
        lower_scale_factor, upper_scale_factor, maximum_delta_sqrt_scale_factor
    )
    first_kicks, drifts, second_kicks = leapfrog_coefficients(
        edges, omega_matter, omega_lambda
    )
    current_forces = forces
    for first_kick, drift, second_kick in zip(
        first_kicks, drifts, second_kicks
    ):
        momenta += first_kick * current_forces
        positions += drift * momenta
        positions %= box_length
        current_forces = rank_ordered_particle_force(positions, box_length)
        momenta += second_kick * current_forces
    return current_forces, edges.size - 1
