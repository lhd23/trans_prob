"""Gaussian first-order Lagrangian initial conditions."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from dynamics import expansion_rate


def file_sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_planck_spectrum(
    path: Path, expected_sha256: str
) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate the supplied redshift-zero three-dimensional power."""

    actual_sha256 = file_sha256(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"spectrum checksum is {actual_sha256}, expected {expected_sha256}"
        )
    table = np.loadtxt(path, comments="#", dtype=np.float64)
    if table.ndim != 2 or table.shape[1] != 2:
        raise ValueError("the spectrum must contain exactly two numeric columns")
    wavenumber = table[:, 0]
    power_3d = table[:, 1]
    if not np.all(np.isfinite(table)):
        raise ValueError("the spectrum contains a non-finite value")
    if np.any(wavenumber <= 0.0) or np.any(np.diff(wavenumber) <= 0.0):
        raise ValueError("spectrum wavenumbers must be positive and increasing")
    if np.any(power_3d <= 0.0):
        raise ValueError("spectrum power must be positive")
    return wavenumber, power_3d


def interpolate_power_logarithmically(
    source_wavenumber: np.ndarray,
    source_power: np.ndarray,
    requested_wavenumber: np.ndarray,
) -> np.ndarray:
    """Interpolate positive power in logarithmic wavenumber and power."""

    requested_wavenumber = np.asarray(requested_wavenumber, dtype=np.float64)
    if np.any(requested_wavenumber < source_wavenumber[0]) or np.any(
        requested_wavenumber > source_wavenumber[-1]
    ):
        raise ValueError("requested modes require forbidden spectrum extrapolation")
    return np.exp(
        np.interp(
            np.log(requested_wavenumber),
            np.log(source_wavenumber),
            np.log(source_power),
        )
    )


def linear_growth_at(
    scale_factor: float, omega_matter: float, omega_lambda: float
) -> tuple[float, float]:
    """Return growth normalized at scale factor one and ``d ln D/d ln a``.

    The integration begins deeply in matter domination with the growing-mode
    conditions ``D=a`` and ``dD/d ln(a)=D``.  Radiation is intentionally absent,
    matching the matter-lambda dynamics used by the particle integrator.
    """

    if scale_factor <= 0.0:
        raise ValueError("growth is required at a positive scale factor")
    early_scale_factor = min(1.0e-5, scale_factor / 100.0)
    final_log_scale_factor = max(0.0, float(np.log(scale_factor)))

    def derivative(log_scale_factor: float, state: np.ndarray) -> np.ndarray:
        current_scale_factor = np.exp(log_scale_factor)
        expansion_squared = (
            omega_matter * current_scale_factor**-3 + omega_lambda
        )
        matter_fraction = (
            omega_matter * current_scale_factor**-3 / expansion_squared
        )
        logarithmic_expansion_derivative = -1.5 * matter_fraction
        growth, logarithmic_velocity = state
        return np.array(
            [
                logarithmic_velocity,
                -(2.0 + logarithmic_expansion_derivative)
                * logarithmic_velocity
                + 1.5 * matter_fraction * growth,
            ],
            dtype=np.float64,
        )

    solution = solve_ivp(
        derivative,
        (np.log(early_scale_factor), final_log_scale_factor),
        np.array([early_scale_factor, early_scale_factor], dtype=np.float64),
        rtol=2.0e-11,
        atol=2.0e-13,
        dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(f"linear growth integration failed: {solution.message}")
    growth_at_scale_factor, velocity_at_scale_factor = solution.sol(
        np.log(scale_factor)
    )
    growth_at_one = solution.sol(0.0)[0]
    normalized_growth = growth_at_scale_factor / growth_at_one
    logarithmic_growth_rate = velocity_at_scale_factor / growth_at_scale_factor
    return float(normalized_growth), float(logarithmic_growth_rate)


def supported_modes(
    particle_count: int,
    physical_box_length: float,
    maximum_wavenumber: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return positive non-Nyquist periodic modes below the configured limit."""

    maximum_mode_number = (particle_count - 1) // 2
    mode_numbers = np.arange(1, maximum_mode_number + 1, dtype=np.int64)
    wavenumbers = 2.0 * np.pi * mode_numbers / physical_box_length
    keep = wavenumbers <= maximum_wavenumber
    if not np.any(keep):
        raise ValueError("the configuration initializes no Fourier modes")
    return mode_numbers[keep], wavenumbers[keep]


def make_1lpt_initial_conditions(
    realization_count: int,
    particle_count: int,
    code_box_length: float,
    physical_box_length: float,
    initial_scale_factor: float,
    maximum_wavenumber: float,
    omega_matter: float,
    omega_lambda: float,
    seed: int,
    spectrum_path: Path,
    spectrum_sha256: str,
    generator: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Displace a quiet lattice with Gaussian growing-mode realizations."""

    source_wavenumber, source_power_3d = load_planck_spectrum(
        spectrum_path, spectrum_sha256
    )
    mode_numbers, wavenumber = supported_modes(
        particle_count, physical_box_length, maximum_wavenumber
    )
    power_3d_z0 = interpolate_power_logarithmically(
        source_wavenumber, source_power_3d, wavenumber
    )
    power_1d_z0 = wavenumber**2 * power_3d_z0 / (2.0 * np.pi)
    power_code_z0 = power_1d_z0 / physical_box_length
    growth, growth_rate = linear_growth_at(
        initial_scale_factor, omega_matter, omega_lambda
    )
    power_code_initial = growth**2 * power_code_z0

    frequency_count = particle_count // 2 + 1
    density_fourier = np.zeros(
        (realization_count, frequency_count), dtype=np.complex128
    )
    if generator is None:
        generator = np.random.default_rng(seed)
    gaussian = generator.standard_normal(
        (realization_count, mode_numbers.size, 2)
    )
    amplitudes = particle_count * np.sqrt(power_code_initial / 2.0)
    density_fourier[:, mode_numbers] = amplitudes[None, :] * (
        gaussian[:, :, 0] + 1j * gaussian[:, :, 1]
    )

    code_wavenumber = 2.0 * np.pi * mode_numbers / code_box_length
    displacement_fourier = np.zeros_like(density_fourier)
    displacement_fourier[:, mode_numbers] = (
        1j * density_fourier[:, mode_numbers] / code_wavenumber[None, :]
    )
    displacement = np.fft.irfft(displacement_fourier, n=particle_count, axis=1)
    spacing = code_box_length / particle_count
    lagrangian_positions = (
        np.arange(particle_count, dtype=np.float64) + 0.5
    ) * spacing
    unwrapped_positions = lagrangian_positions[None, :] + displacement

    neighbor_gaps = np.empty_like(unwrapped_positions)
    neighbor_gaps[:, :-1] = np.diff(unwrapped_positions, axis=1)
    neighbor_gaps[:, -1] = (
        unwrapped_positions[:, 0]
        + code_box_length
        - unwrapped_positions[:, -1]
    )
    minimum_gap = float(np.min(neighbor_gaps))
    if minimum_gap <= 0.0:
        raise ValueError(
            "a first-order Lagrangian configuration has crossed at the "
            "initial scale factor"
        )

    positions = unwrapped_positions % code_box_length
    initial_expansion = float(
        expansion_rate(initial_scale_factor, omega_matter, omega_lambda)
    )
    momenta = (
        initial_scale_factor**2
        * initial_expansion
        * growth_rate
        * displacement
    )
    details = {
        "generator": "numpy.random.Generator(PCG64)",
        "position_distribution": "quiet lattice plus Gaussian 1LPT displacement",
        "momentum_model": "growing mode",
        "growth_factor_at_initial_scale_factor": growth,
        "logarithmic_growth_rate_at_initial_scale_factor": growth_rate,
        "minimum_initial_neighbor_gap_code": minimum_gap,
        "mode_numbers": mode_numbers,
        "wavenumbers_h_over_mpc": wavenumber,
        "power_3d_z0_mpc_over_h_cubed": power_3d_z0,
        "power_1d_z0_mpc_over_h": power_1d_z0,
        "power_code_z0": power_code_z0,
        "power_code_initial": power_code_initial,
        "displacement_code": displacement,
    }
    return positions, momenta, details


def initial_batch(config: dict, count: int, generator: np.random.Generator):
    """Draw the next batch from one continuous random stream."""
    from settings import SPECTRUM_PATH, SPECTRUM_SHA256

    cutoff = config["maximum_initialized_wavenumber_h_over_mpc"]
    return make_1lpt_initial_conditions(
        count, config["particle_count"], 1.0, config["box_length_mpc_over_h"],
        config["initial_scale_factor"], np.inf if cutoff is None else cutoff,
        config["omega_matter"], config["omega_lambda"], config["seed"],
        SPECTRUM_PATH, SPECTRUM_SHA256, generator=generator,
    )
