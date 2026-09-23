"""Paths and settings shared by the two small workflows."""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPECTRUM_PATH = ROOT / "data" / "Pk_Planck18_large.dat"
SPECTRUM_SHA256 = "d4446d8d7674ec1a89f9f4ad5537a8d04c5179d7b5acaa4032976d6c69b3e124"


def load_config(path):
    config = json.loads(Path(path).read_text())
    for key in ("particle_count", "realization_count", "batch_size"):
        value = config[key]
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{key} must be a positive integer")
    if config["particle_count"] < 4 or config["particle_count"] % 2:
        raise ValueError("particle_count must be an even integer of at least four")
    for key in ("box_length_mpc_over_h", "initial_scale_factor",
                "maximum_delta_sqrt_scale_factor"):
        if not np.isfinite(config[key]) or config[key] <= 0:
            raise ValueError(f"{key} must be positive and finite")
    if not (0 < config["omega_matter"] <= 1 and config["omega_lambda"] >= 0
            and np.isclose(config["omega_matter"] + config["omega_lambda"], 1)):
        raise ValueError("the background must be flat matter plus cosmological constant")
    cutoff = config["maximum_initialized_wavenumber_h_over_mpc"]
    if cutoff is not None and (not np.isfinite(cutoff) or cutoff <= 0):
        raise ValueError("the spectrum cutoff must be positive or null")
    return config


def run_directory(config, workflow):
    """Counts and box length propagate to every downstream filename."""
    name = (f"{workflow}_N{config['particle_count']}_R{config['realization_count']}"
            f"_L{config['box_length_mpc_over_h']:g}")
    return ROOT / "data" / name


def numerical_metadata(config):
    return {
        "configuration": config,
        "initial_conditions": "Gaussian first-order Lagrangian growing mode",
        "spectrum_file": "data/Pk_Planck18_large.dat",
        "spectrum_sha256": SPECTRUM_SHA256,
        "spectrum_conversion": "P_1D(k) = k^2 P_3D(k) / (2*pi)",
        "force": "exact periodic rank-ordered equal-mass particle force; no softening",
        "integrator": "scale-factor-integrated kick-drift-kick; canonical momentum",
        "dynamics_precision": "float64",
        "random_generator": "numpy.random.Generator(PCG64)",
    }
