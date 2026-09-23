"""Add the spectrum-based models and draw the five-panel peak plot."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from initial_conditions import linear_growth_at


def zeldovich_correlations(power, wavenumbers, box_length, growth_factors, separations):
    """Evaluate periodic pair conservation for Gaussian displacements."""
    particle_count = 2 * (power.size - 1)
    cell_width = box_length / particle_count
    displacement_power = np.zeros_like(power)
    displacement_power[1:] = power[1:] / wavenumbers[1:]**2
    correlation = particle_count / box_length * np.fft.irfft(displacement_power, n=particle_count)
    relative_variance = np.fft.fftshift(np.maximum(2 * (correlation[0] - correlation), 0))
    lagrangian_separation = (np.arange(particle_count) - particle_count // 2) * cell_width
    result = np.empty((growth_factors.size, separations.size))
    for index, growth in enumerate(growth_factors):
        variance = growth**2 * relative_variance
        valid = variance > np.finfo(float).tiny
        normalization = np.zeros_like(variance)
        normalization[valid] = 1 / np.sqrt(2 * np.pi * variance[valid])
        inverse_variance = np.zeros_like(variance)
        inverse_variance[valid] = 1 / variance[valid]
        for first in range(0, separations.size, 32):
            current = separations[first:first + 32, None]
            distance = ((current - lagrangian_separation + 0.5 * box_length)
                        % box_length - 0.5 * box_length)
            exponent = -0.5 * distance**2 * inverse_variance
            exponent[:, ~valid] = -np.inf
            result[index, first:first + current.shape[0]] = (
                cell_width * np.sum(np.exp(exponent) * normalization, axis=1) - 1)
    return result


def plot(directory):
    metadata = json.loads((directory / "metadata.json").read_text())
    if metadata["workflow"] != "peak":
        raise ValueError("this plot needs a peak simulation")
    config = metadata["configuration"]
    with np.load(directory / "correlations.npz", allow_pickle=False) as data:
        edges = data["radial_edges"]
        times = data["scale_factors"]
        measured = data["correlation"]
        errors = data["standard_error"]
    if times.size != 5:
        raise ValueError("the peak figure has five redshift panels")
    centers = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]
    count = config["particle_count"]
    length = config["box_length_mpc_over_h"]
    wavenumbers = 2 * np.pi * np.fft.rfftfreq(count, d=length / count)
    power = np.zeros(wavenumbers.size)
    with np.load(directory / "spectrum.npz", allow_pickle=False) as spectrum:
        power[spectrum["mode_numbers"]] = spectrum["power_1d_z0"]
    growth = np.array([linear_growth_at(a, config["omega_matter"], config["omega_lambda"])[0]
                       for a in times])
    bin_window = np.sinc(wavenumbers * width / (2 * np.pi))
    linear_grid = count / length * np.fft.irfft(power * bin_window, n=count)
    linear = np.interp(centers, np.arange(count // 2 + 1) * length / count,
                       linear_grid[:count // 2 + 1])
    offsets = np.linspace(-width / 2, width / 2, 9)
    samples = zeldovich_correlations(power, wavenumbers, length, growth,
                                     (centers[:, None] + offsets).ravel())
    samples = samples.reshape(times.size, centers.size, offsets.size)
    # Explicit trapezoidal average works with both NumPy 1 and NumPy 2.
    model = np.sum((samples[:, :, :-1] + samples[:, :, 1:]) * np.diff(offsets) / 2, axis=2) / width
    tables = ROOT / "output/tables" / directory.name
    tables.mkdir(parents=True, exist_ok=True)
    for index in range(times.size):
        np.savetxt(tables / f"correlation_{index + 1}.dat",
                   np.column_stack((centers, measured[index] / growth[index]**2,
                                    model[index] / growth[index]**2, linear,
                                    errors[index] / growth[index]**2)),
                   header="r_mpc_over_h measured_over_D2 model_over_D2 linear standard_error_over_D2")
    amplitude = np.max(np.abs(np.concatenate((measured / growth[:, None]**2,
                                            model / growth[:, None]**2))))
    limit = max(0.02, np.ceil(1.1 * amplitude / 0.01) * 0.01)
    tick_candidates = np.array([1, 2, 5, 10]) * 10.0**np.floor(np.log10(limit / 3))
    tick = tick_candidates[np.argmin(np.abs(tick_candidates - limit / 3))]
    output = ROOT / "output/pdf" / f"{directory.name}.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    variables = (f"TABLES={json.dumps(str(tables))};OUTPUT={json.dumps(str(output))};"
                 f"YMAX={limit:g};YTICK={tick:g};"
                 f"XMIN={edges[0]:g};XMAX={edges[-1]:g};"
                 "array REDSHIFT[5] = [" + ",".join(f"{1/a - 1:g}" for a in times) + "]")
    subprocess.run(["gnuplot", "-e", variables, str(ROOT / "plots/peak.gp")], check=True)
    print(f"wrote {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    plot(args.directory.resolve())
