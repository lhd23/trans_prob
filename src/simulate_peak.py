"""Run the peak ensemble, saving correlation curves rather than particle snapshots."""

import argparse
import json

import numpy as np
from numba import njit, prange

from dynamics import advance_to_scale_factor, rank_ordered_particle_force
from initial_conditions import initial_batch
from settings import ROOT, load_config, numerical_metadata, run_directory


@njit(parallel=True, cache=True)
def correlations_in_batch(positions, edges):
    """Exact forward periodic pair counts, normalized by N^2 times bin width.

    Coordinates and bin edges use a unit box. The continuum normalization
    applies to particles displaced from lattice
    """
    realization_count, particle_count = positions.shape
    result = np.empty((realization_count, edges.size - 1), dtype=np.float64)
    for realization in prange(realization_count):
        ordered = np.sort(np.mod(positions[realization].astype(np.float64), 1.0))
        extended = np.empty(2 * particle_count, dtype=np.float64)
        extended[:particle_count] = ordered
        extended[particle_count:] = ordered + 1.0
        cumulative = np.empty(edges.size, dtype=np.int64)
        for edge_index in range(edges.size):
            pointer = 1
            count = 0
            for first in range(particle_count):
                if pointer < first + 1:
                    pointer = first + 1
                limit = ordered[first] + edges[edge_index]
                while pointer < 2 * particle_count and extended[pointer] < limit:
                    pointer += 1
                count += pointer - first - 1
            cumulative[edge_index] = count
        for radial_bin in range(edges.size - 1):
            continuum_count = particle_count**2 * (edges[radial_bin + 1] - edges[radial_bin])
            result[realization, radial_bin] = (
                (cumulative[radial_bin + 1] - cumulative[radial_bin]) / continuum_count - 1.0)
    return result


def radial_edges(config):
    lower = config["radial_min_mpc_over_h"]
    upper = config["radial_max_mpc_over_h"]
    width = config["radial_bin_width_mpc_over_h"]
    if not (0 < lower < upper < config["box_length_mpc_over_h"] / 2) or width <= 0:
        raise ValueError("radial bins must lie between zero and half the physical box")
    count = int(round((upper - lower) / width))
    if count < 1 or not np.isclose(count * width, upper - lower):
        raise ValueError("the radial range must contain a whole number of bins")
    return np.linspace(lower, upper, count + 1)


def simulate(config, directory):
    edges = radial_edges(config)
    times = 1.0 / (1.0 + np.asarray(config["redshifts"], dtype=float))
    if (config["realization_count"] < 2 or not np.all(np.isfinite(times))
            or times[0] <= config["initial_scale_factor"] or np.any(np.diff(times) <= 0)):
        raise ValueError("need at least two realizations and increasing times after initialization")
    directory.mkdir(parents=True, exist_ok=False)
    generator = np.random.default_rng(config["seed"])
    curve_sum = np.zeros((times.size, edges.size - 1))
    curve_square_sum = np.zeros_like(curve_sum)
    minimum_gap = np.inf
    step_counts = []

    for first in range(0, config["realization_count"], config["batch_size"]):
        count = min(config["batch_size"], config["realization_count"] - first)
        positions, momenta, details = initial_batch(config, count, generator)
        minimum_gap = min(minimum_gap, details["minimum_initial_neighbor_gap_code"])
        if first == 0:
            np.savez(directory / "spectrum.npz", mode_numbers=details["mode_numbers"],
                     wavenumbers=details["wavenumbers_h_over_mpc"],
                     power_1d_z0=details["power_1d_z0_mpc_over_h"])
        forces = rank_ordered_particle_force(positions, 1.0)
        current = config["initial_scale_factor"]
        for snapshot, target in enumerate(times):
            forces, steps = advance_to_scale_factor(
                positions, momenta, forces, current, float(target),
                config["maximum_delta_sqrt_scale_factor"],
                config["omega_matter"], config["omega_lambda"], 1.0)
            current = float(target)
            if first == 0:
                step_counts.append(steps)
            # Convert positions to 32-bit floating-point before pair counting.
            curves = correlations_in_batch(positions.astype(np.float32),
                                           edges / config["box_length_mpc_over_h"])
            curve_sum[snapshot] += curves.sum(axis=0)
            curve_square_sum[snapshot] += (curves**2).sum(axis=0)
        print(f"processed {first + count}/{config['realization_count']} realizations", flush=True)

    count = config["realization_count"]
    measured = curve_sum / count
    variance = np.maximum((curve_square_sum - count * measured**2) / (count - 1), 0)
    np.savez(directory / "correlations.npz", radial_edges=edges,
             scale_factors=times, correlation=measured,
             standard_error=np.sqrt(variance / count))
    metadata = numerical_metadata(config)
    metadata.update(
        workflow="peak", minimum_initial_neighbor_gap_code=minimum_gap,
        integration_steps_between_snapshots=step_counts,
        correlation_estimator="DD / [N^2 * bin_width / L] - 1",
        snapshot_position_precision="float32 before pair counting",
        uncertainty="standard error across independent realizations",
    )
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {directory}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(ROOT / "configs/peak.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    simulate(config, run_directory(config, "peak"))


if __name__ == "__main__":
    main()
