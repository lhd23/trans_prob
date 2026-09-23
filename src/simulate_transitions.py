"""Save selected-pair separations for the histogram and cumulant plots."""

import argparse
import json

import numpy as np

from dynamics import advance_to_scale_factor, rank_ordered_particle_force
from initial_conditions import initial_batch
from settings import ROOT, load_config, numerical_metadata, run_directory
from transitions import measure_selected_signed_separations, select_ulam_pairs


def snapshot_times(config):
    logarithmic = np.geomspace(config["initial_scale_factor"],
                              config["final_scale_factor"],
                              config["logarithmic_snapshot_count"])
    extras = np.asarray(config["extra_scale_factors"], dtype=float)
    times = np.unique(np.concatenate((logarithmic, extras)))
    if (times.size != logarithmic.size + extras.size
            or times[0] != config["initial_scale_factor"]
            or times[-1] != config["final_scale_factor"]
            or np.any(np.diff(times) <= 0)):
        raise ValueError("snapshot times must be distinct and inside the evolution interval")
    return times


def simulate(config, directory):
    times = snapshot_times(config)
    # Checkpoints subdivide integration intervals but are not output snapshots.
    checkpoints = np.asarray(config["integration_checkpoints"], dtype=float)
    if (not np.all(np.isfinite(checkpoints)) or np.any(checkpoints < times[0])
            or np.any(checkpoints > times[-1])):
        raise ValueError("integration checkpoints must be inside the evolution interval")
    integration_times = np.unique(np.concatenate((times, checkpoints)))
    snapshot_indices = {float(time): index for index, time in enumerate(times)}
    targets = np.asarray(config["target_separation_fractions"], dtype=float)
    width = config["source_cell_width_fraction"]
    if (width <= 0 or np.any(targets <= 0) or np.any(targets >= 0.5)
            or np.any(np.diff(targets) >= 0)):
        raise ValueError("source cells need decreasing targets between zero and half a box")
    directory.mkdir(parents=True, exist_ok=False)
    generator = np.random.default_rng(config["seed"])
    counts = np.zeros((targets.size, config["realization_count"]), dtype=np.int64)
    minimum_gap = np.inf
    batch_files = []

    for first in range(0, config["realization_count"], config["batch_size"]):
        count = min(config["batch_size"], config["realization_count"] - first)
        positions, momenta, details = initial_batch(config, count, generator)
        minimum_gap = min(minimum_gap, details["minimum_initial_neighbor_gap_code"])
        if first == 0:
            np.savez(directory / "spectrum.npz",
                     mode_numbers=details["mode_numbers"],
                     wavenumbers=details["wavenumbers_h_over_mpc"],
                     power_1d_z0=details["power_1d_z0_mpc_over_h"])
        # Search one realization at a time, so temporary pair arrays do not
        # scale with the ensemble size. All non-self partners are considered.
        pairs = select_ulam_pairs(positions, targets, width, 1.0, 1)
        counts[:, first:first + count] = pairs["counts_by_target_and_realization"]
        separations = np.empty((times.size, pairs["realization_ids"].size),
                               dtype=np.float32)
        forces = rank_ordered_particle_force(positions, 1.0)
        current = config["initial_scale_factor"]
        steps_in_batch = 0
        for target_time in integration_times:
            if target_time > current:
                forces, steps = advance_to_scale_factor(
                    positions, momenta, forces, current, float(target_time),
                    config["maximum_delta_sqrt_scale_factor"],
                    config["omega_matter"], config["omega_lambda"], 1.0)
                steps_in_batch += steps
                current = float(target_time)
            if float(target_time) in snapshot_indices:
                index = snapshot_indices[float(target_time)]
                stored = measure_selected_signed_separations(positions, pairs, 1.0).astype(np.float32)
                stored[stored >= np.float32(0.5)] -= np.float32(1.0)
                separations[index] = stored
        filename = f"batch_{first:06d}.npz"
        np.savez_compressed(
            directory / filename,
            signed_separations=separations,
            target_offsets=pairs["target_offsets"],
            realization_ids=pairs["realization_ids"] + first,
            first_labels=pairs["first_labels"],
            second_labels=pairs["second_labels"],
            initial_separations=pairs["exact_initial_forward_separations"],
        )
        batch_files.append(filename)
        print(f"saved {first + count}/{config['realization_count']} realizations", flush=True)

    if np.any(counts.sum(axis=1) == 0):
        raise ValueError("a source cell contains no pairs; increase the sampling")
    metadata = numerical_metadata(config)
    metadata.update(
        workflow="transitions", batch_files=batch_files,
        snapshot_scale_factors=times.tolist(),
        pair_counts_by_target_and_realization=counts.tolist(),
        minimum_initial_neighbor_gap_code=minimum_gap,
        integration_step_count=steps_in_batch,
        stored_precision="float32 separations; float64 initial separations",
        conditioning="actual displaced initial forward separation in [q-width/2, q+width/2)",
        stored_orientation="signed periodic interval [-0.5, 0.5)",
    )
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"wrote {directory}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=str, default=str(ROOT / "configs/transitions.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    simulate(config, run_directory(config, "transitions"))


if __name__ == "__main__":
    main()
