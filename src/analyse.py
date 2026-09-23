"""Read transition batches and estimate histograms, cumulants, and trajectories."""

import argparse
import json
from pathlib import Path

import numpy as np


def histogram(samples, edges):
    """Return counts and probability density, normalized by ALL selected pairs.

    The integral is below one if the supplied display bins omit the tails.
    Initial conditioning cells and these destination bins need not match.
    """
    if samples.size == 0:
        raise ValueError("cannot estimate a histogram without selected pairs")
    edges = np.asarray(edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("histogram edges must increase strictly")
    counts, _ = np.histogram(samples, edges)
    return counts, counts / (samples.size * np.diff(edges))


def merge_moments(count, mean, moment2, moment3, samples):
    """Combine population central moments without subtracting large raw moments."""
    batch_count = samples.size
    batch_mean = samples.mean()
    centered = samples - batch_mean
    batch2 = np.sum(centered**2)
    batch3 = np.sum(centered**3)
    total = count + batch_count
    delta = batch_mean - mean
    combined3 = (moment3 + batch3
                 + delta**3 * count * batch_count * (count - batch_count) / total**2
                 + 3 * delta * (count * batch2 - batch_count * moment2) / total)
    combined2 = moment2 + batch2 + delta**2 * count * batch_count / total
    return mean + delta * batch_count / total, combined2, combined3


def analyse(directory, histogram_bins=270):
    metadata = json.loads((directory / "metadata.json").read_text())
    if metadata["workflow"] != "transitions":
        raise ValueError("this analysis needs transition batches")
    config = metadata["configuration"]
    times = np.asarray(metadata["snapshot_scale_factors"])
    targets = np.asarray(config["target_separation_fractions"])
    counts_by_realization = np.asarray(metadata["pair_counts_by_target_and_realization"])
    pair_counts = counts_by_realization.sum(axis=1)
    if histogram_bins < 1 or np.any(pair_counts <= 0):
        raise ValueError("need positive histogram bin and selected-pair counts")
    # Gnuplot draws the histogram-bin centres with histeps, so plateaus are not
    # shifted by half a bin.
    plot_edges = np.linspace(-0.13818682, 0.15, histogram_bins + 1)
    destination_edges = np.linspace(-0.5, 0.5, 1001)
    plot_counts = np.zeros((targets.size, times.size, histogram_bins), dtype=np.int64)
    transition_counts = np.zeros((targets.size, times.size, 1000), dtype=np.int64)
    means = np.zeros((targets.size, times.size))
    moment2 = np.zeros_like(means)
    moment3 = np.zeros_like(means)
    processed = np.zeros(targets.size, dtype=np.int64)
    observed_counts = np.zeros_like(counts_by_realization)
    generator = np.random.default_rng(config["seed"])
    selected = [generator.choice(count, min(75, count), replace=False) for count in pair_counts]
    trajectory_indices = np.flatnonzero(times <= 1.0)
    paths = np.full((targets.size, 75, trajectory_indices.size), np.nan)
    trajectory_ids = np.full((targets.size, 75), -1, dtype=np.int64)
    for row, indices in enumerate(selected):
        trajectory_ids[row, :indices.size] = indices

    for filename in metadata["batch_files"]:
        with np.load(directory / filename, allow_pickle=False) as batch:
            signed = batch["signed_separations"]
            offsets = batch["target_offsets"]
            initial = batch["initial_separations"]
            realizations = batch["realization_ids"]
            first_labels, second_labels = batch["first_labels"], batch["second_labels"]
            if signed.shape != (times.size, int(offsets[-1])):
                raise ValueError(f"invalid separation shape in {filename}")
            if (not np.all(np.isfinite(signed)) or np.any(signed < -0.5)
                    or np.any(signed >= 0.5) or np.any(first_labels == second_labels)):
                raise ValueError(f"invalid signed separation or self-pair in {filename}")
            for row, target in enumerate(targets):
                start, stop = map(int, offsets[row:row + 2])
                count = stop - start
                if count == 0:
                    continue
                half_width = config["source_cell_width_fraction"] / 2
                source = initial[start:stop]
                if np.any(source < target - half_width) or np.any(source >= target + half_width):
                    raise ValueError("a selected pair is outside its displaced initial cell")
                expected = (source + 0.5) % 1 - 0.5
                if not np.allclose(signed[0, start:stop], expected, rtol=0, atol=6e-8):
                    raise ValueError("initial signed separations disagree with selection")
                observed_counts[row] += np.bincount(
                    realizations[start:stop], minlength=config["realization_count"])
                local_indices = selected[row] - processed[row]
                keep = np.flatnonzero((local_indices >= 0) & (local_indices < count))
                for slot in keep:
                    paths[row, slot] = signed[trajectory_indices, start + local_indices[slot]].astype(np.float64) - target
                for snapshot in range(times.size):
                    samples = signed[snapshot, start:stop].astype(np.float64)
                    transition_counts[row, snapshot] += histogram(samples, destination_edges)[0]
                    changes = samples - target
                    plot_counts[row, snapshot] += histogram(changes, plot_edges)[0]
                    means[row, snapshot], moment2[row, snapshot], moment3[row, snapshot] = merge_moments(
                        int(processed[row]), means[row, snapshot], moment2[row, snapshot],
                        moment3[row, snapshot], changes / target)
                processed[row] += count
        print(f"analysed {filename}", flush=True)

    if not np.array_equal(observed_counts, counts_by_realization):
        raise ValueError("stored pair counts disagree with the batch catalogues")
    if not np.all(transition_counts.sum(axis=2) == pair_counts[:, None]):
        raise ValueError("destination histograms lost selected pairs")
    cumulants = np.stack((means, moment2 / pair_counts[:, None],
                         moment3 / pair_counts[:, None]), axis=2)
    np.savez_compressed(
        directory / "analysis.npz", times=times, targets=targets, pair_counts=pair_counts,
        destination_edges=destination_edges, transition_counts=transition_counts,
        transition_probabilities=transition_counts / pair_counts[:, None, None],
        plot_edges=plot_edges, plot_counts=plot_counts,
        density=plot_counts / (pair_counts[:, None, None] * np.diff(plot_edges)),
        cumulants=cumulants, trajectory_times=times[trajectory_indices],
        trajectories=paths, trajectory_pair_indices=trajectory_ids)
    print(f"wrote {directory / 'analysis.npz'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="saved transition run directory")
    parser.add_argument("--histogram-bins", type=int, default=270)
    args = parser.parse_args()
    analyse(args.directory.resolve(), args.histogram_bins)


if __name__ == "__main__":
    main()
