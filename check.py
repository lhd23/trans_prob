"""Small numerical checks; no production simulation or output archives."""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from analyse import histogram, merge_moments
from dynamics import (advance_to_scale_factor, direct_periodic_particle_force,
                      rank_ordered_particle_force)
from initial_conditions import initial_batch
from settings import load_config
from simulate_peak import correlations_in_batch
from simulate_transitions import snapshot_times
from transitions import select_ulam_pairs, measure_selected_signed_separations


def check_dynamics():
    positions = np.random.default_rng(1701).uniform(size=(3, 32))
    force = rank_ordered_particle_force(positions, 1.0)
    np.testing.assert_allclose(force, direct_periodic_particle_force(positions, 1.0), atol=1e-14)
    np.testing.assert_allclose(force, rank_ordered_particle_force((positions + 0.371) % 1, 1), atol=1e-14)
    np.testing.assert_allclose(force.sum(axis=1), 0, atol=1e-14)
    lattice = (np.arange(32) + 0.5) / 32
    displacement = 0.01 * np.sin(2 * np.pi * lattice)
    errors = []
    for step in (0.01, 0.005):
        positions = (lattice + 0.02 * displacement)[None, :]
        momenta = (0.02**1.5 * displacement)[None, :]
        force = rank_ordered_particle_force(positions, 1)
        advance_to_scale_factor(positions, momenta, force, 0.02, 1, step, 1, 0, 1)
        errors.append(np.max(np.abs(positions - lattice - displacement)))
    assert errors[1] < 3.5e-6 and errors[0] / errors[1] > 3.5


def check_pairs():
    positions = np.array([[0.0, 0.125, 0.25, 0.875], [0.0625, 0.1875, 0.3125, 0.9375]])
    pairs = select_ulam_pairs(positions, np.array([0.1875]), 0.125, 1, 1)
    expected = [(r, i, j) for r in range(2) for i in range(4) for j in range(4)
                if i != j and 0.125 <= (positions[r, j] - positions[r, i]) % 1 < 0.25]
    assert list(zip(pairs["realization_ids"], pairs["first_labels"], pairs["second_labels"])) == expected
    moved = (positions + np.array([0, 0.4, -0.3, 0.2])) % 1
    np.testing.assert_array_equal(measure_selected_signed_separations(moved, pairs, 1),
                                  [(moved[r, j] - moved[r, i] + 0.5) % 1 - 0.5 for r, i, j in expected])
    # Exhaustive pair-count reference catches the N(N-1) normalization bug.
    edges = np.array([0.05, 0.15, 0.35, 0.49])
    curves = correlations_in_batch(positions, edges)
    for row, x in enumerate(positions):
        separations = [(x[j] - x[i]) % 1 for i in range(x.size) for j in range(x.size) if i != j]
        counts, _ = np.histogram(separations, edges)
        np.testing.assert_allclose(curves[row], counts / (x.size**2 * np.diff(edges)) - 1)


def check_initial_conditions():
    config = load_config(ROOT / "configs/transitions.json")
    config["particle_count"] = 64
    whole = initial_batch(config, 5, np.random.default_rng(config["seed"]))
    generator = np.random.default_rng(config["seed"])
    first = initial_batch(config, 2, generator)
    second = initial_batch(config, 3, generator)
    for index in (0, 1):
        np.testing.assert_array_equal(whole[index], np.concatenate((first[index], second[index])))
    assert whole[2]["minimum_initial_neighbor_gap_code"] > 0
    times = snapshot_times(config)
    assert times.size == 24
    assert all(np.any(times == a) for a in (0.02, 0.1, 0.3, 0.6, 1, 10))


def check_analysis():
    samples = np.array([-0.5, -0.1, 0, 0.1, 0.49])
    counts, density = histogram(samples, np.array([-0.5, 0, 0.5]))
    np.testing.assert_array_equal(counts, [2, 3])
    assert np.isclose(np.sum(density * 0.5), 1)
    # Cropping the display must not renormalize the retained subset.
    _, density = histogram(samples, np.array([-0.2, 0.2]))
    assert np.isclose(density[0] * 0.4, 3 / 5)
    samples = np.random.default_rng(45).normal(0.3, 2, 101)
    count, mean, moment2, moment3 = 0, 0.0, 0.0, 0.0
    for batch in (samples[:17], samples[17:60], samples[60:]):
        mean, moment2, moment3 = merge_moments(count, mean, moment2, moment3, batch)
        count += batch.size
    np.testing.assert_allclose([mean, moment2/count, moment3/count],
                              [samples.mean(), samples.var(), np.mean((samples-samples.mean())**3)], atol=1e-14)


if __name__ == "__main__":
    check_dynamics()
    check_pairs()
    check_initial_conditions()
    check_analysis()
    print("Force, growing-mode integration, pair selection, labels, seed batching, and estimators passed.")
