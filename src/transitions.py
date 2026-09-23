"""Select pairs by displaced initial separation and follow their labels."""

from __future__ import annotations

import numpy as np


def forward_periodic_separation(
    first_position: np.ndarray,
    second_position: np.ndarray,
    box_length: float,
) -> np.ndarray:
    """Return ``(second - first) mod L`` in the interval ``[0, L)``."""

    return (second_position - first_position) % box_length


def signed_periodic_separation(
    first_position: np.ndarray,
    second_position: np.ndarray,
    box_length: float,
) -> np.ndarray:
    """Return the signed minimum-image separation in ``[-L/2, L/2)``."""

    difference = second_position - first_position
    return (difference + 0.5 * box_length) % box_length - 0.5 * box_length


def select_ulam_pairs(
    initial_positions: np.ndarray,
    target_separations: np.ndarray,
    cell_width: float,
    box_length: float,
    realization_batch_size: int,
) -> dict[str, np.ndarray]:
    """Select every oriented non-self pair in each half-open Ulam cell.

    Membership uses the actual forward periodic initial separation and the
    interval ``[q - width/2, q + width/2)``.  Realizations are batched so the
    exhaustive anchor-by-partner search never creates one ensemble-sized
    three-dimensional temporary array.
    """

    positions = np.asarray(initial_positions, dtype=np.float64)
    targets = np.asarray(target_separations, dtype=np.float64)
    if positions.ndim != 2:
        raise ValueError("initial_positions must have shape (realizations, particles)")
    if targets.ndim != 1 or targets.size == 0:
        raise ValueError("target_separations must be a nonempty vector")
    if cell_width <= 0.0 or realization_batch_size <= 0:
        raise ValueError("cell width and batch size must be positive")
    if np.any(targets - 0.5 * cell_width < 0.0) or np.any(
        targets + 0.5 * cell_width > box_length
    ):
        raise ValueError("every Ulam cell must lie inside the forward interval")

    realization_count, particle_count = positions.shape
    realization_parts: list[np.ndarray] = []
    first_label_parts: list[np.ndarray] = []
    second_label_parts: list[np.ndarray] = []
    separation_parts: list[np.ndarray] = []
    target_index_parts: list[np.ndarray] = []
    counts_by_target_and_realization = np.zeros(
        (targets.size, realization_count), dtype=np.int32
    )
    target_offsets = [0]

    diagonal = np.arange(particle_count)
    for target_index, target in enumerate(targets):
        lower = target - 0.5 * cell_width
        upper = target + 0.5 * cell_width
        target_pair_count = 0
        for start in range(0, realization_count, realization_batch_size):
            stop = min(start + realization_batch_size, realization_count)
            batch = positions[start:stop]
            separations = (batch[:, None, :] - batch[:, :, None]) % box_length
            selected = (separations >= lower) & (separations < upper)
            selected[:, diagonal, diagonal] = False
            local_realization, first_label, second_label = np.nonzero(selected)
            if local_realization.size == 0:
                continue
            realization = local_realization.astype(np.int32) + start
            pair_count = realization.size
            realization_parts.append(realization)
            first_label_parts.append(first_label.astype(np.int32))
            second_label_parts.append(second_label.astype(np.int32))
            separation_parts.append(
                separations[local_realization, first_label, second_label]
            )
            target_index_parts.append(
                np.full(pair_count, target_index, dtype=np.int16)
            )
            counts_by_target_and_realization[target_index] += np.bincount(
                realization, minlength=realization_count
            ).astype(np.int32)
            target_pair_count += pair_count
        target_offsets.append(target_offsets[-1] + target_pair_count)

    if realization_parts:
        realization_ids = np.concatenate(realization_parts)
        first_labels = np.concatenate(first_label_parts)
        second_labels = np.concatenate(second_label_parts)
        exact_initial_separations = np.concatenate(separation_parts)
        target_indices = np.concatenate(target_index_parts)
    else:
        realization_ids = np.empty(0, dtype=np.int32)
        first_labels = np.empty(0, dtype=np.int32)
        second_labels = np.empty(0, dtype=np.int32)
        exact_initial_separations = np.empty(0, dtype=np.float64)
        target_indices = np.empty(0, dtype=np.int16)

    return {
        "realization_ids": realization_ids,
        "first_labels": first_labels,
        "second_labels": second_labels,
        "exact_initial_forward_separations": exact_initial_separations,
        "target_indices": target_indices,
        "target_offsets": np.asarray(target_offsets, dtype=np.int64),
        "counts_by_target": np.sum(
            counts_by_target_and_realization, axis=1, dtype=np.int64
        ),
        "counts_by_target_and_realization": counts_by_target_and_realization,
    }


def measure_selected_signed_separations(
    positions: np.ndarray,
    pairs: dict[str, np.ndarray],
    box_length: float,
) -> np.ndarray:
    """Measure the same ordered labels selected at the initial time."""

    realization = pairs["realization_ids"]
    first = pairs["first_labels"]
    second = pairs["second_labels"]
    return signed_periodic_separation(
        positions[realization, first],
        positions[realization, second],
        box_length,
    )
