#!/usr/bin/env python3
"""Compare learned position-weighted and uniform Hamming kernels by CV."""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict

import numpy as np

from crossmodel_surrogate_recombine import AA_INDEX, spearman
from score_c4_library_hamming_ai import aggregate_observations


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def learned_weights(x: np.ndarray, y: np.ndarray, scale: float) -> np.ndarray:
    if scale == 0:
        return np.ones(x.shape[1])
    total_variance = float(np.var(y)) + 1e-12
    global_mean = float(np.mean(y))
    explained = []
    for column in range(x.shape[1]):
        between = 0.0
        for allele in np.unique(x[:, column]):
            values = y[x[:, column] == allele]
            between += len(values) * (float(np.mean(values)) - global_mean) ** 2
        explained.append(between / (len(y) * total_variance))
    explained_array = np.asarray(explained)
    positive = explained_array[explained_array > 0]
    reference = float(np.median(positive)) if len(positive) else 1.0
    signal = np.clip(explained_array / (3.0 * reference + 1e-12), 0.0, 3.0)
    weights = 1.0 + scale * signal
    return weights / np.mean(weights)


def distances(left: np.ndarray, right: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return ((left[:, None, :] != right[None, :, :]) * weights).sum(axis=2)


def cross_validate(x: np.ndarray, y: np.ndarray, seed: int) -> list[tuple[float, float, float, float]]:
    rng = np.random.default_rng(seed)
    folds = np.arange(len(y)) % min(5, len(y))
    rng.shuffle(folds)
    grids = []
    for scale in (0.0, 0.5, 1.0, 2.0):
        predictions = {
            (gamma, ridge): np.zeros(len(y))
            for gamma in (0.06, 0.10, 0.16, 0.24)
            for ridge in (0.10, 0.30)
        }
        for fold in sorted(set(folds)):
            test = folds == fold
            train = ~test
            weights = learned_weights(x[train], y[train], scale)
            train_distance = distances(x[train], x[train], weights)
            test_distance = distances(x[test], x[train], weights)
            for gamma, ridge in predictions:
                kernel = np.exp(-gamma * train_distance)
                alpha = np.linalg.solve(
                    kernel + ridge * np.eye(int(train.sum())), y[train]
                )
                predictions[(gamma, ridge)][test] = np.exp(-gamma * test_distance) @ alpha
        for (gamma, ridge), predicted in predictions.items():
            grids.append((spearman(y, predicted), scale, gamma, ridge))
    return sorted(grids, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--seed", type=int, default=20260901)
    args = parser.parse_args()

    iptm_observations: dict[str, list[tuple[int, float]]] = defaultdict(list)
    clash_observations: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for path in args.ranking:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            iptm = (
                row.get("mean_iptm") or row.get("af3_mean_iptm")
                or row.get("median_iptm") or row.get("af3_median_iptm")
            )
            if len(sequence) != 386 or not iptm:
                continue
            n_models = int(row.get("n_models") or row.get("af3_n_models") or 5)
            iptm_observations[sequence].append((n_models, float(iptm)))
            zero_clash = (
                row.get("zero_atomic_clash_models")
                or row.get("af3_zero_atomic_clash_models")
                or row.get("zero_clash_models")
            )
            if zero_clash not in (None, ""):
                clash_observations[sequence].append(
                    (n_models, float(zero_clash) / n_models)
                )
            elif row.get("unclashed_fraction") not in (None, ""):
                clash_observations[sequence].append(
                    (n_models, float(row["unclashed_fraction"]))
                )

    outcomes = {
        "iptm": aggregate_observations(iptm_observations, "max_models_mean"),
        "zero_clash": aggregate_observations(clash_observations, "max_models_mean"),
    }
    positions = sorted({
        position
        for labels in outcomes.values()
        for sequence in labels
        for position in range(386)
        if any(other[position] != sequence[position] for other in labels)
    })
    for offset, (name, labels) in enumerate(outcomes.items()):
        sequences = list(labels)
        x = np.asarray([
            [AA_INDEX[sequence[position]] for position in positions]
            for sequence in sequences
        ], dtype=np.uint8)
        y = np.asarray([labels[sequence] for sequence in sequences], dtype=float)
        trials = cross_validate(x, y, args.seed + 7919 * offset)
        best = trials[0]
        uniform = max(trial for trial in trials if trial[1] == 0.0)
        print(
            f"{name}: labels={len(labels)} positions={len(positions)} "
            f"uniform_cv={uniform[0]:.3f} best_cv={best[0]:.3f} "
            f"scale={best[1]:.1f} gamma={best[2]:.2f} ridge={best[3]:.2f}"
        )


if __name__ == "__main__":
    main()
