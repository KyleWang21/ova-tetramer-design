#!/usr/bin/env python3
"""Rank a fixed sequence library with the AF3-labelled Hamming-kernel ensemble."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from crossmodel_surrogate_recombine import (
    AA_INDEX, ensemble_predict, fit_ensemble, spearman,
)


def learned_position_weights(
    train_x: np.ndarray, train_y: np.ndarray, scale: float
) -> np.ndarray:
    if scale == 0:
        return np.ones(train_x.shape[1])
    total_variance = float(np.var(train_y)) + 1e-12
    global_mean = float(np.mean(train_y))
    explained = []
    for column in range(train_x.shape[1]):
        between = 0.0
        for allele in np.unique(train_x[:, column]):
            values = train_y[train_x[:, column] == allele]
            between += len(values) * (float(np.mean(values)) - global_mean) ** 2
        explained.append(between / (len(train_y) * total_variance))
    explained_array = np.asarray(explained)
    positive = explained_array[explained_array > 0]
    reference = float(np.median(positive)) if len(positive) else 1.0
    signal = np.clip(explained_array / (3.0 * reference + 1e-12), 0.0, 3.0)
    weights = 1.0 + scale * signal
    return weights / np.mean(weights)


def weighted_hamming(
    left: np.ndarray, right: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    return ((left[:, None, :] != right[None, :, :]) * weights).sum(axis=2)


@dataclass
class WeightedKernelModel:
    gamma: float
    ridge: float
    train_x: np.ndarray
    weights: np.ndarray
    alpha: np.ndarray

    def predict(self, test_x: np.ndarray, batch: int = 4096) -> np.ndarray:
        output = []
        for start in range(0, len(test_x), batch):
            distance = weighted_hamming(
                test_x[start:start + batch], self.train_x, self.weights
            )
            output.append(np.exp(-self.gamma * distance) @ self.alpha)
        return np.concatenate(output)


def fit_weighted_model(
    train_x: np.ndarray, train_y: np.ndarray, scale: float, gamma: float, ridge: float
) -> WeightedKernelModel:
    weights = learned_position_weights(train_x, train_y, scale)
    kernel = np.exp(-gamma * weighted_hamming(train_x, train_x, weights))
    alpha = np.linalg.solve(kernel + ridge * np.eye(len(train_x)), train_y)
    return WeightedKernelModel(gamma, ridge, train_x, weights, alpha)


def fit_weighted_ensemble(
    train_x: np.ndarray, train_y: np.ndarray, seed: int
) -> tuple[list[WeightedKernelModel], float, list[tuple[float, float, float, float]]]:
    rng = np.random.default_rng(seed)
    folds = np.arange(len(train_x)) % min(5, len(train_x))
    rng.shuffle(folds)
    trials = []
    for scale in (0.0, 0.5, 1.0, 2.0):
        predictions = {
            (gamma, ridge): np.zeros(len(train_y))
            for gamma in (0.06, 0.10, 0.16, 0.24, 0.36)
            for ridge in (0.03, 0.10, 0.30, 1.00)
        }
        for fold in sorted(set(folds)):
            test = folds == fold
            fit = ~test
            weights = learned_position_weights(train_x[fit], train_y[fit], scale)
            train_distance = weighted_hamming(train_x[fit], train_x[fit], weights)
            test_distance = weighted_hamming(train_x[test], train_x[fit], weights)
            for gamma, ridge in predictions:
                kernel = np.exp(-gamma * train_distance)
                alpha = np.linalg.solve(
                    kernel + ridge * np.eye(int(fit.sum())), train_y[fit]
                )
                predictions[(gamma, ridge)][test] = (
                    np.exp(-gamma * test_distance) @ alpha
                )
        for (gamma, ridge), predicted in predictions.items():
            trials.append((spearman(train_y, predicted), scale, gamma, ridge))
    trials.sort(key=lambda item: (
        -(-math.inf if math.isnan(item[0]) else item[0]), item[1], item[2], item[3]
    ))
    chosen = trials[:3]
    return (
        [
            fit_weighted_model(train_x, train_y, scale, gamma, ridge)
            for _, scale, gamma, ridge in chosen
        ],
        chosen[0][0],
        trials,
    )


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: object, default: float = math.nan) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def average_ranks(values: np.ndarray, *, higher: bool) -> np.ndarray:
    order = np.argsort(-values if higher else values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def aggregate_observations(
    observations: dict[str, list[tuple[int, float]]], policy: str
) -> dict[str, float]:
    labels = {}
    for sequence, values in observations.items():
        if policy == "max":
            labels[sequence] = max(score for _, score in values)
            continue
        max_models = max(n_models for n_models, _ in values)
        retained = [score for n_models, score in values if n_models == max_models]
        labels[sequence] = float(np.mean(retained))
    return labels


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--library", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--shortlist", type=int, default=0)
    parser.add_argument(
        "--allow-short-shortlist", action="store_true",
        help="write every distance-qualified sequence when fewer than --shortlist exist",
    )
    parser.add_argument("--minimum-distance", type=int, default=3)
    parser.add_argument(
        "--generation-rank-weight",
        type=float,
        default=0.0,
        help=(
            "Weight assigned to the library's generation_average_rank when "
            "forming the final selection rank.  The remaining weight is "
            "assigned to the 70:30 AF3-ipTM/zero-clash surrogate rank."
        ),
    )
    parser.add_argument(
        "--weighted-iptm",
        action="store_true",
        help="Learn position weights inside each CV fold for the ipTM kernel only.",
    )
    parser.add_argument(
        "--duplicate-label-policy",
        choices=("max_models_mean", "max"),
        default="max_models_mean",
        help=(
            "Resolve repeated AF3 labels for one sequence. max_models_mean "
            "prefers the observation(s) with the largest model count and "
            "averages ties; max reproduces the legacy optimistic behavior."
        ),
    )
    args = parser.parse_args()
    if not 0.0 <= args.generation_rank_weight <= 1.0:
        raise ValueError("--generation-rank-weight must be between 0 and 1")

    observations: dict[str, list[tuple[int, float]]] = defaultdict(list)
    clash_observations: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for path in args.ranking:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            score = (
                row.get("mean_iptm") or row.get("af3_mean_iptm")
                or row.get("median_iptm") or row.get("af3_median_iptm")
            )
            if len(sequence) == 386 and score:
                # Full multi-seed summaries expose n_models.  The strict
                # seed1 summaries are always five samples and lack that
                # column, so five is the correct coverage fallback.
                try:
                    n_models = int(row.get("n_models") or row.get("af3_n_models") or 5)
                except ValueError:
                    n_models = 5
                observations[sequence].append((n_models, float(score)))
                zero_clash = (
                    row.get("zero_atomic_clash_models")
                    or row.get("af3_zero_atomic_clash_models")
                    or row.get("zero_clash_models")
                    or row.get("unclashed_models")
                )
                if zero_clash not in (None, ""):
                    clash_observations[sequence].append(
                        (n_models, float(zero_clash) / n_models)
                    )
                elif row.get("unclashed_fraction") not in (None, ""):
                    clash_observations[sequence].append(
                        (n_models, float(row["unclashed_fraction"]))
                    )
    labels = aggregate_observations(observations, args.duplicate_label_policy)
    clash_labels = aggregate_observations(clash_observations, args.duplicate_label_policy)
    if len(labels) < 50:
        raise ValueError(f"too few AF3-labelled sequences: {len(labels)}")
    library = read_tsv(args.library)
    positions = sorted({
        index for sequence in list(labels) + [row["sequence"] for row in library]
        for index, aa in enumerate(sequence) if any(old[index] != aa for old in labels)
    })
    train_sequences = list(labels)
    train_x = np.asarray([
        [AA_INDEX[sequence[position]] for position in positions] for sequence in train_sequences
    ], dtype=np.uint8)
    train_y = np.asarray([labels[sequence] for sequence in train_sequences], dtype=float)
    test_x = np.asarray([
        [AA_INDEX[row["sequence"][position]] for position in positions] for row in library
    ], dtype=np.uint8)
    if args.weighted_iptm:
        models, cv, _ = fit_weighted_ensemble(train_x, train_y, args.seed)
        iptm_kernel_mode = "learned_position_weighted_hamming"
    else:
        models, cv, _ = fit_ensemble(train_x, train_y, args.seed)
        iptm_kernel_mode = "uniform_hamming"
    predicted, uncertainty = ensemble_predict(models, test_x)
    robust_iptm = predicted - 0.20 * uncertainty
    clash_cv = float("nan")
    predicted_clash = np.ones(len(library), dtype=float)
    clash_uncertainty = np.zeros(len(library), dtype=float)
    if len(clash_labels) >= 30:
        clash_sequences = list(clash_labels)
        clash_x = np.asarray([
            [AA_INDEX[sequence[position]] for position in positions]
            for sequence in clash_sequences
        ], dtype=np.uint8)
        clash_y = np.asarray([clash_labels[sequence] for sequence in clash_sequences], dtype=float)
        clash_models, clash_cv, _ = fit_ensemble(clash_x, clash_y, args.seed + 7919)
        predicted_clash, clash_uncertainty = ensemble_predict(clash_models, test_x)
        predicted_clash = np.clip(predicted_clash, 0.0, 1.0)
    robust_clash = np.clip(predicted_clash - 0.20 * clash_uncertainty, 0.0, 1.0)
    iptm_rank = average_ranks(robust_iptm, higher=True)
    clash_rank = average_ranks(robust_clash, higher=True)
    strict_rank_score = 0.70 * iptm_rank + 0.30 * clash_rank
    generation_values = np.asarray([
        number(row.get("generation_average_rank"), math.inf) for row in library
    ])
    if args.generation_rank_weight:
        if not np.isfinite(generation_values).all():
            raise ValueError(
                "--generation-rank-weight requires finite generation_average_rank "
                "for every library sequence"
            )
        generation_rank = average_ranks(generation_values, higher=False)
    else:
        generation_rank = np.full(len(library), math.nan)
    selection_rank_score = (
        (1.0 - args.generation_rank_weight) * strict_rank_score
        + args.generation_rank_weight
        * np.nan_to_num(generation_rank, nan=0.0)
    )
    for row, mean, sd in zip(library, predicted, uncertainty):
        row["predicted_af3_mean_iptm"] = float(mean)
        row["prediction_sd"] = float(sd)
        row["robust_prediction"] = float(mean - 0.20 * sd)
        row["af3_surrogate_cv_spearman"] = cv
        row["af3_iptm_kernel_mode"] = iptm_kernel_mode
    for row, mean, sd, robust, rank_score, generation_order, selection_score in zip(
        library, predicted_clash, clash_uncertainty, robust_clash, strict_rank_score,
        generation_rank, selection_rank_score,
    ):
        row["predicted_zero_clash_fraction"] = float(mean)
        row["zero_clash_prediction_sd"] = float(sd)
        row["robust_zero_clash_fraction"] = float(robust)
        row["zero_clash_surrogate_cv_spearman"] = clash_cv
        row["af3_strict_surrogate_rank_score"] = float(rank_score)
        row["generation_prior_rank"] = float(generation_order)
        row["selection_rank_score"] = float(selection_score)
    library.sort(key=lambda row: (
        float(row["selection_rank_score"]),
        -float(row["robust_prediction"]), int(row["n_mutations"]), row["sequence"]
    ))
    for rank, row in enumerate(library, 1):
        row["ai_rank"] = rank
    fields = ["ai_rank"] + [field for field in library[0] if field != "ai_rank"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(library)
    with args.out.with_suffix(".fasta").open("w") as handle:
        for row in library:
            edit = row.get("reverted_positions") or row.get("crossover_b_positions", "")
            handle.write(
                f">{row['name']} ai_rank={row['ai_rank']} nmut={row['n_mutations']} "
                f"predicted_af3={float(row['predicted_af3_mean_iptm']):.3f} "
                f"predicted_zero_clash={float(row['predicted_zero_clash_fraction']):.3f} "
                f"edits={edit}\n{row['sequence']}\n"
            )
    if args.shortlist:
        selected: list[dict[str, str]] = []
        for row in library:
            if all(
                sum(a != b for a, b in zip(row["sequence"], old["sequence"]))
                >= args.minimum_distance
                for old in selected
            ):
                selected.append(row)
            if len(selected) >= args.shortlist:
                break
        if len(selected) != args.shortlist and not args.allow_short_shortlist:
            raise ValueError(
                f"only {len(selected)} sequences satisfy minimum distance "
                f"{args.minimum_distance}; requested {args.shortlist}"
            )
        if not selected:
            raise ValueError("no sequences satisfy shortlist distance constraint")
        short_tsv = args.out.with_name(args.out.stem + "_shortlist.tsv")
        with short_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(selected)
        with short_tsv.with_suffix(".fasta").open("w") as handle:
            for rank, row in enumerate(selected, 1):
                handle.write(
                    f">AI-MULTIHIT-{rank:02d} source_ai_rank={row['ai_rank']} "
                    f"nmut={row['n_mutations']} predicted_af3="
                    f"{float(row['predicted_af3_mean_iptm']):.3f} "
                    f"predicted_zero_clash={float(row['predicted_zero_clash_fraction']):.3f}"
                    f"\n{row['sequence']}\n"
                )
    print(
        f"label_observations={sum(map(len, observations.values()))} "
        f"labels={len(labels)} duplicate_sequences="
        f"{sum(len(values) > 1 for values in observations.values())} "
        f"duplicate_policy={args.duplicate_label_policy} "
        f"clash_labels={len(clash_labels)} library={len(library)} "
        f"iptm_kernel={iptm_kernel_mode} cv_spearman={cv:.3f} "
        f"clash_cv_spearman={clash_cv:.3f} "
        f"generation_rank_weight={args.generation_rank_weight:.3f} "
        f"predicted_range={predicted.min():.3f}-{predicted.max():.3f}"
    )
    for row in library[:10]:
        edit = row.get("reverted_positions") or row.get("crossover_b_positions", "")
        print(
            row["ai_rank"], row["name"], row["n_mutations"], edit,
            f"pred={float(row['predicted_af3_mean_iptm']):.3f}",
            f"zero_clash={float(row['predicted_zero_clash_fraction']):.3f}",
        )


if __name__ == "__main__":
    main()
