#!/usr/bin/env python3
"""Re-rank constrained C4 proposals with a cross-validated epistatic ridge model.

The model is deliberately used only as a proposal ranker.  It encodes observed
mutation alleles plus recurrent allele pairs, evaluates them with held-out AF3
labels, and blends them with the existing Hamming-kernel prediction only when
the pairwise signal is independently predictive.  Every selected sequence must
still pass template-free AF3 and the frozen fail-closed screen.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import Counter, defaultdict

import numpy as np

from crossmodel_surrogate_recombine import average_ranks, spearman
from score_c4_library_hamming_ai import aggregate_observations


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def allele_sets(sequences: list[str], minimum_count: int) -> list[tuple[int, str]]:
    counts = Counter(
        (position, amino_acid)
        for sequence in sequences
        for position, amino_acid in enumerate(sequence)
    )
    variable = {
        position
        for position in range(len(sequences[0]))
        if len({sequence[position] for sequence in sequences}) > 1
    }
    return sorted(
        feature for feature, count in counts.items()
        if feature[0] in variable and count >= minimum_count
    )


def main_matrix(sequences: list[str], features: list[tuple[int, str]]) -> np.ndarray:
    return np.asarray([
        [float(sequence[position] == amino_acid) for position, amino_acid in features]
        for sequence in sequences
    ], dtype=np.float64)


def pair_columns(matrix: np.ndarray, features: list[tuple[int, str]], maximum: int) -> list[tuple[int, int]]:
    n = len(matrix)
    candidates: list[tuple[float, int, int]] = []
    for left in range(matrix.shape[1]):
        for right in range(left + 1, matrix.shape[1]):
            if features[left][0] == features[right][0]:
                continue
            count = int(np.sum(matrix[:, left] * matrix[:, right]))
            if 3 <= count <= n - 3:
                # Prefer recurrent, non-constant combinations.
                variance = (count / n) * (1.0 - count / n)
                candidates.append((-variance, left, right))
    candidates.sort()
    return [(left, right) for _, left, right in candidates[:maximum]]


def augment(matrix: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    if not pairs:
        return matrix
    pairwise = np.stack([matrix[:, left] * matrix[:, right] for left, right in pairs], axis=1)
    return np.concatenate([matrix, pairwise], axis=1)


def ridge_predict(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, ridge: float) -> np.ndarray:
    x_mean = train_x.mean(axis=0)
    y_mean = float(train_y.mean())
    centered = train_x - x_mean
    # Dual ridge is substantially cheaper when pair features outnumber labels.
    alpha = np.linalg.solve(centered @ centered.T + ridge * np.eye(len(train_x)), train_y - y_mean)
    beta = centered.T @ alpha
    return (test_x - x_mean) @ beta + y_mean


def random_folds(size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    folds = np.arange(size) % 5
    rng.shuffle(folds)
    return folds


def clustered_folds(sequences: list[str], radius: int = 3) -> tuple[np.ndarray, int]:
    """Keep close sequence neighbors in the same held-out fold."""
    representatives: list[str] = []
    clusters: list[list[int]] = []
    for index, sequence in sorted(enumerate(sequences), key=lambda item: item[1]):
        assigned = False
        for cluster_index, representative in enumerate(representatives):
            if hamming(sequence, representative) <= radius:
                clusters[cluster_index].append(index)
                assigned = True
                break
        if not assigned:
            representatives.append(sequence)
            clusters.append([index])
    folds = np.zeros(len(sequences), dtype=np.int32)
    fold_sizes = [0] * 5
    for cluster in sorted(clusters, key=lambda values: (-len(values), values[0])):
        fold = min(range(5), key=lambda value: (fold_sizes[value], value))
        folds[cluster] = fold
        fold_sizes[fold] += len(cluster)
    if min(fold_sizes) == 0:
        raise ValueError(f"clustered CV produced an empty fold: {fold_sizes}")
    return folds, len(clusters)


def cross_validate(x: np.ndarray, y: np.ndarray, folds: np.ndarray) -> tuple[float, float, np.ndarray]:
    trials: list[tuple[float, float, np.ndarray]] = []
    for ridge in (0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0):
        predicted = np.zeros(len(y), dtype=float)
        for fold in sorted(set(folds)):
            test = folds == fold
            train = ~test
            predicted[test] = ridge_predict(x[train], y[train], x[test], ridge)
        trials.append((spearman(y, predicted), ridge, predicted))
    trials.sort(key=lambda item: (-item[0], item[1]))
    return trials[0]


def kernel_cross_validate(sequences: list[str], y: np.ndarray, folds: np.ndarray) -> tuple[float, float, float]:
    variable_positions = [
        position for position in range(len(sequences[0]))
        if len({sequence[position] for sequence in sequences}) > 1
    ]
    encoded = np.asarray([
        [ord(sequence[position]) for position in variable_positions]
        for sequence in sequences
    ], dtype=np.uint8)
    distance = (encoded[:, None, :] != encoded[None, :, :]).sum(axis=2)
    trials: list[tuple[float, float, float]] = []
    for gamma in (0.06, 0.10, 0.16, 0.24, 0.36):
        for ridge in (0.03, 0.10, 0.30, 1.00):
            predicted = np.zeros(len(y), dtype=float)
            for fold in sorted(set(folds)):
                test = folds == fold
                train = ~test
                kernel = np.exp(-gamma * distance[np.ix_(train, train)])
                alpha = np.linalg.solve(kernel + ridge * np.eye(int(train.sum())), y[train])
                predicted[test] = np.exp(-gamma * distance[np.ix_(test, train)]) @ alpha
            trials.append((spearman(y, predicted), gamma, ridge))
    trials.sort(key=lambda item: (-item[0], item[1], item[2]))
    return trials[0]


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--library", type=pathlib.Path, required=True)
    parser.add_argument("--exclude-tsv", type=pathlib.Path, action="append", default=[])
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, default=32)
    parser.add_argument("--minimum-distance", type=int, default=3)
    parser.add_argument("--maximum-pairs", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument(
        "--library-manifold-only", action="store_true",
        help=(
            "fit the epistatic model only to AF3 labels that match the fixed "
            "background and allowed alleles of the finite input library"
        ),
    )
    args = parser.parse_args()

    observations: dict[str, list[tuple[int, float]]] = defaultdict(list)
    for path in args.ranking:
        for row in read_tsv(path):
            sequence = row.get("sequence", "")
            score = (
                row.get("mean_iptm") or row.get("af3_mean_iptm")
                or row.get("median_iptm") or row.get("af3_median_iptm")
            )
            if len(sequence) != 386 or score in (None, ""):
                continue
            n_models = int(row.get("n_models") or row.get("af3_n_models") or 5)
            observations[sequence].append((n_models, float(score)))
    labels = aggregate_observations(observations, "max_models_mean")

    library = [row for row in read_tsv(args.library) if len(row.get("sequence", "")) == 386]
    if args.library_manifold_only:
        library_sequences_all = [row["sequence"] for row in library]
        variable_positions = [
            position for position in range(386)
            if len({sequence[position] for sequence in library_sequences_all}) > 1
        ]
        variable_set = set(variable_positions)
        template = library_sequences_all[0]
        allowed = {
            position: {sequence[position] for sequence in library_sequences_all}
            for position in variable_positions
        }

        def on_library_manifold(sequence: str) -> bool:
            return all(
                (sequence[position] in allowed[position]) if position in variable_set
                else (sequence[position] == template[position])
                for position in range(386)
            )

        labels = {
            sequence: score for sequence, score in labels.items()
            if on_library_manifold(sequence)
        }
    minimum_labels = 25 if args.library_manifold_only else 50
    if len(labels) < minimum_labels:
        raise ValueError(
            f"expected at least {minimum_labels} AF3-labelled sequences, found {len(labels)}"
        )
    excluded = set(labels)
    for path in args.exclude_tsv:
        excluded.update(
            row.get("sequence", "") for row in read_tsv(path)
            if len(row.get("sequence", "")) == 386
        )
    library = [row for row in library if row["sequence"] not in excluded]
    if len(library) < args.top:
        raise ValueError(f"only {len(library)} novel library sequences")

    train_sequences = list(labels)
    features = allele_sets(train_sequences, minimum_count=2)
    train_main = main_matrix(train_sequences, features)
    pairs = pair_columns(train_main, features, args.maximum_pairs)
    train_x = augment(train_main, pairs)
    train_y = np.asarray([labels[sequence] for sequence in train_sequences], dtype=float)
    random_assignment = random_folds(len(train_y), args.seed)
    clustered_assignment, n_clusters = clustered_folds(train_sequences, radius=3)
    random_cv, random_ridge, _ = cross_validate(train_x, train_y, random_assignment)
    cv, ridge, cv_prediction = cross_validate(train_x, train_y, clustered_assignment)
    kernel_random_cv, kernel_random_gamma, kernel_random_ridge = kernel_cross_validate(
        train_sequences, train_y, random_assignment
    )
    kernel_clustered_cv, kernel_clustered_gamma, kernel_clustered_ridge = kernel_cross_validate(
        train_sequences, train_y, clustered_assignment
    )

    library_sequences = [row["sequence"] for row in library]
    library_x = augment(main_matrix(library_sequences, features), pairs)
    # Five leave-one-fold-out fits provide an uncertainty estimate for proposals.
    fold_predictions = []
    for fold in range(5):
        train = clustered_assignment != fold
        fold_predictions.append(ridge_predict(train_x[train], train_y[train], library_x, ridge))
    fold_predictions_array = np.stack(fold_predictions)
    epistasis_mean = fold_predictions_array.mean(axis=0)
    epistasis_sd = fold_predictions_array.std(axis=0)
    robust_epistasis = epistasis_mean - 0.20 * epistasis_sd

    def kernel_prediction(row: dict[str, str]) -> float:
        for field in ("predicted_af3", "predicted_af3_mean_iptm", "robust_prediction"):
            value = row.get(field)
            if value not in (None, ""):
                return float(value)
        raise ValueError(
            "library is missing a Hamming-kernel prediction column; expected "
            "predicted_af3, predicted_af3_mean_iptm, or robust_prediction"
        )

    kernel = np.asarray([kernel_prediction(row) for row in library], dtype=float)
    if cv >= 0.50:
        epistasis_weight = 0.40
    elif cv >= 0.20:
        epistasis_weight = 0.25
    else:
        epistasis_weight = 0.0
    selection_score = (
        (1.0 - epistasis_weight) * average_ranks(kernel, higher=True)
        + epistasis_weight * average_ranks(robust_epistasis, higher=True)
    )
    order = np.argsort(selection_score)
    selected: list[int] = []
    for minimum in range(args.minimum_distance, 0, -1):
        for index in order:
            index = int(index)
            if index in selected:
                continue
            sequence = library_sequences[index]
            if all(hamming(sequence, library_sequences[old]) >= minimum for old in selected):
                selected.append(index)
            if len(selected) == args.top:
                break
        if len(selected) == args.top:
            break

    args.out.mkdir(parents=True, exist_ok=True)
    selected_rank = {index: rank for rank, index in enumerate(selected, 1)}
    all_rows: list[dict[str, object]] = []
    for index, source_row in enumerate(library):
        row = dict(source_row)
        row.update({
            "epistasis_selected_rank": selected_rank.get(index, 0),
            "epistasis_cv_spearman": cv,
            "epistasis_prediction": epistasis_mean[index],
            "epistasis_prediction_sd": epistasis_sd[index],
            "epistasis_robust_prediction": robust_epistasis[index],
            "epistasis_weight": epistasis_weight,
            "blended_selection_average_rank": selection_score[index],
        })
        all_rows.append(row)
    all_fields: list[str] = []
    for row in all_rows:
        all_fields.extend(field for field in row if field not in all_fields)
    with (args.out / "epistasis_all_library.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, all_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(all_rows)

    rows: list[dict[str, object]] = []
    for rank, index in enumerate(selected, 1):
        source = dict(library[index])
        source.update({
            "epistasis_rank": rank,
            "epistasis_cv_spearman": cv,
            "epistasis_ridge": ridge,
            "epistasis_main_features": len(features),
            "epistasis_pair_features": len(pairs),
            "epistasis_prediction": epistasis_mean[index],
            "epistasis_prediction_sd": epistasis_sd[index],
            "epistasis_used": int(epistasis_weight > 0),
            "epistasis_weight": epistasis_weight,
            "blended_selection_average_rank": selection_score[index],
        })
        rows.append(source)
    fields: list[str] = []
    for row in rows:
        fields.extend(field for field in row if field not in fields)
    with (args.out / "epistasis_reranked_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "epistasis_reranked_candidates.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">EPI-{int(row['epistasis_rank']):03d} source={row['name']} "
                f"nmut={row['n_mutations']} kernel={kernel_prediction(row):.4f} "
                f"epistasis={float(row['epistasis_prediction']):.4f}\n{row['sequence']}\n"
            )
    with (args.out / "epistasis_cv.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "n_labels", "main_features", "pair_features", "cluster_radius",
            "n_sequence_clusters", "clustered_ridge", "clustered_cv_spearman",
            "random_ridge", "random_cv_spearman", "prediction_spearman_check",
            "kernel_clustered_cv_spearman", "kernel_clustered_gamma", "kernel_clustered_ridge",
            "kernel_random_cv_spearman", "kernel_random_gamma", "kernel_random_ridge",
            "epistasis_weight", "used",
        ])
        writer.writerow([
            len(labels), len(features), len(pairs), 3, n_clusters, ridge, cv,
            random_ridge, random_cv, spearman(train_y, cv_prediction),
            kernel_clustered_cv, kernel_clustered_gamma, kernel_clustered_ridge,
            kernel_random_cv, kernel_random_gamma, kernel_random_ridge,
            epistasis_weight, int(epistasis_weight > 0),
        ])
    with (args.out / "epistasis_label_cv_predictions.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sequence", "observed_af3_iptm", "heldout_epistasis_prediction"])
        for sequence, observed, predicted in zip(train_sequences, train_y, cv_prediction):
            writer.writerow([sequence, observed, predicted])
    print(
        f"labels={len(labels)} library={len(library)} main={len(features)} pairs={len(pairs)} "
        f"epistasis_clustered_CV={cv:.3f} random_CV={random_cv:.3f} "
        f"kernel_clustered_CV={kernel_clustered_cv:.3f} random_CV={kernel_random_cv:.3f} "
        f"clusters={n_clusters} manifold_only={int(args.library_manifold_only)} "
        f"ridge={ridge:g} epistasis_weight={epistasis_weight:.2f} selected={len(rows)}"
    )


if __name__ == "__main__":
    main()
