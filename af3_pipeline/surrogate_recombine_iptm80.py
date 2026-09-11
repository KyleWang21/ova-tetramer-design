#!/usr/bin/env python3
"""Rank local C4 sequence recombinations with a small AF3-trained kernel ensemble.

The model is deliberately only a proposal generator.  Labels are completed
template-free AF3 screens; every selected sequence is still re-run through AF3.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib

import numpy as np


NATIVE_CYS = {12, 31, 74, 121, 368, 383}
EPITOPE = range(258, 266)


def read_fasta(path: pathlib.Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3:
        return float("nan")
    return float(np.corrcoef(rankdata(left), rankdata(right))[0, 1])


def hamming(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left[:, None, :] != right[None, :, :]).sum(axis=2)


def kernel_predict(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    gamma: float,
    ridge: float,
) -> np.ndarray:
    train_kernel = np.exp(-gamma * hamming(train_x, train_x))
    alpha = np.linalg.solve(train_kernel + ridge * np.eye(len(train_x)), train_y)
    return np.exp(-gamma * hamming(test_x, train_x)) @ alpha


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", type=pathlib.Path, required=True)
    ap.add_argument("--ranking", type=pathlib.Path, action="append", required=True)
    ap.add_argument("--anchors", type=pathlib.Path, required=True)
    ap.add_argument("--exclude-fasta", type=pathlib.Path, action="append", default=[])
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--proposals", type=int, default=100000)
    ap.add_argument("--top", type=int, default=96)
    ap.add_argument("--seed", type=int, default=20260831)
    args = ap.parse_args()

    reference = read_fasta(args.reference)[0][1]
    anchors = read_fasta(args.anchors)
    if len(reference) != 386 or any(len(seq) != 386 for _, seq in anchors):
        raise SystemExit("reference and anchors must be 386 aa")

    labelled: dict[str, float] = {}
    for path in args.ranking:
        for row in csv.DictReader(path.open(), delimiter="\t"):
            sequence = row.get("sequence", "")
            raw_score = row.get("mean_iptm") or row.get("median_iptm")
            if len(sequence) != 386 or not raw_score:
                continue
            score = float(raw_score)
            labelled[sequence] = max(score, labelled.get(sequence, -math.inf))
    if len(labelled) < 20:
        raise SystemExit(f"too few AF3-labelled sequences: {len(labelled)}")

    excluded = set(labelled)
    for path in args.exclude_fasta:
        excluded.update(sequence for _, sequence in read_fasta(path))

    labelled_items = sorted(labelled.items(), key=lambda item: item[1], reverse=True)
    high = labelled_items[: min(24, len(labelled_items))]
    positions = sorted({
        index for sequence, _ in labelled_items for index, (old, new) in
        enumerate(zip(reference, sequence)) if old != new
    })
    alphabet: dict[int, list[tuple[str, float]]] = {}
    for pos in positions:
        weights: dict[str, float] = {reference[pos]: 0.15}
        for sequence, score in high:
            aa = sequence[pos]
            weights[aa] = weights.get(aa, 0.0) + math.exp(8.0 * (score - 0.70))
        for _, sequence in anchors:
            weights[sequence[pos]] = weights.get(sequence[pos], 0.0) + 3.0
        total = sum(weights.values())
        alphabet[pos] = [(aa, weight / total) for aa, weight in weights.items()]

    rng = np.random.default_rng(args.seed)
    proposals: set[str] = set()
    anchor_sequences = [sequence for _, sequence in anchors]
    top_sequences = [sequence for sequence, _ in high]
    attempts = 0
    max_attempts = max(args.proposals * 100, 100000)
    while len(proposals) < args.proposals and attempts < max_attempts:
        attempts += 1
        if rng.random() < 0.45:
            # Crossover preserves multi-position epistasis observed in an AF3 hit.
            left = list(anchor_sequences[int(rng.integers(len(anchor_sequences)))])
            donor = top_sequences[int(rng.integers(len(top_sequences)))]
            chosen = rng.choice(positions, size=int(rng.integers(2, 8)), replace=False)
            for pos in chosen:
                left[pos] = donor[pos]
        else:
            left = list(anchor_sequences[int(rng.integers(len(anchor_sequences)))])
            chosen = rng.choice(positions, size=int(rng.integers(1, 7)), replace=False)
            for pos in chosen:
                options, probabilities = zip(*alphabet[pos])
                left[pos] = str(rng.choice(options, p=probabilities))
        sequence = "".join(left)
        nmut = sum(a != b for a, b in zip(reference, sequence))
        if not 16 <= nmut < 25 or sequence in excluded:
            continue
        if {i for i, aa in enumerate(sequence, 1) if aa == "C"} != NATIVE_CYS:
            continue
        if any(sequence[i - 1] != reference[i - 1] for i in EPITOPE):
            continue
        proposals.add(sequence)
    if len(proposals) < args.top:
        raise SystemExit(
            f"proposal space exhausted after {attempts} attempts: {len(proposals)} valid"
        )

    aa_order = "ACDEFGHIKLMNPQRSTVWY"
    aa_index = {aa: index for index, aa in enumerate(aa_order)}
    train_sequences = list(labelled)
    train_y = np.asarray([labelled[sequence] for sequence in train_sequences])
    train_x = np.asarray([[aa_index[s[pos]] for pos in positions] for s in train_sequences])
    proposal_sequences = list(proposals)
    proposal_x = np.asarray([[aa_index[s[pos]] for pos in positions] for s in proposal_sequences])

    # Five deterministic folds choose a small Hamming-kernel ridge ensemble.
    folds = np.arange(len(train_x)) % 5
    rng.shuffle(folds)
    trials: list[tuple[float, float, float]] = []
    for gamma in (0.08, 0.12, 0.18, 0.26, 0.38):
        for ridge in (0.01, 0.03, 0.10, 0.30):
            predicted = np.zeros(len(train_y))
            for fold in range(5):
                test = folds == fold
                fit = ~test
                predicted[test] = kernel_predict(
                    train_x[fit], train_y[fit], train_x[test], gamma, ridge
                )
            trials.append((spearman(train_y, predicted), gamma, ridge))
    trials.sort(reverse=True)
    ensemble = trials[:3]
    predictions = np.stack([
        kernel_predict(train_x, train_y, proposal_x, gamma, ridge)
        for _, gamma, ridge in ensemble
    ])
    mean_prediction = predictions.mean(axis=0)
    uncertainty = predictions.std(axis=0)
    order = np.argsort(-(mean_prediction - 0.15 * uncertainty))

    selected: list[int] = []
    for minimum_distance in (3, 2, 1):
        for idx in order:
            if idx in selected:
                continue
            sequence = proposal_sequences[idx]
            if all(sum(a != b for a, b in zip(sequence, proposal_sequences[old])) >= minimum_distance for old in selected):
                selected.append(int(idx))
            if len(selected) >= args.top:
                break
        if len(selected) >= args.top:
            break

    args.out.mkdir(parents=True, exist_ok=True)
    fields = ["rank", "name", "surrogate_rank_score", "ensemble_sd", "n_mutations", "mutations", "nearest_anchor_distance", "sequence"]
    rows: list[dict[str, object]] = []
    for rank, idx in enumerate(selected, 1):
        sequence = proposal_sequences[idx]
        rows.append({
            "rank": rank,
            "name": f"OVA-C4-SURR-{rank:03d}",
            "surrogate_rank_score": float(mean_prediction[idx]),
            "ensemble_sd": float(uncertainty[idx]),
            "n_mutations": sum(a != b for a, b in zip(reference, sequence)),
            "mutations": ",".join(f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b),
            "nearest_anchor_distance": min(sum(a != b for a, b in zip(sequence, anchor)) for anchor in anchor_sequences),
            "sequence": sequence,
        })
    with (args.out / "surrogate_candidates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    with (args.out / "surrogate_candidates.fasta").open("w") as handle:
        for row in rows:
            handle.write(
                f">{row['name']} surrogate_score={row['surrogate_rank_score']:.4f} "
                f"nmut={row['n_mutations']} nearest_anchor={row['nearest_anchor_distance']}\n"
                f"{row['sequence']}\n"
            )
    with (args.out / "cv_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rank", "cv_spearman", "gamma", "ridge"])
        for rank, (correlation, gamma, ridge) in enumerate(trials, 1):
            writer.writerow([rank, correlation, gamma, ridge])
    print(
        f"labels={len(labelled)} positions={len(positions)} proposals={len(proposals)} "
        f"attempts={attempts} "
        f"selected={len(rows)} best_cv_spearman={trials[0][0]:.3f} "
        f"surrogate_score_range={mean_prediction[selected].min():.3f}-{mean_prediction[selected].max():.3f}"
    )


if __name__ == "__main__":
    main()
