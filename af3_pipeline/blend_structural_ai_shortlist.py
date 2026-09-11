#!/usr/bin/env python3
"""Blend frozen generation and AF3/clash surrogate ranks with sequence diversity."""

from __future__ import annotations

import argparse
import csv
import pathlib

import numpy as np


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--shortlist", type=int, default=16)
    parser.add_argument("--minimum-distance", type=int, default=3)
    parser.add_argument("--generation-weight", type=float, default=0.5)
    parser.add_argument("--ai-weight", type=float, default=0.5)
    args = parser.parse_args()

    rows = read_tsv(args.ranking)
    if not rows:
        raise ValueError("empty ranking")
    if args.generation_weight < 0 or args.ai_weight < 0:
        raise ValueError("rank weights must be nonnegative")
    total_weight = args.generation_weight + args.ai_weight
    if total_weight <= 0:
        raise ValueError("at least one rank weight must be positive")

    generation = np.asarray(
        [float(row["generation_average_rank"]) for row in rows], dtype=float
    )
    strict_ai = np.asarray(
        [float(row["af3_strict_surrogate_rank_score"]) for row in rows], dtype=float
    )
    generation_rank = average_ranks(generation)
    ai_rank = average_ranks(strict_ai)
    blended = (
        args.generation_weight * generation_rank + args.ai_weight * ai_rank
    ) / total_weight
    for row, grank, arank, score in zip(rows, generation_rank, ai_rank, blended):
        row["generation_borda_rank"] = float(grank)
        row["strict_ai_borda_rank"] = float(arank)
        row["blended_average_rank"] = float(score)

    rows.sort(key=lambda row: (
        float(row["blended_average_rank"]),
        float(row["af3_strict_surrogate_rank_score"]),
        float(row["generation_average_rank"]),
        int(row["n_mutations"]),
        row["sequence"],
    ))
    selected: list[dict[str, str]] = []
    for minimum in (args.minimum_distance, 1):
        for row in rows:
            if row in selected:
                continue
            if all(hamming(row["sequence"], old["sequence"]) >= minimum for old in selected):
                selected.append(row)
            if len(selected) == args.shortlist:
                break
        if len(selected) == args.shortlist:
            break
    if len(selected) != args.shortlist:
        raise ValueError(
            f"only {len(selected)} diverse sequences available; requested {args.shortlist}"
        )
    for rank, row in enumerate(selected, 1):
        row["blended_rank"] = rank

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["blended_rank"] + [field for field in selected[0] if field != "blended_rank"]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    with args.out.with_suffix(".fasta").open("w") as handle:
        for row in selected:
            handle.write(
                f">{row['name']} blended_rank={row['blended_rank']} "
                f"nmut={row['n_mutations']} predicted_af3="
                f"{float(row['predicted_af3_mean_iptm']):.3f} "
                f"predicted_zero_clash="
                f"{float(row['predicted_zero_clash_fraction']):.3f}\n"
                f"{row['sequence']}\n"
            )
    print(
        f"library={len(rows)} selected={len(selected)} "
        f"minimum_distance={args.minimum_distance} out={args.out}"
    )


if __name__ == "__main__":
    main()
