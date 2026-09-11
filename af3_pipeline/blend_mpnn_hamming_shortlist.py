#!/usr/bin/env python3
"""Blend a generation rank with AF3-labelled Hamming predictions."""

from __future__ import annotations

import argparse
import csv
import pathlib


def average_ranks(rows: list[dict[str, str]], field: str, *, higher: bool) -> list[float]:
    values = [float(row[field]) for row in rows]
    order = sorted(range(len(rows)), key=lambda i: values[i], reverse=higher)
    ranks = [0.0] * len(rows)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2.0
        for index in order[start:end]:
            ranks[index] = rank
        start = end
    return ranks


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, default=16)
    parser.add_argument("--minimum-distance", type=int, default=3)
    parser.add_argument("--fasta-prefix", default="E179-BLEND")
    args = parser.parse_args()

    with args.ranking.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if "af3_strict_surrogate_rank_score" in rows[0]:
        hamming_rank = average_ranks(rows, "af3_strict_surrogate_rank_score", higher=False)
    else:
        hamming_rank = average_ranks(rows, "robust_prediction", higher=True)
    mpnn_rank = average_ranks(rows, "generation_average_rank", higher=False)
    for row, hrank, mrank in zip(rows, hamming_rank, mpnn_rank):
        row["hamming_average_rank"] = hrank
        row["generation_borda_rank"] = mrank
        row["blended_average_rank"] = 0.5 * hrank + 0.5 * mrank
    rows.sort(key=lambda row: (
        float(row["blended_average_rank"]),
        -float(row["robust_prediction"]),
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
            if len(selected) == args.top:
                break
        if len(selected) == args.top:
            break
    if len(selected) != args.top:
        raise ValueError(f"selected only {len(selected)}/{args.top} sequences")
    for rank, row in enumerate(selected, 1):
        row["blended_rank"] = rank
    fields = ["blended_rank"] + [field for field in selected[0] if field != "blended_rank"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
    with args.out.with_suffix(".fasta").open("w") as handle:
        for row in selected:
            handle.write(
                f">{args.fasta_prefix}-{int(row['blended_rank']):02d} "
                f"nmut={row['n_mutations']} "
                f"blend_rank={float(row['blended_average_rank']):.3f} "
                f"pred_af3={float(row['predicted_af3_mean_iptm']):.3f} "
                f"pred_zero_clash={float(row.get('predicted_zero_clash_fraction', 1.0)):.3f}\n"
                f"{row['sequence']}\n"
            )
    print(f"input={len(rows)} selected={len(selected)} out={args.out}")


if __name__ == "__main__":
    main()
