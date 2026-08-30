#!/usr/bin/env python3
"""Plot a candidate Cys-pair geometry in successful vs failed AF3 states."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import statistics

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--positive-seeds", default="1,2")
    ap.add_argument("--title", default="AF3-state-discriminating disulfide geometry")
    args = ap.parse_args()
    positive_seeds = {int(value) for value in args.positive_seeds.split(",")}
    values = {"successful C4 ring": [], "failed double-dimer": []}
    for row in csv.DictReader(args.table.open(), delimiter="\t"):
        seed = int(re.match(r"seed-(\d+)", row["model"]).group(1))
        key = "successful C4 ring" if seed in positive_seeds else "failed double-dimer"
        values[key].append(float(row["mean_cb_distance_a"]))
    labels = list(values)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    parts = ax.violinplot([values[label] for label in labels], showmedians=True, widths=.72)
    for body, color in zip(parts["bodies"], ("#2878b5", "#d95f4b")):
        body.set_facecolor(color); body.set_edgecolor("none"); body.set_alpha(.75)
    ax.axhline(6.0, color="#555555", ls="--", lw=1.2, label="screening geometry ≤6 Å")
    for index, label in enumerate(labels, 1):
        y = values[label]
        annotation_y = max(y) + (4.5 if max(y) < 10 else 1.0)
        ax.text(index, annotation_y, f"n={len(y)}\nmedian={statistics.median(y):.1f} Å",
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks([1, 2], labels)
    ax.set_ylabel("Best four-edge mean Cβ distance for 94→213 (Å)")
    ax.set_ylim(0, max(max(group) for group in values.values()) + 7)
    ax.set_title(args.title)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=.2)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    print(args.out)


if __name__ == "__main__":
    main()
