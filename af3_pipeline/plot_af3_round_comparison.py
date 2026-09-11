#!/usr/bin/env python3
"""Compare strict AF3 seed-1 distributions across active-learning rounds."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--round", action="append", required=True, metavar="LABEL=TSV",
        help="round label and seed1_strict_summary.tsv path",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rounds: list[tuple[str, list[dict[str, str]]]] = []
    for specification in args.round:
        label, raw_path = specification.split("=", 1)
        with Path(raw_path).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            raise ValueError(f"empty AF3 summary: {raw_path}")
        rounds.append((label, rows))

    colors = ["#2f6f9f", "#d1495b", "#2a9d8f", "#8a5fbf"]
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for index, (label, rows) in enumerate(rounds):
        values = np.asarray([float(row["mean_iptm"]) for row in rows])
        passed = np.asarray([bool(int(row["seed1_strict_pass"])) for row in rows])
        order = np.argsort(values)[::-1]
        x = np.arange(1, len(values) + 1)
        color = colors[index % len(colors)]
        axes[0].plot(x, values[order], marker="o", markersize=3.5, lw=1.3,
                     color=color, label=f"{label} ({passed.sum()}/{len(rows)} strict)")
        axes[0].scatter(x[~passed[order]], values[order][~passed[order]], marker="x",
                        s=30, linewidths=1.2, color="black", zorder=3)
        axes[1].bar(index, passed.mean(), color=color, width=0.65)
        axes[1].text(index, passed.mean() + 0.02,
                     f"{passed.sum()}/{len(rows)}\nbest {values[passed].max():.3f}",
                     ha="center", va="bottom", fontsize=9)

    axes[0].axhline(0.80, color="0.45", ls="--", lw=1, label="seed1 gate 0.80")
    axes[0].axhline(0.90, color="#e76f51", ls=":", lw=1.4, label="target 0.90")
    axes[0].set(xlabel="candidate rank within round", ylabel="AF3 seed1 mean ipTM",
                ylim=(0.2, 0.92))
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)
    axes[1].set_xticks(range(len(rounds)), [label for label, _ in rounds])
    axes[1].set(ylabel="strict seed1 pass fraction", ylim=(0, 1.05))
    axes[1].grid(axis="y", alpha=0.2)
    figure.suptitle("AF3-labelled allele active-learning rounds")
    figure.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=220, bbox_inches="tight")
    figure.savefig(args.out.with_suffix(".svg"), bbox_inches="tight")
    print(args.out)


if __name__ == "__main__":
    main()
