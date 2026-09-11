#!/usr/bin/env python3
"""Plot recovered tied-AF2 trajectory candidates and the E257 shortlist."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


COLORS = {"E217": "#0072B2", "E219": "#D55E00", "E221": "#009E73"}


def read(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def source(row: dict[str, str]) -> str:
    text = row.get("recovered_library", "")
    for label in COLORS:
        if label.lower() in text.lower():
            return label
    return "other"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=pathlib.Path, required=True)
    parser.add_argument("--shortlist", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    rows = read(args.library)
    selected = {row["sequence"] for row in read(args.shortlist)}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)
    for label, color in COLORS.items():
        group = [row for row in rows if source(row) == label]
        axes[0].scatter(
            [float(row["af2_iptm"]) for row in group],
            [float(row["predicted_af3_mean_iptm"]) for row in group],
            s=18, alpha=0.45, color=color, label=f"{label} (n={len(group)})",
        )
        axes[1].scatter(
            [int(float(row["n_mutations"])) for row in group],
            [float(row["predicted_af3_mean_iptm"]) for row in group],
            s=18, alpha=0.45, color=color,
        )
    top = [row for row in rows if row["sequence"] in selected]
    axes[0].scatter(
        [float(row["af2_iptm"]) for row in top],
        [float(row["predicted_af3_mean_iptm"]) for row in top],
        s=58, facecolors="none", edgecolors="black", linewidths=1.0,
        label=f"E257 Top32",
    )
    axes[1].scatter(
        [int(float(row["n_mutations"])) for row in top],
        [float(row["predicted_af3_mean_iptm"]) for row in top],
        s=58, facecolors="none", edgecolors="black", linewidths=1.0,
    )
    axes[0].set(
        xlabel="Weakest tied-AF2 interface ipTM",
        ylabel="Cross-validated AF3 surrogate mean ipTM",
        title="Recovered near-discrete trajectory states",
    )
    axes[1].set(
        xlabel="Mutations versus P3-13R",
        ylabel="Cross-validated AF3 surrogate mean ipTM",
        title="Mutation budget and E257 selection",
        xticks=range(12, 25, 2),
    )
    axes[0].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.18)
        ax.set_ylim(0.25, 0.84)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix(".png"), dpi=220)
    fig.savefig(args.out.with_suffix(".pdf"))
    print(f"plotted library={len(rows)} selected={len(top)} to {args.out}")


if __name__ == "__main__":
    main()
