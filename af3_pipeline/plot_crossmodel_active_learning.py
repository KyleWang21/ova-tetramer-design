#!/usr/bin/env python3
"""Plot active-learning proposal trade-offs and selected structural seeds."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=pathlib.Path, required=True)
    parser.add_argument("--selected", type=pathlib.Path, required=True)
    parser.add_argument("--model-usage", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    rows = read_tsv(args.candidates)
    selected = {row["surrogate_name"] for row in read_tsv(args.selected)}
    af3 = np.asarray([float(row["predicted_af3"]) for row in rows])
    protenix = np.asarray([float(row["predicted_protenix"]) for row in rows])
    opendde = np.asarray([float(row["predicted_opendde"]) for row in rows])
    nmut = np.asarray([int(row["n_mutations"]) for row in rows])
    surface = np.asarray([float(row["surface_mutation_fraction"]) for row in rows])
    selected_mask = np.asarray([row["name"] in selected for row in rows])

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.6), constrained_layout=True)
    scatter = axes[0].scatter(
        af3, protenix, c=opendde, cmap="viridis", s=24, alpha=0.75,
        edgecolors="none",
    )
    axes[0].scatter(
        af3[selected_mask], protenix[selected_mask], marker="*", s=150,
        facecolors="#e63946", edgecolors="black", linewidths=0.7, label="E235 seeds",
    )
    axes[0].set(xlabel="Predicted AF3 ipTM", ylabel="Predicted Protenix ipTM",
                title="Cross-model proposal trade-off")
    axes[0].legend(frameon=False, loc="lower right")
    colorbar = fig.colorbar(scatter, ax=axes[0], fraction=0.047, pad=0.03)
    colorbar.set_label("Predicted OpenDDE ipTM")

    axes[1].scatter(nmut, af3, c=surface, cmap="plasma", s=24, alpha=0.75, edgecolors="none")
    axes[1].scatter(
        nmut[selected_mask], af3[selected_mask], marker="*", s=150,
        facecolors="#e63946", edgecolors="black", linewidths=0.7,
    )
    axes[1].axhline(0.90, color="black", linestyle="--", linewidth=0.9, label="target 0.90")
    axes[1].set(xlabel="Mutations vs P3-13R", ylabel="Predicted AF3 ipTM",
                title="Mutation budget and AF3 target", xticks=sorted(set(nmut)))
    axes[1].legend(frameon=False, loc="lower left")

    usage = read_tsv(args.model_usage)
    names = [row["outcome"] for row in usage]
    cv = [float(row["cv_spearman"]) for row in usage]
    colors = ["#2a9d8f" if int(row["used"]) else "#b0b0b0" for row in usage]
    axes[2].barh(names, cv, color=colors)
    axes[2].axvline(0.20, color="black", linestyle="--", linewidth=0.9)
    axes[2].invert_yaxis()
    axes[2].set(xlabel="Five-fold Spearman", title="Surrogate signal quality", xlim=(0, 0.9))
    for index, value in enumerate(cv):
        axes[2].text(value + 0.015, index, f"{value:.2f}", va="center")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    fig.savefig(args.out.with_suffix(".pdf"))
    print(f"wrote {args.out} and {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
