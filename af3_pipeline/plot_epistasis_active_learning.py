#!/usr/bin/env python3
"""Plot held-out epistasis performance and proposal re-ranking diagnostics."""

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
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    args = parser.parse_args()
    labels = read_tsv(args.experiment / "epistasis_label_cv_predictions.tsv")
    library = read_tsv(args.experiment / "epistasis_all_library.tsv")
    cv = read_tsv(args.experiment / "epistasis_cv.tsv")[0]

    observed = np.asarray([float(row["observed_af3_iptm"]) for row in labels])
    predicted = np.asarray([float(row["heldout_epistasis_prediction"]) for row in labels])
    kernel = np.asarray([float(row["predicted_af3"]) for row in library])
    epistasis = np.asarray([float(row["epistasis_robust_prediction"]) for row in library])
    source_rank = np.asarray([float(row["rank"]) for row in library])
    blended_rank = np.asarray([float(row["blended_selection_average_rank"]) for row in library])
    selected = np.asarray([int(row["epistasis_selected_rank"]) > 0 for row in library])

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
    axes[0].scatter(observed, predicted, s=22, alpha=0.65, color="#3569a8")
    lower = min(observed.min(), predicted.min()); upper = max(observed.max(), predicted.max())
    axes[0].plot([lower, upper], [lower, upper], "--", color="0.45", lw=1)
    axes[0].set(
        xlabel="Observed AF3 ipTM", ylabel="Held-out ridge prediction",
        title=(
            "Cluster-blocked five-fold CV "
            f"(Spearman {float(cv['clustered_cv_spearman']):.3f})"
        ),
    )

    axes[1].scatter(kernel[~selected], epistasis[~selected], s=20, alpha=0.35, color="0.55", label="proposal")
    axes[1].scatter(kernel[selected], epistasis[selected], s=34, alpha=0.9, color="#d1495b", label="selected Top32")
    axes[1].set(xlabel="Hamming-kernel AF3 prediction", ylabel="Robust epistasis score",
                title="Independent proposal signals")
    axes[1].legend(frameon=False)

    axes[2].scatter(source_rank[~selected], blended_rank[~selected], s=20, alpha=0.30, color="0.55")
    axes[2].scatter(source_rank[selected], blended_rank[selected], s=34, alpha=0.9, color="#2a9d8f")
    axes[2].set(xlabel="Original cross-model rank", ylabel="Blended average-rank score",
                title="Epistasis rescues non-obvious combinations")
    axes[2].invert_yaxis()
    fig.tight_layout()
    args.experiment.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.experiment / "epistasis_active_learning.png", dpi=220)
    fig.savefig(args.experiment / "epistasis_active_learning.pdf")
    print(args.experiment / "epistasis_active_learning.png")


if __name__ == "__main__":
    main()
