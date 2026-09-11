#!/usr/bin/env python3
"""Plot chain-map repair and live eight-state tied-AF2 progress for E283."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np


PROGRESS = re.compile(r"^(logits|soft|hard) (\d+)/\d+ .* ipTM=([0-9./-]+) nmut=(\d+)")


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    parser.add_argument("--logs", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    mapping = read_tsv(args.experiment / "chain_mapping_audit.tsv")
    progress: dict[tuple[str, int], list[tuple[float, float, int]]] = {}
    for path in sorted(args.logs.glob("*.out")):
        for line in path.read_text().splitlines():
            match = PROGRESS.match(line)
            if not match:
                continue
            stage, iteration_text, values_text, nmut_text = match.groups()
            values = [float(value) for value in values_text.split("/")]
            if len(values) != 8:
                raise RuntimeError(f"expected eight states in {path}: {line}")
            key = (stage, int(iteration_text))
            progress.setdefault(key, []).append((min(values), float(np.mean(values)), int(nmut_text)))
    if not progress:
        raise RuntimeError(f"no progress records found in {args.logs}")

    stage_order = {"logits": 0, "soft": 1, "hard": 2}
    ordered = sorted(progress, key=lambda key: (stage_order[key[0]], key[1]))
    summary_rows: list[dict[str, object]] = []
    for global_iteration, key in enumerate(ordered, 1):
        observations = progress[key]
        if len(observations) != 8:
            continue
        minima = np.asarray([value[0] for value in observations])
        means = np.asarray([value[1] for value in observations])
        summary_rows.append({
            "global_iteration": global_iteration,
            "stage": key[0],
            "stage_iteration": key[1],
            "n_tasks": len(observations),
            "best_task_min_iptm": float(minima.max()),
            "median_task_min_iptm": float(np.median(minima)),
            "best_task_mean_iptm": float(means.max()),
            "median_task_mean_iptm": float(np.median(means)),
        })

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "live_progress_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), constrained_layout=True)
    x = np.arange(1, len(mapping) + 1)
    legacy = np.asarray([float(row["legacy_jaccard"]) for row in mapping])
    selected = np.asarray([float(row["selected_jaccard"]) for row in mapping])
    axes[0].bar(x - 0.18, legacy, 0.36, label="legacy chain labels", color="#b7b7b7")
    axes[0].bar(x + 0.18, selected, 0.36, label="contact-map remapped", color="#2878b5")
    axes[0].axhline(0.70, color="#d62728", linestyle="--", linewidth=1, label="mapping gate 0.70")
    axes[0].set(xlabel="interface state", ylabel="5.3 Å contact-map Jaccard", ylim=(0, 1.08), xticks=x)
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")
    axes[0].set_title("Physical interface consistency")

    iteration = np.asarray([int(row["global_iteration"]) for row in summary_rows])
    axes[1].plot(
        iteration, [float(row["best_task_min_iptm"]) for row in summary_rows],
        marker="o", label="best task: weakest state", color="#d62728",
    )
    axes[1].plot(
        iteration, [float(row["median_task_min_iptm"]) for row in summary_rows],
        marker="o", label="median task: weakest state", color="#ff9896",
    )
    axes[1].plot(
        iteration, [float(row["best_task_mean_iptm"]) for row in summary_rows],
        marker="s", label="best task: 8-state mean", color="#2878b5",
    )
    axes[1].set(
        xlabel="optimization iteration", ylabel="tied-AF2 ipTM proxy",
        ylim=(0, 0.8), xticks=iteration,
    )
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].set_title("E283 live robust-state trajectory")
    figure.suptitle("E283 chain-remapped multibackbone design (not template-free AF3)", fontsize=11)
    figure.savefig(args.out / "e283_chainmapping_progress.png", dpi=220)
    figure.savefig(args.out / "e283_chainmapping_progress.pdf")
    print(f"wrote {args.out} with {len(summary_rows)} complete iterations")


if __name__ == "__main__":
    main()
