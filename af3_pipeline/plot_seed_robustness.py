#!/usr/bin/env python3
"""Combine AF3 model tables and plot per-seed oligomer confidence."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", action="append", required=True,
                    help="model_scores.tsv:SYSTEM_NAME (repeatable)")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--title", required=True)
    args = ap.parse_args()

    rows = []
    for spec in args.source:
        table_name, system = spec.rsplit(":", 1)
        with Path(table_name).open() as handle:
            rows.extend(row for row in csv.DictReader(handle, delimiter="\t")
                        if row["system"] == system)
    grouped = {}
    for row in rows:
        match = re.match(r"seed-(\d+)_", row["model"])
        if match:
            grouped.setdefault(int(match.group(1)), []).append(row)

    summary = []
    for seed, values in sorted(grouped.items()):
        iptm = [float(row["iptm"]) for row in values]
        weak = [float(row["min_incident_iptm"]) for row in values]
        passes = sum(a >= 0.5 and b >= 0.45 for a, b in zip(iptm, weak))
        summary.append({"seed": seed, "models": len(values),
                        "median_iptm": statistics.median(iptm),
                        "median_weakest_chain": statistics.median(weak),
                        "passing_models": passes})
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    seeds = [row["seed"] for row in summary]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(seeds, [row["median_iptm"] for row in summary], "o-", lw=2.2,
            color="#2767c9", label="median ipTM")
    ax.plot(seeds, [row["median_weakest_chain"] for row in summary], "s--", lw=2,
            color="#d23d3d", label="weakest-chain interface")
    ax.axhline(0.5, color="#777777", ls=":", lw=1, label="ipTM 0.5 threshold")
    for row in summary:
        ax.text(row["seed"], 0.04, f"{row['passing_models']}/{row['models']}",
                ha="center", va="bottom", fontsize=9)
    ax.text(0.01, 0.02, "bottom labels: passing models", transform=ax.transAxes,
            fontsize=8, color="#555555")
    ax.set_xticks(seeds)
    ax.set_xlabel("Independent AF3 seed")
    ax.set_ylabel("Interface confidence")
    ax.set_ylim(0, 0.8)
    ax.set_title(args.title)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    print(args.output)


if __name__ == "__main__":
    main()
