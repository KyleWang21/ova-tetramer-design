#!/usr/bin/env python3
"""Plot mean AF3 ipTM against the strict zero-clash sample fraction."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--title", default="AF3 seed1: assembly confidence vs strict geometry")
    args = parser.parse_args()
    with args.summary.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("empty seed1 summary")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(8.2, 5.8))
    strict_names = {
        row["system"] for row in rows if int(row["seed1_strict_pass"]) == 1
    }
    top_failed = sorted(
        (row for row in rows if int(row["seed1_strict_pass"]) == 0),
        key=lambda row: -float(row["mean_iptm"]),
    )[:3]
    label_names = strict_names | {row["system"] for row in top_failed}
    label_index = {
        name: index for index, name in enumerate(sorted(
            label_names,
            key=lambda value: float(next(row["mean_iptm"] for row in rows if row["system"] == value)),
        ))
    }
    for row in rows:
        x = float(row["mean_iptm"])
        y = int(row["zero_atomic_clash_models"]) / 5.0
        passed = int(row["seed1_strict_pass"]) == 1
        color = "#168a55" if passed else ("#d47a1f" if x >= 0.80 else "#8a8f98")
        axis.scatter(x, y, s=58, color=color, edgecolor="white", linewidth=0.6, zorder=3)
        if row["system"] in label_names:
            index = label_index[row["system"]]
            axis.annotate(
                row["system"].removeprefix("ova_body_c4_"), (x, y),
                xytext=(4, 5 + 7 * index), textcoords="offset points", fontsize=8,
            )
    axis.axvline(0.80, color="#777777", linestyle="--", linewidth=1, label="seed1 ipTM gate 0.80")
    axis.axvline(0.90, color="#6a55a3", linestyle=":", linewidth=1.4, label="target mean ipTM 0.90")
    axis.axhline(1.0, color="#444444", linestyle="--", linewidth=1, label="5/5 zero clash")
    axis.set_xlim(min(0.20, min(float(row["mean_iptm"]) for row in rows) - 0.02), 0.92)
    axis.set_ylim(-0.04, 1.08)
    axis.set_xlabel("AF3 seed1 arithmetic mean ipTM (5 samples)")
    axis.set_ylabel("Zero-clash sample fraction")
    axis.set_title(args.title)
    axis.grid(alpha=0.18, linewidth=0.7)
    axis.legend(loc="lower right", fontsize=8, frameon=True)
    figure.tight_layout()
    figure.savefig(args.out, dpi=220)
    figure.savefig(args.out.with_suffix(".pdf"))
    plt.close(figure)
    print(f"wrote {args.out} rows={len(rows)}")


if __name__ == "__main__":
    main()
