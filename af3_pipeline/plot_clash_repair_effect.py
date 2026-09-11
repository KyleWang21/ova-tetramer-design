#!/usr/bin/env python3
"""Compare AF3 seed1 results before and after clash-directed sequence repairs."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    ranking = {int(row["ai_rank"]): row for row in read_tsv(args.ranking)}
    summary = {row["system"]: row for row in read_tsv(args.summary)}
    joined: list[dict[str, object]] = []
    for manifest_row in read_tsv(args.manifest):
        child = summary.get(manifest_row["name"])
        if child is None:
            continue
        ai_rank = int(
            manifest_row["description"].split("source_ai_rank=", 1)[1].split()[0]
        )
        source = ranking[ai_rank]
        joined.append(
            {
                "candidate": manifest_row["name"],
                "repair": source["repair_edits"],
                "parent_iptm": float(source["parent_mean_iptm"]),
                "child_iptm": float(child["mean_iptm"]),
                "parent_clean": int(source["parent_zero_clash_models"]) / 5.0,
                "child_clean": int(child["zero_atomic_clash_models"]) / 5.0,
                "passed": int(child["seed1_strict_pass"]) == 1,
            }
        )
    if not joined:
        raise ValueError("no completed repair candidates could be joined")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table_path = args.out.with_suffix(".tsv")
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(joined[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(joined)

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.8))
    panels = (
        ("parent_iptm", "child_iptm", "Mean AF3 ipTM", (0.72, 0.88)),
        ("parent_clean", "child_clean", "Zero-clash sample fraction", (-0.04, 1.04)),
    )
    for axis, (x_key, y_key, label, limits) in zip(axes, panels):
        for row in joined:
            color = "#168a55" if row["passed"] else "#d47a1f"
            axis.scatter(row[x_key], row[y_key], s=48, color=color,
                         edgecolor="white", linewidth=0.5, zorder=3)
        axis.plot(limits, limits, color="#777777", linestyle="--", linewidth=1)
        axis.set_xlim(*limits)
        axis.set_ylim(*limits)
        axis.set_xlabel(f"Parent {label.lower()}")
        axis.set_ylabel(f"Repaired {label.lower()}")
        axis.grid(alpha=0.18, linewidth=0.7)
    axes[0].axhline(0.80, color="#555555", linestyle=":", linewidth=1)
    axes[1].axhline(1.0, color="#555555", linestyle=":", linewidth=1)
    figure.suptitle(
        f"AF3 clash-directed repair effect ({len(joined)} completed candidates)"
    )
    figure.tight_layout()
    figure.savefig(args.out, dpi=220)
    figure.savefig(args.out.with_suffix(".pdf"))
    plt.close(figure)
    print(f"wrote {args.out} and {table_path} rows={len(joined)}")


if __name__ == "__main__":
    main()
