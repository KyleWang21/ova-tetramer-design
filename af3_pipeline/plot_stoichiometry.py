#!/usr/bin/env python3
"""Plot AF3 chain-count competition from scored experiment tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--source",
        action="append",
        help="system_scores.tsv:SYSTEM_NAME (repeatable)",
    )
    ap.add_argument("--summary-table", type=Path,
                    help="TSV with n_chains, median_iptm, median_min_incident_iptm")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--title", default="AF3 stoichiometry competition")
    args = ap.parse_args()

    records = []
    if args.summary_table:
        with args.summary_table.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                passed = row.get("confidence_pass_models")
                models = row.get("models")
                records.append({"n": int(row["n_chains"]),
                                "iptm": float(row["median_iptm"]),
                                "weak": float(row["median_min_incident_iptm"]),
                                "pass_text": f"{passed}/{models} pass" if passed and models else ""})
    for spec in args.source or []:
        table_name, system = spec.rsplit(":", 1)
        with Path(table_name).open() as handle:
            row = next(row for row in csv.DictReader(handle, delimiter="\t")
                       if row["system"] == system)
        records.append({"n": int(row["n_chains"]), "iptm": float(row["median_iptm"]),
                        "weak": float(row["median_min_incident_iptm"]), "pass_text": ""})
    if not records:
        ap.error("provide --summary-table and/or --source")
    records.sort(key=lambda row: row["n"])

    x = [row["n"] for row in records]
    iptm = [row["iptm"] for row in records]
    weak = [row["weak"] for row in records]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(x, iptm, "o-", lw=2.2, ms=7, label="median ipTM", color="#2266cc")
    ax.plot(x, weak, "s--", lw=2.0, ms=6, label="weakest-chain interface", color="#d13c3c")
    ax.axhline(0.5, color="#777777", lw=1, ls=":", label="0.5 screening threshold")
    for row in records:
        label = f"{row['iptm']:.2f}"
        if row.get("pass_text"):
            label += f"\n{row['pass_text']}"
        ax.annotate(label, (row["n"], row["iptm"]),
                    xytext=(0, 8), textcoords="offset points", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xlabel("Forced AF3 chain count")
    ax.set_ylabel("Interface confidence")
    ax.set_ylim(0, max(0.85, max(iptm + weak) + 0.08))
    ax.set_title(args.title)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    print(args.output)


if __name__ == "__main__":
    main()
