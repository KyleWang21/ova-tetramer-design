#!/usr/bin/env python3
"""Summarize partial AF3 screens and expose mutation/family effects."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def fasta_record(path: pathlib.Path, prefix: str) -> str:
    active = False
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if active:
                break
            active = raw[1:].startswith(prefix)
        elif active:
            chunks.append(raw.strip())
    return "".join(chunks)


def mutations(reference: str, sequence: str) -> list[str]:
    return [
        f"{old}{i}{new}"
        for i, (old, new) in enumerate(zip(reference, sequence), 1)
        if old != new
    ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--reference-name", default="D0-P3-13R")
    ap.add_argument("--hit-threshold", type=float, default=0.78)
    args = ap.parse_args()

    reference = fasta_record(args.reference_fasta, args.reference_name)
    manifest = {
        row["name"]: row
        for row in csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t")
    }
    scores: dict[str, list[float]] = {}
    for path in args.experiment.glob("out_s*/**/*_summary_confidences.json"):
        name = path.name.split("_seed-")[0]
        scores.setdefault(name, []).append(float(json.loads(path.read_text())["iptm"]))

    rows: list[dict[str, object]] = []
    for name, values in scores.items():
        if len(values) != 5:
            continue
        sequence = manifest[name]["sequence"]
        muts = mutations(reference, sequence)
        rows.append({
            "system": name,
            "mean_iptm": float(np.mean(values)),
            "min_iptm": min(values),
            "max_iptm": max(values),
            "n_mutations": len(muts),
            "mutations": ",".join(muts),
            "description": manifest[name]["description"],
            "sequence": sequence,
        })
    rows.sort(key=lambda row: (-float(row["mean_iptm"]), str(row["system"])))
    if not rows:
        raise SystemExit("no complete five-model systems")

    with (args.experiment / "partial_af3_means.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    all_mutations = sorted({m for row in rows for m in str(row["mutations"]).split(",") if m})
    effects: list[dict[str, object]] = []
    global_mean = float(np.mean([float(row["mean_iptm"]) for row in rows]))
    for mutation in all_mutations:
        present = [float(row["mean_iptm"]) for row in rows if mutation in str(row["mutations"]).split(",")]
        absent = [float(row["mean_iptm"]) for row in rows if mutation not in str(row["mutations"]).split(",")]
        effects.append({
            "mutation": mutation,
            "n_present": len(present),
            "mean_present": float(np.mean(present)),
            "mean_absent": float(np.mean(absent)) if absent else global_mean,
            "mean_difference": float(np.mean(present)) - (float(np.mean(absent)) if absent else global_mean),
            "hit_fraction": float(np.mean(np.asarray(present) >= args.hit_threshold)),
        })
    effects.sort(key=lambda row: -float(row["mean_difference"]))
    with (args.experiment / "partial_mutation_effects.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(effects[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(effects)

    hit_rows = [row for row in rows if float(row["mean_iptm"]) >= args.hit_threshold]
    hit_mutations = sorted(
        {m for row in hit_rows for m in str(row["mutations"]).split(",") if m},
        key=lambda m: int(m[1:-1]),
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), constrained_layout=True)
    order = list(reversed(rows))
    axes[0].barh(
        [str(row["system"]).replace("ova_body_c4_", "") for row in order],
        [float(row["mean_iptm"]) for row in order],
        color=["#0f766e" if float(row["mean_iptm"]) > 0.80 else "#94a3b8" for row in order],
    )
    axes[0].axvline(0.80, color="#dc2626", linestyle="--", linewidth=1.5)
    axes[0].set(xlabel="AF3 seed-1 mean ipTM (5 samples)", ylabel="candidate")
    axes[0].set_xlim(0, 0.9)

    if hit_rows and hit_mutations:
        matrix = np.asarray([
            [mutation in str(row["mutations"]).split(",") for mutation in hit_mutations]
            for row in hit_rows
        ], dtype=float)
        axes[1].imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
        axes[1].set_xticks(range(len(hit_mutations)), hit_mutations, rotation=90, fontsize=7)
        axes[1].set_yticks(
            range(len(hit_rows)),
            [f"{str(row['system']).replace('ova_body_c4_', '')} ({float(row['mean_iptm']):.3f})" for row in hit_rows],
        )
        axes[1].set_title(f"mutation patterns for mean ipTM ≥ {args.hit_threshold:.2f}")
    else:
        axes[1].text(0.5, 0.5, "no hits yet", ha="center", va="center")
        axes[1].axis("off")
    fig.suptitle(f"AF3 trajectory screen: {len(rows)} complete candidates")
    fig.savefig(args.experiment / "partial_af3_active_learning.png", dpi=180)
    fig.savefig(args.experiment / "partial_af3_active_learning.svg")

    print(f"complete={len(rows)} hits={sum(float(r['mean_iptm']) > 0.80 for r in rows)}")
    for row in rows[:10]:
        print(row["system"], f"mean={float(row['mean_iptm']):.3f}", f"nmut={row['n_mutations']}")


if __name__ == "__main__":
    main()
