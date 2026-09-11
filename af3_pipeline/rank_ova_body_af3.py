#!/usr/bin/env python3
"""Rank and plot no-module OVA C4 candidates from AF3 model/system scores."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiment", type=pathlib.Path, required=True)
    ap.add_argument("--iptm-min", type=float, default=0.50)
    ap.add_argument("--weak-min", type=float, default=0.45)
    ap.add_argument("--disulfide-table", type=pathlib.Path,
                    help="Optional disulfide_model_scores.tsv; require chemically_clean_ring")
    args = ap.parse_args()
    manifest = {row["name"]: row for row in csv.DictReader((args.experiment / "manifest.tsv").open(), delimiter="\t")}
    systems = list(csv.DictReader((args.experiment / "system_scores.tsv").open(), delimiter="\t"))
    models = list(csv.DictReader((args.experiment / "model_scores.tsv").open(), delimiter="\t"))
    disulfide = {}
    chemical_mode = "none"
    if args.disulfide_table:
        disulfide_rows = list(csv.DictReader(args.disulfide_table.open(), delimiter="\t"))
        disulfide = {(row["system"], row["model"]): row for row in disulfide_rows}
        chemical_mode = (
            "noncovalent" if disulfide_rows
            and all(int(row["expected_ss_bonds"]) == 0 for row in disulfide_rows)
            else "designed_disulfide"
        )
    clean_field = (
        "noncovalent_chemically_clean" if chemical_mode == "noncovalent"
        else "chemically_clean_ring"
    )
    output: list[dict[str, object]] = []
    for system in systems:
        name = system["system"]
        rows = [row for row in models if row["system"] == name]
        confidence_passed = [
            row for row in rows
            if float(row["iptm"]) >= args.iptm_min
            and float(row["min_incident_iptm"]) >= args.weak_min
            and int(row["full_connected"]) == 1
            and float(row["has_clash"]) == 0
        ]
        clean = [row for row in rows
                 if (name, row["model"]) in disulfide
                 and int(disulfide[(name, row["model"])][clean_field]) == 1]
        passed = [row for row in confidence_passed
                  if not disulfide or row in clean]
        description = manifest[name]["description"]
        mutation_text = description.split("|mut=", 1)[1].split("|", 1)[0] if "|mut=" in description else ""
        if mutation_text:
            n_mutations = len(mutation_text.split(","))
        elif "mutations=" in description:
            n_mutations = int(description.split("mutations=", 1)[1].split()[0])
        elif "nmut=" in description:
            n_mutations = int(description.split("nmut=", 1)[1].split()[0])
        else:
            n_mutations = 0
        output.append({
            "rank": 0,
            "system": name,
            "mutations": mutation_text,
            "n_mutations": n_mutations,
            "n_models": len(rows),
            "passing_models": len(passed),
            "pass_fraction": len(passed) / len(rows),
            "confidence_passing_models": len(confidence_passed),
            "chemical_mode": chemical_mode,
            "chemically_clean_models": len(clean) if disulfide else "",
            "chemically_clean_fraction": len(clean) / len(rows) if disulfide else "",
            "mean_iptm": float(np.mean([float(row["iptm"]) for row in rows])),
            "mean_weakest_chain": float(np.mean([float(row["min_incident_iptm"]) for row in rows])),
            "min_iptm": min(float(row["iptm"]) for row in rows),
            "median_iptm": float(system["median_iptm"]),
            "median_weakest_chain": float(system["median_min_incident_iptm"]),
            "median_interface_pae_a": float(system["median_contact_pair_pae"]),
            "median_min_ova_plddt": float(system["median_min_ova_plddt"]),
            "full_connected_fraction": float(system["full_connected_fraction"]),
            "unclashed_fraction": float(system["unclashed_fraction"]),
            "median_interface_jaccard": float(system["median_interface_jaccard"]),
            "best_iptm": float(system["best_iptm"]),
            "best_weakest_chain": float(system["best_min_incident_iptm"]),
            "best_cif": system["best_cif"],
            "sequence": manifest[name]["sequence"],
        })
    output.sort(key=lambda row: (
        int(row["passing_models"]),
        float(row["mean_iptm"]),
        float(row["median_weakest_chain"]),
        float(row["median_iptm"]),
        -float(row["median_interface_pae_a"]),
    ), reverse=True)
    for rank, row in enumerate(output, 1):
        row["rank"] = rank
    with (args.experiment / "body_candidate_ranking.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(output[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(output)

    labels = [str(row["system"]).replace("ova_body_c4_", "BODY-") for row in output]
    iptm = np.array([float(row["median_iptm"]) for row in output])
    weak = np.array([float(row["median_weakest_chain"]) for row in output])
    pass_fraction = np.array([float(row["pass_fraction"]) for row in output])
    x = np.arange(len(output))
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    width = 0.36
    axes[0].bar(x - width / 2, iptm, width, label="median ipTM", color="#1769aa")
    axes[0].bar(x + width / 2, weak, width, label="median weakest-chain ipTM", color="#e07a1f")
    axes[0].axhline(args.iptm_min, color="#1769aa", ls="--", lw=1)
    axes[0].axhline(args.weak_min, color="#e07a1f", ls="--", lw=1)
    axes[0].set_ylim(0, max(0.7, float(max(iptm.max(), weak.max())) + 0.08))
    axes[0].set_ylabel("AF3 confidence")
    axes[0].set_xticks(x, labels)
    axes[0].legend(frameon=False, ncol=2)
    colors = ["#2ca02c" if value >= 0.8 else "#999999" for value in pass_fraction]
    axes[1].bar(x, pass_fraction, color=colors)
    axes[1].axhline(0.8, color="#2ca02c", ls="--", lw=1)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("models passing both gates")
    axes[1].set_xticks(x, labels)
    fig.suptitle("Pure OVA-body C4 AF3 screen (no linker or oligomerization module)")
    figures = args.experiment / "figures"; figures.mkdir(exist_ok=True)
    fig.savefig(figures / "ova_body_c4_screen.png", dpi=200)
    fig.savefig(figures / "ova_body_c4_screen.svg")
    plt.close(fig)
    for row in output:
        print(
            row["rank"], row["system"],
            f"pass={row['passing_models']}/{row['n_models']}",
            f"ipTM={float(row['median_iptm']):.3f}",
            f"weak={float(row['median_weakest_chain']):.3f}",
        )


if __name__ == "__main__":
    main()
