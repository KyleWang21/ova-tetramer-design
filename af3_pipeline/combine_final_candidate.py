#!/usr/bin/env python3
"""Combine pilot and reseed AF3 tables for one final pure-OVA candidate."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import statistics

import matplotlib.pyplot as plt


def selected_rows(specs: list[str]) -> list[dict[str, str]]:
    rows = []
    for spec in specs:
        path_text, system = spec.rsplit(":", 1)
        rows.extend(
            row for row in csv.DictReader(pathlib.Path(path_text).open(), delimiter="\t")
            if row["system"] == system
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-source", action="append", required=True)
    ap.add_argument("--ss-source", action="append", required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--label", default="OVA-BODY-C4-SS3-AI03")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    models = selected_rows(args.model_source)
    ss = selected_rows(args.ss_source)
    noncovalent = bool(ss) and all(int(row["expected_ss_bonds"]) == 0 for row in ss)
    clean_field = "noncovalent_chemically_clean" if noncovalent else "chemically_clean_ring"
    ss_by_model = {row["model"]: row for row in ss}
    if not models or len(models) != len(ss) or set(ss_by_model) != {row["model"] for row in models}:
        raise SystemExit(f"model/disulfide rows do not match: {len(models)} vs {len(ss)}")
    combined = []
    for row in models:
        chem = ss_by_model[row["model"]]
        confidence = (
            float(row["iptm"]) >= .50 and float(row["min_incident_iptm"]) >= .45
            and int(row["full_connected"]) == 1 and float(row["has_clash"]) == 0
        )
        combined.append({
            "seed": int(re.match(r"seed-(\d+)", row["model"]).group(1)),
            "model": row["model"], "iptm": float(row["iptm"]),
            "weakest_chain": float(row["min_incident_iptm"]),
            "interface_pae_a": float(row["contact_pair_pae_median"]),
            "min_ova_plddt": float(row["min_ova_plddt"]),
            "full_connected": int(row["full_connected"]), "has_clash": int(float(row["has_clash"])),
            "confidence_pass": int(confidence),
            "observed_ss_bonds": int(chem["observed_ss_bonds"]),
            "native_c74_c121_bonds": int(chem["native_c74_c121_bonds"]),
            "unexpected_ss_bonds": int(chem["unexpected_ss_bonds"]),
            "chemically_clean": int(chem[clean_field]),
            "joint_pass": int(confidence and int(chem[clean_field])),
            "cif_path": row["cif_path"],
        })
    combined.sort(key=lambda row: (row["seed"], row["model"]))
    with (args.out / "combined_model_scores.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(combined[0]), delimiter="\t"); writer.writeheader(); writer.writerows(combined)
    summary = {
        "candidate": args.label, "seeds": sorted({row["seed"] for row in combined}),
        "n_models": len(combined),
        "median_iptm": statistics.median(row["iptm"] for row in combined),
        "median_weakest_chain": statistics.median(row["weakest_chain"] for row in combined),
        "median_interface_pae_a": statistics.median(row["interface_pae_a"] for row in combined),
        "median_min_ova_plddt": statistics.median(row["min_ova_plddt"] for row in combined),
        "confidence_passing_models": sum(row["confidence_pass"] for row in combined),
        "full_connected_models": sum(row["full_connected"] for row in combined),
        "unclashed_models": sum(not row["has_clash"] for row in combined),
        "chemical_mode": "noncovalent" if noncovalent else "designed_disulfide",
        "chemically_clean_models": sum(row["chemically_clean"] for row in combined),
        "joint_passing_models": sum(row["joint_pass"] for row in combined),
        "unexpected_ss_models": sum(row["unexpected_ss_bonds"] > 0 for row in combined),
    }
    (args.out / "combined_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    per_seed = []
    for seed in summary["seeds"]:
        rows = [row for row in combined if row["seed"] == seed]
        per_seed.append({
            "seed": seed, "models": len(rows),
            "median_iptm": statistics.median(row["iptm"] for row in rows),
            "median_weakest_chain": statistics.median(row["weakest_chain"] for row in rows),
            "confidence_pass": sum(row["confidence_pass"] for row in rows),
            "chemically_clean": sum(row["chemically_clean"] for row in rows),
            "joint_pass": sum(row["joint_pass"] for row in rows),
        })
    with (args.out / "per_seed_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(per_seed[0]), delimiter="\t"); writer.writeheader(); writer.writerows(per_seed)
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.5), constrained_layout=True)
    seeds = [row["seed"] for row in per_seed]
    axes[0].plot(seeds, [row["median_iptm"] for row in per_seed], "o-", lw=2.2, label="median ipTM")
    axes[0].plot(seeds, [row["median_weakest_chain"] for row in per_seed], "s--", lw=2,
                 label="median weakest-chain interface")
    axes[0].axhline(.5, color="#777777", ls=":", lw=1); axes[0].set_ylim(0, .9)
    axes[0].set_ylabel("AF3 confidence"); axes[0].legend(frameon=False); axes[0].grid(axis="y", alpha=.2)
    axes[1].bar([seed - .18 for seed in seeds], [row["confidence_pass"] for row in per_seed], .36,
                color="#2563eb", label="confidence pass / 5")
    clean_label = "zero interchain S-S / 5" if noncovalent else "clean designed S-S / 5"
    axes[1].bar([seed + .18 for seed in seeds], [row["chemically_clean"] for row in per_seed], .36,
                color="#e6a000", label=clean_label)
    axes[1].set_ylim(0, 5.4); axes[1].set_xticks(seeds); axes[1].set_xlabel("Independent AF3 seed")
    axes[1].set_ylabel("models"); axes[1].legend(frameon=False, ncol=2); axes[1].grid(axis="y", alpha=.2)
    fig.suptitle(f"{args.label}: multi-seed AF3 robustness")
    fig.savefig(args.out / "multiseed_robustness.png", dpi=220)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
