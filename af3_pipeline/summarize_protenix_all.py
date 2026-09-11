#!/usr/bin/env python3
"""Join the uniform Protenix-v2 seed31 scores to the final AF3 candidate table."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re

import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def score_name(candidate: str) -> str:
    tag = re.sub(r"[^a-z0-9]+", "_", candidate.lower()).strip("_")
    return f"ova_iptm80_{tag}_protenix_score.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final-candidates", type=pathlib.Path, required=True)
    parser.add_argument("--score-dir", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for final in read_tsv(args.final_candidates):
        score_path = args.score_dir / score_name(final["candidate"])
        if not score_path.exists():
            raise ValueError(f"missing score for {final['candidate']}: {score_path}")
        score = json.loads(score_path.read_text())
        rows.append({
            "rank": int(final["rank"]),
            "candidate": final["candidate"],
            "source_system": final["source_system"],
            "n_mutations": int(final["n_mutations"]),
            "af3_mean_iptm_25models": float(final["mean_iptm"]),
            "protenix_seed": 31,
            "protenix_iptm": score["iptm"],
            "protenix_ptm": score["ptm"],
            "protenix_plddt": score["plddt"],
            "protenix_min_incident_iptm": score["min_incident_iptm"],
            "largest_component": score["largest_component"],
            "full_connected": int(score["full_connected"]),
            "has_clash": int(bool(score["has_clash"])),
            "native_c74_c121_bonds": score["native_c74_c121_bonds"],
            "interchain_ss_bonds": score["interchain_ss_bonds"],
            "contact_edges": ",".join(score["contact_edges"]),
            "score_json": str(score_path),
            "summary_json": score["summary"],
            "cif": score["cif"],
        })

    fields = list(rows[0])
    table = args.out / "all20_protenix_v2_seed31.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    ranked = sorted(rows, key=lambda row: float(row["protenix_iptm"]), reverse=True)
    with (args.out / "all20_protenix_v2_seed31_ranked.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, ["protenix_rank", *fields], delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(ranked, 1):
            writer.writerow({"protenix_rank": rank, **row})

    af3 = [float(row["af3_mean_iptm_25models"]) for row in rows]
    protenix = [float(row["protenix_iptm"]) for row in rows]
    rho, rho_p = spearmanr(af3, protenix)
    pearson, pearson_p = pearsonr(af3, protenix)
    summary = {
        "n_candidates": len(rows),
        "seed": 31,
        "protenix_iptm_min": min(protenix),
        "protenix_iptm_max": max(protenix),
        "protenix_iptm_mean": sum(protenix) / len(protenix),
        "n_iptm_ge_0_75": sum(value >= 0.75 for value in protenix),
        "n_iptm_ge_0_70": sum(value >= 0.70 for value in protenix),
        "n_iptm_ge_0_60": sum(value >= 0.60 for value in protenix),
        "full_connected": sum(int(row["full_connected"]) for row in rows),
        "unclashed": sum(not int(row["has_clash"]) for row in rows),
        "native_ss_complete": sum(int(row["native_c74_c121_bonds"]) == 4 for row in rows),
        "zero_interchain_ss": sum(int(row["interchain_ss_bonds"]) == 0 for row in rows),
        "spearman_af3_vs_protenix": rho,
        "spearman_pvalue": rho_p,
        "pearson_af3_vs_protenix": pearson,
        "pearson_pvalue": pearson_p,
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), gridspec_kw={"width_ratios": [1.45, 1]})
    axes[0].bar(
        [row["candidate"].replace("OVA-C4-IPTM80-", "") for row in ranked],
        [float(row["protenix_iptm"]) for row in ranked], color="#0072B2",
    )
    axes[0].set(ylabel="Protenix-v2 ipTM (seed31)", xlabel="Candidate", ylim=(0.25, 0.86),
                title="Uniform Protenix-v2 validation of all 20 candidates")
    axes[0].tick_params(axis="x", rotation=70)
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].scatter(af3, protenix, color="#D55E00", s=48)
    for row, x, y in zip(rows, af3, protenix):
        axes[1].annotate(row["candidate"].split("-")[-1], (x, y), xytext=(3, 2),
                         textcoords="offset points", fontsize=7)
    axes[1].set(xlabel="AF3 mean ipTM (25 models)", ylabel="Protenix-v2 ipTM (seed31)",
                title=f"AF3 vs Protenix: Spearman ρ={rho:.2f}")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out / "all20_protenix_v2_seed31.png", dpi=220)
    fig.savefig(args.out / "all20_protenix_v2_seed31.svg")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
