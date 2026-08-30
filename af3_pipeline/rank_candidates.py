#!/usr/bin/env python3
"""Rank completed AF3 tetramer candidates with a transparent multi-metric gate."""

from __future__ import annotations

import argparse
import csv
import pathlib


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source = args.experiment / "system_scores.tsv"
    rows = list(csv.DictReader(source.open(), delimiter="\t"))
    ranked = []
    for row in rows:
        if int(row["n_chains"]) != 4:
            continue
        pae = f(row, "median_contact_pair_pae")
        score = (
            0.15 * f(row, "median_iptm")
            + 0.10 * f(row, "median_min_incident_iptm")
            + 0.15 * f(row, "full_connected_fraction")
            + 0.15 * f(row, "module_connected_fraction")
            + 0.10 * min(1.0, f(row, "median_min_ova_plddt") / 100.0)
            + 0.10 * min(1.0, f(row, "median_min_module_plddt") / 100.0)
            + 0.15 * f(row, "median_interface_jaccard")
            + 0.10 * max(0.0, 1.0 - pae / 30.0)
        )
        passed = (
            int(row["complete"]) == 1
            and f(row, "median_iptm") >= 0.20
            and f(row, "median_min_incident_iptm") >= 0.15
            and pae <= 10.0
            and f(row, "full_connected_fraction") >= 0.80
            and f(row, "module_connected_fraction") >= 0.80
            and f(row, "unclashed_fraction") >= 0.90
            and f(row, "median_min_ova_plddt") >= 75.0
            and f(row, "median_min_module_plddt") >= 70.0
            and f(row, "median_interface_jaccard") >= 0.40
        )
        ranked.append({"score": round(score, 4), "pass_gate": int(passed), **row})
    ranked.sort(key=lambda row: (int(row["pass_gate"]), float(row["score"])), reverse=True)
    for index, row in enumerate(ranked, 1):
        row["rank"] = index
    if not ranked:
        raise SystemExit("no four-chain systems scored")
    fields = ["rank", "score", "pass_gate", *[k for k in ranked[0] if k not in {"rank", "score", "pass_gate"}]]
    output = args.experiment / "ranked_candidates.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ranked)
    print("rank\tpass\tscore\tsystem\tipTM\tweak\tPAE\tOVA_pLDDT\tmodule_pLDDT\tmodule_conn\tconsistent")
    for row in ranked:
        print(
            f"{row['rank']}\t{row['pass_gate']}\t{row['score']}\t{row['system']}\t"
            f"{row['median_iptm']}\t{row['median_min_incident_iptm']}\t"
            f"{row['median_contact_pair_pae']}\t{row['median_min_ova_plddt']}\t"
            f"{row['median_min_module_plddt']}\t{row['module_connected_fraction']}\t"
            f"{row['median_interface_jaccard']}"
        )


if __name__ == "__main__":
    main()
