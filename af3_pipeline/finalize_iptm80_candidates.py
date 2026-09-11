#!/usr/bin/env python3
"""Package all strict AF3 ipTM>0.80 hits and plot their seed robustness."""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import shutil
from collections import defaultdict

import matplotlib.pyplot as plt


def read_fasta(path: pathlib.Path) -> str:
    chunks: list[str] = []
    active = False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if active:
                break
            active = True
        elif active:
            chunks.append(line.strip())
    return "".join(chunks)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--seed1-experiment", type=pathlib.Path, action="append", required=True,
        help="one or more seed-1 screens referenced by the validation tables",
    )
    parser.add_argument("--validation", type=pathlib.Path, nargs="+", required=True)
    parser.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--top", type=int, default=0, help="keep only the top N strict hits; 0 keeps all")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    model_dir = args.out / "best_af3_models"
    model_dir.mkdir(exist_ok=True)

    reference = read_fasta(args.reference_fasta)
    hits: list[dict[str, object]] = []
    sources: dict[str, tuple[pathlib.Path, str, pathlib.Path, str]] = {}
    for validation_path in args.validation:
        for row in read_tsv(validation_path):
            if row["passes_final"] != "1":
                continue
            hits.append(dict(row))
            seed1_path = pathlib.Path(row.get("seed1_experiment") or args.seed1_experiment[0])
            multiseed_path = pathlib.Path(
                row.get("multiseed_experiment") or validation_path.parent
            )
            sources[row["sequence"]] = (
                multiseed_path, row["system"], seed1_path, row["seed1_system"]
            )
    hits.sort(key=lambda row: float(row["mean_iptm"]), reverse=True)
    if args.top:
        hits = hits[: args.top]

    seed1_by_system: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for experiment in args.seed1_experiment:
        for row in read_tsv(experiment / "model_scores.tsv"):
            seed1_by_system[(str(experiment.resolve()), row["system"])].append(row)

    final_rows: list[dict[str, object]] = []
    all_models: list[dict[str, object]] = []
    seed_rows: list[dict[str, object]] = []
    mutation_interface_rows: list[dict[str, object]] = []
    fasta: list[str] = []
    for rank, hit in enumerate(hits, 1):
        sequence = str(hit["sequence"])
        nmut = int(hit["n_mutations"])
        label = f"OVA-C4-IPTM80-{nmut:02d}{chr(64 + rank)}"
        multiseed_exp, multiseed_system, seed1_exp, seed1_system = sources[sequence]
        models = [dict(row) for row in seed1_by_system[(str(seed1_exp.resolve()), seed1_system)]]
        models += [dict(row) for row in read_tsv(multiseed_exp / "model_scores.tsv")
                   if row["system"] == multiseed_system]
        if len(models) != 25:
            raise ValueError(f"{label}: expected 25 models, found {len(models)}")
        mutations = [f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference, sequence), 1) if a != b]
        if len(mutations) != nmut:
            raise ValueError(f"{label}: mutation count mismatch")

        by_seed: dict[int, list[float]] = defaultdict(list)
        interface_model_counts: dict[int, int] = defaultdict(int)
        for row in models:
            seed = int(row["model"].split("seed-")[1].split("_")[0])
            by_seed[seed].append(float(row["iptm"]))
            for position in json.loads(row["interface_residues_json"]):
                interface_model_counts[int(position)] += 1
            all_models.append({
                "candidate": label,
                "source_system": seed1_system,
                "seed": seed,
                "model": row["model"],
                "iptm": row["iptm"],
                "min_incident_iptm": row["min_incident_iptm"],
                "contact_pair_pae_median": row["contact_pair_pae_median"],
                "full_connected": row["full_connected"],
                "has_clash": row["has_clash"],
                "cif_path": row["cif_path"],
                "summary_path": row["summary_path"],
            })
        for mutation in mutations:
            position = int(mutation[1:-1])
            mutation_interface_rows.append({
                "candidate": label,
                "mutation": mutation,
                "position": position,
                "interface_models": interface_model_counts[position],
                "interface_recurrence": interface_model_counts[position] / 25,
            })
        for seed in sorted(by_seed):
            values = by_seed[seed]
            seed_rows.append({
                "candidate": label,
                "seed": seed,
                "n_models": len(values),
                "mean_iptm": sum(values) / len(values),
                "min_iptm": min(values),
                "max_iptm": max(values),
            })

        best = max(models, key=lambda row: float(row["iptm"]))
        best_cif = pathlib.Path(best["cif_path"])
        best_summary = pathlib.Path(best["summary_path"])
        shutil.copy2(best_cif, model_dir / f"{label}_best_AF3.cif")
        shutil.copy2(best_summary, model_dir / f"{label}_best_AF3_summary_confidences.json")
        final_rows.append({
            "rank": rank,
            "candidate": label,
            "source_system": seed1_system,
            "n_mutations": nmut,
            "mutations": ",".join(mutations),
            "n_models": 25,
            "mean_iptm": float(hit["mean_iptm"]),
            "median_iptm": float(hit["median_iptm"]),
            "min_iptm": float(hit["min_iptm"]),
            "mean_weakest_chain": float(hit["mean_weakest_chain"]),
            "full_connected_models": hit["full_connected_models"],
            "unclashed_models": hit["unclashed_models"],
            "no_interchain_ss_models": hit["no_interchain_ss_models"],
            "chemically_clean_models": hit["chemically_clean_models"],
            "all_native_c74_c121_models": hit["all_native_c74_c121_models"],
            "siinfecl_intact": int(sequence[257:265] == "SIINFEKL"),
            "best_model_iptm": float(best["iptm"]),
            "sequence": sequence,
        })
        fasta += [f">{label} source={seed1_system} nmut={nmut} mean_AF3_ipTM={float(hit['mean_iptm']):.4f}", sequence]

    final_fields = [
        "rank", "candidate", "source_system", "n_mutations", "mutations", "n_models",
        "mean_iptm", "median_iptm", "min_iptm", "mean_weakest_chain",
        "full_connected_models", "unclashed_models", "no_interchain_ss_models",
        "chemically_clean_models", "all_native_c74_c121_models", "siinfecl_intact",
        "best_model_iptm", "sequence",
    ]
    write_tsv(args.out / "final_candidates.tsv", final_rows, final_fields)
    write_tsv(args.out / "per_seed_summary.tsv", seed_rows,
              ["candidate", "seed", "n_models", "mean_iptm", "min_iptm", "max_iptm"])
    write_tsv(args.out / "per_model_scores.tsv", all_models,
              ["candidate", "source_system", "seed", "model", "iptm", "min_incident_iptm",
               "contact_pair_pae_median", "full_connected", "has_clash", "cif_path", "summary_path"])
    write_tsv(args.out / "mutation_interface_recurrence.tsv", mutation_interface_rows,
              ["candidate", "mutation", "position", "interface_models", "interface_recurrence"])
    (args.out / "final_candidates.fasta").write_text("\n".join(fasta) + "\n")

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for index, row in enumerate(final_rows):
        label = str(row["candidate"])
        points = [seed for seed in seed_rows if seed["candidate"] == label]
        ax.plot([int(p["seed"]) for p in points], [float(p["mean_iptm"]) for p in points],
                marker="o", linewidth=2, label=f"{label} ({row['mean_iptm']:.4f})",
                color=colors[index % len(colors)])
        ax.axhline(float(row["mean_iptm"]), color=colors[index % len(colors)],
                   linewidth=1, alpha=0.25)
    ax.axhline(0.80, color="black", linestyle="--", linewidth=1.5, label="25-model threshold")
    ax.set(xlabel="AF3 random seed", ylabel="Mean ipTM (5 samples)", xticks=[1, 2, 3, 4, 5],
           ylim=(0.70, 0.86), title="Final pure-OVA noncovalent C4 candidates")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(args.out / "af3_seed_robustness.png", dpi=220)
    fig.savefig(args.out / "af3_seed_robustness.svg")
    print(f"wrote {len(final_rows)} candidates to {args.out}")


if __name__ == "__main__":
    main()
