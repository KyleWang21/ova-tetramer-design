#!/usr/bin/env python3
"""Combine an AF3 seed-1 screen with seeds 2--5 validation by sequence."""

from __future__ import annotations

import argparse
import csv
import pathlib
import statistics


def table(path: pathlib.Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t"))


def by_system(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        out.setdefault(row["system"], []).append(row)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--seed1-experiment", type=pathlib.Path, action="append", required=True,
        help="one or more seed-1 screens; sequences must be unique across them",
    )
    ap.add_argument("--multiseed-experiment", type=pathlib.Path, required=True)
    ap.add_argument(
        "--seed1-strict-per-model", type=pathlib.Path, action="append", default=[],
        help="frozen 2.4-A coordinate tables for the seed-1 experiment(s)",
    )
    ap.add_argument(
        "--multiseed-strict-per-model", type=pathlib.Path, action="append", default=[],
        help="frozen 2.4-A coordinate tables for seeds 2--5",
    )
    ap.add_argument("--reference-fasta", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()

    reference = ""
    active = False
    for line in args.reference_fasta.read_text().splitlines():
        if line.startswith(">"):
            if active:
                break
            active = line[1:].startswith("D0-P3-13R")
        elif active:
            reference += line.strip()

    # Some prepared multiseed experiments contain one manifest row per shard,
    # with the same system/sequence repeated for seeds 2--5.  Collapse those
    # rows before combining so one biological candidate yields one summary.
    mm_raw = table(args.multiseed_experiment / "manifest.tsv")
    mm_by_candidate: dict[tuple[str, str], dict[str, str]] = {}
    for row in mm_raw:
        key = (row["name"], row["sequence"])
        mm_by_candidate.setdefault(key, row)
    mm = list(mm_by_candidate.values())
    seed1_source: dict[str, tuple[pathlib.Path, str]] = {}
    score1: dict[tuple[str, str], list[dict[str, str]]] = {}
    chem1: dict[tuple[str, str], list[dict[str, str]]] = {}
    for experiment in args.seed1_experiment:
        key = str(experiment.resolve())
        for row in table(experiment / "manifest.tsv"):
            previous = seed1_source.get(row["sequence"])
            if previous and previous != (experiment, row["name"]):
                raise SystemExit(f"duplicate seed-1 sequence across experiments: {row['name']}")
            seed1_source[row["sequence"]] = (experiment, row["name"])
        for system, rows in by_system(table(experiment / "model_scores.tsv")).items():
            score1[(key, system)] = rows
        for system, rows in by_system(table(experiment / "disulfide_model_scores.tsv")).items():
            chem1[(key, system)] = rows
    scorem = by_system(table(args.multiseed_experiment / "model_scores.tsv"))
    chemm = by_system(table(args.multiseed_experiment / "disulfide_model_scores.tsv"))
    strict1: dict[tuple[str, str], list[dict[str, str]]] = {}
    for path in args.seed1_strict_per_model:
        parent_key = ""
        # Match the table to its owning seed-1 experiment by checking systems.
        rows_by_system = by_system(table(path))
        for experiment in args.seed1_experiment:
            manifest_systems = {row["name"] for row in table(experiment / "manifest.tsv")}
            if manifest_systems & set(rows_by_system):
                parent_key = str(experiment.resolve())
                break
        if not parent_key:
            raise SystemExit(f"cannot match strict seed-1 table to an experiment: {path}")
        for system, rows in rows_by_system.items():
            strict1.setdefault((parent_key, system), []).extend(rows)
    strictm: dict[str, list[dict[str, str]]] = {}
    for path in args.multiseed_strict_per_model:
        for system, rows in by_system(table(path)).items():
            strictm.setdefault(system, []).extend(rows)
    output: list[dict[str, object]] = []
    for row in mm:
        sequence = row["sequence"]
        if sequence not in seed1_source:
            raise SystemExit(f"no seed-1 sequence match for {row['name']}")
        old_experiment, old = seed1_source[sequence]
        source_key = str(old_experiment.resolve())
        scores = score1.get((source_key, old), []) + scorem.get(row["name"], [])
        chemistry = chem1.get((source_key, old), []) + chemm.get(row["name"], [])
        strict = strict1.get((source_key, old), []) + strictm.get(row["name"], [])
        iptm = [float(x["iptm"]) for x in scores]
        weakest = [float(x["min_incident_iptm"]) for x in scores]
        nmut = sum(a != b for a, b in zip(reference, sequence))
        output.append({
            "system": row["name"],
            "multiseed_experiment": str(args.multiseed_experiment),
            "seed1_system": old,
            "seed1_experiment": str(old_experiment),
            "n_mutations": nmut,
            "n_models": len(scores),
            "mean_iptm": statistics.fmean(iptm),
            "median_iptm": statistics.median(iptm),
            "min_iptm": min(iptm),
            "mean_weakest_chain": statistics.fmean(weakest),
            "full_connected_models": sum(int(x["full_connected"]) for x in scores),
            "strict_geometry_models": len(strict),
            "unclashed_models": sum(int(x["atomic_clash_count_2p4"]) == 0 for x in strict),
            "zero_ca_overlap_models": sum(int(x["severe_ca_overlap_count_2p5"]) == 0 for x in strict),
            "valid_c4_network_models": sum(int(x["valid_c4_network"]) for x in strict),
            "two_interfaces_design_ge3_models": sum(
                int(x["two_interfaces_each_design_ge3"]) for x in strict
            ),
            "no_interchain_ss_models": sum(int(x["no_interchain_ss"]) for x in chemistry),
            "chemically_clean_models": sum(int(x["noncovalent_chemically_clean"]) for x in chemistry),
            "all_native_c74_c121_models": sum(int(x["native_c74_c121_bonds"]) == 4 for x in chemistry),
            "passes_final": int(
                len(scores) == 25
                and statistics.fmean(iptm) > 0.80
                and nmut < 25
                and len(strict) == 25
                and all(int(x["atomic_clash_count_2p4"]) == 0 for x in strict)
                and all(int(x["severe_ca_overlap_count_2p5"]) == 0 for x in strict)
                and all(int(x["valid_c4_network"]) for x in strict)
                and all(int(x["native_c74_c121_bonds"]) == 4 for x in strict)
                and all(int(x["interchain_ss_bonds"]) == 0 for x in strict)
                and all(int(x["two_interfaces_each_design_ge3"]) for x in strict)
                and len(chemistry) == 25
                and all(int(x["no_interchain_ss"]) for x in chemistry)
            ),
            "sequence": sequence,
        })
    output.sort(key=lambda x: (-int(x["passes_final"]), -float(x["mean_iptm"])))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(output[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    for row in output:
        print(
            row["system"], f"models={row['n_models']}", f"mean={float(row['mean_iptm']):.4f}",
            f"nmut={row['n_mutations']}", f"pass={row['passes_final']}",
        )


if __name__ == "__main__":
    main()
